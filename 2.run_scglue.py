from __future__ import annotations 
from pathlib import Path
import argparse
import os
import warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import scglue
import networkx as nx
import matplotlib.pyplot as plt
import re
import pytorch_lightning as pl

# ============================================================
# Default hyperparameters
# ============================================================

DEFAULT_N_COMPS      = 50
DEFAULT_MAX_EPOCHS   = 100
DEFAULT_BATCH_SIZE   = 256
DEFAULT_GRAPH_BATCH  = 7500
DEFAULT_PATIENCE     = 20
DEFAULT_LR_PATIENCE  = 15
DEFAULT_ALIGN_BURNIN = 30

# ============================================================
# Argument parser
# ============================================================

parser = argparse.ArgumentParser(description="SCGLUE training pipeline")

parser.add_argument("--n_comps",      type=int, default=DEFAULT_N_COMPS)
parser.add_argument("--max_epochs",   type=int, default=DEFAULT_MAX_EPOCHS)
parser.add_argument("--batch_size",   type=int, default=DEFAULT_BATCH_SIZE)
parser.add_argument("--graph_batch",  type=int, default=DEFAULT_GRAPH_BATCH)
parser.add_argument("--patience",     type=int, default=DEFAULT_PATIENCE)
parser.add_argument("--lr_patience",  type=int, default=DEFAULT_LR_PATIENCE)
parser.add_argument("--align_burnin", type=int, default=DEFAULT_ALIGN_BURNIN)

args = parser.parse_args()

# ============================================================
# Environment
# ============================================================

os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
warnings.filterwarnings("ignore", message="divide by zero")
sc.settings.verbosity = 2

BASE_DIR  = Path(".")
DATA_DIR  = BASE_DIR / "data"
MODEL_DIR = BASE_DIR / "model"
LOG_DIR   = BASE_DIR / "logs"

MODEL_DIR.mkdir(exist_ok=True)
LOG_DIR.mkdir(exist_ok=True)

# ------------------------
# FIXED RUN DIRECTORY
# ------------------------
RUN_DIR = MODEL_DIR / "glue_run"
RUN_DIR.mkdir(parents=True, exist_ok=True)

# ------------------------
# FIXED LOG FILE
# ------------------------
LOG_FILE = LOG_DIR / "glue.log"

def log(msg: str):
    line = f"[LOG] {msg}"
    print(line)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")

log("Start SCGLUE run")

# ============================================================
# Paths
# ============================================================

rna_path  = DATA_DIR / "rna_subsampled.h5ad"
atac_path = DATA_DIR / "atac_subsampled.h5ad"
graph_gml = DATA_DIR / "guidance.graphml.gz"

# ============================================================
# Load data
# ============================================================

log("Loading RNA, ATAC and guidance graph")
rna  = ad.read_h5ad(rna_path)
atac = ad.read_h5ad(atac_path)
G    = nx.read_graphml(graph_gml)
log(f"RNA: {rna.shape}, ATAC: {atac.shape}, Graph: {G.number_of_nodes()} nodes / {G.number_of_edges()} edges")

# ============================================================
# Graph preprocessing
# ============================================================

def add_self_loops_and_attrs(g: nx.Graph):
    for n in list(g.nodes):
        if not g.has_edge(n, n):
            g.add_edge(n, n, weight=1.0, sign=1)
    for _, _, a in g.edges(data=True):
        a.setdefault("weight", 1.0)
        a.setdefault("sign", 1)

add_self_loops_and_attrs(G)

# ============================================================
# Category sanitization
# ============================================================

def sanitize_categories(adata: ad.AnnData):
    for attr in ("obs", "var"):
        df = getattr(adata, attr)
        for c in df.columns:
            if str(df[c].dtype) == "category":
                df[c] = df[c].astype(str)

sanitize_categories(rna)
sanitize_categories(atac)

# ============================================================
# Representation (PCA / LSI)
# ============================================================

def ensure_representation(rna: ad.AnnData, atac: ad.AnnData, n_comps=50):

    if "X_pca" not in rna.obsm:
        if sp.issparse(rna.X):
            rna.X = rna.X.toarray().astype(np.float32)
        sc.pp.highly_variable_genes(rna, n_top_genes=min(4000, rna.n_vars), flavor="seurat_v3")
        sc.pp.pca(rna, n_comps=n_comps, svd_solver="arpack", use_highly_variable=True)
        log(f"RNA PCA computed")
    else:
        log("RNA: Using cached PCA")

    if "X_lsi" not in atac.obsm:
        scglue.data.lsi(atac, n_components=n_comps, n_iter=15)
        log("ATAC LSI computed")
    else:
        log("ATAC: Using cached LSI")

# ============================================================
# Subset ATAC peaks to graph
# ============================================================

def subset_atac_to_graph(atac: ad.AnnData, g: nx.Graph) -> ad.AnnData:
    peaks = list(set(atac.var_names) & set(g.nodes))
    log(f"ATAC peaks before: {atac.n_vars}, after: {len(peaks)}")
    atac_sub = atac[:, peaks].copy()
    if "X_lsi" in atac_sub.obsm:
        del atac_sub.obsm["X_lsi"]
    scglue.data.lsi(atac_sub, n_components=args.n_comps, n_iter=15)
    return atac_sub

atac = subset_atac_to_graph(atac, G)

# ============================================================
# Cell type harmonization
# ============================================================

rna_map = {
    "EN-ET": "EN",
    "EN-IT-1": "EN",
    "EN-IT-2": "EN",
    "IN": "IN",
    "Astro": "Astro",
    "OPC": "OPC",
    "EC": "EC"
}

rna.obs["glue_cell_type"] = rna.obs["H1_annotation"].map(rna_map).fillna("Unknown")

# ATAC mapping rules
rules = {
    "IN": [r"^VIP", r"^SST", r"^PVALB", r"^LAMP5", r"^SNCG", r"^PV_?ChC", r"^SST_?CHODL"],
    "EN": [r"^IT", r"^CT", r"^FOXP2", r"^ERC", r"^PIR", r"^PER", r"^SUB", r"^AMY"],
    "Astro": [r"^ASC"],
    "OPC": [r"^OPC"],
    "EC": [r"^EC"],
}
compiled = {k: [re.compile(p) for p in v] for k, v in rules.items()}

def map_to_glue(label: str) -> str:
    if not isinstance(label, str): return "Unknown"
    for broad, pats in compiled.items():
        for pat in pats:
            if pat.search(label): return broad
    return "Unknown"

atac.obs["glue_cell_type"] = atac.obs["Cell_type (HSC)"].astype(str).map(map_to_glue)

log(str(rna.obs["glue_cell_type"].value_counts()))
log(str(atac.obs["glue_cell_type"].value_counts()))

# ============================================================
# HV features
# ============================================================

if "highly_variable" not in rna.var.columns:
    rna.var["highly_variable"] = True

if "highly_variable" not in atac.var.columns:
    counts = atac.layers["counts_round"] if "counts_round" in atac.layers else (
        atac.layers["counts"] if "counts" in atac.layers else atac.X
    )
    if sp.issparse(counts):
        nnz = np.diff(counts.tocsc().indptr)
    else:
        nnz = (np.asarray(counts) > 0).sum(axis=0)

    rate = nnz / atac.n_obs
    top = np.argpartition(rate, -200_000)[-200_000:]
    mask = np.zeros(atac.n_vars, bool)
    mask[top] = True
    atac.var["highly_variable"] = mask

# ============================================================
# Representation
# ============================================================

ensure_representation(rna, atac, n_comps=args.n_comps)

# ============================================================
# Configure for SCGLUE
# ============================================================

scglue.models.configure_dataset(
    rna, "Normal",
    use_highly_variable=True,
    use_rep="X_pca",
    use_cell_type="glue_cell_type"
)

use_atac_counts = (
    "counts_round" if "counts_round" in atac.layers else
    "counts" if "counts" in atac.layers else None
)

scglue.models.configure_dataset(
    atac, "ZILN",
    use_highly_variable=True,
    use_layer=use_atac_counts,
    use_rep="X_lsi",
    use_cell_type="glue_cell_type"
)

# ============================================================
# Train model
# ============================================================

pl.seed_everything(0)

log("Fitting SCGLUE")
glue = scglue.models.fit_SCGLUE(
    {"rna": rna, "atac": atac},
    G,
    fit_kws={
        "directory": str(RUN_DIR),
        "max_epochs": args.max_epochs,
        "patience": args.patience,
        "reduce_lr_patience": args.lr_patience,
        "graph_batch_size": args.graph_batch,
        "data_batch_size": args.batch_size,
        "align_burnin": args.align_burnin
    }
)
log("Training finished")

# ============================================================
# Save model (NO DATE)
# ============================================================

model_path = RUN_DIR / "glue_model.dill"
glue.save(str(model_path))
log(f"Saved model: {model_path}")

# ============================================================
# Encode embeddings
# ============================================================

rna.obsm["X_glue"]  = glue.encode_data("rna", rna)
atac.obsm["X_glue"] = glue.encode_data("atac", atac)

rna.obs["domain"]  = "RNA"
atac.obs["domain"] = "ATAC"

combined = ad.concat([rna, atac], join="outer")

sc.pp.neighbors(combined, use_rep="X_glue", metric="cosine")
sc.tl.umap(combined)

# ============================================================
# Save UMAP figures (NO DATE)
# ============================================================

fig_rna = RUN_DIR / "umap_rna.png"
sc.pl.umap(
    combined[combined.obs["domain"]=="RNA"],
    color=["glue_cell_type"],
    show=False, frameon=False
)
plt.savefig(fig_rna, dpi=300)
plt.close()
log(f"Saved RNA UMAP: {fig_rna}")

fig_atac = RUN_DIR / "umap_atac.png"
sc.pl.umap(
    combined[combined.obs["domain"]=="ATAC"],
    color=["glue_cell_type"],
    show=False, frameon=False
)
plt.savefig(fig_atac, dpi=300)
plt.close()
log(f"Saved ATAC UMAP: {fig_atac}")

# ============================================================
# Save AnnData (NO DATE)
# ============================================================

rna_out  = RUN_DIR / "rna_with_Xglue.h5ad"
atac_out = RUN_DIR / "atac_with_Xglue.h5ad"

rna.write_h5ad(str(rna_out), compression="gzip")
atac.write_h5ad(str(atac_out), compression="gzip")

log("Saved all outputs. Pipeline completed.")
