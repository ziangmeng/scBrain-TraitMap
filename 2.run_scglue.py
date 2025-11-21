from __future__ import annotations
from pathlib import Path
import argparse
import os
import warnings
import datetime as dt
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
# Default hyperparameters (can be edited directly)
# ============================================================

DEFAULT_N_COMPS      = 50
DEFAULT_MAX_EPOCHS   = 100
DEFAULT_BATCH_SIZE   = 256
DEFAULT_GRAPH_BATCH  = 7500
DEFAULT_PATIENCE     = 20
DEFAULT_LR_PATIENCE  = 15
DEFAULT_ALIGN_BURNIN = 30

# ============================================================
# Argument parser (CLI can override the defaults)
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

TS = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
RUN_DIR = MODEL_DIR / f"glue_run_{TS}"
RUN_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# Logging
# ============================================================

LOG_FILE = LOG_DIR / f"glue_{TS}.log"


def log(msg: str):
    line = f"[{dt.datetime.now().strftime('%H:%M:%S')}] {msg}"
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
# Data loading
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
        if "weight" not in a:
            a["weight"] = 1.0
        if "sign" not in a:
            a["sign"] = 1


add_self_loops_and_attrs(G)

# ============================================================
# Utilities
# ============================================================


def assert_finite_X(adata, name="adata"):
    X = adata.X
    if sp.issparse(X):
        bad = (~np.isfinite(X.data)).sum()
    else:
        bad = (~np.isfinite(np.asarray(X))).sum()
    if bad:
        raise ValueError(f"{name}.X has {bad} non-finite values")
    else:
        log(f"[CHECK] {name}.X all finite")


def sanitize_categories(adata: ad.AnnData):
    for attr in ("obs", "var"):
        df = getattr(adata, attr)
        for c in df.columns:
            if str(df[c].dtype) == "category":
                df[c] = df[c].astype(str)


sanitize_categories(rna)
sanitize_categories(atac)

# ============================================================
# Representation (PCA for RNA, LSI for ATAC)
# ============================================================


def ensure_representation(rna: ad.AnnData, atac: ad.AnnData, n_comps=50):
    # RNA PCA
    if "X_pca" not in rna.obsm or rna.obsm["X_pca"].shape[1] != n_comps:
        if sp.issparse(rna.X):
            rna.X = rna.X.toarray().astype(np.float32)
        assert_finite_X(rna, "RNA")
        if "highly_variable" not in rna.var.columns:
            sc.pp.highly_variable_genes(
                rna,
                n_top_genes=min(4000, rna.n_vars),
                flavor="seurat_v3"
            )
        sc.pp.pca(
            rna,
            n_comps=n_comps,
            svd_solver="arpack",
            use_highly_variable=True
        )
        log(f"RNA: PCA ready (n_comps={n_comps})")
    else:
        log("RNA: use cached X_pca")

    # ATAC LSI
    if "X_lsi" not in atac.obsm or atac.obsm["X_lsi"].shape[1] != n_comps:
        scglue.data.lsi(atac, n_components=n_comps, n_iter=15)
        log(f"ATAC: LSI ready (n_components={n_comps})")
    else:
        log("ATAC: use cached X_lsi")


# ============================================================
# Subset ATAC peaks to graph (memory saving)
# ============================================================


def subset_atac_to_graph(atac: ad.AnnData, g: nx.Graph) -> ad.AnnData:
    peaks_in_graph = list(set(atac.var_names) & set(g.nodes))
    log(f"ATAC peaks before: {atac.n_vars}  after keep-in-graph: {len(peaks_in_graph)}")
    atac_sub = atac[:, peaks_in_graph].copy()
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

if "H1_annotation" in rna.obs.columns:
    rna.obs["glue_cell_type"] = rna.obs["H1_annotation"].map(rna_map).fillna("Unknown")
else:
    raise KeyError("RNA.obs is missing 'H1_annotation'")

rules = {
    "IN": [
        r"^VIP", r"^SST", r"^PVALB", r"^LAMP5", r"^SNCG",
        r"^PV_?ChC", r"^SST_?CHODL",
    ],
    "EN": [
        r"^ITL23", r"^ITL4(?!5)",
        r"^ITL5", r"^ITL6", r"^ITL45", r"^ITL34", r"^ITV1C",
        r"^CT", r"^FOXP2", r"^CTXGA", r"^ICGA", r"^ERC", r"^PIR", r"^PER", r"^SUB", r"^AMY",
    ],
    "Astro": [
        r"^ASCT", r"^ASCNT",
    ],
    "OPC": [
        r"^OPC$",
    ],
    "EC": [
        r"^EC$",
    ],
}
compiled = {k: [re.compile(p) for p in v] for k, v in rules.items()}


def map_to_glue(label: str) -> str:
    if not isinstance(label, str):
        return "Unknown"
    for broad, pats in compiled.items():
        for pat in pats:
            if pat.search(label):
                return broad
    return "Unknown"


if "Cell_type (HSC)" in atac.obs.columns:
    atac.obs["glue_cell_type"] = atac.obs["Cell_type (HSC)"].astype(str).map(map_to_glue)
else:
    raise KeyError("ATAC.obs is missing 'Cell_type (HSC)'")

log("[INFO] RNA glue_cell_type distribution:")
log(str(rna.obs["glue_cell_type"].value_counts(dropna=False)))
log("[INFO] ATAC glue_cell_type distribution:")
log(str(atac.obs["glue_cell_type"].value_counts(dropna=False)))

# ============================================================
# Highly variable features (critical for configure_dataset)
# ============================================================

if "highly_variable" not in rna.var.columns:
    rna.var["highly_variable"] = True
    log("RNA: 'highly_variable' not found, set all to True")

if "highly_variable" not in atac.var.columns:
    prefer_counts = "counts_round" if "counts_round" in atac.layers else (
        "counts" if "counts" in atac.layers else None
    )
    Xc = atac.layers[prefer_counts] if prefer_counts else atac.X
    if sp.issparse(Xc):
        nnz_per_peak = np.diff(Xc.tocsc().indptr)
    else:
        Xc = np.asarray(Xc)
        nnz_per_peak = (Xc > 0).sum(axis=0)
    nonzero_rate = nnz_per_peak / atac.n_obs
    K = min(200_000, atac.n_vars)
    top_idx = np.argpartition(nonzero_rate, -K)[-K:]
    hv_mask = np.zeros(atac.n_vars, dtype=bool)
    hv_mask[top_idx] = True
    atac.var["highly_variable"] = hv_mask
    log(f"ATAC highly variable peaks: {hv_mask.sum()} / {atac.n_vars}")

# ============================================================
# Representation (after HV is ensured)
# ============================================================

ensure_representation(rna, atac, n_comps=args.n_comps)

# ============================================================
# Configure SCGLUE datasets (original API)
# ============================================================

log("Configuring datasets for SCGLUE")

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
if use_atac_counts is not None:
    Xcheck = atac.layers[use_atac_counts]
    if sp.issparse(Xcheck):
        if (Xcheck.data < 0).any():
            raise ValueError("ATAC counts layer has negative values")
    else:
        if (np.asarray(Xcheck) < 0).any():
            raise ValueError("ATAC counts layer has negative values")

scglue.models.configure_dataset(
    atac, "ZILN",
    use_highly_variable=True,
    use_layer=use_atac_counts,
    use_rep="X_lsi",
    use_cell_type="glue_cell_type"
)

# ============================================================
# Training SCGLUE
# ============================================================

pl.seed_everything(0)
pl.utilities.rank_zero_info = print

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
# Save model
# ============================================================

model_path = RUN_DIR / f"glue_model_{TS}.dill"
glue.save(str(model_path))
log(f"Saved model: {model_path}")

# ============================================================
# Encode embeddings and compute UMAP
# ============================================================

log("Encoding X_glue")
rna.obsm["X_glue"]  = glue.encode_data("rna", rna)
atac.obsm["X_glue"] = glue.encode_data("atac", atac)
rna.obs["domain"]  = "RNA"
atac.obs["domain"] = "ATAC"

combined = ad.concat(
    [rna, atac],
    join="outer",
    label="mod",
    keys=["rna", "atac"],
    uns_merge="unique"
)

sc.pp.neighbors(combined, use_rep="X_glue", metric="cosine")
sc.tl.umap(combined, min_dist=0.3, spread=1.0)

rna_only  = combined[combined.obs["domain"] == "RNA"].copy()
atac_only = combined[combined.obs["domain"] == "ATAC"].copy()

fig_rna = RUN_DIR / f"umap_RNA_only_glue_cell_type_{TS}.png"
sc.pl.umap(
    rna_only,
    color=["glue_cell_type"],
    na_color="lightgrey",
    frameon=False,
    size=8,
    legend_loc="on data",
    show=False
)
plt.savefig(fig_rna, dpi=300, bbox_inches="tight")
plt.close()
log(f"Saved UMAP (RNA only): {fig_rna}")

fig_atac = RUN_DIR / f"umap_ATAC_only_glue_cell_type_{TS}.png"
sc.pl.umap(
    atac_only,
    color=["glue_cell_type"],
    na_color="lightgrey",
    frameon=False,
    size=8,
    legend_loc="on data",
    show=False
)
plt.savefig(fig_atac, dpi=300, bbox_inches="tight")
plt.close()
log(f"Saved UMAP (ATAC only): {fig_atac}")

# ============================================================
# Save AnnData with embeddings
# ============================================================

rna_out  = RUN_DIR / f"rna_with_Xglue_{TS}.h5ad"
atac_out = RUN_DIR / f"atac_with_Xglue_{TS}.h5ad"
rna.write_h5ad(str(rna_out), compression="gzip")
atac.write_h5ad(str(atac_out), compression="gzip")
log(f"Saved RNA/ATAC with X_glue: {rna_out} | {atac_out}")

log("All done")
