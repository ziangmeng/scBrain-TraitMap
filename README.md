# Overview

**scBrain-TraitMap** is a streamlined and reproducible framework for integrating **single-cell RNA-seq (scRNA-seq)**, **single-cell ATAC-seq (scATAC-seq)**, and **regulatory guidance graphs** into a unified analytical pipeline. Built on the SCGLUE model, the framework aligns transcriptomic and epigenomic cells into a shared latent space, enabling subtype transfer, peak-level annotation, and downstream trait–cell-type mapping.

![Image text](fig1.png)

The framework provides:

- A script-based, modular pipeline with a clean directory structure (`./data`, `./model`, `./logs`).
- Automated preprocessing, PCA/LSI computation, and harmonized cell-type annotations.
- Cross-modality integration of scRNA-seq and scATAC-seq using SCGLUE with a user-supplied guidance graph.
- KNN-based label transfer from RNA to ATAC cells within the shared embedding.
- Extraction of cell-type–specific accessible peaks for S-LDSC enrichment analysis.
- Easy hyperparameter tuning through both command-line arguments and in-script defaults.
- Reproducible execution with transparent logging and version control.

scBrain-TraitMap is designed for researchers aiming to study cell-type–resolved regulatory landscapes and to map complex-trait associations onto specific neural populations in the human brain.


## 1. Building the Guidance Graph

This notebook generates the *guidance graph* required for SCGLUE-based cross-modality integration.  
The graph encodes regulatory relationships between genes (from scRNA-seq) and chromatin-accessible regions (from scATAC-seq), and serves as the structural prior that anchors the two modalities in a shared latent space.

The workflow includes:

- Loading subsampled scRNA-seq and scATAC-seq datasets (`rna_subsampled.h5ad`, `atac_subsampled.h5ad`)
- Extracting gene coordinates using a GTF annotation (lifted to hg19)
- Parsing ATAC peaks into chromosome-based genomic intervals
- Constructing an RNA-anchored guidance graph using SCGLUE’s regulatory rules
- Sanitizing metadata to ensure H5AD compatibility
- Exporting the graph in both `.graphml.gz` (human-readable) and `.gpickle` (Python-native) formats

This guidance graph is the only structural prior required by the SCGLUE training script (`2.run_scglue.py`), and defines how genes and peaks are linked during cross-modality alignment.

## 2. SCGLUE Cross-Modality Integration (`2.run_scglue.py`)

This script performs the full SCGLUE workflow to integrate subsampled scRNA-seq and scATAC-seq profiles into a shared latent space. It uses the guidance graph generated in Step 1 to anchor genes and chromatin peaks through regulatory relationships. Users may either run the script with its default settings or override any hyperparameter using command-line options.

### Workflow summary

This script implements the following steps:

1. **Load input data**
   - `rna_subsampled.h5ad`
   - `atac_subsampled.h5ad`
   - `guidance.graphml.gz`

2. **Graph preprocessing**
   - Add self-loops to all nodes
   - Standardize edge attributes (`weight`, `sign`)

3. **Metadata cleanup**
   - Validate finite values in `.X`
   - Normalize category fields
   - Ensure the presence of `highly_variable` features

4. **Feature representation**
   - PCA on RNA (default: 50 components)
   - LSI on ATAC (default: 50 components)
   - Optional filtering of ATAC peaks to those in the guidance graph

5. **Cell-type harmonization**
   - Map fine-grained RNA and ATAC labels to broader categories
   - Store harmonized labels in `glue_cell_type`

6. **Configure SCGLUE datasets**
   - RNA modeled using a Normal distribution
   - ATAC modeled using a Zero-Inflated Log-Normal distribution
   - Use PCA/LSI representations and HV features

7. **Train SCGLUE**
   - Fully configurable training schedule:
     - `--max_epochs` (default: 100)
     - `--batch_size` (default: 256)
     - `--align_burnin` (default: 30)
     - `--graph_batch` (default: 7500)
   - All logs are written to `./logs/`

8. **Embedding and visualization**
   - Encode RNA and ATAC into `X_glue`
   - Compute joint UMAP layout
   - Save RNA-only and ATAC-only UMAP images

9. **Save outputs**
   - Trained SCGLUE model (`.dill`)
   - Embedded `.h5ad` files
   - UMAP visualizations
   - Reproducibility log files

### Running the script

Run with default hyperparameters:

```bash
python 2.run_scglue.py
```

Customize hyperparameters:
```bash
python 2.run_scglue.py \
    --n_comps 64 \
    --max_epochs 150 \
    --batch_size 128 \
    --graph_batch 5000 \
    --align_burnin 20
```
### Hyperparameter reference

| Parameter         | CLI Flag          | Default | Description |
|------------------|-------------------|---------|-------------|
| `n_comps`        | `--n_comps`       | 50      | Number of PCA components for RNA and LSI components for ATAC. Controls dimensionality of the initial feature representation. |
| `max_epochs`     | `--max_epochs`    | 100     | Maximum number of training epochs for SCGLUE. |
| `batch_size`     | `--batch_size`    | 256     | Mini-batch size for training RNA/ATAC data batches. |
| `graph_batch`    | `--graph_batch`   | 7500    | Number of graph edges used per graph batch during SCGLUE graph training. |
| `patience`       | `--patience`      | 20      | Number of epochs without improvement before early stopping is triggered. |
| `lr_patience`    | `--lr_patience`   | 15      | Number of stagnant epochs before reducing the learning rate. |
| `align_burnin`   | `--align_burnin`  | 30      | Number of epochs during which the model focuses on learning cell-type structure before enforcing cross-modality alignment. |



Output structure
```bash
model/
    glue_run_YYYYMMDD-HHMMSS/
        glue_model_*.dill
        umap_RNA_only_*.png
        umap_ATAC_only_*.png
        rna_with_Xglue_*.h5ad
        atac_with_Xglue_*.h5ad

logs/
    glue_YYYYMMDD-HHMMSS.log
```
