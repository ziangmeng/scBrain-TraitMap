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
