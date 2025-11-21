# Overview

**scBrain-TraitMap** is a lightweight and reproducible computational framework designed to integrate **single-cell RNA-seq (scRNA-seq)**, **single-cell ATAC-seq (scATAC-seq)**, and **genomic regulatory graphs** to generate a harmonized embedding of brain cell populations. Built upon the SCGLUE model, the pipeline enables the construction of cross-modality embeddings, joint visualization, and downstream trait-mapping analyses in human brain datasets.
![Image text](fig1.png)
The framework provides:

- A fully script-based pipeline with clear input/output structure (`./data`, `./model`, `./logs`).
- Automated preprocessing, feature selection, PCA/LSI computation, and harmonized cell-type assignment.
- Graph-guided integration of RNA and ATAC modalities using SCGLUE.
- Consistent UMAP embeddings for both modalities in shared space.
- Easy customization of training hyperparameters through command-line options.
- Journal-ready reproducibility with transparent logging and version-freezing.

scBrain-TraitMap is intended for researchers studying cell-type-specific regulatory architecture, multimodal single-cell integration, and the mapping of complex trait signals to brain cellular contexts.
