# single-cell-ADOC

This repository contains code for the single-cell RNA-seq analysis of the ADOC (Artificial Dual-channel Organ-on-a-Chip) model of human embryo implantation.

The scripts perform preprocessing, integration, annotation, reference-based transfer, differential expression, and cell–cell communication analyses on 10x Genomics Chromium scRNA-seq data generated from the ADOC platform.

## Study context

The ADOC system is a dual-channel microfluidic model containing organoid-derived endometrial epithelium and primary stromal cells.
It reproduces essential features of the human endometrium, such as epithelial polarization, stromal decidualization, extracellular vesicle release, and hormone-induced receptivity.

The single-cell RNA-seq analysis in this repository focuses on characterizing cellular heterogeneity, transcriptional states, and intercellular communication within this model to elucidate mechanisms of endometrial function and receptivity.

## Code availability

Code used to process the data is published in Zenodo under the DOI: 10.5281/zenodo.17412607 and the GitHub repo: https://github.com/FCS-BioIT/single-cell-ADOC

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17412607.svg)](https://doi.org/10.5281/zenodo.17412607)


## Data availability

Processed data required to reproduce results is published in ArrayExpress under the accession code E-MTAB-15842.

---

## Reproducibility

This repository provides a fully reproducible analysis environment:

- **Pixi**: defines dependencies and R base version. 

---

### 1. Clone the repository

Clone this repository and create the folders input and output in the root folder of the repo

### 2. Activate the Pixi environment

Install the pixi environment using the pixi.toml file, then enter the pixi shell and verify the R installation is in the .pixi folder

```
pixi install
pixi shell -e schandling

which R
```

### 3. Install external dependencies: 

> Install CellChat from github

In R from the pixi shell session, run:

```
devtools::install_github("jinworks/CellChat@623f48f")
```

> Install GenomeInfoDbData (dependency in DE analysis)

In R from the pixi shell session, run:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GenomeInfoDbData", ask = FALSE, update = FALSE)
```


## Analysis workflow

### 1. Raw count matrix preprocessing
Decontamination and filtering of raw matrices using DecontX.

- `src/upstream/decontx_remove_background.R`

### 2. Sample quality control (QC)
Quality control of individual samples and generation of QC plots.

- `src/scrnaseq-preprocess/qc.R`  
- `src/scrnaseq-preprocess/qc_plots.R`  
- `src/scrnaseq-preprocess/utils_qc.R`

### 3. Sample integration
Integration of individual samples into a joint Seurat object.

- `src/scrnaseq-preprocess/sample_integration.R`  
- `src/scrnaseq-preprocess/utils_sample_integration.R`

### 4. Annotation
Cell type annotation and hierarchical classification of clusters.

- `src/annotation/annotation.R`

### 5. Differential expression
Identification of marker genes and pathway enrichment analyses.

- `src/differential_expression/differential_expression.R`  
- `src/differential_expression/de_volcanoplots.R`  
- `src/differential_expression/ora_dotplots_selection.R`  
- `src/differential_expression/ora_gobp.R`  
- `src/differential_expression/ora_kegg.R`
- `src/differential_expression/ora_dotplots_gobp_stromal_paths.R`
- `src/differential_expression/ora_dotplots_gobp.R`
- `src/differential_expression/ora_dotplots_kegg.R`

### 6. Cell–cell communication
Inference of intercellular communication networks using CellChat.

- `src/cellchat/cellchat.R`  
- `src/cellchat/cellchat_grouped.R`  
- `src/cellchat/ccc_infoflow_contrib_heatmaps.R`

### 7. Reference-based integration and similarity analysis
Transfer of cell type labels and evaluation of similarity between the ADOC dataset and a reference dataset using Seurat transfer anchors and cosine similarity.

- `src/convert_h5ad_to_rds.R`
- `src/reference_dataset/load_reference.R`
- `src/reference_dataset/reference_projection.R`

---

## Dependencies
- `envs/pixi.toml`

---

## Contact
Jaime Llera Oyola  
Carlos Simon Foundation  
[jllera@fundacioncarlossimon.com](mailto:jllera@fundacioncarlossimon.com)



