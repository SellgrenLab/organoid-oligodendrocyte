
# Organoid Oligodendrocyte Analysis

This repository contains code and resources for the analysis of single-cell RNA sequencing data from human brain organoids enriched for oligodendrocytes. 

---
## Article
## Human oligodendrocyte progenitor cells mediate synapse elimination through TAM receptor activation
Asimenia Gkogka, Susmita Malwade, Marja Koskuvi, Raj Bose, Sandra Ceccatelli, Jari Koistinaho, Jari Tiihonen, Martin Schalling, Samudyata*, Carl M. Sellgren*


### Abstract
Oligodendrocyte progenitor cells (OPCs) have been implicated in synaptic remodelling in animal models, but the underlying mechanisms and their relevance to human brain development remain unclear. Here, we generate a human multi-lineage forebrain organoid model in which OPCs, together with microglia, form close contacts with synapses and spontaneously internalize synaptic material. Single-nucleus transcriptomic profiling with unbiased cell-cell communication analysis identifies the growth arrest-specific gene 6 (GAS6)-TYRO3, AXL, and MERTK (TAM) receptor axis as a key signalling pathway, with neurons and microglia expressing GAS6 and a subset of OPCs expressing AXL. Further, dose-dependent pharmacological inhibition of TAM receptors demonstrates the importance of AXL, and targeted reduction of AXL expression in OPCs impairs synaptic uptake. These findings reveal a role for GAS6-AXL signalling in driving synaptic internalisation by AXL+ OPCs during early human brain development.

---
## Key Scripts

File	Description

- **`Coexp_modules.R`** - 	Detects and visualizes co-expression modules using WGCNA
- **`crosstalk.R`** - 	Analyzes cell-cell communication networks using CellChat
- **`disease_enrichment.R`** - 	Performs enrichment analysis against disease gene sets
- **`plots.R`** - 	Generates figures for the manuscript
- **`qc_dimred_integration.Rmd`** - 	Performs quality control and dimensionality reduction using Seurat
- **`integrated_velocity.ipynb`** - 	Calculates and visualizes RNA velocity using scVelo

Note that the scripts in the main branch above are to be run independently. To run it linearly access the dev-co branch.

## Citation
[bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2023.09.04.556176v1)
[Nature Communications article](https://www.nature.com/articles/s41467-025-66521-1)



