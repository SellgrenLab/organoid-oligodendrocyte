
# 🧠 Organoid Oligodendrocyte Analysis

This repository contains code and resources for the analysis of single-cell RNA sequencing data from human brain organoids enriched for oligodendrocytes. 

---
## 🧬 Manuscript
## Human oligodendrocyte progenitor cells mediate synapse elimination through TAM receptor activation
Asimenia Gkogka, Susmita Malwade, Marja Koskuvi, Raj Bose, Sandra Ceccatelli, Jari Koistinaho, Jari Tiihonen, Martin Schalling, Samudyata*, Carl M. Sellgren*


### Abstract
Recent work in animal models suggests that oligodendrocyte progenitor cells (OPCs) contribute to elimination of synapses in the developing brain. However, a mechanistic understanding of this process is still missing, and it remains uncertain whether human OPCs also display this feature. In this study, we develop a human multi-lineage forebrain organoid model in which OPCs, alongside microglia, exhibit close interactions with synapses and spontaneously internalize synaptic material. Unbiased cell-cell communication analysis based on single-nucleus transcriptomic profiling predicted GAS6-TAM receptor interaction as a major signaling pathway, with neurons and microglia expressing GAS6 ligand and OPCs expressing the TAM receptor AXL. Dose-dependent pharmacological inhibition of TAM receptors then demonstrated the importance of AXL receptor activation, and specific knockdown of AXL in OPCs resulted in impaired uptake of synaptic material. In summary, these data define a role of GAS6-AXL signaling in promoting OPC-mediated internalization of synaptic structures during early brain development. 

---
## 📜 Key Scripts

File	Description
Coexp_modules.R	Detects and visualizes co-expression modules using WGCNA
crosstalk.R	Analyzes cell-cell communication networks using CellChat
disease_enrichment.R	Performs enrichment analysis against disease gene sets
plots.R	Generates figures for the manuscript
qc_dimred_integration.Rmd	Performs quality control and dimensionality reduction using Seurat
integrated_velocity.ipynb	Calculates and visualizes RNA velocity using scVelo

Note that the scripts in the main branch above are to be run independently. To run it linearly access the dev-co branch.

## 🔧 Dependencies

### 📦 R Packages

The R scripts in the `/code` folder depend on the following packages:

```r
library(WGCNA)
library(biomaRt)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(reshape2)
library(dplyr)
library(readxl)
library(org.Hs.eg.db)
library(clusterProfiler)
library(igraph)
library(ggraph)
library(RColorBrewer)
library(stringr)
library(Seurat)
library(patchwork)
library(Matrix)
library(cowplot)
```
 
### 🐍 Python Packages
Used in the Jupyter notebook integrated_velocity.ipynb:
```
python
scanpy
scvelo
matplotlib
numpy
pandas
anndata
seaborn
```
You can install these with:
```
pip install scanpy scvelo matplotlib numpy pandas anndata seaborn
```



