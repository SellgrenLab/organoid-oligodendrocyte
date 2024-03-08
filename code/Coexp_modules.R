





library(Seurat)
library(WGCNA)
library(hdWGCNA)
library(tidyverse)
library(cowplot)
library(patchwork)
setwd("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/Oligo/scRNAseq/")
# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

# load the snRNA-seq dataset
seurat_obj <- readRDS("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/External datasets/Single Cell Studies/kriegstein human adolescence/Adolescent_Oligodendrocytes_sobj.rds")

p <- DimPlot(seurat_obj, group.by='cell_type', label=TRUE) +
  umap_theme() + ggtitle('Velmeshev et al') + NoLegend()

p

#setup the data for analysis
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "Oligo-organoid" # the name of the hdWGCNA experiment
)

seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj) 
seurat_obj <- RunPCA(seurat_obj)

# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("cell_type", "age"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'pca', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'cell_type' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)

seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "OPCs", # the name of the group of interest in the group.by column
  group.by='cell_type', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)
# Use soft power threshold of 4

#further inspection
#power_table <- GetPowerTable(seurat_obj)
#head(power_table)


# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj, 
  tom_name = 'OPCs' # name of the topoligical overlap matrix written to disk
)

PlotDendrogram(seurat_obj, main='OPCs hdWGCNA Dendrogram')

# need to run ScaleData first or else harmony throws an error:
#seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
  seurat_obj#,
  #group.by.vars="age"
)

# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'cell_type', group_name = 'OPCs'
)

# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "OPCs-M"
)

# plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj, ncol=5)

p

# get the module assignment table:
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')

# show the first 6 columns:
head(modules[,1:6])
# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)

head(hub_df)

saveRDS(seurat_obj, file='./results/hdWGCNA/Velmeshev_OPCs_hdWGCNA_object.rds')

ModuleNetworkPlot(
  seurat_obj,
  outdir = 'ModuleNetworks'
)

