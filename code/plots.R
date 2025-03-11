library(Seurat)

list.files()
oligo <- readRDS("OL_lineage_subclustering_object.rds")
DefaultAssay(oligo) <- "SCT"

library(SeuratExtend)
DimPlot2(oligo, features = c("PDGFRA", "MOG"))

DimPlot2(
  oligo, label=T,
  features = c("PDGFRA"),
  theme = NoAxes()
) + theme_umap_arrows()
ggsave("/Users/susmita.malwade/umap_arrows.pdf", width = 3, height = 2)


pre_OPC_genes <- list(c("OLIG1", "OLIG2", "SOX8", "ASCL1", "EGFR", "DLL3", "HES6", "GAS5", "DLL1", "ZEB1", "DLX1", "DLX2", "DLX5"))
OPC_genes <- list(c("SOX10", "OLIG1", "OLIG2", "NKX2-2", "TRAF4", "NCALD", "LUZP2", "ETV1", "APOD", "PDGFRA"))
oligo <- AddModuleScore(oligo, features = pre_OPC_genes, name = "pre_OPC_genes", nbin = 10)
oligo <- AddModuleScore(oligo, features = OPC_genes, name = "OPC_genes", nbin = 10)


DimPlot2(
  oligo, label=T, cols="A",
  features = c("OLIG1" ),
  theme = NoAxes()
) + theme_umap_arrows()

ggsave("/Users/susmita.malwade/umap_arrows.pdf", width = 4, height = 4)


FeaturePlot3.grid(oligo, features = c("PDGFRA", "AXL"), pt.size = 1.5)
ggsave("/Users/susmita.malwade/AXL_featureplot_grid.pdf", width = 4, height = 4)
