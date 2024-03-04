library(Seurat)
library(EWCE)
library(dplyr)

setwd("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/Oligo/scRNAseq/processed_data/")

oligo <- readRDS("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/Oligo/scRNAseq/processed_data/cleaned/Oligo_batch_integrated_cleaned.rds")

cellchat <- read.csv("./cellchat_clusters.csv")
cellchat <- dplyr::filter(cellchat, cellchat$cellID %in% annot$cell_id)

annot <- oligo@meta.data
annot$cell_id <- rownames(annot)

annot$celltype <- cellchat$oligo.cellchat_clusters[match(cellchat$cellID, annot$cell_id)]
annot$level1class <- annot$celltype
annot$level2class <- annot$celltype
annot <- dplyr::select(annot, cell_id:level2class)
annot <- dplyr::select(annot, -celltype)

oligo_mrna <- list("exp"= oligo[["RNA_CB"]]@data, "annot"= annot)
oligo_mrna$exp <- fix_bad_hgnc_symbols(oligo_mrna$exp, dropNonHGNC = T)
oligo_mrna$exp <- EWCE::drop_uninformative_genes(
exp = oligo_mrna$exp,
input_species = "human",
output_species = "human",
level2annot = oligo_mrna$annot$level2class)

annotLevels <- list(level1class=cortex_mrna$annot$level1class,
                    level2class=cortex_mrna$annot$level2class)

fNames_oligo <- EWCE::generate_celltype_data(
  exp = oligo_mrna$exp,
  annotLevels = annotLevels,
  groupName = "Oligo") 
load(fNames_oligo)[[1]]


leukod <- c("GLA", "FUCA1","GLB1","GALC", "ARSA", "PSAP", "HSD17B4", "ABCD1", "PEX1", "PEX2", "PEX3", "PEX5", "PEX6", "PEX10", "PEX11", "PEX12", "PEX13", "PEX14", "PEX16", "PEX19", "PEX26", "CYP27A1", "AIFM1", "D2HGDH",
          "IHD2", "L2HGDH", "SLC25A1", "ASPA", "GFAP", "LMNB1", "TUBB4A", "POLR3A", "POLR3B", "POLR3C", "POLR3K", "NKX6-2", "SOX10", "DARS1", "EPRS1", "RARS1", "EIF2B1", "EIF2B2", "EIF2B3", "EIF2B4", "EIF2B5", "MMLC1",
          "HEPACAM", "GJA1", "GJC2", "CNP", "CNTNAP1", "MAG", "MAL", "PLP1", "FAM126A", "SLC35B2", "TMEM63A", "TMEM106B", "TMEM163", "TREX1", "RNASEH2A", "SAMHD1", "ADAR1", "IFIH1")

ndd_risk <- read.xlsx('./ndd_risk_genes.xlsx')

hits= leukod

Leuko_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                 sctSpecies = "human",
                                                 genelistSpecies = "human",
                                                 hits = hits, 
                                                 reps = 10000,
                                                 annotLevel = 1)
ewce_plot(Leuko_results$results, mtc_method="BH")

sessionInfo()


