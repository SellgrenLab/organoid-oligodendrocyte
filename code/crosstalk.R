options(future.globals.maxSize = 8000 * 1024^2)
library(CellChat)
# Load the CellChatDB.human.rda file
load("CellChatDB.human.rda")

cellchat1 <- readRDS("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/Oligo/scRNAseq/NN_CodeOcean/Oligo_cellchat.rds")

class(cellchat@meta$labels)

head(cellchat@meta)

combine <- c("Proliferating OPCs", "OPCs", "Pre-OPC-2", "Pre-OPC-1")

cellchat@meta$labels <- as.character(cellchat@meta$labels)
cellchat@meta$labels <- ifelse(cellchat@meta$label %in% combine, "OPCs", cellchat@meta$label)
cellchat@meta$labels <- ifelse(cellchat@meta$labels == "Cortical interneurons", "Neurons-3", cellchat@meta$labels)
table(cellchat@meta$labels)
cellchat@meta$labels <- as.factor(cellchat@meta$labels)

CellChatDB <- CellChatDB.human # use Ce
cellchat@DB <- CellChatDB


#### merge everything

cellchat2 <- cellchat1


combine_o <- c("Proliferating OPCs", "OPCs", "Pre-OPC-2", "Pre-OPC-1", "Myelinating OLs")
combine_n <- c("Neurons-1", "Neurons-2", "Cortical interneurons")
combine_m <- c("Microglia")



cellchat2<- subsetCellChat(
  cellchat2,
  idents.use = c(combine_o, combine_n, combine_m))

#cellchat2@meta$labels <- as.character(cellchat2@meta$labels)
#cellchat2@meta$labels <- ifelse(cellchat2@meta$label %in% combine_o, "OL lineage", cellchat2@meta$label)
#cellchat2@meta$labels <- ifelse(cellchat2@meta$labels %in% combine_n,"Neurons", cellchat2@meta$labels)
#table(cellchat2@meta$labels)
#ellchat2@meta$labels <- as.factor(cellchat2@meta$labels)

#cellchat2@idents <- cellchat2@meta$labels

cellchat2 <- subsetData(cellchat2)
cellchat2 <- identifyOverExpressedGenes(cellchat2)
cellchat2 <- identifyOverExpressedInteractions(cellchat2)

cellchat2 <-  computeCommunProb(cellchat2, type = "triMean") # takes forever and is not parallelized

cellchat2 <- filterCommunication(cellchat2, min.cells = 3)

cellchat2 <- computeCommunProbPathway(cellchat2)

cellchat2 <- aggregateNet(cellchat2)
#execution.time = Sys.time() - ptm
#print(as.numeric(execution.time, units = "secs"))
#groupSize <- as.numeric(table(cellchat2@idents))

pathways.show <- c("GAS", "NRXN", "NCAM")

par(mfrow=c(1,1))
netVisual_aggregate(cellchat2, signaling = pathways.show[1], layout = "chord")

gas_lr <- extractEnrichedLR(
  cellchat1,
  signaling= pathways.show,
  geneLR.return = T,
  enriched.only = F,
  thresh = 0.5)
gas_lr[[1]]


df.net <- subsetCommunication(cellchat2)
View(df.net)
computeAveExpr(cellchat1, features = c("GAS6","AXL"), type =  "truncatedMean", trim = 0.1)



levels(cellchat2@idents)
group.cellType <- levels(cellchat2@idents)
group.cellType[c(1,2,6)] <- "Neurons"
group.cellType[c(3,4,5,8)] <- "OPCs"


names(group.cellType) <- levels(cellchat2@idents)


# Save the output of the plot
pdf("/Users/susmita.malwade/cellchat2_chord_plot.pdf")
netVisual_chord_cell(cellchat2, signaling = pathways.show[1], group = group.cellType, title.name = paste0(pathways.show[1], " signaling network"))
dev.off()



netAnalysis_contribution(cellchat2, signaling = pathways.show[1])
pairLR.CXCL <- extractEnrichedLR(cellchat2, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair


pdf("/Users/susmita.malwade/cellchat2_chord_plot_AXL.pdf")
netVisual_individual(cellchat2, signaling = pathways.show[1],  pairLR.use = LR.show)
dev.off()

LR.show <- pairLR.CXCL[2,] 
pdf("/Users/susmita.malwade/cellchat2_chord_plot_MERTK.pdf")
netVisual_individual(cellchat2, signaling = pathways.show[1],  pairLR.use = LR.show)
dev.off()



netVisual_bubble(cellchat2,signaling = pathways.show[1] ,remove.isolate = FALSE)


saveRDS(cellchat2,"/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/Oligo/scRNAseq/processed_data/Oligo_cellchat2.rds")


pdf("/Users/susmita.malwade/cellchat2_chord_plot_ALL.pdf")
netVisual_chord_gene(cellchat2, sources.use = c(1,2,6), targets.use = c(3,4,5,8,9), slot.name = "netP", legend.pos.x = 15, legend.pos.y = 30,   lab.cex = 0.5)
dev.off()
