library(Seurat)
library(conos)
library(Matrix)

# read .mtx.gz file

oli <- readMM("/Users/susmita.malwade/oli/matrix.mtx.gz")
features <- read.table("/Users/susmita.malwade/oli/features.tsv.gz", header=FALSE, sep="\t", stringsAsFactors=FALSE)
barcodes <- read.table("/Users/susmita.malwade/oli/barcodes.tsv.gz", header=FALSE, sep="\t", stringsAsFactors=FALSE)
meta <- read.table("/Users/susmita.malwade/oli/meta.tsv", header=T, sep="\t", as.is=T, row.names=1)

rownames(oli) <- features$V1
colnames(oli) <- barcodes$V1



olr <- CreateSeuratObject(counts = oli, meta.data = meta)
 rm(oli, features, barcodes, meta)
# read in umap coordinates
umap <- read.table("/Users/susmita.malwade/oli/UMAP.coords.tsv.gz", header=T, sep="\t", as.is=T, row.names=1)
colnames(umap) <- c("UMAP_1", "UMAP_2")
olr[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap), key = "UMAP_", assay = DefaultAssay(olr))

olr_opc <- subset(olr, subset= Cell_Type %in% c("OPC"))
olr_opc[["RNA"]]$data <- olr_opc[["RNA"]]$counts

rm(olr)

olr_opc<- olr_opc %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(npcs=30, verbose=FALSE)

# read in organoid seurat object
olo <- readRDS("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/Oligo/scRNAseq/NN_CodeOcean/OL_lineage_subclustering_object.rds")

# Create a Conos object
panel.preprocessed <- list("Primary_OPC"=olr_opc, "Organoid_OPC"=olo)

rm(olr_opc, olo)

con <- Conos$new(panel.preprocessed, n.cores=3)

qs::qsave(con, "/Users/susmita.malwade/oli/OL_conos_object.qs")

#con <- qs::qread("/Users/susmita.malwade/oli/OL_conos_object.qs")

con$buildGraph(k=30, k.self=5, space='PCA', ncomps=30, n.odgenes=1000, matching.method='mNN', metric='angular', score.component.variance=TRUE, verbose=TRUE)


con$embedGraph(method="UMAP", min.dist=0.01, spread=25, min.prob.lower=1e-3)
#con$findCommunities(method = igraph::walktrap.community, steps=7)
con$plotPanel(clustering='walktrap', size=0.5, use.common.embedding=TRUE)


annot <- panel.preprocessed[[1]]@meta.data$Age_Range
names(annot) <- rownames(panel.preprocessed[[1]]@meta.data)

annot2 <- Idents(panel.preprocessed[[2]])

con$plotGraph(groups = annot2, embedding="UMAP", use.common.embedding=FALSE, size=0.1, title="Age Range", font.size=4)

con$plotPanel(groups=annot, size=0.1, use.common.embedding=TRUE)


new.label.info <- con$propagateLabels(labels = annot, verbose=TRUE)


#con$plotPanel(colors=new.label.info$labels, show.legend=TRUE, legend.title="Uncertainty", legend.pos=c(1, 0), use.common.embedding=TRUE, size=0.1)

con$plotPanel(color=new.label.info$uncertainty, show.legend=FALSE,use.common.embedding=TRUE, size=0.5)

con$plotPanel(groups=new.label.info$labels, show.legend=FALSE, use.common.embedding=TRUE, size=0.5)
pdf("Velm_ol_age_conos.pdf", width=8, height=6)
con$plotPanel(groups=annot, show.legend=FALSE, use.common.embedding=TRUE, size=0.5)
dev.off()

head(new.label.info$label.distribution)


table(new.label.info$labels)

matched_cells <- intersect(names(new.label.info$labels), colnames(panel.preprocessed[[2]]))
panel.preprocessed[[2]]@meta.data$new.labels <- NA
panel.preprocessed[[2]]@meta.data[matched_cells, "new.labels"] <- new.label.info$labels[matched_cells]

library(ggplot2)
library(dplyr)
library(tidyr)

df <- data.frame(
    celltype = Idents(panel.preprocessed[[2]]),
    new_label = panel.preprocessed[[2]]@meta.data$new.labels
)

df <- df %>% 
    filter(!is.na(new_label)) %>%
    group_by(new_label, celltype) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(new_label) %>%
    mutate(prop = count / sum(count))

df$new_label <- factor(df$new_label, levels = c("2nd trimester", "3rd trimester", "0-1 years","2-4 years", "4-10 years","10-20 years", "Adult"))

pdf("celltype_proportion_barplot_rev.pdf", width=8, height=6)
ggplot(df, aes(x = new_label, y = prop, fill = celltype)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("#EFDDCFFF", "#11C2B5FF", "#72E1E1FF" ,"#009474FF", "#DCBE9BFF", "#B0986CFF")) +
    labs(x = "New Label", y = "Proportion", fill = "Cell Type") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("celltype_proportion_barplot.pdf", width=8, height=6)
ggplot(df, aes(x = celltype, y = prop, fill = new_label)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values= paletteer_d("PNWColors::Anemone")) +
    labs(x = "Cell Type", y = "Proportion", fill = "New Label") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

head(new.label.info$labels)

library(paletteer)
label_levels <- levels(df$new_label)
label_colors <- paletteer_d("PNWColors::Anemone")[1:length(label_levels)]
names(label_colors) <- label_levels

# Create a named vector of colors for each cell, named by cell ID
cell_label_colors <- setNames(label_colors[as.character(new.label.info$labels)], names(new.label.info$labels))

# Example usage in con$plotPanel:
# con$plotPanel(groups=new.label.info$labels, colors=label_colors, show.legend=TRUE, use.common.embedding=TRUE, size=0.5)


pdf("conos_panel_plot.pdf", width=8, height=6)
con$plotPanel(groups=new.label.info$labels, colors=cell_label_colors, show.legend=TRUE, use.common.embedding=TRUE, size=0.9)
dev.off()
con$plotPanel(groups=new.label.info$labels, colors=cell_label_colors, show.legend=TRUE,use.common.embedding=TRUE, size=0.9)





# Ensure pseudotime is present in panel.preprocessed[[2]]@meta.data
if (!"pseudotime" %in% colnames(panel.preprocessed[[2]]@meta.data)) {
    stop("pseudotime column not found in panel.preprocessed[[2]]@meta.data")
}

# Prepare data for plotting
plot_df <- data.frame(
    pseudotime = panel.preprocessed[[2]]@meta.data$pseudotime,
    new_label = panel.preprocessed[[2]]@meta.data$new.labels,
    celltype = Idents(panel.preprocessed[[2]]),
    cell_id = rownames(panel.preprocessed[[2]]@meta.data)
)

plot_df <- plot_df %>% filter(!is.na(new_label))

# Set celltype colors
library(paletteer)
celltype_levels <- levels(plot_df$celltype)
celltype_colors <- paletteer_d("fishualize::Acanthurus_olivaceus")[1:length(celltype_levels)]
names(celltype_colors) <- celltype_levels

plot_df$new_label <- factor(plot_df$new_label, levels = c("2nd trimester", "3rd trimester", "0-1 years","2-4 years", "4-10 years","10-20 years", "Adult"))

pdf("pseudotime_by_label_scatterplot.pdf", width=8, height=6)
ggplot(plot_df, aes(x = new_label, y = pseudotime, fill = celltype)) +
    geom_col(alpha=0.7, size=1.2, width=0.3) +
    scale_fill_manual(values=celltype_colors) +
    labs(x = "New Label", y = "Pseudotime", fill = "Cell Type") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + coord_flip()
dev.off()


# Use the same color palette as above
library(paletteer)
label_levels <- levels(df$new_label)
label_colors <- paletteer_d("fishualize::Acanthurus_olivaceus")[1:length(label_levels)]
names(label_colors) <- label_levels

pdf("pseudotime_scatterplot.pdf", width=8, height=6)
ggplot(plot_df, aes(x = pseudotime, y = 0, color = new_label)) +
    geom_point(alpha=0.7, size=1.2, position=position_jitter(height=0.2)) +
    scale_color_manual(values=label_colors) +
    labs(x = "Pseudotime", y = "", color = "New Label") +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
dev.off()


