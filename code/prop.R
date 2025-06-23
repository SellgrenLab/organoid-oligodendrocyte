# Load libraries
library(Seurat)
library(openxlsx)

# Read in the supplementary for subtype classification from : Velmeshev et al. (2023) DOI: 10.1126/science.adf0834
science_data <- read.xlsx("science.adf0834_data_s2.xlsx")
head(science_data)

# Filter for burst pattern and Ex lineage genes
ex_lineage_firstcol <- science_data[grepl("^Ex", science_data$lineage) & science_data$pattern == "burst", "name"]
print(ex_lineage_firstcol)
ex_lineage_firstcol_clean <- sub("_.*", "", ex_lineage_firstcol)
print(ex_lineage_firstcol_clean)


# Score ol for the Ex lineage genes 
ol <- AddModuleScore(ol, features = list(ex_lineage_firstcol_clean), name = "ExLineage", assay = "SCT")
FeaturePlot(ol, features = "ExLineage1")


# Filter for burst pattern and IN lineage genes

in_lineage_firstcol <- science_data[grepl("^IN", science_data$lineage) & science_data$pattern == "burst", "name"]
print(in_lineage_firstcol)
in_lineage_firstcol_clean <- sub("_.*", "", in_lineage_firstcol)
print(in_lineage_firstcol_clean)

# Score ol for the In lineage genes 
ol <- AddModuleScore(ol, features = list(in_lineage_firstcol_clean), name = "InLineage", assay = "SCT")
FeaturePlot(ol, features = "InLineage1")

ol <- AddModuleScore(ol, features = list(ex_lineage_firstcol_clean), name = "ExLineage", assay = "SCT")

colnames(ol@meta.data)

delta_threshold <- 0.01
# Classify neurons based on the difference between excitatory and inhibitory scores
exc_score_col <- "ExLineage1"
inh_score_col <- "InLineage1"
ol$neuron_lineage_type <- ifelse(
  ol[[exc_score_col]] - ol[[inh_score_col]] > delta_threshold,
  "Excitatory",
  ifelse(ol[[inh_score_col]] - ol[[exc_score_col]] > delta_threshold,
         "Inhibitory", "Unclassified")
)

# View proportion of each type
proportions2 <- prop.table(table(ol$neuron_lineage_type))
print(proportions2)
