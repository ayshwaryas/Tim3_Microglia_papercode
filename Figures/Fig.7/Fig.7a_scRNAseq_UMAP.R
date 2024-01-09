library(tidyverse)
library(cowplot)
library(Seurat)

# Load Seurat object
load("R_objects/2023-10-03.scRNAseq_cortex_MG_updated_clusters.RData")

# palette
updated_clusters_pal <- c(
  "#E31A1C", "#FF8B13", "#FDBF6F", "#F5EA5A", 
  "#0C356A", "dodgerblue2", "#89CFF3", "darkseagreen1", "#D2DE32",
  "#A6CF98", "#557C55", "#8E0152", "#EA1179", "#F79BD3",
  rev(brewer.pal(5, "Purples")),
  c("grey", "black")
)

# UMAP split by genotype
DimPlot(scRNAseq_cortex_MG, group.by = "updated_clusters", split.by = "Genotype",
        cols = updated_clusters_pal)
ggsave('figures/UMAP_updated_clusters_by_genotype_v2.pdf', width = 10, height = 4)
ggsave('figures/UMAP_updated_clusters_by_genotype_v2.png', width = 10, height = 4)