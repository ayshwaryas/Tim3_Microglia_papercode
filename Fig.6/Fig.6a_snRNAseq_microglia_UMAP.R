library(Seurat)
library(tidyverse)
library(ggpubfigs)
library(ggtext)

## Seurat object of microglia nuclei
load("R_objects/2022-09-08.nucseq_harmony_MG_recluster_sub_2_3_6_18.RData")

## palette
pal_UMAP <- c(friendly_pal("muted_nine", 3)[c(1, 3, 2)], friendly_pal("nickel_five", 5),
              friendly_pal("bright_seven", 7)[5])

## UMAP
DimPlot(nucseq_harmony_MG_2_3_6_18, group.by = "new_clusters", 
        cols = pal_UMAP, ncol = 2, split.by = "Genotype_labels",
        label = TRUE, repel = TRUE, label.size = 3.5)  +
  labs(title = NULL) + 
  theme(strip.text = element_markdown(),
        panel.spacing.y = unit(0, "mm"))

ggsave("figures/Fig.6_nucseq_figures/Fig.6a_dimplot_MG.png", width = 8, height = 6)
ggsave("figures/Fig.6_nucseq_figures/Fig.6a_dimplot_MG.pdf", width = 8, height = 6)
