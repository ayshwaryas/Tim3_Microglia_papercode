library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(ggtext)
library(rstatix)
library(ggrepel)
library(ggpubr)
library(ggpubfigs)

# Load data  -------------------------------------------------------------------
## Seurat object of microglia nuclei
load("R_objects/2022-09-08.nucseq_harmony_MG_recluster_sub_2_3_6_18.RData")

## palette 
pal <- c("dodgerblue3", "#EF3B2C", friendly_pal("nickel_five", 5)[3],
         friendly_pal("ito_seven", 7)[6])

## Load microglia, BAM, PVM and monocyte signatures from snRNAseq/0.snRNAseq_functions.R
source("snRNAseq/0.snRNAseq_functions.R")

# Extended Data Fig.11-2 a-e: Violin plot of microglia/BAM/PVM/monocyte signatures ------
for(x in sig_list) {
  nucseq_harmony_MG_2_3_6_18 <- AddModuleScore(
    nucseq_harmony_MG_2_3_6_18, list(sig_list[[x]]), name = "UP_in_Clec7a_pos")
  
  VlnPlot(nucseq_harmony_MG_2_3_6_18, paste0(x, 1), cols = pal, 
          group.by = "new_clusters", split.by = "Genotype_labels") +
    labs(title = sig_title[x], x = NULL, y = NULL) +
    theme(legend.text = element_markdown())
  
  filename <- paste0("figures/Fig.6_nucseq_figures/Fig.S11.2_vln_sig_score_", x)  
  ggsave(paste0(filename, ".png"), width = 9, height = 2.7)
  ggsave(paste0(filename, ".pdf"), width = 9, height = 2.7)
}
