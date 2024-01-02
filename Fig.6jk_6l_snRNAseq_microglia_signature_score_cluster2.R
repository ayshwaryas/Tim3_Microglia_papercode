library(Seurat)
library(tidyverse)
library(ggpubfigs)
library(ggtext)
library(cowplot)

## Seurat object of microglia nuclei
load("R_objects/2022-09-08.nucseq_harmony_MG_recluster_sub_2_3_6_18.RData")


## Signatures 
phagosome <- c("Gas7","Lgals3bp","Tac1","Lsp1","Mgat3","Prkch","Gstp1","Ltc4s","Zwint")
alternate_macrophage <- c("Cd163","Fn1","Mgl2","Retnla","Mrc1","Folr2","Selenop","Igf1","Cd36","Psap")

## Loading functions 
source("snRNAseq/0.snRNAseq_functions.R")

# Fig 6k-l: signature score violin plot of alternate macrophage and phagocytosis signatures -------------------
## Alternate macrophage
nucseq_harmony_MG_2_3_6_18 <- AddModuleScore(
  nucseq_harmony_MG_2_3_6_18, list(alternate_macrophage), name = "Alternate_Macrophage")

VlnPlot_custom(nucseq_harmony_MG_2_3_6_18, cluster = 2, 
               feature_name = "Alternate_Macrophage", fig_num = "6k", add_p_value = TRUE,
               split.var = "bimod_genotype", title = "Alternate Macrophage",
               width = 5, height = 3.5)

## Phagosome
nucseq_harmony_MG_2_3_6_18 <- AddModuleScore(
  nucseq_harmony_MG_2_3_6_18, list(phagosome), name = "Phagosome")

VlnPlot_custom(nucseq_harmony_MG_2_3_6_18, cluster = 2, 
               feature_name = "Phagosome", fig_num = "6l", add_p_value = TRUE,
               split.var = "bimod_genotype", title = "Phagocytosis",
               width = 5, height = 3.5, y_pos = c(0.78, 0.98, 0.9))



