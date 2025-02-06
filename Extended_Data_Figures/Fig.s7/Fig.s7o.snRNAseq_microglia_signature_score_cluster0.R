library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(ggtext)
library(rstatix)
library(ggrepel)
library(ggpubr)
library(ggpubfigs)

## Seurat object of microglia nuclei
load("R_objects/2022-09-08.nucseq_harmony_MG_recluster_sub_2_3_6_18.RData")

## Loading functions 
source("snRNAseq/0.snRNAseq_functions.R")

# Extended Data Fig.7o: signature score of down in Tgfbr2cKO signature in cluster 0 ----------------
# Down in Tgfbr2cKO signature top 100 genes down-regulated in Tgfbr2cKO vs control with log2 FC > 2 ()
Top100DEG_DN_in_Tgfbr2KO <- read.csv("results/bulkRNAseq_results_TGFBRII_Lund_2018.csv") %>%
  rename(geneNames = gene_symbol) %>%
  select(geneNames, contains("uG")) %>%
  filter(log2fc.uG < -2) %>%
  arrange(p.ttest.uG) %>%
  dplyr::slice(1:100) %>%
  pull(geneNames)

nucseq_harmony_MG_2_3_6_18 <- AddModuleScore(
  nucseq_harmony_MG_2_3_6_18, list(Top100DEG_DN_in_Tgfbr2KO), name = "DN_in_Tgfbr2KO")

VlnPlot_custom(nucseq_harmony_MG_2_3_6_18, cluster = 0, levs = c("5XFAD", "Tim3_cKO.5XFAD"),
               feature_name = "DN_in_Tgfbr2KO", fig_num = "s7o_ttest", add_p_value = TRUE,
               split.var = "Genotype", title = "Down in <i>Tgfbr2</i><sup>cKO</sup>",
               width = 5, height = 3.5, x_min = 1.115, x_max = 1.335, y_pos = 0.62)
