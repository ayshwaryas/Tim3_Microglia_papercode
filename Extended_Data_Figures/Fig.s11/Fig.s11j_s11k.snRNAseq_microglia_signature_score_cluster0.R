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

# Extended Data Fig.11h: signature score of resting signature in cluster 0 ----------------
## Resting signature: Top 20 genes of resting microglia ordered by fold change (Ellwanger et al. 2021, PMID: 33446504)

## Ellwanger et al. 2021 (PMID: 33446504) Dataset S1B 
pnas_sd01 <- read.csv("data/pnas.2017742118.sd01.csv")

## Select top 20 marker genes for each microglia population
pnas_sig <- sapply(str_remove(colnames(pnas_sd01)[4 + 6 * 0:10], "_Specificity"), function(MG_type) {
  pnas_sd01 %>%
    select(Symbol, contains(MG_type)) %>%
    arrange(desc(!!sym(paste0(MG_type, "_Fold_change")))) %>%
    dplyr::slice(1:20) %>%
    pull(Symbol)
}, simplify = FALSE, USE.NAMES = TRUE) 

nucseq_harmony_MG_2_3_6_18 <- AddModuleScore(
  nucseq_harmony_MG_2_3_6_18, list(pnas_sig$Resting), name = "Resting")

VlnPlot_custom(nucseq_harmony_MG_2_3_6_18, cluster = 0, 
               feature_name = "Resting", fig_num = "s11h",
               split.var = "Genotype", title = "Resting", add_p_value = TRUE,
               width = 5, height = 3.5, x_min = 1.115, x_max = 1.335, y_pos = 1)

# Extended Data Fig.11i: signature score of down in Tgfbr2cKO signature in cluster 0 ----------------
# Down in Tgfbr2cKO signature top 100 genes down-regulated in Tgfbr2cKO vs control with log2 FC > 2 ()
Top100DEG_DN_in_Tgfbr2KO <- read.csv("results/bulkRNAseq_results_TGFBRII_Lund_2018.csv") %>%
  select(geneNames, contains("uG")) %>%
  filter(log2fc.uG > 2) %>%
  arrange(padj.uG) %>%
  dplyr::slice(1:100) %>%
  pull(geneNames)

nucseq_harmony_MG_2_3_6_18 <- AddModuleScore(
  nucseq_harmony_MG_2_3_6_18, list(Top100DEG_DN_in_Tgfbr2KO), name = "DN_in_Tgfbr2KO")

VlnPlot_custom(nucseq_harmony_MG_2_3_6_18, cluster = 0, 
               feature_name = "DN_in_Tgfbr2KO", fig_num = "s11i", add_p_value = TRUE,
               split.var = "Genotype", title = "Down in <i>Tgfbr2</i>^cKO",
               width = 5, height = 3.5, x_min = 1.115, x_max = 1.335, y_pos = 0.62)
