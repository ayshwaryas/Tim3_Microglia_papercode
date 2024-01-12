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

## palette 
pal <- c("dodgerblue3", "#EF3B2C", friendly_pal("nickel_five", 5)[3],
         friendly_pal("ito_seven", 7)[6])

## Top 100 DEGs up/down-regulated in Clec7a+ vs Clec7a- 
Clec7a_sig <- read.csv("results/bulkRNAseq_results_Clec7a_Krasemann_2017.csv") %>%
  arrange(padj) %>% group_by(direction) %>% 
  dplyr::slice(1:100) %>%
  select(tracking_id, direction) %>%
  split(f = .$direction) %>%
  lapply(pull, "tracking_id")

# Extended Data Fig.11f: MGnD signature score ----------------
## MGnD signature: Top 100 DEGs up-regulated in Clec7a+ compared to Clec7a-
nucseq_harmony_MG_2_3_6_18 <- AddModuleScore(
  nucseq_harmony_MG_2_3_6_18, list(Clec7a_sig$up), name = "UP_in_Clec7a_pos")

VlnPlot(nucseq_harmony_MG_2_3_6_18, "UP_in_Clec7a_pos1", cols = pal, 
        group.by = "new_clusters", split.by = "Genotype_labels") +
  labs(title = "MGnD Signature", x = NULL, y = NULL) +
  theme(legend.text = element_markdown())

ggsave("figures/Fig.6_nucseq_figures/Fig.s11f_vln_up_in_clec7a.png", width = 9, height = 2.7)
ggsave("figures/Fig.6_nucseq_figures/Fig.S11f_vln_up_in_clec7a.pdf", width = 9, height = 2.7)

# Extended Data Fig.11g: Homeostasis signature score ----------------
## Homeostasis signature: Top 100 DEGs down-regulated in Clec7a+ compared to Clec7a-
nucseq_harmony_MG_2_3_6_18 <- AddModuleScore(
  nucseq_harmony_MG_2_3_6_18, list(Clec7a_sig$down), name = "DN_in_Clec7a_pos")

VlnPlot(nucseq_harmony_MG_2_3_6_18, "DN_in_Clec7a_pos1", cols = pal, 
        group.by = "new_clusters", split.by = "Genotype_labels") +
  labs(title = "Homeostasis Signature", x = NULL, y = NULL) +
  theme(legend.text = element_markdown())

ggsave("figures/Fig.6_nucseq_figures/Fig.s11g_vln_down_in_clec7a.png", width = 9, height = 2.7)
ggsave("figures/Fig.6_nucseq_figures/Fig.S11g_vln_down_in_clec7a.pdf", width = 9, height = 2.7)
