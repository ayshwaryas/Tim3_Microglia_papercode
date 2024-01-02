library(Seurat)
library(tidyverse)
library(ggpubfigs)
library(ggtext)
library(cowplot)

## Seurat object of microglia nuclei
load("R_objects/2022-09-08.nucseq_harmony_MG_recluster_sub_2_3_6_18.RData")

## Loading functions 
source("snRNAseq/0.snRNAseq_functions.R")

# DEGs comparing P1 and P2 -----------------------------------------------------
DAM_Tim3cKO.5XFAD <- subset(nucseq_harmony_MG_2_3_6_18, seurat_clusters == 2 & Genotype == "Tim3_cKO.5XFAD")

DEG_bimod_TGFB_harmony_2_3_6_18 <- FindMarkers(
  DAM_Tim3cKO.5XFAD, ident.1 = "Tim3_cKO.5XFAD_Top", ident.2 = "Tim3_cKO.5XFAD_Btm",
  min.pct = 0.1, logfc.threshold = 0, test.use = "wilcox", group.by = "bimod_genotype") %>%
  rownames_to_column("gene") %>%
  mutate(direction = ifelse(avg_log2FC > 0, 'up', 'down'))  

save(DEG_bimod_TGFB_harmony_2_3_6_18, file = "results/DEG_bimod_TGFB_harmony_2_3_6_18.RData")

# Get P1 and P2 signatures (top 100 DEGs ordered by avg_log2FC) ----------------
top100_DEGs <- DEG_bimod_TGFB_harmony_2_3_6_18 %>%
  group_by(direction) %>% arrange(desc(abs(avg_log2FC))) %>%
  dplyr::slice(1:100) %>% ungroup() %>% split(f = .$direction) %>%
  lapply(select, gene, avg_log2FC) %>%
  lapply(pull, "gene")


# Violin plot of P1 and P2 signature scores ------------------------------------
nucseq_harmony_MG_2_3_6_18 <- AddModuleScore(nucseq_harmony_MG_2_3_6_18, list(top100_DEGs$up), name = "P1_signature")
nucseq_harmony_MG_2_3_6_18 <- AddModuleScore(nucseq_harmony_MG_2_3_6_18, list(top100_DEGs$down), name = "P2_signature")

VlnPlot_custom(nucseq_harmony_MG_2_3_6_18, feature_name = "P1_signature", fig_num = "6d",
               split.var = "bimod_genotype", title = "P1 Signature")

VlnPlot_custom(nucseq_harmony_MG_2_3_6_18, feature_name = "P2_signature", fig_num = "6e",
               split.var = "bimod_genotype", title = "P2 Signature")
