library(Seurat)
library(tidyverse)
library(ggpubfigs)
library(ggtext)
library(cowplot)
library(RColorBrewer)
library(ggrepel)
library(ComplexHeatmap)

## Seurat object of microglia nuclei
load("R_objects/2022-09-08.nucseq_harmony_MG_recluster_sub_2_3_6_18.RData")

# P1 and P2 signatures (top 100 DEGs ordered by avg_log2FC) ----------------
load("results/DEG_bimod_TGFB_harmony_2_3_6_18.RData")

top100_DEGs <- DEG_bimod_TGFB_harmony_2_3_6_18 %>%
  group_by(direction) %>% arrange(desc(abs(avg_log2FC))) %>%
  dplyr::slice(1:100) %>% ungroup() %>% split(f = .$direction) %>%
  lapply(select, gene, avg_log2FC) %>%
  lapply(pull, "gene")

# Heatmap of the P1 and P2 signature genes -------------------------------------

htmap_df_DAM_Tim3cKO.5XFAD <- FetchData(
  subset(nucseq_harmony_MG_2_3_6_18, seurat_clusters == 2 & Genotype == "Tim3_cKO.5XFAD"),
  c(top100_DEGs$up, top100_DEGs$down, "bimod_genotype")) %>%
  mutate(Population = ifelse(bimod_genotype = "Tim3_cKO.5XFAD_Top", "P1", "P2"))

annot_col <- columnAnnotation(
  Population = htmap_df_DAM_Tim3cKO.5XFAD$Population,
  col = list(Population = setNames(c("#CC79A7", "#009E73"), c("P1", "P2"))),
  annotation_name_gp = gpar(fontsize = 0, fontface = "bold"),
  simple_anno_size = unit(3, "mm")
)

for(x in c("P1", "P2")) {
  htmap_df <- htmap_df_DAM_Tim3cKO.5XFAD %>% 
    select(top100_DEGs[[ifelse(x == "P1", "up", "down")]]) %>% scale() %>% t()
  
  p <- Heatmap(
    htmap_df, name = "Expression",
    col = circlize::colorRamp2(c(-2.5, 0, 2.5), c("magenta", "black", "yellow")),
    top_annotation = annot_col,
    column_split = htmap_df_DAM_Tim3cKO.5XFAD$Population, 
    row_names_gp = gpar(fontsize = 8, fontface = 3), 
    column_title = NULL,
    show_row_dend = FALSE, 
    show_column_dend = FALSE,
    show_column_names = FALSE,
    cluster_columns = FALSE,
    cluster_column_slices = FALSE
    ) 
  
  filename <- paste0("figures/Fig.6_nucseq_figures/Fig.s12_htmap_bimod_", x, ".png")
  png(paste0(filename, ".png"), width = 4, height = 11, res = 400, units = "in")
  draw(p, merge_legend = TRUE)
  dev.off()
  pdf(paste0(filename, ".pdf"), width = 4, height = 11)
  draw(p, merge_legend = TRUE)
  dev.off()
  
}

