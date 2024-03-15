library(tidyverse)
library(ggrepel)
library(pheatmap)
library(ComplexHeatmap)
library(DESeq2)
library(cowplot)
library(gridtext)
library(ggtext)
library(ggpubr)
library(rstatix)

# Read DESeq2 results ------------
load("results/bulkRNAseq_results_ds2_1month.RData")

## Load TPM Matrix & Metadata -----------
tpm <- read.csv("data/expr_mat/Danyang_TPM_matrix.csv", check.names = FALSE)
metadata <- read.csv("data/expr_mat/Danyang_Metadata.csv") 

## Genes of Interest ----------
### MGnD & Homeostasis signature (Top 100 DEGs, Clec7a+ vs Clec7a-)
Clec7a_sig <- read.csv("results/bulkRNAseq_results_Clec7a_Krasemann_2017.csv") %>%
  arrange(padj) %>% group_by(direction) %>% 
  dplyr::slice(1:100) %>%
  pull("tracking_id")

### Tgfbr2 related genes
gene_list_Tgfbr2 <- c("Havcr2", "Axl", "Ly9", "Clec7a", "Cd63", "Cxcl16", 
                      "H2-K1", "H2-D1", "Siglech", "Csf1r")
### Phagocytosis related genes
gene_list_phago <- c("Havcr2", "Cd9", "Cst7", "Cstb", "Cxcl16", "Ly9", "Lyz2", 
                     subset(tpm$gene_symbol, str_detect(tpm$gene_symbol, "H2-")))
### 5XFAD related genes
gene_list_5XFAD <- c("Havcr2", "Axl", "Ccl6", "Cxcl16", "Cd9", "Cd81", "Lyz2",
                     subset(tpm$gene_symbol, str_detect(tpm$gene_symbol, "Cts")))
## KEGG phagosome pathway (mmu04145)
KEGG_phagosome <- read.table("data/signatures/KEGG_phagosome.txt", sep = "\t") %>%
  mutate(gene = str_extract(V2, ".*(?=;)")) %>%
  mutate(details = str_extract(V2, "(?<=; ).*")) %>% select(-V2) %>% pull(gene)



## Genes to highlight in the heatmap
htmap_highlight_genes <- c(
  gene_list_Tgfbr2, gene_list_phago, gene_list_5XFAD, 
  "Cst7", "B2m", "Csf2ra", "Fcrls", Clec7a_sig, KEGG_phagosome)

# Fig. 1h: Heatmap visualization of DEGs (FDR < 0.05) of one month old Havcr2cKO mice ----------

column_title <- c("control", "<i>Havcr2</i><sup>cKO</sup>")

DEG_1mo <- res_ordered %>% filter(padj < 0.05) %>%
  dplyr::select(gene_symbol, direction)

heatmap_meta <- metadata %>%
  mutate(genotype = factor(genotype, c("Tim3_flox", "Tim3_cKO")))

htmap_df <- tpm %>% 
  dplyr::select(gene_symbol, heatmap_meta$Sample_ID) %>%
  inner_join(DEG_1mo, by = "gene_symbol") %>%
  mutate(fontface = ifelse(gene_symbol %in% gene_list, 4, 3)) %>%
  mutate(direction = factor(direction, c("up", "down"))) %>%
  mutate(col = ifelse(direction == "up", "red", "royalblue3"))

htmap_df_scaled <- htmap_df %>%
  select(heatmap_meta$Sample_ID) %>%
  t() %>% scale() %>% t()

htmap_col <- circlize::colorRamp2(
  breaks = c(-2, 0, 2), colors = c("navy", "white", "firebrick2"))

label_df <- htmap_df %>%
  mutate(index = 1:n()) %>%
  filter(gene_symbol %in% gene_list)  

annot_row <- rowAnnotation(
  link = anno_mark(
    at = label_df$index,
    labels = label_df$gene_symbol,
    labels_gp = gpar(fontsize = 10, col = label_df$col, fontface = label_df$fontface)),
  annotation_legend_param = list(title_gp = gpar(fontsize = 10.5, fontface = "bold")))

annot_col <- columnAnnotation(
  genotype = heatmap_meta$genotype,
  sex = heatmap_meta$sex,
  col = list(genotype = setNames(c("#3b58a7", "#ed2e23"), c("Tim3_flox", "Tim3_cKO")),
             sex = setNames(ggpubfigs::friendly_pals$muted_nine[c(8, 4)], c("F", "M"))),
  annotation_name_gp = gpar(fontsize = 0, fontface = "bold"),
  annotation_name_side = "left",
  simple_anno_size = unit(3, "mm"),
  annotation_legend_param = list(
    genotype = list(labels = gt_render(column_title))
  ))

p <- htmap_df_scaled %>%
  `colnames<-`(str_remove(colnames(htmap_df_scaled), "_Tim3.*")) %>%
  Heatmap(
    col = htmap_col,
    row_title = NULL,
    row_split = htmap_df$direction,
    show_row_names = FALSE,
    show_row_dend = FALSE,
    column_title = gt_render(column_title),
    column_names_gp = gpar(fontsize = 8),
    column_names_rot = 45,
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    column_split = heatmap_meta$genotype,
    show_column_names = TRUE,
    show_column_dend = FALSE,
    top_annotation = annot_col,
    right_annotation = annot_row,
    name = "Z-score",
    cluster_columns = FALSE,
    cluster_row_slices = FALSE,
    cluster_column_slices = FALSE,
    row_gap = unit(1.5, "mm"),
    column_gap = unit(1.5, "mm")) 

filename <- "figures/Fig.1h_1l_1month/Fig.1h_htmap_1M_DEGs_padj.05_KEGG_phago"

png(paste0(filename, ".png"), width = 6, height = 7, res = 400, units = "in")
draw(p, heatmap_legend_side = "right", annotation_legend_side = "right",
     legend_grouping = "original", merge_legends = TRUE)
dev.off()

pdf(paste0(filename, ".pdf"), width = 6, height = 7)
draw(p, heatmap_legend_side = "right", annotation_legend_side = "right",
     legend_grouping = "original", merge_legends = TRUE)
dev.off()

