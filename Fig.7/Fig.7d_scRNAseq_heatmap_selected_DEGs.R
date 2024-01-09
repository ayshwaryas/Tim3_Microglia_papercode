library(Seurat)
library(tidyverse)
library(harmony)
library(RColorBrewer)
library(ggtext)
library(parallel)
library(ComplexHeatmap)
library(ggpubr)
library(rstatix)

# Load Data -------------------------------------------------------
## Seurat object
load("R_objects/2023-10-03.scRNAseq_cortex_MG_updated_clusters.RData")

load("results/20231005_scRNAseq_cortex_MG_updated_clusters_DGE_wilcox_min.pct1_named.RData")

## Differential gene expression analysis (Wilcoxon) results, Havcr2icKO;5xFAD vs 5xFAD
load("results/20231005_scRNAseq_cortex_MG_updated_clusters_DGE_wilcox_min.pct1_named.RData")

## DEGs (FDR < 0.05)
scRNAseq_cortex_MG_DEG <- subset(scRNAseq_cortex_MG_DGE_wilcox_min.pct1, FDR < 0.05)
scRNAseq_cortex_MG_DEG_rm_ribo <- subset(scRNAseq_cortex_MG_DGE_wilcox_min.pct1_rm_ribo, FDR < 0.05)


## Genes of interest
gene_list <- list(
  MGnD = data.frame(gene = c("Axl", "Apoe", "Lpl", "Cst7", "Lyz2", "Ctsb", "Ctsz", "Cd9")),
  Homeostatic = data.frame(gene = c("Fcrls", "Cst3", "P2ry12", "Tmem119", "Serinc3", "Siglech", "Hexb")),
  Phagosome = data.frame(gene = c("Calr", "Abca1", "Bin2", "Trem2", "Cdc42", "Camk1d", "Ncf1", "Hspa8")),
  Lysosome = data.frame(gene = c("Cd63", "Cd68", "Ctsb", "Ctsd", "Hexb", "Ctsc", "Cd164")),
  TAM = data.frame(gene = c("Mertk", "Gas6")),
  MHC = data.frame(gene = c("H2-D1", "H2-K1", "H2-M3", "H2-Ab1", "H2-Q7")),
  P1_P2 = data.frame(gene = c("Btg2", "Sp100", "Ski", "Nrp1")),
  TNF = data.frame(gene = c("Nfkbia", "Jun", "Junb")),
  PI3K_Akt = data.frame(gene = c("Sgk1", "Hsp90b1", "Ywhah")),
  Autophagy = data.frame(gene = c("Ubb", "Ubc", "Uba52", "Scoc")),
  IFNA = data.frame(gene = c("Eif2ak2", "B2m", "Bst2", "Cd74", "Tent5a", "Lgals3bp", "Ifitm3", "Ifi44")),
  Other = data.frame(gene = c("Lilrb4a", "Egr1", "Hspa5"))) %>%
  data.table::rbindlist(idcol = "group")  %>%
  mutate(group_simplified = ifelse(group %in% c("MGnD", "Homeostatic", "MHC", "IFNA", "Phagosome", "Lysosome"), group, "Other")) %>%
  mutate(col = case_when(
    group == "MGnD" ~ "brown2",
    group == "Homeostatic" ~ "dodgerblue2", 
    group == "Phagosome" ~ brewer.pal(6, "Dark2")[4], 
    group == "Lysosome" ~ brewer.pal(6, "Dark2")[5],
    group == "MHC" ~ "green",
    group == "IFNA" ~ "orange",
    TRUE ~ "grey15"
  )) %>%
  mutate(face = ifelse(group_simplified == "Other", 3, 4)) %>%
  mutate(group_gene = paste(group, gene, sep = "."))



# Fig.7c Heatmap of selected DEGs ------------------------------------------------------
## color: log2FC comparing Havcr2icKO;5xFAD to 5xFAD
## asterisk: significance of Havcr2icKO;5xFAD to 5xFAD

## Legend 
lgd_signif = Legend(
  pch = c("****", "*** ", "**  ", "*   "), title = "FDR", type = "points", 
  labels = c(" <0.0001", " <0.001", " <0.01", " <0.05"), 
  size = unit(0.06, "mm"), background = "white")

## Data frame (wide format) of log2FC of Havcr2icKO;5xFAD vs. 5xFAD of the genes of interest in each cluster
log2FC_df <- scRNAseq_cortex_MG_DEG_rm_ribo %>%
  data.table::rbindlist(idcol = "cluster") %>%
  mutate(cluster = factor(cluster, names(scRNAseq_cortex_MG_DEG_rm_ribo))) %>%
  filter(gene %in% gene_list$gene) %>%
  arrange(cluster) %>%
  select(gene, cluster, avg_log2FC) %>%
  pivot_wider(names_from = "cluster", values_from = "avg_log2FC") %>%
  arrange(gene) %>%
  right_join(gene_list, by = "gene")  %>%
  mutate(group_gene = factor(group_gene, gene_list$group_gene)) %>%
  arrange(group_gene)

## Data frame (wide format) of the significance of Havcr2icKO;5xFAD vs. 5xFAD of the genes of interest in each cluster
FDR_signif_df <- scRNAseq_cortex_MG_DEG_rm_ribo %>%
  data.table::rbindlist(idcol = "cluster") %>%
  mutate(cluster = factor(cluster, names(scRNAseq_cortex_MG_DEG_rm_ribo))) %>%
  filter(gene %in% gene_list$gene) %>%
  arrange(cluster) %>%
  select(gene, cluster, FDR) %>%
  rstatix::add_significance(p.col = "FDR",
                            symbols = c("****", "***", "**", "*", "")) %>%
  select(-FDR) %>%
  pivot_wider(names_from = "cluster", values_from = "FDR.signif", values_fill = "") %>%
  right_join(gene_list[, c("gene", "group_gene")], by = "gene")  %>%
  mutate(group_gene = factor(group_gene, gene_list$group_gene)) %>%
  arrange(group_gene) %>%
  select(-gene, -group_gene)

htmap_log2FC_df <- log2FC_df %>%
  select(colnames(FDR_signif_df))

all(rownames(FDR_signif_df) == rownames(htmap_log2FC_df))
all(colnames(FDR_signif_df) == colnames(htmap_log2FC_df))


## Heatmap
p_log2FC <- Heatmap(
  htmap_log2FC_df,
  na_col = "grey90",
  border_gp = gpar(col = "grey15"),
  cluster_rows = FALSE, 
  cluster_row_slices = FALSE,
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  column_split = factor(str_remove(colnames(htmap_log2FC_df), "_([0-9]|HMG|DAM)$") %>% str_replace("CC_", "CC_\n"),
                        c("HMG", "DAM", "CC_\nHMG", "CC_\nDAM", "IFN")),
  col = circlize::colorRamp2(c(-1, 0, 0.5), c("blue", "white", "brown2")),
  row_split = factor(str_replace(log2FC_df$group_simplified, "Homeostatic", "Homeo-\nstatic") %>% str_replace("some$", "-\nsome"),
                     c("MGnD", "Homeo-\nstatic", "Phago-\nsome", "Lyso-\nsome", "MHC", "IFNA", "Other")),
  column_names_gp = gpar(fontsize = 10), 
  rect_gp = gpar(col = "white"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(FDR_signif_df[i, j], x, y,
              gp = gpar(fontsize = 13), vjust = 0.8) },
  name = "avg_log2FC",
  row_title_gp = gpar(fontsize = 11, fontface = 2, 
                      col = c("brown2", "dodgerblue2", brewer.pal(6, "Dark2")[4:5], "orange", "#18AAAA", "grey15")),
  column_title_gp = gpar(fontsize = 11, fontface = 2),
  column_names_rot = 45,
  heatmap_legend_param = list(grid_height = unit(6, "mm")),
  row_gap = unit(2, "mm"),
  column_gap = unit(2, "mm"),
  row_names_side = "left",
  row_names_gp = gpar(fontface = 3, fontsize = 10),
  row_labels = log2FC_df$gene
)



png("figures/htmap_log2FC_selected_DEGs_rm_ribo_FDR.png", width = 7, height = 9, unit = "in", res = 400)
draw(p_log2FC, annotation_legend_list = list(lgd_signif), merge_legend = TRUE)
dev.off()

cairo_pdf("figures/htmap_log2FC_selected_DEGs_rm_ribo_FDR.pdf", width = 7, height = 9)
draw(p_log2FC, annotation_legend_list = list(lgd_signif), merge_legend = TRUE)
dev.off()

