library(tidyverse)
library(data.table)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)


# Load Data ---------------------
## 3 month Bulk RNA-seq, TPM Matrix and metadata
outliers_batch1 <- c("KK1", "KK2", "KK7", "KK9", "KK16")
tpm <- read.csv("data/expr_mat/tim3_Kimi_tpm_with_genesymbols_KK.csv")
metadata <- read.csv("data/expr_mat/tim3_Kimi_metadata_KK.csv") %>%
  filter(Sample_ID %in% colnames(tpm)[-(1:2)]) %>%
  mutate(cell_type.genotype = paste(cell_type, genotype))

## 3 month Bulk RNA-seq, Differential expression analysis results
load("results/bulkRNAseq_results_ds1_batch1_3month.RData")

# Genes of interest to be highlighted in the heatmap ---------------------------
### MGnD & Homeostasis signature (Top 100 DEGs, Clec7a+ vs Clec7a-)
Clec7a_sig <- read.csv("data/Oleg_Immunity_2017_DEG_ADpos_vs_neg_all.csv") %>%
  arrange(padj) %>% group_by(direction) %>% 
  dplyr::slice(1:100) %>%
  select(tracking_id, direction)

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


# Generate heatmap -------
htmap_overlap_3sets <- function(overlap, batch, outliers, gene_list, fontsize_row = 8, 
                                height = 8.5, width = 5.5, suffix = "", phago_sig = NULL, padding = 1, link_width = 5,
                                gap = 1) {
  
  meta_sub <- metadata %>%
    filter(Batch == batch & !(Sample_number %in% outliers)) %>%
    mutate(Sample_num = str_extract(Sample_number, "[0-9]+") %>% as.numeric())%>%
    mutate(cell_type = ifelse(cell_type == "phago_plus", "phago+", "phago-")) %>%
    mutate(cell_type = factor(cell_type, c("phago-", "phago+"))) %>%
    mutate(genotype = ifelse(genotype == "control", "control", "Havcr2cKO")) %>%
    mutate(genotype = factor(genotype, c("control", "Havcr2cKO"))) %>%
    arrange(cell_type, genotype)
  
  htmap_df <- tpm %>%
    select(gene_id, gene_name, meta_sub$Sample_ID) %>%
    inner_join(overlap, by = c("gene_id", "gene_name")) %>%
    mutate(fontface = ifelse(gene_name %in% c(gene_list, Clec7a_sig$tracking_id, phago_sig), 4, 3)) %>%
    mutate(col = case_when(!str_detect(dir, "down") ~ "red",
                           !str_detect(dir, "up")  ~ "blue",
                           TRUE ~ "grey15"))
  
  
  htmap_df_scaled <- htmap_df %>% dplyr::select(meta_sub$Sample_ID) %>%
    t() %>% scale() %>% t()
  
  htmap_col <- circlize::colorRamp2(c(-3, 0, 3),
                                    c("navy", "white", "firebrick2"), transparency = 0.01)
  
  annot_col <- columnAnnotation(
    `cell type` = meta_sub$cell_type,
    genotype = meta_sub$genotype,
    col = list(`cell type` = setNames(c("gold", "grey25"), c("phago+", "phago-")),
               genotype = setNames(c("#3b58a7", "red"), c("control", "Havcr2cKO"))),
    annotation_name_gp = gpar(fontsize = 0, fontface = "bold"),
    annotation_name_side = "left",
    simple_anno_size = unit(3, "mm"),
    annotation_legend_param = list(
      `cell type` = list(labels = gt_render(c("phago-", "phago+"))),
      genotype = list(labels = gt_render(c("control", "<i>Havcr2</i><sup>cKO</sup>"))))
  )
  
  label_df <- htmap_df %>%
    mutate(index = 1:n()) %>%
    filter(gene_name %in% c(gene_list, Clec7a_sig$tracking_id, phago_sig)) 
  
  annot_row <- rowAnnotation(
    `phago vs non-phago MG` = htmap_df$direction.phago_sig,
    `Havcr2cKO vs control in phago+` = htmap_df$direction.tim3sig.phagopos,
    `Havcr2cKO vs control in phago-` = htmap_df$direction.tim3sig.phagoneg,
    link = anno_mark(
      at = label_df$index,
      labels = label_df$gene_name,
      link_width = unit(link_width, "mm"), 
      labels_gp = gpar(fontsize = fontsize_row, col = label_df$col,
                       fontface = 4),
      padding = unit(padding, "mm")), 
    col = list(`phago vs non-phago MG` = setNames(c("#F17720", "#0474BA"), c("up", "down")),
               `Havcr2cKO vs control in phago+` = setNames(c("darkorange", "steelblue"), c("up", "down")),
               `Havcr2cKO vs control in phago-` = setNames(c("darkorange", "steelblue"), c("up", "down"))),
    annotation_label = gt_render(c("(1)", "(2)", "(3)", "")),
    annotation_name_rot = 0,
    annotation_name_gp = gpar(fontsize = 9.5, fontface = "bold"),
    annotation_name_side = "top",
    simple_anno_size = unit(3.2, "mm"),
    gap = unit(2, "points"),
    annotation_legend_param = list(
      `phago vs non-phago MG` = list(title = gt_render("(1) phago+ vs phago-<br/>in control"),
                                     title_gp = gpar(fontsize = 10, fontface = "bold", lineheight = 1.15)),
      `Havcr2cKO vs control in phago+` = list(title = gt_render("(2) <i>Havcr2</i><sup>cKO</sup> phago+<br/>vs control phago+"),
                                              title_gp = gpar(fontsize = 10, fontface = "bold", lineheight = 1.15)),
      `Havcr2cKO vs control in phago-` = list(title = gt_render("(3) <i>Havcr2</i><sup>cKO</sup> phago-<br/>vs control phago-"),
                                              title_gp = gpar(fontsize = 10, fontface = "bold", lineheight = 1.15))
    )
  )
  
  p <- Heatmap(
    htmap_df_scaled,
    row_labels = htmap_df$gene_name,
    row_names_gp = gpar(fontsize = fontsize_row, col = htmap_df$col, 
                        fontface = htmap_df$fontface),
    column_names_gp = gpar(fontsize = 7),
    column_labels = meta_sub$Sample_number,
    col = htmap_col,
    column_names_rot = 45,
    column_split = meta_sub$cell_type.genotype,
    top_annotation = annot_col,
    right_annotation = annot_row,
    cluster_rows = TRUE,
    cluster_row_slices = FALSE,
    cluster_columns = TRUE,
    cluster_column_slices = FALSE,
    column_dend_height = unit(4, "mm"),
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    use_raster = TRUE,
    row_title = NULL,
    row_split = htmap_df$dir,
    row_gap = unit(gap, "mm"),
    column_title = NULL,
    name = "Z-score"
  ) 
  
  filename <- paste0("figures/Fig.3g_3month_htmap_overlap_phago_tim3/",
                     "Fig.3g_3month_htmap_overlap_", batch, "_", suffix)
  png(paste0(filename, ".png"), width = width, height = height, 
      res = 400, units = "in")
  draw(p, merge_legends = TRUE)
  dev.off()
  pdf(paste0(filename, ".pdf"), width = width, height = height)
  draw(p, merge_legends = TRUE)
  dev.off()
}

## Top DEGs (order genes by FDR, select the first max(300, n_DEGs) up and down-regulated genes,
## where n_DEGs is the number of DEGs (FDR < 0.1) of a specific direction

### (1) Havcr2cKO vs control in phagocytosing microglia
tim3sig_phagopos_top300 <- results_batch1$Tim3cKOvscontrol_phagopos %>%
  group_by(direction) %>% 
  dplyr::slice(1:max(300, sum(padj < 0.1, na.rm = TRUE))) %>% ungroup %>%
  select(gene_id, gene_name, direction)  
### (2) Havcr2cKO vs control in non-phagocytosing microglia
tim3sig_phagoneg_top300 <- results_batch1$Tim3cKOvscontrol_phagoneg %>%
  group_by(direction) %>% 
  dplyr::slice(1:max(300, sum(padj < 0.1, na.rm = TRUE))) %>% ungroup %>%
  select(gene_id, gene_name, direction)  
### (3) phagocytosing control vs non-phagocytosing control microglia
phagosig_top300 <- results_batch1$phagoposvsneg_control %>% 
  group_by(direction) %>% 
  dplyr::slice(1:max(300, sum(padj < 0.1, na.rm = TRUE))) %>% ungroup %>%
  select(gene_id, gene_name, direction) %>%
  dplyr::rename("direction.phago_sig" = "direction")

## Joining the three DEG dataframes
overlap.phagoSig_ctrl.tim3Sig_phagoPosNeg <- full_join(
  tim3sig_phagopos_top300, tim3sig_phagoneg_top300, by = c('gene_id', 'gene_name'),
  suffix = c(".tim3sig.phagopos", ".tim3sig.phagoneg")) %>%
  full_join(phagosig_top300, by = c('gene_id', 'gene_name')) %>% 
  rowwise() %>%
  mutate(n_na = is.na(direction.tim3sig.phagopos) + is.na(direction.tim3sig.phagoneg) + is.na(direction.phago_sig)) %>%
  filter(n_na < 2) %>% ungroup() %>%
  mutate_at(vars(starts_with("direction")), factor, c("up", "down")) %>%
  mutate(dir = paste(direction.phago_sig, direction.tim3sig.phagopos, direction.tim3sig.phagoneg, sep = ".")) %>%
  mutate_at(vars(starts_with("direction")), 
            list(num = function(x) {case_when(x == "up" ~ 1, x == "down" ~ -1, TRUE ~ 0)})) %>%
  mutate(dir_num = direction.tim3sig.phagopos_num + direction.tim3sig.phagoneg_num + direction.phago_sig_num) %>%
  mutate(dir_num_abs = abs(direction.tim3sig.phagopos_num) + abs(direction.tim3sig.phagoneg_num) + abs(direction.phago_sig_num)) %>%
  arrange(desc(dir_num), desc(dir_num_abs), desc(dir)) %>%
  mutate(dir = factor(dir, unique(.$dir))) 


htmap_overlap_3sets(overlap = overlap.phagoSig_ctrl.tim3Sig_phagoPosNeg,
                    batch = "Batch_1", outliers = outliers_batch1, 
                    gene_list = c(gene_list_phago, "Cd33"), width = 6, height = 9.3, 
                    phago_sig = KEGG_phagosome, 
                    suffix = "tim3.sig_vs_phago.sig_KEGG_phago_posneg", 
                    padding = 0.5, gap = 0.5) 
