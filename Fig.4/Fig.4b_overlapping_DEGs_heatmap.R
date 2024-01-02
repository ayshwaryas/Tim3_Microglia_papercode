library(corrplot)
library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(cowplot)

# Read data -----------------------------------------
## correct gene symbol
correct_gene_symbol <- function(df, var = "gene_symbol") {
  df %>%
    mutate(gene_symbol = case_when(
      gene_symbol == "02-mars" ~ "Mars2",
      str_detect(gene_symbol, "-Mar") ~ paste0("March", gene_symbol),
      str_detect(gene_symbol, "-Sep") ~ paste0("Sept", gene_symbol),
      TRUE ~ gene_symbol),
      gene_symbol = str_remove(gene_symbol, "-Mar|-Sep"))
}

res_order_slice <- function(res, flip = TRUE, slice = TRUE, thres = 0.1) {
  res <- res %>%
    arrange(padj) %>%
    correct_gene_symbol()
  
  if(flip) {
    res <- res %>% 
      mutate(log2FoldChange = -log2FoldChange, 
             direction = ifelse(direction == "up", "down", "up"))
  }
  
  if(slice) {
    res <- res %>% 
      group_by(direction) %>%
      dplyr::slice(1:max(300, sum(padj < thres, na.rm = TRUE))) %>%
      ungroup
  }
  
  res <- res %>% dplyr::select(gene_symbol, log2FoldChange, direction)
  
  return(res)
}


## (1) 1-month-old mice, Havcr2cKO vs Havcr2flox/flox
load("data/DGE_results/2022-01-04.Danyang_DGE_results.RData")
tim3_all <- res_ordered$`1M` %>% correct_gene_symbol() %>%
  select(gene_symbol, log2FoldChange, direction, padj)
tim3_sig <- tim3_all %>% res_order_slice

## (2) 3-month-old mice, phagocytosing control vs non-phagocytosing control microglia
load("data/DGE_results/2022-05-03.Dataset1_DGE_res_ordered.RData")
phago_all <- results_batch1_ordered$`phago+ vs phago- in control` %>%
  dplyr::rename("gene_symbol" = "gene_name") %>%
  correct_gene_symbol() %>%
  select(gene_symbol, log2FoldChange, direction, padj)
phago_sig <- phago_all %>% res_order_slice(flip = FALSE)

## (3) Tgfbr2cKo vs control (Lund et al. 2018, PMID: 29662171)
TGFBRII_all <- read.csv("data/FPKM OB edits - with padj.csv") %>%
  select(geneNames, contains(".uG")) %>%
  dplyr::rename("gene_symbol" = "geneNames", "log2FoldChange" = "log2fc.uG",
                "direction" = "dir.uG", "padj" = "padj.uG") %>%
  filter(rowSums(.[, 2:7]) != 0) %>%
  correct_gene_symbol()  %>%
  select(gene_symbol, log2FoldChange, direction, padj)

TGFBRII_sig <- TGFBRII_all %>% res_order_slice(flip = TRUE, thres = 0.1)

## (4) Clec7a+ vs Clec7a- (Krasemann et al., 2017, PMID: 28930663)
Clec7a_all <- read.csv("data/Oleg_Immunity_2017_DEG_ADpos_vs_neg_all.csv") %>%
  dplyr::rename("gene_symbol" = "tracking_id", "log2FoldChange" = "log2FC") %>%
  correct_gene_symbol() %>%
  select(gene_symbol, log2FoldChange, direction, padj)
Clec7a_sig <- Clec7a_all %>% res_order_slice(flip = FALSE)


## Combine Datasets 
fig5_ls <- list(`Tim3KO` = tim3_sig, `phago+` = phago_sig, 
                `Tgfbr2KO` = TGFBRII_sig, `Clec7a+`= Clec7a_sig)
fig5_ls_all <- list(`Tim3KO` = tim3_all, `phago+` = phago_all, 
                    `Tgfbr2KO` = TGFBRII_all, `Clec7a+`= Clec7a_all)

fig5_df <- fig5_ls %>%
  data.table::rbindlist(idcol = "dataset") %>%
  mutate(direction = factor(direction, c("up", "down"))) %>%
  arrange(dataset, direction)

fig5_df_all <- fig5_ls_all %>%
  data.table::rbindlist(idcol = "dataset") %>%
  mutate(direction = factor(direction, c("up", "down"))) %>%
  arrange(dataset, direction)

gene_list_Tgfbr2 <- c("Havcr2", "Axl", "Ly9", "Clec7a", "Cd63", "Cxcl16", 
                      "H2-K1", "H2-D1", "Siglech", "Csf1r")
gene_list_phago <- c("Havcr2", "Cd9", "Cst7", "Cstb", "Cxcl16", "Ly9", "Lyz2", 
                     unique(subset(fig5_df, str_detect(gene_symbol, "H2\\-"))$gene_symbol))

replace_name <- function(name) {
  data.frame(name = name) %>%
    mutate(name = str_remove(name, "log2FoldChange_")) %>%
    mutate(name = ifelse(str_detect(name, "\\+"),
                         paste0(name, "\nvs. ", str_replace(name, "\\+", "-")),
                         paste0(name, "\n vs. control"))) %>%
    mutate(name = factor(name, .$name)) %>%
    pull(name)
}

# Fig. 4b: Correlation between gene vectors ----------
scale_df <- function(df, exclude_col = 1:2) {
  df_scaled <- t(df[, -exclude_col]) %>% scale() %>% t()
  cbind(df[, exclude_col], df_scaled) %>% as_tibble()
}

annot_param <- list(
  `Tim3KO` = list(title = gt_render("(1) ***Havcr2***<sup>cKO</sup> vs control"),
                  title_gp = gpar(fontsize = 10, fontface = 2, lineheight = 1.15)),
  `phago+` = list(title = gt_render("(2) Phagocytosing vs<br/>Non-phagocytosing MG"),
                  title_gp = gpar(fontsize = 10, fontface = 2, lineheight = 1.15)),
  `Tgfbr2KO` = list(title = gt_render("(3) ***Tgfbr2***<sup>cKO</sup> vs control"),
                    title_gp = gpar(fontsize = 10, fontface = 2, lineheight = 1.15)),
  `Clec7a+` = list(title = gt_render("(4) ***Clec7a***<sup>+</sup> vs ***Clec7a***<sup>-</sup>"),
                   title_gp = gpar(fontsize = 10, fontface = 2, lineheight = 1.15))
)

htmap_corr_alltim3 <- function(fig5F_df = NULL, n_min, suffix = "", gene_list = NULL) {
  
  ## (1) TPM, 1-month-old mice, Havcr2cKO vs Havcr2flox/flox
  tim3_TPM <- read.csv("data/expr_mat/Danyang_TPM_matrix.csv", check.names = FALSE) %>%
    `colnames<-`(c("gene_symbol", colnames(.)[-1])) %>%
    filter(gene_symbol %in% fig5F_df$gene_symbol) %>%
    select(gene_symbol, contains("1M")) %>%
    scale_df(exclude_col = 1) %>%
    dplyr::rename("gene_symbol" = "V1")
  
  
  ## (2) 3-month-old mice, phagocytosing control vs non-phagocytosing control microglia
  ## Metadata
  phago_meta <- read.csv("data/expr_mat/tim3_Kimi_metadata_KK.csv") %>%
    filter(Batch == "Batch_1") %>%
    filter(genotype == "control") %>%
    filter(Sequencing_run != "excluded") %>%
    filter(!Sample_number %in% c("KK1", "KK2", "KK7", "KK9", "KK16")) %>%
    arrange(cell_type)
  
  ## TPM
  phago_TPM <- read.csv("data/expr_mat/tim3_Kimi_tpm_with_genesymbols_KK.csv") %>%
    filter(gene_id %in% results_batch1_ordered$`phago+ vs phago- in control`$gene_id) %>%
    select(gene_symbol = gene_name, phago_meta$Sample_ID) %>%
    filter(gene_symbol %in% fig5F_df$gene_symbol) %>%
    select(gene_symbol, everything())  %>%
    scale_df(exclude_col = 1) %>%
    dplyr::rename("gene_symbol" = "V1")
  
  
  ## (3) FPKM, Tgfbr2cKo vs control (Lund et al. 2018, PMID: 29662171)
  TGFBRII_FPKM <- read.csv("data/FPKM OB edits - with padj.csv") %>%
    dplyr::rename("gene_symbol" = "geneNames") %>% 
    filter(gene_symbol %in% fig5F_df$gene_symbol) %>%
    select(gene_symbol, contains("WT.uG"), contains("KO.uG")) %>%
    scale_df(exclude_col = 1) %>%
    dplyr::rename("gene_symbol" = "V1")
  
  ## (4) FPKM, Clec7a+ vs Clec7a- (Krasemann et al., 2017, PMID: 28930663)
  Clec7a_FPKM <- read.csv("data/Oleg_Immunity_2017_DEG_ADpos_vs_neg_all.csv") %>%
    dplyr::rename("gene_symbol" = "tracking_id") %>% 
    filter(gene_symbol %in% fig5F_df$gene_symbol) %>%
    select(gene_symbol, contains("neg"), contains("pos")) %>%
    scale_df(exclude_col = 1) %>%
    dplyr::rename("gene_symbol" = "V1")
  
  
  htmap_df <- tim3_TPM %>% 
    full_join(phago_TPM, by = "gene_symbol", 
              suffix = c(".tim3KO", "")) %>% 
    full_join(TGFBRII_FPKM, by = "gene_symbol", 
              suffix = c(".phago", ""))%>%
    full_join(Clec7a_FPKM, by = "gene_symbol", 
              suffix = c(".tgfbr2KO", ".clec7a")) %>% 
    mutate_at(2:ncol(.), as.numeric) %>%
    group_by(gene_symbol) %>% dplyr::slice(1) %>% ungroup 
  

  htmap_df_corr <- WGCNA::corAndPvalue(t(htmap_df[, -1]))
  
 
  fig5F_df_wide <- fig5F_df %>%
    select(-log2FoldChange) %>%
    distinct() %>%
    pivot_wider(names_from = "dataset", values_from = "direction") %>%
    select(gene_symbol, Tim3KO, `phago+`, Tgfbr2KO, `Clec7a+`) %>%
    mutate(gene_symbol = factor(gene_symbol, htmap_df$gene_symbol)) %>%
    arrange(gene_symbol) %>%
    mutate(dirs = paste(Tim3KO, `phago+`, Tgfbr2KO, `Clec7a+`, sep = "_")) %>%
    mutate(col = case_when(
      str_detect(dirs, "up")   & !str_detect(dirs, "down") ~ "red",
      str_detect(dirs, "down") & !str_detect(dirs, "up") ~ "blue",
      TRUE ~ "grey15")) %>%
    mutate(index = row_number()) 
  
  label_df <- fig5F_df_wide %>%
    filter(gene_symbol %in% gene_list)
  
  DEG_pal <- setNames(c("gold", "grey25"), factor(c("up", "down"), levels = c("up", "down")))
  row_annot <- rowAnnotation(
    `Tim3KO` = fig5F_df_wide$`Tim3KO`,
    `phago+` = fig5F_df_wide$`phago+`,
    `Tgfbr2KO` = fig5F_df_wide$`Tgfbr2KO`,
    `Clec7a+` = fig5F_df_wide$`Clec7a+`,
    col = list(`Tim3KO` = DEG_pal, `phago+` = DEG_pal,
               `Tgfbr2KO` = DEG_pal, `Clec7a+` = DEG_pal),
    link = anno_mark(at = label_df$index, labels = label_df$gene_symbol,
                     labels_gp = gpar(fontsize = 9.5, fontface = 4, col = label_df$col),
                     padding = unit(0.5, "mm")),
    annotation_label = c("(1)", "(2)", "(3)", "(4)", ""),
    annotation_name_side = "top",
    annotation_name_rot = 0,
    annotation_name_gp = gpar(fontsize = 9.5, fontface = 2),
    annotation_legend_param = annot_param,
    simple_anno_size = unit(3, "mm"),
    gap = unit(2, "points"),
    show_legend = TRUE) 
  
  col_annot <- columnAnnotation(
    `Tim3KO` = fig5F_df_wide$`Tim3KO`,
    `phago+` = fig5F_df_wide$`phago+`,
    `Tgfbr2KO` = fig5F_df_wide$`Tgfbr2KO`,
    `Clec7a+` = fig5F_df_wide$`Clec7a+`,
    col = list(`Tim3KO` = DEG_pal, `phago+` = DEG_pal,
               `Tgfbr2KO` = DEG_pal, `Clec7a+` = DEG_pal),
    annotation_label = c("(1)", "(2)", "(3)", "(4)"),
    annotation_name_gp = gpar(fontsize = 9.5, fontface = 2),
    annotation_legend_param = annot_param,
    simple_anno_size = unit(3, "mm"),
    gap = unit(2, "points"),
    show_legend = FALSE)

  
  cell_size <- 7/ncol(htmap_df_corr$cor)
  filename <- paste0("figures/Fig.4b_all_Tim3KO")
  p <- Heatmap(
    htmap_df_corr$cor, name = "Correlation",
    show_row_names = FALSE,
    show_row_dend = FALSE,
    show_column_names = FALSE,
    show_column_dend = FALSE,
    row_title = NULL,
    column_title = NULL,
    right_annotation = row_annot,
    bottom_annotation = col_annot,
    cluster_rows = TRUE,
    cluster_columns = TRUE
  )
  
  png(paste0(filename, suffix, ".png"), 
      width = 9.5, height = 7.5, units = "in", res = 400)
  draw(p)
  dev.off()
  pdf(paste0(filename, suffix, ".pdf"), 
      width = 9.5, height = 7.5)
  draw(p)
  dev.off()
  p <- draw(p)
  rowOrder <- unname(unlist(row_order(p)))
  rownames(htmap_df_corr$cor) <- htmap_df$gene_symbol
  colnames(htmap_df_corr$cor) <- htmap_df$gene_symbol
  cormat <- htmap_df_corr$cor[rowOrder, rowOrder] %>%
    as.data.frame %>%
    rownames_to_column("gene_symbol") %>%
    left_join(fig5F_df_wide, by = "gene_symbol")
  return(cormat)

}



all_tim3KO_and_atleast3 <- fig5_ls %>%
  data.table::rbindlist(idcol = "dataset") %>%
  arrange(dataset, direction) %>%
  group_by(gene_symbol) %>% mutate(n = n()) %>%
  filter(n >= 3) %>% ungroup %>% select(-n) %>% 
  rbind(cbind(dataset = "Tim3KO", fig5_ls$Tim3KO)) %>%
  mutate(direction = factor(direction, c("up", "down"))) %>%
  distinct()

htmap_corr_all_tim3_order <- htmap_corr_alltim3(
  fig5F_df = all_tim3KO_and_atleast3, suffix = "_and_atleast3",
  gene_list = c(gene_list_Tgfbr2, gene_list_phago, "Sall1", "Apoe", "Cd33"))

## save the orders of the genes
htmap_corr_all_tim3_order %>% 
  select(gene_symbol, names(fig5_ls)) %>%
  left_join(fig5_ls$Tim3KO[, -3], by = "gene_symbol") %>%
  left_join(fig5_ls$`phago+`[, -3], by = "gene_symbol",
            suffix = c(".Tim3KO", "")) %>%
  left_join(fig5_ls$Tgfbr2KO[, -3], by = "gene_symbol",
            suffix = c(".phago", "")) %>%
  left_join(fig5_ls$`Clec7a+`[, -3], by = "gene_symbol",
            suffix = c(".Tgfbr2KO", ".Clec7a")) %>% 
  write.csv("Fig.4b_corr_tim3KO_and_atleast3.csv", row.names = FALSE)
