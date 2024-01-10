
source('bulkRNAseq/0.bulkRNAseq_functions.R')

# Read data -----------------------------------------

## (1) 1-month-old mice, Havcr2cKO vs Havcr2flox/flox
load("results/bulkRNAseq_results_ds2_1month.RData")
tim3_all <- res_ordered$`1M` %>% correct_gene_symbol() %>%
  select(gene_symbol, log2FoldChange, direction, padj)
tim3_DEG <- tim3_all %>% res_order_slice

## (2) 3-month-old mice, phagocytosing control vs non-phagocytosing control microglia
load("results/bulkRNAseq_results_ds1_batch1_3month.RData")
phago_all <- results_batch1$phagoposvsneg_control %>%
  dplyr::rename("gene_symbol" = "gene_name") %>%
  correct_gene_symbol() %>%
  select(gene_symbol, log2FoldChange, direction, padj)
phago_DEG <- phago_all %>% res_order_slice(flip = FALSE)

## (3) Tgfbr2cKo vs control (Lund et al. 2018, PMID: 29662171)
TGFBRII_all <- read.csv("results/bulkRNAseq_results_TGFBRII_Lund_2018.csv") %>%
  select(gene_symbol, contains(".uG")) %>%
  dplyr::rename("log2FoldChange" = "log2fc.uG",
                "direction" = "dir.uG", "padj" = "padj.uG") %>%
  filter(rowSums(.[, 2:7]) != 0) %>%
  correct_gene_symbol()  %>%
  select(gene_symbol, log2FoldChange, direction, padj)

TGFBRII_DEG <- TGFBRII_all %>% res_order_slice(flip = FALSE, thres = 0.1)

## (4) Clec7a+ vs Clec7a- (Krasemann et al., 2017, PMID: 28930663)
Clec7a_all <- read.csv("data/Oleg_Immunity_2017_DEG_ADpos_vs_neg_all.csv") %>%
  dplyr::rename("gene_symbol" = "tracking_id", "log2FoldChange" = "log2FC") %>%
  correct_gene_symbol() %>%
  select(gene_symbol, log2FoldChange, direction, padj)
Clec7a_DEG <- Clec7a_all %>% res_order_slice(flip = FALSE)


## Combine Datasets 
Fig4b_DEG_ls <- list(`Tim3KO` = tim3_DEG, `phago+` = phago_DEG,
                     `Tgfbr2KO` = TGFBRII_DEG, `Clec7a+`= Clec7a_DEG)

Fig4b_DEG_df <- Fig4b_DEG_ls %>%
  data.table::rbindlist(idcol = "dataset") %>%
  mutate(direction = factor(direction, c("up", "down"))) %>%
  arrange(dataset, direction)

gene_list_Tgfbr2 <- c("Havcr2", "Axl", "Ly9", "Clec7a", "Cd63", "Cxcl16", 
                      "H2-K1", "H2-D1", "Siglech", "Csf1r")
gene_list_phago <- c("Havcr2", "Cd9", "Cst7", "Cstb", "Cxcl16", "Ly9", "Lyz2", 
                     unique(subset(Fig4b_DEG_df, str_detect(gene_symbol, "H2\\-"))$gene_symbol))

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

htmap_corr_alltim3 <- function(DEG_df = NULL, n_min, suffix = "", gene_list = NULL) {
  
  ## (1) TPM, 1-month-old mice, Havcr2cKO vs Havcr2flox/flox
  tim3_TPM <- read.csv("data/expr_mat/Danyang_TPM_matrix.csv", check.names = FALSE) %>%
    `colnames<-`(c("gene_symbol", colnames(.)[-1])) %>%
    filter(gene_symbol %in% DEG_df$gene_symbol) %>%
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
    filter(gene_id %in% results_batch1$res_phagoposvsneg_control$gene_id) %>%
    select(gene_symbol = gene_name, phago_meta$Sample_ID) %>%
    filter(gene_symbol %in% DEG_df$gene_symbol) %>%
    select(gene_symbol, everything())  %>%
    scale_df(exclude_col = 1) %>%
    dplyr::rename("gene_symbol" = "V1")
  
  
  ## (3) FPKM, Tgfbr2cKo vs control (Lund et al. 2018, PMID: 29662171)
  TGFBRII_FPKM <- read.csv("data/FPKM OB edits - with padj.csv") %>%
    dplyr::rename("gene_symbol" = "geneNames") %>% 
    filter(gene_symbol %in% DEG_df$gene_symbol) %>%
    select(gene_symbol, contains("WT.uG"), contains("KO.uG")) %>%
    scale_df(exclude_col = 1) %>%
    dplyr::rename("gene_symbol" = "V1")
  
  ## (4) FPKM, Clec7a+ vs Clec7a- (Krasemann et al., 2017, PMID: 28930663)
  Clec7a_FPKM <- read.csv("data/Oleg_Immunity_2017_DEG_ADpos_vs_neg_all.csv") %>%
    dplyr::rename("gene_symbol" = "tracking_id") %>% 
    filter(gene_symbol %in% DEG_df$gene_symbol) %>%
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
  
 
  DEG_df_wide <- DEG_df %>%
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
  
  label_df <- DEG_df_wide %>%
    filter(gene_symbol %in% gene_list)
  
  DEG_pal <- setNames(c("gold", "grey25"), factor(c("up", "down"), levels = c("up", "down")))
  row_annot <- rowAnnotation(
    `Tim3KO` = DEG_df_wide$`Tim3KO`,
    `phago+` = DEG_df_wide$`phago+`,
    `Tgfbr2KO` = DEG_df_wide$`Tgfbr2KO`,
    `Clec7a+` = DEG_df_wide$`Clec7a+`,
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
    `Tim3KO` = DEG_df_wide$`Tim3KO`,
    `phago+` = DEG_df_wide$`phago+`,
    `Tgfbr2KO` = DEG_df_wide$`Tgfbr2KO`,
    `Clec7a+` = DEG_df_wide$`Clec7a+`,
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
    left_join(DEG_df_wide, by = "gene_symbol")
  return(cormat)

}

all_tim3KO_and_atleast3 <- Fig4b_DEG_ls %>%
  data.table::rbindlist(idcol = "dataset") %>%
  arrange(dataset, direction) %>%
  group_by(gene_symbol) %>% mutate(n = n()) %>%
  filter(n >= 3) %>% ungroup %>% select(-n) %>% 
  rbind(cbind(dataset = "Tim3KO", Fig4b_DEG_ls$Tim3KO)) %>%
  mutate(direction = factor(direction, c("up", "down"))) %>%
  distinct()

htmap_corr_all_tim3_order <- htmap_corr_alltim3(
  DEG_df = all_tim3KO_and_atleast3, suffix = "_and_atleast3",
  gene_list = c(gene_list_Tgfbr2, gene_list_phago, "Sall1", "Apoe", "Cd33"))

