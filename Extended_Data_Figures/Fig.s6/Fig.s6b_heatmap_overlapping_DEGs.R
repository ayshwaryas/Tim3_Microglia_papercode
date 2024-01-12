source('bulkRNAseq/0.bulkRNAseq_functions.R')

## Palettes
pal <- c(brewer.pal(9, "Greens")[4], "forestgreen", "steelblue1", "steelblue4",
         "darkgoldenrod1", "orangered", "#FB9A99", "darkred")
pal_text <- c(brewer.pal(9, "Greens")[6], "forestgreen", "steelblue2", "steelblue4",
              "darkgoldenrod3", "orangered", "#fa7573", "darkred")


# Read data --------------------------------------------------------------------

## 1. Differential expression analysis results ---------------------------------
## (1) 1-month-old mice, Havcr2cKO vs Havcr2flox/flox --------------------------
load("results/bulkRNAseq_results_ds2_1month.RData")
tim3_DEG <- res_ordered$`1M` %>% correct_gene_symbol() %>%
  select(gene_symbol, log2FoldChange, direction, padj) %>% 
  res_order_slice(thres = 0.1)

## (2) 3-month-old mice, phagocytosing control vs non-phagocytosing control microglia ---------
load("results/bulkRNAseq_results_ds1_batch1_3month.RData")
phago_DEG <- results_batch1$phagoposvsneg_control %>%
  dplyr::rename("gene_symbol" = "gene_name") %>%
  correct_gene_symbol() %>%
  select(gene_symbol, log2FoldChange, direction, padj) %>% 
  res_order_slice(thres = 0.1)

## (3) Tgfbr2cKo vs control (Lund et al. 2018, PMID: 29662171) -----------------
TGFBRII_DEG <- read.csv("results/bulkRNAseq_results_TGFBRII_Lund_2018.csv") %>%
  select(gene_symbol, contains(".uG")) %>%
  dplyr::rename("log2FoldChange" = "log2fc.uG",
                "direction" = "dir.uG", "padj" = "padj.uG") %>%
  filter(rowSums(.[, 2:7]) != 0) %>%
  correct_gene_symbol()  %>%
  select(gene_symbol, log2FoldChange, direction, padj) %>%
  res_order_slice(thres = 0.1) %>%
  distinct(gene_symbol, direction, .keep_all = TRUE)

## (4) Clec7a+ vs Clec7a- (Krasemann et al., 2017, PMID: 28930663) -------------
Clec7a_DEG <- read.csv("data/Oleg_Immunity_2017_DEG_ADpos_vs_neg_all.csv") %>%
  dplyr::rename("gene_symbol" = "tracking_id", "log2FoldChange" = "log2FC") %>%
  correct_gene_symbol() %>%
  select(gene_symbol, log2FoldChange, direction, padj) %>% 
  res_order_slice(thres = 0.1)



## Direction of change of DEGs shared by at least 3 datasets -------------------

Figs6b_DEG_dir_wide <- list(`Tim3KO` = tim3_DEG, `phago+` = phago_DEG,
                            `Tgfbr2KO` = TGFBRII_DEG, `Clec7a+`= Clec7a_DEG) %>%
  data.table::rbindlist(idcol = "dataset") %>%
  mutate(direction = factor(direction, c("up", "down"))) %>%
  arrange(dataset, direction) %>%
  group_by(gene_symbol) %>%
  filter(n() >= 3) %>% ungroup %>% 
  select(-log2FoldChange) %>%
  pivot_wider(names_from = "dataset", values_from = "direction") %>%
  select(gene_symbol, Tim3KO, `phago+`, Tgfbr2KO, `Clec7a+`)

## 2. Expression matrix and metadata -------------------------------------------

## Function for scaling the expression matrix
scale_df <- function(df, exclude_col = 1:2) {
  df_scaled <- t(df[, -exclude_col]) %>% scale() %>% t()
  cbind(df[, exclude_col], df_scaled) %>% as_tibble()
}

## (1) 1-month-old mice, Havcr2cKO vs Havcr2flox/flox --------------------------
Tim3KO_TPM <- read.csv("data/expr_mat/Danyang_TPM_matrix.csv", check.names = FALSE) %>%
  `colnames<-`(c("gene_symbol", colnames(.)[-1])) %>%
  inner_join(Figs6b_DEG_dir_wide %>% select(gene_symbol, direction = Tim3KO), 
             by = "gene_symbol") %>%
  select(gene_symbol, direction, contains("1M")) %>%
  select(gene_symbol, direction, contains("Tim3_flox"), contains("Tim3_cKO")) %>%
  scale_df() 

## (2) 3-month-old mice, phagocytosing control vs non-phagocytosing control microglia ---------
## metadata
phago_meta <- read.csv("data/expr_mat/tim3_Kimi_metadata_KK.csv") %>%
  filter(Batch == "Batch_1") %>%
  filter(Sequencing_run != "excluded") %>%
  filter(genotype == "control") %>%
  filter(!Sample_number %in% c("KK1", "KK2", "KK7", "KK9", "KK16")) %>%
  arrange(cell_type) %>%
  mutate(Sample_ID_new = paste(Batch, Sample_number, cell_type, sep = "_")) %>%
  mutate(Sample_ID_new = str_replace(Sample_ID_new, "_minus", "-"),
         Sample_ID_new = str_replace(Sample_ID_new, "_plus", "+"))

## expression matrix
phago_TPM <- read.csv("data/expr_mat/tim3_Kimi_tpm_with_genesymbols_KK.csv") %>%
  filter(gene_id %in% results_batch1$phagoposvsneg_control$gene_id) %>%
  select(gene_symbol = gene_name, phago_meta$Sample_ID) %>% 
  `colnames<-`(c("gene_symbol", phago_meta$Sample_ID_new)) %>%
  inner_join(Figs6b_DEG_dir_wide %>% select(gene_symbol, direction = `phago+`), 
             by = "gene_symbol") %>%
  select(gene_symbol, direction, everything())  %>%
  scale_df()


## (3) Tgfbr2cKo vs control (Lund et al. 2018, PMID: 29662171) -----------------
TGFBRII_FPKM <- read.csv("data/FPKM OB edits - with padj.csv") %>%
  dplyr::rename("gene_symbol" = "geneNames") %>% 
  inner_join(Figs6b_DEG_dir_wide %>% select(gene_symbol, direction = Tgfbr2KO), 
             by = "gene_symbol") %>%
  select(gene_symbol, direction, contains("WT.uG"), contains("KO.uG")) %>%
  scale_df() 

## (4) Clec7a+ vs Clec7a- (Krasemann et al., 2017, PMID: 28930663) -------------
Clec7a_FPKM <- read.csv("data/Oleg_Immunity_2017_DEG_ADpos_vs_neg_all.csv") %>%
  dplyr::rename("gene_symbol" = "tracking_id") %>% 
  select(-direction) %>%
  inner_join(Figs6b_DEG_dir_wide %>% select(gene_symbol, direction = `Clec7a+`), 
             by = "gene_symbol") %>% ungroup %>%
  select(gene_symbol, direction, contains("neg"), contains("pos")) %>%
  scale_df()


## Joining the expression matrices ---------------------------------------------
Figs6b_expr_mat <- Tim3KO_TPM %>%
  full_join(phago_TPM, by = "gene_symbol", suffix = c(".tim3KO", "")) %>% 
  full_join(TGFBRII_FPKM, by = "gene_symbol", suffix = c(".phago", "")) %>%
  full_join(Clec7a_FPKM, by = "gene_symbol", suffix = c(".tgfbr2KO", ".clec7a")) %>% 
  distinct() %>%
  mutate(dirs = paste(direction.tim3KO, direction.phago, 
                      direction.tgfbr2KO, direction.clec7a, sep = "_")) %>%
  select(-contains("direction")) %>%
  select(gene_symbol, dirs, 
         contains("Tim3_flox"), contains("phago-"), contains("WT.uG"), contains("neg"), 
         contains("Tim3_cKO"), contains('phago+'), contains("KO.uG"), contains("pos"),
         everything()) 

# Heatmap ----------------------------------------------------------------------



## Main body of the heatmap ----------------------------------------------------
htmap_df <- Figs6b_expr_mat %>%
  group_by(gene_symbol) %>% dplyr::slice(1) %>% ungroup %>%
  separate(dirs, into = c("dir_Tim3KO", "dir_phago", "dir_Tgfbr2KO", "dir_Clec7a"),
           sep = "_", remove = FALSE) %>%
  mutate_at(vars(starts_with("dir_")),
            function(x) {case_when(x == "up" ~ 1, x == "down" ~ -1, TRUE ~ 0)}) %>%
  mutate(dirs_num = dir_Tim3KO + dir_phago + dir_Tgfbr2KO + dir_Clec7a) %>%
  mutate(dirs_num_abs = abs(dir_Tim3KO) + abs(dir_phago) + abs(dir_Tgfbr2KO) + abs(dir_Clec7a)) %>%
  arrange(desc(dirs_num), desc(dirs)) %>%
  mutate(dirs = factor(dirs, unique(.$dirs))) %>%
  select(-starts_with("dir_"), -dirs_num, -dirs_num_abs)



## Row annotations (direction of change of the 4 comparisons) ------------------
## palette
dir_pal <- setNames(c("gold", "grey25"), factor(c("up", "down"), levels = c("up", "down")))

## dataframe for drawing the annotations 
Figs6b_DEG_dir_wide_ordered <- Figs6b_DEG_dir_wide %>%
  mutate(dirs = paste(Tim3KO, `phago+`, Tgfbr2KO, `Clec7a+`, sep = "_")) %>%
  mutate(col = case_when(
    str_detect(dirs, "up")   & !str_detect(dirs, "down") ~ "red",
    str_detect(dirs, "down") & !str_detect(dirs, "up") ~ "blue",
    TRUE ~ "grey15")) %>%
  mutate(gene_symbol = factor(gene_symbol, unique(htmap_df$gene_symbol))) %>%
  arrange(gene_symbol) %>%
  mutate(index = row_number())

## genes to highlight
gene_list_Tgfbr2 <- c("Havcr2", "Axl", "Ly9", "Clec7a", "Cd63", "Cxcl16", 
                      "H2-K1", "H2-D1", "Siglech", "Csf1r")
gene_list_phago <- c("Havcr2", "Cd9", "Cst7", "Cstb", "Cxcl16", "Ly9", "Lyz2", 
                     unique(subset(Figs6b_DEG_dir_wide, str_detect(gene_symbol, "H2\\-"))$gene_symbol))

label_df <- Figs6b_DEG_dir_wide_ordered %>%
  filter(gene_symbol %in% c(gene_list_Tgfbr2, gene_list_phago, "Sall1", "Apoe", "Cd33"))

## legend 
row_annot_legend <- lapply(
  list(`Tim3KO`   = "(1) ***Havcr2***<sup>cKO</sup> vs control",
       `phago+`   = "(2) Phagocytosing vs<br/>Non-phagocytosing MG",
       `Tgfbr2KO` = "(3) ***Tgfbr2***<sup>cKO</sup> vs control",
       `Clec7a+`  = "(4) ***Clec7a***<sup>+</sup> vs ***Clec7a***<sup>-</sup>"),
  function(x) { list(title = gt_render(x),
                     title_gp = gpar(fontsize = 10, fontface = 2, lineheight = 1.15))})

## annotation
row_annot <- rowAnnotation(
  `Tim3KO` = Figs6b_DEG_dir_wide_ordered$`Tim3KO`,
  `phago+` = Figs6b_DEG_dir_wide_ordered$`phago+`,
  `Tgfbr2KO` = Figs6b_DEG_dir_wide_ordered$`Tgfbr2KO`,
  `Clec7a+` = Figs6b_DEG_dir_wide_ordered$`Clec7a+`,
  link = anno_mark(at = label_df$index, labels = label_df$gene_symbol,
                   labels_gp = gpar(fontsize = 9, fontface = 4, col = label_df$col),
                   padding = unit(0.5, "mm"), link_width = unit(6, "mm")),
  col = list(`Tim3KO` = dir_pal, `phago+` = dir_pal,
             `Tgfbr2KO` = dir_pal, `Clec7a+` = dir_pal),
  annotation_label = gt_render(c("(1)", "(2)", "(3)", "(4)", "")),
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontsize = 9.5, fontface = 2),
  annotation_name_rot = 0,
  annotation_legend_param = row_annot_legend,
  simple_anno_size = unit(3, "mm"),
  gap = unit(2, "points")
)

## Column annotations -----------------------------------------------------------

## Column split 
col_split <- colnames(htmap_df)[-(1:2)] %>%
  str_extract("Tim3_flox|Tim3_cKO|phago\\+|phago\\-|WT\\.uG|KO\\.uG|AD.*pos|AD.*neg") %>%
  str_replace("Tim3_cKO", "Tim3KO") %>%
  str_replace("Tim3_flox", "control") %>%
  str_replace("AD.*pos", "Clec7a+") %>%
  str_replace("AD.*neg", "Clec7a-") %>%
  str_replace("WT\\.uG", "WT") %>%
  str_replace("KO\\.uG", "Tgfbr2KO") %>%
  factor(unique(.))

## annotation
col_annot <- columnAnnotation(
  col_names = anno_block(
    gp = gpar(col = NA, fill = NA),
    labels = gt_render(c(
      "control for<br/>***Havcr2***<sup>cKO</sup>", "Non-phago-<br/>cytosing MG",
      "control for<br/>***Tgfbr2***<sup>cKO</sup>", "***Clec7a***<sup>-</sup>",
      "***Havcr2***<sup>cKO</sup>", "Phago-<br/>cytosing MG",
      "***Tgfbr2***<sup>cKO</sup>", "***Clec7a***<sup>+</sup>")),
    labels_gp = gpar(fontsize = 10, fontface = 2, lineheight = 1.15, col = pal_text),
    labels_rot = 45, labels_just ="left", labels_offset = unit(0.1, "npc"), height = unit(8, "mm")
  ),
  genotype = col_split,
  col = list(genotype = setNames(pal, levels(col_split))),
  annotation_name_gp = gpar(fontsize = 0),
  simple_anno_size = unit(3, "mm"),
  show_legend = FALSE)


cell_size <- 7/ncol(htmap_df)
p <-  htmap_df %>% select(-gene_symbol, -dirs) %>%
  Heatmap(
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    top_annotation = col_annot,
    right_annotation = row_annot,
    cluster_columns = TRUE,
    column_split = col_split,
    cluster_column_slices = FALSE,
    row_split = htmap_df$dirs,
    cluster_row_slices = FALSE,
    row_gap = unit(
      c(rep(0, sum(!str_detect(levels(htmap_df$dirs), "down")) - 1), 1,
        rep(0, sum(str_detect(levels(htmap_df$dirs), "up") & str_detect(levels(htmap_df$dirs), "down"))- 1), 1,
        rep(0, sum(!str_detect(levels(htmap_df$dirs), "up")) - 1)),
      "mm"),
    column_gap = unit(c(rep(1, 3), 2, rep(1, 3)), "mm"),
    row_title = NULL,
    column_title = NULL,
    name = "Z-score"
  )
png(paste0("figures/Fig.s6b_overlapping_DEGs_expr_heatmap.png"), width = 10, height = 8.5, units = "in", res = 400)
draw(p, padding = unit(c(2, 2, 10, 2), units = "mm"))
dev.off()
pdf(paste0("figures/Fig.s6b_overlapping_DEGs_expr_heatmap.pdf"), width = 10, height = 8.5)
draw(p)
dev.off()
