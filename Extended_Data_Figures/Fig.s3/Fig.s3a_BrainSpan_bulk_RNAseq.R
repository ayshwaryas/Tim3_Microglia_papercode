library(tidyverse)
library(ComplexHeatmap)

# Load Data ----------------

## BrainSpan bulk RNA-seq dataset (https://www.brainspan.org/)
brainspan_col_meta <- read.csv("data/BrainSpan_genes_matrix_csv/columns_metadata.csv")
brainspan_row_meta <- read.csv("data/BrainSpan_genes_matrix_csv/rows_metadata.csv")
brainspan_expr <- read.csv("data/BrainSpan_genes_matrix_csv/expression_matrix.csv", header = FALSE, row.names = 1)

brainspan_donor_info <- brainspan_col_meta %>%
  mutate(age = factor(age, unique(.$age))) %>%
  group_by(donor_id, donor_name, age, gender) %>%
  dplyr::slice(1) %>% ungroup 

## Immune checkpoints and Tgfb-pathway related genes
checkpoint <- c(
  "HAVCR2", "LAG3", "C10orf54", "PDCD1", "CD274", "CTLA4", "TIGIT", "HAVCR1")


# Fig.s3c: Heatmap visualization the average expressions of selected checkpoint molecules in BrainSpan dataset ---------
for(gene in c(checkpoint)) {
  brainspan_row_gene <- subset(brainspan_row_meta, gene_symbol == gene) 
  if(nrow(brainspan_row_gene) == 0) {return(NULL)}
  
  brainspan_expr_gene <- brainspan_expr[brainspan_row_gene$row_num, ] %>%
    t() %>% `colnames<-`("expr") %>%
    cbind(brainspan_col_meta)
  
  gene_expr_df_mat <- brainspan_expr_gene %>%
    filter(str_detect(structure_name, "cortex")) %>%
    ## remove brain regions with <=5 samples
    group_by(structure_name) %>%
    filter(n() > 5) %>% ungroup %>%
    ## remove donors with samples from <= 5 brain regions
    group_by(donor_id) %>%
    filter(n() > 5) %>% ungroup %>%
    arrange(age, gender) %>%
    ## average gene expression in a certain brain region across samples of the same age 
    group_by(structure_acronym, age) %>%
    summarise(expr = mean(expr, na.rm = TRUE)) %>%
    mutate(age = factor(age, levels(brainspan_donor_info$age))) %>%
    arrange(age) %>%
    pivot_wider(names_from = c("age"), values_from = "expr") %>%
    arrange(structure_acronym) %>%
    column_to_rownames("structure_acronym")

  p_htmap <- log2(gene_expr_df_mat + 1) %>% 
    Heatmap(name = "log2 RPKM",
            col = rev(heat.colors(9)),
            show_column_dend = FALSE,
            cluster_column_slices = FALSE,
            cluster_rows = FALSE,
            row_names_side = "left",
            cluster_columns = FALSE,
            column_names_gp = gpar(fontsize = 10),
            column_names_rot = 45,
            column_names_side = "top",
            show_column_names = TRUE,
            row_names_gp = gpar(fontsize = 9))
  
  width <- ncol(gene_expr_df_mat) * 0.25
  height <- nrow(gene_expr_df_mat) * 0.23
  
  filename <- paste0("figures/brainspan/htmap_bygene_", gene, "_expr_unscaled_avg")
  png(paste0(filename, ".png"), width = width, height = height, res = 400, units = "in")
  draw(p_htmap, column_title = ifelse(gene == "C10orf54", "VSIR", gene),
       column_title_gp = gpar(fontsize = 9, fontface = 2), merge_legends = TRUE)
  dev.off()
  
  pdf(paste0(filename, ".pdf"), width = width, height = height)
  draw(p_htmap, column_title = ifelse(gene == "C10orf54", "VSIR", gene),
       column_title_gp = gpar(fontsize = 9, fontface = 2), merge_legends = TRUE)
  dev.off()
}

