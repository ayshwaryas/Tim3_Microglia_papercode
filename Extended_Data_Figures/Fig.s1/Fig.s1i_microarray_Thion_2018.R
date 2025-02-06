library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggsci)

# Read microarray dataset --------------------
microarray_raw <- read.table("data/microarray/GSE107129_mouse_normalized.txt", sep = "\t", header = TRUE) 
microarray_ref <- readxl::read_xlsx("data/microarray/Thion_2018_mmc1.xlsx", sheet = 1, skip = 2) %>%
  select(IlluminaID, Symbol) %>% distinct()

microarray <- microarray_raw %>%
  select(-starts_with("X")) %>%
  left_join(microarray_ref, by = c("ID_REF" = "IlluminaID")) %>%
  select(ID_REF, Symbol, everything())

meta_microarray <- data.frame(Sample_ID = colnames(microarray)[-c(1:2)]) %>%
  filter(!str_detect(Sample_ID, "A2M.WT.[137]_2$")) %>%
  mutate(Age = str_extract(Sample_ID, "YSM")) %>%
  mutate(Age = ifelse(is.na(Age), str_extract(Sample_ID, "E[0-9]+\\.5|NB|A2M"), Age)) %>%
  mutate(Age = case_when(Age == "NB" ~ "P0",
                         Age == "A2M" ~ "Adult",
                         TRUE ~ Age)) %>%
  mutate(Sample_number = str_extract(Sample_ID, "(?<=\\.)[0-9\\._]+$")) %>%
  mutate(Age = factor(Age, c("YSM", paste0("E", seq(10, 18, 2), ".5"), "P0", "Adult"))) %>%
  arrange(Age, Sample_number)  

microarray_sub <- microarray %>%
  mutate(Symbol = ifelse(Symbol == "4632428N05Rik", "Vsir", Symbol)) %>%
  select(Symbol, meta_microarray$Sample_ID) 

# Genes of interest -------------------------------------------------------
## Well-known MGnD signature genes
MGnD_sig <- c("Axl", "Clec7a", "Trem2", "Itgax", "Ly9", "Cst7")

## Well-known Homeostatic signature genes
Homeostasis_sig <- c("Sall1", "P2ry12", "Tmem119", "Hexb", "Fcrls", "Siglech")

## Top 10 Clec7a+ vs Clec7a- DEGs
## (up-regulated in Clec7a+: MGnD; down-regulated in Clec7a+: Homeostatic)
AD_pos_vs_neg_sig <- read.csv("results/bulkRNAseq_results_Clec7a_Krasemann_2017.csv")
Clec7a_top10 <- AD_pos_vs_neg_sig %>%
  arrange(padj) %>%
  filter(tracking_id %in% microarray$Symbol) %>%
  group_by(direction) %>%
  dplyr::slice(1:10) %>% ungroup %>%
  mutate(type = ifelse(direction == "up", "UP in Clec7a+", "DN in Clec7a+")) %>%
  select(Symbol = tracking_id, type)

MGnD_Homeo <- Clec7a_top10 %>% 
  rbind(data.frame(Symbol = MGnD_sig, type = "MGnD"), 
        data.frame(Symbol = Homeostasis_sig, type = "Homeostasis")) %>%
  mutate(type = factor(type, c("MGnD", "UP in Clec7a+", "Homeostasis", "DN in Clec7a+"))) %>%
  mutate(color = case_when(type == "MGnD" ~ "brown3",
                           type == "UP in Clec7a+" ~ "#FBB4AE",
                           type == "Homeostasis" ~ "grey30",
                           type == "DN in Clec7a+" ~ "grey80"))


# Heatmap -----------------------
Heatmap_custom <- function(gene_list, suffix = "", width = 6.5, height = 4,
                           rownames_col_df = NULL, type_pal = NULL, 
                           cluster_rows = FALSE, fig_num = "1A") {
  htmap_df <- microarray_sub %>% 
    filter(Symbol %in% gene_list)
  
  htmap_df <- as.matrix(htmap_df[, -1]) %>%
    `rownames<-`(htmap_df$Symbol) %>%
    limma::avereps() %>% as.data.frame() %>%
    rownames_to_column("Symbol") %>%
    mutate(Symbol = factor(Symbol, gene_list)) %>%
    arrange(Symbol)
  
  if(is.null(rownames_col_df)) { 
    row_split <- NULL; row_annot <- NULL
  }
  else {
    htmap_df <- htmap_df %>%
      left_join(rownames_col_df, by = "Symbol") 
    row_split <- htmap_df$type
    row_annot <- rowAnnotation(
      type = htmap_df$type,
      col = list(type = setNames(unique(MGnD_Homeo$color), unique(MGnD_Homeo$type))),
      simple_anno_size = unit(4, "mm"),
      show_annotation_name = FALSE
    )
  }
  
  htmap_df_scaled <- htmap_df %>%
    select(meta_microarray$Sample_ID) %>%
    t() %>% scale() %>% t()
  
  htmap_col <- circlize::colorRamp2(
    c(-3, 0, 3), c("navy", "white", "firebrick2"))
  
  col_annot <- columnAnnotation(
    Age = anno_block(
      gp = gpar(fill = c("#FDBF6F", "#FF7F00", "#B2DF8A", "#33A02C", 
                         "#A6CEE3", "#1F78B4", "#08306B", "darkred")),
      height = unit(6, "mm"),
      labels_gp = gpar(col = c(rep("black", 5), rep("white", 3)),
                       fontsize = 11, fontface = "bold"),
      labels = unique(meta_microarray$Age)),
    annotation_name_side = "left"
  )
  
  p <- Heatmap(
    htmap_df_scaled, 
    name = "Z-score",
    col = htmap_col,
    row_title = NULL,
    row_labels = htmap_df$Symbol,
    row_names_gp = gpar(fontface = "italic", fontsize = 12),
    right_annotation = row_annot,
    row_split = row_split,
    row_gap = unit(0, "mm"),
    show_row_dend = FALSE,
    cluster_rows = cluster_rows,
    cluster_row_slices = FALSE,
    cluster_column_slices = FALSE,
    show_column_names = FALSE,
    column_split = meta_microarray$Age,
    column_names_rot = 0,
    column_names_centered = TRUE,
    column_title = NULL,
    column_gap = unit(1, "mm"),
    top_annotation = col_annot,
    show_column_dend = FALSE
  )
  
  filename <- paste0("Fig.", fig_num, "_htmap_", suffix)
  png(paste0("figures/", filename, ".png"), width = width, height = height,
      res = 400, units = "in")
  draw(p)
  dev.off()
  pdf(paste0("figures/", filename, ".pdf"), width = width, height = height)
  draw(p)
  dev.off()
  
  ## Save source data
  p <- draw(p)
  row_order <- row_order(p) %>% unlist
  source_data <- htmap_df %>% 
    dplyr::select("Symbol", meta_microarray$Sample_ID) %>%
    column_to_rownames("Symbol") %>% t() %>% scale() %>%
    as.data.frame() %>%
    dplyr::select(htmap_df$Symbol[row_order])
  write.csv(source_data, paste0("Source_Data/", filename, ".csv"))
}


## Extended Data Figure 1i, expressions of MGnD and Homeostatic-related genes
Heatmap_custom(MGnD_Homeo$Symbol, suffix = "MGnD_Homeo", fig_num = "s3c",
               rownames_col_df = MGnD_Homeo,
               width = 7, height = 7, cluster_rows = TRUE)

