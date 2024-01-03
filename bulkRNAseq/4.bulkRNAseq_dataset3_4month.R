library(dplyr)
library(ggplot2)
library(tibble)
library(stringr)
library(data.table)
library(DESeq2)
library(tidyr)
library(ComplexHeatmap)
library(ggrepel)



# Read data ----------------
meta_ds3_4M <- read.csv("data/expr_mat/GET_metadata.csv") %>%
  mutate(Cohort = str_extract(Cohort, "[0-9]+th"),
         Cohort = factor(Cohort, paste0(8:10, "th"))) %>%
  ## remove outlieres
  filter(!Plate %in% c("C04", "E04") & !str_detect(Plate, "11")) %>%
  mutate(Plate_col = str_extract(Plate, "[0-9]+")) %>%
  mutate(Group = case_when(Plate_col %in% c("06", "08", "09") ~ "Group_1",
                           Plate_col %in% c("03", "10", "12") ~ "Group_2",
                           Plate_col %in% c("01", "02", "04", "05", "07") ~ "Group_3")) %>%
  filter(Genotype != "Tim3_half_cKO.5XFAD") %>%
  mutate(Genotype = factor(Genotype, c("control", "Tim3_cKO", "5XFAD", "Tim3_cKO.5XFAD")))

count_ds3_4M <- read.csv("data/expr_mat/GET_count.csv") %>%
  select(gene_id, gene_symbol, meta_ds3_4M$Sample_ID)


meta_ds3_4M_tbl  <- meta_ds3_4M %>% column_to_rownames("Sample_ID")
count_ds3_4M_tbl <- count_ds3_4M %>% select(-gene_symbol) %>% column_to_rownames("gene_id")


# DESeq2 -------
get_dds <- function(counts_tbl, meta, design = "~ genotype", outliers = NULL, thres = 10){
  
  meta <- meta %>%
    filter(!Plate %in% outliers)
  
  # Construct DESeqDataSet object
  dds <- DESeqDataSetFromMatrix(
    countData = round(counts_tbl[, rownames(meta)]),
    colData = meta,
    design = as.formula(design))
  dds <- DESeq(dds)
  
  # Filter out genes with 0 counts
  keep <- rowSums(counts(dds)) >= thres
  dds <- dds[keep, ]
  dds
}



## DESeq2, adjusting for plate column ----------

## Construct DESeqDataSet
comparisons <- list(
  Tim3cKO_vs_control = c("Tim3_cKO", "control"),
  Tim3cKO.5XFAD_vs_5XFAD = c("Tim3_cKO.5XFAD", "5XFAD"),
  `5XFAD_vs_control` = c("5XFAD", "control"),
  Tim3cKO.5XFAD_vs_Tim3cKO = c("Tim3_cKO.5XFAD", "Tim3_cKO")
)

dds_ds3_F <- lapply(comparisons, function(x) {
  get_dds(counts_tbl = count_ds3_4M_tbl,
          meta = subset(meta_ds3_4M_tbl, Sex == "F" & Genotype %in% x),
          design = "~ Plate_col + Genotype", thres = 20)
})

dds_ds3_M <- lapply(comparisons, function(x) {
  get_dds(counts_tbl = count_ds3_4M_tbl,
          meta = subset(meta_ds3_4M_tbl, Sex == "M" & Genotype %in% x),
          design = "~ Plate_col + Genotype", thres = 20)
})

## Get results table
res_ds3_F <- sapply(names(comparisons), function(x) {
  lfcShrink(dds_ds3_F[[x]], type="apeglm", 
            coef = paste0("Genotype_", paste(comparisons[[x]], collapse = "_vs_")))
}, simplify = FALSE, USE.NAMES = TRUE)

res_ds3_M <- sapply(names(comparisons), function(x) {
  lfcShrink(dds_ds3_M[[x]], type="apeglm", 
            coef = paste0("Genotype_", paste(comparisons[[x]], collapse = "_vs_")))
}, simplify = FALSE, USE.NAMES = TRUE)

## Order results
order_results <- function(res) {
  as.data.frame(res) %>%
    tibble::rownames_to_column("gene_id") %>%
    left_join(count_ds3_4M[, c("gene_id", "gene_symbol")], by = "gene_id") %>%
    dplyr::select(gene_symbol, gene_id, everything()) %>%
    mutate(direction = ifelse(log2FoldChange > 0, 'up', 'down')) %>%
    mutate(direction = factor(direction, c("up", "down"))) %>%
    arrange(desc(direction), padj)
}

results_ds3_F <- lapply(res_ds3_F, order_results)
results_ds3_M <- lapply(res_ds3_M, order_results)

save(results_ds3_F, results_ds3_M, 
     file = "results/bulkRNAseq_results_ds3_4month.RData")

