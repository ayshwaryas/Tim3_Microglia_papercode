library(tidyverse)
library(DESeq2)

# Load count matrix and metadata table ------------------------
Tim3KO_cnts <- read.csv("data/expr_mat/tim3_Kimi_count_with_genesymbols_KK.csv")
Tim3KO_metadata <- read.csv("data/expr_mat/tim3_Kimi_metadata_KK.csv") %>%
  filter(Sample_ID %in% colnames(Tim3KO_cnts)[-(1:2)]) %>%
  mutate(cell_type.genotype = paste(cell_type, genotype, Sample_number, sep = "_")) %>%
  tibble::column_to_rownames("cell_type.genotype") 

# Subset batch 2 samples from metadata and count matrix ---------------
outliers_batch2 <- c("KK1", "KK6", "KK7")

Tim3KO_meta_batch2 <- Tim3KO_metadata %>%
  filter(Batch == "Batch_2" & !(Sample_number %in% outliers_batch2)) %>% 
  mutate(genotype = factor(genotype, c("control", "Tim3.cKO", "5XFAD", "Tim3.cKO.5XFAD"))) %>%
  arrange(genotype) %>%
  `rownames<-`(.$Sample_ID)

Tim3KO_cnts_batch2 <- Tim3KO_cnts %>%
  select(gene_id, Tim3KO_meta_batch2$Sample_ID) %>%
  column_to_rownames("gene_id")

all(colnames(Tim3KO_cnts_batch2) == rownames(Tim3KO_meta_batch2))


# DESeq2 -------------------------------------

## Construct DESeqDataSet
get_dds <- function(counts_tbl, meta, design = "~ genotype", outliers = NULL, thres = 10){
  
  meta <- meta %>%
    filter(!Sample_number %in% outliers)
  
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

comparisons <- list(
  Tim3cKO_vs_control = c("Tim3.cKO", "control"),
  Tim3cKO.5XFAD_vs_5XFAD = c("Tim3.cKO.5XFAD", "5XFAD"),
  `5XFAD_vs_control` = c("5XFAD", "control"),
  Tim3cKO.5XFAD_vs_Tim3cKO = c("Tim3.cKO.5XFAD", "Tim3.cKO")
)

dds_batch2_F <- lapply(comparisons, function(x) {
  get_dds(counts_tbl = Tim3KO_cnts_batch2,
          meta = subset(Tim3KO_meta_batch2, sex == "F" & genotype %in% x),
          design = "~ genotype", outliers = outliers_batch2)
})

dds_batch2_M <- lapply(comparisons, function(x) {
  get_dds(counts_tbl = Tim3KO_cnts_batch2,
          meta = subset(Tim3KO_meta_batch2, sex == "M" & genotype %in% x),
          design = "~ genotype", outliers = outliers_batch2)
})

## Get results table
res_batch2_F <- sapply(names(comparisons), function(x) {
  lfcShrink(dds_batch2_F[[x]], type="apeglm", 
            coef = paste0("genotype_", paste(comparisons[[x]], collapse = "_vs_")))
}, simplify = FALSE, USE.NAMES = TRUE)

res_batch2_M <- sapply(names(comparisons), function(x) {
  lfcShrink(dds_batch2_M[[x]], type="apeglm", 
            coef = paste0("genotype_", paste(comparisons[[x]], collapse = "_vs_")))
}, simplify = FALSE, USE.NAMES = TRUE)

## Order results
res_order <- function(res) {
  as.data.frame(res) %>%
    tibble::rownames_to_column("gene_id") %>%
    left_join(Tim3KO_cnts[, c("gene_id", "gene_name")], by = "gene_id") %>%
    select(gene_id, gene_name, everything()) %>%
    arrange(padj)  %>%
    mutate(direction = ifelse(log2FoldChange > 0, "up", "down")) %>%
    mutate(direction = factor(direction, levels = c("up", "down")))
}

results_batch2_F <- lapply(res_batch2_F, res_order)
results_batch2_M <- lapply(res_batch2_M, res_order)

save(results_batch2_F, results_batch2_M, 
     file = "results/bulkRNAseq_results_ds1_batch2_7month.RData")
