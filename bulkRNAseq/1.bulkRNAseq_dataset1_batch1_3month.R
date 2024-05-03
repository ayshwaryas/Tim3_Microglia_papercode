library(tidyverse)
library(DESeq2)

# Load count matrix and metadata table ------------------------
Tim3KO_cnts <- read.csv("data/expr_mat/tim3_Kimi_count_with_genesymbols_KK.csv")
Tim3KO_metadata <- read.csv("data/expr_mat/tim3_Kimi_metadata_KK.csv") %>%
  filter(Sample_ID %in% colnames(Tim3KO_cnts)[-(1:2)]) %>%
  mutate(cell_type.genotype = paste(cell_type, genotype, Sample_number, sep = "_")) %>%
  tibble::column_to_rownames("cell_type.genotype") 

# Subset batch 1 samples from metadata and count matrix ---------------
outliers_batch1 <- c("KK1", "KK2", "KK7", "KK9", "KK16")

Tim3KO_meta_batch1 <- Tim3KO_metadata %>%
  filter(Batch == "Batch_1" & !(Sample_number %in% outliers_batch1)) %>%
  mutate(Sample_num = str_extract(Sample_number, "[0-9]+") %>% as.numeric())%>%
  mutate(cell_type = ifelse(cell_type == "phago_plus", "phagopos", "phagoneg")) %>%
  mutate(cell_type = factor(cell_type, c("phagoneg", "phagopos"))) %>%
  mutate(genotype = factor(genotype, c("control", "Tim3.cKO"))) %>%
  arrange(cell_type, genotype)

Tim3KO_cnts_batch1 <- Tim3KO_cnts %>%
  select(gene_id, Tim3KO_meta_batch1$Sample_ID) %>%
  column_to_rownames("gene_id") %>%
  `colnames<-`(rownames(Tim3KO_meta_batch1))


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


dds_Tim3cKOvscontrol_phagopos <- get_dds(
  counts_tbl = Tim3KO_cnts_batch1,
  meta = subset(Tim3KO_meta_batch1, cell_type == "phagopos"),
  design = "~ genotype", outliers = outliers_batch1
)

dds_Tim3cKOvscontrol_phagoneg <- get_dds(
  counts_tbl = Tim3KO_cnts_batch1,
  meta = subset(Tim3KO_meta_batch1, cell_type == "phagoneg"),
  design = "~ genotype", outliers = outliers_batch1
)

dds_phagoposvsneg_control <- get_dds(
  counts_tbl = Tim3KO_cnts_batch1,
  meta = subset(Tim3KO_meta_batch1, genotype == "control"),
  design = "~ cell_type", outliers = outliers_batch1
)

dds_phagoposvsneg_Tim3.cKO <- get_dds(
  counts_tbl = Tim3KO_cnts_batch1,
  meta = subset(Tim3KO_meta_batch1, genotype == "Tim3.cKO"),
  design = "~ cell_type", outliers = outliers_batch1
)


## Get results table
res_phagoposvsneg_control  <- lfcShrink(dds_phagoposvsneg_control, coef="cell_type_phagopos_vs_phagoneg", type="apeglm")
res_phagoposvsneg_Tim3.cKO <- lfcShrink(dds_phagoposvsneg_Tim3.cKO, coef="cell_type_phagopos_vs_phagoneg", type="apeglm")
res_Tim3cKOvscontrol_phagopos <- lfcShrink(dds_Tim3cKOvscontrol_phagopos, coef="genotype_Tim3.cKO_vs_control", type="apeglm")
res_Tim3cKOvscontrol_phagoneg <- lfcShrink(dds_Tim3cKOvscontrol_phagoneg, coef="genotype_Tim3.cKO_vs_control", type="apeglm")

res_phagoposvsneg_control_unshrunken  <- results(dds_phagoposvsneg_control, contrast = c("cell_type", "phagopos", "phagoneg"))
res_phagoposvsneg_Tim3.cKO_unshrunken <- results(dds_phagoposvsneg_Tim3.cKO, contrast = c("cell_type", "phagopos", "phagoneg"))
res_Tim3cKOvscontrol_phagopos_unshrunken <- results(dds_Tim3cKOvscontrol_phagopos, contrast = c("genotype", "Tim3.cKO", "control"))
res_Tim3cKOvscontrol_phagoneg_unshrunken <- results(dds_Tim3cKOvscontrol_phagoneg, contrast = c("genotype", "Tim3.cKO", "control"))

## Order results
results_batch1_unshrunken <- list(
  phagoposvsneg_control  = res_phagoposvsneg_control_unshrunken,
  phagoposvsneg_Tim3.cKO = res_phagoposvsneg_Tim3.cKO_unshrunken,
  Tim3cKOvscontrol_phagopos = res_Tim3cKOvscontrol_phagopos_unshrunken,
  Tim3cKOvscontrol_phagoneg = res_Tim3cKOvscontrol_phagoneg_unshrunken) %>%
  lapply(function(res) {
    as.data.frame(res) %>%
      tibble::rownames_to_column("gene_id") %>%
      left_join(Tim3KO_cnts[, c("gene_id", "gene_name")], by = "gene_id") %>%
      select(gene_id, gene_name, everything()) %>%
      arrange(padj)  %>%
      mutate(direction = ifelse(log2FoldChange > 0, "up", "down")) %>%
      mutate(direction = factor(direction, levels = c("up", "down")))
  })

results_batch1 <- list(
  phagoposvsneg_control  = res_phagoposvsneg_control,
  phagoposvsneg_Tim3.cKO = res_phagoposvsneg_Tim3.cKO,
  Tim3cKOvscontrol_phagopos = res_Tim3cKOvscontrol_phagopos,
  Tim3cKOvscontrol_phagoneg = res_Tim3cKOvscontrol_phagoneg) %>%
  lapply(function(res) {
    as.data.frame(res) %>%
      tibble::rownames_to_column("gene_id") %>%
      left_join(Tim3KO_cnts[, c("gene_id", "gene_name")], by = "gene_id") %>%
      select(gene_id, gene_name, everything()) %>%
      arrange(padj)  %>%
      mutate(direction = ifelse(log2FoldChange > 0, "up", "down")) %>%
      mutate(direction = factor(direction, levels = c("up", "down")))
  })

save(results_batch1, results_batch1_unshrunken, file = "results/bulkRNAseq_results_ds1_batch1_3month.RData")
