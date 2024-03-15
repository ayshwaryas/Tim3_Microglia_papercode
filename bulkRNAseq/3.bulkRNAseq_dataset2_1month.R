library(tidyverse)
library(DESeq2)

# Load count matrix and metadata table ------------------------
count_1M_raw <- read.csv("data/expr_mat/Danyang_Count_matrix.csv", check.names = FALSE) %>%
  `colnames<-`(c("gene_symbol", colnames(.)[-1])) 

meta_1M <- read.csv("data/expr_mat/Danyang_Metadata.csv") %>%
  filter(str_detect(Sample_ID, "1M")) %>%
  mutate(sex = ifelse((Age == "1M" & Sample_number %in% c(5, 9)), "F", "M")) %>%
  mutate(Sample_ID_new = paste0(Sample_ID, "_", sex)) %>%
  column_to_rownames("Sample_ID") %>%
  mutate(genotype = factor(genotype, c("Tim3_flox", "Tim3_cKO")))

count_tbl_1M <- count_1M_raw[, c("gene_symbol", rownames(meta_1M))] %>%
  column_to_rownames("gene_symbol")

# DESeq2 -------------------------------------

## Construct DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = round(count_tbl_1M), colData = meta_1M,
                              design= ~ genotype)
dds <- DESeq(dds)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

## Get results table
res <- lfcShrink(dds, coef = "genotype_Tim3_cKO_vs_Tim3_flox", type="apeglm")

## Order results
res_ordered  <- res %>% as.data.frame() %>%
  rownames_to_column("gene_symbol") %>%
  select(gene_symbol, everything()) %>%
  mutate(direction = ifelse(log2FoldChange > 0, "up", "down")) %>%
  mutate(direction = factor(direction, c("up", "down"))) %>%
  arrange(direction, padj) 

save(res_ordered, file = "results/bulkRNAseq_results_ds2_1month.RData")

