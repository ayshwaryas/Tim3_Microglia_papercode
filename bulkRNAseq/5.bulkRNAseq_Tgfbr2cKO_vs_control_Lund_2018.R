library(tidyverse)
library(edgeR)

TGFBRII_all <- read.csv("data/FPKM OB edits.csv") %>%
  filter(!is.na(WT.uG)) %>%
  select(gene_symbol = geneNames, contains("WT.uG"), contains("KO.uG")) 

TGFBRII <- TGFBRII_all[rowSums(TGFBRII_all[, -1]) != 0, ]

for(i in 1:nrow(TGFBRII)) {
  WT.uG <- as.numeric(TGFBRII[i, 2:4])
  KO.uG <- as.numeric(TGFBRII[i, 5:7])
  
  TGFBRII$p.ttest.uG[i]  <- t.test(WT.uG, KO.uG)$p.value
  TGFBRII$log2fc.uG[i]  <- log2(mean(KO.uG) + 1) - log2(mean(WT.uG) + 1)
}

TGFBRII$padj.uG  <- p.adjust(TGFBRII$p.ttest.uG, method = "BH")
TGFBRII$dir.uG  <- ifelse(TGFBRII$log2fc.uG > 0, "up", "down")


write.csv(TGFBRII, "results/bulkRNAseq_results_TGFBRII_Lund_2018.csv", row.names = FALSE)
