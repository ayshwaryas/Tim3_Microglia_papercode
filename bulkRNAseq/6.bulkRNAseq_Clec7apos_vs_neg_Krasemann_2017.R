library(tidyverse)

## Krasemann et al. 2017 (PMID: 28930663) Table S2
Krasemann_2017_TableS2 <- readxl::read_xlsx(
  "data/NIHMS901481-supplement-3.xlsx", 
  sheet = "Fig.2E RNAseq Clec7a+ in APPPS1", skip = 1) %>%
  select(-contains("WT"))

## Filter genes (total RPKM > 0.1)
ADpos_vs_neg <- Krasemann_2017_TableS2[rowSums(Krasemann_2017_TableS2[, -1]) > 0.1, ] %>%
  mutate(pval = NA, log2FC = NA)

for(i in 1:nrow(ADpos_vs_neg)) {
  ADpos <- as.numeric(select(ADpos_vs_neg, contains("pos"))[i, ])
  ADneg <- as.numeric(select(ADpos_vs_neg, contains("neg"))[i, ])
  
  ADpos_vs_neg$pval[i]  <- t.test(ADpos, ADneg)$p.value
  ADpos_vs_neg$log2FC[i]  <- log2(mean(ADpos) + 1) - log2(mean(ADneg) + 1)
}

ADpos_vs_neg$padj  <- p.adjust(ADpos_vs_neg$pval, method = "BH")
ADpos_vs_neg$direction  <- ifelse(ADpos_vs_neg$log2FC > 0, "up", "down")

write.csv(ADpos_vs_neg, "results/bulkRNAseq_results_Clec7a_Krasemann_2017.csv", row.names = FALSE)


