library(tidyverse)

## Krasemann et al. 2017 (PMID: 28930663) Table S2
Krasemann_2017_TableS2 <- readxl::read_xlsx(
  "data/NIHMS901481-supplement-3.xlsx", 
  sheet = "Fig.2E RNAseq Clec7a+ in APPPS1", skip = 1) %>%
  select(-contains("WT"))

## Filter genes (total RPKM > 0.1)
Clec7apos_vs_neg <- Krasemann_2017_TableS2[rowSums(Krasemann_2017_TableS2[, -1]) > 0.1, ] %>%
  mutate(pval = NA, log2FC = NA)

for(i in 1:nrow(Clec7apos_vs_neg)) {
  Clec7apos <- as.numeric(select(Clec7apos_vs_neg, contains("pos"))[i, ])
  Clec7aneg <- as.numeric(select(Clec7apos_vs_neg, contains("neg"))[i, ])
  
  Clec7apos_vs_neg$pval[i]  <- t.test(Clec7apos, Clec7aneg)$p.value
  Clec7apos_vs_neg$log2FC[i]  <- log2(mean(Clec7apos) + 1) - log2(mean(Clec7aneg) + 1)
}

Clec7apos_vs_neg$padj  <- p.adjust(Clec7apos_vs_neg$pval, method = "BH")
Clec7apos_vs_neg$direction  <- ifelse(Clec7apos_vs_neg$log2FC > 0, "up", "down")

write.csv(Clec7apos_vs_neg, "results/bulkRNAseq_results_Clec7a_Krasemann_2017.csv", row.names = FALSE)


