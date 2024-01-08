library(Seurat)
library(tidyverse)
library(RColorBrewer)

# Load Seurat object -----------------------------------------------------------
load("R_objects/2023-10-03.scRNAseq_cortex_MG_updated_clusters.RData")

# Clec7a+ signature (Top 100 DEGs of Clec7a+ vs Clec7a-) -----------------------
## Up in Clec7a+: MGnD signature
## Down in Clec7a+: Homeostasis signature

Clec7a_sig <- read.csv("data/Oleg_Immunity_2017_DEG_ADpos_vs_neg_all.csv") %>%
  arrange(padj) %>% group_by(direction) %>% 
  dplyr::slice(1:100) %>%
  split(f = .$direction) %>%
  lapply(pull, tracking_id)

# Compute signature score ------------------------------------------------------
scRNAseq_cortex_MG <- AddModuleScore(scRNAseq_cortex_MG, list(Clec7a_sig$up), name = "MGnD")
scRNAseq_cortex_MG <- AddModuleScore(scRNAseq_cortex_MG, list(Clec7a_sig$down), name = "Homeostasis")


# Fig.s14a: Feature plot of signature score ------------------------------------
FeaturePlot(scRNAseq_cortex_MG, c("MGnD1", "Homeostasis1"),
            label = TRUE, label.size = 4, repel = TRUE, cols = rev(brewer.pal(11, "Spectral"))) 
ggsave("figures/cortex_updated_clusters/featureplot_MG_updated_clusters_clec7a.png", width = 14, height = 5.5)
ggsave("figures/cortex_updated_clusters/featureplot_MG_updated_clusters_clec7a.pdf", width = 14, height = 5.5)
