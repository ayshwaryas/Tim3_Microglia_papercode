library(Seurat)
library(tidyverse)
library(RColorBrewer)

# Load Seurat object -----------------------------------------------------------
load("R_objects/2023-10-03.scRNAseq_cortex_MG_updated_clusters.RData")

# Clec7a+ signature (Top 100 DEGs of Clec7a+ vs Clec7a-) -----------------------
## Up in Clec7a+: MGnD signature
## Down in Clec7a+: Homeostasis signature

Clec7a_sig <- read.csv("results/bulkRNAseq_results_Clec7a_Krasemann_2017.csv") %>%
  arrange(padj) %>% group_by(direction) %>% 
  dplyr::slice(1:100) %>%
  split(f = .$direction) %>%
  lapply(pull, tracking_id)

# IFN signature: Top 20 genes of IFN-rich microglia ordered by fold change (Ellwanger et al. 2021, PMID: 33446504)
pnas_sd01 <- read.csv("data/pnas.2017742118.sd01.csv")

## Select top 20 marker genes for each microglia population
pnas_sig <- sapply(str_remove(colnames(pnas_sd01)[4 + 6 * 0:10], "_Specificity"), function(MG_type) {
  pnas_sd01 %>%
    select(Symbol, contains(MG_type)) %>%
    arrange(desc(!!sym(paste0(MG_type, "_Fold_change")))) %>%
    dplyr::slice(1:20) %>%
    pull(Symbol)
}, simplify = FALSE, USE.NAMES = TRUE) 


# Compute signature score ------------------------------------------------------
scRNAseq_cortex_MG <- AddModuleScore(scRNAseq_cortex_MG, list(Clec7a_sig$up), name = "MGnD")
scRNAseq_cortex_MG <- AddModuleScore(scRNAseq_cortex_MG, list(Clec7a_sig$down), name = "Homeostasis")

scRNAseq_cortex_MG <- AddModuleScore(scRNAseq_cortex_MG, list(pnas_sig$IFNR), name = "IFNR")

# Fig.s14a: Feature plot of signature score ------------------------------------
FeaturePlot(scRNAseq_cortex_MG, c("MGnD1", "Homeostasis1"),
            label = TRUE, label.size = 4, repel = TRUE, cols = rev(brewer.pal(11, "Spectral"))) 
ggsave("figures/cortex_updated_clusters/featureplot_MG_updated_clusters_clec7a.png", width = 14, height = 5.5)
ggsave("figures/cortex_updated_clusters/featureplot_MG_updated_clusters_clec7a.pdf", width = 14, height = 5.5)


FeaturePlot(scRNAseq_cortex_MG, c("IFNR1"),
            label = TRUE, label.size = 4, repel = TRUE, cols = rev(brewer.pal(11, "Spectral"))) 
ggsave("figures/cortex_updated_clusters/featureplot_MG_updated_clusters_ifnr.png", width = 7, height = 5.5)
ggsave("figures/cortex_updated_clusters/featureplot_MG_updated_clusters_ifnr.pdf", width = 7, height = 5.5)
