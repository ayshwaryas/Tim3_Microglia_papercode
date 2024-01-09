library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(ggtext)
library(parallel)
library(ComplexHeatmap)
library(ggpubr)
library(rstatix)
library(effsize)

# 1. Load data -----------------------------------
## Load Seurat object -------------------------
load("R_objects/2023-10-03.scRNAseq_cortex_MG_updated_clusters.RData")

## Markers of microglia populations in Ellwanger et al. 2021 (PMID: 33446504) -------
pnas_sd01 <- read.csv("data/pnas.2017742118.sd01.csv") ## Dataset S1B

## Select top 20 marker genes for each microglia population
pnas_sig <- sapply(str_remove(colnames(pnas_sd01)[4 + 6 * 0:10], "_Specificity"), function(MG_type) {
  pnas_sd01 %>%
    select(Symbol, contains(MG_type)) %>%
    arrange(desc(!!sym(paste0(MG_type, "_Fold_change")))) %>%
    dplyr::slice(1:20) %>%
    pull(Symbol)
}, simplify = FALSE, USE.NAMES = TRUE) 

# 2. Compute signature score --------------------------------------------------------
for(x in names(pnas_sig)) {
  scRNAseq_cortex_MG <- AddModuleScore(
    scRNAseq_cortex_MG, list(pnas_sig[[x]]), name = x)  }

# 3. T-test on signatures scores, Havcr2icKO;5xFAD vs 5xFAD in each cluster ---------
PNAS_ttest <- sapply(names(all_sigs)[13:23], function(x) {
  FetchData(scRNAseq_cortex_MG, c(paste0(x, 1), "updated_clusters", "Genotype")) %>%
    mutate(Genotype = factor(Genotype, c("Tim3icKO_5XFAD", "5XFAD"))) %>%
    split(f = .$updated_clusters) %>%
    lapply(function(df) {
      res <- t.test(as.formula(paste0(paste0(x, 1), "~", "Genotype")), data = df)
      data.frame(cluster = df$updated_clusters[1], 
                 statistic = res$statistic, p_value = res$p.value)}) %>%
    do.call(rbind, .) %>%
    mutate(padj = p.adjust(.$p_value, method = "BH"))
}, USE.NAMES = TRUE, simplify = FALSE)

## Data frame of T statistics in wide format
htmap_PNAS_tstat_df <- PNAS_ttest %>%
  data.table::rbindlist(idcol = "signature") %>%
  select(signature, cluster, statistic) %>%
  pivot_wider(names_from = cluster, values_from = statistic) %>%
  column_to_rownames("signature")

## Data frame of Direction of change (up/down-regulated in Havcr2icKO;5xFAD vs 5xFAD) in wide format
htmap_PNAS_dir_df <- PNAS_ttest %>%
  data.table::rbindlist(idcol = "signature") %>%
  mutate(direction = ifelse(statistic > 0, ">0 (up in Tim3icKO 5XFAD)", "<0")) %>%
  mutate(direction = factor(direction, c(">0 (up in Tim3icKO 5XFAD)", "<0"))) %>%
  select(signature, cluster, direction) %>%
  pivot_wider(names_from = cluster, values_from = direction) %>%
  column_to_rownames("signature")

## Data frame of FDR in wide format
htmap_PNAS_padj_df <- PNAS_ttest %>%
  data.table::rbindlist(idcol = "signature") %>%
  mutate(padj = ifelse(padj < 0.05, "grey15", "grey95")) %>%
  select(signature, cluster, padj) %>%
  pivot_wider(names_from = cluster, values_from = padj) %>%
  column_to_rownames("signature")

# 4. Fig. 7c: Heatmap of t-test results  ---------------------------------------
## T-test comparing signature scores of Havcr2icKO;5xFAD to 5xFAD in each cluster
## color: directions of change
## text: T statistics
## text color: black if FDR < 0.05; white if FDR >= 0.05

## Legend
lgd = Legend(labels = c("FDR<0.05", "FDR>=0.05"), title = NULL, type = "points", 
             pch = c("a", "a"), legend_gp = gpar(col = c("grey15", "grey95")), 
             size = unit(1, "mm"), background = "grey80")

p_PNAS_htmap <- Heatmap(
  htmap_PNAS_dir_df,
  col = list(`>0 (up in Tim3icKO 5XFAD)` = "#FBB4AE", `<0` = "lightsteelblue"),
  rect_gp = gpar(col = "white"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(ifelse(is.na(htmap_PNAS_tstat_df[i, j]), "",
                     formatC(htmap_PNAS_tstat_df[i, j], digits = 1, format = "f")), x, y, 
              gp = gpar(fontsize = 8, col = htmap_PNAS_padj_df[i, j],
                        fontface = ifelse(htmap_PNAS_padj_df[i, j] == "grey15", 2, 1))) },
  name = "t-statistic",
  column_names_rot = 45,
  cluster_rows = FALSE,
  cluster_column_slices = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10), 
  column_split = factor(str_remove(levels(PNAS_ttest[[1]]$cluster), "_([0-6]|HMG|DAM)$"),
                        c("HMG", "DAM", "CC_HMG", "CC_DAM", "IFN")),
  column_title_gp = gpar(fontsize = 12,fontface = 2)
)


png("figures/PNAS_ttest_htmap_RdBu.png", width = 9, height = 4, res = 400, units = "in")
draw(p_PNAS_htmap, padding = unit(c(0, 6, 1, 1), "mm"),
     annotation_legend_list = list(lgd), merge_legend = TRUE)
dev.off()

cairo_pdf("figures/PNAS_ttest_htmap_RdBu.pdf", width = 9, height = 4)
draw(p_PNAS_htmap, padding = unit(c(0, 6, 1, 1), "mm"),
     annotation_legend_list = list(lgd), merge_legend = TRUE)
dev.off()
