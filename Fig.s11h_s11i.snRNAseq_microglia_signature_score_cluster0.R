library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(ggtext)
library(rstatix)
library(ggrepel)
library(ggpubr)
library(ggpubfigs)

## Seurat object of microglia nuclei
load("R_objects/2022-09-08.nucseq_harmony_MG_recluster_sub_2_3_6_18.RData")

## palette 
pal <- c("dodgerblue3", "#EF3B2C", friendly_pal("nickel_five", 5)[3],
         friendly_pal("ito_seven", 7)[6])

## Function for drawing signature score violing plot --------------------
VlnPlot_custom <- function(
    SeuratData, feature_name, cluster = NULL, fig_num = '6',
    group.var = "new_clusters", split.by = "Genotype", title = "", 
    width = 10, height = 3, suffix = "", add_p_value = FALSE, 
    x_min = c(1, 1, 1.18), x_max = c(1.18, 1.36, 1.36), add1 = TRUE,
    y_pos = c(0.78, 0.9, 0.68), ylim = 0, hjust = 0.5, angle = 0,
    folder = "Fig.6_nucseq_figures", 
    levs = c("5XFAD", "<i>Havcr2</i><sup>icKO</sup> 5XFAD"))  {
  if(!is.null(cluster)) {
    SeuratData <- subset(SeuratData, new_clusters == cluster)
  }
  
  suffix <- paste0(feature_name, suffix)
  
  feature <- paste0(feature_name, ifelse(add1, "1", ""))
  p <- VlnPlot(SeuratData, features = feature, group.by = group.var, split.by = split.var) +
    labs(x = NULL, y = NULL, title = title) +
    theme_cowplot(font_size = 12) +
    scale_fill_manual(values = pal) +
    theme(legend.text = element_markdown(size = 12),
          plot.title = element_markdown(hjust = 0.5),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12, angle = angle, hjust= hjust))
  
  if(add_p_value) {
    p <- add_p_value(p, var = feature_name, y_pos = y_pos, ylim = ylim, 
                     levs = levs, x_min = x_min, x_max = x_max)
  }
  
  if(!is.null(cluster)) {
    suffix <- paste0(suffix, "_sub", cluster)  
  }
  
  filename <- paste0("figures/", folder, "/Fig.", fig_num, "_vln_", suffix)
  ggsave(paste0(filename, ".png"), p, width = width, height = height, dpi = 400)
  ggsave(paste0(filename, ".pdf"), p, width = width, height = height)
  
}

# Extended Data Fig.11h: signature score of resting signature in cluster 0 ----------------
## Resting signature: Top 20 genes of resting microglia ordered by fold change (Ellwanger et al. 2021, PMID: 33446504)

## Ellwanger et al. 2021 (PMID: 33446504) Dataset S1B 
pnas_sd01 <- read.csv("data/pnas.2017742118.sd01.csv")

## Select top 20 marker genes for each microglia population
pnas_sig <- sapply(str_remove(colnames(pnas_sd01)[4 + 6 * 0:10], "_Specificity"), function(MG_type) {
  pnas_sd01 %>%
    select(Symbol, contains(MG_type)) %>%
    arrange(desc(!!sym(paste0(MG_type, "_Fold_change")))) %>%
    dplyr::slice(1:20) %>%
    pull(Symbol)
}, simplify = FALSE, USE.NAMES = TRUE) 

nucseq_harmony_MG_2_3_6_18 <- AddModuleScore(
  nucseq_harmony_MG_2_3_6_18, list(pnas_sig$Resting), name = "Resting")

VlnPlot_custom(nucseq_harmony_MG_2_3_6_18, cluster = 0, 
               feature_name = "Resting", fig_num = "s11h",
               split.var = "Genotype", title = "Resting", add_p_value = TRUE,
               width = 5, height = 3.5, x_min = 1.115, x_max = 1.335, y_pos = 1)

# Extended Data Fig.11i: signature score of down in Tgfbr2cKO signature in cluster 0 ----------------
# Down in Tgfbr2cKO signature top 100 genes down-regulated in Tgfbr2cKO vs control with log2 FC > 2 ()
Top100DEG_DN_in_Tgfbr2KO <- read.csv("data/FPKM OB edits - with padj.csv") %>%
  select(geneNames, contains("uG")) %>%
  filter(log2fc.uG > 2) %>%
  arrange(padj.uG) %>%
  dplyr::slice(1:100) %>%
  pull(geneNames)

nucseq_harmony_MG_2_3_6_18 <- AddModuleScore(
  nucseq_harmony_MG_2_3_6_18, list(Top100DEG_DN_in_Tgfbr2KO), name = "DN_in_Tgfbr2KO")

VlnPlot_custom(nucseq_harmony_MG_2_3_6_18, cluster = 0, 
               feature_name = "DN_in_Tgfbr2KO", fig_num = "s11i", add_p_value = TRUE,
               split.var = "Genotype_labels", title = "Down in <i>Tgfbr2</i>^cKO",
               width = 5, height = 3.5, x_min = 1.115, x_max = 1.335, y_pos = 0.62)
