library(Seurat)
library(tidyverse)
library(ggtext)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(cowplot)

setwd("/broad/kuchroolab/kimi_microglia/nucseq_20220216")

## palette
c25 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00","black", "gold1", 
         "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", "gray70", "khaki2", 
         "maroon", "orchid1","deeppink1","blue1","steelblue4","darkturquoise", "green3", 
         "yellow4", "yellow2", "#FEC44F", "brown4","darkcyan", "darkgoldenrod2", 
         "darkolivegreen4", "cornsilk", "aquamarine1", "darkseagreen1", "honeydew3", "salmon", 
         "slateblue3", "thistle2", "slategray2", "springgreen", "rosybrown1", "plum2", 
         "chocolate", "blanchedalmond", "brown1", "aliceblue", "cadetblue1", "coral1", 
         "firebrick1", "darkslategrey", "darkslateblue", "#745745", "#a01b68", "#ffc372",
         "#fff0f0", "#ebd4d4", "#eeeeee", "#de4463", "#e8ffc1", "#19d3da", "#ee6f57", "#ff9642",
         "#646464", "#0072ce", "#dc9cbf", "#ED8E7C", "#F1ECC3", "#D9DD6B", "#7C83FD", "#FDD2BF",
         "#E98580", "#C6B4CE", "#FFDADA", "#B5EAEA", "#BB8760", "#C9E4C5", "#FAEBE0", "#F7DBF0",
         "#66DE93", "#FFC074", "#B2B8A3", "#F2F4C3", "#D1D9D9", "#A7D0CD", "#114E60", "#28B5B5")


## Genotype levels and labels to use in the figures
genotype_breaks <- c("control", "Tim3_cKO", "5XFAD", "Tim3_cKO.5XFAD")
genotype_labels <- c("control", "<i>Havcr2</i><sup>icKO</sup>", "5XFAD", "<i>Havcr2</i><sup>icKO</sup> 5XFAD")


# Add Azimuth annotations to the Seurat object ------------------------------

## Load Seurat object
load("results/doubletfinder/nucseq_harmony_PC20_res1_doubletFinder.RData")  ## All cells
load("results/doubletfinder/nucseq_harmony_subcluster_doubletFinder.RData") ## Microglia/PVM clusters

## Load Azimuth output 
azimuth_pred <- read.delim(paste0("results/azimuth/azimuth_pred_harmony.tsv"), row.names = 1) 

all(colnames(nucseq_harmony) == rownames(azimuth_pred))
nucseq_harmony <- AddMetaData(object = nucseq_harmony, metadata = azimuth_pred)

# Subset the Seurat object to remove doublets --------------------
## All cells
nucseq_harmony_rm_doublets <- subset(nucseq_harmony, DF.classifications_0.25_0.01_3433 == "Singlet")

## Microglia/PVM clusters
nucseq_harmony_subcluster_rm_doublets <- subset(nucseq_harmony_subcluster, DF.classifications_0.25_0.01_3433 == "Singlet")


# Finalize annotations ---------------
## Split L2/3 IT into 3 groups 
nucseq_harmony_rm_doublets_L23_IT_raw <- subset(nucseq_harmony_rm_doublets, predicted.subclass == "L2/3 IT")
nucseq_harmony_rm_doublets_L23_IT <- RunUMAP(nucseq_harmony_rm_doublets_L23_IT_raw,
                                             reduction = "harmony", dims = 1:20)
DimPlot(nucseq_harmony_rm_doublets_L23_IT, cols = c25)
ggsave("figures/nucseq_harmony_rm_doublets_L23_IT.png")

nucseq_harmony_rm_doublets_L23_IT_group <- FetchData(
  nucseq_harmony_rm_doublets_L23_IT, c("UMAP_1", "UMAP_2", "seurat_clusters")) %>%
  cbind(Cells = Cells(nucseq_harmony_rm_doublets_L23_IT)) %>%
  mutate(group = case_when(UMAP_1 < -9 ~ "L2/3 IT_2", UMAP_1 > 0 ~ "L2/3 IT_1", TRUE ~ "L2/3 IT_3"))

nucseq_harmony_rm_doublets_L23_IT$group <- nucseq_harmony_rm_doublets_L23_IT_group$group

DimPlot(nucseq_harmony_rm_doublets_L23_IT, group.by = "group", cols = c25)
ggsave("figures/nucseq_harmony_rm_doublets_L23_IT_group.png")


annot_levs <- c("Microglia", "Microglia-2", "BAM-PVM", 
                "Astro", "Endo", 
                "L2/3 IT_1", "L2/3 IT_2", "L2/3 IT_3", "L5 ET", "L5 IT", "L5/6 NP", "L6b", "L6 CT", "L6 IT", 
                "Oligo", "OPC", "Peri", "VLMC", "Cycling", "Meis2", 
                "Lamp5", "Pvalb", "Sncg", "Sst", "Sst Chodl", "Vip")
annot_new <- nucseq_harmony_rm_doublets@meta.data %>%
  rownames_to_column("Cells") %>%
  mutate(annot = case_when(
    seurat_clusters == 27 ~ "Cycling",
    Cells %in% Cells(subset(nucseq_harmony_subcluster_rm_doublets, seurat_clusters %in% c(2, 3, 6, 18))) ~ "Microglia",
    Cells %in% Cells(subset(nucseq_harmony_subcluster_rm_doublets, seurat_clusters == 12)) ~ "BAM-PVM",
    Cells %in% Cells(subset(nucseq_harmony_subcluster_rm_doublets, seurat_clusters == 15)) ~ "BAM",
    TRUE ~ predicted.subclass)) %>%
  left_join(nucseq_harmony_rm_doublets_L23_IT_group[, c("Cells", "group")], by = "Cells") %>%
  mutate(annot_new = ifelse(!is.na(group), group, annot)) %>%
  mutate(annot_new = ifelse(annot_new == "L6 IT Car3", "L6 IT", annot_new)) %>%
  mutate(annot_new = ifelse(annot_new == "BAM", "BAM-PVM", annot_new)) %>%
  mutate(annot_new = ifelse(annot_new == "Micro-PVM", "Microglia-2", annot_new)) %>%
  mutate(annot_new = factor(annot_new, annot_levs)) 

nucseq_harmony_rm_doublets$annot_new <- annot_new$annot_new

save(nucseq_harmony_rm_doublets, file = "R_objects/nucseq_harmony_rm_doublets_annotated.RData")