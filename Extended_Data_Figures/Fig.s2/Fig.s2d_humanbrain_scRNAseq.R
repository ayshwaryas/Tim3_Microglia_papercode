library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(cowplot)

setwd("/broad/kuchroolab/kimi_microglia/manuscript_code")

## Load Seurat object and updated it from v3 to v4
load("data/human_brain/2019-08-08.brain_nuclei_res1.RObj")
human_brain <- UpdateSeuratObject(exp4)

## Load and map the cluster annotations
human_brain_clusters <- read.table("data/human_brain/clusters.txt", header = TRUE)

human_brain$celltype <- plyr::mapvalues(
  x = human_brain$res.1, from = human_brain_clusters$cluster,
  to = human_brain_clusters$celltype)
human_brain_genes <- rownames(human_brain)

save(human_brain, file = "data/human_brain/2022-10-25.brain_nuclei_res1.RData")

# Dotplot of checkpoint molecules and Tgfb pathway related genes
checkpoint_Tgfb <- c("HAVCR2", "LAG3", "C10orf54","PDCD1", "CD274", "CTLA4",
                     "TGFB1", "TGFBR1", "TGFBR2", 
                     "SMAD2", "SMAD3", "SMAD4")

DotPlot(human_brain, features = checkpoint_Tgfb, group.by = "celltype",
        scale = TRUE, dot.scale = 8) +
  coord_flip() +
  scale_x_discrete(limits = rev) +
  scale_color_gradientn(colours = rev(brewer.pal(11, "RdYlBu"))) +
  labs(y = NULL, x = NULL) +
  theme(axis.text.x = element_text(hjust = 1, angle = 45),
        axis.text.y = element_text(face = 3))

ggsave("figures/human_brain_co_inhibit_scaled.png", width = 9.5, height = 4.5, dpi = 400, units = "in")
ggsave("figures/human_brain_co_inhibit_scaled.pdf", width = 9.5, height = 4.5)
