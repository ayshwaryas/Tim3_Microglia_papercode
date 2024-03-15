library(Seurat)
library(tidyverse)
library(ggpubfigs)
library(ggtext)
library(cowplot)
library(RColorBrewer)
library(ggrepel)

## Seurat object of microglia nuclei
load("R_objects/2022-09-08.nucseq_harmony_MG_recluster_sub_2_3_6_18.RData")

## Function for drawing percentage barplot 
perc_bar <- function(SeuratObject, group_by = "new_clusters", fill_by = "orig.ident", pal = NULL, title = NULL,
                     width = 7.5, height = 3.5) {
  SeuratObject@meta.data %>%
    group_by_at(.vars = c(group_by, fill_by)) %>%
    summarise(n = n()) %>%
    ggplot(aes(x = .data[[group_by]], y = n, fill = .data[[fill_by]]))+ 
    geom_bar(position="fill", stat="identity", alpha = 0.8, width = 0.75) +
    scale_fill_manual(values = pal, name = title)+
    scale_y_continuous(expand = c(0, 0)) +
    labs(y = NULL, x = NULL) +
    theme_cowplot(font_size = 14) +
    theme(panel.grid = element_blank(),
          legend.text = element_markdown(),
          aspect.ratio = 2/3,
          axis.text.x = element_markdown(hjust = 1, angle = 45))
  
  filename <- paste("figures/Fig.6_nucseq_figures/Fig.s11_perc_bar", group_by, "by", fill_by, sep = "_")
  ggsave(paste0(filename, ".png"), width = width, height = height, dpi = 400)
  ggsave(paste0(filename, ".pdf"), width = width, height = height, dpi = 400)
}

## Create labels for mouse ID
mouse_id <- nucseq_harmony_MG_2_3_6_18@meta.data %>%
  select(Genotype_labels, orig.ident) %>%
  group_by(Genotype_labels, orig.ident) %>%
  dplyr::slice(1) %>%
  arrange(Genotype_labels, orig.ident) %>%
  group_by(Genotype_labels) %>%
  mutate(Mouse_ID = paste(Genotype_labels, "Mouse", row_number(), sep = "_")) %>%
  ungroup() %>%
  mutate(Mouse_ID = factor(Mouse_ID, levels = .$Mouse_ID))

meta <- nucseq_harmony_MG_2_3_6_18@meta.data %>%
  left_join(mouse_id, by = c("Genotype_labels", "orig.ident"))

nucseq_harmony_MG_2_3_6_18$Mouse_ID <- meta$Mouse_ID

# Extended Data Fig.13a: Percentage barplot by cluster, colored by mouse ID -------------------------------
perc_bar(nucseq_harmony_MG_2_3_6_18, fill_by = "Mouse_ID", 
         pal =  c("dodgerblue1", "dodgerblue4", "#FCAE91", "#DE2D26", 
                  brewer.pal(6, "Purples")[c(2, 4, 6)], "gold1", "#E69F00", "saddlebrown"),
         title = "Mouse_ID")

# Extended Data Fig.13b: Percentage barplot by cluster, colored by genotype -------------------------------
perc_bar(nucseq_harmony_MG_2_3_6_18, fill_by = "Genotype_labels", 
         pal =  c("dodgerblue3", "#EF3B2C", friendly_pal("nickel_five", 5)[3],
                 friendly_pal("ito_seven", 7)[6]),
         title = "Genotype", width = 6)

# Extended Data Fig.13c: Percentage barplot by Mouse_ID, colored by cluster -------------------------------
perc_bar(nucseq_harmony_MG_2_3_6_18, group_by = "Mouse_ID", fill_by = "new_clusters",
         pal =  c(friendly_pal("muted_nine", 3)[c(1, 3, 2)],
                  friendly_pal("nickel_five", 5), friendly_pal("bright_seven", 7)[5]),
         title = "Cluster", width = 5.6, height = 4.5)

