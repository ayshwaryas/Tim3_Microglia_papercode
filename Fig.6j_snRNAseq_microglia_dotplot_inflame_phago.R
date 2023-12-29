library(Seurat)
library(tidyverse)
library(ggpubfigs)
library(ggtext)
library(cowplot)

## Seurat object of microglia nuclei
load("R_objects/2022-09-08.nucseq_harmony_MG_recluster_sub_2_3_6_18.RData")

## P1, P2 marker genes that are known to have inflammatory and phagocytic properties
Fig6j_genes <- c("Tlr4", "Il21r", "Mlxipl", "Stat1",                             # Pro-inflammatory
                 "Il10ra", "Nrp1", "Arpc1b",                                     # Anti-inflammatory
                 "Tmem59", "Gas7", "Mgat3", "Lgals3bp", "Tac1", "Lsp1", "Prkch", # Phagocytosis
                 "Zwint", "Ltc4s")                                               # Additional genes

## Dotplot of inflammatory and phagocytosis genes
bimod_genotype_breaks <- levels(nucseq_harmony_MG_2_3_6_18$bimod_genotype)[-4]
bimod_genotype_labels <- c("control", "<i>Havcr2</i><sup>icKO</sup>", "5XFAD", 
                           "<i>Havcr2</i><sup>icKO</sup> 5XFAD;P1", "<i>Havcr2</i><sup>icKO</sup> 5XFAD;P2")

DotPlot(subset(nucseq_harmony_MG_2_3_6_18, new_clusters == 2),
               features = Fig6j_genes, group.by = "bimod_genotype", scale = TRUE) +
  scale_y_discrete(breaks = bimod_genotype_breaks,
                   labels = bimod_genotype_labels,
                   limits = rev) +
  scale_color_gradientn(colors = rev(brewer.pal(11, "RdBu"))) +
  scale_radius(range = c(1.2, 12)) +
  guides(color = guide_colorbar(order = 1, barheight = 0.7, barwidth = 7, 
                                title = "Average Expression", title.vjust = 1), 
         size = guide_legend(title = "Percent Expression", nrow = 1, label.position = "bottom"), 
         radius = guide_legend(order = 2, title = "Percent Expression", nrow = 1, label.position = "bottom")) +
  labs(x = NULL, y = NULL) + 
  theme_cowplot(font_size = 14) +
  theme(axis.text.x = element_text(face = 3, size = 13, angle = 30, hjust = 1),
        axis.text.y = element_markdown(size = 13),
        legend.position = "bottom",
        legend.box.just = "bottom",
        legend.justification = c(0.5, 0)) 

ggsave("figures/Fig.6_nucseq_figures/Fig.6j_dotplot_inflammatory_phago_bimod_scaled.png", width = 9, height = 4)
ggsave("figures/Fig.6_nucseq_figures/Fig.6j_dotplot_inflammatory_phago_bimod_scaled.pdf", width = 9, height = 4)
