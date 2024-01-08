library(Seurat)
library(tidyverse)
library(RColorBrewer)

## Seurat object of microglia nuclei
load("R_objects/2022-09-08.nucseq_harmony_MG_recluster_sub_2_3_6_18.RData")

## FindAllMarkers results on the microglia clusters
nucseq_harmony_MG_2_3_6_18 <- SetIdent(nucseq_harmony_MG_2_3_6_18, value = "new_clusters")
nucseq_harmony_MG_2_3_6_18_markers <- FindAllMarkers(
  nucseq_harmony_MG_2_3_6_18, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

save(nucseq_harmony_MG_2_3_6_18_markers, file = "results/2022-09-09.nucseq_harmony_MG_2_3_6_18_markers.RData")

## Dot plot of the top 5 markers ordered by FDR
top5 <- nucseq_harmony_MG_2_3_6_18_markers %>% group_by(cluster) %>% 
  arrange(p_val) %>% dplyr::slice(1:5) 

DotPlot(nucseq_harmony_MG_2_3_6_18, features = top5$gene, col.min = -2.2, col.max = 2.2,
        group.by = "new_clusters") +
  scale_color_gradientn(colors = rev(brewer.pal(11, 'RdYlBu'))) +
  scale_radius(range = c(0.1, 7)) +
  coord_flip() + scale_x_discrete(limits = rev) +
  labs(x = NULL, y = NULL)+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(face = 3, size = 12))
ggsave("figures/Fig.6_nucseq_figures/Fig.6b_dotplot_top5.png", height = 10, width = 8, dpi = 400)
ggsave("figures/Fig.6_nucseq_figures/Fig.6b_dotplot_top5.pdf", height = 10, width = 8)
