library(tidyverse)
library(Seurat)
library(RColorBrewer)

# Load Seurat object -----------------------------------------------------------
load("R_objects/2023-10-03.scRNAseq_cortex_MG_updated_clusters.RData")

# Get cluster markers using FindMarkers ----------------------------------------
scRNAseq_cortex_MG_markers <- FindAllMarkers(scRNAseq_cortex_MG, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
save(scRNAseq_cortex_MG_markers, file = "results/2023-10-04.scRNAseq_cortex_MG_markers_updated_clusters.RData")

# Fig.s14b: dotplot of top 10 markers orderd by log2FC --------------------------
top10_markers_log2FC <- scRNAseq_cortex_MG_markers %>%
  split(f = .$cluster) %>%
  lapply(function(df) {
    df %>% mutate(FDR = p.adjust(.$p_val, method = "BH")) }) %>%
  do.call(rbind, .) %>%
  filter(FDR < 0.05) %>% 
  arrange(cluster, desc(avg_log2FC)) %>%
  group_by(cluster) %>% 
  dplyr::slice(1:10) %>% ungroup

DotPlot(scRNAseq_cortex_MG, features = unique(top10_markers_log2FC$gene), 
        group.by = "updated_clusters", scale = TRUE) +
  scale_color_gradientn(colors = rev(brewer.pal(11, "RdBu"))) +
  scale_radius(range = c(0.1, 7)) +
  coord_flip() + 
  scale_x_discrete(limits = rev) +
  labs(x = NULL, y = NULL)+
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(face = 3, size = 12))

filename <- "figures/cortex_updated_clusters/top10_markers_updated_clusters_order_by_log2FC_scaled"
ggsave(paste0(filename, ".png"), height = 9, width = 20, dpi = 400)
ggsave(paste0(filename, ".pdf"), height = 9, width = 20, dpi = 400)
