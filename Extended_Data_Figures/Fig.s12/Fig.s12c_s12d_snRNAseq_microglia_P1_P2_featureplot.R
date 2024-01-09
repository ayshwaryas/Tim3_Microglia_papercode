library(Seurat)
library(tidyverse)
library(ggpubfigs)
library(ggtext)
library(cowplot)
library(RColorBrewer)
library(ggrepel)

## Seurat object of microglia nuclei
load("R_objects/2022-09-08.nucseq_harmony_MG_recluster_sub_2_3_6_18.RData")

# P1 and P2 signatures (top 100 DEGs ordered by avg_log2FC) ----------------
load("results/DEG_bimod_TGFB_harmony_2_3_6_18.RData")

top100_DEGs <- DEG_bimod_TGFB_harmony_2_3_6_18 %>%
  group_by(direction) %>% arrange(desc(abs(avg_log2FC))) %>%
  dplyr::slice(1:100) %>% ungroup() %>% split(f = .$direction) %>%
  lapply(select, gene, avg_log2FC) %>%
  lapply(pull, "gene")


# Featureplot of P1 and P2 signature scores split by genotype ------------------------------------
nucseq_harmony_MG_2_3_6_18 <- AddModuleScore(nucseq_harmony_MG_2_3_6_18, list(top100_DEGs$up), name = "P1_signature")
nucseq_harmony_MG_2_3_6_18 <- AddModuleScore(nucseq_harmony_MG_2_3_6_18, list(top100_DEGs$down), name = "P2_signature")

for(x in c("P1_signature", "P2_signature")) {
  p_data <- FetchData(nucseq_harmony_MG_2_3_6_18, c(
    "UMAP_1", "UMAP_2", "new_clusters", "Genotype_labels", paste0(x, "1")))
  
  p_data_labels <- p_data %>%
    group_by(new_clusters, Genotype_labels) %>%
    summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
  
  ggplot(p_data, aes_string(x = "UMAP_1", y = "UMAP_2", color = paste0(x, 1))) +
    geom_point(shape = 16, size = 1) +
    geom_text_repel(data = p_data_labels, aes(label = new_clusters), 
                    size = 3.5, color = "black", seed = 42) +
    scale_color_gradientn(colors = rev(brewer.pal(11, "RdYlBu")), name = NULL) +
    facet_wrap(~ Genotype_labels) +
    theme_cowplot() +
    labs(title = str_to_title(str_replace(x, "_", " "))) +
    theme(strip.text = element_markdown(face = "bold"),
          strip.background = element_blank(),
          plot.title = element_text(hjust = 0.5),
          panel.spacing.y = unit(0, "mm"))
  filename <- paste0("figures/Fig.6_nucseq_figures/Fig.s12_featureplot_", x, "_by_genotype")
  ggsave(paste0(filename, ".png"), width = 6.5, height = 5.5)
  ggsave(paste0(filename, ".pdf"), width = 6.5, height = 5.5)
}