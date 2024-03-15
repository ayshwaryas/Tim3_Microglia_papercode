library(tidyverse)
library(Seurat)
library(RColorBrewer)

# Load Seurat object -----------------------------------------------------------
load("R_objects/2023-10-03.scRNAseq_cortex_MG_updated_clusters.RData")


# Fig.s14c: percentage barplot by genotype -------------------------------------
scRNAseq_cortex_MG@meta.data %>%
  group_by(updated_clusters, Genotype) %>%
  summarise(n = n()) %>%
  ggplot(aes_string(x = updated_clusters, y = "n", fill = Genotype)) + 
  geom_bar(position="fill", stat="identity", alpha = 0.8) +
  geom_hline(yintercept = 0.5, color = "red", lty = 2) +
  scale_fill_manual(values = c("#785EF0", "#E69F00"), drop = TRUE)+
  coord_cartesian(expand = FALSE) +
  labs(y = NULL) +
  theme_bw(base_size = 12) +
  guides(fill = guide_legend(keywidth = 0.7, keyheight = 0.7)) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(hjust = 1, angle = 45))
ggsave('figures/cortex_updated_clusters/perc_bar_updated_clusters_by_genotype.png', width = 7, height = 3)
ggsave('figures/cortex_updated_clusters/perc_bar_updated_clusters_by_genotype.pdf', width = 7, height = 3)
