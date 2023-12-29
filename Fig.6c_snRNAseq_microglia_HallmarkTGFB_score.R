library(Seurat)
library(tidyverse)
library(ggpubfigs)
library(ggtext)
library(egg)
library(cowplot)

## Seurat object of microglia nuclei
load("R_objects/2022-09-08.nucseq_harmony_MG_recluster_sub_2_3_6_18.RData")


## Denisty of Hallmark TGFB score of 5XFAD and Havcr2cKO 5XFAD nuclei in cluster 2
TGFB_sig_score_df <- FetchData(subset(nucseq_harmony_MG_2_3_6_18, seurat_clusters == 2), 
                               c("Genotype", "bimod_genotype", "Hallmark_TGFB1")) %>%
  filter(str_detect(Genotype, "5XFAD")) %>%
  mutate(label = ifelse(Genotype == "5XFAD", "5XFAD", "<i>Havcr2</i><sup>icKO</sup> 5XFAD"))
 
p_dens <- TGFB_sig_score_df %>%
  ggplot(aes(x = Hallmark_TGFB1, color = label, fill = label)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = -0.01, lty = 2, size = 0.4, color = "grey15") +
  scale_y_continuous(expand = expansion(add = c(0.1, 0.8))) +
  scale_color_manual(values = c("#785EF0", "#E69F00")) +
  scale_fill_manual(values = c("#785EF0", "#E69F00")) +
  coord_cartesian(clip = "off") +
  theme_cowplot(font_size = 12) +
  labs(x = NULL, title = "Hallmark TGFB") +
  theme(plot.title = element_text(hjust = 0.5, margin = margin(0, 5, 3, 5)),
        legend.text = element_markdown())

p_rug <- TGFB_sig_score_df %>%
  ggplot(aes(x = Hallmark_TGFB1, color = label, fill = label)) +
  geom_rug(data = subset(TGFB_sig_score_df, Genotype == "5XFAD"),
           size = 0.3, length = unit(0.4, "npc"), sides = "tr")  +
  geom_rug(data = subset(TGFB_sig_score_df, Genotype == "Tim3_cKO.5XFAD"),
           size = 0.3, length = unit(0.4, "npc"))  +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_manual(values = c("#785EF0", "#E69F00")) +
  scale_fill_manual(values = c("#785EF0", "#E69F00")) +
  coord_cartesian(clip = "off") +
  theme_cowplot() +
  labs(x = NULL, y = NULL, title = NULL) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_markdown(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none",
        plot.margin = margin(29, 5, -39, 5))

g <- egg::ggarrange(p_rug, p_dens, ncol = 1, heights = c(0.11, 1))

ggsave(paste0("figures/Fig.6_nucseq_figures/Fig.6c_Hallmark_TGFB_density.png"), 
       g, width = 5, height = 2.3, dpi = 400)
ggsave(paste0("figures/Fig.6_nucseq_figures/Fig.6c_Hallmark_TGFB_density.pdf"), 
       g, width = 5, height = 2.3)
