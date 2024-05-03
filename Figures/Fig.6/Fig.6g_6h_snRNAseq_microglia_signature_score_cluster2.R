library(Seurat)
library(tidyverse)
library(ggpubfigs)
library(ggtext)
library(cowplot)
library(rstatix)

## Seurat object of microglia nuclei
load("R_objects/2022-09-08.nucseq_harmony_MG_recluster_sub_2_3_6_18.RData")

## Signatures 
phagosome <- c("Gas7","Lgals3bp","Tac1","Lsp1","Mgat3","Prkch","Gstp1","Ltc4s","Zwint")
alternate_macrophage <- c("Cd163","Fn1","Mgl2","Retnla","Mrc1","Folr2","Selenop","Igf1","Cd36","Psap")

## Palette
bimod_genotype_breaks <- c("control", "Tim3_cKO", "5XFAD", "Tim3_cKO.5XFAD_Top", "Tim3_cKO.5XFAD_Btm")
bimod_genotype_labels <- c("control", "<i>Havcr2</i><sup>icKO</sup>", "5XFAD",
                           "<i>Havcr2</i><sup>icKO</sup> 5XFAD;P1", "<i>Havcr2</i><sup>icKO</sup> 5XFAD;P2")

pal <- c("dodgerblue3", "#EF3B2C", friendly_pal("nickel_five", 5)[3],
         friendly_pal("ito_seven", 7)[c(6, 4:3)])


# Fig 6g-h: signature score violin plot of alternate macrophage and phagocytosis signatures -------------------
## Alternate macrophage -------------------------------------------------------
nucseq_harmony_MG_2_3_6_18 <- AddModuleScore(
  nucseq_harmony_MG_2_3_6_18, list(phagosome), name = "Alternate_Macrophage")

score_df_alternate <- FetchData(subset(nucseq_harmony_MG_2_3_6_18, new_clusters == 2), 
                                c("bimod_genotype", "Alternate_Macrophage1")) %>%
  filter(str_detect(bimod_genotype, "5XFAD")) %>%
  mutate(bimod_genotype = factor(bimod_genotype, c("5XFAD", "Tim3_cKO.5XFAD_Btm", "Tim3_cKO.5XFAD_Top")))
  
## Anova
anova_res_alternate <- aov(Alternate_Macrophage1 ~ bimod_genotype, data = score_df_alternate)
  
## Tukey's HSD
TukeyHSD_res_alternate <- TukeyHSD(anova_res_alternate)$bimod_genotype %>%
  as.data.frame() %>%
  rownames_to_column("comparison")

stat.test_alternate <- TukeyHSD_res_alternate %>%
  separate("comparison", sep = "-", into = c("group1", "group2")) %>%
  rename("padj" = "p adj") %>%
  add_significance("padj") %>%
  mutate(padj.signif = ifelse(padj >= 0.05, paste("p =", round(padj, 3)), padj.signif)) %>%
  mutate(size = ifelse(padj < 0.05, 6, 2.5),
         vjust = ifelse(padj < 0.05, 0.5, -0.1)) 

write.csv(score_df_alternate, "/broad/kuchroolab/kimi_microglia/manuscript_code/Source_Data/Fig.6g_vln_alternate_sub2_score.csv")
write.csv(stat.test_alternate, "/broad/kuchroolab/kimi_microglia/manuscript_code/Source_Data/Fig.6g_vln_alternate_sub2_TukeyHSD.csv")

VlnPlot(subset(nucseq_harmony_MG_2_3_6_18, new_clusters == 2), "Alternate_Macrophage1",
        group.by = "bimod_genotype") +
  scale_fill_manual(values = pal[-4], breaks = bimod_genotype_breaks,
                    labels = bimod_genotype_labels, name = "Genotype") +
  stat_pvalue_manual(stat.test_alternate, size = 3.5, fontface = 2, inherit.aes = FALSE,
                     y.position = c(0.8, 0.7, 0.55) + 0.1, label.size = "size",
                     bracket.size = 0.5, label = "padj.signif", 
                     tip.length = 0.015, vjust = "vjust") +
  cowplot::theme_cowplot(font_size = 12) +
  scale_y_continuous(expand = expansion(add = 0.0, 0.06)) +
  labs(x = NULL, y = "Signature Score", title = "Alternate") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(),
        panel.spacing.x = unit(1.5, "mm"),
        plot.title = element_text(hjust = 0.5))

ggsave("figures/Fig.6g_vln_alternate_sub2_TukeyHSD.png", width = 5, height = 3.5)
ggsave("figures/Fig.6g_vln_alternate_sub2_TukeyHSD.pdf", width = 5, height = 3.5)

## Phagosome --------------------------------------------------------
nucseq_harmony_MG_2_3_6_18 <- AddModuleScore(
  nucseq_harmony_MG_2_3_6_18, list(phagosome), name = "Phagosome")

score_df_phagosome <- FetchData(subset(nucseq_harmony_MG_2_3_6_18, new_clusters == 2),
                                c("bimod_genotype", "Phagosome1")) %>%
  filter(str_detect(bimod_genotype, "5XFAD")) %>%
  mutate(bimod_genotype = factor(bimod_genotype, c("5XFAD", "Tim3_cKO.5XFAD_Btm", "Tim3_cKO.5XFAD_Top")))

## Anova
anova_res_phagosome <- aov(Phagosome1 ~ bimod_genotype, data = score_df_phagosome)

## Tukey's HSD
TukeyHSD_res_phagosome <- TukeyHSD(anova_res_phagosome)$bimod_genotype %>%
  as.data.frame() %>%
  rownames_to_column("comparison")

stat.test_phagosome <- TukeyHSD_res_phagosome %>%
  separate("comparison", sep = "-", into = c("group1", "group2")) %>%
  rename("padj" = "p adj") %>%
  add_significance("padj") %>%
  mutate(padj.signif = ifelse(padj >= 0.05, paste("p =", round(padj, 3)), padj.signif)) %>%
  mutate(size = ifelse(padj < 0.05, 6, 2.5),
         vjust = ifelse(padj < 0.05, 0.5, -0.1)) 

write.csv(score_df_phagosome, "/broad/kuchroolab/kimi_microglia/manuscript_code/Source_Data/Fig.6h_vln_phagosome_sub2_score.csv")
write.csv(stat.test_phagosome, "/broad/kuchroolab/kimi_microglia/manuscript_code/Source_Data/Fig.6h_vln_phagosome_sub2_TukeyHSD.csv")

VlnPlot(subset(nucseq_harmony_MG_2_3_6_18, new_clusters == 2), "Phagosome1",
        group.by = "bimod_genotype") +
  scale_fill_manual(values = pal[-4], breaks = bimod_genotype_breaks,
                    labels = bimod_genotype_labels, name = "Genotype") +
  stat_pvalue_manual(stat.test_phagosome, size = 3.5, fontface = 2, inherit.aes = FALSE,
                     y.position = c(0.88, 0.7, 0.8) + 0.1, label.size = "size",
                     bracket.size = 0.5, label = "padj.signif", 
                     tip.length = 0.015, vjust = "vjust") +
  cowplot::theme_cowplot(font_size = 12) +
  scale_y_continuous(expand = expansion(add = 0, 0.05)) +
  labs(x = NULL, y = "Signature Score", title = "Phagocytosis") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(),
        panel.spacing.x = unit(1.5, "mm"),
        plot.title = element_text(hjust = 0.5))

ggsave("figures/Fig.6h_vln_phagosome_sub2_TukeyHSD.png", width = 5, height = 3.5)
ggsave("figures/Fig.6h_vln_phagosome_sub2_TukeyHSD.pdf", width = 5, height = 3.5)


