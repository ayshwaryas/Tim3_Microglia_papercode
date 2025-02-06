library(tidyverse)
library(rstatix)
library(ggtext)
library(ggpubr)
library(Seurat)

## Seurat object of microglia nuclei
load("R_objects/2022-09-08.nucseq_harmony_MG_recluster_sub_2_3_6_18.RData")

perc_df <- nucseq_harmony_MG_2_3_6_18@meta.data %>%
  filter(str_detect(Genotype, "5XFAD")) %>%
  group_by(new_clusters, orig.ident, Genotype) %>%
  summarise(n = n()) %>%
  select(new_clusters, orig.ident, Genotype, n) %>%
  group_by(orig.ident) %>%
  mutate(n_total = sum(n)) %>%
  mutate(perc = n/n_total) %>%
  group_by(new_clusters, Genotype) %>%
  mutate(perc_avg = mean(perc))

prop_ttest <- perc_df %>%
  filter(str_detect(Genotype, "5XFAD")) %>% 
  mutate(Genotype = factor(Genotype, c("Tim3_cKO.5XFAD", "5XFAD"))) %>%
  split(f = .$new_clusters) %>%
  lapply(function(df) {
    ttest <- t.test(perc ~ Genotype, data = df)
    data.frame(tstat = -ttest$statistic,
               pval = ttest$p.value) }) %>%
  data.table::rbindlist(idcol = "new_clusters") %>%
  filter(new_clusters %in% as.character(0:4)) %>%
  mutate(padj = p.adjust(.$pval, method = "bonferroni")) 

stat.test <- prop_ttest %>%
  add_significance("padj") %>%
  mutate(group1 = "5XFAD", group2 = "Tim3_cKO.5XFAD") %>%
  mutate(padj.signif = ifelse(padj >= 0.05, paste("p =", round(padj, 3)), padj.signif)) %>%
  mutate(size = ifelse(padj < 0.05, 5, 3),
         vjust = ifelse(padj < 0.05, 0.5, -0.1)) %>%
  left_join(perc_df %>% filter(str_detect(Genotype, "5XFAD")) %>%group_by(new_clusters) %>% 
              summarise(y.position = max(perc)+ 0.03), by = "new_clusters") 

perc_df %>%
  filter(str_detect(Genotype, "5XFAD")) %>% 
  filter(new_clusters %in% as.character(0:4)) %>%
  ggplot(aes(x = Genotype, y = perc, color = Genotype, shape = Genotype)) +
  geom_col(aes(y = perc_avg), position = "dodge", fill = NA, linewidth = 0.5) + 
  geom_point(position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.7, seed = 42)) +
  stat_pvalue_manual(stat.test, size = 3.5, fontface = 2, inherit.aes = FALSE,
                     bracket.size = 0.5, label = "padj.signif", 
                     tip.length = 0.015, vjust = -0.5) +
  scale_shape_discrete(labels = c("5XFAD", "<i>Havcr2</i><sup>icKO</sup>;5XFAD")) +
  scale_color_manual(values = c("#785EF0", "#E69F00"),
                     labels = c("5XFAD", "<i>Havcr2</i><sup>icKO</sup>;5XFAD")) +
  labs(x = NULL, y = "Proportion") +
  facet_grid(~new_clusters, switch = "x") +
  cowplot::theme_cowplot(font_size = 12) +
  scale_y_continuous(expand = expansion(add = c(0.01, 0.04))) +
  theme(legend.text = element_markdown(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank())

ggsave("figures/nucseq_proportion_ttest_bonferroni.png", width = 6, height = 3)
ggsave("figures/nucseq_proportion_ttest_bonferroni.pdf", width = 6, height = 3)