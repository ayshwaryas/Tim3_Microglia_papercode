library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(ggtext)
library(rstatix)
library(ggrepel)
library(ggpubr)
library(ggpubfigs)
library(Matrix)

## Seurat object of microglia nuclei
load("R_objects/2022-09-08.nucseq_harmony_MG_recluster_sub_2_3_6_18.RData")


## Subset 5XFAD and Havcr2icKO 5XFAD nuclei from cluster 2/DAM/MGnD
DAM_5XFAD <- subset(nucseq_harmony_MG_2_3_6_18, seurat_clusters == 2 & 
                 Genotype %in% c("5XFAD", "Tim3_cKO.5XFAD"))

groups <- DAM_5XFAD@meta.data[, c("orig.ident", "bimod_genotype")]
cnts_raw <- DAM_5XFAD@assays$RNA@counts

## Aggregate raw counts for pseudobulk analysis
pb <- aggregate(t(cnts_raw), groupings = groups, fun = "sum") %>% t()

## Normalize aggregated counts
pb_total <- apply(pb, 2, sum)
pb_norm <- sweep(pb, 2, pb_total, "/") * 1E6

## Exclude genes with 0 counts
pb_norm_nonzero <- pb_norm[apply(pb_norm, 1, sum) != 0, ]

## Compute Spearman correlation
corr_pb_log <- cor(log2(pb_norm_nonzero + 1), method = "spearman")
corr_pb_df <- corr_pb_log
corr_pb_df[lower.tri(corr_pb_df, diag = TRUE)] <- NA


# Fig. 6i: Boxplot of pairwise Spearman correlation of pseudobulk expression between conditions ----------

corr_pb_df_fig <- corr_pb_df %>% as.data.frame() %>%
  rownames_to_column("id1") %>%
  pivot_longer(2:ncol(.), values_to = "rho", names_to = "id2") %>%
  filter(!is.na(rho)) %>%
  mutate(Sample_ID1 = str_extract(id1, "Sample[0-9]+")) %>%
  mutate(Sample_ID2 = str_extract(id2, "Sample[0-9]+")) %>%
  mutate(Genotype1 = str_remove(id1, "Sample[0-9]+_")) %>%
  mutate(Genotype2 = str_remove(id2, "Sample[0-9]+_")) %>%
  filter(Genotype1  != Genotype2) %>%
  mutate(Contrast_id = paste(Sample_ID1, Sample_ID2, sep = " vs "))  %>%
  mutate(Contrast = paste(Genotype1, Genotype2, sep = " vs ")) %>%
  mutate(Contrast = case_when(Contrast == "5XFAD vs Tim3_cKO.5XFAD_Top" ~ "Tim3_cKO.5XFAD_Top vs 5XFAD",
                              Contrast == "5XFAD vs Tim3_cKO.5XFAD_Btm" ~ "Tim3_cKO.5XFAD_Btm vs 5XFAD",
                              Contrast == "Tim3_cKO.5XFAD_Btm vs Tim3_cKO.5XFAD_Top" ~ "Tim3_cKO.5XFAD_Top vs Tim3_cKO.5XFAD_Btm",
                              TRUE ~ Contrast)) %>%
  mutate(Contrast = str_replace(Contrast, " vs ", "\nvs ")) %>%
  mutate(Contrast = factor(Contrast, c("Tim3_cKO.5XFAD_Top\nvs 5XFAD",
                                       "Tim3_cKO.5XFAD_Btm\nvs 5XFAD",
                                       "Tim3_cKO.5XFAD_Top\nvs Tim3_cKO.5XFAD_Btm"))) %>%
  group_by(Contrast) %>% mutate(rho_med = median(rho)) %>% ungroup

## T-test
stat.test <- corr_pb_df_fig %>% t_test(rho ~ Contrast) %>%
  add_significance("p") %>%
  mutate(p.signif = ifelse(p >= 0.05, paste("p =", round(p, 3)), p.signif))  %>%
  mutate(size = ifelse(p < 0.05, 5, 3),
         vjust = ifelse(p < 0.05, 0.5, -0.1)) 

genotype_breaks <- c("Tim3_cKO.5XFAD_Top\nvs 5XFAD",
                     "Tim3_cKO.5XFAD_Btm\nvs 5XFAD",
                     "Tim3_cKO.5XFAD_Top\nvs Tim3_cKO.5XFAD_Btm")
genotype_labels <- c("<i>Havcr2<sup>icKO</sup></i> 5XFAD;P1<br/>vs 5XFAD",
                     "<i>Havcr2<sup>icKO</sup></i> 5XFAD;P2<br/>vs 5XFAD",
                     "<i>Havcr2<sup>icKO</sup></i> 5XFAD;P1<br/>vs <i>Havcr2<sup> icKO</sup></i> 5XFAD;P2")

## Boxplot
corr_pb_df_fig %>%
  ggplot(aes(x = Contrast, y = rho)) +
  geom_boxplot(aes(fill = rho_med), width = 0.4, alpha = 0.8, size = 0.6, color = "grey15") +
  geom_point(size = 2, shape = 16, color = "grey15", 
             position = position_jitter(width = 0.1, seed = 42))  +
  scale_fill_gradientn(colors = brewer.pal(9, "YlOrBr")[3:7]) +
  scale_x_discrete(breaks = genotype_breaks, labels = genotype_labels,
                   expand = expansion(add = c(0.5, 0.5))) +
  stat_pvalue_manual(stat.test, size = 5, fontface = 2,
                     y.position = c(0.645, 0.66, 0.63), 
                     bracket.size = 0.5, label = "p.signif", 
                     tip.length = 0.015, vjust = 0.5) +
  theme_cowplot(font_size = 12) +
  labs(x = NULL, y = "Spearman Correlation") +
  theme(axis.text.x = element_markdown(lineheight = 1.2, size = 11, angle = 15, halign = 0.5, 
                                       hjust = 0.8, vjust = 0.8, margin = margin(l = -5)),
        axis.text.y = element_text(size = 11),
        legend.position = "none", 
        plot.margin = margin(5, 10, 5, 5))

ggsave("figures/Fig.6_nucseq_figures/Fig.6i_pseudobulk_corr_stripchart.png", height = 5, width = 4.3)
ggsave("figures/Fig.6_nucseq_figures/Fig.6i_pseudobulk_corr_stripchart.pdf", height = 5, width = 4.3)
