library(Seurat)
library(tidyverse)
library(cowplot)
library(ggrepel)
library(ggtext)

# Load Data -----------------------------------------------------------
## Differential gene expression analysis (Wilcoxon) results, Havcr2icKO;5xFAD vs 5xFAD
load("results/20231005_scRNAseq_cortex_MG_updated_clusters_DGE_wilcox_min.pct1_named.RData")

## Clec7a+ vs Clec7a-
Clec7a_all <- read.csv("results/bulkRNAseq_results_Clec7a_Krasemann_2017.csv") %>%
  dplyr::rename("gene_symbol" = "tracking_id", "log2FoldChange" = "log2FC")  %>%
  select(gene_symbol, log2FoldChange, padj)

# Scatter plot of log2 FC between (1) Clec7a- vs Clec7a+ (Y-axis) (2) Havcr2icKO;5xFAD vs 5xFAD (X-axis) ------------

## Genes to highlight 
selected_DEGs <- c("Axl", "Gas6", "Apoe", "Lilrb4a", "Lpl", "Cst7", "Lyz2", "Npc2", 
                   "Cd9", "Cst3", "Fcrls", "P2ry12", "Tmem119", 
                   "Serinc3", "Siglech", "Hexb", "Mertk", "Ctsb", "Ctsz")

MGnD <- c("Axl", "Apoe", "Lpl", "Cst7", "Lyz2", "Ctsb", "Ctsz", "Cd9")
Homeostatic <- c("Fcrls", "Cst3", "P2ry12", "Tmem119", "Serinc3", "Siglech", "Hexb")

## Fig.s15a: HMG_0 ------------------------------------
AllGenes_joined_HMG_0 <- scRNAseq_cortex_MG_DGE_wilcox_min.pct1$HMG_0 %>%  
  inner_join(Clec7a_all, by = c("gene" = "gene_symbol"))
DEG_joined_HMG_0 <- scRNAseq_cortex_MG_DGE_wilcox_min.pct1$HMG_0 %>%  
  filter(FDR < 0.05) %>%
  inner_join(subset(Clec7a_all, padj < 0.05), by = c("gene" = "gene_symbol"))  %>%
  mutate(nudge_x = ifelse(log2FoldChange < 0, 0.2, -0.2)) %>%
  mutate(col = case_when(gene %in% MGnD ~ "brown2",
                         gene %in% Homeostatic ~ "dodgerblue3",
                         TRUE ~ "grey15"))

corr_Clec7a_HMG_0 <- cor.test(
  AllGenes_joined_HMG_0$avg_log2FC, -AllGenes_joined_HMG_0$log2FoldChange, 
  method = "spearman")

rho_HMG_0 <- round(unname(corr_Clec7a_HMG_0$estimate), 2)
p_val_HMG_0 <- corr_Clec7a_HMG_0$p.value
p_val_HMG_0 <- paste0("= ", round(as.numeric(str_remove(p_val_HMG_0, "e.*$")), 2),
                      " \u00D7 10<sup>", str_remove(p_val_HMG_0, "^.*e") %>%
                        str_replace("-0", "-"), "</sup>")

AllGenes_joined_HMG_0 %>%
  ggplot(aes(x = -log2FoldChange, y = avg_log2FC)) +
  geom_point(size = 0.3, color = "grey80") +
  geom_vline(xintercept = 0, linewidth = 0.4, lty = 2, color = "grey30") +
  geom_hline(yintercept = 0, linewidth = 0.4, lty = 2, color = "grey30") +
  geom_point(data = DEG_joined_HMG_0, size = 0.8, color = "grey15") +
  geom_text_repel(
    data = subset(AllGenes_joined_HMG_0, avg_log2FC < -0.45 | gene %in% selected_DEGs) %>%
      mutate(col = ifelse(avg_log2FC >= -0.45, NA, "grey40")),
    aes(label = gene, color = col), 
    max.overlaps = Inf, segment.size = 0.4, 
    segment.alpha = 0.7, min.segment.length = 0.2, seed = 1, 
    box.padding = 0.2, show.legend = FALSE, 
    fontface = 3) +
  geom_label_repel(
    data = subset(DEG_joined_HMG_0, gene %in% c(selected_DEGs)),
    aes(label = gene, col = col),
    nudge_x = subset(DEG_joined_HMG_0, gene %in% c(selected_DEGs))$nudge_x,
    max.overlaps = Inf, segment.size = 0.4, size = 4,
    segment.alpha = 0.7, label.padding = 0.2, min.segment.length = 0.2, seed = 1, 
    box.padding = 0.2, show.legend = FALSE, label.size = 0.4, 
    fontface = 4) +
  geom_richtext(data = data.frame(log2FoldChange = -Inf, avg_log2FC = -Inf,
                                  label = paste0("r = ", rho_HMG_0, "<br/>p ", p_val_HMG_0)), 
                aes(label = label), hjust = 0, vjust = 1) +
  geom_point(data = subset(DEG_joined_HMG_0, gene %in% selected_DEGs), aes(color = col), size = 1.5) +
  labs(x = "log<sub>2</sub>(<i>Clec7a</i><sup>-</sup>/<i>Clec7a</i><sup>+</sup>)", 
       y = "log<sub>2</sub>(<i>Havcr2</i><sup>icKO</sup>5XFAD/5XFAD), HMG_0") +
  scale_y_continuous(limits = c(-3, 0.25)) +
  scale_color_identity() +
  coord_flip(ylim = c(-0.45, 0.25)) +
  theme_cowplot() +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_markdown())

ggsave("figures/cortex_updated_clusters/scatter_Clec7a_vs_Tim3icKO_HMG_0.png", width = 6, height = 5)
ggsave("figures/cortex_updated_clusters/scatter_Clec7a_vs_Tim3icKO_HMG_0.pdf", width = 6, height = 5)


## Fig.s15b: DAM_0 ------------------------------------
AllGenes_joined_DAM_0 <- scRNAseq_cortex_MG_DGE_wilcox_min.pct1$DAM_0 %>%  
  inner_join(Clec7a_all, by = c("gene" = "gene_symbol"))
DEG_joined_DAM_0 <- scRNAseq_cortex_MG_DGE_wilcox_min.pct1$DAM_0 %>%  
  filter(FDR < 0.05) %>%
  inner_join(subset(Clec7a_all, padj < 0.05), by = c("gene" = "gene_symbol"))  %>%
  mutate(nudge_x = ifelse(log2FoldChange > 0, 0.2, -0.2))  %>%
  mutate(col = case_when(gene %in% MGnD ~ "brown2",
                         gene %in% Homeostatic ~ "dodgerblue3",
                         TRUE ~ "grey15"))

corr_Clec7a_DAM_0 <- cor.test(
  AllGenes_joined_DAM_0$avg_log2FC, -AllGenes_joined_DAM_0$log2FoldChange, 
  method = "spearman")
rho_DAM_0 <- round(unname(corr_Clec7a_DAM_0$estimate), 2)
p_val_DAM_0 <- corr_Clec7a_DAM_0$p.value
p_val_DAM_0 <- paste0("= ", round(as.numeric(str_remove(p_val_DAM_0, "e.*$")), 2),
                      " \u00D7 10<sup>", str_remove(p_val_DAM_0, "^.*e") %>%
                        str_replace("-0", "-"), "</sup>")


AllGenes_joined_DAM_0 %>%
  ggplot(aes(x = -log2FoldChange, y = avg_log2FC)) +
  geom_point(size = 0.3, color = "grey80") +
  geom_vline(xintercept = 0, linewidth = 0.4, lty = 2, color = "grey30") +
  geom_hline(yintercept = 0, linewidth = 0.4, lty = 2, color = "grey30") +
  geom_point(data = DEG_joined_DAM_0, size = 0.8, color = "grey15") +
  geom_text_repel(
    data = subset(AllGenes_joined_DAM_0, avg_log2FC < -0.45 | gene %in% selected_DEGs) %>%
      mutate(col = ifelse(avg_log2FC >= -0.45, NA, "grey40")),
    aes(label = gene, color = col), 
    max.overlaps = Inf, segment.size = 0.4, 
    segment.alpha = 0.7, min.segment.length = 0.2, seed = 1, 
    box.padding = 0.2, show.legend = FALSE, 
    fontface = 3) +
  geom_label_repel(
    data = subset(DEG_joined_DAM_0, gene %in% c(selected_DEGs) ),
    aes(label = gene, color = col), 
    nudge_x = subset(DEG_joined_DAM_0, gene %in% c(selected_DEGs))$nudge_x,
    max.overlaps = Inf, segment.size = 0.4, 
    segment.alpha = 0.7, label.padding = 0.2, min.segment.length = 0.2, seed = 1, 
    box.padding = 0.2, show.legend = FALSE, label.size = 0.4, 
    fontface = 4, size = 4) +
  geom_point(data = subset(DEG_joined_DAM_0, gene %in% selected_DEGs), 
             aes(color = col), size = 1.5) +
  geom_richtext(data = data.frame(log2FoldChange = -Inf, avg_log2FC = -Inf,
                                  label = paste0("r = ", rho_DAM_0, "<br/>p ", p_val_DAM_0)), 
                aes(label = label), hjust = 0, vjust = 1) +
  labs(x = "log<sub>2</sub>(<i>Clec7a</i><sup>-</sup>/<i>Clec7a</i><sup>+</sup>)", 
       y = "log<sub>2</sub>(<i>Havcr2</i><sup>icKO</sup>5XFAD/5XFAD), DAM_0") +
  scale_y_continuous(limits = c(-3, 0.42)) +
  scale_color_identity() +
  coord_flip(ylim = c(-0.45, 0.42)) +
  theme_cowplot() +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_markdown())

ggsave("figures/cortex_updated_clusters/scatter_Clec7a_vs_Tim3icKO_DAM_0.png", width = 6, height = 5)
ggsave("figures/cortex_updated_clusters/scatter_Clec7a_vs_Tim3icKO_DAM_0.pdf", width = 6, height = 5)

