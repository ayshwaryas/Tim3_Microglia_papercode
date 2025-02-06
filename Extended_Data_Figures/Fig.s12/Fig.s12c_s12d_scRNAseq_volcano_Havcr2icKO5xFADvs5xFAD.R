library(tidyverse)
library(cowplot)
library(ggrepel)
library(ggtext)

# Load Data --------------------------------------------------------------------
## Differential gene expression analysis (Wilcoxon) results, Havcr2icKO;5xFAD vs 5xFAD
load("results/20231005_scRNAseq_cortex_MG_updated_clusters_DGE_wilcox_min.pct1_named.RData")

## Genes to highlight
gene_list <- list(
  MGnD = data.frame(gene = c("Axl", "Apoe", "Lpl", "Cst7", "Lyz2", "Ctsb", "Ctsz", "Cd9")),
  Homeostatic = data.frame(gene = c("Fcrls", "Cst3", "P2ry12", "Tmem119", "Serinc3", "Siglech", "Hexb")),
  TAM = data.frame(gene = c("Mertk", "Tyro3", "Gas6")),
  MHC = data.frame(gene = c("H2-D1", "H2-K1", "H2-M3", "H2-Ab1", "H2-Q7")),
  P1_P2 = data.frame(gene = c("Btg2", "Sp100", "Ski", "Nrp1")),
  TNF = data.frame(gene = c("Nfkbia", "Jun", "Junb")),
  PI3K_Akt = data.frame(gene = c("Sgk1", "Hsp90b1", "Ywhah")),
  Autophagy = data.frame(gene = c("Ubb", "Ubc", "Uba52", "Scoc")),
  IFNA = data.frame(gene = c("Eif2ak2", "B2m", "Bst2", "Cd74", "Tent5a", "Lgals3bp", "Ifitm3", "Ifi44")),
  Other = data.frame(gene = c("Lilrb4a", "Egr1", "Hspa5", "Calr", "Ctsd"))) %>%
  data.table::rbindlist(idcol = "group")  %>%
  mutate(group_simplified = ifelse(group %in% c("MGnD", "Homeostatic", "MHC", "IFNA"), group, "Other")) %>%
  mutate(col = case_when(
    group == "MGnD" ~ "brown2",
    group == "Homeostatic" ~ "dodgerblue2", 
    group == "MHC" ~ "green",
    group == "IFNA" ~ "orange",
    TRUE ~ "grey15"
  )) %>%
  mutate(face = ifelse(group_simplified == "Other", 3, 4))


# Volcano plot of Havcr2icKO;5xFAD vs 5xFAD ------------------------------------
## (1) Fig.s12c: HMG_0 ---------------------------------------------------------

label_df_HMG_0 <- scRNAseq_cortex_MG_DGE_wilcox_min.pct1$HMG_0 %>%
  inner_join(gene_list, by = "gene") %>%
  filter(FDR < 0.05 | gene == "Havcr2")  %>%
  mutate(nudge_x = ifelse(avg_log2FC > 0, 0.02, -0.02),
         nudge_y = ifelse(gene == "H2-Ab1", 0.5, 0))
scRNAseq_cortex_MG_DGE_wilcox_min.pct1$HMG_0 %>%
  ggplot(aes(x = avg_log2FC, y = -log10(FDR))) +
  geom_point(color = "grey80", size = 0.8) +
  geom_hline(yintercept = -log10(0.05), lty = 2, color = "grey40") +
  geom_text_repel(data = label_df_HMG_0, aes(label = gene, color = group_simplified),
                  max.overlaps = Inf, size = 4.5, 
                  nudge_x = label_df_HMG_0$nudge_x, nudge_y = label_df_HMG_0$nudge_y,
                  fontface = label_df_HMG_0$face,
                  segment.size = 0.4, segment.alpha = 0.7, min.segment.length = 0.2, 
                  seed = 1, show.legend = FALSE) +
  geom_point(data = label_df_HMG_0, aes(color = group_simplified), size = 1.5) +
  scale_x_continuous(limits = c(-3, NA)) +
  scale_y_continuous(limits = c(0, Inf), expand = expansion(add = c(2, 0)),
                     sec.axis = dup_axis(breaks = -log10(0.05), labels = "FDR=0.05"), name = NULL) +
  scale_color_manual(values = c("brown2", "dodgerblue2", "orange", "#18AAAA", "grey15"),
                     breaks = c("MGnD", "Homeostatic", "MHC", "IFNA", "Other"),
                     name = NULL) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  coord_cartesian(xlim = c(-0.4, NA), ylim = c(0, 120)) +
  cowplot::theme_cowplot() +
  theme(legend.box.margin = margin(l = -70),
        axis.line.y.right = element_blank())
ggsave("figures/volcano_HMG_0.png", width = 6, height = 5.5)
ggsave("figures/volcano_HMG_0.pdf", width = 6, height = 5.5)

## (2) Fig.s12d: DAM_0 ---------------------------------------------------------
label_df_DAM_0 <- scRNAseq_cortex_MG_DGE_wilcox_min.pct1$DAM_0 %>%
  inner_join(gene_list, by = "gene") %>%
  filter(FDR < 0.05 | gene == "Havcr2")  %>%
  mutate(nudge_x = ifelse(avg_log2FC > 0, 0.02, -0.02)) %>%
  mutate(nudge_x = ifelse(gene == "Cst7", 0.03, nudge_x))

scRNAseq_cortex_MG_DGE_wilcox_min.pct1$DAM_0 %>%
  ggplot(aes(x = avg_log2FC, y = -log10(FDR))) +
  geom_point(color = "grey80", size = 0.8) +
  geom_hline(yintercept = -log10(0.05), lty = 2, color = "grey40") +
  geom_text_repel(data = label_df_DAM_0, aes(label = gene, color = group_simplified),
                  max.overlaps = Inf, size = 4.5, nudge_x = label_df_DAM_0$nudge_x,
                  fontface = label_df_DAM_0$face, box.padding = 0.3,
                  segment.size = 0.4, segment.alpha = 0.7, min.segment.length = 0.2,
                  seed = 1, show.legend = FALSE) +
  geom_point(data = label_df_DAM_0, aes(color = group_simplified), size = 1.5) +
  scale_x_continuous(limits = c(-3, NA)) +
  scale_y_continuous(limits = c(0, Inf), expand = expansion(add = c(1, 0)),
                     sec.axis = dup_axis(breaks = -log10(0.05), labels = "FDR=0.05"), name = NULL) +
  scale_color_manual(values = c("brown2", "dodgerblue2", "orange", "#18AAAA", "grey15"),
                     breaks = c("MGnD", "Homeostatic", "MHC", "IFNA", "Other"),
                     name = NULL) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  coord_cartesian(xlim = c(-0.4, 0.3), ylim = c(0, 35)) +
  cowplot::theme_cowplot() +
  theme(legend.box.margin = margin(l = -70),
        axis.line.y.right = element_blank())
ggsave("figures/volcano_DAM_0.png", width = 6, height = 5.5)
ggsave("figures/volcano_DAM_0.pdf", width = 6, height = 5.5)

