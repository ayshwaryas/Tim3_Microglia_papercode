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


# P1/P2 vs Hallmark TGFB signature score, split by genotype ------------------------------------
## Compute P1/P2 Signature score
nucseq_harmony_MG_2_3_6_18 <- AddModuleScore(
  nucseq_harmony_MG_2_3_6_18, list(top100_DEGs$up), name = "P1_signature")
nucseq_harmony_MG_2_3_6_18 <- AddModuleScore(
  nucseq_harmony_MG_2_3_6_18, list(top100_DEGs$down), name = "P2_signature")

for(x in c("P1", "P2")) {
  
  ## metadata of cluster 2/DAM/MGnD
  meta_sub2 <- subset(nucseq_harmony_MG_2_3_6_18, new_clusters == 2)@meta.data
  
  ## Spearman correlation between P1/P2 and Hallmark TGFB score
  cor_df <- meta_sub2 %>%
    filter(str_detect(Genotype, "5XFAD")) %>%
    split(.$Genotype_labels) %>%
    lapply(function(df) {
      test <- cor.test(df[[paste0(x, "_signature", 1)]], df$"Hallmark_TGFB1", 
                       method = "spearman")
      data.frame(rho = test$estimate, p_val = test$p.value)
    }) %>%
    data.table::rbindlist(idcol = "Genotype_labels") %>%
    mutate(Genotype_labels = factor(Genotype_labels, levels(nucseq_harmony_MG_2_3_6_18$Genotype_labels))) %>%
    mutate(cor_label = paste0("r = ", round(rho, 3), "\np = ", round(p_val, 3)))
  
  meta_sub2 %>%
    ggplot(aes_string(x = "Hallmark_TGFB1", y = paste0(y, "_signature1"))) +
    geom_point(aes(color = bimod_genotype), shape = 16, size = 1) +
    geom_text(data = cor_df, aes(x = Inf, y = Inf, label = cor_label),
              hjust = 1, vjust = 1.2, size = 3.5, color = "red") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    scale_color_manual(name = "Genotype", values = c("#785EF0", "#CC79A7", "#009E73"),
                       labels = c("5XFAD", paste0("<i>Havcr2</i><sup>icKO</sup> 5XFAD;", c("P1", "P2")))) +
    guides(color = guide_legend(override.aes = list(size = 1.6))) +
    facet_wrap(~ Genotype_labels) +
    theme_cowplot() +
    labs(x = "Hallmark TGFB", y = paste(x, "Signature"), 
         title = paste(x, "Signature vs Hallmark TGFB")) +
    theme(plot.title = element_text(hjust = 0.5),
          strip.text.x = element_markdown())
  filename <- paste0("figures/Fig.6_nucseq_figures/Fig.S_", x, "_vs_Hallmark_TGFB")
  ggsave(paste0(filename, ".png"),  width = 6.7, height = 5)
  ggsave(paste0(filename, ".pdf"),  width = 6.7, height = 5)
}

