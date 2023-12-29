library(Seurat)
library(tidyverse)
library(ggpubfigs)
library(ggtext)
library(cowplot)

## Seurat object of microglia nuclei
load("R_objects/2022-09-08.nucseq_harmony_MG_recluster_sub_2_3_6_18.RData")

## palette 
pal <- c("dodgerblue3", "#EF3B2C", friendly_pal("nickel_five", 5)[3],
         friendly_pal("ito_seven", 7)[c(6, 4:3)])

# DEGs comparing P1 and P2 -----------------------------------------------------
DAM_Tim3cKO.5XFAD <- subset(nucseq_harmony_MG_2_3_6_18, seurat_clusters == 2 & Genotype == "Tim3_cKO.5XFAD")

DEG_bimod_TGFB_harmony_2_3_6_18 <- FindMarkers(
  DAM_Tim3cKO.5XFAD, ident.1 = "Tim3_cKO.5XFAD_Top", ident.2 = "Tim3_cKO.5XFAD_Btm",
  min.pct = 0.1, logfc.threshold = 0, test.use = "wilcox", group.by = "bimod_genotype") %>%
  rownames_to_column("gene") %>%
  mutate(direction = ifelse(avg_log2FC > 0, 'up', 'down'))  

save(DEG_bimod_TGFB_harmony_2_3_6_18, file = "results/DEG_bimod_TGFB_harmony_2_3_6_18.RData")

# Get P1 and P2 signatures (top 100 DEGs ordered by avg_log2FC) ----------------
top100_DEGs <- DEG_bimod_TGFB_harmony_2_3_6_18 %>%
  group_by(direction) %>% arrange(desc(abs(avg_log2FC))) %>%
  dplyr::slice(1:100) %>% ungroup() %>% split(f = .$direction) %>%
  lapply(select, gene, avg_log2FC) %>%
  lapply(pull, "gene")


# Violin plot of P1 and P2 signature scores ------------------------------------
nucseq_harmony_MG_2_3_6_18 <- AddModuleScore(nucseq_harmony_MG_2_3_6_18, list(top100_DEGs$up), name = "P1_signature")
nucseq_harmony_MG_2_3_6_18 <- AddModuleScore(nucseq_harmony_MG_2_3_6_18, list(top100_DEGs$down), name = "P2_signature")

bimod_genotype_breaks <- levels(nucseq_harmony_MG_2_3_6_18$bimod_genotype)
bimod_genotype_labels <- c("control", "<i>Havcr2</i><sup>icKO</sup>", "5XFAD",
                           "<i>Havcr2</i><sup>icKO</sup> 5XFAD",
                           "<i>Havcr2</i><sup>icKO</sup> 5XFAD;P1",
                           "<i>Havcr2</i><sup>icKO</sup> 5XFAD;P2")

for(i in c("P1", "P2")) {
  VlnPlot(nucseq_harmony_MG_2_3_6_18, features = paste0(i, "_signature1"), 
          split.by = "bimod_genotype", group.by = "new_clusters") +
    labs(x = NULL, y = NULL, title = paste0(i, " signature")) +
    scale_fill_manual(breaks = bimod_genotype_breaks, 
                      labels = bimod_genotype_labels, 
                      values = pal) +
    theme_cowplot(font_size = 12) +
    theme(legend.text = element_markdown(size = 12),
          plot.title = element_text(hjust = 0.5))
}
