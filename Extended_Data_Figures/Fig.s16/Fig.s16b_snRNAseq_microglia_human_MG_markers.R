library(Seurat)
library(tidyverse)
library(harmony)
library(ggpubfigs)
library(ggpubr)
library(ggtext)
library(rstatix)

## Seurat object of microglia nuclei
load("R_objects/2022-09-08.nucseq_harmony_MG_recluster_sub_2_3_6_18.RData")

## Genotype labels and color palette
bimod_genotype_breaks <- c("control", "Tim3_cKO", "5XFAD", "Tim3_cKO.5XFAD_Top", "Tim3_cKO.5XFAD_Btm")
bimod_genotype_labels <- c("control", "<i>Havcr2</i><sup>icKO</sup>", "5XFAD",
                           "<i>Havcr2</i><sup>icKO</sup> 5XFAD;P1", "<i>Havcr2</i><sup>icKO</sup> 5XFAD;P2")

pal <- c("dodgerblue3", "#EF3B2C", friendly_pal("nickel_five", 5)[3],
         friendly_pal("ito_seven", 7)[c(6, 4:3)])

## Function to convert human gene symbols to mouse gene symbols
mouse_human_genes_raw <- read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
convert_human_to_mouse <- function(gene_list){
  mouse_human_genes_raw %>%
    filter(DB.Class.Key %in% subset(mouse_human_genes_raw, Symbol %in% gene_list )$DB.Class.Key) %>% 
    filter(Common.Organism.Name == "mouse, laboratory") %>%
    mutate(Symbol = ifelse(is.na(Symbol), str_to_title(Symbol), Symbol)) %>%
    pull(Symbol) %>% unique
}


## Sun et al 2023 Top 25 Markers
Sun_2023_StateMarkers <- readxl::read_xlsx("data/Sun_et_al_2023_TableS1_HumanAD.xlsx", sheet = "Page 2.StateMarkers")
Sun_2023_StateMarkers_top25_ls <- Sun_2023_StateMarkers %>%
  split(f = .$microgliaState) %>%
  lapply(function(df) {
    df %>% mutate(FDR = p.adjust(.$p_val, method = "BH")) %>%
      filter(FDR < 0.05) %>%
      arrange(desc(avg_log2FC)) %>%
      slice(1:25) %>%
      pull(gene) %>% convert_human_to_mouse
  })


for(x in names(Sun_2023_StateMarkers_top25_ls)) {
  nucseq_harmony_MG_2_3_6_18 <- AddModuleScore(
    nucseq_harmony_MG_2_3_6_18, list(Sun_2023_StateMarkers_top25_ls[[x]]), name = x, assay = "RNA")  }

## Anova test with Tukey's HSD correction on MG8 score in 5XFAD, P1, P2 in cluster 2 (DAM)
Sun_2023_Top25Markers_score_sub2 <- FetchData(subset(nucseq_harmony_MG_2_3_6_18, new_clusters == 2),
                                              c("bimod_genotype", paste0(names(Sun_2023_StateMarkers_top25_ls), 1))) %>%
  pivot_longer(2:ncol(.), names_to = "Signature", values_to = "Score") %>%
  mutate(Signature = str_remove(Signature, "1$")) %>%
  mutate(Signature = factor(Signature, paste0("MG", 0:12))) 

## Bartlett Test to examine whether 3 groups have different variance
## MG8: p = 0.1325 ==> fail to reject H0 (equal variance)
Sun_2023_Top25Markers_score_sub2 %>%
  filter(str_detect(bimod_genotype, "5XFAD")) %>%
  mutate(bimod_genotype = factor(bimod_genotype, c("5XFAD", "Tim3_cKO.5XFAD_Btm", "Tim3_cKO.5XFAD_Top"))) %>%
  mutate(Signature = as.character(Signature)) %>%
  split(f = .$Signature) %>%
  lapply(function(df) {
    test_res <- bartlett.test(Score ~ bimod_genotype, data = df) 
    data.frame(Bartlett_K = test_res$statistic, p_val = test_res$p.value)
  }) %>%
  data.table::rbindlist(idcol = "Signature")


TukeyHSD_res <- lapply(unique(Sun_2023_Top25Markers_score_sub2$Signature), function(x) {
  score_df <- Sun_2023_Top25Markers_score_sub2 %>%
    filter(Signature == x & str_detect(bimod_genotype, "5XFAD")) %>%
    mutate(bimod_genotype = factor(bimod_genotype, c("5XFAD", "Tim3_cKO.5XFAD_Btm", "Tim3_cKO.5XFAD_Top")))
  
  ## Anova
  anova_res <- aov(Score ~ bimod_genotype, data = score_df)
  
  ## Tukey's HSD
  TukeyHSD_res <- TukeyHSD(anova_res)$bimod_genotype %>%
    as.data.frame() %>%
    rownames_to_column("comparison") %>%
    cbind(Signature = x, .)
})%>%
  do.call(rbind, .) %>%
  mutate(Signature = factor(Signature, paste0("MG", 0:12))) %>%
  arrange(Signature)

write.csv(TukeyHSD_res, "results/nucseq_Sun_2023_MG0to12_score_sub2_TukeyHSD.csv", row.names = FALSE)

stat.test <- TukeyHSD_res %>% filter(Signature == "MG8") %>%
  separate("comparison", sep = "-", into = c("group1", "group2")) %>%
  rename("padj" = "p adj") %>%
  add_significance("padj") %>%
  mutate(padj.signif = ifelse(padj >= 0.05, paste("p =", round(padj, 3)), padj.signif)) %>%
  mutate(size = ifelse(padj < 0.05, 6, 2.5),
         vjust = ifelse(padj < 0.05, 0.5, -0.1)) %>%
  mutate(Signature = factor("MG8", paste0("MG", 0:12)))

## Violin plot of MG0-MG12 signature scores in cluster 2, split by genotype 
Sun_2023_Top25Markers_score_sub2 %>%
  ggplot(aes(x = bimod_genotype, y = Score, fill = bimod_genotype)) +
  geom_violin(scale = "width", linewidth = 0.4) +
  geom_point(size = 0.2, position = position_jitter(width = 0.3)) +
  scale_fill_manual(values = pal[-4], breaks = bimod_genotype_breaks,
                    labels = bimod_genotype_labels, name = "Genotype") +
  facet_grid(~ Signature) +
  stat_pvalue_manual(stat.test, size = 3.5, fontface = 2, inherit.aes = FALSE,
                     y.position = c(0.38, 0.55, 0.28) + 0.1, label.size = "size",
                     bracket.size = 0.5, label = "padj.signif", 
                     tip.length = 0.015, vjust = "vjust") +
  cowplot::theme_cowplot(font_size = 12) +
  labs(x = NULL, y = "Signature Score") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown(),
        panel.spacing.x = unit(1.5, "mm"))

ggsave('figures/Fig.s16b_nucseq_Sun_2023_MG0to12_Top25Markers_sub2.png', width = 10, height = 3)
ggsave('figures/Fig.s16b_nucseq_Sun_2023_MG0to12_Top25Markers_sub2.pdf', width = 10, height = 3)

