library(Seurat)
library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(ggtext)
library(ComplexHeatmap)

# Signatures -------------------------------------
## P1/P2 signatures: top 100 DEGs of P1 vs P2 ordered by avg_log2FC) ----------------
load("results/DEG_bimod_TGFB_harmony_2_3_6_18.RData")

top100_DEGs <- DEG_bimod_TGFB_harmony_2_3_6_18 %>%
  group_by(direction) %>% arrange(desc(abs(avg_log2FC))) %>%
  dplyr::slice(1:100) %>% ungroup() %>% split(f = .$direction) %>%
  lapply(select, gene, avg_log2FC) %>%
  lapply(pull, "gene")

## Hallmark TGFB signature
Hallmark_TGFB <- read.table("data/signatures/HALLMARK_TGF_BETA_SIGNALING.txt", header = FALSE)$V1

# Load  microglia from public AD datasets & compute signature score ---------------------------
## (1) GSE140510 (7M snRNA-seq) ---------------------------
load("R_objects/GSE140510_ADmouse_7M_MG.RData")

## Hallmark TGFB signature
GSE140510_MG_obj <- AddModuleScore(GSE140510_MG_obj, list(Hallmark_TGFB), name = "Hallmark_TGFB")

## P1/P2 Signature score
GSE140510_MG_obj <- AddModuleScore(GSE140510_MG_obj, list(top100_DEGs$up), name = "P1_signature")
GSE140510_MG_obj <- AddModuleScore(GSE140510_MG_obj, list(top100_DEGs$down), name = "P2_signature")

## (2) GSE98969 (6M scRNA-seq) ---------------------------
load("R_objects/GSE98969_ADmouse_6M.RData")

## Hallmark TGFB signature
GSE98969_obj <- AddModuleScore(GSE98969_obj, list(Hallmark_TGFB), name = "Hallmark_TGFB")

## P1/P2 Signature score
GSE98969_obj <- AddModuleScore(GSE98969_obj, list(top100_DEGs$up), name = "P1_signature")
GSE98969_obj <- AddModuleScore(GSE98969_obj, list(top100_DEGs$down), name = "P2_signature")

# P1/P2 vs Hallmark TGFB signature score, split by genotype ------------------------------------

draw_scatter <- function(obj, sub = 1, suffix = "", title = "") {
  for(x in c("P1", "P2")) {
    
    ## metadata of cluster 2/DAM/MGnD
    meta <- subset(obj, seurat_clusters == sub)@meta.data %>%
      filter(Genotype == "5XFAD")
    
    ## Spearman correlation between P1/P2 and Hallmark TGFB score
    ttest <- cor.test(meta[[paste0(x, "_signature", 1)]], 
                      meta$"Hallmark_TGFB1", method = "spearman")
    
    cor_df <- data.frame(rho = ttest$estimate, p_val = ttest$p.value, Genotype = "5XFAD") %>%
      mutate(cor_label = paste0("r = ", round(rho, 3), "\np = ", round(p_val, 3)))
    
    meta %>%
      ggplot(aes_string(x = "Hallmark_TGFB1", y = paste0(y, "_signature1"))) +
      geom_point(shape = 16, size = 1) +
      geom_text(data = cor_df, aes(x = Inf, y = Inf, label = cor_label),
                hjust = 1, vjust = 1.2, size = 3.5, color = "red") +
      scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
      guides(color = guide_legend(override.aes = list(size = 1.6))) +
      facet_wrap(~ Genotype) +
      theme_cowplot() +
      labs(x = "Hallmark TGFB", y = paste(x, "Signature"), 
           title = title) +
      theme(plot.title = element_text(hjust = 0.5),
            strip.text.x = element_markdown())
    filename <- paste0("figures/Fig.S_", x, "_vs_Hallmark_TGFB_", suffix, "_sub", sub)
    ggsave(paste0(filename, ".png"),  width = 6, height = 4.5)
    ggsave(paste0(filename, ".pdf"),  width = 6, height = 4.5)
  }
}

## GSE98969, 6M, Cluster 2 (DAMs)
draw_scatter(GSE98969_obj, sub = 2, suffix = "GSE98969", 
             title = "GSE98969, 6M, Cluster 2 (DAMs)")
## GSE98969, 6M, Cluster 3 (DAMs)
draw_scatter(GSE98969_obj, sub = 3, suffix = "GSE98969", 
             title = "GSE98969, 6M, Cluster 3 (DAMs)")
## GSE140510, 7M, Cluster 1 (DAMs)
draw_scatter(GSE140510_MG_obj, sub = 1, suffix = "GSE140510", 
             title = "GSE140510, 7M, Cluster 1 (DAMs)")
