library(Seurat)
library(tidyverse)
library(cowplot)
library(ggrepel)
library(ggtext)

# Load Seurat object -----------------------------------------------------------
load("R_objects/2023-10-03.scRNAseq_cortex_MG_updated_clusters.RData")


# Differential gene expression analysis (Wilcoxon), IFN_DAM vs IFN_HMG ---------
scRNAseq_cortex_MG_IFN_DAMvsHMG <- FindMarkers(
  scRNAseq_cortex_MG,
  ident.1 = "IFN_DAM", ident.2 = "IFN_HMG", group.by = "updated_clusters", 
  logfc.threshold = 0, min.pct = 0.1, test.use = "wilcox") %>%
  rownames_to_column("gene") %>%
  filter(!str_detect(gene, "^Rp[ls]")) %>%             ## remove ribosomal genes
  mutate(FDR = p.adjust(.$p_val, method = "BH")) %>%
  mutate(direction = ifelse(avg_log2FC > 0, "up", "down")) %>%
  mutate(direction = factor(direction, c('up', 'down'))) 


# Late pre-MGnD vs Early pre-MGnD, Yin et al. 2023 (PMID: 37291336) Supplementary Table 2 ---------
## Female 
pre_MGnD_F <- readxl::read_xlsx("data/Yin_Herron_Supplementary_Table_2.xlsx",
                                sheet = "5 clusters_combined_f_Fig. 2e", skip = 3,
                                .name_repair = "minimal") %>% select(1:6) %>%
  filter(`Early pre-MGnD` != 0 | `Late pre-MGnD` != 0) %>%
  mutate(log2FoldChange = log2(`Late pre-MGnD`+ 1) - log2(`Early pre-MGnD` + 1))

## Male
pre_MGnD_M <- readxl::read_xlsx("data/Yin_Herron_Supplementary_Table_2.xlsx",
                                sheet = "5 clusters_combined_m_ED_Fig.6k", skip = 3,
                                .name_repair = "minimal") %>% select(1:6) %>%
  filter(`Early-pre-MGnD` != 0 | `Late-pre-MGnD` != 0) %>%
  mutate(log2FoldChange = log2(`Late-pre-MGnD`+ 1) - log2(`Early-pre-MGnD` + 1))


# Scatter plot of log2 FC between (1) late vs early pre-MGnD (Y-axis) and (2) IFN_DAM vs IFN_HMG (X-axis) ------------

## Joining the results table of comparison (1) and (2)
AllGenes_joined_F <- scRNAseq_cortex_MG_IFN_DAMvsHMG %>%  
  inner_join(pre_MGnD_F, by = c("gene" = "gene_id"))
AllGenes_joined_M <- scRNAseq_cortex_MG_IFN_DAMvsHMG_All[[x]] %>%  
  inner_join(pre_MGnD_M, by = c("gene" = "gene_id"))

AllGenes_joined <- list(`F` = AllGenes_joined_F, `M` = AllGenes_joined_M)

## Computing the Spearman correlation between log2 FC
corr_ls <- lapply(AllGenes_joined, function(df) {
  cor.test(df$log2FoldChange, df$avg_log2FC, method = "spearman")
})


## Scatter plot
## Genes to highlight
selected_DEGs <- c("Axl", "Gas6", "Apoe", "Lilrb4a", "Lpl", "Cst7", "Lyz2", "Npc2", 
                   "Cd9", "Cst3", "Fcrls", "P2ry12", "Tmem119", 
                   "Serinc3", "Siglech", "Hexb", "Mertk", "Ctsb", "Ctsz")

MGnD <- c("Axl", "Apoe", "Lpl", "Cst7", "Lyz2", "Ctsb", "Ctsz", "Cd9")
Homeostatic <- c("Fcrls", "Cst3", "P2ry12", "Tmem119", "Serinc3", "Siglech", "Hexb")

for(sex in c("F", "M")) {
  df <- AllGenes_joined[[sex]]
  corr <- corr_ls[[sex]]
  
  rho <- round(unname(corr$estimate), 2)
  p_val <- corr$p.value
  
  if(p_val < 10E-10) {p_val <- "< 10<sup>-10</sup>"}
  else {p_val <- paste0("= ", round(as.numeric(str_remove(p_val, "e.*$")), 3),
                        " \u00D7 10<sup>", str_remove(p_val, "^.*e"), "</sup>")}
  
  label_df <- subset(df, gene %in% selected_DEGs & FDR < 0.05) %>% 
    select(gene, avg_log2FC, log2FoldChange) %>%
    mutate(color = case_when(
      gene %in% Homeostatic ~ "dodgerblue3",
      gene %in% MGnD ~ "brown2", 
      TRUE ~ "grey15"))
  
  df %>% 
    ggplot(aes(x = log2FoldChange, y = avg_log2FC)) +
    geom_point(size = 0.4, color = "grey") +
    geom_vline(xintercept = 0, linewidth = 0.4, lty = 2, color = "grey30") +
    geom_hline(yintercept = 0, linewidth = 0.4, lty = 2, color = "grey30") +
    geom_point(data = subset(df, FDR < 0.05), size = 0.4, color = "grey15") +
    geom_point(data = label_df, aes(color = color), size = 1.5) +
    geom_label_repel(
      data = label_df,
      aes(label = gene, color = color), 
      max.overlaps = Inf, segment.size = 0.4, 
      segment.alpha = 0.7, label.padding = 0.2, min.segment.length = 0.2, seed = 1, 
      box.padding = 0.3, show.legend = FALSE, label.size = 0.4, 
      fontface = 4) +
    geom_richtext(data = data.frame(log2FoldChange = Inf, avg_log2FC = -Inf,
                                    label = paste0("r = ", rho, "<br/>p ", p_val)), 
                  aes(label = label), hjust = 0, vjust = 1) +
    labs(x = paste0("log<sub>2</sub>(late pre-MGnD/early pre-MGnD), ", sex), 
         y = paste0("log<sub>2</sub>(IFN_DAM/IFN_HMG)")) +
    scale_color_identity() +
    coord_flip() +
    theme_cowplot() +
    theme(axis.title.x = element_markdown(),
          axis.title.y = element_markdown())
  ggsave(paste0("figures/figure_panels/scatter_preMGnD_", sex, "_vs_IFN_DAMvsHMG_All.png"), width = 6, height = 5)
  ggsave(paste0("figures/figure_panels/scatter_preMGnD_", sex, "_vs_IFN_DAMvsHMG_All.pdf"), width = 6, height = 5)
  
  
}
