library(tidyverse)
library(ggrepel)
library(pheatmap)
library(ComplexHeatmap)
library(DESeq2)
library(cowplot)
library(gridtext)
library(ggtext)

# Read DESeq2 results ------------
load("results/bulkRNAseq_results_ds2_1month.RData")

## Load TPM Matrix & Metadata -----------
tpm <- read.csv("data/expr_mat/Danyang_TPM_matrix.csv", check.names = FALSE)
metadata <- read.csv("data/expr_mat/Danyang_Metadata.csv") 

# Genes of Interest ----------
## MGnD & Homeostasis signature (Top 100 DEGs, Clec7a+ vs Clec7a-)
Clec7a_sig <- read.csv("results/bulkRNAseq_results_Clec7a_Krasemann_2017.csv") %>%
  arrange(padj) %>% group_by(direction) %>% 
  dplyr::slice(1:100) %>%
  pull(tracking_id)

## Tgfbr2 related genes
gene_list_Tgfbr2 <- c("Havcr2", "Axl", "Ly9", "Clec7a", "Cd63", "Cxcl16", 
                      "H2-K1", "H2-D1", "Siglech", "Csf1r")
## Phagocytosis related genes
gene_list_phago <- c("Havcr2", "Cd9", "Cst7", "Cstb", "Cxcl16", "Ly9", "Lyz2", 
                     subset(tpm$gene_symbol, str_detect(tpm$gene_symbol, "H2-")))
## 5XFAD related genes
gene_list_5XFAD <- c("Havcr2", "Axl", "Ccl6", "Cxcl16", "Cd9", "Cd81", "Lyz2",
                     subset(tpm$gene_symbol, str_detect(tpm$gene_symbol, "Cts")))
## KEGG phagosome pathway (mmu04145)
KEGG_phagosome <- read.table("data/signatures/KEGG_phagosome.txt", sep = "\t") %>%
  mutate(gene = str_extract(V2, ".*(?=;)")) %>%
  mutate(details = str_extract(V2, "(?<=; ).*")) %>% select(-V2) %>% pull(gene)

## All Genes of interest
volcano_gene_list <- c("Axl", "Ly9", "Clec7a", "Cd63", "Cxcl16", "H2-K1", 
                       "H2-D1", "Siglech", "Csf1r", 
                       Clec7a_sig, gene_list_phago, gene_list_5XFAD, gene_list_Tgfbr2,
                       "Cst7", "B2m", "Csf2ra", "Fcrls", "Itgax")

# Fig.1i Volcano plot of differential gene expression analysis results (DESeq2) in of one month old Havcr2cKO mice -----------

lim_x = 3   ## x-axis limit
lim_y = 10  ## y-axis limit

## top 10 DEGs and the genes cut-offed by the x-axis limits
top_DEGs <- res %>%
  filter(padj < 0.05) %>%
  filter(abs(log2FoldChange) > 1) %>%
  group_by(direction) %>%
  dplyr::slice(1:10) %>%
  rbind(subset(res_ordered, abs(log2FoldChange) > lim_x)) %>%
  distinct() 

## Genes to highlight
text_df <- res_ordered %>%
  filter(padj < 0.05 & gene_symbol %in% volcano_gene_list) %>%
  rbind(top_DEGs) %>%
  distinct() %>%
  mutate(pt.color = ifelse(log2FoldChange > 0, "red", "#194fe3")) %>%
  mutate(color = ifelse(log2FoldChange > 0, "#ed2e23", "#3b58a7")) %>%
  mutate(color = ifelse(log2FoldChange > lim_x, "grey15", color)) %>%
  mutate(hjust = case_when(-log10(padj) > 5 ~ 0.4, 
                           log2FoldChange > 0 ~ 0, 
                           log2FoldChange < 0 ~ 0.7),
         vjust = ifelse(-log10(padj) > lim_y, 1, 0.5),
         nudge_x = case_when(log2FoldChange > 0 ~ 0.1, TRUE ~ -0.1)) %>%
  mutate(nudge_y = ifelse(padj > 0.035 & -log10(padj) < 5, 0.25, 0))

p <- res_ordered %>%
  filter(!is.na(padj)) %>% 
  ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
  geom_hline(yintercept = -log10(0.05), size = 0.4, linetype = 2, color = "grey40") +
  geom_vline(xintercept = c(1, -1), size = 0.4, linetype = 2, color = "grey40") +
  geom_point(size = 0.4, alpha = 0.8, shape = 16, color = "grey65") +
  geom_point(data = text_df, aes(color = pt.color), size = 1.5, shape = 16) +
  geom_text_repel(data = text_df, aes(label = gene_symbol, color = color, 
                                      hjust = hjust, vjust = vjust),
                  size = 5, fontface = 3, 
                  nudge_x = text_df$nudge_x, nudge_y = text_df$nudge_y, max.overlaps = Inf, 
                  segment.size = 0.35, segment.alpha = 0.85, segment.color = "grey35",
                  min.segment.length = 0.25, seed = 42, box.padding = 0.15) +
  scale_color_identity() +
  scale_x_continuous(breaks = c(1, -1, seq(-round(lim_x), round(lim_x), step)),
                     labels = c(1, -1, seq(-round(lim_x), round(lim_x), step)),
                     limits = range(res$log2FoldChange)) +
  scale_y_continuous(sec.axis = dup_axis(breaks = -log10(0.05), 
                                         label = paste0("padj=", 0.05), name = NULL),
                     expand = expansion(mult = c(0.02, 0.05))) +
  coord_cartesian(xlim = c(-lim_x, lim_x), ylim = c(-0.02 * lim_y, lim_y), expand = F) +
  labs(subtitle = "control", title = "<i>Havcr2</i><sup>cKO</sup>",
       y = expression(-log[10](padj)), x = expression(log[2]~"Fold Change"))+
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(),
        legend.key = element_rect(color = "grey25"),
        legend.background = element_blank(),
        legend.box.margin = margin(0, 0, 0, -50),
        axis.text = element_text(size = 13),
        axis.text.y.right = element_text(size = 14),
        axis.title = element_text(size = 16),
        plot.title = element_markdown(face = "bold", color = "red", hjust = 1, size = 16,
                                      margin = margin(5, 5, -10, 5)),
        plot.subtitle = element_markdown(face = "bold", color = "#194fe3", size = 16)) 

filename <- "figures/Fig.1h_1l_1month/Fig.1i_1M_volcano_DEG_p.0.05"
ggsave(paste0(filename, ".png"), p, width = 7.5, height = 6.5, dpi = 500)
ggsave(paste0(filename, ".pdf"), p, width = 7.5, height = 6.5)


