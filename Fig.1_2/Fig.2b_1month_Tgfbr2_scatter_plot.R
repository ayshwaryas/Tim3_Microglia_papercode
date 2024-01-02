library(tidyverse)
library(data.table)
library(DESeq2)
library(edgeR)
library(ggrepel)
library(ggtext)
library(ComplexHeatmap)
library(circlize)
library(cowplot)

# Load data ---------
## 1 month, Havcr2cKO vs Havcr2flox/flox ----------
### Expression Matrices & Metadata 
count <- read.csv("data/expr_mat/Danyang_Count_matrix.csv", check.names = FALSE)
tpm <- read.csv("data/expr_mat/Danyang_TPM_matrix.csv", check.names = FALSE) 
metadata <- read.csv("data/expr_mat/Danyang_Metadata.csv") 

### Differential gene expression analysis results
load("data/DGE_results/2023-01-30.Danyang_DGE_results.RData")
tim3_all <- res_ordered %>%
  dplyr::select(gene_symbol, log2FoldChange, direction, padj)

## TGFBRII KO vs control RPKM -------------
TGFBRII <- read.csv("data/FPKM OB edits - with padj.csv") %>%
  mutate(log2fc.uG = - log2fc.uG)  %>%
  mutate(dir.uG = ifelse(dir.uG == "up", "down", "up")) %>%
  dplyr::rename("gene_symbol" = "geneNames")

keep_TGFBRII <- filterByExpr(TGFBRII[, 2:7], group = c(rep("WT", 3), rep("KO", 3)),
                             min.count = 2.5, min.prop = 0.8, min.total.count = 5)

TGFBRII_all <- TGFBRII[keep_TGFBRII, ] %>% 
  dplyr::select(gene_symbol, log2FoldChange = log2fc.uG,
                direction = dir.uG, padj = padj.uG) 

## Genes to highlight -----------
## Well knwon MGnD and Homeostatic microglia genes
gene_list <- c("Axl", "Ly9", "Clec7a", "Cd63", "Cxcl16", "H2-K1", 
               "H2-D1", "Siglech", "Csf1r")
### Immune checkpoint molecules and Tgfb pathway related genes
checkpoint_tgfb <- c("Havcr2", "Lag3", "Vsir", "Tgfb1", "Tgfbr1", "Tgfbr2", 
                         "Smad2", "Smad3", "Smad4")
### Top 100 Clec7a+ vs Clec7a- DEGs (up-regulated in Clec7a+: MGnD; down-regulated: Homeostatic)
Clec7a_sig <- read.csv("data/Oleg_Immunity_2017_DEG_ADpos_vs_neg_all.csv") %>%
  arrange(padj) %>% group_by(direction) %>% 
  dplyr::slice(1:100)

# Scatter plot ------------
## Function to calculate spearman correlation
lm_eqn <- function(df) {
  cor <- cor.test(-df$log2FoldChange.Tim3KO, -df$log2FoldChange.Tgfbr2KO, method = "spearman",exact=FALSE)
  R <- unname(cor$estimate)
  p <- formatC(cor$p.value, format = "e", digits = 2) %>% 
    str_replace("e-", " \u00D7 10<sup>-") %>% paste0("</sup>")
  paste0("R = ", round(R, 3), "<br/>p = ", p)
}

## Function to create signature plot
overlap_scatter <- function(RES, DDS, suffix, lim_x = NA, lim_y = 8,
                            add_line = TRUE, width = 8, height = 7, 
                            highlight_genes = NULL) {
  count_sub <- assay(DDS)
  keep <- filterByExpr(count_sub, group = colData(DDS)$genotype,
                       min.count = 5, min.prop = 0.8, min.total.count = 10)
  res_all <- RES %>% 
    filter(gene_symbol %in% rownames(count_sub[keep, ])) 
  
  res_sig_atleast300 <- RES %>%
    group_by(direction) %>%
    arrange(padj) %>% mutate(index = 1:n()) %>%
    dplyr::slice(1:max(sum(padj < 0.1, na.rm = TRUE), 300)) %>%
    filter(gene_symbol %in% rownames(count_sub[keep, ])) %>%
    dplyr::select(gene_symbol, log2FoldChange, padj, direction, index)
  
  TGFBRII_atleast300 <- TGFBRII_all %>%
    group_by(direction) %>%
    arrange(padj) %>% mutate(index = 1:n()) %>%
    dplyr::slice(1:max(sum(padj < 0.1, na.rm = TRUE), 300)) %>%
    filter(gene_symbol %in% TGFBRII_all$gene_symbol)
  
  overlap_all <- inner_join(res_all, TGFBRII_all,
                            by = "gene_symbol", suffix = c(".Tim3KO", ".Tgfbr2KO")) 
  
  
  overlap_atleast300 <- res_sig_atleast300 %>%
    inner_join(TGFBRII_atleast300, by = "gene_symbol", suffix = c(".Tim3KO", ".Tgfbr2KO")) %>%
    mutate(col = case_when(!gene_symbol %in% c(Clec7a_sig$tracking_id, gene_list, highlight_genes, checkpoint_tgfb) ~ "grey15",
                           log2FoldChange.Tim3KO > 0 & log2FoldChange.Tgfbr2KO > 0 ~ "red",
                           log2FoldChange.Tim3KO < 0 & log2FoldChange.Tgfbr2KO < 0 ~ "dodgerblue3",
                           TRUE ~ "black"))
  
  stat_df_all <- data.frame(x = -lim_x + 0.1,  y = -10 + 0.2, label = lm_eqn(overlap_all))
  
  label_df <- overlap_atleast300 %>%
    filter(col != "grey15") %>%
    mutate(nudge_x = ifelse(log2FoldChange.Tim3KO > 0 | abs(log2FoldChange.Tim3KO) >= lim_x, 0.2, -0.2),
           nudge_y = ifelse(log2FoldChange.Tgfbr2KO > 0, 0.2, -0.2))
  
  text_df <- overlap_all %>%
    filter(abs(log2FoldChange.Tim3KO) > lim_x | abs(log2FoldChange.Tgfbr2KO) >= lim_y |
             gene_symbol %in% checkpoint_tgfb) %>%
    filter(!gene_symbol %in% label_df$gene_symbol) %>%
    mutate(col = ifelse(gene_symbol %in% checkpoint_tgfb, NA, "grey40")) %>%
    mutate(nudge_x = case_when(
      (gene_symbol %in% checkpoint_tgfb) & (log2FoldChange.Tim3KO > 0) ~ 0.2,
      (gene_symbol %in% checkpoint_tgfb) & (log2FoldChange.Tim3KO <= 0) ~ -0.1,
      log2FoldChange.Tim3KO < -lim_x ~ 0.5,
      log2FoldChange.Tim3KO > lim_x ~ -0.4,
      TRUE ~ 0)) %>%
    mutate(nudge_y = case_when(
      (gene_symbol %in% checkpoint_tgfb) & (log2FoldChange.Tgfbr2KO > 0) ~ 0.1,
      (gene_symbol %in% checkpoint_tgfb) & (log2FoldChange.Tgfbr2KO <= -0.4) ~ -0.5,
      (gene_symbol %in% checkpoint_tgfb) & (log2FoldChange.Tgfbr2KO <= 0) ~ 0,
      log2FoldChange.Tgfbr2KO < 2 ~ 0.2,
      log2FoldChange.Tgfbr2KO >= 2 ~ -0.2)) %>%
    mutate(fontface = ifelse(gene_symbol %in% checkpoint_tgfb, 4, 3))
  
  p <- overlap_all %>%
    ggplot(aes(x = log2FoldChange.Tim3KO, y = log2FoldChange.Tgfbr2KO)) +
    geom_point(color = "grey80", size = 0.4, shape = 16) +
    geom_vline(xintercept = 0, size = 0.4, lty = 2, color = "grey30") +
    geom_hline(yintercept = 0, size = 0.4, lty = 2, color = "grey30") +
    geom_smooth(method = "lm", se = TRUE, formula = y ~ x, fullrange = TRUE,
                size = 0.5, color = "goldenrod3", fill = "goldenrod3", alpha = 0.2) +
    geom_point(data = overlap_atleast300, color = "grey15", size = 0.6) +
    geom_point(data = subset(overlap_atleast300, col != "grey15"), 
               aes(color = col), size = 1.5) +
    geom_point(data = subset(text_df, gene_symbol %in% checkpoint_tgfb), 
               aes(color = col), size = 0.6, shape = 16) +
    geom_label_repel(data = label_df, aes(label = gene_symbol, col = col), size = 4.5, 
                     max.overlaps = Inf, segment.size = 0.4, nudge_x = label_df$nudge_x, nudge_y = label_df$nudge_y,
                     segment.alpha = 0.7, label.padding = 0.2, min.segment.length = 0.2, seed = 1, 
                     box.padding = 0.2, show.legend = FALSE, label.size = 0.4, 
                     fontface = 4) +
    geom_text_repel(data = text_df, aes(label = gene_symbol, col = col), size = 4.5,
                    max.overlaps = Inf, segment.size = 0.4, box.padding = 0.2, 
                    segment.alpha = 0.7, min.segment.length = 0.1, seed = 1, 
                    show.legend = FALSE, fontface = text_df$fontface, 
                    nudge_x = text_df$nudge_x, nudge_y = text_df$nudge_y) +
    labs(x = "log<sub>2</sub>(<i>Havcr2</i><sup>cKO</sup>/control)", 
         y = "log<sub>2</sub>(<i>Tgfbr2</i><sup>cKO</sup>/control)") +
    scale_x_continuous(limits = range(overlap_all$log2FoldChange.Tim3KO, -lim_x, lim_x),
                       breaks = seq(-lim_x, lim_x, length.out = 5),
                       labels = seq(-lim_x, lim_x, length.out = 5)) +
    scale_color_identity()+ 
    coord_cartesian(xlim = c(-lim_x, lim_x), ylim = c(-10, 10), expand = FALSE)  +
    guides(fill = guide_legend(override.aes = list(color = "grey70")),
           size = guide_legend(override.aes = list(size = 0.8))) +
    geom_richtext(data = stat_df_all, aes(x = x, y = y, label = label),
                  hjust = 0, vjust = 0, size = 5, inherit.aes = FALSE) +
    theme_cowplot(font_size = 14) + 
    theme(legend.title = element_text(size = 14, face = "bold"),
          axis.text = element_markdown(size = 13),
          axis.title.x = element_markdown(hjust = 0.5, size = 14),
          axis.title.y = element_markdown(hjust = 0.5, size = 14)) 
  
  filename <- paste0("figures/Fig.2b_1month_scatter/Fig.2b_overlap_scatter_TGFBRIIKO_", suffix)
  ggsave(paste0(filename, ".png"), p, width = width, height = height, dpi = 400)
  ggsave(paste0(filename, ".pdf"), p, width = width, height = height)
  
  return(list(overlap_atleast300 = overlap_atleast300, p = p,
              capped_genes = subset(overlap_all, abs(log2FoldChange.Tim3KO) > lim_x)))
}

scatter_1M <- overlap_scatter(tim3_all, DDS = dds, suffix = "1M", 
                              lim_x = 2.5, width = 7, height = 7)

