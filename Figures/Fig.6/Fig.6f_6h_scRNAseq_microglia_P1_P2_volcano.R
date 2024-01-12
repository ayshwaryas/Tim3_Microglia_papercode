library(tidyverse)
library(ggrepel)
library(DESeq2)
library(cowplot)
library(gridtext)
library(ggtext)
library(ggpubr)
library(rstatix)

## Palette
pal <- c("#CC79A7", "#009E73")


## Function for correcting gene symbols
correct_gene_symbol <- function(df, var = "gene_symbol") {
  df %>%
    mutate(gene_symbol = case_when(
      gene_symbol == "02-mars" ~ "Mars2",
      str_detect(gene_symbol, "-Mar") ~ paste0("March", gene_symbol),
      str_detect(gene_symbol, "-Sep") ~ paste0("Sept", gene_symbol),
      TRUE ~ gene_symbol),
      gene_symbol = str_remove(gene_symbol, "-Mar|-Sep"))
}

## Function for generating volcano plots
volcano_plot <- function(res, genes, width = 7, height = 7, step = 1,
                         position = c(0.25, 0.95), suffix = "", fontsize = 3, fig_num = "Fig.6F",
                         thres = 0.05, lfc_thres = 1, lim_x = NA, lim_y = NA, show_text = TRUE,
                         legend.box.margin = margin(0, 0, 0, -50),
                         genes_extra = "", geno = c("control", "<i>Havcr2</i>^cKO")) {
  if(length(lim_x) == 1) {lim_x <- c(-lim_x, lim_x)}
  highlight_genes <- res %>%
    filter(padj < thres) %>%
    filter(gene_symbol %in% unlist(genes))
  
  text_df <- highlight_genes %>%
    mutate(pt.color = ifelse(gene_symbol %in% genes$up, "P1", "P2")) %>% 
    mutate(hjust = ifelse(log2FoldChange > 0, 0, 1),
           nudge_x = ifelse(log2FoldChange > 0, 0.1, -0.1)) 
  
  if(floor(lfc_thres) != lfc_thres) {x_text <- paste0("log2(", 2^lfc_thres, ")")}
  else {x_text <- lfc_thres}
  
  p <- res %>%
    filter(!is.na(padj)) %>% 
    ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(size = 0.4, alpha = 0.8, shape = 16, color = "grey82") +
    geom_hline(yintercept = -log10(thres), size = 0.4, linetype = 2, color = "grey50") +
    geom_vline(xintercept = c(lfc_thres, -lfc_thres), size = 0.4, linetype = 2, color = "grey50") +
    geom_point(data = text_df, aes(color = pt.color, shape = pt.color),
               size = 1.6) +
    scale_x_continuous(limits = range(res$log2FoldChange)) +
    scale_y_continuous(sec.axis = dup_axis(breaks = -log10(thres), 
                                           label = paste0("padj=", thres), name = NULL),
                       expand = expansion(mult = c(0.02, 0.05))) +
    scale_color_manual(name = NULL, values = pal) +
    scale_shape_manual(name = NULL, values = c(16, 17)) +
    coord_cartesian(xlim = lim_x, ylim = c(-0.02 * lim_y, lim_y), expand = F) +
    labs(subtitle = geno[1], title = geno[2]) +
    theme_cowplot(font_size = 11) +
    theme(legend.position = position,
          legend.key = element_rect(color = "grey25"),
          legend.background = element_blank(),
          legend.title = element_text(size = 10.5, lineheight = 1.05, face = "bold"),
          legend.text = element_markdown(size = 10.5),
          axis.line = element_blank(),
          axis.text = element_text(size = 11),
          panel.border = element_rect(size = 0.7, color = "grey15"),
          plot.background = element_blank(),
          legend.box.margin = legend.box.margin,
          
          plot.title = element_markdown(face = "bold", hjust = 1, size = 12,
                                        margin = margin(5, 5, -11, 5)),
          plot.subtitle = element_markdown(face = "bold", size = 12)) 
  if(show_text) {
    p <- p + geom_text_repel(
      data = text_df,
      aes(label = gene_symbol, color = pt.color),
      size = fontsize, fontface = 3, 
      nudge_x =  text_df$nudge_x, show.legend = FALSE,
      max.overlaps = Inf, point.padding = 0.2,
      segment.size = 0.35, segment.alpha = 0.85, segment.color = "grey35",
      min.segment.length = 0.25, seed = 42, box.padding = 0.2
    )
  }
  filename <- paste0("figures/Fig.6/Fig.", fig_num, "_volcano_", suffix)
  ggsave(paste0(filename, ".png"), p, width = width, height = height, dpi = 500)
  ggsave(paste0(filename, ".pdf"), p, width = width, height = height)
}


# P1 & P2 signatures (top 100 DEGs of P1 vs P2 ordered by avg_log2FC) ----------------
load("results/DEG_bimod_TGFB_harmony_2_3_6_18.RData")
top100_DEGs <- DEG_bimod_TGFB_harmony_2_3_6_18 %>%
  group_by(direction) %>% arrange(desc(abs(avg_log2FC))) %>%
  dplyr::slice(1:100) %>% ungroup() %>% split(f = .$direction) %>%
  lapply(select, gene, avg_log2FC) %>%
  lapply(pull, "gene")

# Load Bulk RNA-seq differential gene expression analysis results --------------
## (1) Phagocytosing control vs non-phagocytosing control microglia from 3-month-old mice ------
load("data/DGE_results/2022-05-03.Dataset1_DGE_res_ordered.RData")
phago_all <- results_batch1_ordered$`phago+ vs phago- in control` %>%
  dplyr::rename("gene_symbol" = "gene_name") 

## (2) Clec7a+ vs Clec7a- ------
Clec7a_all <- read.csv("results/bulkRNAseq_results_Clec7a_Krasemann_2017.csv") %>%
  dplyr::rename("gene_symbol" = "tracking_id", "log2FoldChange" = "log2FC")

## (3) Tgfbr2cKO vs control ------
TGFBRII_all <- read.csv("results/bulkRNAseq_results_TGFBRII_Lund_2018.csv")  %>%
  select(geneNames, contains("uG")) %>%
  filter(rowSums(.[, 2:7]) != 0) %>%
  mutate(padj = p.adjust(.$p.ttest.uG, method = "BH")) %>%
  dplyr::rename("gene_symbol" = "geneNames",  "log2FoldChange" = "log2fc.uG") %>%
  correct_gene_symbol()


# Volcano plot --------------
## (1) Fig.6f: Phagocytosing control vs non-phagocytosing control microglia from 3-month-old mice -----
volcano_plot(phago_all, genes = top100_DEGs, 
             width = 4.5, height = 4, thres = 0.05, show_text = FALSE,
             position = "right", suffix = "nucseq_bimod_phago_pos_vs_neg",
             fig_num = "6F", legend.box.margin = margin(l = -55, r = 10),
             lfc_thres = Inf, lim_x = 3, lim_y = 10, 
             geno = c('Non-Phago MG', 'Phago MG'))

## (2) Fig.6g: Clec7a+ vs Clec7a- --------------------
volcano_plot(Clec7a_all, genes  = top100_DEGs,
             width = 4.5, height = 4, thres = 0.05, show_text = FALSE,
             legend.box.margin = margin(l = -55, r = 10),
             position = "right", suffix = "nucseq_bimod_Clec7a_pos_vs_neg", fig_num = "6G",
             lfc_thres = Inf, lim_x = 4, lim_y = 4, 
             geno = c('<i>Clec7a</i>^-', '<i>Clec7a</i>^+'))

## (3) Fig. 6h: Tgfbr2cKO vs WT --------------------
volcano_plot(TGFBRII_all, genes = top100_DEGs,
             width = 4.5, height = 4, thres = 0.05, show_text = FALSE,
             position = "right", suffix = "nucseq_bimod_Tgfbr2KO_vs_WT", fig_num = "6H",
             lfc_thres = Inf, lim_x = 8, lim_y = 3, 
             geno = c('control', '<i>Tgfbr2</i>^KO'),
             legend.box.margin = margin(l = -55, t = -50, r = 10))


# Spearman correlation between log2 FC of P1 vs P2 in snRNA-seq and bulk RNA-seq comparisons --------------------

cor_test_bulk_bimod <- function(bimod_res, bulk_res) {
  overlap <- inner_join(bulk_res, bimod_res, by = c("gene_symbol" = "gene"))
  
  corr <- cor.test(overlap$log2FoldChange, overlap$avg_log2FC, method = "spearman", exact = FALSE)
  data.frame(rho = corr$estimate, p = corr$p.value)
}

corr_res <- list(phago = phago_all,
                 Clec7a = Clec7a_all, 
                 Tgfbr2 = TGFBRII_all) %>%
  lapply(cor_test_bulk_bimod, bimod_res = DEG_bimod_TGFB_harmony_2_3_6_18) %>%
  data.table::rbindlist(idcol = 'bulk_ds') 
