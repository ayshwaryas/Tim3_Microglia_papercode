library(tidyverse)
library(ggrepel)
library(pheatmap)
library(ComplexHeatmap)
library(DESeq2)
library(cowplot)
library(gridtext)
library(ggtext)
library(ggpubr)
library(rstatix)

# Read DESeq2 results ------------
load("results/bulkRNAseq_results_ds2_1month.RData")

## Load TPM Matrix & Metadata -----------
tpm <- read.csv("data/expr_mat/Danyang_TPM_matrix.csv", check.names = FALSE)
metadata <- read.csv("data/expr_mat/Danyang_Metadata.csv") 

## Signatures ----------
### MGnD & Homeostasic microglia signature (Top 100 DEGs, Clec7a+ vs Clec7a-)
Clec7a_sig <- read.csv("results/bulkRNAseq_results_Clec7a_Krasemann_2017.csv") %>%
  arrange(padj) %>% group_by(direction) %>% 
  dplyr::slice(1:100) %>%
  select(tracking_id, direction)

## KEGG phagosome pathway (mmu04145)
KEGG_phagosome <- read.table("data/signatures/KEGG_phagosome.txt", sep = "\t") %>%
  mutate(gene = str_extract(V2, ".*(?=;)")) %>%
  mutate(details = str_extract(V2, "(?<=; ).*")) %>% select(-V2) %>% pull(gene)


# Fig.1j-1l: Signature score boxplot of 1 month old mice, Havcr2cKO vs control ------------
## Exclude genes with average TPM < 1 across replicates
gene_set_filter <- function(gene_list) {
  tpm %>%
    filter(gene_symbol %in% gene_list) %>%
    tidyr::pivot_longer(cols = 2:ncol(.), values_to = "TPM", names_to = "Sample_ID") %>%
    left_join(metadata, by = "Sample_ID") %>%
    group_by(genotype, gene_symbol) %>%
    summarise(avg = mean(TPM)) %>%
    filter(avg >= 1) %>%
    pull(gene_symbol) %>% 
    unique() 
}

MGnD_filtered <- gene_set_filter(subset(Clec7a_sig, direction == "up")$tracking_id)
Homeostasis_filtered <- gene_set_filter(subset(Clec7a_sig, direction == "down")$tracking_id)
KEGG_phagosome_filtered <- gene_set_filter(KEGG_phagosome)

## Function to generate signature score boxplot 
sig_score <- function(gene_list, title = "", suffix,legend_pos = "right") {
  
  ## Calculate signature score
  score_df <- tpm %>%
    filter(gene_symbol %in% gene_list) %>%
    tibble::column_to_rownames("gene_symbol") %>%
    colMeans() %>% as.data.frame() %>%
    tibble::rownames_to_column("Sample_ID") %>%
    dplyr::rename("score" = ".") %>% 
    left_join(metadata, by = "Sample_ID") %>%
    mutate(genotype = factor(genotype, c("Tim3_flox", "Tim3_cKO")))
  
  ## T-test on signature score
  stat.test <- score_df %>% t_test(score ~ genotype) %>%
    add_significance("p") %>%
    mutate(p.signif = ifelse(p >= 0.05, paste("p =", round(p, 3)), p.signif))
  
  boxplot <- score_df %>%
    ggplot(aes(x = genotype, y = log2(score), color = genotype,
               group = interaction(Age, genotype))) +
    geom_boxplot(size = 0.6, width = 0.3, outlier.alpha = 0) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
              fill = "white", alpha = 0.025, color = NA)+
    geom_point(aes(shape = sex), 
               size = 1.5, position = position_dodge(width = 0.3)) +
    stat_pvalue_manual(stat.test, size = ifelse(stat.test$p < 0.05, 5, 3),
                       y.position = max(log2(score_df$score)) * 1.002, 
                       bracket.size = 0.6, label = "p.signif", 
                       tip.length = 0, vjust = ifelse(stat.test$p < 0.05, 0.5, -0.1)) +
    theme_cowplot(font_size = 9) +
    scale_color_manual(values = c("#3b58a7", "red"),
                       breaks = c("Tim3_flox", "Tim3_cKO"), 
                       labels = c("control", "<i>Havcr2</i><sup>cKO</sup>")) +
    scale_x_discrete(breaks = c("Tim3_flox", "Tim3_cKO"), 
                     labels = c("control", "<i>Havcr2</i><sup>cKO</sup>"),
                     expand = expansion(add = c(0.4, 0.4))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.06))) +
    guides(color = guide_legend(order = 1),
           shape = guide_legend(order = 2)) +
    labs(x = NULL, y = "Log2 Signature Score", title = title) +
    theme(axis.text.x = element_markdown(vjust = 0),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 9.5, face = "bold", lineheight = 1),
          legend.text = element_markdown(),
          legend.position = legend_pos)

  return(list(score_df = score_df, boxplot = boxplot, p_val_ttest = stat.test))
}

## Fig 1j: MGnD signature score boxplot
MGnD_score <- sig_score(
  MGnD_filtered, title = "MGnD\nSignature Genes",
  suffix = "clec7a_up", legend_pos = "none")

## Fig 1k: Homeostasis microglia signature score boxplot
Homeostasis_score <- sig_score(
  Homeostasis_filtered, title = "Homeostasis\nSignature Genes",
  suffix = "clec7a_dn", legend_pos = "none")

## Fig 1l: KEGG phagosome signature score boxplot
KEGG_phagosome_score <- sig_score(
  KEGG_phagosome_filtered, title = "Phagosome\nSignature Genes",
  suffix = "KEGG Phagosome")


## Combine the three boxplots
sig_score_boxplots <- egg::ggarrange(
  MGnD_score$boxplot, Homeostasis_score$boxplot,
  KEGG_phagosome_score$boxplot, nrow = 1)
ggsave("figures/Fig.1h_1l_1month/Fig.1j_1l_1M_signature_boxplot.png",
       sig_score_boxplots, width = 5.5, height = 3, dpi = 400)
ggsave("figures/Fig.1h_1l_1month/Fig.1j_1l_1M_signature_boxplot.pdf",
       sig_score_boxplots, width = 5.5, height = 3)

MGnD_score$score_df %>%
  select(Sample_ID, genotype, sex, score) %>% 
  write.csv("Source_Data/Fig.1j_1M_signature_boxplot_MGnD.csv")
Homeostasis_score$score_df %>%
  select(Sample_ID, genotype, sex, score) %>% 
  write.csv("Source_Data/Fig.1k_1M_signature_boxplot_Homeostatic.csv")
KEGG_phagosome_score$score_df %>%
  select(Sample_ID, genotype, sex, score) %>% 
  write.csv("Source_Data/Fig.1l_1M_signature_boxplot_KEGG_phagosome.csv")
