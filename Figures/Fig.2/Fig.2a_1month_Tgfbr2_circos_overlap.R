source('bulkRNAseq/0.bulkRNAseq_functions.R')
ptions(scipen=999)

## palette
pal_text <- as.list(
  setNames(c("darkgoldenrod2", "orangered", "#FB9A99", "darkred", ggpubfigs::friendly_pals$nickel_five[5],
             brewer.pal(9, "Greens")[5], "forestgreen", "steelblue1", "steelblue4", "grey15"),
           c(paste0(c("Havcr2cKO", "phago MG", "Tgfbr2cKO", "Clec7a+", "5XFAD"), ".up"),
             paste0(c("Havcr2cKO", "phago MG", "Tgfbr2cKO", "Clec7a+", "5XFAD"), ".down"))))

## genes to highlight
gene_list_Tgfbr2 <- c("Havcr2", "Axl", "Ly9", "Clec7a", "Cd63", "Cxcl16", 
                      "H2-K1", "H2-D1", "Siglech", "Csf1r")

## MGnD/Homeostasis signatures (Top 100 DEGs of up/down-regulated in Clec7a+ vs Clec7a-)
Clec7a_sig <- read.csv("results/bulkRNAseq_results_Clec7a_Krasemann_2017.csv") %>%
  arrange(padj) %>% group_by(direction) %>% dplyr::slice(1:100) %>% pull(tracking_id)


# Load differential expression analysis results --------------------------------

## (1) 1-month-old mice, Havcr2cKO vs Havcr2flox/flox microglia
load("results/bulkRNAseq_results_ds2_1month.RData")
ds2_tim3_sig <- res_ordered %>% res_order_slice(thres = 0.1)

## (2) Tgfbr2KO vs control
TGFBRII_all <- read.csv("results/bulkRNAseq_results_TGFBRII_Lund_2018.csv") %>%
  select(gene_symbol, contains(".uG")) %>%
  dplyr::rename("log2FoldChange" = "log2fc.uG",
                "direction" = "dir.uG", "padj" = "padj.uG") %>%
  filter(!is.na(padj)) %>%
  correct_gene_symbol()  %>%
  select(gene_symbol, log2FoldChange, direction, padj)

TGFBRII_sig <- TGFBRII_all %>% res_order_slice(thres = 0.1)


# Premutation test on the overlapping genes ------------------------------------
## Gene sets G1 and G2 have g1 and g2 DEGs, respectively
## Denote the actual number of overlaps between g1 and g2 as n_obs
## The permutation test p-values were computed by repeating the following steps 10,000 times
## (1) Randomly select n1 genes from G1, n2 genes from G2 
## (2) Count the number of overlap between the two sets of selected genes, denote as n_i, i=1, 2, ..., 10,000.
## If the majority of n_i > n_obs, then we used the right-tailed p-value (the sum of I{n_i>n_pos} divided by 10,000); 
## otherwise we used the left-tailed p-value (the sum of I{n_i<=n_pos} divided by 10,000); 
## The significance level of the permutation test was p<0.025. 

overlap_ds2_tim3_Tgfbr2_sig <- perm_test(
  set1 = ds2_tim3_sig, set2 = TGFBRII_sig,
  all_genes1 = res_ordered$gene_symbol, 
  all_genes2 = TGFBRII_all$gene_symbol,
  name1 = "Tim3cKO", name2 = "Tgfbr2KO",
  prefix1 = "DS2_1M_", prefix2 = "", 
  n_perm = 10000
)

write.csv(overlap_ds2_tim3_Tgfbr2_sig$p_overlap_perm,
          "Source_Data/Fig.2a_perm_test.csv")

# Fig. 2a: Circos plot comparison of up- and down-regulated DEGs ---------------
## Comparison 1: Havcr2cKO vs. Havcr2 flox/flox from 1-month-old mice  
## Comparison 2: Tgfbr2cKO vs. control 
circo_overlap(res1 = ds2_tim3_sig, res2 = TGFBRII_sig, Fig_num = "2a",
              name1 = "Havcr2cKO", name2 = "Tgfbr2cKO", 
              suffix = "ds2_tim3_Tgfbr2KO_sub", size = 7.5, gene_width = 25,
              num_cex = 0.9, gene_cex = 0.75, gene_width_short = 3, 
              big_gap = 5, small_gap = 0.4, 
              gene_list = c(gene_list_Tgfbr2, Clec7a_sig, "H2-Q4", "Hexb"),
              degree = 100, step = 300, show_selected_genes = TRUE, circo_pal = pal_text,
              overlap_perm_test = overlap_ds2_tim3_Tgfbr2_sig, 
              nudge_x = c(0, 0, 50, -50), niceFacing = FALSE)

