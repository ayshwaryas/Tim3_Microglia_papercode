source('bulkRNAseq/0.bulkRNAseq_functions.R')
ptions(scipen=999)

## palette
pal_text <- as.list(
  setNames(c("darkgoldenrod2", "orangered", "#FB9A99", "darkred", ggpubfigs::friendly_pals$nickel_five[5],
             brewer.pal(9, "Greens")[5], "forestgreen", "steelblue1", "steelblue4", "grey15"),
           c(paste0(c("Havcr2cKO", "phago MG", "Tgfbr2cKO", "Clec7a+", "5XFAD"), ".up"),
             paste0(c("Havcr2cKO", "phago MG", "Tgfbr2cKO", "Clec7a+", "5XFAD"), ".down"))))

## genes to highlight
gene_list_phago <- c("Havcr2", "Cd9", "Cst7", "Cstb", "Cxcl16", "Ly9", "Lyz2", "B2m", 
                     "H2-M3", "H2-M2", "H2-M10.2", "H2-Oa", "H2-M10.1", "H2-M5", 
                     "H2-Q4", "H2-Aa", "H2-M10.6", "H2-M10.5", "H2-M1", "H2-M11", 
                     "H2-DMb2", "H2-DMa", "H2-Ob", "H2-M10.4", "H2-T24", "H2-T3", 
                     "H2-T22", "H2-M10.3", "H2-Q7", "H2-Eb1", "H2-K1", "H2-M9", 
                     "H2-T23", "H2-Q10", "H2-Eb2", "H2-Q6", "H2-D1", "H2-Ab1", 
                     "H2-Ke6", "H2-Q1", "H2-DMb1", "H2-Q2")

## KEGG Phagosome signature
KEGG_phagosome <- read.table("data/signatures/KEGG_phagosome.txt", sep = "\t") %>%
  mutate(gene = str_extract(V2, ".*(?=;)")) %>%
  mutate(details = str_extract(V2, "(?<=; ).*")) %>% select(-V2) %>% pull(gene)

## MGnD/Homeostasis signatures (Top 100 DEGs of up/down-regulated in Clec7a+ vs Clec7a-)
Clec7a_sig <- read.csv("results/bulkRNAseq_results_Clec7a_Krasemann_2017.csv") %>%
  arrange(padj) %>% group_by(direction) %>% dplyr::slice(1:100) %>% pull(tracking_id)


# Load differential expression analysis results -----------------------------------------
load("results/bulkRNAseq_results_ds1_batch1_3month.RData")

## (1) 3-month-old mice, control phagocyosing vs control non-phagocytosing microglia
ds1_batch1_tim3_sig <- results_batch1$Tim3cKOvscontrol_phagoneg %>%
  dplyr::rename("gene_symbol" = "gene_name") %>% 
  res_order_slice(thres = 0.1) 

## (2) 3-month-old mice, Havcr2cKO non-phagocytosing vs control non-phagocytosing microglia
ds1_batch1_phago_sig <- results_batch1$phagoposvsneg_control %>%
  dplyr::rename("gene_symbol" = "gene_name") %>% 
  res_order_slice(thres = 0.1) 


# Premutation test on the overlapping genes -----------------------------------
## Gene sets G1 and G2 have g1 and g2 DEGs, respectively
## Denote the actual number of overlaps between g1 and g2 as n_obs
## The permutation test p-values were computed by repeating the following steps 10,000 times
## (1) Randomly select n1 genes from G1, n2 genes from G2 
## (2) Count the number of overlap between the two sets of selected genes, denote as n_i, i=1, 2, ..., 10,000.
## If the majority of n_i > n_obs, then we used the right-tailed p-value (the sum of I{n_i>n_pos} divided by 10,000); 
## otherwise we used the left-tailed p-value (the sum of I{n_i<=n_pos} divided by 10,000); 
## The significance level of the permutation test was p<0.025. 

overlap_ds1_tim3_phago_sig <- perm_test(
  set1 = ds1_batch1_tim3_sig, set2 = ds1_batch1_phago_sig,
  all_genes1 = results_batch1$Tim3cKOvscontrol_phagoneg$gene_name,
  all_genes2 = results_batch1$phagoposvsneg_control$gene_name,
  name1 = "Tim3cKO", name2 = "phago+", 
  prefix1 = "DS1_Batch1_", prefix2 = "DS1_Batch1_",
  n_perm = 10000
)

# Fig. 3i: Circos plot comparison of up- and down-regulated DEGs from 3-month-old mice ------------------------------
## Comparison 1: Havcr2cKO vs control non-phagocytosing microglia
## Comparison 2: control phagocyosing vs control non-phagocytosing microglia

circo_overlap(res1 = ds1_batch1_tim3_sig,
              res2 = ds1_batch1_phago_sig,
              name1 = "Havcr2cKO", name2 = "phago MG", 
              suffix = "atleast300_KEGG_phago_sub", size = 8.5, 
              gene_width = 28, gene_width_short = 1.5,
              num_cex = 0.9, gene_cex = 0.75, p_cex = 0.8,
              big_gap = 3, small_gap = 0.15, circo_pal = pal_text, 
              nudge_x = c(0, 0, 50, -50),
              degree = 100, gene_list = c(gene_list_phago, Clec7a_sig),
              step = 300, phago_sig = KEGG_phagosome, 
              show_selected_genes = TRUE,
              overlap_perm_test = overlap_ds1_tim3_phago_sig,
              conflict_facing = "inside")

