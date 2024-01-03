source('bulkRNAseq/0.bulkRNAseq_functions.R')
ptions(scipen=999)

## palette
pal_text <- as.list(
  setNames(c("darkgoldenrod2", "orangered", "#FB9A99", "darkred", ggpubfigs::friendly_pals$nickel_five[5],
             brewer.pal(9, "Greens")[5], "forestgreen", "steelblue1", "steelblue4", "grey15"),
           c(paste0(c("Havcr2cKO", "phago MG", "Tgfbr2cKO", "Clec7a+", "5XFAD"), ".up"),
             paste0(c("Havcr2cKO", "phago MG", "Tgfbr2cKO", "Clec7a+", "5XFAD"), ".down"))))

## genes to highlight
gene_list_5XFAD <- c("Havcr2", "Axl", "Ccl6", "Cxcl16", "Cd9", "Cd81", "Lyz2",
                     "Ctse", "Ctsd", "Ctsz", "Ctsa", "Ctsq", "Cts7", "Cts6", "Ctsl", 
                     "Ctsb", "Ctsw", "Ctso", "Ctsk", "Ctsc", "Ctss", "Ctsg", "Ctsj", 
                     "Ctsr", "Ctsll3", "Cts8", "Cts3", "Ctsm") 

## MGnD/Homeostasis signatures (Top 100 DEGs of up/down-regulated in Clec7a+ vs Clec7a-)
Clec7a_sig <- read.csv("results/bulkRNAseq_results_Clec7a_Krasemann_2017.csv") %>%
  arrange(padj) %>% group_by(direction) %>% dplyr::slice(1:100) %>% pull(tracking_id)


# Load differential expression analysis results --------------------------------
load("results/bulkRNAseq_results_ds3_4month.RData")

## (1) 4-month-old female mice, Havcr2cKO vs control -------
ds3_tim3_F_DEGs <- results_ds3_F$Tim3cKO_vs_control %>%
  dplyr::rename("gene_symbol" = "gene_name") %>%
  res_order_slice(flip = FALSE, thres = 0.1)

## (2) 4-month-old female mice, 5XFAD vs control -------
ds3_5XFAD_F_DEGs <- results_ds3_F$`5XFAD_vs_control` %>%
  dplyr::rename("gene_symbol" = "gene_name") %>%
  res_order_slice(flip = FALSE, thres = 0.1)


# Premutation test on the overlapping genes -----------------------------------
## Gene sets G1 and G2 have n_obs overlapping genes
## The permutation test p-values were computed by repeating the following steps 10,000 times
## (1) Randomly select n1 genes from G1, n2 genes from G2 
## (2) Count the number of overlap between the two sets of selected genes, denote as n_i, i=1, 2, ..., 10,000.
## If the majority of n_i > n_obs, then we used the right-tailed p-value (the sum of I{n_i>n_pos} divided by 10,000); 
## otherwise we used the left-tailed p-value (the sum of I{n_i<=n_pos} divided by 10,000); 
## The significance level of the permutation test was p<0.025. 

overlap_ds3_tim3_5XFAD_F_sig <- perm_test(
  set1 = ds3_tim3_F_DEGs, set2 = ds3_5XFAD_F_DEGs,
  all_genes1 = results_ds3_F$Tim3cKO_vs_control$gene_name, 
  all_genes2 = results_ds3_F$`5XFAD_vs_control`$gene_name,
  name1 = "Tim3cKO", name2 = "5XFAD", n_perm = 10000,
  prefix1 = paste0("DS3_F_"),
  prefix2 = paste0("DS3_F_")
)

# Extended Data Fig. 11a: Circos plot comparison of up- and down-regulated DEGs from 4-month-old mice ------------------------------
## Comparison 1: Havcr2cKO vs control
## Comparison 2: 5XFAD vs control

circo_overlap(res1 = ds3_tim3_F_DEGs, res2 = ds3_5XFAD_F_DEGs,
              name1 = "Havcr2cKO", name2 = "5XFAD", circo_pal = pal_text,
              overlap_perm_test = overlap_ds3_tim3_5XFAD_F_sig,
              gene_cex = 0.7, small_gap = 1, suffix = "top300_F",
              size = 7.5, gene_width = 6, degree = 90,
              gene_list = setdiff(c(gene_list_5XFAD, Clec7a_sig), "Apoe"))
