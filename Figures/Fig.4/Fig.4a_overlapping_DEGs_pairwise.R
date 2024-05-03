library(corrplot)
library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(cowplot)
library(limma)
library(ggtext)

# Load functions for the analysis
source('bulkRNAseq/0.bulkRNAseq_functions.R')

# Palette 
pal <- c("darkgoldenrod1", "orangered", "#FB9A99", "darkred",
         brewer.pal(9, "Greens")[4], "forestgreen", "steelblue1", "steelblue4")
pal_text <- c("darkgoldenrod2", "orangered", "#FB9A99", "darkred",
              brewer.pal(9, "Greens")[5], "forestgreen", "steelblue1", "steelblue4")


# Read data -----------------------------------------


## (1) 1-month-old mice, Havcr2cKO vs Havcr2flox/flox
load("results/bulkRNAseq_results_ds2_1month.RData")
tim3_all <- res_ordered %>% correct_gene_symbol()  %>%
  select(gene_symbol, log2FoldChange, direction, padj)
tim3_DEG <- tim3_all %>% res_order_slice(thres = 0.1)

## (2) 3-month-old mice, phagocytosing control vs non-phagocytosing control microglia
load("results/bulkRNAseq_results_ds1_batch1_3month.RData")
phago_all <- results_batch1$phagoposvsneg_control %>%
  dplyr::rename("gene_symbol" = "gene_name") %>%
  correct_gene_symbol() %>%
  select(gene_symbol, log2FoldChange, direction, padj)
phago_DEG <- phago_all %>% res_order_slice(thres = 0.1)

## (3) Tgfbr2cKo vs control (Lund et al. 2018, PMID: 29662171)
TGFBRII_all <- read.csv("results/bulkRNAseq_results_TGFBRII_Lund_2018.csv") %>%
  select(gene_symbol, contains(".uG")) %>%
  dplyr::rename("log2FoldChange" = "log2fc.uG",
                "direction" = "dir.uG", "padj" = "padj.uG") %>%
  filter(rowSums(.[, 2:7]) != 0) %>%
  correct_gene_symbol()  %>%
  select(gene_symbol, log2FoldChange, direction, padj)

TGFBRII_DEG <- TGFBRII_all %>% res_order_slice(thres = 0.1)

## (4) Clec7a+ vs Clec7a- (Krasemann et al., 2017, PMID: 28930663)
Clec7a_all <- read.csv("results/bulkRNAseq_results_Clec7a_Krasemann_2017.csv") %>%
  dplyr::rename("gene_symbol" = "tracking_id", "log2FoldChange" = "log2FC") %>%
  correct_gene_symbol() %>%
  select(gene_symbol, log2FoldChange, direction, padj)
Clec7a_DEG <- Clec7a_all %>% res_order_slice(thres = 0.1)

## DEGs of the 4 comparisons
DEG_ls <-  list(`Havcr2cKO` = tim3_DEG,
                `Phagocytosing MG` = phago_DEG,
                `Tgfbr2cKO` = TGFBRII_DEG,
                `Clec7a+` = Clec7a_DEG)
## Full gene lists of the 4 comparisons
all_genes_ls <- list(`Havcr2cKO` = tim3_all,
                     `Phagocytosing MG` = phago_all,
                     `Tgfbr2cKO` = TGFBRII_all,
                     `Clec7a+` = Clec7a_all)

# Permutation test (pairwise) ---------
## Gene sets G1 and G2 have g1 and g2 DEGs, respectively
## Denote the actual number of overlaps between g1 and g2 as n_obs
## The permutation test p-values were computed by repeating the following steps 10,000 times
## (1) Randomly select n1 genes from G1, n2 genes from G2 
## (2) Count the number of overlap between the two sets of selected genes, denote as n_i, i=1, 2, ..., 10,000.
## If the majority of n_i > n_obs, then we used the right-tailed p-value (the sum of I{n_i>n_pos} divided by 10,000); 
## otherwise we used the left-tailed p-value (the sum of I{n_i<=n_pos} divided by 10,000); 
## The significance level of the permutation test was p<0.025. 

perm_test <- function(set1, set2, all_genes1, all_genes2, name1, name2, n_perm = 10000) {
  n_overlap <- inner_join(set1, set2, by = "gene_symbol",
                          suffix = c(".set1", ".set2")) %>%
    mutate(direction.set1 = paste(direction.set1, "in", name1),
           direction.set2 = paste(direction.set2, "in", name2)) %>%
    mutate(direction.set1 = forcats::fct_rev(direction.set1),
           direction.set2 = forcats::fct_rev(direction.set2)) %>%
    group_by(direction.set1, direction.set2) %>%
    summarise(n_overlap = n())
  
  set1_ls <- split(set1, set1$direction)
  set2_ls <- split(set2, set2$direction)
  
  set.seed(42)
  n_overlap_perm <- replicate(n_perm, {
    set1_perm <- lapply(set1_ls, function(set1_sub) {
      sample(all_genes1, size = nrow(set1_sub), replace = FALSE)
    })
    set2_perm <- lapply(set2_ls, function(set2_sub) {
      sample(all_genes2, size = nrow(set2_sub), replace = FALSE)
    })
    
    mapply(
      FUN = function(x, y) {length(intersect(x, y))},
      x = list(up.up = set1_perm$up, down.up = set1_perm$down, 
               up.down = set1_perm$up, down.down = set1_perm$down),
      y = list(up.up = set2_perm$up, down.up = set2_perm$up,
               up.down = set2_perm$down, down.down = set2_perm$down),
      USE.NAMES = TRUE) }) %>% 
    as.data.frame() %>% rownames_to_column("dirs") %>%
    separate(dirs, c("direction.set1", "direction.set2"), sep = "\\.") %>%
    mutate(direction.set1 = paste(direction.set1, "in", name1),
           direction.set2 = paste(direction.set2, "in", name2)) %>%
    mutate(direction.set1 = forcats::fct_rev(direction.set1),
           direction.set2 = forcats::fct_rev(direction.set2)) %>%
    pivot_longer(cols = 3:ncol(.), names_to = "index", values_to = "n_overlap_perm")
  
  p_overlap_perm <- n_overlap_perm %>% 
    left_join(n_overlap, by = c("direction.set1", "direction.set2")) %>%
    group_by(direction.set1, direction.set2) %>%
    summarise(p.perm = sum(n_overlap_perm > n_overlap)/n_perm) %>%
    mutate(p.perm_1sided = ifelse(p.perm > 0.5, 1 - p.perm, p.perm)) %>%
    as.data.frame() 
  
  return(list(n_overlap = n_overlap, p_overlap_perm = p_overlap_perm))
}

replace_compare_names <- function(x) {
  x <- str_replace(x, "cKO", "<sup>cKO</sup>") %>%
    str_replace("up in (?=Phago|Clec7a)", "<span style='color:white'>up in ") %>%
    str_replace("down in (?=Phago|Clec7a)", "<span style='color:white'>down in ") %>%
    str_replace("in ", "in <br/>") %>%
    str_replace("\\+", "<sup>+</sup>") %>%
    str_replace("Havcr2", "<i>Havcr2</i>") %>%
    str_replace("Tgfbr2", "<i>Tgfbr2</i>") %>%
    str_replace("Clec7a", "<i>Clec7a</i>") %>%
    str_replace("Phagocytosing MG", "phago<sup>+</sup></span>")
}


## Pairs of comparisons
compare_ls <- data.frame(x = rep(1:4, 4), y = c(sapply(1:4, rep, 4))) %>%
  filter(x >= y)

## Permutation test for pairwise overlapping genes
perm_test_pairwise <- sapply(1:nrow(compare_ls), function(i) {
  x <- compare_ls$x[i]; y <- compare_ls$y[i]
  perm_test(
    set1 = DEG_ls[[x]], set2 = DEG_ls[[y]],
    all_genes1 = all_genes_ls[[x]]$gene_symbol, all_genes2 = all_genes_ls[[y]]$gene_symbol,
    name1 = names(DEG_ls)[x], name2 = names(DEG_ls)[y],
    n_perm = 10000 )
}, simplify = FALSE)

save(perm_test_pairwise, file = "results/2023-12-28.perm_test_pairwise_padj.1.RData")

# load("results/2023-12-28.perm_test_pairwise_padj.1.RData")

## No. of overlapping genes of each pair of the comparisons
perm_test_pairwise_n_df <- lapply(perm_test_pairwise, "[[", "n_overlap") %>%
  do.call(rbind, .)

## Permutation test p-value of each pair of the comparisons
perm_test_pairwise_p_df <- lapply(perm_test_pairwise, "[[", "p_overlap_perm") %>%
  do.call(rbind, .)

perm_test_pairwise_p_df %>%
  full_join(perm_test_pairwise_n_df, by = c("direction.set1", "direction.set2")) %>%
  write.csv("Source_Data/Fig.4a_perm_test.csv")

perm_test_pairwise_df <- perm_test_pairwise_p_df %>%
  full_join(perm_test_pairwise_n_df, by = c("direction.set1", "direction.set2")) %>%
  mutate_at(vars(starts_with("direction")), replace_compare_names) %>%
  mutate_at(vars(starts_with("direction")), factor, 
            replace_compare_names(paste0(
              c(sapply(c("up in ", "down in "), rep, 4)),
              rep(names(all_genes_ls), 2)))) %>%
  select(direction.set1, direction.set2, n_overlap, p = p.perm, p.1sided = p.perm_1sided)

# Fig. 4a Dotplot displaying the overlaps between each pair of DEGs up- and down-regulated in the below 4 comparisons -----------------------------------
## Circle cize and text: no. of overlapping genes 
## Color: permutation test p-value

## Function for changing the colors of the strips in ggplot
strip_col <- function(p, cols) {
  g <- ggplot_gtable(ggplot_build(p))
  
  strips <- which(grepl('strip-', g$layout$name))
  
  strips_pal <- rep(cols, 2)
  
  for (i in seq_along(strips)) {
    k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
    l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
    g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- strips_pal[i]
  }
  return(g)
}

df <- perm_test_pairwise_df %>%
  `colnames<-`(c("direction.set2", "direction.set1", "n_overlap", "p", "p.1sided")) %>%
  select(direction.set1, direction.set2, n_overlap, p, p.1sided) %>%
  rbind(perm_test_pairwise_df) %>%
  distinct() %>%
  mutate(fct.set1 = as.numeric(direction.set1),
         fct.set2 = as.numeric(direction.set2)) %>%
  filter(fct.set1 <= fct.set2) %>%
  mutate(fill = ifelse(p > 0.5, "#B3CDE3", "#FBB4AE")) %>%
  mutate(fill = ifelse(fct.set1 == fct.set2, "white", fill)) %>%
  mutate(p_cat = case_when(p.1sided < 1E-5 ~ "< 0.00001",
                           p.1sided < 1E-3 ~ "< 0.001",
                           p.1sided < 1E-2 ~ "< 0.01",
                           p.1sided < 0.025 ~ "< 0.025",
                           TRUE ~ "ns")) %>%
  mutate(p_cat = ifelse(p > 0.5 & p_cat != "ns", paste0(p_cat, " (dissimilar)"), p_cat)) %>%
  mutate(p_cat = factor(p_cat, c("< 0.00001", "< 0.001", "< 0.01", "< 0.025",
                                 "< 0.025 (dissimilar)", "< 0.01 (dissimilar)", "< 0.00001 (dissimilar)", "ns")))

p <- df %>%
  mutate(n_overlap = ifelse(fct.set1 == fct.set2, NA, n_overlap)) %>%
  ggplot(aes(x = direction.set1, y = direction.set2)) +
  geom_rect(xmin = 0.35, xmax = 1.6, ymin = 0.35, ymax = 1.6,
            alpha = 0.2, fill = "white",
            size = 0.6, color = "grey55") +
  geom_point(aes(size = n_overlap, fill = p_cat), shape = 21, color = "grey40") +
  geom_text(data = subset(df, !(is.na(p) | direction.set1 == direction.set2)), 
            aes(label = paste("n =", n_overlap)), size = 3, vjust = 3.8, color = "grey15", fontface = 3) +
  geom_text(data = subset(df, direction.set1 == direction.set2) %>% arrange(direction.set1), 
            aes(label = paste("N =", n_overlap)),
            color = pal_text, size = 4, fontface = 4) +
  geom_segment(data = subset(df, is.na(p)),
               x = 0, xend = 2, y = 2, yend = 0,
               color = "grey55", size = 0.3) +
  scale_y_discrete(limits = rev) +
  theme_cowplot(font_size = 12) +
  guides(fill = guide_legend(override.aes = list(size = 3.5), order = 1),
         size = guide_legend(override.aes = list(fill = "grey15"), order = 2)) + 
  scale_size(range = c(2, 12), breaks = seq(50, 600, 100), name = "# Overlaps") +
  facet_grid(direction.set2 ~ direction.set1, scales = "free", switch = "both") +
  scale_fill_manual(values = c(rev(brewer.pal(9, "YlOrRd")[c(1, 3, 5, 8)]), 
                               brewer.pal(9, "Blues")[c(2, 5, 9)], "grey60"),
                    name = "p", drop = FALSE) +
  theme(strip.text.x = element_markdown(face = "bold", lineheight = 1.15, margin = margin(3, 0, 3, 0)),
        strip.text.y.left = element_markdown(face = "bold", lineheight = 1.15, margin = margin(0, 3, 0, 3)),
        strip.background.x = element_rect(color = "grey35", size = 0.6),
        strip.background.y = element_rect(color = "grey35", size = 0.6),
        panel.border = element_blank(),
        aspect.ratio = 1, axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.box = "horizontal",
        legend.position = c(0.55, 0.82),
        legend.text = element_text(size = 11),
        panel.spacing = unit(0.001, "mm"),
        axis.line = element_blank()) +
  labs(x = NULL, y = NULL)

g <- strip_col(p, cols = pal)

ggsave("figures/Fig.4a_perm_test_pairwise.png", g, width = 7, height = 7, dpi = 400)
ggsave("figures/Fig.4a_perm_test_pairwise.pdf", g,  width = 7, height = 7)

