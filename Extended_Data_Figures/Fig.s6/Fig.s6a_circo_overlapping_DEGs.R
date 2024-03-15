source('bulkRNAseq/0.bulkRNAseq_functions.R')

## Palettes
pal <- c("darkgoldenrod1", "orangered", "#FB9A99", "darkred",
         brewer.pal(9, "Greens")[4], "forestgreen", "steelblue1", "steelblue4")
pal_text <- c("darkgoldenrod3", "orangered", "#fa7573", "darkred",
              brewer.pal(9, "Greens")[6], "forestgreen", "steelblue2", "steelblue4")


# Read data --------------------------------------------------------------------

## 1. Differential expression analysis results ---------------------------------
## (1) 1-month-old mice, Havcr2cKO vs Havcr2flox/flox --------------------------
load("results/bulkRNAseq_results_ds2_1month.RData")
tim3_all <- res_ordered %>% correct_gene_symbol() %>%
  select(gene_symbol, log2FoldChange, direction, padj) 
tim3_DEG <- res_order_slice(tim3_all, thres = 0.05)

## (2) 3-month-old mice, phagocytosing control vs non-phagocytosing control microglia ---------
load("results/bulkRNAseq_results_ds1_batch1_3month.RData")
phago_all <- results_batch1$phagoposvsneg_control %>%
  dplyr::rename("gene_symbol" = "gene_name") %>%
  correct_gene_symbol() %>%
  select(gene_symbol, log2FoldChange, direction, padj)
phago_DEG <- res_order_slice(phago_all, thres = 0.05)

## (3) Tgfbr2cKo vs control (Lund et al. 2018, PMID: 29662171) -----------------
TGFBRII_all <- read.csv("results/bulkRNAseq_results_TGFBRII_Lund_2018.csv") %>%
  select(gene_symbol, contains(".uG")) %>%
  dplyr::rename("log2FoldChange" = "log2fc.uG",
                "direction" = "dir.uG", "padj" = "padj.uG") %>%
  filter(rowSums(.[, 2:7]) != 0) %>%
  correct_gene_symbol()  %>%
  select(gene_symbol, log2FoldChange, direction, padj)  %>%
  distinct(gene_symbol, direction, .keep_all = TRUE) 

TGFBRII_DEG <- res_order_slice(TGFBRII_all, thres = 0.05)

## (4) Clec7a+ vs Clec7a- (Krasemann et al., 2017, PMID: 28930663) -------------
Clec7a_all <- read.csv("results/bulkRNAseq_results_Clec7a_Krasemann_2017.csv") %>%
  dplyr::rename("gene_symbol" = "tracking_id", "log2FoldChange" = "log2FC") %>%
  correct_gene_symbol() %>%
  select(gene_symbol, log2FoldChange, direction, padj) 

Clec7a_DEG <- res_order_slice(Clec7a_all, thres = 0.05)

## Combining list of DEGs ------------------------------------------------------
all_genes_ls <- list(`Havcr2cKO` = tim3_all$gene_symbol,
                     `Phagocytosing MG` = phago_all$gene_symbol,
                     `Tgfbr2cKO` = TGFBRII_all$gene_symbol,
                     `Clec7a+` = Clec7a_all$gene_symbol)

DEG_ls <- list(`Havcr2cKO` = tim3_DEG,
               `Phagocytosing MG` = phago_DEG,
               `Tgfbr2cKO` = TGFBRII_DEG,
               `Clec7a+` = Clec7a_DEG)


## Direction of change of DEGs shared by at least 3 datasets -------------------

Figs6a_DEG_dir <- DEG_ls %>%
  data.table::rbindlist(idcol = "dataset") %>%
  mutate(direction = factor(direction, c("up", "down"))) %>%
  arrange(dataset, direction) 

gene_list_Tgfbr2 <- c("Havcr2", "Axl", "Ly9", "Clec7a", "Cd63", "Cxcl16", 
                      "H2-K1", "H2-D1", "Siglech", "Csf1r")
gene_list_phago <- c("Havcr2", "Cd9", "Cst7", "Cstb", "Cxcl16", "Ly9", "Lyz2", 
                     unique(subset(Figs6a_DEG_dir, str_detect(gene_symbol, "H2\\-"))$gene_symbol))

# Permutation test -------------------------------------------------------------

DEG_set_ls <- lapply(DEG_ls, function(x) {split(x, x$direction)})

n_overlap <- DEG_ls %>%
  data.table::rbindlist(idcol = "dataset") %>%
  group_by(gene_symbol) %>%
  filter(n() >= 3) %>% select(-log2FoldChange) %>%
  pivot_wider(names_from = "dataset", values_from = "direction") %>%
  mutate_at(2:5, factor) %>%
  mutate_at(2:5, function(x) { factor(x, str_sort(levels(x), decreasing = TRUE))}) %>%
  group_by(`Havcr2cKO`, `Phagocytosing MG`, `Tgfbr2cKO`, `Clec7a+`) %>%
  summarise(n_overlap = n()) %>%
  mutate(dir = paste(Havcr2cKO, `Phagocytosing MG`, Tgfbr2cKO, `Clec7a+`, sep = "_")) %>%
  ungroup()


set.seed(42)
n_perm <- 10000 # number of permutations
n_overlap_perm <- replicate(n_perm, {
  set_perm <- sapply(names(DEG_set_ls), function(x) {
    sample_up <- data.frame(gene_symbol = sample(all_genes_ls[[x]], size = nrow(DEG_set_ls[[x]]$up), replace = FALSE))
    sample_dn <- data.frame(gene_symbol = sample(all_genes_ls[[x]], size = nrow(DEG_set_ls[[x]]$down), replace = FALSE))
    list(up = sample_up, down = sample_dn) %>% 
      data.table::rbindlist(idcol = 'direction')}, 
    USE.NAMES = TRUE, simplify = FALSE) %>%
    data.table::rbindlist(idcol = "dataset") %>%
    mutate(dataset = factor(dataset, names(DEG_ls))) %>%
    mutate(direction = factor(direction, c("up", "down"))) 
  
  set_perm_summ <- set_perm %>%
    group_by(gene_symbol, dataset) %>%
    mutate(index = 1:n()) %>%
    group_by(dataset, direction, index) %>%
    pivot_wider(names_from = "dataset", values_from  = "direction") %>%
    mutate(dir = paste(Havcr2cKO, `Phagocytosing MG`, Tgfbr2cKO, `Clec7a+`, sep = "_")) %>%
    filter(str_count(dir, "NA") <= 1) %>%
    group_by(dir) %>%
    summarise(n_overlap_perm = n()) 
  
  return(set_perm_summ) },
  simplify = FALSE)

p_overlap_perm <- n_overlap_perm %>% 
  data.table::rbindlist(idcol = "iter") %>%
  inner_join(n_overlap %>% select(dir, n_overlap), by = "dir") %>%
  group_by(dir) %>%
  summarise(p.perm = sum(n_overlap_perm > n_overlap)/n_perm) %>%
  mutate(p.perm = ifelse(p.perm > 0.5, 1 - p.perm, p.perm))

# save(n_overlap_perm, p_overlap_perm, file = "results/Fig.s6a_perm_test_results.05.RData")

# Circos plot ------------------------------------------------------------------

Figs6a_circos <- function(
    DEG_overlaps = NULL, p_overlap_perm = NULL, highlight_genes = NULL,
    big.gap = 1.5, small.gap = 0.15, p.perm_thres = 0.05,
    dir_cex = 0.7, gene_cex = c(0.6, 0.55), p_perm_cex = c(0.8, 0.8), label_cex = 0.5,
    gene_width = c(60, 60), track.height = 0.07,
    y_p_perm = 1.5, y_dir = -0.5, degree = 100, 
    size = 10, suffix = "", fig_num = "s6a") {
  
  names <- paste0(rep(names(DEG_ls), 2), c(rep(".up", 4), rep(".down", 4)))
  
  DEG_overlaps_circo <- DEG_overlaps %>%
    group_by(gene_symbol) %>% mutate(n = n()) %>%
    filter(n >= 3) %>% ungroup %>%
    mutate(from = paste(dataset, direction, sep = "."),
           to = gene_symbol) %>%
    mutate(width = ifelse(gene_symbol %in% c(highlight_genes), gene_width[1], gene_width[2])) %>%
    mutate(value = 1, value2 = width / n * 3) %>%
    select(from, to, value, value2)
  
  count_ds <- DEG_overlaps  %>%
    mutate(from = paste(dataset, direction, sep = ".")) %>% dplyr::count(from) %>%
    full_join(dplyr::count(DEG_overlaps_circo, from), by = "from",
              suffix = c(".total", ".overlap")) %>%
    mutate(n.overlap = replace_na(n.overlap, 0)) %>%
    mutate(to = from, value = (n.total - n.overlap)/2, value2 = value) %>% 
    select(from, to, value, value2)
  
  group_genes <- DEG_overlaps_circo %>%
    mutate(from = factor(from, names)) %>% arrange(from) %>% 
    mutate(from_num = as.numeric(from)) %>%
    mutate(dir_num = ifelse(str_detect(from, "up"), 1, -1)) %>%
    separate(from, sep = "\\.", into = c("dataset", "direction")) %>%
    select(to, dataset, direction, from_num, dir_num) %>% 
    pivot_wider(names_from = "dataset", values_from = c("direction", "dir_num", "from_num")) %>% 
    unite(group, starts_with("direction_")) %>% 
    mutate(n = str_count(group, "up") * 1 - str_count(group, "down")) %>%
    mutate(dir_num = rowSums(.[, 3:6], na.rm = TRUE)) %>%
    unite(group_num, starts_with("from_num_"), sep = "\n") %>%
    arrange(sign(n), desc(abs(n)), group_num) %>%
    left_join(p_overlap_perm, by = c("group" = "dir")) %>%
    filter(p.perm < p.perm_thres) %>%
    select(-starts_with("dir_num_"), -p.perm) %>%
    arrange(sign(n), desc(abs(n)), group_num)
  
  ## grouping the genes
  group <- c(str_extract(names[1:4], "\\.up|\\.down"), 
             group_genes[group_genes$n > 0, ]$group,
             str_extract(names[8:5], "\\.up|\\.down"),
             group_genes[group_genes$n <= 0, ]$group)
  group <- structure(
    factor(group, unique(group)),  
    names = c(names[1:4], group_genes[group_genes$n > 0, ]$to,
              names[8:5], group_genes[group_genes$n <= 0, ]$to))
  
  grid.col.ds <- setNames(pal, names)
  grid.col <- c(grid.col.ds, setNames(rep("grey", nrow(group_genes)), group_genes$to))
  grid.col.text.ds <- setNames(pal_text, names)
  
  circo_df <- rbind(subset(DEG_overlaps_circo, to %in% group_genes$to), count_ds) 
  
  ## add significance levels to permutation test p-values
  p_overlap_perm_circo <- p_overlap_perm %>% 
    rstatix::add_significance("p.perm",
                              cutpoints = c(0, 1e-05, 1e-04, 0.001, 0.01, 0.025, 1),
                              symbols = c("*****", "****", "***", "**", "*", "ns")) 
  
  pdf(paste0("figures/Fig.", fig_num, "_circos", suffix,  ".pdf"), width = size, height = size)
  circos.clear()
  circos.par(clock.wise = TRUE, start.degree = degree, message = FALSE)
  chordDiagram(circo_df, grid.col = grid.col,
               link.visible = circo_df[[1]] != circo_df[[2]],
               self.link = 2, 
               big.gap = big.gap, small.gap = small.gap,
               order = names(group), group = group,
               link.sort = TRUE, link.decreasing = TRUE,
               transparency = 0.05,
               annotationTrack = c("grid"), 
               preAllocateTracks = list(
                 ## tracks showing the gene symbols
                 list(track.height = track.height),
                 ## 4 tracks indicating if the gene is DE in a specific comparison
                 list(track.height = 0.015, track.margin = c(0, 0.005), bg.border = "white"),
                 list(track.height = 0.015, track.margin = c(0, 0.005), bg.border = "white"),
                 list(track.height = 0.015, track.margin = c(0, 0.005), bg.border = "white"),
                 list(track.height = 0.015, track.margin = c(0, 0.005), bg.border = "white")
               )
  )
  
  ## gene symbols
  for(si in group_genes$to) {
    xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
    ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
    dir <- group_genes$group[group_genes$to == si]
    
    if (si %in% highlight_genes & !str_detect(dir, "down")) {
      col = "red"; face = 4; cex = gene_cex[1]}
    else if (si %in% highlight_genes & !str_detect(dir, "up")) {
      col = "blue"; face = 4; cex = gene_cex[1]}
    else {col = "grey50"; face = 3; cex = gene_cex[2]}
    
    circos.text(mean(xlim), ylim[1], si, sector.index = si, track.index = 1, 
                facing = "clockwise", adj = c(0, 0.5), col = col,
                cex = cex, niceFacing = TRUE, font = face)
  }
  
  ## comparisons (e.g. Up in Havcr2cKO)
  for(si in unique(circo_df$from)) {
    
    circos.axis(h = 0, major.at = seq(0, 3000, 500),
                labels.cex = label_cex, labels.facing = "clockwise", 
                labels.niceFacing = FALSE,
                sector.index = si, track.index = 5)
    
    si_label <- str_remove(si, "\\.up|\\.down")
    si_dir <- paste0(str_to_title(str_extract(si, "up|down")), " in")
    
    xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
    ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
    
    if(!str_detect(si_label, "Phagocytosing")) {
      if(si_dir == "Up in")  {x = mean(xlim) + c(-150, 150)}
      if(si_dir == "Down in")  {x = mean(xlim) + c(-150, 150)}
      if(si == "Havcr2cKO.up") {x = x + 100}
      if(si == "Tgfbr2cKO.down") {x = x + c(0, 70)}
      si_sup <- str_extract(si_label, "cKO|\\+|\\-")
      si_base <- str_remove(si_label, "cKO|\\+|\\-")
      label_dir <- substitute(bold(dir), list(dir = si_dir)) 
      label <- c(label_dir, substitute(bolditalic(x)^bold(y), list(x = si_base, y = si_sup))) %>%
        as.expression()
    }
    else { x = mean(xlim); y = 0.5; label <- paste(si_dir, "\nPhago MG") }
    
    circos.text(x, y_dir, label, adj = 0,
                sector.index = si, track.index = 1, col = grid.col.text.ds[si],
                facing = "clockwise", font = 2, 
                cex = dir_cex, niceFacing = FALSE)
  }
  
  for(si in unique(group_genes$group)) {
    si_ls <- unlist(str_split(si, "_"))
    for(i in which(si_ls != "NA")) {
      highlight.sector(group_genes[group_genes$group == si, ]$to,
                       track.index = i + 1, col = grid.col.ds[paste(names(DEG_ls)[i], si_ls[i], sep = ".")],
                       cex = 0.1, text.col = "white", niceFacing = TRUE)}
  }
  
  group_ls <- split(group, group)
  group_mid <- group_ls %>%
    lapply(function(x) {
      index <- ceiling(length(x)/2)
      if(nchar(as.character(names(x[index]))) > 10 & index >= 3) {index <- index + 2}
      names(x[index])}) 
  
  
  for(x in names(group_mid)) {
    if(x %in% c(".up", ".down")) { next }
    si <- group_mid[[x]]
    xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
    ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
    
    label <- subset(p_overlap_perm_circo, dir == x)$p.perm.signif
    if(length(group_ls[[x]]) == 2) {x = mean(xlim)}
    circos.text(mean(xlim), y_p_perm, label, sector.index = si, track.index = 1, adj = 0,
                cex = ifelse(label == "ns", p_perm_cex[2], p_perm_cex[1]), niceFacing = TRUE, font = 2, facing = "clockwise")
  }
  
  
  dev.off()
} 

Figs6a_circos(DEG_overlaps = Figs6a_DEG_dir,  p_overlap_perm = p_overlap_perm,
              highlight_genes = c(gene_list_phago, gene_list_Tgfbr2, "Apoe", "Sall1", "Cd33"),
              dir_cex = 0.77, y_dir = -0.3, 
              p.perm_thres = 1, fig_num = "s6a")