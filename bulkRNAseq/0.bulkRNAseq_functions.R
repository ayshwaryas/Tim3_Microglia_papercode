library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(cowplot)

## Functions for correcting gene symbols ---------------------------------------
correct_gene_symbol <- function(df, var = "gene_symbol") {
  df %>%
    mutate(gene_symbol = case_when(
      gene_symbol == "02-mars" ~ "Mars2",
      str_detect(gene_symbol, "-Mar") ~ paste0("March", gene_symbol),
      str_detect(gene_symbol, "-Sep") ~ paste0("Sept", gene_symbol),
      TRUE ~ gene_symbol),
      gene_symbol = str_remove(gene_symbol, "-Mar|-Sep"))
}

## Function for subsetting differential expression analysis results ------------
res_order_slice <- function(res, slice = TRUE, thres = 0.1) {
  res <- res %>%
    arrange(padj) %>%
    correct_gene_symbol()

  if(slice) {
    res <- res %>% 
      group_by(direction) %>%
      dplyr::slice(1:max(300, sum(padj < thres, na.rm = TRUE))) %>%
      ungroup
  }
  
  res <- res %>% dplyr::select(gene_symbol, log2FoldChange, direction)
  
  return(res)
}


## Function for permutation test comparing 2 sets of DEGs ----------------------
perm_test <- function(set1, set2, all_genes1, all_genes2, name1, name2, 
                      prefix1 = "", prefix2 = "", n_perm = 10000) {
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
    # mutate(p.perm = ifelse(p.perm > 0.5, 1 - p.perm, p.perm)) %>%
    as.data.frame() %>%
    pivot_wider(names_from = "direction.set1", values_from = "p.perm")
  
  n_bins <- min(max(n_overlap_perm$n_overlap_perm, na.rm = TRUE), 15)
  hist_overlap_perm <- n_overlap_perm %>%
    ggplot(aes(x = n_overlap_perm)) +
    geom_bar(stat = "count", size = 0.4) +
    geom_vline(data = n_overlap, aes(xintercept = n_overlap),
               color = "red", lty = 2) + 
    ggh4x::facet_grid2(direction.set2 ~ direction.set1, 
                       scales = "free", independent = "all") +
    theme_cowplot(font_size = 11)
  
  ggsave(paste0("figures/perm_test/perm_test_", prefix1, name1, "_", prefix2, name2, ".png"), 
         hist_overlap_perm, width = 5, height = 4, dpi = 400)
  return(list(n_overlap = n_overlap, p_overlap_perm = p_overlap_perm,
              hist_overlap_perm = hist_overlap_perm))
}


## Functions for drawing circos plots -------------------------------------

circo_overlap <- function(res1, res2, name1, name2, niceFacing = TRUE, 
                          circo_pal = NULL, gene_width_short = 2, nudge_x = rep(0, 4),
                          suffix = "", size, gene_width = 10, gene_cex = 0.6, num_cex = 0.75, p_cex = 0.85,
                          degree = -250, big_gap = 5, small_gap = 1, gene_list = "", step = 150,
                          phago_sig = NULL, show_selected_genes = FALSE, break_p_label = FALSE,
                          overlap_perm_test = NULL, conflict_facing = "outside") {
  res1 <- res1 %>% select(gene_symbol, direction)
  res2 <- res2 %>% select(gene_symbol, direction)
  
  
  names <- c(paste0(name2, c(".up", ".down")), paste0(name1, c(".up", ".down")))
  names <- factor(names, names)
  
  overlap_genes <- res1 %>%
    inner_join(res2, by = "gene_symbol", suffix = c(".1", ".2")) %>%
    mutate(direction.1 = ifelse(is.na(direction.1), NA, paste0(name1, ".", direction.1)),
           direction.2 = ifelse(is.na(direction.2), NA, paste0(name2, ".", direction.2))) %>%
    pivot_longer(cols = 2:3, names_to = "group", values_to = "from") %>%
    dplyr::rename("to" = "gene_symbol")%>%
    select(-group) %>% 
    mutate(from = str_replace(from, "2$", name2)) %>%
    mutate(value = 1, value2 = gene_width) %>% 
    select(from, to, value, value2)
  
  if(show_selected_genes) {
    overlap_genes <- overlap_genes %>%
      mutate(value2 = ifelse(to %in% c(gene_list, phago_sig), value2, gene_width_short))
  }
  
  count_ds <- rbind(
    dplyr::count(res1, direction) %>% mutate(direction = paste0(name1, ".", direction)),
    dplyr::count(res2, direction) %>% mutate(direction = paste0(name2, ".", direction))) %>%
    dplyr::rename("from" = "direction") %>% 
    full_join(dplyr::count(overlap_genes, from), by = "from", suffix = c(".total", ".overlap")) %>% 
    mutate(n.overlap = replace_na(n.overlap, 0)) %>%
    mutate(to = from, value = (n.total - n.overlap)/2, value2 = value) %>% 
    select(from, to, value, value2)
  
  group_genes <- overlap_genes %>%
    group_by(to) %>%
    summarise(group = paste(from, collapse = "_")) %>%
    mutate(group = str_remove_all(group, paste0(name1, "\\.|", name2, "\\."))) %>%
    arrange(group)
  
  
  group <- structure(c(rep("dirs", 4), as.character(group_genes$group)), 
                     names = c(as.character(names), group_genes$to)) %>% 
    factor(c(unique(as.character(group_genes$group)), "dirs"))
  
  group_mid <- split(group, group) %>%
    lapply(function(x) {x[ceiling(length(x)/2)]}) %>%
    lapply(names) 
  
  if(is.null(circo_pal)) {grid.col.ds <- setNames(c("#f7970c", "#4a4e4d", "brown1", "#0275d8"), names)}
  else {grid.col.ds <- unlist(circo_pal[as.character(names)])}
  grid.col <- c(grid.col.ds, setNames(rep("grey", nrow(group_genes)), group_genes$to))
  
  circo_df <- rbind(overlap_genes, count_ds) 
  if(show_selected_genes) {genes_display <- intersect(group_genes$to, c(gene_list, phago_sig))} 
  else { genes_display <- group_genes$to }
  
  par(bg="white")
  
  dev.copy(pdf, str_replace_all(paste0("figures/circos/circos_", name1, "_", name2, "_", suffix, ".pdf"),
                                "\\\n| ", "_") %>% str_replace_all("_+", "_") %>% str_remove_all(","), 
           height = size, width = size)
  
  circos.clear()
  circos.par(clock.wise = FALSE, start.degree = degree, message = FALSE)
  chordDiagram(circo_df, grid.col = grid.col,
               link.visible = circo_df[[1]] != circo_df[[2]],
               self.link = 2, 
               big.gap = big_gap, small.gap = small_gap,
               order = names(group),
               group = group,
               link.sort = TRUE, 
               transparency = 0.1,
               annotationTrack = c("grid"), 
               preAllocateTracks = list(track.height = max(strwidth(names(group)))))
  
  
  for(si in genes_display) {
    xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
    ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
    if(si %in% phago_sig) {col = "#00bf46"; face = 4}
    else if(si %in% intersect(subset(group_genes, group == "up_up")$to, gene_list)) {col = "brown1"; face = 4}
    else if(si %in% intersect(subset(group_genes, group == "down_down")$to, gene_list)) {col = "dodgerblue3"; face = 4}
    else if (si %in% gene_list) {col = "black"; face = 4}
    else {col = "grey50"; face = 3}
    
    circos.text(mean(xlim), ylim[1], si, sector.index = si, track.index = 1, 
                facing = "clockwise", adj = c(0, 0.5), col = col,
                cex = gene_cex, niceFacing = TRUE, font = face)
  }
  
  perm_test_res <- overlap_perm_test$p_overlap_perm %>%
    pivot_longer(2:3, names_to = "direction.set1", values_to = "p") %>%
    mutate(dir = paste(str_extract(direction.set1, c("up|down")),
                       str_extract(direction.set2, c("up|down")), sep = "_"))
  
  p_facing <- list(up_up = "outside", up_down = conflict_facing,
                   down_up = conflict_facing, down_down = "inside")
  nudge_x <- setNames(nudge_x, names)
  
  for(x in names(group_mid)[-length(group_mid)]) {
    si <- group_mid[[x]]
    xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
    ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
    
    p <- subset(perm_test_res, dir == x)$p
    label <- paste0("p = ", round(ifelse(p > 0.5, 1 - p, p), 4), 
                    ifelse(p > 0.5, " *", ""))
    if(p == 0 | p == 1) {
      label <- substitute(
        paste(bold(x)^bold(y), bold(z)), list(x = "p < 10", y = "-5", z = ifelse(p > 0.5, " *", ""))) %>% as.expression()}
    if(sum(group == x) < 4 & break_p_label) {
      label <- str_replace(label, "= ",  "=\n")
    }
    
    circos.text(mean(xlim), 0.55, label, sector.index = si, track.index = 1,
                facing = p_facing[[x]], cex = p_cex, niceFacing = FALSE, font = 2)
  }
  
  for(si in unique(circo_df$from)) {
    si_dir <- paste0(str_to_title(str_extract(si, "up|down")), " in")
    si_name <- ifelse(str_detect(si, name1), name1, name2)
    
    xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
    ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
    
    
    if(str_detect(si_name, "cKO|\\+|\\-")) {
      si_sup <- str_extract(si_name, "cKO|\\+|\\-")
      si_base <- str_remove(si_name, "cKO|\\+|\\-")
      y = c(0.37, 0.26)
      label_dir <- substitute(bold(dir), list(dir = si_dir)) 
      label <- c(label_dir, substitute(bolditalic(x)^bold(y), list(x = si_base, y = si_sup))) %>%
        as.expression()}
    else { y = 0.3; label <- paste0(si_dir, "\n", si_name) }
    circos.text(
      mean(xlim) + nudge_x[si], y, label,
      sector.index = si, track.index = 1, col = grid.col.ds[si],
      facing = "inside", font = 2,
      cex=num_cex, niceFacing = niceFacing
    )
    
    
    circos.axis(h = 0, major.at = seq(0, 3000, step),
                labels.cex = 0.5, labels.facing = "outside",
                sector.index = si, track.index = 1)
  }
  
  dev.off()
}

