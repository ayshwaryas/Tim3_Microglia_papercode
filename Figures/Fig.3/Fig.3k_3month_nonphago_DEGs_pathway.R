library(tidyverse)
library(enrichR)
library(ggtext)

## DESeq2 results, 3 month old
load("data/DGE_results/2022-05-03.Dataset1_DGE_res_ordered.RData")

# EnrichR ---------

## KEGG pathways with disease pathways removed
## Disease pathways: pathways under section 6.1 to 6.10 with identifier number beginning with 049-054
## (https://www.kegg.jp/kegg/pathway.html)
KEGG_rm_disease <- read.table("data/KEGG_rm_disease.txt", sep = "\t", header = TRUE)

db <- "KEGG_2019_Mouse"
enrich <- function(res, rm_ribo = FALSE) {
  if(rm_ribo) {
    res <- res %>% filter(!str_detect(gene_name, "(Rp[ls]|Mrp[ls])"))
  }
  if(length(res$gene_name) <= 5) {enr_res <- NULL} 
  else {
    enr_res <- enrichr(genes = res$gene_name, databases = db)
    enr_res$KEGG_2019_Mouse <- enr_res$KEGG_2019_Mouse %>%
      filter(!Term %in% KEGG_rm_disease$Term) %>%
      filter(!str_detect(toupper(Term), "CANCER")) %>%
      filter(!str_detect(toupper(Term), "DISEASE")) %>%
      mutate(Adjusted.P.value = p.adjust(.$P.value, method = "BH"))
    enr_res <- enr_res %>%
      data.table::rbindlist(idcol = "database") %>%
      arrange(database, Adjusted.P.value)
  }
  
  return(enr_res)
}

## DEGs (FDR < 0.1), Havcr2cKO non-phagocytosing vs control non-phagocytosing microglia from 3 month old mice
tim3_sig_nonphago <- results_batch1_ordered$`control vs Tim3.cKO in phago-` %>%
  mutate(direction = ifelse(direction == "up", "down", "up")) %>%
  group_by(direction) %>% filter(padj < 0.1) %>%
  select(gene_name, padj, direction) %>% 
  split(f = .$direction) 

## Pathway enrichment analysis on DEGs (removing ribosomal genes) using enrichR
enriched_tim3_nonphago_sig.1_no_ribo <- tim3_sig_nonphago %>%
  lapply(enrich, rm_ribo = TRUE) %>% ## remove rimosomal genes
  data.table::rbindlist(idcol = "direction")

# Plot pathways -------------------------------
pal <- colorRampPalette(brewer.pal(9, "Reds"))(20)[c(20, 12, 6, 3)]
plotPathway <- function(enr_res, db = "MSigDB_Hallmark_2020", width = 5, legend_pos = "right",
                        title = "", title_prefix = "", subfolder = "", fig_num = "3k") {
  df <- enr_res %>%
    filter(database == db) %>%
    filter(Adjusted.P.value < 0.1) %>%
    mutate(direction = paste(str_to_title(direction), "in <i>Havcr2</i>^cKO")) %>%
    mutate(text_col = ifelse(Adjusted.P.value < 0.001, "white", "black")) %>%
    mutate(fill_col = case_when(Adjusted.P.value < 0.001 ~ "< 0.001",
                                Adjusted.P.value < 0.01  ~ "< 0.01",
                                Adjusted.P.value < 0.05  ~ "< 0.05",
                                Adjusted.P.value < 0.1   ~ "< 0.1")) %>%
    mutate(fill_col = factor(fill_col, paste0('< ', c(0.001, 0.01, 0.05, 0.1))))
  
  
  if(nrow(df) == 0) {return(NULL)}
  
  p <- df %>%
    mutate(direction = fct_rev(direction)) %>%
    arrange(direction, Adjusted.P.value) %>%
    ungroup() %>%  
    mutate(index = n()-row_number() + 1) %>%  
    ggplot(aes(x = -log10(Adjusted.P.value), y = factor(index))) +
    geom_vline(xintercept = 1:floor(max(-log10(df$Adjusted.P.value))), alpha = 0.4,
               color = "grey30", lwd = 0.2, lty = "21") +
    geom_col(aes(fill = fill_col)) +
    geom_text(aes(x = max(-log10(df$Adjusted.P.value)) * 0.02,
                  label = Term, color = text_col), size = 3.5, hjust = 0) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
    scale_color_identity() +
    scale_fill_manual(breaks = c("< 0.001", "< 0.01", "< 0.05", "< 0.1"), 
                      values = pal, name = "padj", drop = FALSE) +
    labs(x = expression(-log[10](padj)), y = NULL, 
         title = paste0(title_prefix, str_remove_all(db, "_[0-9]+|_Mouse"), " Pathways"))+
    guides(fill = guide_legend(keywidth = 0.65, keyheight = 0.65)) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
          axis.title.x = element_text(size = 11),
          axis.text.y = element_blank(),
          strip.text.y = element_markdown(),
          axis.ticks.x = element_line(size = 0.25),
          axis.ticks.length.x = unit(0.05, "cm"),
          axis.ticks.y = element_blank(),
          legend.position = legend_pos) + 
    facet_grid(direction ~., scale= "free_y", space = "free_y")
  
  
  height <- case_when(nrow(df) == 1 ~ 1.2,
                      nrow(df) == 2 ~ 1.5,
                      nrow(df) == 3 ~ 1.8,
                      nrow(df) == 4 ~ 1.9,
                      nrow(df) == 5 ~ 2,
                      nrow(df) <= 10 ~ nrow(df) * 0.35,
                      TRUE ~ nrow(df) * 0.3)
  
  filename <- paste0("figures/Fig.", fig_num, "_pathway", subfolder,
                     "/Fig.", fig_num, "_enr_", 
                     str_replace_all(title, ": |, | ", "_"), "_", db)
  ggsave(paste0(filename, ".png"), p, width = width, height = height, dpi = 400)
  ggsave(paste0(filename, ".pdf"), p, width = width, height = height)
}



plotPathway(enriched_tim3_nonphago_sig.1_no_ribo , db = db,
            title = "tim3_nonphago_sig_padj.1_no_ribo")

