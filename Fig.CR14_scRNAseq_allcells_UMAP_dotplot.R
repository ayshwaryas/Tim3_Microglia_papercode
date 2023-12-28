library(Seurat)
library(tidyverse)
library(ggtext)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(cowplot)

setwd("/broad/kuchroolab/kimi_microglia/nucseq_20220216")
load("R_objects/nucseq_harmony_rm_doublets_annotated.RData")

## DimPlot --------------
annot_pal <- c("plum2", "deeppink1", "maroon4", 
               "grey75", "grey20", 
               "#E31A1C", "#770000", "#FEE5D9", 
               "#9ECAE1", "dodgerblue2", "#6A3D9A", "gold1", "pink2", "orange",  
               "yellow2", "yellow4", "darkorange2", "darkgoldenrod4", "darkturquoise", "springgreen", 
               "cornsilk2", "slategray2", "darkseagreen1", "aquamarine1", "green3", "darkcyan")

DimPlot(nucseq_harmony_rm_doublets, group.by = "annot_new",
        label = TRUE, label.size = 4, repel = TRUE, 
        cols = annot_pal, pt.size = 0.1) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3))) +
  labs(title = NULL)

ggsave("figures/dimplot_final_annot.png", width = 10, height = 7)
ggsave("figures/dimplot_final_annot.pdf", width = 10, height = 7)

## Dotplot of Havcr2 and Cx3cr1 ----------------------
nucseq_harmony_rm_doublets$Genotype_labels <- as.character(nucseq_harmony_rm_doublets$Genotype)
nucseq_harmony_rm_doublets$Genotype_labels[nucseq_harmony_rm_doublets$Genotype == "Tim3_cKO"] <- "<i>Havcr2</i><sup>icKO</sup>"
nucseq_harmony_rm_doublets$Genotype_labels[nucseq_harmony_rm_doublets$Genotype == "Tim3_cKO.5XFAD"] <- "<i>Havcr2</i><sup>icKO</sup> 5XFAD"

nucseq_harmony_rm_doublets$group <- paste(as.character(nucseq_harmony_rm_doublets$annot_new), 
                                          nucseq_harmony_rm_doublets$Genotype_labels, sep = "__")


parallel::mclapply(
  c("Havcr2", "Cx3cr1"),
  function(x) {
    p <- DotPlot(nucseq_harmony_rm_doublets, features = x, group.by = "group", scale = FALSE) +
      scale_color_gradientn(colors = rev(brewer.pal(11, "RdBu")))
    
    p$data$Genotype <- str_extract(p$data$id, "(?<=__).*") %>%
      factor(c("control", "<i>Havcr2</i><sup>icKO</sup>", 
               "5XFAD", "<i>Havcr2</i><sup>icKO</sup> 5XFAD"))
    p$data$id <- str_extract(p$data$id, ".*(?=__)") %>% factor(annot_levs)
    p + facet_grid(~Genotype) + 
      scale_y_discrete(limits = rev) +
      labs(x = NULL, y = NULL, title = x) +
      guides(color = guide_colorbar(order = 1, title = "Average Expression"),
             size = guide_legend(order = 2, title = "Percent Expression")) +
      theme_cowplot(font_size = 12) +
      theme(strip.background = element_rect(fill = NA),
            strip.clip = "off",
            plot.title = element_text(hjust = 0.5),
            strip.text = element_markdown(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.spacing = unit(0, "mm"))
    ggsave(paste0("figures/dotplot_final_annot_rm_doublets_", x, ".png"),
           width = 6, height = 5.5)
    ggsave(paste0("figures/dotplot_final_annot_rm_doublets_", x, ".pdf"),
           width = 6, height = 5.5)
  }, mc.cores = 10)
