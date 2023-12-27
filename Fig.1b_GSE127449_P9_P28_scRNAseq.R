library(Seurat)
library(tidyverse)
library(parallel)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(cowplot)
library(harmony)
library(parallel)

## ggplot2 theme 
theme_custom <- theme_cowplot(font_size = 11) +
  theme(axis.text.y = element_text(face = 3),
        legend.title = element_markdown(size = 12.5),
        legend.box = "horizontal")
theme_set(theme_custom)

## color palettes
c25 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00","black", "gold1", 
         "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", "gray70", "khaki2", 
         "maroon", "orchid1","deeppink1","blue1","steelblue4","darkturquoise", "green1", 
         "yellow4", "yellow3", "darkorange4", "brown","darkcyan", "darkgoldenrod1", 
         "darkolivegreen4", "cornsilk", "aquamarine1", "darkseagreen1", "honeydew3", "salmon", 
         "slateblue3", "thistle2", "slategray2", "springgreen", "rosybrown1", "plum2", 
         "chocolate", "blanchedalmond", "brown1", "aliceblue", "cadetblue1", "coral1", 
         "firebrick1", "darkslategrey", "darkslateblue", "#745745", "#a01b68", "#ffc372",
         "#fff0f0", "#ebd4d4", "#eeeeee", "#de4463", "#e8ffc1", "#19d3da", "#ee6f57", "#ff9642",
         "#646464", "#0072ce", "#dc9cbf", "#ED8E7C", "#F1ECC3", "#D9DD6B", "#7C83FD", "#FDD2BF",
         "#E98580", "#C6B4CE", "#FFDADA", "#B5EAEA", "#BB8760", "#C9E4C5", "#FAEBE0", "#F7DBF0",
         "#66DE93", "#FFC074", "#B2B8A3", "#F2F4C3", "#D1D9D9", "#A7D0CD", "#114E60", "#28B5B5")

pal_c <- c(colorRampPalette(c("#011233", "blue4", dichromat::colorschemes$DarkRedtoBlue.18[1]))(4), 
           colorRampPalette(c(dichromat::colorschemes$DarkRedtoBlue.18[2:5]))(5), 
           colorRampPalette(dichromat::colorschemes$DarkRedtoBlue.18[6:10])(6),
           colorRampPalette(brewer.pal(9, "YlOrRd")[1:4])(5),
           colorRampPalette(brewer.pal(9, "YlOrRd")[-(1:4)])(6))

## Border-associated microglia (BAM) signature
BAM_sig <- c("C3ar1", "Nrp1", "Cd63", "Ms4a7", "Adrb2", "Ms4a6b", "Ms4a6c",
             "Ifnar1", "Ptger4", "Lifr", "Clec12a", "Zfp36l2", "Ehd4", "Zfp36l1",
             "Myo5a", "Swap70", "Rasa4", "B3galnt1", "St8sia4", "H1f0", "Sesn1", 
             "Apoe", "Cd14", "Sdc4", "Aoah", "Pla2g7", "Tgfbi", "Ier2", "Maf")

## Immune checkpoints, and TGFb pathway-related molecules
checkpoint_tgfb <-  list(
  Checkpoint = c("Havcr2", "Lag3", "Vsir", "Pdcd1", "Cd274", "Ctla4"),
  Tgfb = c("Tgfb1", "Tgfbr1", "Tgfbr2"),
  Smad = c("Smad2", "Smad3", "Smad4"))

## Function for snRNA-seq processing
scrna_process <- function(seurat_obj, npc, res, normalize = TRUE, 
                          run_harmony = FALSE, group_var = "orig.ident") {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  seurat_obj <- subset(seurat_obj, subset = nCount_RNA>= 1000 & nFeature_RNA>= 200 & percent.mt < 80)
  
  if(normalize) {
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  }
  
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj)) %>% 
    RunPCA(features = VariableFeatures(object = seurat_obj), npc = npc) 
  
  if(run_harmony) {
    seurat_obj <- RunHarmony(seurat_obj, group.by.vars = group_var,
                             max.iter.harmony = 25)
    seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:npc) 
    seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:npc) 
    seurat_obj <- FindClusters(seurat_obj, resolution = res)
  }
  else {
    seurat_obj <- seurat_obj %>%
      FindNeighbors(dims = 1:npc) %>% 
      FindClusters(resolution = res) %>% 
      RunUMAP(dims = 1:npc)
  }
  return(seurat_obj)
}


# Load Data --------------
## Metadata
GSE127449_meta <- read.table("data/GSE127449_P9_P28/GSE127449/GSE127449_microglia.barcode.tsv",
                             sep = "\t", header = FALSE, col.names = "barcode") %>%
  mutate(dataset = "GSE127449") %>%
  mutate(age = str_extract(toupper(barcode), "P9|P28")) %>%
  mutate(age = factor(age, c("P9", "P28")))
  mutate(tissue = "cortex") %>%
  mutate(genotype = str_extract(barcode, "ctx|het|ko")) %>%
  mutate(genotype = case_when(genotype == "ctx"~ "WT",
                              genotype == "het" ~ "ST2_het",
                              genotype == "KO" ~ "ST2_KO")) %>%
  mutate(strain = "B6") %>% mutate(sex = "F") %>%
  tibble::column_to_rownames("barcode")

## Count matrix
GSE127449_count <- Read10X("data/GSE127449_P9_P28/GSE127449")

## Create Seurat object
GSE127449_seurat <- CreateSeuratObject(
  counts = GSE127449_count, assay = "RNA", meta.data = GSE127449_meta) %>%
  subset(genotype == "WT")

# Processing the scRNA-seq dataset --------------
GSE127449_seurat_processed <- scrna_process(
  GSE127449_seurat, npc = 50, res = 1, normalize = TRUE, run_harmony = FALSE)

## save(GSE127449_seurat_processed, file = "data/GSE127449_P9_P28/GSE127449_seurat_processed.RData")

# Identify & exclude BAM cluster --------------
## (Cluster 9 and 12 are BAM clusters)
GSE127449_seurat <- AddModuleScore(GSE127449_seurat, list(BAM_sig), name = "BAM_signature")

FeaturePlot(GSE127449_seurat, features = "BAM_signature1", 
            label = TRUE, repel = TRUE, label.size = 3) +
  scale_colour_gradientn(colours = rev( brewer.pal(11, "Spectral")))
ggsave("figures/Fig.1b_GSE127449_P9_P28/featureplot_BAM_signature.png",
       width = 7, height = 5, dpi = 400)

VlnPlot(GSE127449_seurat, features = "BAM_signature1",
        cols = c25, pt.size = 0.2) + NoLegend()
ggsave(paste0("figures/Fig.1b_GSE127449_P9_P28/vlnplot_BAM_signature.png"),
       width = 7, height = 3, dpi = 400)


# Generate DotPlot (Fig. 1b) ----------------
## Rmove BAM clusters
BAM_clusters_P9_P28 <- c(9, 12)
MG_clusters_P9_P28 <- setdiff(unique(GSE127449_seurat_processed$seurat_clusters), 
                              c(BAM_clusters_P9_P28))

GSE127449_seurat_MG <- subset(GSE127449_seurat_processed, 
                              seurat_clusters %in% MG_clusters_P9_P28)

## Dotplot of immune checkpoints, and TGFb pathway-related molecules

p_ls <- lapply(names(checkpoint_tgfb), function(type) {
  if(type == "Smad") {scale_breaks = seq(0, 12, 4)}
  else {scale_breaks = seq(0, 75, 25)}
  p <- DotPlot(GSE127449_seurat_MG, features = checkpoint_tgfb[[type]], scale = FALSE, 
          group.by = "age", cols = c("blue", "red"), scale.by = "radius") + 
    coord_flip() + 
    labs(x = NULL, y = NULL, radius = "Average Expression") + 
    scale_x_discrete(limits = rev) +
    scale_color_gradientn(name = "Average Expression",colors = pal_c, trans = "log10", na.value = "grey80") +
    scale_radius(range = c(2, 9), name = "Pct Expr", limits = c(1E-99, NA), breaks = scale_breaks) +
    guides(color = guide_colorbar(order = 1, barheight = 4.5, barwidth = 0.7), 
           radius = guide_legend(order = 2, title = "Pct Expr"),
           size = guide_legend(title = "Pct Expr"))
  if(type != "Smad") {
    p <- p + theme(
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(), 
      axis.line.x = element_blank(),
      plot.margin = margin(5, 5, 0, 5))
  }
  p
})

g <- egg::ggarrange(
  p_ls[[1]], 
  p_ls[[2]] + theme(legend.box.margin = margin(t = -10)), 
  p_ls[[3]] + theme(legend.box.margin = margin(t = 10)), 
  ncol = 1, heights = c(1, 0.6, 0.6))

filename <- paste0("figures/Fig.1b_GSE127449_P9_P28/Fig.1b_GSE127449_P9_P28_dotplot_")
ggsave(paste0(filename, ".png"), g, width = 4.2, height = 4.5, dpi = 400, units = "in")
ggsave(paste0(filename, ".pdf"), g, width = 4.2, height = 4.5)


