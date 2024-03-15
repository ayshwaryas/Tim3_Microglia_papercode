library(tidyverse)
library(Seurat)
library(harmony)
library(parallel)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(cowplot)
library(ggpubfigs)
library(DoubletFinder)
library(ggtext)


# palette -------------------------
c25 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00","black", "gold1", 
         "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", "gray70", "khaki2", 
         "maroon", "orchid1","deeppink1","blue1","steelblue4","darkturquoise", "green3", 
         "yellow4", "yellow2", "#FEC44F", "brown4","darkcyan", "darkgoldenrod2", 
         "darkolivegreen4", "cornsilk", "aquamarine1", "darkseagreen1", "honeydew3", "salmon", 
         "slateblue3", "thistle2", "slategray2", "springgreen", "rosybrown1", "plum2", 
         "chocolate", "blanchedalmond", "brown1", "aliceblue", "cadetblue1", "coral1", 
         "firebrick1", "darkslategrey", "darkslateblue", "#745745", "#a01b68", "#ffc372",
         "#fff0f0", "#ebd4d4", "#eeeeee", "#de4463", "#e8ffc1", "#19d3da", "#ee6f57", "#ff9642",
         "#646464", "#0072ce", "#dc9cbf", "#ED8E7C", "#F1ECC3", "#D9DD6B", "#7C83FD", "#FDD2BF",
         "#E98580", "#C6B4CE", "#FFDADA", "#B5EAEA", "#BB8760", "#C9E4C5", "#FAEBE0", "#F7DBF0",
         "#66DE93", "#FFC074", "#B2B8A3", "#F2F4C3", "#D1D9D9", "#A7D0CD", "#114E60", "#28B5B5")


# Genotype levels and labels to use in the figures ----------------------
genotype_breaks <- c("control", "Tim3_cKO", "5XFAD", "Tim3_cKO.5XFAD")
genotype_labels <- c("control", "<i>Havcr2</i><sup>icKO</sup>", "5XFAD", "<i>Havcr2</i><sup>icKO</sup> 5XFAD")

# Function for processing snRNA-seq data ----------------------------
scrna_process <- function(seurat_obj, npc, res, normalize = TRUE, 
                          run_harmony = FALSE, max.iter.harmony = 25) {
  if(normalize) {
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  }
  
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj)) %>% 
    RunPCA(features = VariableFeatures(object = seurat_obj), npc = npc) 
  
  if(run_harmony) {
    seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "orig.ident", 
                             max.iter.harmony = max.iter.harmony)
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


# Signature list ---------------------------------------------
## * Microglia: Known microglial genes, Van Hove et al. 2019, Monaghan et al. 2019, Zeisel et al. 2018
## * Border-associated macrophage (BAM): Van Hove et al. 2019
## * Perivascular macrophages (PVMs): Zeisel et al. 2018, Yang et al. 2019
## * Monocyte/monocyte-derived antigen-presenting cells: Monaghan et al. 2019
## * Infiltrating monocyte: https://www.biorxiv.org/content/10.1101/2021.05.30.446342v1

sig_list <- list(
  MG_sanity_check = c("Havcr2", "Tmem119", "Hexb", "P2ry12", "Sall1", "Mafb", "Spi1",
                      "Cx3cr1", "Trem2", "C1qc", "Cst7", "Apoe", "Axl", "Itgax", "Lilrb4",
                      "Ly9", "Il1b", "Tnf", "Lag3", "Cd274", "Pdcd1"),
  Microglia_Monaghan_2019 = c("Tmem119", "Sall1", "Hexb", "Slc2a5", "P2ry12",
                              "Siglech", "Trem2", "Bhlhe41", "Gpr34", "Serpine2"),
  Microglia_MouseBrain = c("P2ry12", "Ccr6", "Tmem119", "Casp8", "Ccl4", "Tnf"),
  Microglia_VanHove_2019 = c("Sall1", "Adgrg1", "Il1a", "Golm1", "Adamts1", "Selplg", "Sparc", 
                             "Gem", "Serpine2", "P2ry12", "Fscn1", "Rtn4rl1", "Atp6v0a2", 
                             "Arhgap5", "Siglech", "Ctsl", "Smad7", "P2ry13", "Tmem119", 
                             "Fam102b", "Smap2", "Csf3r", "Soga1", "Kcnk6", "Wsb1", "Rgs1", "Ldhb", 
                             "Gpr34", "Glul", "Cst3", "Lhfpl2", "Epb41l2", "Gpr155", "Rasal3", "Dtx4", 
                             "Srgap2", "Frmd4a", "Lpcat2", "Spata13", "Qk", "Entpd1", "Hexb", "Vsir", 
                             "Tanc2", "Bin1", "Trem2", "Olfml3", "Ctsf", "Elmo1", "Asph", "Gna15", "Cmtm6", 
                             "Herpud1", "Brd2", "Pde3b", "Tgfbr1", "Daglb", "Nfkbia", "Tmem173", "Basp1", "Arsb", 
                             "Inpp5d", "Fmnl3", "Slc29a3", "Sipa1", "Cx3cr1", "Snx18", "F11r", "Itgb5", "Fcrls"),
  BAM_signature = c("C3ar1", "Nrp1", "Cd63", "Ms4a7", "Adrb2", "Ms4a6b", "Ms4a6c", "Ifnar1",
                    "Ptger4", "Lifr", "Clec12a", "Zfp36l2", "Ehd4", "Zfp36l1", "Myo5a", 
                    "Swap70", "Rasa4", "B3galnt1", "St8sia4", "H1f0", "Sesn1", "Apoe", 
                    "Cd14", "Sdc4", "Aoah", "Pla2g7", "Tgfbi", "Ier2", "Maf"),
  PVM = c("Cd163", "Mrc1", "Lyve1", "Pf4", "Cd74", "Cxcl2"),
  Monocyte_and_derived_APCs = c("Ly6c2", "Ccr2", "Cd44", "Fcgr1", "Plac8", "Nr4a1"),
  Infiltrating_monocyte = c("Nupr1", "Adgre4", "Plac8", "Ifitm6", "Smpdl3a",
                            "Vim", "Gpx1", "Napsa", "S100a6", "S100a4", "S100a11",
                            "Crip1", "Ear2", "Fabp4", "Chil3", "Ccr2")
)

sig_title <- setNames(c("Microglia Sanity Check", "Microglia (Monaghan et al. 2019)",
                        "Microglia (Mouse Brain)", "Microglia (Van Hove et al. 2019)"),
                      "BAM Signature", "PVM Signature",
                      "Monocyte/Monocyte-Derived APCs", "Infiltrating Monocyte",
                      names(sig_list))

# Function for adding p-values to signature score violin plot -------------------------------
add_p_value <- function(p, var, y_pos = c(0.78, 0.9, 0.68), ylim = 0, 
                        levs = c("5XFAD", "Tim3_cKO.5XFAD_Top", "Tim3_cKO.5XFAD_Btm"),
                        x_min = c(1, 1, 1.18), x_max = c(1.18, 1.36, 1.36)) {
  stat.test <- p$data %>% 
    filter(str_detect(split, "5XFAD")) %>%
    mutate(split = factor(split, levs)) %>%
    t_test(as.formula(paste0(var, "1 ~ split"))) %>%
    add_xy_position() %>%
    add_significance("p") %>%
    mutate(p.signif = ifelse(p >= 0.05, paste("p =", round(p, 3)), p.signif))  %>%
    mutate(size = ifelse(p < 0.05, 5, 2),
           vjust = ifelse(p < 0.05, 0.5, -0.1),
           x_min = x_min, x_max = x_max,
           y_pos = y_pos) 
  
  p + stat_pvalue_manual(
    stat.test, fontface = 2, 
    xmin = "x_min", xmax = "x_max",
    y.position = "y_pos",
    bracket.size = 0.5, label = "p.signif", 
    tip.length = 0.015, vjust = "vjust", inherit.aes = FALSE) +
    scale_y_continuous(limits = c(NA, max(c(y_pos, ylim))))
}

# Function to draw dotplot of mouse brain markers ----------------------------------------
## (1) Top 10 markers of human brain dataset (https://www.biorxiv.org/content/10.1101/2021.05.30.446342v1)
human_brain_markers <- readRDS("data/markers/2019-08-09.exp4_res1_markers.RDS") %>%
  as.data.frame() %>% filter(avg_logFC > 0.25) 

human_brain_clusters <- read.table("data/markers/clusters.txt", 
                                   sep = "\t", header = TRUE)  %>%
  mutate(cluster = factor(cluster, 0:max(cluster))) %>%
  mutate(celltype = case_when(cluster == 13 ~ "Excitatory-1",
                              cluster == 15 ~ "Excitatory-2",
                              TRUE ~ celltype)) %>%
  arrange(celltype) %>%
  mutate(celltype_index = 1:n())

## (2) Marker genes from Zeisel et. al. 2015 (PMID: 25700174)
Zeisel <- read.csv("data/markers/Zeisel_2015_marker.csv") %>%
  mutate(Genes = strsplit(Genes, ", "))%>%
  tidyr::unnest(Genes)

## (3) Marker genes from Di Bella et. al. 2021 (PMID: 34163074)
Di_Bella <- read.csv("data/markers/Di_Bella_2021_marker.csv") %>%
  mutate(celltype = subtype_abbr) %>%
  mutate(Genes = strsplit(Genes, ", ")) %>%
  tidyr::unnest(Genes)

human_brain_top10 <- human_brain_markers %>% 
  filter(p_val_adj < 0.1) %>%
  inner_join(data.frame(gene_lower = rownames(nucseq_harmony),
                        gene_upper = toupper(rownames(nucseq_harmony))), by = c("gene" = "gene_upper")) %>%
  left_join(human_brain_clusters, by = "cluster") %>%
  dplyr::group_by(celltype) %>% 
  arrange(desc(avg_logFC)) %>%
  slice(1:10) %>% mutate(index = 1:n()) %>% ungroup %>%
  select(cluster, celltype_index, celltype, gene_lower, index)

## Function to draw dotplot 
draw_dotplot <- function(df, group = FALSE, type = "", suffix = "", width = 15, height = 12,
                         rot = FALSE) {
  p <- ggplot(df, aes(x = id, y = features.plot)) +
    geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
    scale_color_gradient(low = "grey", high = "blue") +
    scale_size(range = c(0, 4), breaks = seq(0, 100, 25),
               limits = c(0, 100 - 1E-6)) + 
    cowplot::theme_cowplot(font_size = 10) +
    guides(color = guide_colorbar(order = 1),
           size = guide_legend(order = 2)) + 
    labs(x = NULL, y = NULL) +
    theme(strip.text = element_text(size = 7),
          panel.spacing = unit(1, "mm"))
  
  if(rot) { p <- p + coord_flip() + 
    facet_grid( ~ celltype, switch = "x", scale = "free_x", space = "free_x") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_x_discrete(limits = rev)
  }
  else {p <- p + facet_grid(celltype ~ ., switch = "y", scale = "free_y", space = "free_y")}
  ggsave(paste0("figures/DotPlot_harmony/dotplot_", type, 
                ifelse(group, paste0(" - ", df$group[1]), ""), suffix, ".png"), p,
         width = width, height = height)
}

## Dotplot of marker genes
dotplot_markers <- function(SeuratData, suffix = "", width = 15, height = 12, rot = FALSE) {
  
  
  ## (1) Top 10 markers of human brain dataset (https://www.biorxiv.org/content/10.1101/2021.05.30.446342v1)
  dotplot_human_brain_top10 <- DotPlot(SeuratData, features = unique(human_brain_top10$gene_lower))
  dotplot_human_brain_top10$data %>%
    full_join(human_brain_top10, by = c("features.plot" = "gene_lower")) %>%
    mutate(group = case_when(celltype_index <= 9 ~ 1, 
                             str_detect(celltype, "^L[0-9]") ~ 2, 
                             TRUE ~ 3)) %>%
    split(f = .$group) %>% 
    lapply(draw_dotplot, type = "human_brain", group = TRUE,
           suffix = suffix, width = width, height = height, rot = rot)
  
  ## (2) Marker genes from Zeisel et. al. 2015 (PMID: 25700174)
  dotplot_Zeisel <- DotPlot(SeuratData, features = unique(Zeisel$Genes))
  dotplot_Zeisel$data %>%
    full_join(Zeisel, by = c("features.plot" = "Genes")) %>%
    draw_dotplot(type = "Zeisel_2015", suffix = suffix,  width = width)
  
  ## (3) Marker genes from Di Bella et. al. 2021 (PMID: 34163074)
  dotplot_Di_Bella <- DotPlot(SeuratData, features = unique(Di_Bella$Genes))
  dotplot_Di_Bella$data %>%
    full_join(Di_Bella %>% mutate(subtype_abbr = factor(subtype_abbr, unique(subtype_abbr))), 
              by = c("features.plot" = "Genes"))  %>%
    draw_dotplot(type = "Di_Bella_2021", suffix = suffix, width = width, height = height, rot = rot)
}


# Function for drawing signature score violin plot ----------------------------------
VlnPlot_custom <- function(
    SeuratData, feature_name, cluster = NULL, fig_num = '6',
    group.var = "new_clusters", split.var = c("bimod_genotype", "Genotype"), title = "", log = FALSE,
    width = 10, height = 3, suffix = "", add_p_value = FALSE, 
    x_min = c(1, 1, 1.18), x_max = c(1.18, 1.36, 1.36), add1 = TRUE,
    y_pos = c(0.78, 0.9, 0.68), ylim = 0, pal_custom = NULL, hjust = 0.5, angle = 0,
    folder = "Fig.6_nucseq_figures", levs = c("5XFAD", "Tim3_cKO.5XFAD_Top", "Tim3_cKO.5XFAD_Btm"))  {
  
  ## Palette
  pal <- c("dodgerblue3", "#EF3B2C", friendly_pal("nickel_five", 5)[3],
           friendly_pal("ito_seven", 7)[c(6, 4:3)])
  
  if(!is.null(cluster)) {
    SeuratData <- subset(SeuratData, new_clusters == cluster)
  }
  
  if(split.var == "bimod_genotype") {
    breaks <- levels(SeuratData$bimod_genotype)
    labels <- c("control", "<i>Havcr2</i><sup>icKO</sup>", "5XFAD",
                "<i>Havcr2</i><sup>icKO</sup> 5XFAD",
                "<i>Havcr2</i><sup>icKO</sup> 5XFAD;P1",
                "<i>Havcr2</i><sup>icKO</sup> 5XFAD;P2")
    cols <- pal
    
    if(!is.null(cluster))  {
      if(cluster == 2) {
        breaks <- breaks[-4]
        labels <- labels[-4]
        cols <- cols[-4]
      }
    }
  }
  else if(split.var == "Genotype") {
    breaks <- c("control", "Tim3_cKO", "5XFAD", "Tim3_cKO.5XFAD")
    labels <- c("control", "<i>Havcr2</i><sup>icKO</sup>",
                "5XFAD", "<i>Havcr2</i><sup>icKO</sup> 5XFAD")
    cols <- pal[1:4]
  }
  else {
    breaks <- levels(pull(SeuratData@meta.data, split.var))
    labels <- breaks
    cols <- pal_custom
  }
  
  suffix <- paste0(feature_name, ifelse(split.var == "bimod_genotype","_bimod", ""), suffix)
  
  feature <- paste0(feature_name, ifelse(add1, "1", ""))
  p <- VlnPlot(SeuratData, features = feature, log = log,
               split.by = split.var, group.by = group.var) +
    labs(x = NULL, y = NULL, title = title) +
    theme_cowplot(font_size = 12) +
    scale_fill_manual(breaks = breaks, values = cols, labels = labels, limits = force) +
    theme(legend.text = element_markdown(size = 12),
          legend.position = "right",
          plot.title = element_markdown(hjust = 0.5),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12, angle = angle, hjust= hjust))
  
  if(add_p_value) {
    p <- add_p_value(p, var = feature_name,  y_pos = y_pos, ylim = ylim, 
                     levs = levs, x_min = x_min, x_max = x_max)
  }
  
  if(!is.null(cluster)) {
    suffix <- paste0(suffix, "_sub", cluster)
    
  }
  filename <- paste0("figures/", folder, "/Fig.", fig_num, "_vln_", suffix)
  ggsave(paste0(filename, ".png"), p,  width = width, height = height, dpi = 400)
  ggsave(paste0(filename, ".pdf"), p, width = width, height = height)
  
}