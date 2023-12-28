library(Seurat)
library(tidyverse)
library(parallel)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(harmony)
library(cowplot)
library(ggpubfigs)
library(DoubletFinder)


# 1. Loading 10X files, creating and merging Seurat objects ------------------------------

nucseq <- sapply(paste0("Sample", 4:13), function(x) {
  sample <- Read10X(paste0("raw_data/", x))
  CreateSeuratObject(counts = sample, project = x)
}, simplify = FALSE, USE.NAMES = TRUE)

nucseq <- merge(x = nucseq[[1]], y = nucseq[-1],
                add.cell.ids = paste0("Sample", 4:13))
save(nucseq, file = "R_objects/nucseq_merged.RData")


# 2. QC --------------------------------
nucseq[["percent.mt"]] <- PercentageFeatureSet(nucseq, pattern = "^mt-")
nucseq <- subset( nucseq, subset = nCount_RNA>= 1000 & nFeature_RNA>= 200 & percent.mt < 80)

VlnPlot(nucseq, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# 3. Processing the snRNA-seq data (with Harmony batch correction) ------------------

npc <- 20; res <- 1
nucseq_harmony <- NormalizeData(nucseq, verbose = FALSE)
nucseq_harmony <- FindVariableFeatures(nucseq_harmony, selection.method = "vst", nfeatures = 2000)
nucseq_harmony <- ScaleData(nucseq_harmony, features = rownames(nucseq_harmony)) %>% 
  RunPCA(features = VariableFeatures(object = nucseq_harmony), npc = npc)

nucseq_harmony <- RunHarmony(nucseq_harmony, group.by.vars = "orig.ident", max.iter.harmony = 25)
nucseq_harmony <- RunUMAP(nucseq_harmony, reduction = "harmony", dims = 1:npc) 
nucseq_harmony <- FindNeighbors(nucseq_harmony, reduction = "harmony", dims = 1:npc) 
nucseq_harmony <- FindClusters(nucseq_harmony, resolution = res)

# 4. Add metadata --------------------------------------------
meta_tbl <- read.csv("data/nucseq_metadata.csv") %>%
  mutate(Genotype = factor(Genotype, c("control", "Tim3_cKO", "5XFAD", "Tim3_cKO.5XFAD")))

nucseq_harmony$Genotype <- plyr::mapvalues(
  x = nucseq_harmony$orig.ident, from = meta_tbl$SampleID, to = meta_tbl$Genotype)
nucseq_harmony$Batch <- plyr::mapvalues(
  x = nucseq_harmony$orig.ident, from = meta_tbl$SampleID,
  to = meta_tbl$Experiment)


save(nucseq_harmony, file = "R_objects/nucseq_harmony_dim20_res1.RData")

# 5. Identify microglia clusters ----------------------------
## ==> Cluster 6, 11, 32, 35, 41 are potentially microglia or PVM

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

dotplot_markers(nucseq_harmony)

# 6. Combine and subcluster clusters that are potentially microglia and PVM (6, 11, 32, 35, 41) ----------------------------
## (1) Subclustering ------------------------
nucseq_harmony_subcluster <- subset(nucseq_harmony, seurat_clusters %in% c(6, 11, 32, 35, 41))
nucseq_harmony_subcluster$orig.cluster <- nucseq_harmony$seurat_clusters

npc <- 20; res <- 1
nucseq_harmony_subcluster <- nucseq_harmony_subcluster(nucseq_harmony_subcluster, selection.method = "vst", nfeatures = 2000)
nucseq_harmony_subcluster <- ScaleData(nucseq_harmony, features = rownames(nucseq_harmony_subcluster)) %>% 
  RunPCA(features = VariableFeatures(object = nucseq_harmony_subcluster), npc = npc)

nucseq_harmony_subcluster <- RunHarmony(nucseq_harmony_subcluster, group.by.vars = "orig.ident", max.iter.harmony = 25)
nucseq_harmony_subcluster <- RunUMAP(nucseq_harmony_subcluster, reduction = "harmony", dims = 1:npc) 
nucseq_harmony_subcluster <- FindNeighbors(nucseq_harmony_subcluster, reduction = "harmony", dims = 1:npc) 
nucseq_harmony_subcluster <- FindClusters(nucseq_harmony_subcluster, resolution = res)

save(nucseq_harmony_subcluster, file = "R_objects/nucseq_markers_subcluster_harmony.RData")

## (2) Scoring microglia, PVM, BAM, monocyte signatures ------------------------
### Subluster 2, 3, 6, 18 => Microglia
### Sluster 12 => BAM-PVM
### Sluster 15 => BAM

### Signature list
### * Microglia: Known microglial genes, Van Hove et al. 2019, Monaghan et al. 2019, Zeisel et al. 2018
### * Border-associated microglia (BAM): Van Hove et al. 2019
### * Perivascular macrophages (PVMs): Zeisel et al. 2018, Yang et al. 2019
### * Monocyte/monocyte-derived antigen-presenting cells: Monaghan et al. 2019
### * Infiltrating monocyte: https://www.biorxiv.org/content/10.1101/2021.05.30.446342v1

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

parallel::mclapply(names(sig_list), function(x){
  nucseq_harmony_subcluster <- AddModuleScore(
    nucseq_harmony_subcluster, list(sig_list[[x]]), name = x)
  
  VlnPlot(nucseq_harmony_subcluster, features = paste0(x, 1), group.by = "seurat_clusters", 
          cols = c25, pt.size = 0.1) +
    labs(title = sig_title[x], x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  ggsave(paste0('figures/harmony_subcluster/vln_', x, "_harmony_subcluster.png"), width = 7, height = 3.2)
  ggsave(paste0('figures/harmony_subcluster/vln_', x, "_harmony_subcluster.pdf"), width = 7, height = 3.2)
})
