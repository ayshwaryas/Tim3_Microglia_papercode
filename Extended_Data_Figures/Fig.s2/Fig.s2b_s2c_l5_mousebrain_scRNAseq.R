library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(scales)
library(parallel)
library(ggtext)

load("data/mouse_brain/l5_seurat.RData")
load("data/mouse_brain/l5_PVM_MGL.RData")

## Color Palette
pal <- c(brewer.pal(9, "Reds")[c(3, 6, 9)], brewer.pal(6, "Blues")[c(4, 6)], 
         brewer.pal(8, "Pastel2"), brewer.pal(9, "Pastel1")) 

## Figure theme
theme_custom <- cowplot::theme_cowplot(font_size = 11) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank())
theme_set(theme_custom)

## Fucntion for processing scRNA-seq
scrna_process <- function(seurat_obj, npc, res, normalize = TRUE, 
                          run_harmony = FALSE, harmony = TRUE) {
  if(normalize) {
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  }
  
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj)) %>% 
    RunPCA(features = VariableFeatures(object = seurat_obj), npc = npc) 
  
  if(run_harmony) {
    seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "orig.ident")
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


# Load scRNA-seq dataset from Zeisel et al 2018 (PMID: 30096314) ---------------------
## Read the loom file and convert it to a Seurat object
l5.immune <- SeuratDisk::Connect(filename = "data/mouse_brain/l5_all.loom", mode = "r")
l5.seurat <- as.Seurat(l5.immune)
save(l5.seurat, file = "data/mouse_brain/l5_seurat_raw.RData")
l5.immune$close_all()

# Processing the scRNA-seq dataset -----------

## QC
l5.seurat[["percent.mt"]] <- PercentageFeatureSet(l5.seurat, pattern = "^mt-")
l5.seurat <- subset(l5.seurat, subset = nCount_RNA>= 1000 & nFeature_RNA>= 200 & percent.mt < 80)

## Normalization, scaling, dimensional reduction and clustering (nPC = 50, resolution = 0.1)
l5.seurat <- scrna_process(l5.seurat, npc = 50, res = 0.1, normalize = TRUE)
save(l5.seurat, file = "data/mouse_brain/l5_seurat.RData")

## Relevel clusters 
ClusterName_MGL_PVM <- c("MGL1", "MGL2", "MGL3", "PVM1", "PVM2")
celltype_labels <- FetchData(l5.seurat, c("Subclass", "ClusterName", "TaxonomyRank3", "TaxonomyRank4")) %>%
  mutate(label = ifelse(ClusterName %in% ClusterName_MGL_PVM, ClusterName, Subclass))

fct_celltype_labels <- c(ClusterName_MGL_PVM, setdiff(unique(celltype_labels$label), ClusterName_MGL_PVM))
l5.seurat$celltype_labels <- factor(celltype_labels$label, fct_celltype_labels)

## Subset microglia & PVM clusters
l5_PVM_MGL <- subset(l5.seurat, ClusterName %in% c("PVM1", "PVM2", "MGL3", "MGL2", "MGL1"))
l5_PVM_MGL <- scrna_process(l5_PVM_MGL, npc = 50, res = 0.1, normalize = TRUE)
l5_PVM_MGL$ClusterName <- factor(l5_PVM_MGL$ClusterName, ClusterName_MGL_PVM)
l5_PVM_MGL <- SetIdent(l5_PVM_MGL, value = "ClusterName")
save(l5_PVM_MGL, file = "data/mouse_brain/l5_PVM_MGL.RData")

# Extended Data Fig. 2b: Dotplot of immune checkpoint molecules -----------------
scale_breaks <- function(x) {
  breaks_extend <- extended_breaks(n = 4)(x)
  if(round(breaks_extend[2]) != breaks_extend[2]) {
    c(0.1, breaks_extend[-1])
  }
  else{c(1, breaks_extend[-1])}
}

DotPlot_custom <- function(SeuratData, genes, group_var, short_class = FALSE, 
                           pt.size = 8, font.size = 8, 
                           width = 7, height = 7, break_ls = NULL) {
  p_split <- lapply(genes, function(x) {
    p_sub <- DotPlot(SeuratData, features = x, group.by = group_var,
                     cols = c("blue", "red"), scale = FALSE) +
      coord_flip() +
      labs(x = NULL, y = NULL) +
      guides(color = guide_colorbar(order = 1, barheight = 3, 
                                    barwidth = 0.5, title = "Avg Expr"),
             size = guide_legend(order = 2, title = "Pct Expr")) +
      theme_cowplot(font_size = font.size) +
      theme(legend.box = "horizontal",
            axis.text.y = element_text(size = 9))
    
    if(x == tail(genes, 1)) {
      p_sub <- p_sub + 
        theme(axis.text.x = element_text(hjust = 1, angle = 45))}
    else {
      p_sub <- p_sub + 
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank())}  
  })
  
  p <- egg::ggarrange(plots = p_split, ncol = 1)
  
  filename <- paste0("figures/Fig.s2b_s2c_l5_mousebrain/Fig.s2b_DotPlot_by_", group_var, "_split")
  ggsave(paste0(filename, ".png"), p, width = width, height = height, dpi = 400)
  ggsave(paste0(filename, ".pdf"), p, width = width, height = height)
  
  return(p_split$data)
}

Fig2b_dotplot <- DotPlot_custom(
  l5.seurat, 
  genes = c("Havcr2", "Lag3", "Vsir", "Pdcd1", "Cd274", "Ctla4"), 
  group_var = "celltype_labels",
  pt.size = 7, font.size = 10, 
  height = 7.5, width = 7,
  break_ls = list(c(1, 1:3 * 4), c(1, 1:3 * 10),
                  c(1, 1:3 * 20), c(0.1, 1:3 * 0.6),
                  c(0.1, 1:3 * 0.6), c(0.01, 1:3 * 0.05))
)

Fig2b_dotplot %>%
  write.csv("Source_Data/Fig.s2c_dotplot_checkpoint_Tgfb.csv")

# Extended Data Fig. 2c: microglia and PVM clusters -----------------
## (1) Dimplot colored by ClusterNames ---------------------
pal <- c(brewer.pal(9, "Reds")[c(3, 6, 9)], brewer.pal(6, "Blues")[c(4, 6)], 
         brewer.pal(8, "Pastel2"), brewer.pal(9, "Pastel1"))

p <- DimPlot(l5_PVM_MGL, reduction = "umap", group.by = "ClusterName", 
             cols = pal[1:5], label = FALSE) + NoLegend() + labs(title = NULL)
LabelClusters(p, id = "ClusterName",  size = 4, repel = TRUE, box = TRUE,
              color = c("black", rep("white", 3), "black"), label.size = 0.5,
              fontface = "bold", max.overlaps = Inf, seed = 42,
              nudge_x = c(rep(-0.4, 3), rep(0.4, 2)))

ggsave("figures/Fig.s2b_s2c_l5_mousebrain/Fig.s2c_DimPlot_PVM_MGL_ClusterName.png", width = 6, height = 4.5, dpi = 400)
ggsave("figures/Fig.s2b_s2c_l5_mousebrain/Fig.s2c_DimPlot_PVM_MGL_ClusterName.pdf", width = 6, height = 4.5)

## (2) Microglia dissociation signature score -------------

## Compute microglia dissociation signature score (Marsh et al, 2022, PMID: 35260865)
Diss_MG_sig <- read.csv("data/mouse_brain/Marsh_2022_TableS4.csv") %>%
  select(-contains("X")) %>%
  `colnames<-`(str_replace(colnames(.), "\\.\\.", "\\.")) %>%
  pull(Micro.Myeloid.Shared.Act.Score) %>%
  setdiff("")

l5_PVM_MGL <- AddModuleScore(l5_PVM_MGL, list(Diss_MG_sig), name = "Microglia_Dissociation")

## Feature plot
FeaturePlot(l5_PVM_MGL, features = "Microglia_Dissociation1", 
            cols = rev(brewer.pal(11, "Spectral")), label = TRUE, label.size = 4, repel = TRUE) +  
  labs(x = NULL, y = NULL, title = "Microglia Dissociation") + NoLegend()

ggsave("figures/Fig.s2b_s2c_l5_mousebrain//Fig.s2c_featureplot_Microglia_Dissociation_PVM_MGL.png",
       width = 6, height = 4.5, dpi = 400)
ggsave("figures/Fig.s2b_s2c_l5_mousebrain//Fig.s2c_featureplot_Microglia_Dissociation_PVM_MGL.pdf",
       width = 6, height = 4.5, dpi = 400)

## Violin plot
VlnPlot(l5_PVM_MGL, features = "Microglia_Dissociation1", 
        group.by = "ClusterName", pt.size = 0.1, cols = pal) +  
  labs(x = NULL, y = NULL, title = "Microglia Dissociation") + NoLegend()
ggsave("figures/Fig.s2b_s2c_l5_mousebrain/Fig.s2c_vln_Microglia_Dissociation_PVM_MGL.png",
       width = 6, height = 4.5, dpi = 400)
ggsave("figures/Fig.s2b_s2c_l5_mousebrain/Fig.s2c_vln_Microglia_Dissociation_PVM_MGL.pdf",
       width = 6, height = 4.5, dpi = 400)

