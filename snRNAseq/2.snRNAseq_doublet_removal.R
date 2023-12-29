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

load("R_objects/nucseq_harmony_dim20_res1.RData")
load("R_objects/nucseq_harmony_subcluster.RData")

# Doublet removal using DoubletFinder ---------------
## (1) pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_nucseq <- paramSweep_v3(nucseq_harmony, PCs = 1:20, sct = FALSE)
sweep.stats_nucseq <- summarizeSweep(sweep.res.list_nucseq, GT = FALSE)
bcmvn_nucseq <- find.pK(sweep.stats_nucseq)

ggplot(bcmvn_nucseq, aes(x = pK, y = BCmetric)) +
  geom_point() +
  geom_line(group = 1) +
  theme_cowplot(font_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave('results/doubletfinder/choose_pK.png', width = 8, height = 4)

## (2) Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(nucseq_harmony@meta.data$seurat_clusters)
nExp_poi <- round(0.054 * nrow(nucseq_harmony@meta.data))  ## Assuming 5.4% doublet formation rate
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


## (3) Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
nucseq_harmony <- doubletFinder_v3(nucseq_harmony, PCs = 1:20, pN = 0.25, pK = 0.01, 
                                   nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nucseq_harmony <- doubletFinder_v3(nucseq_harmony, PCs = 1:20, pN = 0.25, pK = 0.01, 
                                   nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.01_3594", sct = FALSE)

### UMAP colored by doublets (red) and singlets (grey)
p <- DimPlot(nucseq_harmony, group.by = "DF.classifications_0.25_0.01_3594", 
             cols = c("red", "grey80"), pt.size = 0.05)
p$data$seurat_clusters <- nucseq_harmony$seurat_clusters
p <- LabelClusters(p, id = "seurat_clusters", size = 3.5, repel = TRUE)
ggsave('results/doubletfinder/dimplot_DF.classifications_0.25_0.01_3433_1st_run.png', 
       p, width = 10, height = 7, dpi = 400)

### UMAP colored by probability of doublet
p <- DimPlot(nucseq_harmony, group.by = "DF.classifications_0.25_0.01_3433",
             cols = c("red", "grey80"), pt.size = 0.05)
p$data$seurat_clusters <- nucseq_harmony$seurat_clusters
p <- LabelClusters(p, id = "seurat_clusters", size = 3.5, repel = TRUE)
ggsave('results/doubletfinder/dimplot_DF.classifications_0.25_0.01_3433_2nd_run.png', 
       p, width = 10, height = 7, dpi = 400)

save(nucseq_harmony, file = "results/doubletfinder/nucseq_harmony_PC20_res1_doubletFinder.RData")

# Match DoubletFinder results to microglia/PVM subclusters ------------------------
MG_subcluster_doubletfinder <- nucseq_harmony_subcluster@meta.data %>%
  rownames_to_column("Cells") %>%
  left_join(nucseq_harmony@meta.data %>% rownames_to_column("Cells") %>%
              select(Cells, pANN_0.25_0.01_3594, DF.classifications_0.25_0.01_3594, DF.classifications_0.25_0.01_3433),
            by = "Cells") %>%
  column_to_rownames("Cells")
nucseq_harmony_subcluster <- AddMetaData(nucseq_harmony_subcluster, MG_subcluster_doubletfinder)

save(nucseq_harmony_subcluster, "results/doubletfinder/nucseq_harmony_subcluster_doubletFinder.RData")
