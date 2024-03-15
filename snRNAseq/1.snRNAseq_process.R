source("snRNAseq/0.snRNAseq_functions.R")

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
nucseq_harmony <- scrna_process(nucseq, run_harmony = TRUE, npc = 20, res = 1)

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
dotplot_markers(nucseq_harmony)

# 6. Combine and subcluster clusters that are potentially microglia and PVM (6, 11, 32, 35, 41) ----------------------------
## (1) Subclustering ------------------------
nucseq_harmony_subcluster <- subset(nucseq_harmony, seurat_clusters %in% c(6, 11, 32, 35, 41))
nucseq_harmony_subcluster$orig.cluster <- nucseq_harmony$seurat_clusters

nucseq_harmony_subcluster <- scrna_process(
  nucseq_harmony_subcluster, normalize = FALSE,
  run_harmony = TRUE, npc = 20, res = 1)

save(nucseq_harmony_subcluster, file = "R_objects/nucseq_markers_subcluster_harmony.RData")

## (2) Scoring microglia, PVM, BAM, monocyte signatures ------------------------
### Subluster 2, 3, 6, 18 => Microglia
### Sluster 12 => BAM-PVM
### Sluster 15 => BAM


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
