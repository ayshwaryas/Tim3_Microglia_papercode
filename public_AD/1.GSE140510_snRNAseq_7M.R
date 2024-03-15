source("snRNAseq/0.snRNAseq_functions.R")

# 1. Creating and merging Seurat objects ---------------------------------
# system("tar -xvf data/GSE140510_RAW.tar -C data/GSE140510")
sample_names <- list.files("data/GSE140510") %>%
  str_remove("(_(barcodes|features).tsv.gz)|_matrix.mtx.gz") %>%
  unique

for(i in sample_names) {
  system(paste0("mkdir data/GSE140510/", i))
  system(paste0("mv data/GSE140510/", i, "_barcodes.tsv.gz data/GSE140510/", i, "/barcodes.tsv.gz"))
  system(paste0("mv data/GSE140510/", i, "_features.tsv.gz data/GSE140510/", i, "/features.tsv.gz"))
  system(paste0("mv data/GSE140510/", i, "_matrix.mtx.gz data/GSE140510/", i, "/matrix.mtx.gz"))
}

GSE140510_obj <- sapply(sample_names, function(i) {
  expr <- Read10X(paste0("data/GSE140510/", i))
  CreateSeuratObject(counts = expr, project = i)
}, simplify = FALSE, USE.NAMES = TRUE)


GSE140510_obj <- merge(
  x = GSE140510_obj$GSM4173504_WT_1, 
  y = c(GSE140510_obj$GSM4173505_WT_2, GSE140510_obj$GSM4173506_WT_3,
        GSE140510_obj$GSM4173510_WT_5XFAD_1, GSE140510_obj$GSM4173511_WT_5XFAD_2,
        GSE140510_obj$GSM4173512_WT_5XFAD_3),
  add.cell.ids = c("GSM4173504_WT_1", "GSM4173505_WT_2", "GSM4173506_WT_3",
                   "GSM4173510_WT_5XFAD_1", "GSM4173511_WT_5XFAD_2", "GSM4173512_WT_5XFAD_3"),
  project = "GSE140510"
)

# 2. QC ----------------------------------
GSE140510_obj[["percent.mt"]] <- PercentageFeatureSet(GSE140510_obj, pattern = "^mt-")
GSE140510_obj <- subset(GSE140510_obj, subset = nCount_RNA>= 1000 & nFeature_RNA>= 200 & percent.mt < 80)

# 3. Processing  ----------------------------------
GSE140510_obj <- scrna_process(GSE140510_obj, npc = 50, res = 1, normalize = TRUE)
GSE140510_obj$orig.ident <- str_remove(GSE140510_obj$orig.ident, "GSM[0-9]+_")
GSE140510_meta <- data.frame(
  orig.ident = GSE140510_obj$orig.ident) %>%
  mutate(Genotype = str_extract(orig.ident, "5XFAD")) %>%
  mutate(Genotype = ifelse(is.na(Genotype), 'WT', '5XFAD')) %>%
  mutate(Genotype = factor(Genotype, c("WT", "5XFAD"))) 

GSE140510_obj <- AddMetaData(GSE140510_obj, GSE140510_meta)

save(GSE140510_obj, file = "R_objects/GSE140510_Admouse_7M.RData")


# 4. Identify microglia clusters ----------------------------
## ==> Cluster 14, 19, 36, 40, 44, 45 are potentially microglia or PVM
dotplot_markers(GSE140510_obj, width = 12, height = 10, suffix = "GSE140510")

# 5. Subsetting and subclustering microglia clusters (14, 19, 36, 40, 44, 45) ------------------------------------
GSE140510_MG_obj <- subset(GSE140510_obj, seurat_clusters %in% c(14, 19, 36, 40, 44, 45))
GSE140510_MG_obj$orig.clusters <- GSE140510_MG_obj$seurat_clusters
GSE140510_MG_obj$orig.clusters <- factor(GSE140510_MG_obj$orig.clusters, c(14, 19, 36, 40, 44, 45))
GSE140510_MG_obj <- scrna_process(GSE140510_MG_obj, npc = 50, res = 1, normalize = FALSE)

save(GSE140510_MG_obj, file = "R_objects/GSE140510_ADmouse_7M_MG.RData")

# 6. Identify DAM/MGnD cluster -----------------------------
## ==> Subcluster 1 is DAM/MGnD
## MGnD and Homeostasis signature (Top 100 DEGs up/down in Clec7a+ vs Clec7a-)

Clec7a_sig <- read.csv("results/bulkRNAseq_results_Clec7a_Krasemann_2017.csv") %>%
  arrange(padj) %>% group_by(direction) %>% 
  dplyr::slice(1:100) %>%
  split(f = .$direction) %>%
  lapply(pull, tracking_id)

## Compute MGnD and Homeostasis signature score
GSE140510_MG_obj <- AddModuleScore(GSE140510_MG_obj, list(Clec7a_sig$up), name = "MGnD")
GSE140510_MG_obj <- AddModuleScore(GSE140510_MG_obj, list(Clec7a_sig$down), name = "Homeostasis")

for(x in c("MGnD", "Homeostasis")) {
  FeaturePlot(GSE140510_MG_obj, features = paste0(x, 1), cols = rev(brewer.pal(11, "Spectral")),
              label = TRUE, label.size = 3, repel = TRUE) + labs(title = x)
  ggsave(paste0("figures/GSE140510/sig_score_GSE140510_MG_", x, ".png"), width = 5, height = 3.5)
  
  VlnPlot(GSE140510_MG_obj, features = paste0(x, 1), cols = c25) + NoLegend()
  ggsave(paste0("figures/GSE140510/vln_GSE140510_MG_", x, ".png"), width = 4, height = 2.5)
  
  VlnPlot(GSE140510_MG_obj, features = paste0(x, 1), cols = c25, split.by = "Genotype") 
  ggsave(paste0("figures/GSE140510/vln_GSE140510_MG_", x, "_by_genotype.png"), width = 5.5, height = 2.5)
}
