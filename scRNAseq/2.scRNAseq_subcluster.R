source("scRNAseq/0.scRNAseq_functions.R")

load("R_objects/scRNAseq_processed_no_intron_cortex_only.RData")


# Subclustering cluster 6 ------------------------
scRNAseq_cortex_sub6 <- subset(scRNAseq_cortex, seurat_clusters == 6)
scRNAseq_cortex_sub6 <-  scrna_process(scRNAseq_cortex_sub6, npc = 50, res = 1, normalize = FALSE)
save(scRNAseq_cortex_sub6, file = "R_objects/scRNAseq_processed_no_intron_cortex_only_sub6.RData")

# Subclustering cluster 7 ------------------------
scRNAseq_cortex_sub7 <- subset(scRNAseq_cortex, seurat_clusters == 7)
scRNAseq_cortex_sub7 <- scrna_process(scRNAseq_cortex_sub7, npc = 50, res = 0.4, normalize = FALSE)
save(scRNAseq_cortex_sub7, file = "R_objects/scRNAseq_processed_no_intron_cortex_only_sub7.RData")


# Subclustering cluster 8 ---------------------------------------
scRNAseq_cortex_sub8 <- subset(scRNAseq_cortex, seurat_clusters == 8)
scRNAseq_cortex_sub8 <- scrna_process(scRNAseq_cortex_sub8, npc = 50, res = 0.4, normalize = FALSE)
save(scRNAseq_cortex_sub8, file = "R_objects/scRNAseq_processed_no_intron_cortex_only_sub8.RData")


# Update cluster labels including the subcluster labels ------------------------------
new_clusters <- scRNAseq_cortex@meta.data %>%
  rownames_to_column("Cells") %>%
  select(Cells, seurat_clusters) %>%
  left_join(scRNAseq_cortex_sub6@meta.data %>% rownames_to_column("Cells") %>%
              select(Cells, seurat_clusters),
            suffix = c("", ".sub6")) %>%
  left_join(scRNAseq_cortex_sub7@meta.data %>% rownames_to_column("Cells") %>% 
              select(Cells, seurat_clusters),
            suffix = c("", ".sub7")) %>%
  left_join(scRNAseq_cortex_sub8@meta.data %>% rownames_to_column("Cells") %>% 
              select(Cells, seurat_clusters),
            suffix = c("", ".sub8")) %>%
  mutate(new_clusters = case_when(
    seurat_clusters == 6 ~ seurat_clusters.sub6,
    seurat_clusters == 7 ~ seurat_clusters.sub7,
    seurat_clusters == 8 ~ seurat_clusters.sub8,
    TRUE ~ seurat_clusters
  ))  %>%
  select(Cells, new_clusters) %>%
  column_to_rownames("Cells")


scRNAseq_cortex <- AddMetaData(scRNAseq_cortex, new_clusters)
save(scRNAseq_cortex, file = "R_objects/scRNAseq_processed_no_intron_cortex_only.RData")






