source("snRNAseq/0.snRNAseq_functions.R")

# Add Azimuth annotations to the Seurat object ------------------------------

## Load Seurat object
load("results/doubletfinder/nucseq_harmony_PC20_res1_doubletFinder.RData")  ## All cells
load("results/doubletfinder/nucseq_harmony_subcluster_doubletFinder.RData") ## Microglia/PVM clusters

## Load Azimuth output 
azimuth_pred <- read.delim(paste0("results/azimuth/azimuth_pred_harmony.tsv"), row.names = 1) 

all(colnames(nucseq_harmony) == rownames(azimuth_pred))
nucseq_harmony <- AddMetaData(object = nucseq_harmony, metadata = azimuth_pred)

# Subset the Seurat object to remove doublets --------------------
## All cells
nucseq_harmony_rm_doublets <- subset(nucseq_harmony, DF.classifications_0.25_0.01_3433 == "Singlet")

## Microglia/PVM clusters
nucseq_harmony_subcluster_rm_doublets <- subset(nucseq_harmony_subcluster, DF.classifications_0.25_0.01_3433 == "Singlet")


# Finalize annotations ---------------
## Split L2/3 IT into 3 groups 
nucseq_harmony_rm_doublets_L23_IT_raw <- subset(nucseq_harmony_rm_doublets, predicted.subclass == "L2/3 IT")
nucseq_harmony_rm_doublets_L23_IT <- RunUMAP(nucseq_harmony_rm_doublets_L23_IT_raw,
                                             reduction = "harmony", dims = 1:20)
DimPlot(nucseq_harmony_rm_doublets_L23_IT, cols = c25)
ggsave("figures/nucseq_harmony_rm_doublets_L23_IT.png")

nucseq_harmony_rm_doublets_L23_IT_group <- FetchData(
  nucseq_harmony_rm_doublets_L23_IT, c("UMAP_1", "UMAP_2", "seurat_clusters")) %>%
  cbind(Cells = Cells(nucseq_harmony_rm_doublets_L23_IT)) %>%
  mutate(group = case_when(UMAP_1 < -9 ~ "L2/3 IT_2", UMAP_1 > 0 ~ "L2/3 IT_1", TRUE ~ "L2/3 IT_3"))

nucseq_harmony_rm_doublets_L23_IT$group <- nucseq_harmony_rm_doublets_L23_IT_group$group

DimPlot(nucseq_harmony_rm_doublets_L23_IT, group.by = "group", cols = c25)
ggsave("figures/nucseq_harmony_rm_doublets_L23_IT_group.png")


annot_levs <- c("Microglia", "Microglia-2", "BAM-PVM", 
                "Astro", "Endo", 
                "L2/3 IT_1", "L2/3 IT_2", "L2/3 IT_3", "L5 ET", "L5 IT", "L5/6 NP", "L6b", "L6 CT", "L6 IT", 
                "Oligo", "OPC", "Peri", "VLMC", "Cycling", "Meis2", 
                "Lamp5", "Pvalb", "Sncg", "Sst", "Sst Chodl", "Vip")
annot_new <- nucseq_harmony_rm_doublets@meta.data %>%
  rownames_to_column("Cells") %>%
  mutate(annot = case_when(
    seurat_clusters == 27 ~ "Cycling",
    Cells %in% Cells(subset(nucseq_harmony_subcluster_rm_doublets, seurat_clusters %in% c(2, 3, 6, 18))) ~ "Microglia",
    Cells %in% Cells(subset(nucseq_harmony_subcluster_rm_doublets, seurat_clusters == 12)) ~ "BAM-PVM",
    Cells %in% Cells(subset(nucseq_harmony_subcluster_rm_doublets, seurat_clusters == 15)) ~ "BAM",
    TRUE ~ predicted.subclass)) %>%
  left_join(nucseq_harmony_rm_doublets_L23_IT_group[, c("Cells", "group")], by = "Cells") %>%
  mutate(annot_new = ifelse(!is.na(group), group, annot)) %>%
  mutate(annot_new = ifelse(annot_new == "L6 IT Car3", "L6 IT", annot_new)) %>%
  mutate(annot_new = ifelse(annot_new == "BAM", "BAM-PVM", annot_new)) %>%
  mutate(annot_new = ifelse(annot_new == "Micro-PVM", "Microglia-2", annot_new)) %>%
  mutate(annot_new = factor(annot_new, annot_levs)) 

nucseq_harmony_rm_doublets$annot_new <- annot_new$annot_new

save(nucseq_harmony_rm_doublets, file = "R_objects/nucseq_harmony_rm_doublets_annotated.RData")