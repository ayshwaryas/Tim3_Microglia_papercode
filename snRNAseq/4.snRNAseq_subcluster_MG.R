source("snRNAseq/0.snRNAseq_functions.R")

# Combining and clustering subcluster 2, 3, 6, 18 ---------------------
load("results/doubletfinder/nucseq_harmony_subcluster_doubletFinder.RData") ## Microglia/PVM clusters
nucseq_harmony_subcluster_rm_doublets <- subset(nucseq_harmony_subcluster, DF.classifications_0.25_0.01_3433 == "Singlet")


nucseq_harmony_MG_2_3_6_18 <- subset(nucseq_harmony_subcluster_rm_doublets, seurat_clusters %in% c(2, 3, 6, 18))
nucseq_harmony_MG_2_3_6_18$orig.subcluster <- nucseq_harmony_MG_2_3_6_18$seurat_clusters


nucseq_harmony_MG_2_3_6_18 <- scrna_process(
  nucseq_harmony_MG_2_3_6_18, run_harmony = TRUE, normalize = FALSE, npc = 50, res = 1)


# Subclustering cluster 5 ------------------------------
nucseq_harmony_MG_2_3_6_18_sub5 <- subset(nucseq_harmony_MG_2_3_6_18, seurat_clusters == 5)
nucseq_harmony_MG_2_3_6_18_sub5 <- scrna_process(
  nucseq_harmony_MG_2_3_6_18_sub5, run_harmony = FALSE, 
  npc = 50, res = 1, normalize = FALSE)

nucseq_harmony_MG_2_3_6_18_sub5$seurat_clusters <- paste0(
  "5_", as.character(nucseq_harmony_MG_2_3_6_18_sub5$seurat_clusters))

new_clusters <- nucseq_harmony_MG_2_3_6_18@meta.data %>%
  rownames_to_column("Cells") %>%
  left_join(nucseq_harmony_MG_2_3_6_18_sub5 %>% rownames_to_column("Cells"),
            suffix = c("", ".sub5"), by = "Cells") %>%
  mutate(new_clusters = ifelse(
    is.na(seurat_clusters.sub5), as.character(seurat_clusters),
    seurat_clusters.sub5)) %>%
  mutate(new_clusters = factor(new_clusters)) %>%
  mutate(Genotype_labels = str_replace(as.character(Genotype), "\\.", " ")) %>%
  mutate(Genotype_labels = str_replace(Genotype_labels, "Tim3_cKO", "<i>Havcr2</i><sup>icKO</sup>")) %>%
  mutate(Genotype_labels = factor(Genotype_labels, levels = c("control", "<i>Havcr2</i><sup>icKO</sup>", 
                                                              "5XFAD", "<i>Havcr2</i><sup>icKO</sup> 5XFAD") ))


nucseq_harmony_MG_2_3_6_18$new_clusters <- new_clusters$new_clusters
nucseq_harmony_MG_2_3_6_18$Genotype_labels <- new_clusters$Genotype_labels

save(nucseq_harmony_MG_2_3_6_18, file = "R_objects/2022-09-08.nucseq_harmony_MG_recluster_sub_2_3_6_18.RData")

# Split Havcr2cKO 5XFAD nuclei in cluster 2 into P1 and P2 based on Hallmark TGFB signature score --------------------
Hallmark_TGFB <- read.table("data/signatures/HALLMARK_TGF_BETA_SIGNALING.txt", header = FALSE)$V1
nucseq_harmony_MG_2_3_6_18 <- AddModuleScore(nucseq_harmony_MG_2_3_6_18, list(Hallmark_TGFB), name = "Hallmark_TGFB")

sub2 <- nucseq_harmony_MG_2_3_6_18@meta.data %>%
  rownames_to_column("Cells") %>%
  mutate(bimod_genotype = case_when(
    seurat_clusters == 2 & Genotype == "Tim3_cKO.5XFAD" & Hallmark_TGFB1 > 0 ~ "Tim3_cKO.5XFAD_Top",
    seurat_clusters == 2 & Genotype == "Tim3_cKO.5XFAD" & Hallmark_TGFB1 < 0 ~ "Tim3_cKO.5XFAD_Btm",
    TRUE ~ as.character(Genotype))) %>%
  mutate(bimod_genotype = factor(bimod_genotype, c("control", "Tim3_cKO", "5XFAD", "Tim3_cKO.5XFAD",
                                                   "Tim3_cKO.5XFAD_Top", "Tim3_cKO.5XFAD_Btm"))) %>%
  select(Cells, Hallmark_TGFB1, bimod_genotype) %>%
  column_to_rownames("Cells")
nucseq_harmony_MG_2_3_6_18 <- AddMetaData(nucseq_harmony_MG_2_3_6_18, sub2)


# Getting marker genes of microglia clusters using FindAllMarkers --------------------
nucseq_harmony_MG_2_3_6_18 <- SetIdent(nucseq_harmony_MG_2_3_6_18, value = "new_clusters")
nucseq_harmony_MG_2_3_6_18_markers <- FindAllMarkers(
  nucseq_harmony_MG_2_3_6_18, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

save(nucseq_harmony_MG_2_3_6_18_markers, file = "results/2022-09-09.nucseq_harmony_MG_2_3_6_18_markers.RData")
