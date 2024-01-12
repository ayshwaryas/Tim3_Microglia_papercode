source("scRNAseq/0.scRNAseq_functions.R")

load("R_objects/scRNAseq_processed_no_intron_cortex_only.RData")

# 1. Remove cluster 5 (low-quality) and 9 (PVM) -----------------------------------
scRNAseq_cortex_MG <- subset(scRNAseq_cortex, seurat_clusters %in% c(0:4, 6:8, 10))
scRNAseq_cortex_MG <- RunUMAP(scRNAseq_cortex_MG, dims = 1:50)
scRNAseq_cortex_MG <- SetIdent(scRNAseq_cortex_MG, value = "new_clusters")

# 2. Identify MGnD, Homeostasis and Cycling clusters -------------------------------

## (1) Identify Homeostasis/MGnD clusters using Clec7a+ signatures (Top 100 DEGs of Clec7a+ vs Clec7a-) --------------
Clec7a_sig <- read.csv("results/bulkRNAseq_results_Clec7a_Krasemann_2017.csv") %>%
  arrange(padj) %>% group_by(direction) %>% 
  dplyr::slice(1:100) %>%
  split(f = .$direction) %>%
  lapply(pull, tracking_id)

scRNAseq_cortex_MG <- AddModuleScore(scRNAseq_cortex_MG, list(Clec7a_sig$up), name = "MGnD")
scRNAseq_cortex_MG <- AddModuleScore(scRNAseq_cortex_MG, list(Clec7a_sig$down), name = "Homeostasis")

## (2) Identify Cycling clusters ----------------------------------------------------
## Convert human gene symbols to mouse gene symbols
mouse_human_genes_raw <- read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
convert_human_to_mouse <- function(gene_list){
  mouse_human_genes_raw %>%
    filter(DB.Class.Key %in% subset(mouse_human_genes_raw, Symbol %in% gene_list )$DB.Class.Key) %>% 
    filter(Common.Organism.Name == "mouse, laboratory") %>%
    mutate(Symbol = ifelse(is.na(Symbol), str_to_title(Symbol), Symbol)) %>%
    pull(Symbol) %>% unique
}

s.genes_mouse <- convert_human_to_mouse(cc.genes.updated.2019$s.genes)
g2m.genes_mouse <- convert_human_to_mouse(cc.genes.updated.2019$g2m.genes)

scRNAseq_cortex_MG <- CellCycleScoring(scRNAseq_cortex_MG, s.features = s.genes_mouse, g2m.features = g2m.genes_mouse)


## (3) Identify IFN clusters ----------------------------------------------------
## IFN signature: Top 20 genes of IFN-rich microglia ordered by fold change (Ellwanger et al. 2021, PMID: 33446504)
## Ellwanger et al. 2021 (PMID: 33446504) Dataset S1B 
pnas_sd01 <- read.csv("data/pnas.2017742118.sd01.csv")

## Select top 20 marker genes for each microglia population
pnas_sig <- sapply(str_remove(colnames(pnas_sd01)[4 + 6 * 0:10], "_Specificity"), function(MG_type) {
  pnas_sd01 %>%
    select(Symbol, contains(MG_type)) %>%
    arrange(desc(!!sym(paste0(MG_type, "_Fold_change")))) %>%
    dplyr::slice(1:20) %>%
    pull(Symbol)
}, simplify = FALSE, USE.NAMES = TRUE) 

scRNAseq_cortex_MG <- AddModuleScore(
  scRNAseq_cortex_MG, list(pnas_sig$IFNR), name = "IFNR")


## Visualization of signature score
FeaturePlot(scRNAseq_cortex_MG, c("MGnD1", "Homeostasis1", "S.Score", "G2M.Score", "IFNR1"),
            label = TRUE, label.size = 4, repel = TRUE, cols = rev(RColorBrewer::brewer.pal(11, "Spectral"))) 
ggsave("figures/cortex_combine_clusters/featureplot_clec7a_cycling.png", width = 12, height = 12)

VlnPlot(scRNAseq_cortex_MG, c("MGnD1", "Homeostasis1", "S.Score", "G2M.Score", "IFNR1"),
        cols = c25, ncol = 1) & theme(axis.title = element_blank())
ggsave("figures/cortex_combine_clusters/vln_clec7a_cycling.png", width = 7, height = 11)

VlnPlot(scRNAseq_cortex_MG, c("MGnD1", "Homeostasis1", "S.Score", "G2M.Score", "IFNR1"),
        cols =  c("#785EF0", "#E69F00"), ncol = 1, split.by = "Genotype") &
  theme(axis.title = element_blank())
ggsave("figures/cortex_combine_clusters/vln_clec7a_cycling_split.png", width = 9, height = 11)

# 3. Reclustering -----------------

## (1) Combine and recluster Homeostatic & IFN clusters ------------------------
scRNAseq_cortex_Homeo_recluster <- subset(scRNAseq_cortex, new_clusters %in% c(0, 1, 4, "6_0","6_2", "7_0", "7_1"))
scRNAseq_cortex_Homeo_recluster <- scrna_process(scRNAseq_cortex_Homeo_recluster, npc = 50, res = 0.6, normalize = FALSE)


## (2) Combine and recluster DAM clusters --------------------------------------
### Initial clustering
scRNAseq_cortex_DAM_recluster <- subset(scRNAseq_cortex, new_clusters %in% c(2, 3, "7_0", "7_1", "8_0", 10))
scRNAseq_cortex_DAM_recluster <- scrna_process(scRNAseq_cortex_DAM_recluster, npc = 50, res = 0.6, normalize = FALSE)

### Combine subcluster 5 with cluster 7_0, 7_1 and recluster
scRNAseq_cortex_DAM_recluster_sub5 <- subset(
  scRNAseq_cortex_DAM_recluster, new_clusters %in% c("7_0", "7_1") | RNA_snn_res.0.6 == 5)
scRNAseq_cortex_DAM_recluster_sub5 <- scrna_process(scRNAseq_cortex_DAM_recluster_sub5, npc = 50, res = 0.4, normalize = FALSE)


# 4. Assign annotations --------------------------------------------------------
updated_clusters <- scRNAseq_cortex_MG@meta.data %>%
  rownames_to_column("Cells") %>%
  left_join(rownames_to_column(scRNAseq_cortex_Homeo_recluster@meta.data, "Cells") %>%
              select(Cells, cluster.HMG = RNA_snn_res.0.6), by = "Cells") %>%
  left_join(rownames_to_column(scRNAseq_cortex_DAM_recluster@meta.data, "Cells") %>%
              select(Cells, cluster.DAM = RNA_snn_res.0.6), by = "Cells") %>%
  left_join(rownames_to_column(scRNAseq_cortex_DAM_recluster_sub5@meta.data, "Cells") %>%
              select(Cells, cluster.old7.new5 = RNA_snn_res.0.4), by = "Cells") %>% 
  mutate(updated_clusters = case_when(
    ## (1) Homeostatic & IFN ---------------
    cluster.old7.new5 == 0 ~ "HMG_2",
    cluster.HMG == 3 ~ "HMG_1",
    cluster.HMG == 8 ~ "HMG_3",
    cluster.HMG == 4 ~ "IFN_HMG",
    cluster.HMG == 5 ~ "IFN_DAM",
    cluster.HMG == 7 | new_clusters %in% c("6_0", "6_2") ~ "CC_HMG_0",
    cluster.HMG %in% c(0, 1, 2) | new_clusters %in% c(0, 1, 4)  ~ "HMG_0",
    ## (2) DAM -----------------------------
    cluster.old7.new5 == 1 ~ "DAM_6",
    cluster.DAM  == 4 | new_clusters == "8_0" ~ "DAM_4",
    new_clusters == 3 ~ "DAM_0",
    cluster.DAM  == 6 | new_clusters == 10 ~ "DAM_5",
    cluster.DAM  == 1 ~ "DAM_1",
    cluster.DAM  == 7 ~ "DAM_3",
    cluster.DAM %in% c(2, 3) | new_clusters == 2 ~ "DAM_2",
    ## (3) Cycling ------------------------
    new_clusters == "6_1" ~ "CC_DAM_0",
    new_clusters == "6_3" ~ "CC_HMG_1",
    new_clusters == "6_4" ~ "CC_DAM_1",
    new_clusters == "6_5" ~ "CC_HMG_2",
    new_clusters == "6_6" ~ "CC_DAM_2",
    new_clusters == "8_1" ~ "CC_DAM_3",
    new_clusters == "8_2" ~ "CC_DAM_4")) %>%
  mutate(updated_clusters = factor(
    updated_clusters, c(paste0("HMG_", 0:3), paste0("DAM_", 0:6),
                        paste0("CC_HMG_", 0:2), paste0("CC_DAM_", 0:4), 
                        "IFN_HMG", "IFN_DAM"))) %>%
  select(Cells, updated_clusters) %>%
  column_to_rownames("Cells")

scRNAseq_cortex_MG <- AddMetaData(scRNAseq_cortex_MG, updated_clusters)
scRNAseq_cortex_MG <- SetIdent(scRNAseq_cortex_MG, value = "updated_clusters")

save(scRNAseq_cortex_MG, file = "R_objects/2023-10-03.scRNAseq_cortex_MG_updated_clusters.RData")
