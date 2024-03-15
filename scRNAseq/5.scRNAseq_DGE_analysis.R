source("snRNAseq/0.snRNAseq_functions.R")

# Function for performing differential gene expression analysis
scrna_DGE <- function(SeuratData, logfc.threshold = 0, min.pct = 0.1, 
                      test.use = "wilcox", rm_ribosome = FALSE) {
  if(rm_ribosome) {
    features <- rownames(SeuratData)[!str_detect(rownames(SeuratData), "^(Rp[ls]|Mrp[ls])")]
  }
  else { features <- NULL }
  clusters <- levels(SeuratData$updated_clusters)
  results <- mclapply(clusters, function(x) {
    res <- FindMarkers(subset(SeuratData, updated_clusters == x), features = features,
                       ident.1 = "Tim3icKO_5XFAD", ident.2 = "5XFAD", group.by = "Genotype", 
                       logfc.threshold = logfc.threshold, min.pct = min.pct, test.use = test.use)
    if(nrow(res) == 0) { return(NULL) }
    else { res %>% rownames_to_column("gene") %>%
        mutate(direction = ifelse(avg_log2FC > 0, "up", "down")) %>%
        mutate(direction = factor(direction, c("up", "down"))) %>%
        mutate(FDR = p.adjust(.$p_val, method = "BH"))}},
    mc.cores = length(clusters)) 
  names(results) <- clusters
  return(results)
}

## All genes, Havcr2icKO;5xFAD vs 5xFAD (Wilcoxon)
scRNAseq_cortex_MG_DGE_wilcox_min.pct1 <- scrna_DGE(scRNAseq_cortex_MG)

## Remove ribosonal genes,  Havcr2icKO;5xFAD vs 5xFAD (Wilcoxon)
scRNAseq_cortex_MG_DGE_wilcox_min.pct1_rm_ribo <- scrna_DGE(scRNAseq_cortex_MG, rm_ribosome = TRUE)


save(scRNAseq_cortex_MG_DGE_wilcox_min.pct1, scRNAseq_cortex_MG_DGE_wilcox_min.pct1_rm_ribo,
     file = "results/20231005_scRNAseq_cortex_MG_updated_clusters_DGE_wilcox_min.pct1_named.RData")
