source("snRNAseq/0.snRNAseq_functions.R")

load("R_objects/scRNAseq_processed_no_intron_cortex_only.RData")

# pK Identification (no ground-truth) ------------------------------------------
sweep.res.list_scRNAseq <- paramSweep_v3(scRNAseq_cortex, PCs = 1:50, sct = FALSE)
sweep.stats_scRNAseq <- summarizeSweep(sweep.res.list_scRNAseq, GT = FALSE)
bcmvn_scRNAseq <- find.pK(sweep.stats_scRNAseq)


# Homotypic Doublet Proportion Estimate ----------------------------------------
homotypic.prop <- modelHomotypic(scRNAseq_cortex@meta.data$new_clusters)
nExp_poi <- round(0.076 * nrow(scRNAseq_cortex@meta.data))  ## Assuming 7.6% doublet formation rate
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder with varying classification stringencies -------------------
scRNAseq_cortex <- doubletFinder_v3(scRNAseq_cortex, PCs = 1:50, pN = 0.25, pK = 0.005, 
                                    nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
scRNAseq_cortex <- doubletFinder_v3(scRNAseq_cortex, PCs = 1:50, pN = 0.25, pK = 0.005, 
                                    nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_3289", sct = FALSE)

save(scRNAseq_cortex, file = "results/doubletfinder/scRNAseq_processed_no_intron_cortex_only_doubletFinder.RData")


# Visualization ----------------------------------------------------------------
DimPlot(subset(scRNAseq_cortex, DF.classifications_0.25_0.005_2713 == "Singlet"), group.by = "new_clusters",
        cols = c25, pt.size = 0.05, label = TRUE, repel = TRUE, label.size = 4)
ggsave('figures/doubletfinder/dimplot_rm_doublet_0.25_0.005_2713.png', 
       width = 7.5, height = 5, dpi = 400)

p <- FeaturePlot(scRNAseq_cortex, feature = "pANN_0.25_0.005_3289",
                 cols = rev(brewer.pal(11, "Spectral")), pt.size = 0.05)
p$data$new_clusters <- scRNAseq_cortex$new_clusters
p <- LabelClusters(p, id = "new_clusters", size = 4, repel = TRUE)
ggsave('figures/doubletfinder/feature_pANN_0.25_0.005_3289.png', 
       width = 7.5, height = 5, dpi = 400)