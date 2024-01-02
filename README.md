## Immune checkpoint molecule Tim-3 regulates microglial function and the development of Alzheimer’s disease pathology

### Processing and annotating snRNA-seq data

* `snRNAseq/0.snRNAseq_functions.R`: functions, color palettes and gene signatures used in the analysis
* `snRNAseq/1.snRNAseq_process.R`: QC, processing the snRNA-seq data; combining and subclustering microglia and PVM clusters
* `snRNAseq/2.snRNAseq_doublet_removal.R`: doublet removal using DoubletFinder
* `snRNAseq/3.snRNAseq_annotation.R`: annotating the snRNA-seq clusters
* `snRNAseq/4.snRNAseq_subcluster_MG.R`: further subclustering microglia cells; using FindAllMarkers to get gene expression markers for each cluster; splitting *Havcr2*<sup>icKO</sup> 5XFAD nuclei in cluster 2 (DAM/MGnD) into subpopulations P1 and P2 based on [MSigDB Hallmark TGFβ signaling](https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/HALLMARK_TGF_BETA_SIGNALING.html) signature score

### Fig. 1 and Extended Fig. 3b

* **Fig. 1a, Extended Fig. 3b** (`Fig.1a_s3b_microarray_Thion_2018.R`): Heatmap visualization of the developmental alterations of gene expressions in microglia in a published microarray dataset ([Thion et al. 2018, PMID: 29275859](https://pubmed.ncbi.nlm.nih.gov/29275859/))
    * **Fig. 1a**: Immune checkpoints and TGFβ pathway-related molecules 
    * **Extended Data Fig. 3b**: MGnD and homeostasis associated genes
    
* **Fig. 1b** (`Fig.1b_GSE127449_P9_P28_scRNAseq.R`): Dotplot visualization of gene expression at postnatal stages (P9 and P28) in a published scRNA-seq dataset ([He et al. 2018, PMID: 34982959](https://pubmed.ncbi.nlm.nih.gov/34982959/); [GSE127449](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127449))

* **Fig. 1h** (`Fig.1h_1month_heatmap_DEGs.R`): Heatmap visualization of DEGs (FDR < 0.05) of 1-month-old *Havcr2*<sup>cKO</sup> mice (n = 4 (3 males, 1 female)) compared to *Havcr2*<sup>flox/flox</sup> mice (n = 5 (4 males, 1 female)). Rows represent genes and columns are biological replicates. Gene expression is row-normalized across replicates. 

* **Fig. 1i** (`Fig.1i_1month_volcanoplot.R`): Volcano plot of differential gene expression analysis performed by DESeq2 comparing 1-month-old *Havcr2*<sup>cKO</sup> to *Havcr2*<sup>flox/flox</sup> mice

* **Fig. 1j-l** (`Fig.1j_1l_1month_signature_score.R`): Boxplot visualization of the the scores of MGnD (j), homeostasis (k), and KEGG phagosome (l) signatures in 1-month-old *Havcr2*<sup>cKO</sup> compared to *Havcr2*<sup>flox/flox</sup> mice. The Y-axis represents log2 transformed average TPM, and significance was calculated using the Student's two-tailed t-test.

### Fig. 2

* **Fig. 2b** (`Fig.2b_1month_Tgfbr2_scatter_plot.R`): Scatter plot of genes based on expressional difference represented by log2-transformed fold changes in *Havcr2*<sup>cKO</sup> (X-axis) and *Tgfbr2*<sup>cKO</sup> (Y-axis) compared to control microglia

### Fig. 3

* **Fig. 3g** (`Fig.3g_3month_heatmap_overlapping_DEGs.R`): Heatmap visualization of top DEGs shared by at least two of the following three comparisons (3-month-old mice)
    (1) control phagocytosing versus control non-phagocytosing microglia
    (2) *Havcr2*<sup>cKO</sup> phagocytosing versus control phagocytosing microglia
    (3) *Havcr2*<sup>cKO</sup> non-phagocytosing versus control non-phagocytosing microglia
* **Fig. 3k** (`Fig.3k_3month_nonphago_DEGs_pathway.R`):  Pathway analysis of DEGs in non-phagocytosing microglia from 3-month-old *Havcr2*<sup>cKO</sup> mice compared to control mice. Disease pathways (pathways under section 6.1-6.10 from https://www.kegg.jp/kegg/pathway.html) and ribosomal genes were excluded from the analysis.

### Fig. 4

* **Fig. 4a** (`Fig.4a_overlapping_DEGs_pairwise.R`): Number of overlapped genes and permutation test p-values between each pair of DEGs up- and down-regulated in *Havcr2*<sup>cKO</sup>, phagocytosing, *Tgfbr2*<sup>cKO</sup>, and *Clec7a*<sup>+</sup> microglia compared to control microglia
* **Fig. 4b** (`Fig.4b_overlapping_DEGs_heatmap.R`): Heatmap visualization of correlation between vectors of the expression levels of DEGs of *Havcr2*<sup>cKO</sup> microglia compared to control microglia, as well as the DEGs shared by all three other comparisons (phagocytosing, *Tgfbr2*<sup>cKO</sup>, *Clec7a*<sup>+</sup> microglia compared to their corresponding control groups).

### Fig. 6: snRNA-seq data

* **Fig. 6a** (`Fig.6a_snRNAseq_microglia_UMAP.R`): UMAP visualization of microglia clusters split by genotype
* **Fig. 6b** (`Fig.6b_snRNAseq_microglia_top5markers.R`): Dotplot visualization of top 5 genes distinguishing each cell cluster from Fig. 6a
* **Fig. 6c** (`Fig.6c_snRNAseq_microglia_HallmarkTGFB_score.R`): Probability density curves of the signature score for the Hallmark TGFβ pathway in the 5xFAD (purple) and *Havcr2*<sup>icKO</sup>:5xFAD (yellow) phenotype. Individual cells of each genotype are represented by the bars on top. 
* **Fig. 6d-e** (`Fig.6d_6e_snRNAseq_microglia_P1_P2_score.R`): Violin plot of signature scores (Y-axis) identifying *Havcr2*<sup>icKO</sup>:5xFAD populations P1 (d) and P2 (e) across all genotypes (color) and microglia clusters 
* **Fig. 6f-h**: Projection of P1(magenta) and P2 (green) marker genes on the volcano plot of DEGs in different microglial perturbations: (f) Phagocytosis assay, (g) MGnD, (h) Tgfbr2cKO. Only marker genes which are also significant in each perturbation are indicated. i, Spearman correlation and significance of log2 fold changes of all expressed genes in each perturbation condition (f-h) with those in P1 vs P2. 
* **Fig. 6i** (`Fig.6i_snRNAseq_microglia_DAM_pseudobulk_corr.R`): Boxplot of pairwise spearman correlation coefficient (Y-axis) of whole gene expression pseudobulk microglial profiles between conditions (X-axis). Each point represents the spearman correlation coefficient between profiles of two biological replicates. Significance is computed by t-test.
* **Fig. 6j** (`Fig.6j_snRNAseq_microglia_dotplot_inflame_phago.R`): Dotplot representation of P1 and P2 marker genes that are significant DEGs in perturbations (f-h) and known anti-inflammatory, pro-inflammatory and phagocytic properties
* **Fig. 6k-l** (`Fig.6jk_6l_snRNAseq_microglia_signature_score_cluster2.R`): Violin plot of alternate macrophage (l) and phagocytic (m) signature genes among genotypes in cluster 2/MGnD/DAM. Significances are computed using Wilcoxon test.

### Extended Data Fig. 11: snRNA-seq data

* **Fig. 11c-e** (`Fig.s11c_s11e.snRNAseq_microglia_pct_barplot.R`): Barplot visualization of proportions (Y-axis) of individual mouse (c), and genotype (d) per cluster (X-axis); and cluster per mouse (X-axis) (e)
* **Fig. 11f-g** (`Fig.s11f_s11g.snRNAseq_microglia_signature_score.R`): Violin plot of the MGnD (f) and homeostasis (g) signature score (Y-axis) across microglia clusters (X-axis) and genotypes. 
* **Fig. 11h-i** (`Fig.s11h_s11i.snRNAseq_microglia_signature_score_cluster0.R`): Violin plot of scores for homeostatic (h) and *Tgfbr2*<sup>cKO</sup> (i) signature genes among genotypes in cluster 0. Significances are computed using Wilcoxon test.

### Extended Data Fig. 12: snRNA-seq data

* **Fig. 12a-b** (`Fig.s12a_s12b_snRNAseq_microglia_P1_P2_heatmap.R`): Heatmap of top DEGs distinguishing P1 (a) and P2 (b) subpopulations in snRNA-seq cluster 2/MGnD/DAM
* **Fig. 12c-d** (`Fig.s12c_s12d_snRNAseq_microglia_P1_P2_featureplot.R`): UMAP Visualization of microglial snRNA-seq clusters colored by P1 (c) and P2 (d) signature scores split by genotype

### Extended Data Fig. X 

* `Fig.CR14_snRNAseq_allcells_UMAP_dotplot.R`: snRNA-seq (1) UMAP of all cells colored by annotations (2) Dotplot of Havcr2 and Cx3cr1 expression in each cell type, split by genotype
