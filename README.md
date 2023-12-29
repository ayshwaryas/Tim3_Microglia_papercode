## Immune checkpoint molecule Tim-3 regulates microglial function and the development of Alzheimer’s disease pathology

### Processing and annotating snRNA-seq data

* `snRNAseq/snRNAseq_process.R`: QC, processing the snRNA-seq data; combine and subclustering microglia/PVM clusters
* `snRNAseq/snRNAseq_doublet_removal.R`: doublet removal using DoubletFinder
* `snRNAseq/snRNAseq_annotation.R`: annotating the snRNA-seq clusters

### Fig. 1 and Extended Fig. 3b

* `Fig.1a_s3b_microarray_Thion_2018.R`: Heatmap visualization of the developmental alterations of gene expressions in microglia in a published microarray dataset ([Thion et al. 2018, PMID: 29275859](https://pubmed.ncbi.nlm.nih.gov/29275859/))
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

* **Fig. 4a** (`Fig.4a_overlapping_DEGs_pairwise.R`): Number of overlapped genes and permutation test p-values between each pair of DEGs up- and down-regulated in *Havcr2*<sup>cKO</sup>, phagocytosing, *Tgfbr2*<sup>cKO</sup, and *Clec7a*<sup>+</sup> microglia compared to control microglia
* **Fig. 4b** (`Fig.4b_overlapping_DEGs_heatmap.R`): Heatmap visualization of correlation between vectors of the expression levels of DEGs of *Havcr2*<sup>cKO</sup> microglia compared to control microglia, as well as the DEGs shared by all three other comparisons (phagocytosing, *Tgfbr2*<sup>cKO</sup>, *Clec7a*<sup>+</sup> microglia compared to their corresponding control groups).

### Extended Data Fig.X 

* `Fig.CR14_snRNAseq_allcells_UMAP_dotplot.R`: snRNA-seq (1) UMAP of all cells colored by annotations (2) Dotplot of Havcr2 and Cx3cr1 expression in each cell type, split by genotype
