## Immune checkpoint molecule Tim-3 regulates microglial function and the development of Alzheimer’s disease pathology

### Analysis of bulk RNA-seq data (📁[bulkRNAseq](/bulkRNAseq))
* Functions used in the analysis: `0.bulkRNAseq_functions.R`
* In-house generated datasets
    * **Dataset 1 Batch 1 (3-month-old mice)** (`1.bulkRNAseq_dataset1_batch1_3month.R`): control phagocytosing vs. <i>Havcr2</i><sup>cKO</sup> phagocytosing vs. control non-phagocytosing vs. <i>Havcr2</i><sup>cKO</sup> non-phagocytosing microglia
    * **Dataset 1 Batch 2 (7-month-old mice)** (`2.bulkRNAseq_dataset1_batch2_7month.R`): control vs. <i>Havcr2</i><sup>cKO</sup> vs. 5xFAD vs. <i>Havcr2</i><sup>cKO</sup>;5xFAD microglia
    * **Dataset 2 (1-month-old mice)** (`3.bulkRNAseq_dataset2_1month.R`): <i>Havcr2</i><sup>flox/flox</sup> vs. <i>Havcr2</i><sup>cKO</sup>
    * **Dataset 3 (4-month-old mice)** (`4.bulkRNAseq_dataset3_4month.R`): control vs. <i>Havcr2</i><sup>cKO</sup> vs. 5xFAD vs. <i>Havcr2</i><sup>cKO</sup>;5xFAD microglia
* Public datasets
    * [**Lund et al. 2018 (PMID: 29662171)**](https://pubmed.ncbi.nlm.nih.gov/29662171/) (`5.bulkRNAseq_Tgfbr2cKO_vs_control_Lund_2018.R`):  <i>Tgfbr2</i><sup>cKO</sup> vs. control microglia
    * [**Krasemann et al. 2017 (PMID: 28930663)**](https://pubmed.ncbi.nlm.nih.gov/28930663/) (`6.bulkRNAseq_Clec7apos_vs_neg_Krasemann_2017.R`):  <i>Clec7a</i><sup>+</sup> vs. <i>Clec7a</i><sup>-</sup> microglia

### Processing and annotating snRNA-seq data (📁[snRNAseq](/snRNAseq))

* `0.snRNAseq_functions.R`: functions, color palettes and gene signatures used in the analysis
* `1.snRNAseq_process.R`: QC, processing the snRNA-seq data; combining and subclustering microglia and PVM clusters
* `2.snRNAseq_doublet_removal.R`: doublet removal using DoubletFinder
* `3.snRNAseq_annotation.R`: annotating the snRNA-seq clusters
* `4.snRNAseq_subcluster_MG.R`: further subclustering microglia cells; splitting *Havcr2*<sup>icKO</sup> 5XFAD nuclei in cluster 2 (DAM/MGnD) into subpopulations P1 and P2 based on [MSigDB Hallmark TGFβ signaling](https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/HALLMARK_TGF_BETA_SIGNALING.html) signature score


### Processing and annotating scRNA-seq data (📁[scRNAseq](/scRNAseq))

* `0.scRNAseq_functions.R`: functions, color palettes and gene signatures used in the analysis
* `1.scRNAseq_process.R`: QC and processing of the scRNA-seq data
* `2.scRNAseq_subcluster.R`: subclustering cluster 6, 7, 8 obtained in 1.
* `3.scRNAseq_DoubletFinder.R`: doublet identification using DoubletFinder
* `4.scRNAseq_annotation.R`: identify MGnD, Homeostasis, interferon-rich and cycling populations based on signature score; combining and reclustering the cells based on the population characteristics; finalizing annotations
* `5.scRNAseq_DGE_analysis.R`: differential gene expression analysis (Wilcoxon) comparing <i>Havcr2</i><sup>icKO</sup>;5xFAD and 5xFAD cells in each cluster using the `FindMarkers` function

### Processing and annotating public AD sc/snRNA-seq datasets (📁[public_AD](/public_AD))
* `1.GSE140510_snRNAseq_7M.R`: processing 7-month snRNA-seq dataset from GSE140510
* `2.GSE98969_scRNAseq_6M.R`: processing 6-month scRNA-seq dataset from GSE98969

### [Fig. 1](/Figures/Fig.1)
* **Fig. 1a** (`Fig.1a_microarray_Thion_2018.R`): Heatmap visualization of the developmental alterations of gene expressions of immune checkpoints and TGFβ pathway-related molecules in microglia in a published microarray dataset ([Thion et al. 2018, PMID: 29275859](https://pubmed.ncbi.nlm.nih.gov/29275859/))
* **Fig. 1b** (`Fig.1b_GSE127449_P9_P28_scRNAseq.R`): Dotplot visualization of gene expression at postnatal stages (P9 and P28) in a published scRNA-seq dataset ([He et al. 2018, PMID: 34982959](https://pubmed.ncbi.nlm.nih.gov/34982959/); [GSE127449](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127449))
* **Fig. 1h** (`Fig.1h_1month_heatmap_DEGs.R`): Heatmap visualization of DEGs (FDR < 0.05) of 1-month-old *Havcr2*<sup>cKO</sup> mice (n = 4 (3 males, 1 female)) compared to *Havcr2*<sup>flox/flox</sup> mice (n = 5 (4 males, 1 female)). Rows represent genes and columns are biological replicates. Gene expression is row-normalized across replicates. 
* **Fig. 1i** (`Fig.1i_1month_volcanoplot.R`): Volcano plot of differential gene expression analysis performed by DESeq2 comparing 1-month-old *Havcr2*<sup>cKO</sup> to *Havcr2*<sup>flox/flox</sup> mice
* **Fig. 1j-l** (`Fig.1j_1l_1month_signature_score.R`): Boxplot visualization of the the scores of MGnD (j), homeostasis (k), and KEGG phagosome (l) signatures in 1-month-old *Havcr2*<sup>cKO</sup> compared to *Havcr2*<sup>flox/flox</sup> mice. The Y-axis represents log2 transformed average TPM, and significance was calculated using the Student's two-tailed t-test.

### [Fig. 2](/Figures/Fig.2)
* **Fig. 2a** (`Fig.2a_1month_Tgfbr2_circos_overlap.R`): Circos plot comparison of up- and down-regulated DEGs in 
    (1) <i>Havcr2</i><sup>cKO</sup> vs. control microglia from 1-month-old mice
    (2) <i>Tgfbr2</i><sup>cKO</sup> vs. control microglia
* **Fig. 2b** (`Fig.2b_1month_Tgfbr2_scatter_plot.R`): Scatter plot of genes based on expressional difference represented by log2-transformed fold changes in *Havcr2*<sup>cKO</sup> (X-axis) and *Tgfbr2*<sup>cKO</sup> (Y-axis) compared to control microglia
* **Fig. 2e** (`Fig.2e_1month_DEG_homer_motif.Rmd`): Transcription factor motif enrichment analysis using DEGs comparing *Havcr2*<sup>cKO</sup> to control microglia from 1-month-old mice

### [Fig. 3](/Figures/Fig.3)
* **Fig. 3g** (`Fig.3g_3month_heatmap_overlapping_DEGs.R`): Heatmap visualization of top DEGs shared by at least two of the following three comparisons (3-month-old mice)
    (1) control phagocytosing versus control non-phagocytosing microglia
    (2) *Havcr2*<sup>cKO</sup> phagocytosing versus control phagocytosing microglia
    (3) *Havcr2*<sup>cKO</sup> non-phagocytosing versus control non-phagocytosing microglia
* **Fig. 3k** (`Fig.3k_3month_nonphago_DEGs_pathway.R`):  Pathway analysis of DEGs in non-phagocytosing microglia from 3-month-old *Havcr2*<sup>cKO</sup> mice compared to control mice. Disease pathways (pathways under section 6.1-6.10 from https://www.kegg.jp/kegg/pathway.html) and ribosomal genes were excluded from the analysis.
* **Fig. 3i** (`Fig.3i_3month_nonphago_vs_Tgfbr2_circos_overlap.R`): Circos plot comparison of up- and down-regulated DEGs from 3-month-old mice in below comparisons 
    (1) <i>Havcr2</i><sup>cKO</sup> non-phagocytosing vs. control non-phagocytosing
    (2) control phagocytosing vs. control non-phagocytosing microglia

### [Fig. 4](/Figures/Fig.4)
* **Fig. 4a** (`Fig.4a_overlapping_DEGs_pairwise.R`): Number of overlapped genes and permutation test p-values between each pair of DEGs up- and down-regulated in *Havcr2*<sup>cKO</sup>, phagocytosing, *Tgfbr2*<sup>cKO</sup>, and *Clec7a*<sup>+</sup> microglia compared to control microglia
* **Fig. 4b** (`Fig.4b_overlapping_DEGs_corr_hetamap.R`): Heatmap visualization of correlation between vectors of the expression levels of DEGs of *Havcr2*<sup>cKO</sup> microglia compared to control microglia, as well as the DEGs shared by all three other comparisons (phagocytosing, *Tgfbr2*<sup>cKO</sup>, *Clec7a*<sup>+</sup> microglia compared to their corresponding control groups).

### [Fig. 6](/Figures/Fig.6): snRNA-seq data
* **Fig. 6a** (`Fig.6a_snRNAseq_microglia_UMAP.R`): UMAP visualization of microglia snRNA-seq clusters split by genotype
* **Fig. 6b** (`Fig.6b_snRNAseq_microglia_HallmarkTGFB_score.R`): Probability density curves of the signature score for the Hallmark TGFβ pathway in the 5xFAD (purple) and *Havcr2*<sup>icKO</sup>:5xFAD (yellow) phenotype. Individual cells of each genotype are represented by the bars on top. 
* **Fig. 6c-e**: Projection of P1(magenta) and P2 (green) marker genes on the volcano plot of DEGs in different microglial perturbations: (c) Phagocytosis assay, (d) MGnD, (e) <i>Tgfbr2</i><sup>cKO</sup>. Only marker genes which are also significant in each perturbation are indicated. i, Spearman correlation and significance of log2 fold changes of all expressed genes in each perturbation condition (f-h) with those in P1 vs P2. 
* **Fig. 6f** (`Fig.6f_snRNAseq_microglia_dotplot_inflame_phago.R`): Dotplot representation of P1 and P2 marker genes that are significant DEGs in perturbations (f-h) and known anti-inflammatory, pro-inflammatory and phagocytic properties
* **Fig. 6g-h** (`Fig.6g_6h_snRNAseq_microglia_signature_score_cluster2.R`): Violin plot of alternate macrophage (l) and phagocytic (m) signature genes among genotypes in cluster 2/MGnD/DAM. Significances are computed using one-way ANOVA with Tukey’s HSD correction.

### [Fig. 7](/Figures/Fig.7): scRNA-seq data
* **Fig. 7a** (`Fig.7a_scRNAseq_UMAP.R`): UMAP visualization of microglia scRNA-seq clusters split by genotype
* **Fig. 7b** (`Fig.7b_scRNAseq_signature_score_hedgesg.R`): Heatmap visualization of the Hedges' g effect size (color) and significance (text) comparing the signature scores in <i>Havcr2</i><sup>icKO</sup>;5xFAD and 5xFAD cells in each microglia cluster. Significance is computed by t-test and adjusted using the Benjamini-Hochberg procedures. Signatures include:
    * cGAS-STING signature: Top 100 up- and down-regulated DEGs ordered by p-value comparing Cgas<sup>WT/R241E</sup> to Cgas<sup>WT/WT</sup> mice from [Gulen et al. 2023 (PMID: 37532932)](https://www.nature.com/articles/s41586-023-06373-1)
    * [MSigDB Hallmark interferon alpha signaling](https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/HALLMARK_INTERFERON_ALPHA_RESPONSE.html)
    * [MSigDB Hallmark TNFA signaling](https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/HALLMARK_TNFA_SIGNALING_VIA_NFKB.html) and [KEGG TNF signaling](https://www.genome.jp/dbget-bin/www_bget?path:mmu04668) pathways
    * [KEGG phagosome](https://www.genome.jp/dbget-bin/www_bget?path:mmu04145) and [GO phagocytosis](https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/GOBP_PHAGOCYTOSIS.html) pathways, excluding MHC genes
    * [KEGG lysosome pathway](https://www.genome.jp/dbget-bin/www_bget?path:mmu04142)
    * [MSigDB Hallmark PI3K/AKT/mTOR signaling](https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/HALLMARK_PI3K_AKT_MTOR_SIGNALING.html) and [KEGG PI3K-Akt signaling](https://www.genome.jp/dbget-bin/www_bget?path:mmu04151) pathways
    * [MSigDB Hallmark MTORC1 signaling pathway](https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/HALLMARK_MTORC1_SIGNALING.html)
    * [MSigDB Hallmark hypoxia signaling pathway](https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/HALLMARK_HYPOXIA.html)
* **Fig. 7c** (`Fig.7c_scRNAseq_signature_score_ttest.R`): Heatmap visualization of T-test results comparing the scores of microglia population markers from [Ellwanger et al. 2021 (PMID: 33446504)](https://pubmed.ncbi.nlm.nih.gov/33446504/) in <i>Havcr2</i><sup>icKO</sup>;5xFAD and 5xFAD cells in each microglia cluster. T-statistics are displayed in text (black if FDR < 0.05; white if FDR $\geq$ 0.05). Colors represent whether the average signature scores are higher in <i>Havcr2</i><sup>icKO</sup>;5xFAD or 5xFAD (pink: higher in <i>Havcr2</i><sup>icKO</sup>;5xFAD; blue: higher in 5xFAD). 
* **Fig. 7d** (`Fig.7d_scRNAseq_heatmap_selected_DEGs.R`): Heatmap visualization of log2 fold-change and significances of selected DEGs comparing <i>Havcr2</i><sup>icKO</sup>;5xFAD to 5xFAD in each cluster.

### [Extended Data Fig. 2](/Extended_Data_Figures/Fig.s2)
* **Extended Data Fig. 2b-c** (`Fig.s2b_s2c_l5_mousebrain_scRNAseq.R`): analysis of a public scRNA-seq dataset of mouse nervous system ([Zeisel et al. 2018, PMID: 30096314](https://pubmed.ncbi.nlm.nih.gov/30096314/))
    * **Extended Data Fig. 2b**: The expression of immune checkpoint genes, including Havcr2, Lag3, and Vsir in mouse microglia and other cell populations in the central and peripheral nervous system
    * **Extended Data Fig. 2c**: Dissociation signature score in microglia and perivascular macrophages
* **Extended Data Fig. 2d** (`Fig.s2d_humanbrain_scRNAseq.R`): The expression of immune checkpoint genes and TGFβ pathway-related genes in human microglia and other cell populations in the brain ([Gaublomme et al. 2019, PMID: 31266958](https://pubmed.ncbi.nlm.nih.gov/31266958/))

### [Extended Data Fig. 3](/Extended_Data_Figures/Fig.s3)
* **Extended Data Fig. 3a** (`Fig.s3a_BrainSpan_bulk_RNAseq.R`): Heatmap visualization the developmental alterations of the expressions of selected checkpoint molecules in multiple brain regions in a published bulk RNA-seq dataset ([Miller et al. 2014, PMID: 24695229](https://pubmed.ncbi.nlm.nih.gov/24695229/)). The RPKM values were averaged across samples of the same age and then log2-transformed for visual clarity
* **Extended Data Fig. 3c** (`Fig.s3c_microarray_Thion_2018.R`): Heatmap visualization of the developmental alterations of the expressions of MGnD and homeostasis associated genes in microglia in a published microarray dataset ([Thion et al. 2018, PMID: 29275859](https://pubmed.ncbi.nlm.nih.gov/29275859/))

### [Extended Data Fig. 6](/Extended_Data_Figures/Fig.s6)
* **Extended Data Fig. 6a** (`Fig.s6a_circo_overlapping_DEGs.R`): Circos plot comparison of DEGs up- and down-regulated in *Havcr2*<sup>cKO</sup>, phagocytosing, *Tgfbr2*<sup>cKO</sup>, and *Clec7a*<sup>+</sup> microglia compared to control microglia. Only islands with DEGs shared by at least 3 of the 4 comparisons were included. Permutation test p-values were displayed as ns P > 0.025, \*P ≤ 0.025, \*\*P ≤ 0.01, \*\*\*P ≤ 0.001, \*\*\*\*P ≤ 0.0001, and \*\*\*\*\*P ≤ 0.00001.
* **Extended Data Fig. 6b** (`Fig.s6b_heatmap_overlapping_DEGs.R`): Heatmap visualization of the expressions of DEGs shared by at least 3 of the 4 aforementioned comparisons. 

### [Extended Data Fig. 12](/Extended_Data_Figures/Fig.s12)
* **Extended Data Fig. 12a-b** (`Fig.s12a_4month_F_tim3cKO_5XFAD_circos_overlap.R`, `Fig.s12b_7month_F_tim3cKO_5XFAD_circos_overlap.R`): Circos plot comparison of up- and down-regulated DEGs from 4- (a) and 7-month-old (b) female mice in the below comparisons
    (1) <i>Havcr2</i><sup>cKO</sup> vs. control microglia 
    (2) 5xFAD vs. control microglia

### [Extended Data Fig. 13](/Extended_Data_Figures/Fig.s13): snRNA-seq data
* **Extended Data Fig. 13a-c** (`Fig.s11c_s11e.snRNAseq_microglia_pct_barplot.R`): Barplot visualization of proportions (Y-axis) of individual mouse (a), and genotype (b) per cluster (X-axis); and cluster per mouse (X-axis) (c)
* **Fig. 13d** (`Fig.13d_snRNAseq_microglia_top5markers.R`): Dotplot visualization of top 5 genes distinguishing each cell cluster from Fig. 6a
* **Extended Data Fig. 13e-f** (`Fig.s13e_s13f_snRNAseq_allcells_UMAP_dotplot.R`): snRNA-seq UMAP of all cells colored by annotations (e) and Dotplot of Havcr2 and Cx3cr1 expression in each cell type, split by genotype (f)
* **Extended Data Fig. 13g-h** (`Fig.s13g_s13h.snRNAseq_microglia_signature_score.R`): Violin plot of the MGnD (g) and homeostasis (h) signature score (Y-axis) across microglia clusters (X-axis) and genotypes. 
* **Extended Data Fig. 13i** (`Fig.s13i.snRNAseq_microglia_proportion_test.R`): Proportions of microglial cluster 0-4 in 5xFAD and *Havcr2*<sup>icKO</sup>:5xFAD (n = 3 for each).
* **Extended Data Fig. 13j-k** (`Fig.s13j_s13k.snRNAseq_microglia_signature_score_cluster0.R`): Violin plot of scores for homeostatic (j) and *Tgfbr2*<sup>cKO</sup> (k) signature genes among genotypes in cluster 0. Significances are computed using Wilcoxon test.

### [Extended Data Fig. 14](/Extended_Data_Figures/Fig.s12): snRNA-seq data
* **Extended Data Fig. 14a-b** (`Fig.s14a_s14b_snRNAseq_microglia_P1_P2_heatmap.R`): Heatmap of top DEGs distinguishing P1 (a) and P2 (b) subpopulations in snRNA-seq cluster 2/MGnD/DAM
* **Extended Data Fig. 14c-d** (`Fig.s14c_s14d_snRNAseq_microglia_P1_P2_featureplot.R`): UMAP visualization of microglial snRNA-seq clusters colored by P1 (c) and P2 (d) signature scores split by genotype


### [Extended Data Fig. 15](/Extended_Data_Figures/Fig.s15): snRNA-seq data
* **Extended Data Fig. 15 a-h** (`Fig.s15a_s15h.snRNAseq_microglia_BAM_PVM_monocyte_signature_score.R`): Violin plot of border-associated macrophage (BAM, a), perivascular macrophage (PVM, b), monocyte (c-d) and microglia (e-h) signature scores across microglia clusters (X-axis) and genotypes. 
* **Extended Data Fig. s15i-j** (`Fig.s15i_s15j_snRNAseq_microglia_P1_P2_score.R`): Violin plot of signature scores (Y-axis) identifying *Havcr2*<sup>icKO</sup>:5xFAD populations P1 (i) and P2 (j) across all genotypes (color) and microglia clusters 

### [Extended Data Fig. 16](/Extended_Data_Figures/Fig.s16): snRNA-seq data
* **Extended Data Fig. 16a** (`Fig.s16a_snRNAseq_microglia_DAM_pseudobulk_corr.R`): Boxplot of pairwise spearman correlation coefficient (Y-axis) of whole gene expression pseudobulk microglial profiles between conditions (X-axis). Each point represents the spearman correlation coefficient between profiles of two biological replicates. Significance is computed by t-test.
* **Extended Data Fig. 16b** (`Fig.s16b_snRNAseq_microglia_human_MG_markers.R`): Twelve human microglial signatures curated using samples from patients with AD ([Sun et al. 2023, PMID: 37774678](https://pubmed.ncbi.nlm.nih.gov/37774678/)). These signature scores were plotted for cluster 2 in our snRNAseq data. Significance were computed using one-way ANOVA with Tukey's HSD correction.
* **Extended Data Fig. 16c-d** (`Fig.s16c_s16d_snRNAseq_microglia_P1_P2_vs_Hallmark_TGFB_sub2.R`): Scatter plot of P1 (Y-axis, c) or P2 (Y-axis, d), and Hallmark TGFB signature scores among 5xFAD and <i>Havcr2</i><sup>icKO</sup>;5xFAD microglia in current study
* **Extended Data Fig. 16e-g** (`Fig.s16e_s16g_public_snRNAseq_P1_P2_vs_Hallmark_TGFB_sub2.R`): Scatterplot of P1 (Y-axis) or P2 (Y-axis), and Hallmark TGFB signature scores among 5xFAD microglia in MGnD clusters in public datasets (GSE98969 cluster 2 (e), GSE98969 cluster 3 (f) , GSE140510 cluster 1 (g))


### [Extended Data Fig. 18](/Extended_Data_Figures/Fig.s18): scRNA-seq data
* **Extended Data Fig. 18a** (`Fig.s18a_scRNAseq_Clec7a_signature.R`): UMAP of the scRNA-seq data colored by MGnD (a) and Homeostatis (b) signature scores
* **Extended Data Fig. 18b** (`Fig.s18b_scRNAseq_dotplot_top10markers_by_log2FC.R`): Dotplot visualization of top 10 genes distinguishing each cell cluster
* **Extended Data Fig. 18c** (`Fig.s18c_scRNAseq_perc_barplot_by_genotype.R`): Barplot visualization of proportions (Y-axis) of genotype per cluster
* **Extended Data Fig. 18d-e** (`Fig.s18d_s18e_scRNAseq_scatter_IFNDAMvsHMG_preMGnD.R`): Scatter plot of genes based on expressional difference represented by log2-transformed fold changes of
    * IFN_DAM compared to IFN_HMG cluster in the scRNA-seq data (X-axis) and
    * late IFN-responsive pre-MGnD compared to early IFN-responsive pre-MGnD clusters from male (d) and female (e) mice defined in [Yin et al. 2023 (PMID: 37291336)](https://pubmed.ncbi.nlm.nih.gov/37291336/)

### [Extended Data Fig. 19](Fig.s19): scRNA-seq data
* **Extended Data Fig. 19a-b** (`Fig.s19a_s19b_scRNAseq_scatter_Havcr2iKO5xFADvs5XFAD_preMGnD.R`): Scatter plot of genes based on expressional difference represented by log2-transformed fold changes of
    * <i>Havcr2</i><sup>icKO</sup>;5xFAD vs 5xFAD in cluster HMG_0 (a) and DAM_0 (b)
    * <i>Clec7a</i><sup>-</sup> vs <i>Clec7a</i><sup>+</sup> (Y-axis)
* **Extended Data Fig. 19c-d** (`Fig.s19c_s19d_scRNAseq_volcano_Havcr2icKO5xFADvs5xFAD.R`): volcano plot of differential gene expression analysis comparing <i>Havcr2</i><sup>icKO</sup>;5xFAD and 5xFAD in cluster HMG_0 (c) and DAM_0 (d)
