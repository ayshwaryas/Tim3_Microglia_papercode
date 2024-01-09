library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(ggtext)
library(parallel)
library(ComplexHeatmap)
library(ggpubr)
library(rstatix)
library(effsize)

# 1. Load data -----------------------------------
## Load Seurat object -------------------------
load("R_objects/2023-10-03.scRNAseq_cortex_MG_updated_clusters.RData")

## Signatures ---------------------------------

## Function for converting human gene symbol to mouse gene symbol
mouse_human_genes_raw <- read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
convert_human_to_mouse <- function(gene_list){
  mouse_human_genes_raw %>%
    filter(DB.Class.Key %in% subset(mouse_human_genes_raw, Symbol %in% gene_list )$DB.Class.Key) %>% 
    filter(Common.Organism.Name == "mouse, laboratory") %>%
    filter(!is.na(Symbol)) %>%
    pull(Symbol) %>% unique
}


## mouse MHCI genes
mouse_MHCI_genes <- mouse_human_genes_raw %>%
  filter(DB.Class.Key %in% subset(
    mouse_human_genes_raw, str_detect(Symbol, "HLA-[ABCEFGKL].*") )$DB.Class.Key) %>% 
  filter(Common.Organism.Name == "mouse, laboratory") %>%
  filter(!is.na(Symbol)) %>%
  pull(Symbol) %>% unique

## mouse MHCII genes
mouse_MHCII_genes <- mouse_human_genes_raw %>%
  filter(DB.Class.Key %in% subset(
    mouse_human_genes_raw,
    Symbol %in% c("HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1",
                  "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2", 
                  "HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5"))$DB.Class.Key) %>%
  filter(Common.Organism.Name == "mouse, laboratory") %>%
  filter(!is.na(Symbol)) %>%
  pull(Symbol) %>% unique()

### (1) cGAS-STING signatures ------------------------------------
### cGAS-STING signatures: Top 100 DEGs of CgasWT/R241E vs. CgasWT/WT mice, ordered by p-value
#### Gulen et al. 2023 (PMID: 37532932) Supplementary Table 6
cGAS_TableS6_sheets <- readxl::excel_sheets("data/cGAS_TableS6.xlsx")
cGAS_TableS6 <- lapply(cGAS_TableS6_sheets, function(x) {
  readxl::read_xlsx("data/cGAS_TableS6.xlsx", sheet = x) }) %>%
  `names<-`(str_replace_all(cGAS_TableS6_sheets, "\\-| ", "_") %>%
              str_replace("\\+", "."))

cGAS_signatures <- lapply(names(cGAS_TableS6), function(x) {
  cGAS_TableS6[[x]] %>%
    `colnames<-`(c("gene", colnames(.)[-1])) %>%
    mutate(FDR = p.adjust(.$p_val, method = "BH")) %>%
    mutate(direction = ifelse(avg_log2FC > 0, paste0(x, ".up"), paste0(x, ".down"))) %>%
    filter(FDR < 0.05) %>%
    group_by(direction) %>%
    arrange(p_val) %>%
    slice(1:100)}) %>%
  do.call(rbind, .) %>%
  split(f = .$direction) %>%
  lapply(pull, gene)


### (2) Interferon alpha response pathway ----------------------------
# MSigDB Hallmark interferon alpha signaling pathway (https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/HALLMARK_INTERFERON_ALPHA_RESPONSE.html)
Hallmark_IFNA <- read.table("data/signatures/HALLMARK_INTERFERON_ALPHA_RESPONSE.txt")$V1

### (3) TNF signaling pathways ---------------------------------------
# MSigDB Hallmark TNFA signaling pathway (https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/HALLMARK_TNFA_SIGNALING_VIA_NFKB.html)
Hallmark_TNFA <- read.table("data/signatures/HALLMARK_TNFA_SIGNALING_VIA_NFKB.txt")$V1

### KEGG TNF signaling pathway (mmu04668)
KEGG_TNF <- read.table("data/signatures/KEGG_TNF_Signaling.txt", sep = "\t")$V2 %>% str_extract("^.*(?=;)")

### (4) Phagosome pathways removing MHC genes ------------------------------
# KEGG phagosome pathway (mmu04145)
KEGG_Phagosome_rm_MHC <- read.table("data/signatures/KEGG_phagosome.txt", sep = "\t")$V2 %>%
  str_extract("^.*(?=;)") %>%
  setdiff(c(mouse_MHCI_genes, mouse_MHCII_genes)) 

# GO phagocytosis pathway (https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/GOBP_PHAGOCYTOSIS.html)
GO_Phagosome_rm_MHC <- read.table("data/signatures/GO_phagocytosis.txt")$V1 %>%
  setdiff(c(mouse_MHCI_genes, mouse_MHCII_genes))

### (5) Lysosome pathway --------------------------------------------------------
# KEGG lysosome pathway (mmu04142)
KEGG_Lysosome <- read.table("data/KEGG_lysosome.txt", sep = "\t")$V2 %>% str_extract("^.*(?=;)")

### (6) PI3K AKT pathway -------------------------------------------------------
# MSigDB Hallmark PI3K/AKT/mTOR signaling pathway (https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/HALLMARK_PI3K_AKT_MTOR_SIGNALING.html)
Hallmark_PI3K_AKT_MTOR <- read.table("data/HALLMARK_PI3K_AKT_MTOR_SIGNALING.txt")$V1

# KEGG PI3K-Akt signaling pathway (mmu04151)
KEGG_PI3K_AKT <- read.delim("data/KEGG_PI3K_AKT.txt", header = FALSE)$V2 %>% str_extract("^.*(?=;)")

### (7) MTORC1 signaling pathway -----------------------------------------------
# MSigDB Hallmark MTORC1 signaling pathway (https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/HALLMARK_MTORC1_SIGNALING.html)
Hallmark_MTORC1 <- read.table("data/HALLMARK_MTORC1_SIGNALING.txt")$V1

### (8) Hypoxia signaling pathway -----------------------------------------------
# MSigDB Hallmark Hypoxia signaling pathway (https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/HALLMARK_HYPOXIA.html)
Hallmark_Hypoxia <-  read.table("data/HALLMARK_HYPOXIA.txt")$V1



# 2. Compute signature score ------------------------------------------------------

## Compile signatuers into a list
signatures <- list(
  cGAS_STING_Up = cGAS_signatures$All_MG.up,
  cGAS_STING_Down = cGAS_signatures$All_MG.down,
  Hallmark_IFNA = Hallmark_IFNA,
  Hallmark_TNFA = Hallmark_TNFA,
  KEGG_TNF = KEGG_TNF,
  KEGG_Phagosome_rm_MHC = KEGG_Phagosome_rm_MHC,
  GO_Phagosome_rm_MHC = GO_Phagosome_rm_MHC,
  KEGG_Lysosome = KEGG_Lysosome,
  Hallmark_PI3K_AKT_MTOR = Hallmark_PI3K_AKT_MTOR,
  KEGG_PI3K_AKT = KEGG_PI3K_AKT,
  Hallmark_MTORC1 = Hallmark_MTORC1,
  Hallmark_Hypoxia = Hallmark_Hypoxia
)

for(x in names(signatures)) {
  scRNAseq_cortex_MG <- AddModuleScore(
    scRNAseq_cortex_MG, list(signatures[[x]]), name = x)  }


# 3. T-test on signatures scores, Havcr2icKO;5xFAD vs 5xFAD in each cluster ----------------
signatures_ttest <- sapply(names(signatures), function(x) {
  FetchData(scRNAseq_cortex_MG, c(paste0(x, 1), "updated_clusters", "Genotype")) %>%
    mutate(Genotype = factor(Genotype, c("Tim3icKO_5XFAD", "5XFAD"))) %>%
    split(f = .$updated_clusters) %>%
    lapply(function(df) {
      res <- t.test(as.formula(paste0(paste0(x, 1), "~", "Genotype")), data = df)
      data.frame(cluster = df$updated_clusters[1], 
                 statistic = res$statistic, p_value = res$p.value)}) %>%
    do.call(rbind, .) %>%
    mutate(padj = p.adjust(.$p_value, method = "BH"))
}, USE.NAMES = TRUE, simplify = FALSE)

## Convert FDRs to significance symbols and convert the data frame to wide format
signatures_ttest_signif_wide <- signatures_ttest %>%
  data.table::rbindlist(idcol = "signature") %>% 
  rstatix::add_significance(p.col = "padj",
                            symbols = c("****", "***", "**", "*", "")) %>%
  select(signature, cluster, padj.signif) %>%
  pivot_wider(names_from = cluster, values_from = padj.signif) %>%
  filter(!str_detect(signature, "Clec7a")) %>%
  mutate(signature = str_remove(signature, "_Signaling")) %>%
  column_to_rownames("signature")

# 4. Compute Hedges' g effect size on signatures scores, Havcr2icKO;5xFAD vs 5xFAD in each cluster ----------------
signatures_hedges_g <- sapply(names(signatures), function(x) {
  FetchData(scRNAseq_cortex_MG, c(paste0(x, 1), "updated_clusters", "Genotype")) %>%
    mutate(Genotype = factor(Genotype, c("Tim3icKO_5XFAD", "5XFAD"))) %>%
    split(f = .$updated_clusters) %>%
    lapply(function(df) {
      res <- cohen.d(f = df$Genotype, d = df[[paste0(x, 1)]], hedges.correction=TRUE)
      data.frame(cluster = df$updated_clusters[1], estimate = res$estimate)}) %>%
    do.call(rbind, .)
}, USE.NAMES = TRUE, simplify = FALSE)

## convert to wide format
signatures_hedges_g_wide <- signatures_hedges_g %>%
  data.table::rbindlist(idcol = "signature") %>%
  select(signature, cluster, estimate) %>%
  pivot_wider(names_from = cluster, values_from = estimate) %>%
  mutate(signature = str_remove(signature, "_Signaling")) %>%
  filter(!str_detect(signature, "Clec7a")) %>%
  column_to_rownames("signature")


all(rownames(signatures_hedges_g_wide) == rownames(signatures_ttest_signif_wide))


# Fig. 7b: Heatmap displaying signature score comparison between genotype (Havcr2icKO;5xFAD vs 5xFAD) ---------------
## Color: Hedges'g effect size
## Text:  t-test p-values adjusted using Benjamini Hochberg method

## Legend showing the significance levels
lgd_signif = Legend(
  pch = c("****", "*** ", "**  ", "*   "), title = "FDR-adjusted\nt-test p-value", type = "points", 
  labels = c(" <0.0001", " <0.001", " <0.01", " <0.05"), 
  size = unit(0.02, "mm"), background = "white")


p_hedges_g_htmap <- Heatmap(
  signatures_hedges_g_wide,
  col = circlize::colorRamp2(
    c(-0.56, -0.3, 0, 0.3, 0.56), 
    c("#2171B5", "#9ECAE1", "white", "#FC9272", "#CB181D")),
  border_gp = gpar(col = "black"),
  rect_gp = gpar(col = "white"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(signatures_ttest_signif_wide[i, j], x, y,
              gp = gpar(fontsize = 13), vjust = 0.8) },
  name = "Hedges' g\nEffect Size",
  column_names_rot = 45,
  cluster_rows = FALSE,
  cluster_column_slices = FALSE,
  cluster_columns = FALSE,
  column_gap = unit(1.5, "mm"),
  row_labels = str_replace(rownames(signatures_hedges_g_wide), "rm", "w/o"),
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10), 
  column_split = factor(str_remove(levels(signatures_hedges_g[[1]]$cluster), "_([0-6]|HMG|DAM)$"),
                        c("HMG", "DAM", "CC_HMG", "CC_DAM", "IFN")),
  column_title_gp = gpar(fontsize = 12,fontface = 2)
)

png("figures/signatures_hedges_g_htmap_RdBu.png", width = 10, height = 4, res = 400, units = "in")
draw(p_hedges_g_htmap, padding = unit(c(0, 6, 1, 1), "mm"),
     annotation_legend_list = list(lgd_signif), merge_legend = TRUE)
dev.off()

cairo_pdf("figures/signatures_hedges_g_htmap_RdBu.pdf", width = 10, height = 4)
draw(p_hedges_g_htmap, padding = unit(c(0, 6, 1, 1), "mm"),
     annotation_legend_list = list(lgd_signif), merge_legend = TRUE)
dev.off()



