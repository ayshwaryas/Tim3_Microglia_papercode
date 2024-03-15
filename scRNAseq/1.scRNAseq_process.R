source("scRNAseq/0.scRNAseq_functions.R")

# 1. Load Data & Create Seurat object --------------
path_10X <- "raw_data_no_intron"
samples <- list.files(path_10X)

scRNAseq_cortex <- lapply(samples, function(i) {
  obj <- Read10X(paste0(path_10X, i))
  obj <- CreateSeuratObject(counts = obj, project = i)
})

scRNAseq_cortex <- merge(
  x = scRNAseq_cortex[[1]], y = scRNAseq_cortex[-1], 
  add.cell.ids = samples, project = "scRNAseq_20230206")

## Add metadata
metadata <- readxl::read_xlsx("SK-5GF3.xlsx") %>%
  filter(!is.na(`Sample ID`))

metadata <- metadata %>% 
  mutate(Genotype = str_extract(Alias...3, "Tim3icKO_5XFAD|5XFAD")) %>%
  select(Sample_ID = "Alias...4", Genotype, Gender) %>%
  mutate(Tissue = str_extract(Sample_ID, "cortex|hippocampus")) %>%
  mutate(Mouse_ID = str_extract(Sample_ID, "(cKO|con)(1|2)")) %>%
  mutate(Mouse_ID = factor(Mouse_ID, c("con1", "con2", "cKO1", "cKO2")))

scRNAseq_cortex_meta <- scRNAseq_cortex@meta.data %>%
  rownames_to_column("Cells") %>%
  select(Cells, orig.ident) %>%
  left_join(metadata, by = c("orig.ident" = "Sample_ID")) %>%
  column_to_rownames("Cells") %>%
  select(-orig.ident)

scRNAseq_cortex <- AddMetaData(scRNAseq_cortex, scRNAseq_cortex_meta)


# 2. QC -------------------
scRNAseq_cortex[["percent.mt"]] <- PercentageFeatureSet(scRNAseq_cortex, pattern = "^mt-")
scRNAseq_cortex <- subset(scRNAseq_cortex, subset = nCount_RNA>= 1000 & nFeature_RNA>= 200 & percent.mt < 80)

# 3. Processing the scRNA-seq data ------------------
scRNAseq_cortex <- scrna_process(scRNAseq_cortex, npc = 50, res = 0.4, normalize = TRUE)
save(scRNAseq_cortex, file = "R_objects/scRNAseq_processed_no_intron_cortex_only.RData")

## ==> Cluster 5 is low-quality
VlnPlot(scRNAseq_cortex, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# 4. Scoring microglia, PVM, BAM, monocyte signatures ------------------------
## ==> Cluster 9 is potentially PVM
parallel::mclapply(names(sig_list), function(x){
  scRNAseq_cortex <- AddModuleScore(
    scRNAseq_cortex, list(sig_list[[x]]), name = x)
  
  VlnPlot(scRNAseq_cortex, features = paste0(x, 1), group.by = "seurat_clusters", 
          cols = c25, pt.size = 0.1) +
    labs(title = sig_title[x], x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  ggsave(paste0('figures/signatures_cortex/vln_', x, ".png"), width = 6, height = 3.2)
  ggsave(paste0('figures/signatures_cortex/vln_', x, ".pdf"), width = 6, height = 3.2)
  
  VlnPlot(scRNAseq_cortex, features = paste0(x, 1), split.by = "Genotype",
          group.by = "seurat_clusters", pt.size = 0.1) +
    labs(title = sig_title[x], x = NULL, y = NULL) +
    scale_fill_manual(breaks = c("5XFAD", "Tim3icKO_5XFAD"),
                      labels = c("5XFAD", "<i>Havcr2</i><sup>icKO</sup> 5XFAD"),
                      values = c("#785EF0", "#E69F00")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.text = element_markdown())
  ggsave(paste0('figures/signatures_cortex/vln_', x, "_by_genotype.png"), width = 9, height = 3.2)
  ggsave(paste0('figures/signatures_cortex/vln_', x, "_by_genotype.pdf"), width = 9, height = 3.2)
  
})

