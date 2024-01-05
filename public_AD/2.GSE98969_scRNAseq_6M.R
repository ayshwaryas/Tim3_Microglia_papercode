source("snRNAseq/0.snRNAseq_functions.R")

# 1. Creating and merging Seurat objects ---------------------------------
GSE98969_expr_design <- read.delim("data/GSE98969/GSE98969_experimental_design_f.txt", sep = "\t", skip = 16) %>%
  mutate(Cell = paste(Amp_batch_ID, Well_ID, sep = "_")) %>%
  mutate(Genotype = ifelse(Mouse_ID == "C57BL/6", "WT", Mouse_ID)) %>%
  mutate(Genotype =  factor(Genotype, c("WT", "Trem2KO", "5XFAD", "Trem2KO_5XFAD", "SOD1"))) %>%
  mutate(Tissue = str_extract(Batch_desc, "cerebellum|cortex")) 

sample_name_6M <- unique(subset(GSE98969_expr_design, Seq_batch_ID == "SB126")$Amp_batch_ID)

GSE98969_obj <- lapply(sample_name_6M, function(i) {
  expr <- read.table(paste0("data/GSE98969/expr/", i))
  CreateSeuratObject(counts = expr, project = str_extract(i,  "AB[0-9]+"))
}) %>% `names<-`(names(sample_name_6M))

GSE98969_obj <- merge(
  x = GSE98969_obj[[1]],
  y = GSE98969_obj[-1],
  add.cell.ids = names(sample_name_6M),
  project = "GSE98969"
)

# 2. QC ----------------------------------------------------------------------
GSE98969_obj[["percent.mt"]] <- PercentageFeatureSet(GSE98969_obj, pattern = "^mt", assay = "RNA")
GSE98969_obj <- subset(GSE98969_obj, subset = nCount_RNA>= 1000 & nFeature_RNA>= 200 & percent.mt < 80)

# 3. Processing --------------------------------------------------------------
GSE98969_obj <- scrna_process(GSE98969_obj, npc = 50, res = 1, normalize = TRUE)


GSE98969_meta <- data.frame(Cell = Cells(GSE98969_obj)) %>%
  left_join(GSE98969_expr_design, by = "Cell") %>%
  column_to_rownames("Cell")
GSE98969_obj <- AddMetaData(GSE98969_obj, GSE98969_meta)
GSE98969_obj$Mouse_ID <- str_extract(GSE98969_obj$Batch_desc, "mouse[1-2]") %>%
  paste(GSE98969_obj$Genotype, ., sep = "_") %>%
  factor(c("WT_mouse1", "WT_mouse2", "5XFAD_mouse1", "5XFAD_mouse2"))

save(GSE98969_obj, file = "R_objects/GSE98969_ADmouse_6M.RData")


# 4. Identify microglia clusters ----------------------------
dotplot_markers(GSE98969_obj, width = 5, height = 10, folder = "GSE98969")


# 5. Identify DAM/MGnD clusters -----------------------------
## ==> Subcluster 2 and 3 are DAM/MGnD
## MGnD and Homeostasis signature (Top 100 DEGs up/down in Clec7a+ vs Clec7a-)

Clec7a_sig <- read.csv("data/Oleg_Immunity_2017_DEG_ADpos_vs_neg_all.csv") %>%
  arrange(padj) %>% group_by(direction) %>% 
  dplyr::slice(1:100) %>%
  split(f = .$direction) %>%
  lapply(pull, tracking_id)

## Compute MGnD and Homeostasis signature score
GSE98969_obj <- AddModuleScore(GSE98969_obj, list(Clec7a_sig$up), name = "MGnD")
GSE98969_obj <- AddModuleScore(GSE98969_obj, list(Clec7a_sig$down), name = "Homeostasis")

for(x in c("MGnD", "Homeostasis")) {
  FeaturePlot(GSE98969_obj, features = paste0(x, 1), cols = rev(brewer.pal(11, "Spectral")),
              label = TRUE, label.size = 3, repel = TRUE) + labs(title = x)
  ggsave(paste0("figures/GSE98969/sig_score_GSE98969_MG_", x, ".png"), width = 5, height = 3.5)
  
  VlnPlot(GSE98969_obj, features = paste0(x, 1), cols = c25) + NoLegend()
  ggsave(paste0("figures/GSE98969/vln_GSE98969_MG_", x, ".png"), width = 4, height = 2.5)
  
  VlnPlot(GSE98969_obj, features = paste0(x, 1), cols = c25, split.by = "Genotype") 
  ggsave(paste0("figures/GSE98969/vln_GSE989690_MG_", x, "_by_genotype.png"), width = 5.5, height = 2.5)
}

