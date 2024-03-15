library(tidyverse)
library(Seurat)
library(harmony)
library(parallel)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(cowplot)
library(ggpubfigs)
library(DoubletFinder)
library(ggtext)

# palette ----------------------------------------------------------------------
c25 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00","black", "gold1", 
         "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", "gray70", "khaki2", 
         "maroon", "orchid1","deeppink1","blue1","steelblue4","darkturquoise", "green3", 
         "yellow4", "yellow2", "#FEC44F", "brown4","darkcyan", "darkgoldenrod2", 
         "darkolivegreen4", "cornsilk", "aquamarine1", "darkseagreen1", "honeydew3", "salmon", 
         "slateblue3", "thistle2", "slategray2", "springgreen", "rosybrown1", "plum2", 
         "chocolate", "blanchedalmond", "brown1", "aliceblue", "cadetblue1", "coral1", 
         "firebrick1", "darkslategrey", "darkslateblue", "#745745", "#a01b68", "#ffc372",
         "#fff0f0", "#ebd4d4", "#eeeeee", "#de4463", "#e8ffc1", "#19d3da", "#ee6f57", "#ff9642",
         "#646464", "#0072ce", "#dc9cbf", "#ED8E7C", "#F1ECC3", "#D9DD6B", "#7C83FD", "#FDD2BF",
         "#E98580", "#C6B4CE", "#FFDADA", "#B5EAEA", "#BB8760", "#C9E4C5", "#FAEBE0", "#F7DBF0",
         "#66DE93", "#FFC074", "#B2B8A3", "#F2F4C3", "#D1D9D9", "#A7D0CD", "#114E60", "#28B5B5")

# Figure theme -----------------------------------------------------------------
theme_custom <- cowplot::theme_cowplot(font_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank())
theme_set(theme_custom)

# Function for processing snRNA-seq data ---------------------------------------
scrna_process <- function(seurat_obj, npc, res, normalize = TRUE, 
                          run_harmony = FALSE, max.iter.harmony = 25) {
  if(normalize) {
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  }
  
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj)) %>% 
    RunPCA(features = VariableFeatures(object = seurat_obj), npc = npc) 
  
  if(run_harmony) {
    seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "orig.ident", 
                             max.iter.harmony = max.iter.harmony)
    seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:npc) 
    seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:npc) 
    seurat_obj <- FindClusters(seurat_obj, resolution = res)
  }
  else {
    seurat_obj <- seurat_obj %>%
      FindNeighbors(dims = 1:npc) %>% 
      FindClusters(resolution = res) %>% 
      RunUMAP(dims = 1:npc)
  }
  return(seurat_obj)
}


# Signature list ---------------------------------------------------------------
## * Microglia: Known microglial genes, Van Hove et al. 2019, Monaghan et al. 2019, Zeisel et al. 2018
## * Border-associated macrophage (BAM): Van Hove et al. 2019
## * Perivascular macrophages (PVMs): Zeisel et al. 2018, Yang et al. 2019
## * Monocyte/monocyte-derived antigen-presenting cells: Monaghan et al. 2019
## * Infiltrating monocyte: https://www.biorxiv.org/content/10.1101/2021.05.30.446342v1

sig_list <- list(
  MG_sanity_check = c("Havcr2", "Tmem119", "Hexb", "P2ry12", "Sall1", "Mafb", "Spi1",
                      "Cx3cr1", "Trem2", "C1qc", "Cst7", "Apoe", "Axl", "Itgax", "Lilrb4",
                      "Ly9", "Il1b", "Tnf", "Lag3", "Cd274", "Pdcd1"),
  Microglia_Monaghan_2019 = c("Tmem119", "Sall1", "Hexb", "Slc2a5", "P2ry12",
                              "Siglech", "Trem2", "Bhlhe41", "Gpr34", "Serpine2"),
  Microglia_MouseBrain = c("P2ry12", "Ccr6", "Tmem119", "Casp8", "Ccl4", "Tnf"),
  Microglia_VanHove_2019 = c("Sall1", "Adgrg1", "Il1a", "Golm1", "Adamts1", "Selplg", "Sparc", 
                             "Gem", "Serpine2", "P2ry12", "Fscn1", "Rtn4rl1", "Atp6v0a2", 
                             "Arhgap5", "Siglech", "Ctsl", "Smad7", "P2ry13", "Tmem119", 
                             "Fam102b", "Smap2", "Csf3r", "Soga1", "Kcnk6", "Wsb1", "Rgs1", "Ldhb", 
                             "Gpr34", "Glul", "Cst3", "Lhfpl2", "Epb41l2", "Gpr155", "Rasal3", "Dtx4", 
                             "Srgap2", "Frmd4a", "Lpcat2", "Spata13", "Qk", "Entpd1", "Hexb", "Vsir", 
                             "Tanc2", "Bin1", "Trem2", "Olfml3", "Ctsf", "Elmo1", "Asph", "Gna15", "Cmtm6", 
                             "Herpud1", "Brd2", "Pde3b", "Tgfbr1", "Daglb", "Nfkbia", "Tmem173", "Basp1", "Arsb", 
                             "Inpp5d", "Fmnl3", "Slc29a3", "Sipa1", "Cx3cr1", "Snx18", "F11r", "Itgb5", "Fcrls"),
  BAM_signature = c("C3ar1", "Nrp1", "Cd63", "Ms4a7", "Adrb2", "Ms4a6b", "Ms4a6c", "Ifnar1",
                    "Ptger4", "Lifr", "Clec12a", "Zfp36l2", "Ehd4", "Zfp36l1", "Myo5a", 
                    "Swap70", "Rasa4", "B3galnt1", "St8sia4", "H1f0", "Sesn1", "Apoe", 
                    "Cd14", "Sdc4", "Aoah", "Pla2g7", "Tgfbi", "Ier2", "Maf"),
  PVM = c("Cd163", "Mrc1", "Lyve1", "Pf4", "Cd74", "Cxcl2"),
  Monocyte_and_derived_APCs = c("Ly6c2", "Ccr2", "Cd44", "Fcgr1", "Plac8", "Nr4a1"),
  Infiltrating_monocyte = c("Nupr1", "Adgre4", "Plac8", "Ifitm6", "Smpdl3a",
                            "Vim", "Gpx1", "Napsa", "S100a6", "S100a4", "S100a11",
                            "Crip1", "Ear2", "Fabp4", "Chil3", "Ccr2")
)

sig_title <- setNames(c("Microglia Sanity Check", "Microglia (Monaghan et al. 2019)",
                        "Microglia (Mouse Brain)", "Microglia (Van Hove et al. 2019)"),
                      "BAM Signature", "PVM Signature",
                      "Monocyte/Monocyte-Derived APCs", "Infiltrating Monocyte",
                      names(sig_list))

