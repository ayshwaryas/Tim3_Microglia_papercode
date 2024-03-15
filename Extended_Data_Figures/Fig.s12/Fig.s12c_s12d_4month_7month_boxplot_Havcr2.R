library(tidyverse)

# Load data --------------------------------------------------------------------
## Dataset 1 Batch 2, 7 month
outliers_batch2 <- c("KK1", "KK6", "KK7")

tpm_ds1_batch2 <- read.csv("data/expr_mat/tim3_Kimi_tpm_with_genesymbols_KK.csv") %>% 
  select(gene_name, rownames(meta_ds1_batch2))

meta_ds1_batch2 <- read.csv("data/expr_mat/tim3_Kimi_metadata_KK.csv") %>%
  filter(Batch == "Batch_2") %>%
  filter(Sample_ID %in% colnames(tpm_ds1_batch2)[-(1:2)]) %>%
  filter(Batch == "Batch_2" & !(Sample_number %in% outliers_batch2)) %>%
  mutate(genotype = factor(genotype, c("control", "Tim3.cKO", "5XFAD", "Tim3.cKO.5XFAD")))

tpm_ds1_batch2_long <- tpm_ds1_batch2 %>% 
  pivot_longer(2:ncol(.), values_to = "TPM", names_to = "Sample_ID") %>%
  left_join(meta_ds1_batch2, by = "Sample_ID")

## Dataset 3, 4 month
metadata_ds3 <- read.csv("data/expr_mat/GET_metadata.csv") %>%
  mutate(Cohort = str_extract(Cohort, "[0-9]+th"),
         Cohort = factor(Cohort, paste0(8:10, "th"))) %>%
  filter(!Plate %in% c("C04", "E04") & !str_detect(Plate, "11")) %>%
  mutate(Plate_col = str_extract(Plate, "[0-9]+")) %>%
  filter(Genotype != "Tim3_half_cKO.5XFAD") %>%
  mutate(Group = case_when(Plate_col %in% c("06", "08", "09") ~ "Group_1",
                           Plate_col %in% c("03", "10", "12") ~ "Group_2",
                           Plate_col %in% c("01", "02", "04", "05", "07") ~ "Group_3")) %>%
  mutate(genotype = Genotype, sex = Sex)  %>%
  mutate(genotype = factor(genotype, c("control", "Tim3_cKO", "5XFAD", "Tim3_cKO.5XFAD")))

tpm_ds3 <- read.csv("data/expr_mat/GET_tpm.csv") %>%
  select(gene_name = gene_symbol, metadata_ds3$Sample_ID)


tpm_ds3_long <- tpm_ds3 %>% pivot_longer(2:ncol(.), values_to = "TPM", names_to = "Sample_ID") %>%
  left_join(metadata_ds3, by = "Sample_ID")

# Boxplot of Havcr2 expression ------------------------------------------------------
boxplot_genes <- function(genes, tpm, x_var = "sex", pal = brewer.pal(4, "Set2"), shape_var = "sex",
                          x_breaks = c("F", "M"), suffix = "", ncol = 3, width = 7, height = 7, log = TRUE) {
  if(log) {tpm <- tpm %>% mutate(expr = log2(TPM + 1))}
  else {tpm <- tpm %>% mutate(expr = TPM)}
  tpm %>%
    filter(gene_name %in% genes) %>%
    ggplot(aes(y = expr, color = genotype)) +
    geom_boxplot(aes_string(x = x_var), width = 0.7, outlier.alpha = 0, size = 0.4) +
    geom_point(aes_string(x = x_var, shape = shape_var, group = "genotype"), size = 1,
               position = position_dodge(0.7)) +
    facet_wrap(~gene_name, scales = "free_y", ncol = ncol) +
    scale_color_manual(values = pal) +
    scale_x_discrete(labels = x_breaks) +
    labs(x = NULL, y = ifelse(log, "Log2(TPM + 1)", "TPM")) +
    theme_cowplot(font_size = 10)
  
  ggsave(paste0("figures/boxplot_", suffix, ifelse(log, "_log2", ""), ".png"), 
         dpi = 400, width = width, height = height)
  ggsave(paste0("figures/boxplot_", suffix, ifelse(log, "_log2", ""), ".pdf"), 
         width = width, height = height)
}

## palette
pal <- c("dodgerblue3", "#EF3B2C", friendly_pal("nickel_five", 5)[3],
         friendly_pal("ito_seven", 7)[c(6, 4:3)])

## Havcr2 expression in Dataset 3, 4 month
boxplot_genes("Havcr2", tpm = tpm_ds3_long, ncol = 1, suffix = "ds3_4M_Havcr2", 
              height = 2.5, width = 4, pal = pal, log = FALSE)
## Havcr2 expression in Dataset 1 Batch 2, 7 month
boxplot_genes("Havcr2", tpm = tpm_ds1_batch2_long, ncol = 1, suffix = "ds1_batch2_7M_Havcr2", 
              height = 2.5, width = 4, pal = pal, log = FALSE)

