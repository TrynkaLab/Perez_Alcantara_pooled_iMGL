
#!/usr/bin/Rscript
library(tidyverse)
library(dplyr)
library(SummarizedExperiment)
library(rlang)
library(limma)
library(edgeR)
library(sessioninfo)

#-------------------------------------------------------------------------------
#        5.2 DGE analysis on permuted cell line phagocytosis phenotypes
#-------------------------------------------------------------------------------
#  Code to run within-treatment DGE on permuted cell line phagocytosis 
#  estimates to assess the significance of phagocytosis associations with
#  gene expr.
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

## Set working dir
setwd("/lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_Daianna/")

## Define output data and plot dir 
outdir = paste(getwd(), "output_data", "05_Diff_expr_phagocytosis_permutations", "02_permuted_phagocytosis_DGE", sep = "/")
plotdir = paste(getwd(), "plots", "05_Diff_expr_phagocytosis_permutations", "02_permuted_phagocytosis_DGE", sep = "/")
dir.create(outdir, recursive = T)
dir.create(plotdir, recursive = T)

## Input dir
input_dir01 = paste(getwd(), "output_data", "05_Diff_expr_phagocytosis_permutations", "01_run_phagocytosis_DGE", sep = "/")
## Job array index (i)
args = commandArgs(trailingOnly=TRUE)
index = args[1]

## load rse with phagocytosis data
load(paste0(input_dir01, "/rse_non_prolif_filt_expr_phago_IFN.Rdata"), verbose = T)
load(paste0(input_dir01, "/rse_non_prolif_filt_expr_phago_LPS.Rdata"), verbose = T)
load(paste0(input_dir01, "/rse_non_prolif_filt_expr_phago_untreated.Rdata"), verbose = T)


dge_permutation_results <- list()

## Run DGE within each treatment for each permutation i
for(treatment in c("IFN", "LPS", "untreated")){
  
  rse <- eval(parse_expr(paste("rse_non_prolif_filt_expr_phago", treatment, sep = "_")))
  
  ## Add TMM norm factors for voom normalization
  rse_norm <- calcNormFactors(rse, method = "TMM")

  ## DGE for permuted scaled log-fractions 
  coef <- paste0("permutation_phago_mean_scaled_log_fraction_", index)
  formula <- as.formula(paste("~", coef, "+ Sex + phagocytosis_mean_min_prop + genotype_PC1 + genotype_PC2"))
  model = model.matrix(formula, data = rse_norm$samples) %>%  as.data.frame()
  
  v = voom(rse_norm, design = model, plot = TRUE)
  
  ## Estimate intra-pool corr in gene expr
  cor = duplicateCorrelation(v, design = model, block = rse_norm$samples$pool)
  
  ## Re-compute voom weights based on within-pool corr
  v2 = voom(rse_norm, design = model, plot = TRUE, block = rse_norm$samples$pool, correlation = cor$consensus)
  ## Estimate corr based on corrected weights
  cor2 = duplicateCorrelation(v2, design = model, block = rse_norm$samples$pool)

  ## Fit linear model
  fit = lmFit(v2, design = model, block = rse_norm$samples$pool, correlation = cor2$consensus)
  eBGene = eBayes(fit)

  ## Extract results for coeff of interest
  top_genes = topTable(eBGene, coef = coef, p.value = 1, number = nrow(rse), sort.by = "none")
  
  ## DEGs?
  de_genes_per <- subset(top_genes, adj.P.Val<0.05)
  message("DGE on permuted results:")
  print(paste(dim(de_genes_per)[1], "DEGs:", 
              dim(subset(de_genes_per, logFC>0))[1], "up-regulated and", 
              dim(subset(de_genes_per, logFC<0))[1], "down-regulated"))
  
  dge_permutation_results[[treatment]] <- top_genes
}

save(dge_permutation_results, file = paste0(outdir, "/dge_results_permutation_", index, ".Rdata"))







## Reproducibility info
options(width = 120)
session_info()

# ─ Session info ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.1 (2023-06-16)
# os       Ubuntu 22.04.4 LTS
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_GB.utf8
# ctype    en_GB.utf8
# tz       UTC
# date     2024-12-27
# pandoc   3.1.12.3 @ /opt/view/bin/pandoc
# 
# ─ Packages ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [1] CRAN (R 4.3.1)
# Biobase              * 2.62.0    2023-10-24 [1] Bioconductor
# BiocGenerics         * 0.48.1    2023-11-01 [1] Bioconductor
# bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.3.1)
# cli                    3.6.2     2023-12-11 [1] CRAN (R 4.3.1)
# colorspace             2.1-0     2023-01-23 [1] CRAN (R 4.3.1)
# crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.1)
# DelayedArray           0.28.0    2023-10-24 [1] Bioconductor
# dplyr                * 1.1.4     2023-11-17 [1] CRAN (R 4.3.1)
# edgeR                * 4.0.16    2024-02-18 [1] Bioconductor 3.18 (R 4.3.1)
# fansi                  1.0.6     2023-12-08 [1] CRAN (R 4.3.1)
# forcats              * 1.0.0     2023-01-29 [1] CRAN (R 4.3.1)
# generics               0.1.3     2022-07-05 [1] CRAN (R 4.3.1)
# GenomeInfoDb         * 1.38.8    2024-03-15 [1] Bioconductor 3.18 (R 4.3.1)
# GenomeInfoDbData       1.2.11    2024-11-06 [1] Bioconductor
# GenomicRanges        * 1.54.1    2023-10-29 [1] Bioconductor
# ggplot2              * 3.5.1     2024-04-23 [1] CRAN (R 4.3.1)
# glue                   1.7.0     2024-01-09 [1] CRAN (R 4.3.1)
# gtable                 0.3.4     2023-08-21 [1] CRAN (R 4.3.1)
# hms                    1.1.3     2023-03-21 [1] CRAN (R 4.3.1)
# IRanges              * 2.36.0    2023-10-24 [1] Bioconductor
# lattice                0.22-6    2024-03-20 [1] CRAN (R 4.3.1)
# lifecycle              1.0.4     2023-11-07 [1] CRAN (R 4.3.1)
# limma                * 3.58.1    2023-10-31 [1] Bioconductor
# locfit                 1.5-9.9   2024-03-01 [1] CRAN (R 4.3.1)
# lubridate            * 1.9.3     2023-09-27 [1] CRAN (R 4.3.1)
# magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.3.1)
# Matrix                 1.6-5     2024-01-11 [1] CRAN (R 4.3.1)
# MatrixGenerics       * 1.14.0    2023-10-24 [1] Bioconductor
# matrixStats          * 1.2.0     2023-12-11 [1] CRAN (R 4.3.1)
# munsell                0.5.1     2024-04-01 [1] CRAN (R 4.3.1)
# pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.3.1)
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.1)
# purrr                * 1.0.2     2023-08-10 [1] CRAN (R 4.3.1)
# R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.1)
# Rcpp                   1.0.12    2024-01-09 [1] CRAN (R 4.3.1)
# RCurl                  1.98-1.14 2024-01-09 [1] CRAN (R 4.3.1)
# readr                * 2.1.5     2024-01-10 [1] CRAN (R 4.3.1)
# rlang                * 1.1.3     2024-01-10 [1] CRAN (R 4.3.1)
# S4Arrays               1.2.1     2024-03-04 [1] Bioconductor 3.18 (R 4.3.1)
# S4Vectors            * 0.40.2    2023-11-23 [1] Bioconductor 3.18 (R 4.3.1)
# scales                 1.3.0     2023-11-28 [1] CRAN (R 4.3.1)
# sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.3.1)
# SparseArray            1.2.4     2024-02-11 [1] Bioconductor 3.18 (R 4.3.1)
# statmod                1.5.0     2023-01-06 [1] CRAN (R 4.3.1)
# stringi                1.8.3     2023-12-11 [1] CRAN (R 4.3.1)
# stringr              * 1.5.1     2023-11-14 [1] CRAN (R 4.3.1)
# SummarizedExperiment * 1.32.0    2023-10-24 [1] Bioconductor
# tibble               * 3.2.1     2023-03-20 [1] CRAN (R 4.3.1)
# tidyr                * 1.3.1     2024-01-24 [1] CRAN (R 4.3.1)
# tidyselect             1.2.1     2024-03-11 [1] CRAN (R 4.3.1)
# tidyverse            * 2.0.0     2023-02-22 [1] CRAN (R 4.3.1)
# timechange             0.3.0     2024-01-18 [1] CRAN (R 4.3.1)
# tzdb                   0.4.0     2023-05-12 [1] CRAN (R 4.3.1)
# utf8                   1.2.4     2023-10-22 [1] CRAN (R 4.3.1)
# vctrs                  0.6.5     2023-12-01 [1] CRAN (R 4.3.1)
# withr                  3.0.0     2024-01-16 [1] CRAN (R 4.3.1)
# XVector                0.42.0    2023-10-24 [1] Bioconductor
# zlibbioc               1.48.2    2024-03-13 [1] Bioconductor 3.18 (R 4.3.1)

