
#!/usr/bin/Rscript
library(tidyverse)
library(dplyr)
library(SummarizedExperiment)
library(rlang)
library(limma)
library(edgeR)
library(sessioninfo)

#-------------------------------------------------------------------------------
#        7.2 DGE analysis on permuted cell line proliferation phenotypes
#-------------------------------------------------------------------------------
#  Code to run within-treatment DGE on permuted cell line proliferation 
#  fractions to assess the significance of the associations between 
#  proliferation from preMac -> microglia and gene expr.
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

## Set working dir
setwd("/lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_Daianna/")

## Define output data and plot dir 
outdir = paste(getwd(), "output_data", "07_Diff_expr_proliferation_permutations", "02_permuted_proliferation_DGE", sep = "/")
plotdir = paste(getwd(), "plots", "07_Diff_expr_proliferation_permutations", "02_permuted_proliferation_DGE", sep = "/")
dir.create(outdir, recursive = T)
dir.create(plotdir, recursive = T)

## Input dir
input_dir01 = paste(getwd(), "output_data", "07_Diff_expr_proliferation_permutations", "01_run_proliferation_DGE", sep = "/")
## Job array index (i)
args = commandArgs(trailingOnly = TRUE)
index = args[1]

## load rse with migration data
load(paste0(input_dir01, "/rse_non_prolif_filt_expr_proliferation_IFN.Rdata"), verbose = T)
load(paste0(input_dir01, "/rse_non_prolif_filt_expr_proliferation_LPS.Rdata"), verbose = T)
load(paste0(input_dir01, "/rse_non_prolif_filt_expr_proliferation_untreated.Rdata"), verbose = T)


dge_permutation_results <- list()

## Run DGE within each treatment for each permutation i
for(treatment in c("IFN", "LPS", "untreated")){
  
  rse <- eval(parse_expr(paste("rse_non_prolif_filt_expr_proliferation", treatment, sep = "_")))
  
  ## Add TMM norm factors for voom normalization
  rse_norm <- calcNormFactors(rse, method = "TMM")
  
  ## DGE for permuted scaled log-fractions 
  coef <- paste0("permutation_proliferation_mean_scaled_log_fraction_", index)
  formula <- as.formula(paste("~", coef, "+ Sex + proliferation_mean_min_prop + genotype_PC1 + genotype_PC2"))
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

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.1 (2023-06-16)
# os       Ubuntu 22.04.4 LTS
# system   x86_64, linux-gnu
# ui       RStudio
# language (EN)
# collate  en_GB.UTF-8
# ctype    en_GB.UTF-8
# tz       Europe/Belfast
# date     2025-01-17
# rstudio  2024.04.0+735 Chocolate Cosmos (server)
# pandoc   3.1.12.3 @ /opt/view/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version    date (UTC) lib source
# abind                  1.4-5      2016-07-21 [1] CRAN (R 4.3.1)
# aod                    1.3.3      2023-12-13 [1] CRAN (R 4.3.1)
# backports              1.4.1      2021-12-13 [1] CRAN (R 4.3.1)
# Biobase              * 2.62.0     2023-10-24 [1] Bioconductor
# BiocGenerics         * 0.48.1     2023-11-01 [1] Bioconductor
# BiocParallel         * 1.36.0     2023-10-24 [1] Bioconductor
# bitops                 1.0-7      2021-04-24 [1] CRAN (R 4.3.1)
# boot                   1.3-30     2024-02-26 [1] CRAN (R 4.3.1)
# broom                  1.0.5      2023-06-09 [1] CRAN (R 4.3.1)
# Cairo                  1.6-2      2023-11-28 [1] CRAN (R 4.3.1)
# caTools                1.18.2     2021-03-28 [1] CRAN (R 4.3.1)
# circlize               0.4.16     2024-02-20 [1] CRAN (R 4.3.1)
# cli                    3.6.2      2023-12-11 [1] CRAN (R 4.3.1)
# clue                   0.3-65     2023-09-23 [1] CRAN (R 4.3.1)
# cluster                2.1.6      2023-12-01 [1] CRAN (R 4.3.1)
# codetools              0.2-20     2024-03-31 [1] CRAN (R 4.3.1)
# colorspace             2.1-0      2023-01-23 [1] CRAN (R 4.3.1)
# ComplexHeatmap       * 2.18.0     2023-10-24 [1] Bioconductor
# corpcor                1.6.10     2021-09-16 [1] CRAN (R 4.3.1)
# cowplot              * 1.1.3      2024-01-22 [1] CRAN (R 4.3.1)
# crayon                 1.5.2      2022-09-29 [1] CRAN (R 4.3.1)
# DelayedArray           0.28.0     2023-10-24 [1] Bioconductor
# digest                 0.6.35     2024-03-11 [1] CRAN (R 4.3.1)
# doParallel             1.0.17     2022-02-07 [1] CRAN (R 4.3.1)
# dplyr                * 1.1.4      2023-11-17 [1] CRAN (R 4.3.1)
# edgeR                * 4.0.16     2024-02-18 [1] Bioconductor 3.18 (R 4.3.1)
# EnvStats               2.8.1      2023-08-22 [1] CRAN (R 4.3.1)
# fANCOVA                0.6-1      2020-11-13 [1] CRAN (R 4.3.1)
# fansi                  1.0.6      2023-12-08 [1] CRAN (R 4.3.1)
# farver                 2.1.1      2022-07-06 [1] CRAN (R 4.3.1)
# forcats              * 1.0.0      2023-01-29 [1] CRAN (R 4.3.1)
# foreach                1.5.2      2022-02-02 [1] CRAN (R 4.3.1)
# generics               0.1.3      2022-07-05 [1] CRAN (R 4.3.1)
# GenomeInfoDb         * 1.38.8     2024-03-15 [1] Bioconductor 3.18 (R 4.3.1)
# GenomeInfoDbData       1.2.11     2024-11-06 [1] Bioconductor
# GenomicRanges        * 1.54.1     2023-10-29 [1] Bioconductor
# GetoptLong             1.0.5      2020-12-15 [1] CRAN (R 4.3.1)
# ggplot2              * 3.5.1      2024-04-23 [1] CRAN (R 4.3.1)
# ggrepel                0.9.5      2024-01-10 [1] CRAN (R 4.3.1)
# GlobalOptions          0.1.2      2020-06-10 [1] CRAN (R 4.3.1)
# glue                   1.7.0      2024-01-09 [1] CRAN (R 4.3.1)
# gplots                 3.1.3.1    2024-02-02 [1] CRAN (R 4.3.1)
# gtable                 0.3.4      2023-08-21 [1] CRAN (R 4.3.1)
# gtools                 3.9.5      2023-11-20 [1] CRAN (R 4.3.1)
# hms                    1.1.3      2023-03-21 [1] CRAN (R 4.3.1)
# IRanges              * 2.36.0     2023-10-24 [1] Bioconductor
# iterators              1.0.14     2022-02-05 [1] CRAN (R 4.3.1)
# KernSmooth             2.23-22    2023-07-10 [1] CRAN (R 4.3.1)
# labeling               0.4.3      2023-08-29 [1] CRAN (R 4.3.1)
# lattice                0.22-6     2024-03-20 [1] CRAN (R 4.3.1)
# lifecycle              1.0.4      2023-11-07 [1] CRAN (R 4.3.1)
# limma                * 3.58.1     2023-10-31 [1] Bioconductor
# lme4                   1.1-35.2   2024-03-28 [1] CRAN (R 4.3.1)
# lmerTest               3.1-3      2020-10-23 [1] CRAN (R 4.3.1)
# locfit                 1.5-9.9    2024-03-01 [1] CRAN (R 4.3.1)
# lubridate            * 1.9.3      2023-09-27 [1] CRAN (R 4.3.1)
# magick                 2.8.3      2024-02-18 [1] CRAN (R 4.3.1)
# magrittr               2.0.3      2022-03-30 [1] CRAN (R 4.3.1)
# MASS                   7.3-60.0.1 2024-01-13 [1] CRAN (R 4.3.1)
# Matrix                 1.6-5      2024-01-11 [1] CRAN (R 4.3.1)
# MatrixGenerics       * 1.14.0     2023-10-24 [1] Bioconductor
# matrixStats          * 1.2.0      2023-12-11 [1] CRAN (R 4.3.1)
# minqa                  1.2.6      2023-09-11 [1] CRAN (R 4.3.1)
# munsell                0.5.1      2024-04-01 [1] CRAN (R 4.3.1)
# mvtnorm                1.2-4      2023-11-27 [1] CRAN (R 4.3.1)
# nlme                   3.1-164    2023-11-27 [1] CRAN (R 4.3.1)
# nloptr                 2.0.3      2022-05-26 [1] CRAN (R 4.3.1)
# numDeriv               2016.8-1.1 2019-06-06 [1] CRAN (R 4.3.1)
# pbkrtest               0.5.2      2023-01-19 [1] CRAN (R 4.3.1)
# pillar                 1.9.0      2023-03-22 [1] CRAN (R 4.3.1)
# pkgconfig              2.0.3      2019-09-22 [1] CRAN (R 4.3.1)
# plyr                   1.8.9      2023-10-02 [1] CRAN (R 4.3.1)
# png                    0.1-8      2022-11-29 [1] CRAN (R 4.3.1)
# purrr                * 1.0.2      2023-08-10 [1] CRAN (R 4.3.1)
# R6                     2.5.1      2021-08-19 [1] CRAN (R 4.3.1)
# ragg                   1.3.0      2024-03-13 [1] CRAN (R 4.3.1)
# rbibutils              2.2.16     2023-10-25 [1] CRAN (R 4.3.1)
# RColorBrewer           1.1-3      2022-04-03 [1] CRAN (R 4.3.1)
# Rcpp                   1.0.12     2024-01-09 [1] CRAN (R 4.3.1)
# RCurl                  1.98-1.14  2024-01-09 [1] CRAN (R 4.3.1)
# Rdpack                 2.6        2023-11-08 [1] CRAN (R 4.3.1)
# readr                * 2.1.5      2024-01-10 [1] CRAN (R 4.3.1)
# remaCor                0.0.18     2024-02-08 [1] CRAN (R 4.3.1)
# reshape2               1.4.4      2020-04-09 [1] CRAN (R 4.3.1)
# RhpcBLASctl            0.23-42    2023-02-11 [1] CRAN (R 4.3.1)
# rjson                  0.2.21     2022-01-09 [1] CRAN (R 4.3.1)
# rlang                * 1.1.3      2024-01-10 [1] CRAN (R 4.3.1)
# rstudioapi             0.16.0     2024-03-24 [1] CRAN (R 4.3.1)
# S4Arrays               1.2.1      2024-03-04 [1] Bioconductor 3.18 (R 4.3.1)
# S4Vectors            * 0.40.2     2023-11-23 [1] Bioconductor 3.18 (R 4.3.1)
# scales                 1.3.0      2023-11-28 [1] CRAN (R 4.3.1)
# sessioninfo          * 1.2.2      2021-12-06 [1] CRAN (R 4.3.1)
# shape                  1.4.6.1    2024-02-23 [1] CRAN (R 4.3.1)
# SparseArray            1.2.4      2024-02-11 [1] Bioconductor 3.18 (R 4.3.1)
# statmod                1.5.0      2023-01-06 [1] CRAN (R 4.3.1)
# stringi                1.8.3      2023-12-11 [1] CRAN (R 4.3.1)
# stringr              * 1.5.1      2023-11-14 [1] CRAN (R 4.3.1)
# SummarizedExperiment * 1.32.0     2023-10-24 [1] Bioconductor
# systemfonts            1.0.6      2024-03-07 [1] CRAN (R 4.3.1)
# textshaping            0.3.7      2023-10-09 [1] CRAN (R 4.3.1)
# tibble               * 3.2.1      2023-03-20 [1] CRAN (R 4.3.1)
# tidyr                * 1.3.1      2024-01-24 [1] CRAN (R 4.3.1)
# tidyselect             1.2.1      2024-03-11 [1] CRAN (R 4.3.1)
# tidyverse            * 2.0.0      2023-02-22 [1] CRAN (R 4.3.1)
# timechange             0.3.0      2024-01-18 [1] CRAN (R 4.3.1)
# tzdb                   0.4.0      2023-05-12 [1] CRAN (R 4.3.1)
# utf8                   1.2.4      2023-10-22 [1] CRAN (R 4.3.1)
# variancePartition    * 1.32.5     2024-02-16 [1] Bioconductor 3.18 (R 4.3.1)
# vctrs                  0.6.5      2023-12-01 [1] CRAN (R 4.3.1)
# withr                  3.0.0      2024-01-16 [1] CRAN (R 4.3.1)
# XVector                0.42.0     2023-10-24 [1] Bioconductor
# zlibbioc               1.48.2     2024-03-13 [1] Bioconductor 3.18 (R 4.3.1)
# 
# [1] /opt/view/rlib/R/library
# [2] /software/hgi/installs/softpack/rstudio/.spack/linux-ubuntu22.04-x86_64_v3/gcc-11.4.0/r-4.3.1-bfwldrk76z6f52upk47zepliekn7ayqz/rlib/R/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
