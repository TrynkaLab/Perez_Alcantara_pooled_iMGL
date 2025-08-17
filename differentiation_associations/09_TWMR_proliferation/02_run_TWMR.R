
library(tidyverse)
library(dplyr)
library(sessioninfo)

#-------------------------------------------------------------------------------
#             9.2 Run multivariable MR analysis per focal gene*
#-------------------------------------------------------------------------------
#  Code to find the multivariate causal effect of each focal gene on the 
#  proliferation from preMac -> microglia in IFN/LPS/untreated based on the   
#  inverse-variance weighted method for eQTL and GWAS summary statistics.
#  * Note: code based on Marta's code, which in turn is based on the authors 
#          implementation in https://github.com/eleporcu/TWMR/blob/master/MR.R.
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

## Define dirs
input_dir = paste("..", "..", "input_data", "09_TWMR_proliferation", "02_run_TWMR", sep = "/")
plotdir = paste("..", "..", "plots", "09_TWMR_proliferation", "02_run_TWMR", sep = "/")

dir.create(input_dir, recursive = T, showWarnings = F)
dir.create(plotdir, recursive = T, showWarnings = F)

## Test
# focal_gene <-  "FLOT1"
# treatment <- "IFN"

for (treatment in c("IFN", "LPS", "untreated")){
  
  ## Input dirs
  input_dir01out <- paste("..", "..", "output_data", "09_TWMR_proliferation", "01_IVs_and_exposures_selection", treatment, sep = "/") 
  input_dir01in <-paste("..", "..", "input_data", "09_TWMR_proliferation", "01_IVs_and_exposures_selection", sep = "/")
  
  ## Output dir
  outdir = paste("..", "..", "output_data", "09_TWMR_proliferation", "02_run_TWMR", treatment, sep = "/")
  dir.create(outdir, recursive = T, showWarnings = F)
  
  ## Confirm all genes tested in eQTL file have a TWMR model
  nominal_eqtl_results = readr::read_delim(paste0(input_dir01in, "/nominal_eqtl_results_", treatment, ".csv")) %>% 
    filter(! is.na(gene_symbol))

  length(list.files(input_dir01out)) == length(unique(nominal_eqtl_results$gene_symbol))
  # TRUE
  
  ## All genes 
  all_genes <- str_split_fixed(list.files(input_dir01out), "_", 2)[,1]
  
  for(focal_gene in all_genes){
    
    ## Skip genes already read 
    if(focal_gene %in% str_split_fixed(list.files(outdir), "_", 2)[,1]){
      message("Gene already read.")
      next()
    }
    
    ## Model for focal gene
    load(paste0(input_dir01out, "/", focal_gene, "_input_matrices.Rdata"), verbose = T)
    
    # * E (nxk): matrix of effects of the n independent IVs on k uncorr gene expression exposures
    # * G (1xn): vector of effects of the n independent IVs on proliferation 
    # * C (nxn): pair-wise LD matrix between the n independent IVs 
    # * Ngwas: number of participants in GWAS
    # * N_eQTLs: number of participants in eQTL
    # * n_for_LD: N used for LD calculation (from 1000 Genomes)
 
    E = inputs$E
    G = inputs$G
    C = inputs$C
    N_eQTLs = inputs$N_eQTLs
    Ngwas = inputs$Ngwas
    n_for_LD = inputs$n_for_LD
    
    
    if(is.null(names(inputs))){
      E = inputs[[1]]
      G = inputs[[2]]
      C = inputs[[3]]
      N_eQTLs = inputs[[4]]
      Ngwas = inputs[[5]]
      n_for_LD = inputs[[6]] 
    }
    
    
    ## Sanity checks!!!! :)
    ## Check same num of IVs in inputs
    if(! (nrow(E) == length(G) & nrow(E) == nrow(C) & nrow(C) == ncol(C))){
      message(paste0("Incorrect dimensions of input matrices for ", focal_gene, "."))
      stop()
    }
    ## Check same IVs included
    if(! identical(rownames(E), names(G))){
      message(paste0("Not same IVs included in E and G for ", focal_gene, "."))
      stop()
    }
    if(! (identical(rownames(E), rownames(C))) | ! (identical(rownames(E), colnames(C))) ){
      message(paste0("Not same IVs included in E and C for ", focal_gene, "."))
      stop()
    }
    ## Check focal gene is col 1 in E
    if(colnames(E)[1] != focal_gene){
      message(paste0("Put ", focal_gene, " in 1st column of E."))
      stop()
    }
    

    ## Genes with all zero effects?
    x = colSums(abs(E))
    zeroGenes = which(x == 0)
    if (length(zeroGenes) > 0) {
      E = E[, -zeroGenes]
    }
    E = as.matrix(E)
    
    
    ## Calculate alpha (equation 2): 
    
    ## Inverse of C: SNP-SNP variance-covariance matrix
    C_inv = solve(C)
    ## E'C^-1E: indirect/direct SNP effects on each gene pair
    EtCinvE = t(E) %*% C_inv %*% E
    H = (1 - 1 / sqrt(n_for_LD)) * EtCinvE + (1 / sqrt(n_for_LD)) * diag(nrow(EtCinvE)) # ?
    ## E'C^-1G: aggregated SNP effects on gene and prolif
    EtCinvG = t(E) %*% C_inv %*% G
    ## Gene causal effects
    alpha = solve(H) %*% EtCinvG
    alpha = as.vector(alpha)
    
    
    # Calculate standard error (equation 3)
    GCG_inv = solve(H)
    
    df_dg = GCG_inv %*% t(E) %*% C_inv
    df_dG = (GCG_inv %x% (t(G) %*% C_inv %*% ((E %*% GCG_inv %*% t(E)) %*% C_inv + diag(nrow(E))
    ))) +
      ((-t(G) %*% C_inv %*% E %*% GCG_inv) %x% (GCG_inv %*% t(E) %*% C_inv))
    J = cbind(df_dG, df_dg)
    
    SEs = c(rep(1 / sqrt(N_eQTLs), length(E[1, ]) * length(E[, 1])), rep(1 / sqrt(Ngwas), length(G)))
    R = diag(ncol(E) + 1)
    Sigma = (SEs %*% t(SEs)) * (C %x% R)
    V = J %*% Sigma %*% t(J)
    se = sqrt(V[1, 1])
    
    
    ## Calculate focal gene causal effect Z-statistic and p-value
    Z = alpha[1] / se
    pval = 2 * pnorm(abs(Z), lower.tail = FALSE)
    
    ## Save outputs
    Nsnps = nrow(E)
    Ngene = ncol(E)
    
    gene_results <- list(
      "alpha" = alpha[1],
      "Z" = Z,
      "SE" = se,
      "P" = pval,
      "Nsnps" = Nsnps,
      "Ngene" = Ngene
    )
    
    save(gene_results, file = paste0(outdir, "/", focal_gene, "_results.Rdata"))
  }
}


## Check all genes have results
length(list.files(outdir)) == length(list.files(input_dir01out))
# [1] TRUE







# Reproducibility info
# session_info()

# ─ Session info ──────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.1 (2023-06-16)
# os       Ubuntu 22.04.4 LTS
# system   x86_64, linux-gnu
# ui       RStudio
# language (EN)
# collate  en_GB.UTF-8
# ctype    en_GB.UTF-8
# tz       Europe/Belfast
# date     2025-03-19
# rstudio  2024.04.0+735 Chocolate Cosmos (server)
# pandoc   3.1.12.3 @ /opt/view/bin/ (via rmarkdown)
# 
# ─ Packages ──────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package          * version     date (UTC) lib source
# AnnotationDbi      1.64.1      2023-11-03 [1] Bioconductor
# Biobase            2.62.0      2023-10-24 [1] Bioconductor
# BiocFileCache      2.10.2      2024-03-27 [1] Bioconductor 3.18 (R 4.3.1)
# BiocGenerics       0.48.1      2023-11-01 [1] Bioconductor
# biomaRt          * 2.58.2      2024-01-30 [1] Bioconductor 3.18 (R 4.3.1)
# Biostrings         2.70.3      2024-03-13 [1] Bioconductor 3.18 (R 4.3.1)
# bit                4.0.5       2022-11-15 [1] CRAN (R 4.3.1)
# bit64              4.0.5       2020-08-30 [1] CRAN (R 4.3.1)
# bitops             1.0-7       2021-04-24 [1] CRAN (R 4.3.1)
# blob               1.2.4       2023-03-17 [1] CRAN (R 4.3.1)
# cachem             1.0.8       2023-05-01 [1] CRAN (R 4.3.1)
# cli                3.6.2       2023-12-11 [1] CRAN (R 4.3.1)
# colorspace         2.1-0       2023-01-23 [1] CRAN (R 4.3.1)
# cowplot          * 1.1.3       2024-01-22 [1] CRAN (R 4.3.1)
# crayon             1.5.2       2022-09-29 [1] CRAN (R 4.3.1)
# curl               5.2.1       2024-03-01 [1] CRAN (R 4.3.1)
# DBI                1.2.2       2024-02-16 [1] CRAN (R 4.3.1)
# dbplyr             2.5.0       2024-03-19 [1] CRAN (R 4.3.1)
# digest             0.6.35      2024-03-11 [1] CRAN (R 4.3.1)
# dplyr            * 1.1.4       2023-11-17 [1] CRAN (R 4.3.1)
# evaluate           0.23        2023-11-01 [1] CRAN (R 4.3.1)
# fansi              1.0.6       2023-12-08 [1] CRAN (R 4.3.1)
# fastmap            1.1.1       2023-02-24 [1] CRAN (R 4.3.1)
# filelock           1.0.3       2023-12-11 [1] CRAN (R 4.3.1)
# forcats          * 1.0.0       2023-01-29 [1] CRAN (R 4.3.1)
# generics           0.1.3       2022-07-05 [1] CRAN (R 4.3.1)
# GenomeInfoDb       1.38.8      2024-03-15 [1] Bioconductor 3.18 (R 4.3.1)
# GenomeInfoDbData   1.2.11      2025-02-07 [1] Bioconductor
# ggplot2          * 3.5.1       2024-04-23 [1] CRAN (R 4.3.1)
# glue               1.7.0       2024-01-09 [1] CRAN (R 4.3.1)
# gtable             0.3.4       2023-08-21 [1] CRAN (R 4.3.1)
# here               1.0.1       2020-12-13 [1] CRAN (R 4.3.1)
# hms                1.1.3       2023-03-21 [1] CRAN (R 4.3.1)
# htmltools          0.5.8       2024-03-25 [1] CRAN (R 4.3.1)
# httr             * 1.4.7       2023-08-15 [1] CRAN (R 4.3.1)
# IRanges            2.36.0      2023-10-24 [1] Bioconductor
# jsonlite           1.8.8       2023-12-04 [1] CRAN (R 4.3.1)
# KEGGREST           1.42.0      2023-10-24 [1] Bioconductor
# knitr              1.45        2023-10-30 [1] CRAN (R 4.3.1)
# lattice            0.22-6      2024-03-20 [1] CRAN (R 4.3.1)
# lifecycle          1.0.4       2023-11-07 [1] CRAN (R 4.3.1)
# lubridate        * 1.9.3       2023-09-27 [1] CRAN (R 4.3.1)
# magrittr           2.0.3       2022-03-30 [1] CRAN (R 4.3.1)
# Matrix             1.6-5       2024-01-11 [1] CRAN (R 4.3.1)
# memoise            2.0.1       2021-11-26 [1] CRAN (R 4.3.1)
# munsell            0.5.1       2024-04-01 [1] CRAN (R 4.3.1)
# pillar             1.9.0       2023-03-22 [1] CRAN (R 4.3.1)
# pkgconfig          2.0.3       2019-09-22 [1] CRAN (R 4.3.1)
# png                0.1-8       2022-11-29 [1] CRAN (R 4.3.1)
# prettyunits        1.2.0       2023-09-24 [1] CRAN (R 4.3.1)
# progress           1.2.3       2023-12-06 [1] CRAN (R 4.3.1)
# purrr            * 1.0.2       2023-08-10 [1] CRAN (R 4.3.1)
# R6                 2.5.1       2021-08-19 [1] CRAN (R 4.3.1)
# rappdirs           0.3.3       2021-01-31 [1] CRAN (R 4.3.1)
# Rcpp               1.0.12      2024-01-09 [1] CRAN (R 4.3.1)
# RCurl              1.98-1.14   2024-01-09 [1] CRAN (R 4.3.1)
# readr            * 2.1.5       2024-01-10 [1] CRAN (R 4.3.1)
# reticulate         1.35.0      2024-01-31 [1] CRAN (R 4.3.1)
# rlang              1.1.3       2024-01-10 [1] CRAN (R 4.3.1)
# rmarkdown          2.26        2024-03-05 [1] CRAN (R 4.3.1)
# rprojroot          2.0.4       2023-11-05 [1] CRAN (R 4.3.1)
# RSQLite            2.3.6       2024-03-31 [1] CRAN (R 4.3.1)
# rstudioapi         0.16.0      2024-03-24 [1] CRAN (R 4.3.1)
# S4Vectors          0.40.2      2023-11-23 [1] Bioconductor 3.18 (R 4.3.1)
# scales             1.3.0       2023-11-28 [1] CRAN (R 4.3.1)
# sessioninfo      * 1.2.2       2021-12-06 [1] CRAN (R 4.3.1)
# stringi            1.8.3       2023-12-11 [1] CRAN (R 4.3.1)
# stringr          * 1.5.1       2023-11-14 [1] CRAN (R 4.3.1)
# tibble           * 3.2.1       2023-03-20 [1] CRAN (R 4.3.1)
# tidyr            * 1.3.1       2024-01-24 [1] CRAN (R 4.3.1)
# tidyselect         1.2.1       2024-03-11 [1] CRAN (R 4.3.1)
# tidyverse        * 2.0.0       2023-02-22 [1] CRAN (R 4.3.1)
# timechange         0.3.0       2024-01-18 [1] CRAN (R 4.3.1)
# tzdb               0.4.0       2023-05-12 [1] CRAN (R 4.3.1)
# utf8               1.2.4       2023-10-22 [1] CRAN (R 4.3.1)
# vctrs              0.6.5       2023-12-01 [1] CRAN (R 4.3.1)
# withr              3.0.0       2024-01-16 [1] CRAN (R 4.3.1)
# xfun               0.43        2024-03-25 [1] CRAN (R 4.3.1)
# XML                3.99-0.16.1 2024-01-22 [1] CRAN (R 4.3.1)
# xml2               1.3.6       2023-12-04 [1] CRAN (R 4.3.1)
# XVector            0.42.0      2023-10-24 [1] Bioconductor
# yaml               2.3.8       2023-12-11 [1] CRAN (R 4.3.1)
# zlibbioc           1.48.2      2024-03-13 [1] Bioconductor 3.18 (R 4.3.1)
# 
# [1] /opt/view/rlib/R/library
# [2] /software/hgi/installs/softpack/rstudio/.spack/linux-ubuntu22.04-x86_64_v3/gcc-11.4.0/r-4.3.1-bfwldrk76z6f52upk47zepliekn7ayqz/rlib/R/library
# 
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
