
library(tidyverse)
library(dplyr)
library(SummarizedExperiment)
library(rlang)
library(edgeR)
library(sessioninfo)


################################################################################
##                6. Permutation analysis for migration DGE
################################################################################

#-------------------------------------------------------------------------------
#                     6.0 Data exploration and processing
#-------------------------------------------------------------------------------
#  Code to explore, process, and add cell line migration data to sample metadata
#  to run DGE on permuted phenotypes in the next script.
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

## Set working dir
setwd("/lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_Daianna/")

## Define output data and plot dir
outdir = paste(getwd(), "output_data", "06_Diff_expr_migration_permutations", "00_Data_exploration_processing", sep = "/")
plotdir = paste(getwd(), "plots", "06_Diff_expr_migration_permutations", "00_Data_exploration_processing", sep = "/")
dir.create(outdir, recursive = T)
dir.create(plotdir, recursive = T)

## Input dirs
input_dir05 = paste(getwd(), "output_data", "05_Diff_expr_phagocytosis_permutations", "00_Data_exploration_processing", sep = "/")

## Load complete rse (with genotype PCs)
load(paste0(input_dir05, "/rse_complete_with_phagocitosis_data.Rdata"), verbose = T)


##############################################
##  Add migration data to sample metadata
##############################################

line_prop_changes = read_csv("../OTAR2065_phenotypic_QTLs/data/results/migration/1.check_line_proportions/line_prop_changes_per_well.csv") %>% as.data.frame()
dim(line_prop_changes)
# [1] 1657   14
dim(unique(line_prop_changes))
# [1] 1657   14

head(line_prop_changes, 3)
#         line well log_fraction_mean prop_unadjusted_max_value prop_unadjusted_min_value   pool condition treatment replicate  donor
# 1 Arlene-003  366        0.08296595                 0.1583728                 0.1455170 pool15      +C5a untreated      rep1 Arlene
# 2 Arlene-003  367        0.08246210                 0.1593949                 0.1465331 pool15      +C5a untreated      rep2 Arlene
# 3 Arlene-003  370        0.00842465                 0.1473069                 0.1460452 pool15      +C5a       IFN      rep1 Arlene
#   centering_value log_fraction_centered_media_only scaled_log_fraction    sex
# 1      -0.3913577                        0.4743236           0.2831531 Female
# 2      -0.3913577                        0.4738198           0.4125250 Female
# 3       0.1424393                       -0.1340146          -0.2482179 Female

## Num of lines, donors, pools, and wells
length(unique(line_prop_changes$line))
# [1] 216
length(unique(line_prop_changes$donor))
# [1] 209
length(unique(line_prop_changes$pool))
# [1] 15
length(unique(line_prop_changes$well))
# [1] 79

## All in +C5a (fractions already normalized to -C5a)
table(line_prop_changes$condition)
# +C5a 
# 1657 

## Samples = line x pool x treatment x well
dim(unique(line_prop_changes[, c("line", "pool", "treatment", "well")]))
# [1] 1657    4

## 1-3 replicates per sample
line_prop_changes$line_pool_treatment <- paste(line_prop_changes$line, line_prop_changes$pool, line_prop_changes$treatment)
table(line_prop_changes$line_pool_treatment) %>%  table
# 1   2   3 
# 53 592 140 

## Subset to samples with min prop > 0.005
line_prop_changes_filt <- subset(line_prop_changes, prop_unadjusted_min_value > 0.005)
dim(line_prop_changes_filt)
# [1] 1221   15

## Lines left
length(unique(line_prop_changes_filt$line))
# [1] 175

## Mean of scaled log-fractions and min props across replicates 
line_prop_changes_across_replicates <- line_prop_changes_filt %>% group_by(line, pool, treatment, line_pool_treatment) %>% 
                                            summarise(migration_mean_scaled_log_fraction = mean(scaled_log_fraction), 
                                                      migration_mean_min_prop = mean(prop_unadjusted_min_value)) %>% as.data.frame()
  
## % samples in rse with migration fractions
length(which(rse$line_pool_treatment %in% line_prop_changes_across_replicates$line_pool_treatment)) / length(rse$line_pool_treatment) * 100
# [1]  53.97208

## num lines in rse with migration fractions
rse$line[which(rse$line %in% line_prop_changes_across_replicates$line)] %>% unique %>% length()
# 175 / 250

## Add 
metadata <- colData(rse) %>% as.data.frame() %>% left_join(., line_prop_changes_across_replicates)
rse$migration_mean_scaled_log_fraction <- metadata$migration_mean_scaled_log_fraction
rse$migration_mean_min_prop <- metadata$migration_mean_min_prop


######################################
##  Final rse for migration DGE 
###################################### 

## All pseudobulk samples 
dim(rse)
# [1] 18648  2077

## Non-prolif cluster
rse_non_prolif <- rse[, which(rse$proliferation_status == "Not_proliferating")]
dim(rse_non_prolif)
# [1] 18648   1244

## - samples with <100 cells
rse_non_prolif_filt <- rse_non_prolif[, rse_non_prolif$ncells >=100]
dim(rse_non_prolif_filt)
# [1] 18648   781

## Final samples (migration data)
rse_non_prolif_filt_migration <- rse_non_prolif_filt[, which(!is.na(rse_non_prolif_filt$migration_mean_scaled_log_fraction))]
dim(rse_non_prolif_filt_migration)
# [1] 18648   533

length(unique(rse_non_prolif_filt_migration$line))
# [1] 166


## Samples x treatment
rse_non_prolif_filt_migration_IFN <- rse_non_prolif_filt_migration[, which(rse_non_prolif_filt_migration$treatment == "IFN")]
ncol(rse_non_prolif_filt_migration_IFN)
# [1] 123

rse_non_prolif_filt_migration_LPS <- rse_non_prolif_filt_migration[, which(rse_non_prolif_filt_migration$treatment == "LPS")]
ncol(rse_non_prolif_filt_migration_LPS)
# [1] 199

rse_non_prolif_filt_migration_untreated <- rse_non_prolif_filt_migration[, which(rse_non_prolif_filt_migration$treatment == "untreated")]
ncol(rse_non_prolif_filt_migration_untreated)
# [1] 211


## Num of replicates per sample x treatment in final rse
line_prop_changes_in_rse <- subset(line_prop_changes_filt, line_pool_treatment %in% rse_non_prolif_filt_migration$line_pool_treatment)
by(line_prop_changes_in_rse$line_pool_treatment, line_prop_changes_in_rse$treatment, function(x){table(table(x))})
# line_prop_changes_in_rse$treatment: IFN
#  1   2 
# 29  94 
# ----------------------------------------------------------------------------------------------------- 
# line_prop_changes_in_rse$treatment: LPS
#  1   2   3 
#  2 161  36 
# ----------------------------------------------------------------------------------------------------- 
# line_prop_changes_in_rse$treatment: untreated
# 1   2   3 
# 6 174  31 


## Add 100 permutations of migration scaled log-fractions x treatment
set.seed(01152025)
colData(rse_non_prolif_filt_migration_IFN)[, paste0("permutation_migration_mean_scaled_log_fraction_", 1:100)] <- replicate(n = 100, 
                                                                                                                    expr = sample(rse_non_prolif_filt_migration_IFN$migration_mean_scaled_log_fraction, replace = FALSE))

colData(rse_non_prolif_filt_migration_LPS)[, paste0("permutation_migration_mean_scaled_log_fraction_", 1:100)] <- replicate(n = 100, 
                                                                                                                    expr = sample(rse_non_prolif_filt_migration_LPS$migration_mean_scaled_log_fraction, replace = FALSE))

colData(rse_non_prolif_filt_migration_untreated)[, paste0("permutation_migration_mean_scaled_log_fraction_", 1:100)] <- replicate(n = 100, 
                                                                                                                          expr = sample(rse_non_prolif_filt_migration_untreated$migration_mean_scaled_log_fraction, replace = FALSE))

## Explore lines per treatment 
by(rse_non_prolif_filt_migration$line, rse_non_prolif_filt_migration$treatment, function(x){length(unique(x))})
# rse_non_prolif_filt_migration$treatment: IFN
# [1] 100
# ----------------------------------------------------------------------------------------------------- 
# rse_non_prolif_filt_migration$treatment: LPS
# [1] 156
# ----------------------------------------------------------------------------------------------------- 
# rse_non_prolif_filt_migration$treatment: untreated
# [1] 157


save(rse, file = paste0(outdir, "/rse_complete_with_migration_data.Rdata"))
save(rse_non_prolif_filt_migration, file = paste0(outdir, "/rse_non_prolif_filt_migration.Rdata"))
save(rse_non_prolif_filt_migration_IFN, file = paste0(outdir, "/rse_non_prolif_filt_migration_IFN.Rdata"))
save(rse_non_prolif_filt_migration_LPS, file = paste0(outdir, "/rse_non_prolif_filt_migration_LPS.Rdata"))
save(rse_non_prolif_filt_migration_untreated, file = paste0(outdir, "/rse_non_prolif_filt_migration_untreated.Rdata"))







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
# date     2025-01-15
# rstudio  2024.04.0+735 Chocolate Cosmos (server)
# pandoc   3.1.12.3 @ /opt/view/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [1] CRAN (R 4.3.1)
# Biobase              * 2.62.0    2023-10-24 [1] Bioconductor
# BiocGenerics         * 0.48.1    2023-11-01 [1] Bioconductor
# bit                    4.0.5     2022-11-15 [1] CRAN (R 4.3.1)
# bit64                  4.0.5     2020-08-30 [1] CRAN (R 4.3.1)
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
# rstudioapi             0.16.0    2024-03-24 [1] CRAN (R 4.3.1)
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
# vroom                  1.6.5     2023-12-05 [1] CRAN (R 4.3.1)
# withr                  3.0.0     2024-01-16 [1] CRAN (R 4.3.1)
# XVector                0.42.0    2023-10-24 [1] Bioconductor
# zlibbioc               1.48.2    2024-03-13 [1] Bioconductor 3.18 (R 4.3.1)
# 
# [1] /opt/view/rlib/R/library
# [2] /software/hgi/installs/softpack/rstudio/.spack/linux-ubuntu22.04-x86_64_v3/gcc-11.4.0/r-4.3.1-bfwldrk76z6f52upk47zepliekn7ayqz/rlib/R/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────