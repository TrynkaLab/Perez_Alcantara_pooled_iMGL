
library(tidyverse)
library(dplyr)
library(SummarizedExperiment)
library(rlang)
library(sessioninfo)


################################################################################
##                7. Permutation analysis for proliferation DGE
################################################################################

#-------------------------------------------------------------------------------
#                     7.0 Data exploration and processing
#-------------------------------------------------------------------------------
#  Code to explore, process, and add cell line proliferation data to sample
#  metadata to run DGE on permuted phenotypes in the next script.
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

## Set working dir
setwd("/lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_Daianna/")

## Define output data and plot dir
outdir = paste(getwd(), "output_data", "07_Diff_expr_proliferation_permutations", "00_Data_exploration_processing", sep = "/")
plotdir = paste(getwd(), "plots", "07_Diff_expr_proliferation_permutations", "00_Data_exploration_processing", sep = "/")
dir.create(outdir, recursive = T)
dir.create(plotdir, recursive = T)

## Input dirs
input_dir06 = paste(getwd(), "output_data", "06_Diff_expr_migration_permutations", "00_Data_exploration_processing", sep = "/")

## Load complete rse (with genotype PCs)
load(paste0(input_dir06, "/rse_complete_with_migration_data.Rdata"), verbose = T)

## Load prolif log-fractions per donor x pool:
#   from preMac (diff ages) -> microglia (in IFN/LPS/untreated):
line_prop_changes_microglia_vs_premac <- as.data.frame(read_csv("/lustre/scratch123/hgi/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/data/results/1.2.scale_proliferation/line_prop_changes_microglia_premac.csv"))

  
#################################################
##   Add proliferation data to sample metadata
#################################################

########### (young - old) preMac -> microglia (in IFN/LPS/untreated) ###########

## Samples = line x pool x batch x preMAC_age x treatment 
dim(line_prop_changes_microglia_vs_premac)
# [1] 2344   12
dim(unique(line_prop_changes_microglia_vs_premac[, c("line", "pool", "differentiation", "preMAC_age", "treatment")]))
# [1] 2344    5

## Treatments
table(line_prop_changes_microglia_vs_premac$treatment)
# IFNg       LPS untreated 
#  823       957       564 
line_prop_changes_microglia_vs_premac$treatment  <- 
        line_prop_changes_microglia_vs_premac$treatment %>%  replace(., which(. == "IFNg"), "IFN")

## Samples at the required level
line_prop_changes_microglia_vs_premac$line_pool_treatment <- paste(line_prop_changes_microglia_vs_premac$line, 
                                                                   line_prop_changes_microglia_vs_premac$pool, 
                                                                   line_prop_changes_microglia_vs_premac$treatment)
length(unique(line_prop_changes_microglia_vs_premac$line_pool_treatment ))
# [1] 1091

## preMac ages
range(line_prop_changes_microglia_vs_premac$preMAC_age)
# [1] 28 60

## Num of entries per sample
table(line_prop_changes_microglia_vs_premac$line_pool_treatment) %>%  table
#   1   2   3   4 
# 260 430 380  21 

## Aggregate prolif estimates across batches and preMac ages
line_prop_changes_microglia_vs_premac_agg_per_sample <- line_prop_changes_microglia_vs_premac %>% 
                                                            group_by(line_pool_treatment, line, pool, treatment) %>% 
                                                            summarize(proliferation_mean_scaled_log_fraction = mean(scaled_log_fraction),
                                                                      proliferation_mean_min_prop = mean(prop_unadjusted_min_value))

## % of samples in rse with prolif data
length(which(rse$line_pool_treatment %in% line_prop_changes_microglia_vs_premac_agg_per_sample$line_pool_treatment)) / length(rse$line_pool_treatment) * 100
# [1] 86.27829 (1792 / 2077)

## Add to sample metadata 
complete_metadata <- left_join(as.data.frame(colData(rse)), line_prop_changes_microglia_vs_premac_agg_per_sample)
rse$proliferation_mean_scaled_log_fraction <- complete_metadata$proliferation_mean_scaled_log_fraction
rse$proliferation_mean_min_prop <- complete_metadata$proliferation_mean_min_prop

save(rse, file = paste0(outdir, "/rse_complete_with_proliferation_data.Rdata"))


## Non-prolif cluster
rse_non_prolif <- rse[, which(rse$proliferation_status == "Not_proliferating")]
dim(rse_non_prolif)
# [1] 18648   1244

## - samples with <100 cells
rse_non_prolif_filt <- rse_non_prolif[, rse_non_prolif$ncells >=100]
dim(rse_non_prolif_filt)
# [1] 18648   781

## Final samples (prolif data)
rse_non_prolif_filt_proliferation <- rse_non_prolif_filt[, which(!is.na(rse_non_prolif_filt$proliferation_mean_scaled_log_fraction))]
dim(rse_non_prolif_filt_proliferation)
# [1] 18648   757

length(unique(rse_non_prolif_filt_proliferation$line))
# [1] 193


## Samples x treatment
rse_non_prolif_filt_proliferation_IFN <- rse_non_prolif_filt_proliferation[, which(rse_non_prolif_filt_proliferation$treatment == "IFN")]
ncol(rse_non_prolif_filt_proliferation_IFN)
# [1] 257

rse_non_prolif_filt_proliferation_LPS <- rse_non_prolif_filt_proliferation[, which(rse_non_prolif_filt_proliferation$treatment == "LPS")]
ncol(rse_non_prolif_filt_proliferation_LPS)
# [1] 245

rse_non_prolif_filt_proliferation_untreated <- rse_non_prolif_filt_proliferation[, which(rse_non_prolif_filt_proliferation$treatment == "untreated")]
ncol(rse_non_prolif_filt_proliferation_untreated)
# [1] 255

## Explore lines per treatment 
by(rse_non_prolif_filt_proliferation$line, rse_non_prolif_filt_proliferation$treatment, function(x){length(unique(x))})
# rse_non_prolif_filt_proliferation$treatment: IFN
# [1] 191
# ----------------------------------------------------------------------------------------------------------------- 
# rse_non_prolif_filt_proliferation$treatment: LPS
# [1] 188
# ----------------------------------------------------------------------------------------------------------------- 
# rse_non_prolif_filt_proliferation$treatment: untreated
# [1] 183


## Add 100 permutations of prolif scaled log-fractions x treatment
set.seed(01172025)
colData(rse_non_prolif_filt_proliferation_IFN)[, paste0("permutation_proliferation_mean_scaled_log_fraction_", 1:100)] <- replicate(n = 100, 
                                                                                                                                    expr = sample(rse_non_prolif_filt_proliferation_IFN$proliferation_mean_scaled_log_fraction, replace = FALSE))

colData(rse_non_prolif_filt_proliferation_LPS)[, paste0("permutation_proliferation_mean_scaled_log_fraction_", 1:100)] <- replicate(n = 100, 
                                                                                                                                    expr = sample(rse_non_prolif_filt_proliferation_LPS$proliferation_mean_scaled_log_fraction, replace = FALSE))

colData(rse_non_prolif_filt_proliferation_untreated)[, paste0("permutation_proliferation_mean_scaled_log_fraction_", 1:100)] <- replicate(n = 100, 
                                                                                                                                          expr = sample(rse_non_prolif_filt_proliferation_untreated$proliferation_mean_scaled_log_fraction, replace = FALSE))

save(rse_non_prolif_filt_proliferation, file = paste0(outdir, "/rse_non_prolif_filt_proliferation.Rdata"))
save(rse_non_prolif_filt_proliferation_IFN, file = paste0(outdir, "/rse_non_prolif_filt_proliferation_IFN.Rdata"))
save(rse_non_prolif_filt_proliferation_LPS, file = paste0(outdir, "/rse_non_prolif_filt_proliferation_LPS.Rdata"))
save(rse_non_prolif_filt_proliferation_untreated, file = paste0(outdir, "/rse_non_prolif_filt_proliferation_untreated.Rdata"))







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
# pandoc   3.1.12.3 @ /opt/view/bin/ (via rmarkdown)
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
# digest                 0.6.35    2024-03-11 [1] CRAN (R 4.3.1)
# dplyr                * 1.1.4     2023-11-17 [1] CRAN (R 4.3.1)
# evaluate               0.23      2023-11-01 [1] CRAN (R 4.3.1)
# fansi                  1.0.6     2023-12-08 [1] CRAN (R 4.3.1)
# fastmap                1.1.1     2023-02-24 [1] CRAN (R 4.3.1)
# forcats              * 1.0.0     2023-01-29 [1] CRAN (R 4.3.1)
# generics               0.1.3     2022-07-05 [1] CRAN (R 4.3.1)
# GenomeInfoDb         * 1.38.8    2024-03-15 [1] Bioconductor 3.18 (R 4.3.1)
# GenomeInfoDbData       1.2.11    2024-11-06 [1] Bioconductor
# GenomicRanges        * 1.54.1    2023-10-29 [1] Bioconductor
# ggplot2              * 3.5.1     2024-04-23 [1] CRAN (R 4.3.1)
# ggrepel                0.9.5     2024-01-10 [1] CRAN (R 4.3.1)
# glue                   1.7.0     2024-01-09 [1] CRAN (R 4.3.1)
# gtable                 0.3.4     2023-08-21 [1] CRAN (R 4.3.1)
# hms                    1.1.3     2023-03-21 [1] CRAN (R 4.3.1)
# htmltools              0.5.8     2024-03-25 [1] CRAN (R 4.3.1)
# IRanges              * 2.36.0    2023-10-24 [1] Bioconductor
# jtools                 2.2.2     2023-07-11 [1] CRAN (R 4.3.1)
# knitr                  1.45      2023-10-30 [1] CRAN (R 4.3.1)
# lattice                0.22-6    2024-03-20 [1] CRAN (R 4.3.1)
# lifecycle              1.0.4     2023-11-07 [1] CRAN (R 4.3.1)
# lubridate            * 1.9.3     2023-09-27 [1] CRAN (R 4.3.1)
# magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.3.1)
# Matrix                 1.6-5     2024-01-11 [1] CRAN (R 4.3.1)
# MatrixGenerics       * 1.14.0    2023-10-24 [1] Bioconductor
# matrixStats          * 1.2.0     2023-12-11 [1] CRAN (R 4.3.1)
# munsell                0.5.1     2024-04-01 [1] CRAN (R 4.3.1)
# pander                 0.6.5     2022-03-18 [1] CRAN (R 4.3.1)
# pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.3.1)
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.1)
# purrr                * 1.0.2     2023-08-10 [1] CRAN (R 4.3.1)
# R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.1)
# Rcpp                   1.0.12    2024-01-09 [1] CRAN (R 4.3.1)
# RCurl                  1.98-1.14 2024-01-09 [1] CRAN (R 4.3.1)
# readr                * 2.1.5     2024-01-10 [1] CRAN (R 4.3.1)
# rlang                * 1.1.3     2024-01-10 [1] CRAN (R 4.3.1)
# rmarkdown              2.26      2024-03-05 [1] CRAN (R 4.3.1)
# rstudioapi             0.16.0    2024-03-24 [1] CRAN (R 4.3.1)
# S4Arrays               1.2.1     2024-03-04 [1] Bioconductor 3.18 (R 4.3.1)
# S4Vectors            * 0.40.2    2023-11-23 [1] Bioconductor 3.18 (R 4.3.1)
# scales                 1.3.0     2023-11-28 [1] CRAN (R 4.3.1)
# sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.3.1)
# SparseArray            1.2.4     2024-02-11 [1] Bioconductor 3.18 (R 4.3.1)
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
# xfun                   0.43      2024-03-25 [1] CRAN (R 4.3.1)
# XVector                0.42.0    2023-10-24 [1] Bioconductor
# yaml                   2.3.8     2023-12-11 [1] CRAN (R 4.3.1)
# zlibbioc               1.48.2    2024-03-13 [1] Bioconductor 3.18 (R 4.3.1)
# 
# [1] /opt/view/rlib/R/library
# [2] /software/hgi/installs/softpack/rstudio/.spack/linux-ubuntu22.04-x86_64_v3/gcc-11.4.0/r-4.3.1-bfwldrk76z6f52upk47zepliekn7ayqz/rlib/R/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────


