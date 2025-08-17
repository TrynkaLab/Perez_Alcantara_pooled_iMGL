
library(tidyverse)
library(dplyr)
library(rlang)
library(sessioninfo)


################################################################################
##      8. Genome Wide Association Study for proliferation phenotype
################################################################################

#-------------------------------------------------------------------------------
#                     8.0 Data exploration and processing
#-------------------------------------------------------------------------------
#  Code to explore and process genotype and proliferation data, preparing input 
#  data for proliferation GWAS for each variant chunk.
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

## Define dirs
outdir = paste("..", "..", "output_data", "08_GWAS_proliferation", "00_Data_exploration_processing", sep = "/")
plotdir = paste("..", "..", "plots", "08_GWAS_proliferation", "00_Data_exploration_processing", sep = "/")
logs_dir = paste("logs", "00_Data_exploration_processing", sep = "/")
err_dir = paste("errors", "00_Data_exploration_processing", sep = "/")
dir.create(outdir, recursive = T)
dir.create(plotdir, recursive = T)
dir.create(logs_dir, recursive = T)
dir.create(err_dir, recursive = T)

## Input dirs 
genotype_input_dir = "/lustre/scratch125/humgen/projects_v2/otar2065/OTAR2065_phenotypic_QTLs/data/genotypes/full_genotype/"
input_dir04 = paste("..", "..", "output_data", "04_Burden_test_proliferation", "01_Burden_tests", sep = "/")


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#                          8.0.1 Prepare input data
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

## Total number of chunks
grep("genotype_minor_allele_dosage_", list.files(genotype_input_dir)) %>% length
# [1] 582

## Test with 1st chunk
# chunk = 1
# genotype <- read_csv(genotype_input_dir) %>%
#   as.data.frame()

args = commandArgs(trailingOnly = TRUE)

if (length(args)<6) {
  stop("You need to detail 1 genotype input and 5 output file paths.n", call. = FALSE)
} else if (length(args) == 6) {
  genotype_path = args[1]
  line_prop_changes_premac_vs_iPSC_genotype_out_path = args[2]
  line_prop_changes_old_vs_young_premac_genotype_out_path = args[3]
  line_prop_changes_microglia_vs_premac_IFN_genotype_out_path = args[4]
  line_prop_changes_microglia_vs_premac_LPS_genotype_out_path = args[5]
  line_prop_changes_microglia_vs_premac_untreated_genotype_out_path = args[6]
}

## Extract chunk
chunk = genotype_path %>% strsplit(., ".csv") %>% .[[1]] %>% strsplit(., "_") %>% .[[1]] %>% .[length(.)]

## 9999 variants per chunk
genotype <- read_csv(genotype_path)
dim(genotype)
# [1] 9999  251

## Doses:
genotype[, -251] %>% as.vector() %>% unlist() %>% unique()
# [1] 1.0 0.0 0.5

## All SNPs? No

## Duplicated doses (not independent variants) 
which(duplicated(genotype[, -251])) %>% length() 
# [1] 5136      (example for chunk 1)
genotype[1:2, c(1:10, 251)] # (example for chunk 1)
# Gilma-009 Keoni-014 Mindy-016 Dexter-006 Fiona-010 Hector-011 Qiana-022 Imani-012 Nestor-017 Arlene-003             rn
#         1         1         1        0.5       0.5          1       0.5         1           1          1  1_817186_G_A
#         1         1         1        0.5       0.5          1       0.5         1           1          1  1_817341_A_G

## Recode allele dosages
genotype[genotype == 1] <- 2
genotype[genotype == 0.5] <- 1


################################################################################
##                    iPSC -> young preMac proliferation 
################################################################################

load(paste0(input_dir04, "/line_prop_changes_premac_vs_iPSC.Rdata"), verbose = T)

## Mean prolif estimate per line across pools
## Num pools each line is present in
line_prop_changes_premac_vs_iPSC_mean <- line_prop_changes_premac_vs_iPSC %>% 
                                            group_by(line) %>% 
                                            mutate(mean_scaled_log_fraction = mean(scaled_log_fraction),
                                                   npool = length(unique(pool))) %>% 
                                            ungroup() %>% 
                                            dplyr::select(line, mean_scaled_log_fraction, npool, sex, paste0("PC", 1:10)) %>% 
                                            unique()
## Num lines
dim(line_prop_changes_premac_vs_iPSC_mean)
# [1] 228   14

## Subset to lines with prolif estimates
if(length(which(! line_prop_changes_premac_vs_iPSC_mean$line %in% colnames(genotype))) > 0){
  stop("Not all lines with estimates of proliferation from iPSC to preMac have genotype data.")
}
genotype_prolif_premac_vs_iPSC <- genotype[, which(colnames(genotype) %in% 
                                                     c("rn", line_prop_changes_premac_vs_iPSC_mean$line))]
message(paste("There are", ncol(genotype_prolif_premac_vs_iPSC)-1, "lines with genotype info and estimates of prolif from iPSC -> preMac"))
# There are 228 lines with genotype info and estimates of prolif from iPSC -> preMac

## (almost) no variants with MAF<5% 
genotype_prolif_premac_vs_iPSC_MAFs <- genotype_prolif_premac_vs_iPSC %>% 
          mutate(MAF = genotype_prolif_premac_vs_iPSC[, -229] %>% rowSums(.) / (ncol(genotype_prolif_premac_vs_iPSC[, -229]) *2))

range(genotype_prolif_premac_vs_iPSC_MAFs$MAF)
# [1]  0.04385965 0.95614035

## Merge prolif + cov + genotype
genotype_prolif_premac_vs_iPSC_t <- t(genotype_prolif_premac_vs_iPSC[,-229]) %>% as.data.frame()
colnames(genotype_prolif_premac_vs_iPSC_t) <- genotype_prolif_premac_vs_iPSC$rn
genotype_prolif_premac_vs_iPSC_t$line <- rownames(genotype_prolif_premac_vs_iPSC_t)

line_prop_changes_premac_vs_iPSC_genotype <- line_prop_changes_premac_vs_iPSC_mean %>% 
                    inner_join(., genotype_prolif_premac_vs_iPSC_t)

## Make sure all variants in chunk are included
if (chunk != 582) {
  if(nrow(line_prop_changes_premac_vs_iPSC_genotype) == 228 & 
     ncol(line_prop_changes_premac_vs_iPSC_genotype) == 9999 + 14){
    save(line_prop_changes_premac_vs_iPSC_genotype, file = line_prop_changes_premac_vs_iPSC_genotype_out_path)
  }
  else{
    stop(paste0("Not all variants or samples are present in input data for chunk ", chunk, "."))
  }
} else{
  if(nrow(line_prop_changes_premac_vs_iPSC_genotype) == 228 &
     ncol(line_prop_changes_premac_vs_iPSC_genotype) == 1888 + 14){
    save(line_prop_changes_premac_vs_iPSC_genotype, file = line_prop_changes_premac_vs_iPSC_genotype_out_path)
  }
  else{
    stop("Not all variants or samples are present in input data for last chunk.")
  }
}


################################################################################
##                 young preMac -> old preMac proliferation
################################################################################

load(paste0(input_dir04, "/line_prop_changes_old_vs_young_premac.Rdata"), verbose = T)

line_prop_changes_old_vs_young_premac_mean <- line_prop_changes_old_vs_young_premac %>% 
  group_by(line) %>% 
  mutate(mean_scaled_log_fraction = mean(scaled_log_fraction),
         npool = length(unique(pool))) %>% 
  ungroup() %>% 
  dplyr::select(line, mean_scaled_log_fraction, npool, sex, paste0("PC", 1:10)) %>% 
  unique()

## Num lines
dim(line_prop_changes_old_vs_young_premac_mean)
# [1] 210   14

## Subset to lines with prolif estimates
if(length(which(! line_prop_changes_old_vs_young_premac_mean$line %in% colnames(genotype))) > 0){
  stop("Not all lines with estimates of proliferation from young to old preMac have genotype data.")
}
genotype_prolif_old_vs_young_premac <- genotype[, which(colnames(genotype) %in% 
                                                     c("rn", line_prop_changes_old_vs_young_premac_mean$line))]
message(paste("There are", ncol(genotype_prolif_old_vs_young_premac)-1, "lines with genotype info and estimates of prolif from young -> old preMac"))
# There are 210 lines with genotype info and estimates of prolif from young -> old preMac

## Almost no variants with MAF<5% 
genotype_prolif_old_vs_young_premac_MAFs <- genotype_prolif_old_vs_young_premac %>% 
  mutate(MAF = genotype_prolif_old_vs_young_premac[, -211] %>% rowSums(.) / (ncol(genotype_prolif_old_vs_young_premac[, -211]) *2))

range(genotype_prolif_old_vs_young_premac_MAFs$MAF)
# [1] 0.04285714 0.94285714

## Merge prolif + cov + genotype
genotype_prolif_old_vs_young_premac_t <- t(genotype_prolif_old_vs_young_premac[,-211]) %>% as.data.frame()
colnames(genotype_prolif_old_vs_young_premac_t) <- genotype_prolif_old_vs_young_premac$rn
genotype_prolif_old_vs_young_premac_t$line <- rownames(genotype_prolif_old_vs_young_premac_t)

line_prop_changes_old_vs_young_premac_genotype <- line_prop_changes_old_vs_young_premac_mean %>% 
  inner_join(., genotype_prolif_old_vs_young_premac_t)

if (chunk != 582) {
  if(nrow(line_prop_changes_old_vs_young_premac_genotype) == 210 &
     ncol(line_prop_changes_old_vs_young_premac_genotype) == 9999 + 14){
    save(line_prop_changes_old_vs_young_premac_genotype, file = line_prop_changes_old_vs_young_premac_genotype_out_path)
  }
  else{
    stop(paste0("Not all variants or samples are present in input data for chunk ", chunk, "."))
  }
} else{
  if(nrow(line_prop_changes_old_vs_young_premac_genotype) == 210 &
     ncol(line_prop_changes_old_vs_young_premac_genotype) == 1888 + 14){
    save(line_prop_changes_old_vs_young_premac_genotype, file = line_prop_changes_old_vs_young_premac_genotype_out_path)
  }
  else{
    stop("Not all variants or samples are present in input data for last chunk.")
  }
}


################################################################################
##           preMac -> microglia (in IFN/LPS/untreated) proliferation
################################################################################

load(paste0(input_dir04, "/line_prop_changes_microglia_vs_premac.Rdata"), verbose = T)

#############################   Microglia in IFN   #############################

line_prop_changes_microglia_vs_premac_mean_IFN <- line_prop_changes_microglia_vs_premac %>% 
  filter(treatment == "IFNg") %>% 
  group_by(line) %>% 
  mutate(mean_scaled_log_fraction = mean(scaled_log_fraction),
         npool = length(unique(pool))) %>% 
  ungroup() %>% 
  dplyr::select(line, mean_scaled_log_fraction, npool, sex, paste0("PC", 1:10)) %>% 
  unique()

## Num lines
dim(line_prop_changes_microglia_vs_premac_mean_IFN)
# [1] 230   14

## Subset to lines with prolif estimates
if(length(which(! line_prop_changes_microglia_vs_premac_mean_IFN$line %in% colnames(genotype))) > 0){
  stop("Not all lines with estimates of proliferation from preMac to microglia in IFN have genotype data.")
}
genotype_prolif_microglia_vs_premac <- genotype[, which(colnames(genotype) %in% 
                                                          c("rn", line_prop_changes_microglia_vs_premac_mean_IFN$line))]
message(paste("There are", ncol(genotype_prolif_microglia_vs_premac)-1, "lines with genotype info and estimates of prolif from preMac -> microglia in IFN"))
# There are 230 lines with genotype info and estimates of prolif from preMac -> microglia in IFN

## No variants with MAF<5% 
genotype_prolif_microglia_vs_premac_MAFs <- genotype_prolif_microglia_vs_premac %>% 
  mutate(MAF = genotype_prolif_microglia_vs_premac[, -231] %>% rowSums(.) / (ncol(genotype_prolif_microglia_vs_premac[, -231]) *2))

range(genotype_prolif_microglia_vs_premac_MAFs$MAF)
# [1] 0.0500000 0.9456522

## Merge prolif + cov + genotype
genotype_prolif_microglia_vs_premac_t <- t(genotype_prolif_microglia_vs_premac[,-231]) %>% as.data.frame()
colnames(genotype_prolif_microglia_vs_premac_t) <- genotype_prolif_microglia_vs_premac$rn
genotype_prolif_microglia_vs_premac_t$line <- rownames(genotype_prolif_microglia_vs_premac_t)

line_prop_changes_microglia_vs_premac_IFN_genotype <- line_prop_changes_microglia_vs_premac_mean_IFN %>% 
  inner_join(., genotype_prolif_microglia_vs_premac_t)

if (chunk != 582) {
  if(nrow(line_prop_changes_microglia_vs_premac_IFN_genotype) == 230 & 
     ncol(line_prop_changes_microglia_vs_premac_IFN_genotype) == 9999 + 14){
    save(line_prop_changes_microglia_vs_premac_IFN_genotype, file = line_prop_changes_microglia_vs_premac_IFN_genotype_out_path)
  }
  else{
    stop(paste0("Not all variants or samples are present in input data for chunk ", chunk, "."))
  }
} else{
  if(nrow(line_prop_changes_microglia_vs_premac_IFN_genotype) == 230 &
     ncol(line_prop_changes_microglia_vs_premac_IFN_genotype) == 1888 + 14){
    save(line_prop_changes_microglia_vs_premac_IFN_genotype, file = line_prop_changes_microglia_vs_premac_IFN_genotype_out_path)
  }
  else{
    stop("Not all variants or samples are present in input data for last chunk.")
  }
}


#############################   Microglia in LPS   #############################

line_prop_changes_microglia_vs_premac_mean_LPS <- line_prop_changes_microglia_vs_premac %>% 
  filter(treatment == "LPS") %>% 
  group_by(line) %>% 
  mutate(mean_scaled_log_fraction = mean(scaled_log_fraction),
         npool = length(unique(pool))) %>% 
  ungroup() %>% 
  dplyr::select(line, mean_scaled_log_fraction, npool, sex, paste0("PC", 1:10)) %>% 
  unique()

## Num lines
dim(line_prop_changes_microglia_vs_premac_mean_LPS)
# [1] 229   14

## Subset to lines with prolif estimates
if(length(which(! line_prop_changes_microglia_vs_premac_mean_LPS$line %in% colnames(genotype))) > 0){
  stop("Not all lines with estimates of proliferation from preMac to microglia in LPS have genotype data.")
}
genotype_prolif_microglia_vs_premac <- genotype[, which(colnames(genotype) %in% 
                                                          c("rn", line_prop_changes_microglia_vs_premac_mean_LPS$line))]
message(paste("There are", ncol(genotype_prolif_microglia_vs_premac)-1, "lines with genotype info and estimates of prolif from preMac -> microglia in LPS"))
# There are 229 lines with genotype info and estimates of prolif from preMac -> microglia in LPS

## No variants with MAF<5% 
genotype_prolif_microglia_vs_premac_MAFs <- genotype_prolif_microglia_vs_premac %>% 
  mutate(MAF = genotype_prolif_microglia_vs_premac[, -230] %>% rowSums(.) / (ncol(genotype_prolif_microglia_vs_premac[, -230]) *2))

range(genotype_prolif_microglia_vs_premac_MAFs$MAF)
# [1] 0.05021834 0.94541485

## Merge prolif + cov + genotype
genotype_prolif_microglia_vs_premac_t <- t(genotype_prolif_microglia_vs_premac[,-230]) %>% as.data.frame()
colnames(genotype_prolif_microglia_vs_premac_t) <- genotype_prolif_microglia_vs_premac$rn
genotype_prolif_microglia_vs_premac_t$line <- rownames(genotype_prolif_microglia_vs_premac_t)

line_prop_changes_microglia_vs_premac_LPS_genotype <- line_prop_changes_microglia_vs_premac_mean_LPS %>% 
  inner_join(., genotype_prolif_microglia_vs_premac_t)

if (chunk != 582) {
  if(nrow(line_prop_changes_microglia_vs_premac_LPS_genotype) == 229 & 
     ncol(line_prop_changes_microglia_vs_premac_LPS_genotype) == 9999 + 14){
    save(line_prop_changes_microglia_vs_premac_LPS_genotype, file = line_prop_changes_microglia_vs_premac_LPS_genotype_out_path)
  }
  else{
    stop(paste0("Not all variants or samples are present in input data for chunk ", chunk, "."))
  }
} else{
  if(nrow(line_prop_changes_microglia_vs_premac_LPS_genotype) == 229 & 
     ncol(line_prop_changes_microglia_vs_premac_LPS_genotype) == 1888 + 14){
    save(line_prop_changes_microglia_vs_premac_LPS_genotype, file = line_prop_changes_microglia_vs_premac_LPS_genotype_out_path)
  }
  else{
    stop("Not all variants or samples are present in input data for last chunk.")
  }
}


##########################   Microglia in untreated   ##########################

line_prop_changes_microglia_vs_premac_mean_untreated <- line_prop_changes_microglia_vs_premac %>% 
  filter(treatment == "untreated") %>% 
  group_by(line) %>% 
  mutate(mean_scaled_log_fraction = mean(scaled_log_fraction),
         npool = length(unique(pool))) %>% 
  ungroup() %>% 
  dplyr::select(line, mean_scaled_log_fraction, npool, sex, paste0("PC", 1:10)) %>% 
  unique()

## Num lines
dim(line_prop_changes_microglia_vs_premac_mean_untreated)
# [1] 209   14

## Subset to lines with prolif estimates
if(length(which(! line_prop_changes_microglia_vs_premac_mean_untreated$line %in% colnames(genotype))) > 0){
  stop("Not all lines with estimates of proliferation from preMac to microglia in untreated have genotype data.")
}
genotype_prolif_microglia_vs_premac <- genotype[, which(colnames(genotype) %in% 
                                                          c("rn", line_prop_changes_microglia_vs_premac_mean_untreated$line))]
message(paste("There are", ncol(genotype_prolif_microglia_vs_premac)-1, "lines with genotype info and estimates of prolif from preMac -> microglia in untreated"))
# There are 209 lines with genotype info and estimates of prolif from preMac -> microglia in untreated

## MAF<5% 
genotype_prolif_microglia_vs_premac_MAFs <- genotype_prolif_microglia_vs_premac %>% 
  mutate(MAF = genotype_prolif_microglia_vs_premac[, -210] %>% rowSums(.) / (ncol(genotype_prolif_microglia_vs_premac[, -210]) *2))

range(genotype_prolif_microglia_vs_premac_MAFs$MAF)
# [1] 0.04545455 0.94258373

## Merge prolif + cov + genotype
genotype_prolif_microglia_vs_premac_t <- t(genotype_prolif_microglia_vs_premac[,-210]) %>% as.data.frame()
colnames(genotype_prolif_microglia_vs_premac_t) <- genotype_prolif_microglia_vs_premac$rn
genotype_prolif_microglia_vs_premac_t$line <- rownames(genotype_prolif_microglia_vs_premac_t)

line_prop_changes_microglia_vs_premac_untreated_genotype <- line_prop_changes_microglia_vs_premac_mean_untreated %>% 
  inner_join(., genotype_prolif_microglia_vs_premac_t)

if (chunk != 582) {
  if(nrow(line_prop_changes_microglia_vs_premac_untreated_genotype) == 209 &
     ncol(line_prop_changes_microglia_vs_premac_untreated_genotype) == 9999 + 14){
    save(line_prop_changes_microglia_vs_premac_untreated_genotype, file = line_prop_changes_microglia_vs_premac_untreated_genotype_out_path)
  }
  else{
    stop(paste0("Not all variants or samples are present in input data for chunk ", chunk, "."))
  }
} else{
  if(nrow(line_prop_changes_microglia_vs_premac_untreated_genotype) == 209 &
     ncol(line_prop_changes_microglia_vs_premac_untreated_genotype) == 1888 + 14){
    save(line_prop_changes_microglia_vs_premac_untreated_genotype, file = line_prop_changes_microglia_vs_premac_untreated_genotype_out_path)
  }
  else{
    stop("Not all variants or samples are present in input data for last chunk.")
  }
}







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
# date     2025-02-05
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
# bit                    4.0.5      2022-11-15 [1] CRAN (R 4.3.1)
# bit64                  4.0.5      2020-08-30 [1] CRAN (R 4.3.1)
# bitops                 1.0-7      2021-04-24 [1] CRAN (R 4.3.1)
# boot                   1.3-30     2024-02-26 [1] CRAN (R 4.3.1)
# broom                  1.0.5      2023-06-09 [1] CRAN (R 4.3.1)
# caTools                1.18.2     2021-03-28 [1] CRAN (R 4.3.1)
# cli                    3.6.2      2023-12-11 [1] CRAN (R 4.3.1)
# codetools              0.2-20     2024-03-31 [1] CRAN (R 4.3.1)
# colorspace             2.1-0      2023-01-23 [1] CRAN (R 4.3.1)
# corpcor                1.6.10     2021-09-16 [1] CRAN (R 4.3.1)
# crayon                 1.5.2      2022-09-29 [1] CRAN (R 4.3.1)
# DelayedArray           0.28.0     2023-10-24 [1] Bioconductor
# digest                 0.6.35     2024-03-11 [1] CRAN (R 4.3.1)
# dplyr                * 1.1.4      2023-11-17 [1] CRAN (R 4.3.1)
# edgeR                * 4.0.16     2024-02-18 [1] Bioconductor 3.18 (R 4.3.1)
# EnvStats               2.8.1      2023-08-22 [1] CRAN (R 4.3.1)
# fANCOVA                0.6-1      2020-11-13 [1] CRAN (R 4.3.1)
# fansi                  1.0.6      2023-12-08 [1] CRAN (R 4.3.1)
# forcats              * 1.0.0      2023-01-29 [1] CRAN (R 4.3.1)
# generics               0.1.3      2022-07-05 [1] CRAN (R 4.3.1)
# GenomeInfoDb         * 1.38.8     2024-03-15 [1] Bioconductor 3.18 (R 4.3.1)
# GenomeInfoDbData       1.2.11     2024-11-06 [1] Bioconductor
# GenomicRanges        * 1.54.1     2023-10-29 [1] Bioconductor
# ggplot2              * 3.5.1      2024-04-23 [1] CRAN (R 4.3.1)
# glue                   1.7.0      2024-01-09 [1] CRAN (R 4.3.1)
# gplots                 3.1.3.1    2024-02-02 [1] CRAN (R 4.3.1)
# gtable                 0.3.4      2023-08-21 [1] CRAN (R 4.3.1)
# gtools                 3.9.5      2023-11-20 [1] CRAN (R 4.3.1)
# hms                    1.1.3      2023-03-21 [1] CRAN (R 4.3.1)
# IRanges              * 2.36.0     2023-10-24 [1] Bioconductor
# iterators              1.0.14     2022-02-05 [1] CRAN (R 4.3.1)
# jtools                 2.2.2      2023-07-11 [1] CRAN (R 4.3.1)
# KernSmooth             2.23-22    2023-07-10 [1] CRAN (R 4.3.1)
# lattice                0.22-6     2024-03-20 [1] CRAN (R 4.3.1)
# lifecycle              1.0.4      2023-11-07 [1] CRAN (R 4.3.1)
# limma                * 3.58.1     2023-10-31 [1] Bioconductor
# lme4                   1.1-35.2   2024-03-28 [1] CRAN (R 4.3.1)
# lmerTest               3.1-3      2020-10-23 [1] CRAN (R 4.3.1)
# locfit                 1.5-9.9    2024-03-01 [1] CRAN (R 4.3.1)
# lubridate            * 1.9.3      2023-09-27 [1] CRAN (R 4.3.1)
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
# pander                 0.6.5      2022-03-18 [1] CRAN (R 4.3.1)
# pbkrtest               0.5.2      2023-01-19 [1] CRAN (R 4.3.1)
# pillar                 1.9.0      2023-03-22 [1] CRAN (R 4.3.1)
# pkgconfig              2.0.3      2019-09-22 [1] CRAN (R 4.3.1)
# plyr                   1.8.9      2023-10-02 [1] CRAN (R 4.3.1)
# purrr                * 1.0.2      2023-08-10 [1] CRAN (R 4.3.1)
# R6                     2.5.1      2021-08-19 [1] CRAN (R 4.3.1)
# rbibutils              2.2.16     2023-10-25 [1] CRAN (R 4.3.1)
# Rcpp                   1.0.12     2024-01-09 [1] CRAN (R 4.3.1)
# RCurl                  1.98-1.14  2024-01-09 [1] CRAN (R 4.3.1)
# Rdpack                 2.6        2023-11-08 [1] CRAN (R 4.3.1)
# readr                * 2.1.5      2024-01-10 [1] CRAN (R 4.3.1)
# remaCor                0.0.18     2024-02-08 [1] CRAN (R 4.3.1)
# reshape2               1.4.4      2020-04-09 [1] CRAN (R 4.3.1)
# RhpcBLASctl            0.23-42    2023-02-11 [1] CRAN (R 4.3.1)
# rlang                * 1.1.3      2024-01-10 [1] CRAN (R 4.3.1)
# rstudioapi             0.16.0     2024-03-24 [1] CRAN (R 4.3.1)
# S4Arrays               1.2.1      2024-03-04 [1] Bioconductor 3.18 (R 4.3.1)
# S4Vectors            * 0.40.2     2023-11-23 [1] Bioconductor 3.18 (R 4.3.1)
# scales                 1.3.0      2023-11-28 [1] CRAN (R 4.3.1)
# sessioninfo          * 1.2.2      2021-12-06 [1] CRAN (R 4.3.1)
# SparseArray            1.2.4      2024-02-11 [1] Bioconductor 3.18 (R 4.3.1)
# statmod                1.5.0      2023-01-06 [1] CRAN (R 4.3.1)
# stringi                1.8.3      2023-12-11 [1] CRAN (R 4.3.1)
# stringr              * 1.5.1      2023-11-14 [1] CRAN (R 4.3.1)
# SummarizedExperiment * 1.32.0     2023-10-24 [1] Bioconductor
# tibble               * 3.2.1      2023-03-20 [1] CRAN (R 4.3.1)
# tidyr                * 1.3.1      2024-01-24 [1] CRAN (R 4.3.1)
# tidyselect             1.2.1      2024-03-11 [1] CRAN (R 4.3.1)
# tidyverse            * 2.0.0      2023-02-22 [1] CRAN (R 4.3.1)
# timechange             0.3.0      2024-01-18 [1] CRAN (R 4.3.1)
# tzdb                   0.4.0      2023-05-12 [1] CRAN (R 4.3.1)
# utf8                   1.2.4      2023-10-22 [1] CRAN (R 4.3.1)
# variancePartition    * 1.32.5     2024-02-16 [1] Bioconductor 3.18 (R 4.3.1)
# vctrs                  0.6.5      2023-12-01 [1] CRAN (R 4.3.1)
# vroom                  1.6.5      2023-12-05 [1] CRAN (R 4.3.1)
# withr                  3.0.0      2024-01-16 [1] CRAN (R 4.3.1)
# XVector                0.42.0     2023-10-24 [1] Bioconductor
# zlibbioc               1.48.2     2024-03-13 [1] Bioconductor 3.18 (R 4.3.1)
# 
# [1] /opt/view/rlib/R/library
# [2] /software/hgi/installs/softpack/rstudio/.spack/linux-ubuntu22.04-x86_64_v3/gcc-11.4.0/r-4.3.1-bfwldrk76z6f52upk47zepliekn7ayqz/rlib/R/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────


