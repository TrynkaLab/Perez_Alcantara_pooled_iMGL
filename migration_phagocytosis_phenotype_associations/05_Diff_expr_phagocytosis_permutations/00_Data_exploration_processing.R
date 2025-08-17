
library(tidyverse)
library(dplyr)
library(SummarizedExperiment)
library(rlang)
library(edgeR)
library(variancePartition)
library(sessioninfo)


################################################################################
##               5. Permutation analysis for phagocytosis DGE
################################################################################

#-------------------------------------------------------------------------------
#                     5.0 Data exploration and processing
#-------------------------------------------------------------------------------
#  Code to explore and process gene expression and cell line phenotypic 
#  metadata to run DGE on permuted phenotypes in the next script.
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

## Set working dir
setwd("/lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_Daianna/")

## Define output data and plot dir
outdir = paste(getwd(), "output_data", "05_Diff_expr_phagocytosis_permutations", "00_Data_exploration_processing", sep = "/")
plotdir = paste(getwd(), "plots", "05_Diff_expr_phagocytosis_permutations", "00_Data_exploration_processing", sep = "/")
dir.create(outdir, recursive = T)
dir.create(plotdir, recursive = T)

## Input dirs
input_dir0 = paste(getwd(), "output_data", "01_Diff_expr_PRS_permutations", "00_Data_exploration_processing", sep = "/")

## color dict
load(paste0(input_dir0, "/var_colors.Rdata"), verbose = T)

## Load rse
load(paste0(input_dir0, "/rse_pseudobulk.Rdata"), verbose = T)


######################################### 
##  Add genotype PCs to sample metadata
######################################### 

## PCs
genotype_pcs = read.table("../OTAR2065_sc_eQTL/data/genotype/plink_genotypes/all_pools.all_donors.genotype.MAF05.eigenvec")[,-1] 
rownames(genotype_pcs) <- genotype_pcs[,1]
colnames(genotype_pcs) <- c("line", paste0("genotype_PC", 1:(ncol(genotype_pcs)-1)))

## Add first 5 PCs to sample metadata
five_genotype_pcs = genotype_pcs[,c(1:6)]

## All lines in sample metadata with PCs 
length(which(! unique(rse$line) %in% five_genotype_pcs$line))
# [1] 0

colData(rse)[, c(colnames(five_genotype_pcs))] <- as.data.frame(colData(rse)) %>% 
                                                    left_join(five_genotype_pcs, by = "line") %>% 
                                                    select(colnames(five_genotype_pcs))


##############################################
##  Add phagocytosis data to sample metadata
##############################################

line_prop_changes = read_csv("../OTAR2065_phenotypic_QTLs/data/results/phagocytosis/1.check_line_proportions/line_prop_changes_per_well.csv") %>% as.data.frame()
dim(line_prop_changes)
# [1] 1821   14

head(line_prop_changes, 3)
#         line  well log_fraction_mean log_fraction_median log_fraction_max prop_unadjusted_max_value prop_unadjusted_min_value   pool
# 1 Arlene-003 15A-2        -0.5259682          -0.5277118       -0.6243050                0.01064571               0.007093153 pool15
# 2 Arlene-003 15A-3        -0.4582498          -0.4594732       -0.5403332                0.01047626               0.007362023 pool15
# 3 Arlene-003 15B-2        -0.4851264          -0.4806216       -0.6551755                0.01109790               0.007571185 pool15
#      condition treatment replicate  donor scaled_log_fraction    sex
# 1 Cargo_Well_A untreated      rep1 Arlene          0.16688201 Female
# 2 Cargo_Well_B untreated      rep2 Arlene          0.25865586 Female
# 3 Cargo_Well_A       IFN      rep1 Arlene          0.09037676 Female

## Data per line x pool x treatment x replica (well)
dim(unique(line_prop_changes[, c("line", "pool", "treatment", "replicate")]))
# [1] 1821    4
dim(unique(line_prop_changes[, c("line", "pool", "treatment", "well")]))
# [1] 1821    4

## Num of unique lines
length(unique(line_prop_changes$line))
# [1] 202

## Num of pools
length(unique(line_prop_changes$pool))
# [1] 14

## Num of treatments
length(unique(line_prop_changes$treatment))
# [1] 3

## Num of wells
length(unique(line_prop_changes$well))
# [1] 93

## Num of conditions
unique(line_prop_changes$condition)
# [1] "Cargo_Well_A" "Cargo_Well_B" "Cargo_Well_C"

## Num of replicates per line x pool x treatment
line_prop_changes$line_pool_treatment <- paste(line_prop_changes$line, line_prop_changes$pool, line_prop_changes$treatment)
table(line_prop_changes$line_pool_treatment) %>% table
#  1   2   3 
# 72 627 165 

## Num of replicates per sample per treatment 
by(line_prop_changes$line_pool_treatment, line_prop_changes$treatment, function(x){table(table(x))})

# line_prop_changes$treatment: IFN
# 1   2   3 
# 18 212  65 
# ------------------------------------------------------------------------------------------------------- 
#   line_prop_changes$treatment: LPS
# 
# 1   2   3 
# 35 189  64 
# ------------------------------------------------------------------------------------------------------- 
#   line_prop_changes$treatment: untreated
# 
# 1   2   3 
# 19 226  36 

## Check all scaled and unscaled log-fractions are valid
which(is.na(line_prop_changes$log_fraction_mean))
# integer(0)
which(line_prop_changes$log_fraction_mean %in% c(NA,Inf,-Inf))
# integer(0)
which(is.na(line_prop_changes$scaled_log_fraction))
# integer(0)
which(line_prop_changes$scaled_log_fraction %in% c(NA,Inf,-Inf))
# integer(0)

## Subset to samples with min prop > 0.005
line_prop_changes_filt <- subset(line_prop_changes, prop_unadjusted_min_value > 0.005)
dim(line_prop_changes_filt)
# [1] 1267   15

## Lines left
length(unique(line_prop_changes_filt$line))
# [1] 159

## Mean scaled log-fractions and min proportions across replicates
mean_line_phenotype_across_replicates <-  line_prop_changes_filt %>%  dplyr::group_by(line, pool, treatment, line_pool_treatment) %>% 
                                                                      dplyr::summarise(scaled_log_fraction_mean_across_replicates = mean(scaled_log_fraction), 
                                                                                       prop_unadjusted_min_value_mean_across_replicates = mean(prop_unadjusted_min_value)) %>% 
                                                                      as.data.frame()
                                                                    
## % of samples in rse with phago data
rse$line_pool_treatment <- paste(rse$line, rse$pool, rse$treatment)
length(which(rse$line_pool_treatment %in% mean_line_phenotype_across_replicates$line_pool_treatment)) / length(rse$line_pool_treatment) * 100
# [1] 50.93885

## Num of lines in rse with phago data
length(which(unique(rse$line) %in% mean_line_phenotype_across_replicates$line))
# 158 / 250

## Add line phagocytosis estimates to sample metadata
colData(rse)$phagocytosis_mean_scaled_log_fraction <- as.data.frame(colData(rse)) %>% 
                                                  left_join(mean_line_phenotype_across_replicates, by = c("line_pool_treatment")) %>% 
                                                  select(scaled_log_fraction_mean_across_replicates) %>% unname() %>% unlist()

colData(rse)$phagocytosis_mean_min_prop <- as.data.frame(colData(rse)) %>% 
  left_join(mean_line_phenotype_across_replicates, by = c("line_pool_treatment")) %>% 
  select(prop_unadjusted_min_value_mean_across_replicates) %>% unname() %>% unlist()


######################################
##  Final rse for phagocytosis DGE 
###################################### 

## All pseudobulk samples (prolif x treatment x pool x cell line)
dim(rse)
# [1] 18648  2077

## Subset to non-prolif cluster
rse_non_prolif <- rse[, which(rse$proliferation_status == "Not_proliferating")]
dim(rse_non_prolif)
# [1] 18648   1244

## Unique cell lines
length(unique(rse_non_prolif$line))
# [1] 250

## Discard samples with <100 cells
rse_non_prolif_filt <- rse_non_prolif[, rse_non_prolif$ncells >=100]
dim(rse_non_prolif_filt)
# [1] 18648   781

length(unique(rse_non_prolif_filt$line))
# [1] 194

## Final sample set (samples with PCs and phago data)
rse_non_prolif_filt_phago <- rse_non_prolif_filt[, which(!is.na(rse_non_prolif_filt$phagocytosis_mean_scaled_log_fraction))]
dim(rse_non_prolif_filt_phago)
# [1] 18648   537

length(unique(rse_non_prolif_filt_phago$line))
# [1] 154

## Check data are complete
which(is.na(rse_non_prolif_filt_phago$phagocytosis_mean_min_prop))
# integer(0)
apply(colData(rse_non_prolif_filt_phago)[, paste0("genotype_PC", 1:5)], 2, function(x){which(!is.double(x))})
# integer(0)


## Num of replicates per sample per treatment in rse
line_prop_changes_in_rse <- subset(line_prop_changes_filt, line_pool_treatment %in% rse_non_prolif_filt_phago$line_pool_treatment)
by(line_prop_changes_in_rse$line_pool_treatment, line_prop_changes_in_rse$treatment, function(x){table(table(x))})

# line_prop_changes_in_rse$treatment: IFN
# 
#  1   2   3 
# 12 130  38 
# ------------------------------------------------------------------------------------------------------- 
#   line_prop_changes_in_rse$treatment: LPS
# 
#  1   2   3 
# 21 122  39 
# ------------------------------------------------------------------------------------------------------- 
#   line_prop_changes_in_rse$treatment: untreated
# 
# 1   2   3 
# 6 144  25 


## Add 100 permutations of phagocytosis scaled log-fractions for each treatment
set.seed(12262024)
rse_non_prolif_filt_phago_IFN <- rse_non_prolif_filt_phago[, which(rse_non_prolif_filt_phago$treatment == "IFN")]
colData(rse_non_prolif_filt_phago_IFN)[, paste0("permutation_phago_mean_scaled_log_fraction_", 1:100)] <- replicate(n = 100, 
                                                                                 expr = sample(rse_non_prolif_filt_phago_IFN$phagocytosis_mean_scaled_log_fraction, replace = FALSE))

rse_non_prolif_filt_phago_LPS <- rse_non_prolif_filt_phago[, which(rse_non_prolif_filt_phago$treatment == "LPS")]
colData(rse_non_prolif_filt_phago_LPS)[, paste0("permutation_phago_mean_scaled_log_fraction_", 1:100)] <- replicate(n = 100, 
                                                                                                             expr = sample(rse_non_prolif_filt_phago_LPS$phagocytosis_mean_scaled_log_fraction, replace = FALSE))

rse_non_prolif_filt_phago_untreated <- rse_non_prolif_filt_phago[, which(rse_non_prolif_filt_phago$treatment == "untreated")]
colData(rse_non_prolif_filt_phago_untreated)[, paste0("permutation_phago_mean_scaled_log_fraction_", 1:100)] <- replicate(n = 100, 
                                                                                                             expr = sample(rse_non_prolif_filt_phago_untreated$phagocytosis_mean_scaled_log_fraction, replace = FALSE))

## Explore lines per treatment 
by(rse_non_prolif_filt_phago$line, rse_non_prolif_filt_phago$treatment, function(x){length(unique(x))})
# rse_non_prolif_filt_phago$treatment: IFN
# [1] 139
# ------------------------------------------------------------------------------------------------------- 
#   rse_non_prolif_filt_phago$treatment: LPS
# [1] 146
# ------------------------------------------------------------------------------------------------------- 
#   rse_non_prolif_filt_phago$treatment: untreated
# [1] 140


save(rse, file = paste0(outdir, "/rse_complete_with_phagocitosis_data.Rdata"))
save(rse_non_prolif_filt_phago, file = paste0(outdir, "/rse_non_prolif_filt_phago.Rdata"))
save(rse_non_prolif_filt_phago_IFN, file = paste0(outdir, "/rse_non_prolif_filt_phago_IFN.Rdata"))
dim(rse_non_prolif_filt_phago_IFN)
# [1] 18648   180
save(rse_non_prolif_filt_phago_LPS, file = paste0(outdir, "/rse_non_prolif_filt_phago_LPS.Rdata"))
dim(rse_non_prolif_filt_phago_LPS)
# [1] 18648   182
save(rse_non_prolif_filt_phago_untreated, file = paste0(outdir, "/rse_non_prolif_filt_phago_untreated.Rdata"))
dim(rse_non_prolif_filt_phago_untreated)
# [1] 18648   175


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
##       5.0.1 Gene expression variance partition analysis per treatment 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

for(treatment in c("IFN", "LPS", "untreated")){
  
  rse_treatment <- eval(parse_expr(paste0("rse_non_prolif_filt_phago_", treatment)))
  
  ## Filter lowly-expressed genes
  expressed_genes <- rowSums(edgeR::cpm(rse_treatment) > 0.1) >= floor(0.3*ncol(rse_treatment))
  message(paste0(sum(expressed_genes), " genes out of ", nrow(rse), " with >0.1 CPM in at least 30% of samples"))
  rse_expr = rse_treatment[expressed_genes, ]  
  
  ## Log-normalize counts with TMM 
  assays(rse_expr)$logcounts <- edgeR::cpm(calcNormFactors(rse_expr, method = "TMM"), 
                                                log = TRUE, prior.count = 0.5)
  ## Covariates
  formula <- ~ phagocytosis_mean_scaled_log_fraction + phagocytosis_mean_min_prop + genotype_PC1 + genotype_PC2 + (1|Sex) + (1|pool)
  
  ## Discard genes with var = 0
  genes_var_zero <- which(apply(assays(rse_expr)$logcounts, 1, var)==0)
  if (length(genes_var_zero)>0){
    rse_expr <- rse_expr[-genes_var_zero, ]
  }
  
  ## Fit model per expressed gene and extract FVE for each covariate
  varPart<- fitExtractVarPartModel(assays(rse_expr)$logcounts, formula, colData(rse_expr))
  dim(varPart)
  # [1] 17413     7
  
  # Sort variables by median FVE
  vp <- sortCols(varPart)
  p <- plotVarPart(vp, col = var_colors) + theme(axis.text.x = element_text(size = 9))
  ggsave(filename= paste0(plotdir, "/VarPart_", treatment,".pdf"),
         p, width = 23, height = 12, units = "cm")
  
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
# date     2024-12-18
# rstudio  2024.04.0+735 Chocolate Cosmos (server)
# pandoc   3.1.12.3 @ /opt/view/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version    date (UTC) lib source
# abind                  1.4-5      2016-07-21 [1] CRAN (R 4.3.1)
# backports              1.4.1      2021-12-13 [1] CRAN (R 4.3.1)
# Biobase              * 2.62.0     2023-10-24 [1] Bioconductor
# BiocGenerics         * 0.48.1     2023-11-01 [1] Bioconductor
# bit                    4.0.5      2022-11-15 [1] CRAN (R 4.3.1)
# bit64                  4.0.5      2020-08-30 [1] CRAN (R 4.3.1)
# bitops                 1.0-7      2021-04-24 [1] CRAN (R 4.3.1)
# boot                   1.3-30     2024-02-26 [1] CRAN (R 4.3.1)
# broom                  1.0.5      2023-06-09 [1] CRAN (R 4.3.1)
# car                    3.1-2      2023-03-30 [1] CRAN (R 4.3.1)
# carData                3.0-5      2022-01-06 [1] CRAN (R 4.3.1)
# circlize               0.4.16     2024-02-20 [1] CRAN (R 4.3.1)
# cli                    3.6.2      2023-12-11 [1] CRAN (R 4.3.1)
# clue                   0.3-65     2023-09-23 [1] CRAN (R 4.3.1)
# cluster                2.1.6      2023-12-01 [1] CRAN (R 4.3.1)
# codetools              0.2-20     2024-03-31 [1] CRAN (R 4.3.1)
# colorspace             2.1-0      2023-01-23 [1] CRAN (R 4.3.1)
# ComplexHeatmap         2.18.0     2023-10-24 [1] Bioconductor
# crayon                 1.5.2      2022-09-29 [1] CRAN (R 4.3.1)
# DelayedArray           0.28.0     2023-10-24 [1] Bioconductor
# digest                 0.6.35     2024-03-11 [1] CRAN (R 4.3.1)
# doParallel             1.0.17     2022-02-07 [1] CRAN (R 4.3.1)
# dplyr                * 1.1.4      2023-11-17 [1] CRAN (R 4.3.1)
# edgeR                * 4.0.16     2024-02-18 [1] Bioconductor 3.18 (R 4.3.1)
# fansi                  1.0.6      2023-12-08 [1] CRAN (R 4.3.1)
# forcats              * 1.0.0      2023-01-29 [1] CRAN (R 4.3.1)
# foreach                1.5.2      2022-02-02 [1] CRAN (R 4.3.1)
# generics               0.1.3      2022-07-05 [1] CRAN (R 4.3.1)
# GenomeInfoDb         * 1.38.8     2024-03-15 [1] Bioconductor 3.18 (R 4.3.1)
# GenomeInfoDbData       1.2.11     2024-11-06 [1] Bioconductor
# GenomicRanges        * 1.54.1     2023-10-29 [1] Bioconductor
# GetoptLong             1.0.5      2020-12-15 [1] CRAN (R 4.3.1)
# ggplot2              * 3.5.1      2024-04-23 [1] CRAN (R 4.3.1)
# ggpubr                 0.6.0      2023-02-10 [1] CRAN (R 4.3.1)
# ggsignif               0.6.4      2022-10-13 [1] CRAN (R 4.3.1)
# GlobalOptions          0.1.2      2020-06-10 [1] CRAN (R 4.3.1)
# glue                   1.7.0      2024-01-09 [1] CRAN (R 4.3.1)
# gtable                 0.3.4      2023-08-21 [1] CRAN (R 4.3.1)
# hms                    1.1.3      2023-03-21 [1] CRAN (R 4.3.1)
# IRanges              * 2.36.0     2023-10-24 [1] Bioconductor
# iterators              1.0.14     2022-02-05 [1] CRAN (R 4.3.1)
# jtools                 2.2.2      2023-07-11 [1] CRAN (R 4.3.1)
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
# nlme                   3.1-164    2023-11-27 [1] CRAN (R 4.3.1)
# nloptr                 2.0.3      2022-05-26 [1] CRAN (R 4.3.1)
# numDeriv               2016.8-1.1 2019-06-06 [1] CRAN (R 4.3.1)
# pander                 0.6.5      2022-03-18 [1] CRAN (R 4.3.1)
# pillar                 1.9.0      2023-03-22 [1] CRAN (R 4.3.1)
# pkgconfig              2.0.3      2019-09-22 [1] CRAN (R 4.3.1)
# pkgload                1.3.4      2024-01-16 [1] CRAN (R 4.3.1)
# png                    0.1-8      2022-11-29 [1] CRAN (R 4.3.1)
# purrr                * 1.0.2      2023-08-10 [1] CRAN (R 4.3.1)
# R6                     2.5.1      2021-08-19 [1] CRAN (R 4.3.1)
# RColorBrewer           1.1-3      2022-04-03 [1] CRAN (R 4.3.1)
# Rcpp                   1.0.12     2024-01-09 [1] CRAN (R 4.3.1)
# RCurl                  1.98-1.14  2024-01-09 [1] CRAN (R 4.3.1)
# readr                * 2.1.5      2024-01-10 [1] CRAN (R 4.3.1)
# rjson                  0.2.21     2022-01-09 [1] CRAN (R 4.3.1)
# rlang                  1.1.3      2024-01-10 [1] CRAN (R 4.3.1)
# rstatix                0.7.2      2023-02-01 [1] CRAN (R 4.3.1)
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
# tibble               * 3.2.1      2023-03-20 [1] CRAN (R 4.3.1)
# tidyr                * 1.3.1      2024-01-24 [1] CRAN (R 4.3.1)
# tidyselect             1.2.1      2024-03-11 [1] CRAN (R 4.3.1)
# tidyverse            * 2.0.0      2023-02-22 [1] CRAN (R 4.3.1)
# timechange             0.3.0      2024-01-18 [1] CRAN (R 4.3.1)
# tzdb                   0.4.0      2023-05-12 [1] CRAN (R 4.3.1)
# utf8                   1.2.4      2023-10-22 [1] CRAN (R 4.3.1)
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


