
library(tidyverse)
library(dplyr)
library(rlang)
library(sessioninfo)

#-------------------------------------------------------------------------------
#             8.1 Testing for genotype-proliferation association 
#-------------------------------------------------------------------------------
#  Code to fit linear models to each proliferation phenotype for each 
#  variant in a chunk.
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

outdir = paste("..", "..", "output_data", "08_GWAS_proliferation", "01_Association_testing", sep = "/")
plotdir = paste("..", "..", "plots", "08_GWAS_proliferation", "01_Association_testing", sep = "/")
logs_dir = paste("logs", "01_Association_testing", sep = "/")
err_dir = paste("errors", "01_Association_testing", sep = "/")
dir.create(outdir, recursive = T)
dir.create(plotdir, recursive = T)
dir.create(logs_dir, recursive = T)
dir.create(err_dir, recursive = T)


args = commandArgs(trailingOnly = TRUE)

if (length(args)<10) {
  stop("You need to detail 5 input and 5 output file paths.n", call. = FALSE)
} else if (length(args) == 10) {
  line_prop_changes_premac_vs_iPSC_genotype_input_path = args[1]
  line_prop_changes_old_vs_young_premac_genotype_input_path = args[2]
  line_prop_changes_microglia_vs_premac_IFN_genotype_input_path = args[3]
  line_prop_changes_microglia_vs_premac_LPS_genotype_input_path = args[4]
  line_prop_changes_microglia_vs_premac_untreated_genotype_input_path = args[5]

  ## Out files for each chunk
  GWAS_premac_vs_iPSC_out_path = args[6]
  GWAS_old_vs_young_premac_out_path = args[7]
  GWAS_microglia_vs_premac_IFN_out_path = args[8]
  GWAS_microglia_vs_premac_LPS_out_path = args[9]
  GWAS_microglia_vs_premac_untreated_out_path = args[10]
}

## Extract chunk
chunk = line_prop_changes_premac_vs_iPSC_genotype_input_path %>% strsplit(., ".Rdata") %>% .[[1]] %>% strsplit(., "_") %>% .[[1]] %>% .[length(.)]

# args[1] = "../../output_data/08_GWAS_proliferation/00_Data_exploration_processing/premac_vs_iPSC_data/line_prop_changes_premac_vs_iPSC_genotype_4231.Rdata"
# args[2] = "../../output_data/08_GWAS_proliferation/00_Data_exploration_processing/old_vs_young_premac_data/line_prop_changes_old_vs_young_premac_genotype_4231.Rdata"
# args[3] = "../../output_data/08_GWAS_proliferation/00_Data_exploration_processing/microglia_vs_premac_IFN_data/line_prop_changes_microglia_vs_premac_IFN_genotype_4231.Rdata"
# args[4] = "../../output_data/08_GWAS_proliferation/00_Data_exploration_processing/microglia_vs_premac_LPS_data/line_prop_changes_microglia_vs_premac_LPS_genotype_4231.Rdata"
# args[5] = "../../output_data/08_GWAS_proliferation/00_Data_exploration_processing/microglia_vs_premac_untreated_data/line_prop_changes_microglia_vs_premac_untreated_genotype_4231.Rdata"
# 
# args[6] = "../../output_data/08_GWAS_proliferation/01_Association_testing/premac_vs_iPSC_results/GWAS_premac_vs_iPSC_results_4231.Rdata"
# args[7] = "../../output_data/08_GWAS_proliferation/01_Association_testing/old_vs_young_premac_results/GWAS_old_vs_young_premac_results_4231.Rdata"
# args[8] = "../../output_data/08_GWAS_proliferation/01_Association_testing/microglia_vs_premac_IFN_results/GWAS_microglia_vs_premac_IFN_results_4231.Rdata"
# args[9] = "../../output_data/08_GWAS_proliferation/01_Association_testing/microglia_vs_premac_LPS_results/GWAS_microglia_vs_premac_LPS_results_4231.Rdata"
# args[10] = "../../output_data/08_GWAS_proliferation/01_Association_testing/microglia_vs_premac_untreated_results/GWAS_microglia_vs_premac_untreated_results_4231.Rdata"


## Fit lm for each variant
fit_prolif_lm <- function(prolif_phenotype){
  
  input_data <- eval(parse(text = paste0("line_prop_changes_", prolif_phenotype, "_genotype")))
  out_path <- eval(parse(text = paste0("GWAS_", prolif_phenotype, "_out_path")))
  results <- list()
  
  ## Expected num of samples and variants in input
  num_samples <- c("premac_vs_iPSC" = 228, 
                   "old_vs_young_premac" = 210,
                   "microglia_vs_premac_IFN" = 230,
                   "microglia_vs_premac_LPS" = 229,
                   "microglia_vs_premac_untreated" = 209)
  
  num_samples_pheno = num_samples[prolif_phenotype]
  if(chunk != 582){
    num_variants = 9999 + 14
  }else{
    num_variants = 1888 + 14
  }
  
  ## Check input data were correctly generated
  if(nrow(input_data) == num_samples_pheno & ncol(input_data) == num_variants){
    
    ## Fit lm model for each variant in chunk
    for(i in 1:(ncol(input_data) - 14) ){
      
      input_data$allele_dosage <- input_data[, colnames(input_data)[i + 14]] %>% unlist() %>% as.double()
    
      model <- mean_scaled_log_fraction ~ allele_dosage + sex + npool + PC1 + PC2 + PC3 + PC4 + PC5
      
      head(model.matrix(model, data = input_data))
      # (Intercept) allele_dosage sexMale npool       PC1        PC2         PC3         PC4        PC5
      # 1           1             2       0     1 0.1157790  0.0886567 -0.02629670 -0.00873555 -0.1311270
      # 2           1             2       0     1 0.0601301  0.0703484  0.00865988 -0.02587890 -0.0715146
      # 3           1             2       0     1 0.0757558 -0.0274609 -0.06032140  0.14498100 -0.0952725
      # 4           1             1       1     1 0.0244120  0.1056540 -0.01054090 -0.13276500 -0.0490484
      # 5           1             1       0     1 0.0977595  0.0318457 -0.04368480 -0.00747156 -0.1009210
      # 6           1             2       0     2 0.1111900  0.1342830 -0.03189210 -0.06852020 -0.0559980
      
      fit = lm(model, data = input_data)
      sum_res = summary(fit)
      results[[colnames(input_data)[i + 14]]] = list(
        coefficients = sum_res$coefficients,
        residuals =  sum_res$residuals
      )
    }
    
    save(results, file = out_path)
  }
  else{
    stop(paste0("Input data not properly generated for chunk ", chunk, "."))
  }
}


################################################################################
##                    iPSC -> young preMac proliferation 
################################################################################

load(line_prop_changes_premac_vs_iPSC_genotype_input_path, verbose = T) 
# Loading objects:
#   line_prop_changes_premac_vs_iPSC_genotype

prolif_phenotype <- "premac_vs_iPSC"
fit_prolif_lm(prolif_phenotype)


################################################################################
##                 young preMac -> old preMac proliferation
################################################################################

load(line_prop_changes_old_vs_young_premac_genotype_input_path, verbose = T)
# Loading objects:
#   line_prop_changes_old_vs_young_premac_genotype

prolif_phenotype <- "old_vs_young_premac"
fit_prolif_lm(prolif_phenotype)


################################################################################
##           preMac -> microglia (in IFN/LPS/untreated) proliferation
################################################################################

#############################   Microglia in IFN   #############################

load(line_prop_changes_microglia_vs_premac_IFN_genotype_input_path, verbose = T)
# Loading objects:
#   line_prop_changes_microglia_vs_premac_IFN_genotype

prolif_phenotype <- "microglia_vs_premac_IFN"
fit_prolif_lm(prolif_phenotype)

#############################   Microglia in LPS   #############################

load(line_prop_changes_microglia_vs_premac_LPS_genotype_input_path, verbose = T)
# Loading objects:
#   line_prop_changes_microglia_vs_premac_LPS_genotype

prolif_phenotype <- "microglia_vs_premac_LPS"
fit_prolif_lm(prolif_phenotype)

##########################   Microglia in untreated   ##########################

load(line_prop_changes_microglia_vs_premac_untreated_genotype_input_path, verbose = T)
# Loading objects:
#   line_prop_changes_microglia_vs_premac_untreated_genotype

prolif_phenotype <- "microglia_vs_premac_untreated"
fit_prolif_lm(prolif_phenotype)







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
# date     2025-02-10
# rstudio  2024.04.0+735 Chocolate Cosmos (server)
# pandoc   3.1.12.3 @ /opt/view/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package     * version date (UTC) lib source
# cli           3.6.2   2023-12-11 [1] CRAN (R 4.3.1)
# colorspace    2.1-0   2023-01-23 [1] CRAN (R 4.3.1)
# dplyr       * 1.1.4   2023-11-17 [1] CRAN (R 4.3.1)
# fansi         1.0.6   2023-12-08 [1] CRAN (R 4.3.1)
# forcats     * 1.0.0   2023-01-29 [1] CRAN (R 4.3.1)
# generics      0.1.3   2022-07-05 [1] CRAN (R 4.3.1)
# ggplot2     * 3.5.1   2024-04-23 [1] CRAN (R 4.3.1)
# glue          1.7.0   2024-01-09 [1] CRAN (R 4.3.1)
# gtable        0.3.4   2023-08-21 [1] CRAN (R 4.3.1)
# here          1.0.1   2020-12-13 [1] CRAN (R 4.3.1)
# hms           1.1.3   2023-03-21 [1] CRAN (R 4.3.1)
# jsonlite      1.8.8   2023-12-04 [1] CRAN (R 4.3.1)
# lattice       0.22-6  2024-03-20 [1] CRAN (R 4.3.1)
# lifecycle     1.0.4   2023-11-07 [1] CRAN (R 4.3.1)
# lubridate   * 1.9.3   2023-09-27 [1] CRAN (R 4.3.1)
# magrittr      2.0.3   2022-03-30 [1] CRAN (R 4.3.1)
# Matrix        1.6-5   2024-01-11 [1] CRAN (R 4.3.1)
# munsell       0.5.1   2024-04-01 [1] CRAN (R 4.3.1)
# pillar        1.9.0   2023-03-22 [1] CRAN (R 4.3.1)
# pkgconfig     2.0.3   2019-09-22 [1] CRAN (R 4.3.1)
# png           0.1-8   2022-11-29 [1] CRAN (R 4.3.1)
# purrr       * 1.0.2   2023-08-10 [1] CRAN (R 4.3.1)
# R6            2.5.1   2021-08-19 [1] CRAN (R 4.3.1)
# rappdirs      0.3.3   2021-01-31 [1] CRAN (R 4.3.1)
# Rcpp          1.0.12  2024-01-09 [1] CRAN (R 4.3.1)
# readr       * 2.1.5   2024-01-10 [1] CRAN (R 4.3.1)
# reticulate    1.35.0  2024-01-31 [1] CRAN (R 4.3.1)
# rlang         1.1.3   2024-01-10 [1] CRAN (R 4.3.1)
# rprojroot     2.0.4   2023-11-05 [1] CRAN (R 4.3.1)
# rstudioapi    0.16.0  2024-03-24 [1] CRAN (R 4.3.1)
# scales        1.3.0   2023-11-28 [1] CRAN (R 4.3.1)
# sessioninfo * 1.2.2   2021-12-06 [1] CRAN (R 4.3.1)
# stringi       1.8.3   2023-12-11 [1] CRAN (R 4.3.1)
# stringr     * 1.5.1   2023-11-14 [1] CRAN (R 4.3.1)
# tibble      * 3.2.1   2023-03-20 [1] CRAN (R 4.3.1)
# tidyr       * 1.3.1   2024-01-24 [1] CRAN (R 4.3.1)
# tidyselect    1.2.1   2024-03-11 [1] CRAN (R 4.3.1)
# tidyverse   * 2.0.0   2023-02-22 [1] CRAN (R 4.3.1)
# timechange    0.3.0   2024-01-18 [1] CRAN (R 4.3.1)
# tzdb          0.4.0   2023-05-12 [1] CRAN (R 4.3.1)
# utf8          1.2.4   2023-10-22 [1] CRAN (R 4.3.1)
# vctrs         0.6.5   2023-12-01 [1] CRAN (R 4.3.1)
# withr         3.0.0   2024-01-16 [1] CRAN (R 4.3.1)
# 
# [1] /opt/view/rlib/R/library
# [2] /software/hgi/installs/softpack/rstudio/.spack/linux-ubuntu22.04-x86_64_v3/gcc-11.4.0/r-4.3.1-bfwldrk76z6f52upk47zepliekn7ayqz/rlib/R/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────



