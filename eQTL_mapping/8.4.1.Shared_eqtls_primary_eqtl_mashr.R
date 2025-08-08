# checking shared/specific eQTLs with mashr
# after loading module containing R libraries with softpack
# module load HGI/softpack/groups/otar2065/otar2065_test2/1
# based on https://cran.r-project.org/web/packages/mashr/vignettes/intro_mash.html
# https://cran.r-project.org/web/packages/mashr/vignettes/eQTL_outline.html
# and https://cran.r-project.org/web/packages/mashr/vignettes/mash_sampling.html
library(tidyverse)
library(mashr)
source("./functions.R")

set.seed(123)
outdir = "../../data/results/8.4.1.Shared_eqtls_primary_eqtl_mashr"
treatment_cols =  c(untreated = "#8D918B", IFN = "#3A5683", LPS = "#F8766D")

####### do as 8.4.Shared__eqtls_mashr.R but with the external primary eQTLs from macromap
