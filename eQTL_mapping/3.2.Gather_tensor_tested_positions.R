# Gathering all unique positions tested to then subset VCF for eQTL inspection
# this saves time in 4.Process_eQTL_results.R
library(tidyverse)

# read tensorQTL files

tensorqtl = list()
tensorqtl[[1]] = readr::read_delim("../../data/results/tensorqtl/63/sum_sizefactorsNorm_log2_scaled_centered_untreated_Not_proliferating_common_500kb_window_tensorQTL_nominal.txt")
tensorqtl[[2]] = readr::read_delim("../../data/results/tensorqtl/73/sum_sizefactorsNorm_log2_scaled_centered_IFN_Not_proliferating_common_500kb_window_tensorQTL_nominal.txt")
tensorqtl[[3]] = readr::read_delim("../../data/results/tensorqtl/84/sum_sizefactorsNorm_log2_scaled_centered_LPS_Not_proliferating_common_500kb_window_tensorQTL_nominal.txt")

tensorqtl = do.call("rbind",tensorqtl)
# gather unique positions

tensorqtl %>%
  dplyr::distinct(variant_id) %>%
  tidyr::separate_wider_delim(cols = variant_id,names =  c("CHR", "POS", "REF","ALT"), delim = "_") %>%
  dplyr::select(c("CHR", "POS")) %>%
  dplyr::mutate(CHR = as.numeric(CHR),POS = as.numeric(POS)) %>%
  dplyr::arrange(CHR,POS) %>%
  write.table(.,"../../data/genotype/tensorQTL_vars_tested.txt",sep = "\t",
              col.names = FALSE, row.names = FALSE, quote = FALSE)


