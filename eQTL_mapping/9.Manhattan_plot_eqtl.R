# Manhattan plot eQTL
.libPaths(c("/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/",
            "/software/teamtrynka/conda/otar2065/lib/R/library",
            "/software/teamtrynka/ma23/R4.1/libs"))
library(tidyverse)
source("./functions.R")

output_path = "../../data/results/4.Inspect_eQTL_results/"

eQTL_nominal_untreated=readr::read_tsv("../../data/results/tensorqtl/60/sum_sizefactorsNorm_log2_scaled_centered_Not_proliferating_untreated_common_500kb_window_tensorQTL_nominal.txt")
# eQTL_nominal_LPS=readr::read_tsv("../../data/results/tensorqtl/60/")
# eQTL_nominal_IFN=readr::read_tsv("../../data/results/tensorqtl/60/")
eQTL_sign = readr::read_csv("../../../OTAR2065_sc_eQTL/data/results/4.Inspect_eQTL_results/tensorQTL_variant_gene_60PCs.csv")
eQTL_sign = eQTL_sign %>%
  dplyr::filter(group == "60_Not_proliferating_untreated" & qval < 0.05)

test = eQTL_nominal_untreated %>%
  dplyr::mutate(CHR=as.numeric(stringr::str_split_fixed(string = variant_id,pattern = "_",n = 4)[,1]),
                BP=as.numeric(stringr::str_split_fixed(string = variant_id,pattern = "_",n = 4)[,2])) %>%
  dplyr::rename(P=pval_nominal,SNP=variant_id) %>%
  dplyr::group_by(SNP) %>%
  # select smallest p-value for each SNP
  dplyr::filter(P == min(P)) %>%
  dplyr::ungroup() %>%
  distinct()
# subset highest significance
test_sign = test %>%
  dplyr::filter(SNP %in% unique(eQTL_sign$variant_id)) %>%
  dplyr::filter(P < 5e-8) # something weird is happening - probably nominal and annotated were not generated at asame time
# I shouldn't need extra pval filters, pvalues don't exactly match
# test_sign2 = test %>%
#   dplyr::filter(P < 5e-38)
# subsample those non-significant variants so the plot doesn't take ages
 
test_notsign = test %>%
  dplyr::filter(!(SNP %in% unique(eQTL_sign$variant_id))) %>%
  dplyr::sample_n(500000) 

test_subset = rbind(test_sign,test_notsign,test_sign2)
 myManhattan( test_subset, graph.title = "eQTL: untreated", 
              suggestiveline = FALSE,
                genomewideline = FALSE, 
              even.facet = T,
              highlight = unique(test_sign$SNP),
              highlight.col = "red",
                 chrom.lab = c(as.character(1:22))) 
# usual saving doesn't work, make sure it's proper ggplot
# png(paste0(output_path,"eqtl_manhattan_untreated.png"),
#     width = 15, height = 7)
# plot(p)
# dev.off()
