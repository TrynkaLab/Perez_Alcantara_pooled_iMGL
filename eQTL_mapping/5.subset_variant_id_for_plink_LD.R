# LD annotation by plink requires variant IDs
# subsetting BIM files to interesting (moderately significant) positions we got with tensorQTL
.libPaths(c('/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/',"/software/teamtrynka/ma23/R4.1/libs",
            "/software/R-4.1.0/lib/R/library"))
library(stringr)
library(tidyverse)

# Load the bim file
bim_file = readr::read_table("../../data/genotype/plink_genotypes/all_pools.genotype.MAF05.filtered.bim", 
                             col_names = c("chr", "variant_id", "gd", "pos", "a1", "a2"),
                             col_types = cols(
                               chr = col_character(), # there's an X
                               variant_id = col_character(),
                               gd = col_double(),
                               pos = col_double(),
                               a1 = col_character(),
                               a2 = col_character()
                             ))

# tensorQTL file - all results

tensorqtl_all = readr::read_csv("../../data/results/4.Inspect_eQTL_results/tensorQTL_variant_gene_60PCs.csv")

#########
#### filtering tensor variants #######
tensorqtl_pfiltered = tensorqtl_all %>%
  dplyr::filter(pval_beta < 0.20) #### filter threshold in pval_beta
reformated_bim = bim_file %>%
  mutate(tensor_variant_id = paste(chr,pos,a2,a1, sep = "_"), # tensor variants with respect to bim file are chr_pos_A2_A1 (weird)
         wrong_id = paste(chr,pos,a1,a2, sep = "_")) # to check
#########
#######

sum(tensorqtl_pfiltered$variant_id %in% reformated_bim$tensor_variant_id) == nrow(tensorqtl_pfiltered) # TRUE: all should be present
# if not, investigate because the bim file loaded for eQTL analysis might have not be the correct one
sum(tensorqtl_pfiltered$variant_id %in% reformated_bim$wrong_id) # none present
length(unique(tensorqtl_pfiltered$variant_id)) # 15378
length(unique(tensorqtl_all$variant_id)) # 53382
sum(reformated_bim$tensor_variant_id %in% unique(tensorqtl_pfiltered$variant_id)) # 15378 = 0 missing 
sum(reformated_bim$tensor_variant_id %in% unique(tensorqtl_all$variant_id)) # 6466076 - 53383 = many missing 
# missing variants in tensor QTL results probably not reported because of the windows tested around genes and other tensorQTL filters
# numbers of reported gene-variant pairs are different in each cluster-treatment result, but don't vary between PCs
# so probably just depends on genes tested
sum(reformated_bim$wrong_id %in% unique(tensorqtl_pfiltered$variant_id)) # 0 

reformated_filtered_bim = reformated_bim %>%
  dplyr::filter(tensor_variant_id %in% unique(tensorqtl_pfiltered$variant_id))
rm(reformated_bim)
gc()

# save variant ids from filtered bim file
write_lines(reformated_filtered_bim$variant_id,
            "../../data/genotype/plink_genotypes/tensor_variants_to_subset.txt")
