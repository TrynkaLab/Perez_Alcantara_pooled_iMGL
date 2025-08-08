# compare new and old genotype
.libPaths(c('/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/',
            "/software/teamtrynka/ma23/R4.1/libs",
            "/software/R-4.1.0/lib/R/library"
)) 
library(tidyverse)
library(patchwork)
options(stringsAsFactors=F)

new_genotype = readr::read_tsv("../../data/all_pools_no_curn_iukl/merged_phased_imputed_UK10K_1KG_HRC_hg38_Bellenguez_checked.txt") %>%
  tidyr::pivot_longer(cols =!c(CHROM,POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT),
                      names_to = "donor",values_to = "genotype") %>%
  dplyr::mutate(genotype = str_split(genotype, ":", simplify = TRUE)[,1]) %>%
  dplyr::mutate(allele1 = as.numeric(str_split_i(genotype, "\\|", i=1)),
                allele2=as.numeric(str_split_i(genotype, "\\|", i=2))) %>%
  dplyr::mutate(gtCount = allele1 + allele2) %>%
  dplyr::select(CHROM,POS,REF,ALT,INFO,donor,gtCount) %>%
  dplyr::mutate(type = "new")
old_genotype = readr::read_tsv("../old_hipsci_imputation/check_AD_old.vcf") %>%
  tidyr::pivot_longer(cols =!c(CHROM,POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT),
                      names_to = "donor",values_to = "genotype") %>%
  dplyr::mutate(genotype = str_split(genotype, ":", simplify = TRUE)[,1]) %>%
  dplyr::mutate(allele1 = as.numeric(str_split_i(genotype, "\\|", i=1)),
                allele2=as.numeric(str_split_i(genotype, "\\|", i=2))) %>%
  dplyr::mutate(gtCount = allele1 + allele2) %>%
  dplyr::select(CHROM,POS,REF,ALT,INFO,donor,gtCount) %>%
  dplyr::mutate(type = "old")
kaur_genotype = readr::read_tsv("/lustre/scratch123/hgi/projects/otar2065/resources/hiPSCi_Alasoo_reimputed_GRCh38_2021/data/check_AD_kaur.vcf") %>%
  tidyr::pivot_longer(cols =!c(CHROM,POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT),
                      names_to = "donor",values_to = "genotype") %>%
  dplyr::mutate(genotype = str_split(genotype, ":", simplify = TRUE)[,1]) %>%
  dplyr::mutate(allele1 = as.numeric(str_split_i(genotype, "\\|", i=1)),
                allele2=as.numeric(str_split_i(genotype, "\\|", i=2)),
                CHROM=paste0("chr",CHROM),
                donor = str_split_i(donor, "-", i=2)) %>%
  dplyr::mutate(gtCount = allele1 + allele2) %>%
  dplyr::select(CHROM,POS,REF,ALT,INFO,donor,gtCount) %>%
  dplyr::mutate(type = "kaur") %>%
  dplyr::filter(donor %in% unique(new_genotype$donor))
  
coxy_kaur = kaur_genotype %>%
  dplyr::select(c("CHROM","POS",contains("coxy")))



# I thiink kaur has extra genotype alts

HRC_new = readr::read_tsv("../../data/all_pools_no_curn_iukl/HRC/check_AD_HRC.vcf") %>%
  tidyr::pivot_longer(cols =!c(CHROM,POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT),
                      names_to = "donor",values_to = "genotype") %>%
  dplyr::mutate(genotype = str_split(genotype, ":", simplify = TRUE)[,1]) %>%
  dplyr::mutate(allele1 = as.numeric(str_split_i(genotype, "\\|", i=1)),
                allele2=as.numeric(str_split_i(genotype, "\\|", i=2))) %>%
  dplyr::mutate(gtCount = allele1 + allele2) %>%
  dplyr::select(CHROM,POS,REF,ALT,INFO,donor,gtCount) %>%
  dplyr::mutate(type = "HRC_new")
tenK_new = readr::read_tsv("../../data/all_pools_no_curn_iukl/1K_10K/check_AD_1K_10K.vcf") %>%
  tidyr::pivot_longer(cols =!c(CHROM,POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT),
                      names_to = "donor",values_to = "genotype") %>%
  dplyr::mutate(genotype = str_split(genotype, ":", simplify = TRUE)[,1]) %>%
  dplyr::mutate(allele1 = as.numeric(str_split_i(genotype, "\\|", i=1)),
                allele2=as.numeric(str_split_i(genotype, "\\|", i=2))) %>%
  dplyr::mutate(gtCount = allele1 + allele2) %>%
  dplyr::select(CHROM,POS,REF,ALT,INFO,donor,gtCount) %>%
  dplyr::mutate(type = "tenK_new")

merged_long = rbind(new_genotype,old_genotype,kaur_genotype,HRC_new,tenK_new) %>%
  dplyr::arrange(CHROM,POS,type)

merged_wide = merged_long %>%
  tidyr::pivot_wider(names_from = c("donor"),values_from = "gtCount") %>%
  dplyr::arrange(CHROM,POS)

# how many shared positions are identical
identical_tally = list()
# new vs old
identical_tally[["new_old"]] = merged_long %>%
  dplyr::filter(type %in% c("new","old")) %>%
  dplyr::group_by(CHROM,POS,REF,ALT,donor) %>%
  dplyr::filter(n()>1) %>%
  dplyr::summarise(are_gt_identical = all(gtCount == first(gtCount))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(CHROM,POS,REF,ALT) %>%
  dplyr::count(are_gt_identical) %>%
  dplyr::mutate(prop_identical = n/sum(n)) %>%
  dplyr::filter(are_gt_identical == TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(.,merged_long[merged_long$type=="new",c("CHROM",  "POS", "REF",   "ALT",   "INFO")]) %>% # adding new imputation info
  dplyr::distinct() %>%
  dplyr::mutate(comparison = "new_old")
  
# new vs kaur
identical_tally[["new_kaur"]] = merged_long %>%
  dplyr::filter(type %in% c("new","kaur")) %>%
  dplyr::group_by(CHROM,POS,REF,ALT,donor) %>%
  dplyr::filter(n()>1) %>%
  dplyr::summarise(are_gt_identical = all(gtCount == first(gtCount))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(CHROM,POS,REF,ALT) %>%
  dplyr::count(are_gt_identical) %>%
  dplyr::mutate(prop_identical = n/sum(n)) %>%
  dplyr::filter(are_gt_identical == TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(.,merged_long[merged_long$type=="new",c("CHROM",  "POS", "REF",   "ALT",   "INFO")]) %>% # adding new imputation info
  dplyr::distinct() %>%
  dplyr::mutate(comparison = "new_kaur")

# old vs kaur
identical_tally[["old_kaur"]] = merged_long %>%
  dplyr::filter(type %in% c("old","kaur")) %>%
  dplyr::group_by(CHROM,POS,REF,ALT,donor) %>%
  dplyr::filter(n()>1) %>%
  dplyr::summarise(are_gt_identical = all(gtCount == first(gtCount))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(CHROM,POS,REF,ALT) %>%
  dplyr::count(are_gt_identical) %>%
  dplyr::mutate(prop_identical = n/sum(n)) %>%
  dplyr::filter(are_gt_identical == TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(.,merged_long[merged_long$type=="new",c("CHROM",  "POS", "REF",   "ALT",   "INFO")]) %>% # adding new imputation info
  dplyr::distinct() %>%
  dplyr::mutate(comparison = "old_kaur")
# HRC new vs kaur
identical_tally[["HRC_new_kaur"]] =merged_long %>%
  dplyr::filter(type %in% c("HRC_new","kaur")) %>%
  dplyr::group_by(CHROM,POS,REF,ALT,donor) %>%
  dplyr::filter(n()>1) %>%
  dplyr::summarise(are_gt_identical = all(gtCount == first(gtCount))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(CHROM,POS,REF,ALT) %>%
  dplyr::count(are_gt_identical) %>%
  dplyr::mutate(prop_identical = n/sum(n)) %>%
  dplyr::filter(are_gt_identical == TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(.,merged_long[merged_long$type=="new",c("CHROM",  "POS", "REF",   "ALT",   "INFO")]) %>% # adding new imputation info
  dplyr::distinct() %>%
  dplyr::mutate(comparison = "HRC_new_kaur")


# tenk new vs kaur
identical_tally[["tenK_new_kaur"]] = merged_long %>%
  dplyr::filter(type %in% c("tenK_new","kaur")) %>%
  dplyr::group_by(CHROM,POS,REF,ALT,donor) %>%
  dplyr::filter(n()>1) %>%
  dplyr::summarise(are_gt_identical = all(gtCount == first(gtCount))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(CHROM,POS,REF,ALT) %>%
  dplyr::count(are_gt_identical) %>%
  dplyr::mutate(prop_identical = n/sum(n)) %>%
  dplyr::filter(are_gt_identical == TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(.,merged_long[merged_long$type=="new",c("CHROM",  "POS", "REF",   "ALT",   "INFO")]) %>% # adding new imputation info
  dplyr::distinct() %>%
  dplyr::mutate(comparison = "tenk_new_kaur")
identical_tally_wide = do.call("rbind",identical_tally) %>%
  tidyr::pivot_wider(names_from = c("comparison"),values_from = "prop_identical") %>%
  dplyr::arrange(CHROM,POS)

## prop of identical positions is correlated with imputation R2 within each dataset compared
summary(identical_tally_wide$new_old)
summary(identical_tally_wide$new_kaur)
summary(identical_tally_wide$old_kaur) # highest prop of identical positions
summary(identical_tally_wide$HRC_new_kaur)
summary(identical_tally_wide$tenk_new_kaur)
