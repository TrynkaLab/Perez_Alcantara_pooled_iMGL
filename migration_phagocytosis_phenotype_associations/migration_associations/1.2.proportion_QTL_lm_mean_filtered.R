# line proportion QTL pipeline start

library(patchwork)
library(tidyverse)
library(lme4) # mixed models
library(lmerTest) # p-values of lme4
#library("emmeans") # emmeans is needed for follow-up tests 
# library("multcomp") # for advanced control for multiple testing/Type 1 errors.
library(ggbeeswarm)
library(ggpubr)
library(corrplot)
library(psych)
# library(ComplexHeatmap)
# library(cowplot)
# library(gridExtra)
# library(biglm)
# library(tictoc)
# library(foreach)
# library(doParallel)
# library(qvalue)
# library(microbenchmark)
source("../functions.R")

args = commandArgs(trailingOnly=TRUE)


if (length(args)<5) {
  stop("You need to detail 4 input and 1 output file paths.n",
       call. = FALSE)
} else if (length(args) == 5) {
  sample_info_path = args[1]
  line_prop_changes_path = args[2]
  genotype_path = args[3]
  genotype_pc_path=args[4]
  output_csv=args[5]
}

#output_csv = "../../../data/results/migration/1.2.proportion_QTL_lm_mean_filtered/prop_QTL_1.csv"
output_path = dirname(output_csv)
dir.create(output_path,recursive = T)

# to test
# sample_info_path = "../../../data/all_pools_migration_sample_info.csv"
# line_prop_changes_path = "../../../data/results/migration/1.check_line_proportions/line_prop_changes_averages_variances_1pct_filtered.csv"
# genotype_path = "../../../data/genotypes/full_genotype/genotype_minor_allele_dosage_1.csv"
# genotype_pc_path = "../../../../OTAR2065_sc_eQTL/data/genotype/plink_genotypes/all_pools.no_outliers.genotype.MAF05.eigenvec"

sample_info = read.csv(sample_info_path)
line_prop_changes = read.csv(line_prop_changes_path)

message("There are ", 
length(unique(line_prop_changes$line)), " lines in the analysis") # 158
# # load minor allele dosage file for all lines
# # 
genotype = readr::read_csv(genotype_path) %>%
  dplyr::relocate(rn)

### main ###

# ggplot(line_prop_changes,aes(x=log(mean_max_prop_rep),y = log(var_outcome_rep))) + 
#   geom_point() + theme_minimal()
# negative relationship - lines with larger proportions have smaller variances in their outcomes
# we already new this

# treat lines across pools as replicates, do not average across pools
line_prop_changes = line_prop_changes %>%
  dplyr::select(line,treatment,mean_outcome_pool,n_pools,sex) %>%
  distinct() 

message("There are ", length(unique(line_prop_changes$line)), " lines in the analysis after Na and Inf filters")

summary(line_prop_changes$mean_outcome_pool)
hist(line_prop_changes$mean_outcome_pool)

genotype = genotype %>% 
  dplyr::select(colnames(genotype)[colnames(genotype) %in% c("rn",unique(line_prop_changes$line))]) 
length(colnames(genotype)) - 1
message("There are ", length(colnames(genotype)) - 1, # 150
        " lines in the analysis after matching with the genotype file (removing clones).")

###### there are some missing lines - the individuals we don't want to check because
# their ancestry is uncertain

# remove genotypes that don't pass filters
# more than 20% of samples in different bins ( for 53 lines that's 11)
# to_subset = genotype %>% 
#   dplyr::select(!rn) %>% 
#   dplyr::mutate(hom_0 = rowSums(. == 0), hom_1 = rowSums(. == 1), het = rowSums(. == 0.5), 
#                 threshold = ceiling(0.2 * ncol(.)),rn = genotype$rn) %>%
#   dplyr::filter(hom_0 > threshold & hom_1 > threshold & het > threshold) %>%
#   .$rn

genotype_subset = genotype %>%
  dplyr::ungroup() %>%
  #dplyr::filter(rn %in% to_subset) %>%
  tidyr::pivot_longer(cols = !rn, values_to = "alt_allele_dosage", names_to = "line") %>%
  dplyr::mutate(alt_allele_dosage = dplyr::case_when(
    alt_allele_dosage == 0.5 ~ 1,
    alt_allele_dosage == 1   ~ 2,
    TRUE          ~ alt_allele_dosage  # Recode so var of allele dosage is 1, Keep 0 unchanged
  ))

length(unique(genotype_subset$rn)) # number of SNPs after filtering

# remove genotypes under 5% MAF
# doing the filtering again here because some lines have been lost after the NA filtering
# probably should do the filtering of the NAs in the pheno file before and then do the QC again with plink?

snps_to_keep = genotype_subset %>%
  dplyr::mutate(alt_allele_dosage = as.factor(alt_allele_dosage)) %>%
  count(rn,alt_allele_dosage) %>%
  dplyr::group_by(rn) %>%
  dplyr::filter(min(n) >= (0.05 * length(unique(genotype_subset$line)))) %>%
  ungroup() %>%
  dplyr::select(rn) %>%
  distinct() %>%
  unlist()


genotype_subset = genotype_subset %>%
  dplyr::filter(rn %in% snps_to_keep) %>%
  dplyr::mutate(rn = str_replace(rn, "chr", "")) 

length(unique(genotype_subset$rn)) # number of SNPs after filtering ~ 400


# tensorqtl swaps alleles so that the minor allele is always the less abundant one - checking in eQTL results
# because many are swapped compared to the original VCF
# this might take a while 
# 
# eqtl = eqtl %>%
#   dplyr::mutate(swapped_tensor_alleles = map_chr(variant_id, .f = swap_REF_ALT_alelle_names))
# 
# eqtl = eqtl %>%
#   dplyr::filter(variant_id %in% unique(genotype_subset$rn) | swapped_tensor_alleles %in% unique(genotype_subset$rn))
# 
# if(nrow(eqtl)!=0){
#   
#   
#   genotype_subset_eqtl = genotype_subset %>%
#     dplyr::filter(rn %in% eqtl$variant_id | rn %in% eqtl$swapped_tensor_alleles)
# } else {
#   message("There are no rows in the eqtl file - saving empty dataframe")
#   res_df = data.frame(coef_untreated=NA,p_untreated=NA,se_untreated = NA,
#                       coef_IFN = NA,p_IFN = NA,se_IFN = NA,
#                       coef_LPS = NA,p_LPS = NA,se_LPS = NA,
#                       snp=NA)
#   write.csv(res_df,output_csv,quote = FALSE,row.names = FALSE)
# }

# not subsetting - dirty fix to check
genotype_subset_eqtl = genotype_subset
if(exists("genotype_subset_eqtl")) {
  if (nrow(genotype_subset_eqtl) != 0) {
    length(unique(genotype_subset_eqtl$rn)) # number of SNPs after filtering
    # if subset to eqtl: we've reduced set from 10k to 3.5k approx, so to 35%
    # if not,
    
    ## add allele dosage, also genotype PCs
    message("Adding genotype PCs")
    message("Genotype pc path ", genotype_pc_path)
    
    full_genotype_pcs = read.table(genotype_pc_path)
    rownames(full_genotype_pcs) = ifelse(
      full_genotype_pcs$V1 == full_genotype_pcs$V2,
      yes = full_genotype_pcs$V1,
      no = paste(full_genotype_pcs$V1, full_genotype_pcs$V2, sep = "_")
    )
    full_genotype_pcs = full_genotype_pcs[c(2:7)] %>%
      dplyr::rename(
        line = V2,
        genotypePC1 = V3,
        genotypePC2 = V4,
        genotypePC3 = V5,
        genotypePC4 = V6,
        genotypePC5 = V7
      )
    
    line_prop_changes_IFN = line_prop_changes %>%
      dplyr::filter(treatment == "IFN") %>%
      dplyr::inner_join(., genotype_subset_eqtl) %>%
      dplyr::distinct() %>%
      dplyr::inner_join(., full_genotype_pcs) %>%
      dplyr::distinct()
    
    line_prop_changes_untreated = line_prop_changes %>%
      dplyr::filter(treatment == "untreated") %>%
      dplyr::inner_join(., genotype_subset_eqtl) %>%
      dplyr::distinct() %>%
      dplyr::inner_join(., full_genotype_pcs) %>%
      dplyr::distinct()
    
    line_prop_changes_LPS = line_prop_changes %>%
      dplyr::filter(treatment == "LPS") %>%
      dplyr::inner_join(., genotype_subset_eqtl) %>%
      dplyr::distinct() %>%
      dplyr::inner_join(., full_genotype_pcs) %>%
      dplyr::distinct()
    
    # List to store the results
    res_untreated = list()
    res_IFN = list()
    res_LPS = list()
    message("Fitting")
    
    for (i in 1:length(unique(genotype_subset_eqtl$rn))) {
      # fit linear mixed model first with random effects, then fit residuals in normal regression
      # this is pre-adjustment for pool and donor effects
      snp =  unique(genotype_subset_eqtl$rn)[i]
      
      
      # subset shared donors to fit random effects
      # explanation random slope vs random intercept: https://www.bristol.ac.uk/cmm/learning/videos/random-slopes.html
      # random pool intercept and only donor in slope?
      
      test = line_prop_changes_untreated[line_prop_changes_untreated$rn == snp, ] %>%
        dplyr::filter(!is.na(sex)) # remove missing data causing problems
      
      my_fit = lm(
        mean_outcome_pool ~ sex  + n_pools + genotypePC1 + genotypePC2 + genotypePC3 + genotypePC4 + genotypePC5 + alt_allele_dosage   ,
        data = test)
      
      sum_res = summary(my_fit)
      res_untreated[[i]] = list(
        coefficients = sum_res$coefficients,
        residuals =  sum_res$residuals,
        call =  sum_res$call      )
     
      test = line_prop_changes_IFN[line_prop_changes_IFN$rn == snp, ] %>%
        dplyr::group_by(line) %>%
        dplyr::filter(!is.na(sex)) # remove missing data causing problems
      # nlme implementation is too slow with just 2 replicates
      my_fit = lm(
        mean_outcome_pool ~ sex  +  n_pools + genotypePC1 + genotypePC2 + genotypePC3 + genotypePC4 + genotypePC5 + 
          alt_allele_dosage ,
        data = test
      )
      sum_res = summary(my_fit)
      res_IFN[[i]] = list(
        coefficients = sum_res$coefficients,
        residuals =  sum_res$residuals,
        call =  sum_res$call
      )
      
      test = line_prop_changes_LPS[line_prop_changes_LPS$rn == snp, ] %>%

        dplyr::filter(!is.na(sex)) # remove missing data causing problems
      my_fit = lm(
        mean_outcome_pool ~ sex  +  n_pools + genotypePC1 + genotypePC2 + genotypePC3 + genotypePC4 + genotypePC5 + 
          alt_allele_dosage ,
        data = test
      )
      
      sum_res = summary(my_fit)
      res_LPS[[i]] = list(
        coefficients = sum_res$coefficients,
        residuals =  sum_res$residuals,
        call =  sum_res$call
      )
     
      
      if (i %% 100 == 0)  {
        print(paste0("i is ", i))
      }
    }
    
   
    p_untreated = list()
    coef_untreated = list()
    se_untreated = list()
    p_IFN = list()
    coef_IFN = list()
    se_IFN = list()
    p_LPS = list()
    coef_LPS = list()
    se_LPS = list()
    for (i in 1:length(res_untreated)) {
      if ("alt_allele_dosage" %in% rownames(res_untreated[[i]]$coefficients)) {
        p_untreated[[i]] = res_untreated[[i]]$coefficients["alt_allele_dosage", "Pr(>|t|)"] # pval
        coef_untreated[[i]] = res_untreated[[i]]$coefficients["alt_allele_dosage", "Estimate"] # coefficient
        se_untreated[[i]] = res_untreated[[i]]$coefficients["alt_allele_dosage", "Std. Error"] # std error
        
      }
      if (!"alt_allele_dosage" %in% rownames(res_untreated[[i]]$coefficients)) {
        p_untreated[[i]] = NA
        coef_untreated[[i]] = NA
        se_untreated[[i]] = NA
        
      }
      if ("alt_allele_dosage" %in% rownames(res_IFN[[i]]$coefficients)) {
        p_IFN[[i]] = res_IFN[[i]]$coefficients["alt_allele_dosage", "Pr(>|t|)"] # pval
        coef_IFN[[i]] = res_IFN[[i]]$coefficients["alt_allele_dosage", "Estimate"] # coefficient
        se_IFN[[i]] = res_IFN[[i]]$coefficients["alt_allele_dosage", "Std. Error"] # std error
      }
      if (!"alt_allele_dosage" %in% rownames(res_IFN[[i]]$coefficients)) {
        p_IFN[[i]] = NA # pval
        coef_IFN[[i]] = NA
        se_IFN[[i]] = NA
      }
      if ("alt_allele_dosage" %in% rownames(res_LPS[[i]]$coefficients)) {
        p_LPS[[i]] = res_LPS[[i]]$coefficients["alt_allele_dosage", "Pr(>|t|)"] # pval
        coef_LPS[[i]] = res_LPS[[i]]$coefficients["alt_allele_dosage", "Estimate"] # coefficient
        se_LPS[[i]] = res_LPS[[i]]$coefficients["alt_allele_dosage", "Std. Error"] # std error
      }
      if (!"alt_allele_dosage" %in% rownames(res_LPS[[i]]$coefficients)) {
        p_LPS[[i]] = NA
        coef_LPS[[i]] = NA
        se_LPS[[i]] = NA
      }
      if (i %% 100 == 0)  {
        print(paste0("i is ", i))
      }
    }
    p_untreated = unlist(p_untreated)
    coef_untreated = unlist(coef_untreated)
    se_untreated = unlist(se_untreated)
    # hist(p_untreated)
    # hist(coef_untreated)
    # hist(se_untreated)
    
    p_IFN = unlist(p_IFN)
    coef_IFN = unlist(coef_IFN)
    se_IFN = unlist(se_IFN)
    p_LPS = unlist(p_LPS)
    coef_LPS = unlist(coef_LPS)
    se_LPS = unlist(se_LPS)
    
    res_df = data.frame(
      coef_untreated = coef_untreated,
      p_untreated = p_untreated,
      se_untreated = se_untreated,
      coef_IFN = coef_IFN,
      p_IFN = p_IFN,
      se_IFN = se_IFN,
      coef_LPS = coef_LPS,
      p_LPS = p_LPS,
      se_LPS = se_LPS,
      snp = unique(genotype_subset_eqtl$rn)
    ) %>%
      distinct()
    
    
    write.csv(res_df,
              output_csv,
              quote = FALSE,
              row.names = FALSE)
  } else{
    # for empty snp lists - not close to any gene tested in eQTL
    res_df = data.frame(
      coef_untreated = NA,
      p_untreated = NA,
      se_untreated = NA,
      coef_IFN = NA,
      p_IFN = NA,
      se_IFN = NA,
      coef_LPS = NA,
      p_LPS = NA,
      se_LPS = NA,
      snp = NA
    )
    write.csv(res_df,
              output_csv,
              quote = FALSE,
              row.names = FALSE)
  }
}else{
  message("No overlapping SNPs - Exiting")
}
