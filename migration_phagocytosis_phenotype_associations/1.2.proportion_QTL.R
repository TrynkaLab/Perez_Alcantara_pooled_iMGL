# Donor proportion QTL pipeline start
# 3 pools with donor mix, migration assay with different conditions
# 

.libPaths(c("/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/",
            "/software/teamtrynka/conda/otar2065/lib/R/library",
            "/software/teamtrynka/ma23/R4.1/libs"))
library(patchwork)
library(tidyverse)
# library(afex) # mixed models
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
library(foreach)
library(doParallel)
library(qvalue)
source("./functions.R")

args = commandArgs(trailingOnly=TRUE)


if (length(args)<5) {
  stop("You need to detail 4 input and 1 output file paths.n",
       call. = FALSE)
} else if (length(args) == 5) {
  sample_info_path = args[1]
  donor_prop_changes_path = args[2]
  genotype_path = args[3]
  eqtl_path=args[4]
  output_path=args[5]
}


dir.create(file.path("../../data/results/5.proportion_QTL"),recursive = T)

# to test
# sample_info = read.csv("../../data/all_pools_migration_sample_info.csv")
# donor_prop_changes = read.csv("../../data/results/4.check_donor_proportions/line_prop_changes_averages_variances_1pct_filtered.csv")
# # load minor allele dosage file for all donors
# # 
# genotype = readr::read_csv("../../data/genotypes/full_genotype/genotype_minor_allele_dosage_1.csv") %>%
#   dplyr::relocate(rn)
# load nominal eQTL results and subset to those there - 500kb around expressed genes 
# have been combined for all treatments and clusters
# eqtl = read.csv("../../../OTAR2065_sc_eQTL/data/results/tensorqtl/60/sum_sizefactorsNorm_log2_scaled_centered_unique_nominal_snps.csv")

sample_info = read.csv(sample_info_path)
donor_prop_changes = read.csv(donor_prop_changes_path)
length(unique(donor_prop_changes$line))
# load minor allele dosage file for all donors
# 
genotype = readr::read_csv(genotype_path) %>%
  dplyr::relocate(rn)
# load significant eQTL results and subset to those there

eqtl = read.csv(eqtl_path)

### main ###

# retain only samples with treatments present in all three pools (Media only and C5a - chemoattractant)
donor_prop_changes = donor_prop_changes %>%
  dplyr::filter(condition %in% c("media_only", "3nM_C5a")) %>%
  unique()

# Retain C5a samples for the moment to check large migration effects
# How to incorporate the Media-only migration rate in this analysis?

# gather averages, for more than one pool say that pool is "shared"
# remove mean duplicates
donor_prop_changes = donor_prop_changes %>%
  dplyr::filter(!mean_outcome_pool %in% c(NA,Inf,-Inf)) %>%
  dplyr::mutate(pool= case_when(n_pools>1 ~ "shared",
                                n_pools==1 ~pool)) %>%
  dplyr::select(mean_outcome_pool,var_outcome_pool,line,condition,treatment,pool,n_pools,sex) %>%
  distinct() %>%
  dplyr::mutate(var_outcome_pool= case_when(is.na(var_outcome_pool) ~ 1,
                                            !is.na(var_outcome_pool) ~ var_outcome_pool))

length(unique(donor_prop_changes$line)) # 78 lines

# filter to remove extreme proportions

donor_prop_changes = donor_prop_changes %>%
  dplyr::filter(condition == "3nM_C5a")


# linear regression

# retain donors present in the pools 
genotype = genotype %>% 
  dplyr::select(colnames(genotype)[colnames(genotype) %in% c("rn",unique(donor_prop_changes$line))]) 

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

# remove genotypes under 10% MAF
snps_to_keep = genotype_subset %>%
  dplyr::mutate(alt_allele_dosage = as.factor(alt_allele_dosage)) %>%
  count(rn,alt_allele_dosage) %>%
  dplyr::group_by(rn) %>%
  dplyr::filter(min(n) >= (0.1 * length(unique(genotype_subset$line)))) %>%
  ungroup() %>%
  dplyr::select(rn) %>%
  distinct() %>%
  unlist()


genotype_subset = genotype_subset %>%
  dplyr::filter(rn %in% snps_to_keep) %>%
  dplyr::mutate(rn = str_replace(rn, "chr", "")) 

length(unique(genotype_subset$rn)) # number of SNPs after filtering


# tensorqtl swaps alleles so that the minor allele is always the less abundant one - checking in eQTL results
# because many are swapped compared to the original VCF
# this might take a while 

eqtl = eqtl %>%
  dplyr::mutate(swapped_tensor_alleles = map_chr(variant_id, .f = swap_REF_ALT_alelle_names))

eqtl = eqtl %>% 
  dplyr::filter(variant_id %in% unique(genotype_subset$rn) | swapped_tensor_alleles %in% unique(genotype_subset$rn))

if(nrow(eqtl)!=0){
  
  
  genotype_subset_eqtl = genotype_subset %>%
    dplyr::filter(rn %in% eqtl$variant_id | rn %in% eqtl$swapped_tensor_alleles)
  
  length(unique(genotype_subset_eqtl$rn)) # number of SNPs after filtering
  # we've reduced set from 10k to 3.5k approx, so to 35%
  
  ## add allele dosage, also genotype PCs
  message("Adding genotype PCs")
  
  full_genotype_pcs = read.table("../../../OTAR2065_sc_eQTL/data/genotype/plink_genotypes/all_pools.genotype.MAF05.eigenvec")
  rownames(full_genotype_pcs) = ifelse(full_genotype_pcs$V1 ==full_genotype_pcs$V2,
                                       yes = full_genotype_pcs$V1,
                                       no = paste(full_genotype_pcs$V1,full_genotype_pcs$V2,sep = "_"))
  full_genotype_pcs = full_genotype_pcs[c(2:7)] %>%
    dplyr::rename(line=V2,genotypePC1 = V3,genotypePC2 = V4,genotypePC3 = V5,genotypePC4 = V6,genotypePC5 = V7)
  
  
  line_prop_changes_IFN = donor_prop_changes %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(.,genotype_subset_eqtl) %>%
    dplyr::filter(treatment == "IFN") %>%
    dplyr::distinct() %>%
    dplyr::inner_join(.,full_genotype_pcs) %>%
    dplyr::distinct()
  
  
  line_prop_changes_untreated = donor_prop_changes %>%
    dplyr::ungroup() %>%
    dplyr::left_join(.,genotype_subset_eqtl) %>%
    dplyr::filter(treatment == "Untreated")  %>%
    dplyr::distinct() %>%
    dplyr::inner_join(.,full_genotype_pcs) %>%
    dplyr::distinct()
  
  line_prop_changes_LPS = donor_prop_changes %>%
    dplyr::ungroup() %>%
    dplyr::left_join(.,genotype_subset_eqtl) %>%
    dplyr::filter(treatment == "LPS") %>%
    dplyr::distinct() %>%
    dplyr::inner_join(.,full_genotype_pcs) %>%
    dplyr::distinct()
  
  # test mixed model
  # checking distribution again
  # 
  # p1 = rbind(subset_IFN,subset_untreated) %>%
  #   dplyr::select(well:pool) %>%
  #   unique() %>%
  #   ggplot(., aes(bottom_top_fraction)) +
  #   geom_histogram(bins = 10) +
  #   geom_vline(xintercept = 1,col = "red") +
  #   facet_wrap(vars(treatment), scales = "free_x") + 
  #   theme_bw()
  # p2 = rbind(subset_IFN,subset_untreated) %>%
  #   dplyr::select(well:pool) %>%
  #   unique() %>%
  #   ggplot(., aes(log(bottom_top_fraction))) +
  #   geom_histogram(bins = 10) +
  #   geom_vline(xintercept = 0,col = "red") +
  #   facet_wrap(vars(treatment), scales = "free_x") + 
  #   theme_bw()
  # 
  # png("../../data/QTL_results/1.proportion_QTL/histogram_fractions_C5a_postfilter.png",
  #     width = 9, height = 4, units = "in", res = 400)
  # p1 / p2 + patchwork::plot_annotation(title = "Distribution of top/bottom fraction",
  #                                      subtitle = "With C5a. As-is (top) or logged (bottom). Red line = no change.")
  # 
  # dev.off()
  
  # subset of genotype first
  
  
  # List to store the results
  
  res_IFN = list()
  # res_IFN_afex = list()
  res_untreated = list()
  res_LPS=list()
  
  registerDoParallel(cores=30)
  # parallelized
  # res_IFN = list()
  # res_untreated = list()
  # res_LPS=list()
  
  res_untreated<- foreach(i = 1:length(unique(genotype_subset_eqtl$rn)), .combine = append, .packages = c("data.table")) %dopar%
    {      
      snp =  unique(genotype_subset_eqtl$rn)[i]
      my_fit <- lm("mean_outcome_pool ~ pool + sex + genotypePC1 + genotypePC2 + genotypePC3 + genotypePC4 + genotypePC5 + alt_allele_dosage", 
                   data = line_prop_changes_untreated[line_prop_changes_untreated$rn == snp, ])     
      sum = summary(my_fit)
      par_results <- list(
        coefs =sum$coefficients[,1],
        se = sum$coefficients[,2],
        fit=my_fit
      )
      par_results <- list(par_results)
      names(par_results) <- snp
      return(par_results)
    }
  
  res_IFN <- foreach(i = 1:length(unique(genotype_subset_eqtl$rn)), .combine = append, .packages = c("data.table")) %dopar%
    {      
      snp =  unique(genotype_subset_eqtl$rn)[i]
      # my_fit <- lm("mean_outcome_pool ~ pool + sex + var_outcome_pool + genotypePC1 + genotypePC2 + genotypePC3 + genotypePC4 + genotypePC5 + alt_allele_dosage", 
      #              data = line_prop_changes_IFN[line_prop_changes_IFN$rn == snp, ]) 
      my_fit <- lm("mean_outcome_pool ~ pool + sex + genotypePC1 + genotypePC2 + genotypePC3 + genotypePC4 + genotypePC5 + alt_allele_dosage", 
                   data = line_prop_changes_IFN[line_prop_changes_IFN$rn == snp, ])    
      sum = summary(my_fit)
      par_results <- list(
        coefs =sum$coefficients[,1],
        se = sum$coefficients[,2],
        fit=my_fit
      )
      par_results <- list(par_results)
      names(par_results) <- snp
      return(par_results)
    }
  
  
  res_LPS <- foreach(i = 1:length(unique(genotype_subset_eqtl$rn)), .combine = append, .packages = c("data.table")) %dopar%
    {      
      snp =  unique(genotype_subset_eqtl$rn)[i]
      my_fit <- lm("mean_outcome_pool ~ pool + sex + genotypePC1 + genotypePC2 + genotypePC3 + genotypePC4 + genotypePC5 + alt_allele_dosage", 
                   data = line_prop_changes_LPS[line_prop_changes_LPS$rn == snp, ])     
      sum = summary(my_fit)
      par_results <- list(
        coefs =sum$coefficients[,1],
        se = sum$coefficients[,2],
        fit=my_fit
      )
      par_results <- list(par_results)
      names(par_results) <- snp
      return(par_results)
    }
  
  stopImplicitCluster()
  # 5 min for ~5k snps and 2 cores
  # there are singularities in the result - can't fit line and pool at the same time
  
  # sometimes bigglm throws error: Error in coef.bigqr(object$qr) : 
  #NA/NaN/Inf in foreign function call (arg 5)
  # likely due to seriously imbalanced categories
  # biglm does not throw this error
  #with(line_prop_changes_IFN[line_prop_changes_IFN$rn == snp,],table(alt_allele_dosage,pool))
  
  # save.image( file = "~/proportion_allele_effect.RData")
  # load("~/proportion_allele_effect.RData")
  
  
  # check significant SNPs and plot them
  # how to use the media only vs C5a effect?
  # should alt allele dosage be a factor?
  # res_IFN_afex[[snp]]$anova_table # effect of allele dosage vs reduced model w. intercept
  # res_IFN_afex[[snp]]$full_model #fit of full model
  # summary(res_IFN_afex[[snp]])$varcor 
  #summary(res_IFN[[snp]]$fit)
  #anova(res_IFN[[snp]])
  # extract p-value for allelic dosage coefficient
  # p_IfN_afex  = sapply(res_IFN_afex, FUN = function(x) x$anova_table$`Pr(>F)`[8])
  # print example
  message("The anova structure looks like: ")
  
  print(anova(res_untreated[[1]]$fit))
  print(anova(res_untreated[[1]]$fit)$`Pr(>F)`[8])
  print(summary(res_untreated[[1]]$fit))
  
  p_IFN = sapply(res_IFN,FUN = function(x) anova(x$fit)$`Pr(>F)`[8]) # always double-check you're taking the correct row
  p_untreated = sapply(res_untreated,FUN = function(x) anova(x$fit)$`Pr(>F)`[8])
  p_LPS= sapply(res_LPS,FUN = function(x) anova(x$fit)$`Pr(>F)`[8])
  
  # coefficients
  coef_IFN = sapply(res_IFN,FUN = function(x) x$coef["alt_allele_dosage"]) # always double-check you're taking the correct row
  names(coef_IFN) = names(res_IFN)
  coef_untreated = sapply(res_untreated,FUN = function(x) x$coef["alt_allele_dosage"])
  names(coef_untreated) = names(res_untreated)
  coef_LPS= sapply(res_LPS,FUN = function(x) x$coef["alt_allele_dosage"])
  names(coef_LPS) = names(res_LPS)
  
  # standard errors
  se_IFN = sapply(res_IFN,FUN = function(x) x$se["alt_allele_dosage"]) # always double-check you're taking the correct row
  names(se_IFN) = names(res_IFN)
  se_untreated = sapply(res_untreated,FUN = function(x) x$se["alt_allele_dosage"])
  names(se_untreated) = names(res_untreated)
  se_LPS= sapply(res_LPS,FUN = function(x) x$se["alt_allele_dosage"])
  names(se_LPS) = names(res_LPS)
  
  
  res_df = data.frame(coef_untreated=coef_untreated,p_untreated=p_untreated,se_untreated = se_untreated,
                      coef_IFN = coef_IFN,p_IFN = p_IFN,se_IFN = se_IFN,
                      coef_LPS = coef_LPS,p_LPS = p_LPS,se_LPS = se_LPS) %>%
    tibble::rownames_to_column(var = "snp") %>%
    distinct()
  write.csv(res_df,output_path,quote = FALSE,row.names = FALSE)
}else{
  # for empty snp lists - not close to any gene tested in eQTL
  res_df = data.frame(coef_untreated=NA,p_untreated=NA,se_untreated = NA,
                      coef_IFN = NA,p_IFN = NA,se_IFN = NA,
                      coef_LPS = NA,p_LPS = NA,se_LPS = NA,
                      snp=NA)
  write.csv(res_df,output_path,quote = FALSE,row.names = FALSE)
}
