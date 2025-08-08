# line proportion QTL pipeline start
# 3 pools with line mix, migration assay with different conditions
# 

.libPaths(c("/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/",
            "/software/teamtrynka/conda/otar2065/lib/R/library",
            "/software/teamtrynka/ma23/R4.1/libs"))
library(patchwork)
library(tidyverse)
library(lme4) # mixed models
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
library(microbenchmark)
source("../functions.R")

args = commandArgs(trailingOnly=TRUE)


if (length(args)<5) {
  stop("You need to detail 4 input and 1 output file paths.n",
       call. = FALSE)
} else if (length(args) == 5) {
  sample_info_path = args[1]
  line_prop_changes_path = args[2]
  genotype_path = args[3]
  eqtl_path=args[4]
  output_path=args[5]
}


#dir.create(file.path("../../../data/results/1.2.proportion_QTL_mixed_models"),recursive = T)

# to test
#sample_info = read.csv("../../../data/all_pools_migration_sample_info.csv")
#line_prop_changes = read.csv("../../../data/results/1.check_line_proportions/line_prop_changes_per_well.csv")
# load minor allele dosage file for all lines
#
#genotype = readr::read_csv("../../../data/genotypes/full_genotype/genotype_minor_allele_dosage_1.csv") %>%
 # dplyr::relocate(rn)
# load nominal eQTL results and subset to those there - 500kb around expressed genes
# have been combined for all treatments and clusters
#eqtl = read.csv("../../../../OTAR2065_sc_eQTL/data/results/tensorqtl/60/sum_sizefactorsNorm_log2_scaled_centered_unique_nominal_snps.csv")

# sample_info = read.csv(sample_info_path)
# line_prop_changes = read.csv(line_prop_changes_path)
# length(unique(line_prop_changes$line))
# # load minor allele dosage file for all lines
# # 
# genotype = readr::read_csv(genotype_path) %>%
#   dplyr::relocate(rn)
# # load significant eQTL results and subset to those there
# 
# eqtl = read.csv(eqtl_path)

### main ###

# retain only samples with treatments present in all three pools (Media only and C5a - chemoattractant)
line_prop_changes = line_prop_changes %>%
  dplyr::filter(condition %in% c("media_only", "3nM_C5a")) %>%
  unique()

# Retain C5a samples for the moment to check large migration effects
# How to incorporate the Media-only migration rate in this analysis?

# ggplot(line_prop_changes,aes(x=log(mean_max_prop_rep),y = log(var_outcome_rep))) + 
#   geom_point() + theme_minimal()
# negative relationship - lines with larger proportions have smaller variances in their outcomes
# we already new this

# treat lines across pools as replicates, do not average across pools
line_prop_changes = line_prop_changes %>%
  # dplyr::filter(prop_unadjusted_max_value > 0.005) %>%
  dplyr::filter(!log_fraction_mean %in% c(NA,Inf,-Inf)) %>%
  dplyr::select(log_fraction_mean,replicate,line,condition,treatment,pool,sex) %>%
  distinct() 

length(unique(line_prop_changes$line)) # 182 lines without prop filter, 
# with prop filter 1% = 102, with 0.5% = 133

##### just take those with chemoattractant for the time being
line_prop_changes = line_prop_changes %>%
  dplyr::filter(condition == "3nM_C5a")

# checking effects of line and pool on variability
# scatter plots, fraction on the y-axis for both, points grouped by donor for one and pool 
line_var = line_prop_changes %>%
  dplyr::group_by(line) %>%
  dplyr::summarise(fraction_var = var(log_fraction_mean)) 
pool_var = line_prop_changes %>%
  dplyr::group_by(pool) %>%
  dplyr::summarise(fraction_var = var(log_fraction_mean)) 

summary(line_var)
summary(pool_var)
# can't do scatter plot grouping by donor and by pool using all donors because
# donor and pool are highly correlated
# this is an argument to include only pool as random effect
# and then donor would be covered by the genotype PC?

# PCA assessment
pca_data = line_prop_changes %>%
  dplyr::filter(line %in% c("hegp_3","aowh_2") & treatment == "Untreated") %>%
  tidyr::pivot_wider(names_from = replicate, values_from = log_fraction_mean) %>%
  tidyr::drop_na()

metadata = pca_data %>%
  dplyr::select(!c(rep1,rep2,rep3)) %>%
  as.data.frame()

rownames(metadata) = paste(metadata$line,metadata$condition,metadata$treatment,metadata$pool, sep = "_")

pca_data = pca_data %>%
  dplyr::mutate(new_names = paste(line,condition,treatment,pool, sep = "_")) %>%
  dplyr::select(new_names,rep1,rep2,rep3) %>%
  as.data.frame() 

rownames(pca_data) = pca_data$new_names
pca_data = pca_data[-1]
pca_data = t(pca_data)

pca_res = pca_modified(pca_data, metadata = metadata)

biplot(pca_res,colby = "line")
biplot(pca_res,colby = "pool")
biplot(pca_res,colby = "treatment")

# linear regression - mixed models
# explicitly pre-correct for pool effects, then take residuals

# retain lines present in the pools 
###### remove iukl_1 and curn_3 for the time being because I don't have them in the genotype file yet #####
line_prop_changes = line_prop_changes %>%
  dplyr::filter(!line %in% c("iukl_1","curn_3"))
length(unique(line_prop_changes$line))

genotype = genotype %>% 
  dplyr::select(colnames(genotype)[colnames(genotype) %in% c("rn",unique(line_prop_changes$line))]) 
length(colnames(genotype)) - 1
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

# eqtl = eqtl %>%
#   dplyr::mutate(swapped_tensor_alleles = map_chr(variant_id, .f = swap_REF_ALT_alelle_names))
# 
# eqtl = eqtl %>% 
#   dplyr::filter(variant_id %in% unique(genotype_subset$rn) | swapped_tensor_alleles %in% unique(genotype_subset$rn))
# 
# if(nrow(eqtl)!=0){
#   
# 
# genotype_subset_eqtl = genotype_subset %>%
#   dplyr::filter(rn %in% eqtl$variant_id | rn %in% eqtl$swapped_tensor_alleles)

genotype_subset_eqtl = genotype_subset # let's try testing everything, not subsetting to eqtl

if(nrow(genotype_subset_eqtl)!=0){
  length(unique(genotype_subset_eqtl$rn)) # number of SNPs after filtering
  # if subset to eqtl: we've reduced set from 10k to 3.5k approx, so to 35%
  # if not, 
  
  ## add allele dosage, also genotype PCs
  message("Adding genotype PCs")
  
  full_genotype_pcs = read.table("../../../../OTAR2065_sc_eQTL/data/genotype/plink_genotypes/all_pools.genotype.MAF05.eigenvec")
  rownames(full_genotype_pcs) = ifelse(full_genotype_pcs$V1 ==full_genotype_pcs$V2,
                                       yes = full_genotype_pcs$V1,
                                       no = paste(full_genotype_pcs$V1,full_genotype_pcs$V2,sep = "_"))
  full_genotype_pcs = full_genotype_pcs[c(2:7)] %>%
    dplyr::rename(line=V2,genotypePC1 = V3,genotypePC2 = V4,genotypePC3 = V5,genotypePC4 = V6,genotypePC5 = V7)
  
  line_prop_changes_IFN = line_prop_changes %>%
    dplyr::filter(treatment == "IFN") %>%
    dplyr::group_by(line, pool) %>%
    dplyr::filter(n() >2) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(.,genotype_subset_eqtl) %>%
    dplyr::distinct() %>%
    dplyr::inner_join(.,full_genotype_pcs) %>%
    dplyr::distinct() 
  
  line_prop_changes_untreated = line_prop_changes %>%
    dplyr::filter(treatment == "Untreated") %>%
    dplyr::group_by(line, pool) %>%
    dplyr::filter(n() >2) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(.,genotype_subset_eqtl) %>%
    dplyr::distinct() %>%
    dplyr::inner_join(.,full_genotype_pcs) %>%
    dplyr::distinct() 
  
  line_prop_changes_LPS = line_prop_changes %>%
    dplyr::filter(treatment == "LPS") %>%
    dplyr::group_by(line, pool) %>%
    dplyr::filter(n() >2) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(.,genotype_subset_eqtl) %>%
    dplyr::distinct() %>%
    dplyr::inner_join(.,full_genotype_pcs) %>%
    dplyr::distinct() 
  
  # List to store the results
  res_untreated = list()
  res_IFN = list()
  res_LPS=list()
  
  # Is there greater variability in the pool or in the donor component?
  
  result <- microbenchmark(
    
    for(i in 1:length(unique(genotype_subset_eqtl$rn))){
      # fit linear mixed model first with random effects, then fit residuals in normal regression
      # this is pre-adjustment for pool and donor effects
      snp =  unique(genotype_subset_eqtl$rn)[i]
      
      
      # subset shared donors to fit random effects
      # explanation random slope vs random intercept: https://www.bristol.ac.uk/cmm/learning/videos/random-slopes.html
      # random pool intercept and only donor in slope?
      
      test = line_prop_changes_untreated[line_prop_changes_untreated$rn == snp, ] %>%
        dplyr::group_by(line, replicate) %>%
        # dplyr::filter(n() >1) %>% # shared donors
        dplyr::ungroup() %>%
        dplyr::filter(!is.na(sex)) # remove missing data causing problems
      
      
      # log_fraction_mean ~ line +  (1+pool|pool ) also gives singular fit
      # log_fraction_mean ~ sex  +  genotypePC1 + genotypePC2 + genotypePC3 + genotypePC4 + genotypePC5 + alt_allele_dosage  + (1|pool ) + (0+line|pool) also gives singular fit [(0+line|pool) means random slope with no covariance]
      my_fit = lmerTest::lmer(log_fraction_mean ~ sex  +  genotypePC1 + genotypePC2 + genotypePC3 + genotypePC4 + genotypePC5 + alt_allele_dosage + (1|line) + (1|pool) , # log_fraction_mean ~ (...) +  (1+line+pool|pool ) gives singular fit
                              data = test,
                              control = lmerControl(optCtrl=list(maxfun=5000) ))  # higher values don't improve convergence
      myfit_2= lmerTest::lmer(log_fraction_mean ~ sex  +  genotypePC1 + genotypePC2 + genotypePC3 + genotypePC4 + genotypePC5 + alt_allele_dosage + (1+pool|pool) , # log_fraction_mean ~ (...) +  (1+pool|pool ) gives singular fit
                              data = test,
                              control = lmerControl(optCtrl=list(maxfun=5000) ),
                              verbose = 0)  # higher values don't improve convergence
      
      anova_res = anova(my_fit, myfit_2, refit = FALSE)
      p_anova_res[[i]]= anova_res$`Pr(>Chisq)`[2] # there is no significant difference between the model fits (p==1), so we must decide which one is more informative
      # log_fraction_mean ~ sex  +  genotypePC1 + genotypePC2 + genotypePC3 + genotypePC4 + genotypePC5 + alt_allele_dosage + (1|line) + (1|pool) is not singular and gives the lowest REML criterion at convergence
      # can't fit line or pool as random slope - not enough degrees of freedom
      # see https://www.learn-mlms.com/07-module-7.html
      sum_res = summary(my_fit)
      res_untreated[[i]] = list(
        coefficients = sum_res$coefficients,
        varcor =  sum_res$varcor,
        residuals =  sum_res$residuals,
        call =  sum_res$call,
        AICtab =  sum_res$AICtab)
      # testing assumptions
      # https://sscc.wisc.edu/sscc/pubs/MM/MM_DiagInfer.html
      # plot(residuals(my_fit))
      # ggplot(data.frame(x1=test$genotypePC1,pearson=residuals(my_fit,type="pearson")),
      #        aes(x=x1,y=pearson)) +
      #   geom_point() +
      #   theme_bw()
      # 
      # ggplot(data.frame(x2=test$sex,pearson=residuals(my_fit,type="pearson")),
      #        aes(x=x2,y=pearson)) +
      #   geom_point() +
      #   theme_bw()
      # ggplot(data.frame(x1=test$genotypePC1,pearson=residuals(my_fit,type="pearson")),
      #        aes(x=x1,y=pearson)) +
      #   geom_point() +
      #   theme_bw()
      # 
      # ggplot(data.frame(x2=test$alt_allele_dosage,pearson=residuals(my_fit,type="pearson")),
      #        aes(x=x2,y=pearson)) +
      #   geom_point() +
      #   theme_bw()
      # qqnorm(residuals(my_fit))
      test = line_prop_changes_IFN[line_prop_changes_IFN$rn == snp, ] %>%
        dplyr::group_by(line, replicate) %>%
        # dplyr::filter(n() >1) %>% # shared donors
        dplyr::ungroup() %>%
        dplyr::filter(!is.na(sex)) # remove missing data causing problems
      # nlme implementation is too slow with just 2 replicates
      my_fit = lmerTest::lmer(log_fraction_mean ~ sex  +  genotypePC1 + genotypePC2 + genotypePC3 + genotypePC4 + genotypePC5 + alt_allele_dosage + (1|line) + (1|pool), 
                              data = test,
                              control = lmerControl(optCtrl=list(maxfun=10000) ))
      sum_res = summary(my_fit)
      res_IFN[[i]] = list(
        coefficients = sum_res$coefficients,
        varcor =  sum_res$varcor,
        residuals =  sum_res$residuals,
        call =  sum_res$call,
        AICtab =  sum_res$AICtab)
      
      test = line_prop_changes_LPS[line_prop_changes_LPS$rn == snp, ] %>%
        dplyr::group_by(line, replicate) %>%
        # dplyr::filter(n() >1) %>% # shared donors
        dplyr::ungroup() %>%
        dplyr::filter(!is.na(sex)) # remove missing data causing problems
      
      my_fit = lmerTest::lmer(log_fraction_mean ~ sex  +  genotypePC1 + genotypePC2 + genotypePC3 + genotypePC4 + genotypePC5 + alt_allele_dosage + (1|line) + (1|pool), # log_fraction_mean ~ line + pool +  (1+line+pool|pool ) gives singular fit
                              data = test,
                              control = lmerControl(optCtrl=list(maxfun=10000) ))
      
      sum_res = summary(my_fit)
      res_LPS[[i]] = list(
        coefficients = sum_res$coefficients,
        varcor =  sum_res$varcor,
        residuals =  sum_res$residuals,
        call =  sum_res$call,
        AICtab =  sum_res$AICtab)
      
      # nlme::lme ( fixed=log_fraction_mean ~ line, 
      #             random=~1+line|line,
      #             data = test,
      #             control =  nlme::lmeControl(opt='optim'))
      # this doesn't converge
      
      
      if(i%% 1000 ==0)  { print(paste0("i is ",i))}
    }
    
    ,times=1)
  print(result) # around 9h for 3 treatments (20h for singular fits)
  # nlme is slower - takes long time to converge because there are only two replicates
  
  # save these as rds
  #saveRDS(object = res_untreated,paste0(output_path,"/res_untreated.rds"))
  #saveRDS(object = res_IFN,paste0(output_path,"/res_IFN.rds"))
  #saveRDS(object = res_LPS,paste0(output_path,"/res_LPS.rds"))
  #saveRDS(object = p_anova_res,paste0(output_path,"/p_anova_res.rds")) # all 1s
  # identify per chunk of snps
  
  # ggplot(test, aes(x=alt_allele_dosage, y=log_fraction_mean, col = pool)) + 
  #   geom_point() + theme_minimal()
  # 5 min for ~5k snps and 2 cores
  # there are singularities in the result - can't fit line and pool at the same time
  p_untreated = list()
  coef_untreated = list()
  se_untreated = list()
  p_IFN = list()
  coef_IFN = list()
  se_IFN= list()
  p_LPS = list()
  coef_LPS = list()
  se_LPS= list()
  for(i in 1:length(res_untreated)){
    p_untreated[[i]] = res_untreated[[i]]$coefficients["alt_allele_dosage","Pr(>|t|)"] # pval
    coef_untreated[[i]] = res_untreated[[i]]$coefficients["alt_allele_dosage","Estimate"] # coefficient
    se_untreated[[i]] = res_untreated[[i]]$coefficients["alt_allele_dosage","Std. Error"] # std error
    
    p_IFN[[i]] = res_IFN[[i]]$coefficients["alt_allele_dosage","Pr(>|t|)"] # pval
    coef_IFN[[i]] = res_IFN[[i]]$coefficients["alt_allele_dosage","Estimate"] # coefficient
    se_IFN[[i]] = res_IFN[[i]]$coefficients["alt_allele_dosage","Std. Error"] # std error
    
    p_LPS[[i]] = res_LPS[[i]]$coefficients["alt_allele_dosage","Pr(>|t|)"] # pval
    coef_LPS[[i]] = res_LPS[[i]]$coefficients["alt_allele_dosage","Estimate"] # coefficient
    se_LPS[[i]] = res_LPS[[i]]$coefficients["alt_allele_dosage","Std. Error"] # std error
    
    if(i%% 1000 ==0)  { print(paste0("i is ",i))}
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
  
  
  res_df = data.frame(coef_untreated=coef_untreated,p_untreated=p_untreated,se_untreated = se_untreated,
                      coef_IFN = coef_IFN,p_IFN = p_IFN,se_IFN = se_IFN,
                      coef_LPS = coef_LPS,p_LPS = p_LPS,se_LPS = se_LPS, snp = unique(genotype_subset_eqtl$rn)) %>%
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
