library(patchwork)
library(tidyverse)
source("./functions.R")

input_dir = "../../data/results/1.alluvial_plots"
output_dir = "../../data/results/2.associations_preMAC_microglia_proliferation"
dir.create(output_dir, recursive = TRUE)
data = read.table(paste0(input_dir,
                         "/pools2-11_13-15_changing_props_iPSC_preMacs_microglia_WGS_sc.txt"),
                  header = TRUE)

# read in proportion error estimates
error_estimates = read.table("../../../OTAR2065_phenotypic_QTLs/data/w/error_approximations_generic_pool.txt", header = TRUE) %>%
  dplyr::filter(coverage ==5 & genotype =="old")

# Perform the inner join based on the closest proportion values
data_with_error =  data %>%
  dplyr::rowwise() %>%
  dplyr::mutate(closest_error_proportion = find_closest_value(prop, error_estimates$w_real)) %>%
  dplyr::ungroup() %>%
  dplyr::inner_join(error_estimates, by = c("closest_error_proportion" = "w_real")) %>%
dplyr::mutate(prop_adjusted_mean = dplyr::case_when(prop <0.1 ~ prop - mean_wdif,
                                                    prop >=0.1 ~ prop + mean_wdif),
              prop_adjusted_median = dplyr::case_when(prop <0.1 ~ prop - median_wdif,
                                                      prop >=0.1 ~ prop + median_wdif),
              prop_adjusted_max = dplyr::case_when(prop <0.1 ~ prop - max_wdif,
                                                   prop >=0.1 ~ prop + max_wdif))
# if any estimate becomes negative, make it zero:
data_with_error[data_with_error$prop_adjusted_mean<0, "prop_adjusted_mean"] = 0
data_with_error[data_with_error$prop_adjusted_median<0, "prop_adjusted_median"] = 0
data_with_error[data_with_error$prop_adjusted_max<0, "prop_adjusted_max"] = 0

# need to distinguish between those failed at point of preMACS and at point of microglia
# If Inf, then they were 0 at preMAC
# If 0, they disappear at microglia
efficiency = prolif_microglia_premac(data_with_error)
 summary(efficiency)
 efficiency = efficiency %>%
   dplyr::filter(!scaled_proportion %in% c(Inf,NA,NaN)) %>%
   dplyr::filter(line!="sh5y5y") %>% # not interested in phagocytosis line
   dplyr::mutate(sequencing = case_when(stringr::str_detect(sample,"sc") ~ "single_cell",
                 stringr::str_detect(sample,"phago") ~ "WGS",
                 stringr::str_detect(sample,"migr") ~ "WGS"),
                 treatment = case_when(stringr::str_detect(sample,"untreated") ~ "untreated",
                                       stringr::str_detect(sample,"IFN") ~ "IFN",
                                       stringr::str_detect(sample,"LPS") ~ "LPS")) 
 summary(efficiency$scaled_proportion)
 hist(efficiency$scaled_proportion, breaks = 100)
 # where do the NAs come from?
 # this is not normal at all
 hist(log10(efficiency$scaled_proportion),breaks = 40)
 # still very extreme values
 # model without transforming via GLMM?

 # diagnostic plots
 
 # Run a diagnostic lot of the group variances vs group
 # means (genotype x nutrient x clipping grouping).  Code
 # used to produce the plot :
 # https://github.com/QCBSRworkshops/workshop07/blob/main/pres-fr/data/glmm_e.r

 # without logging
 
 summary_tibble = efficiency %>%
   dplyr::group_by(line,treatment,pool) %>%
   dplyr::summarise(means =mean(scaled_proportion),
                    vars = var(scaled_proportion))
 # dplyr::filter(!log_means %in% c(-Inf,Inf,NA,NaN))
 vars = summary_tibble$vars
 means = summary_tibble$means
 
 # Quasi-Poisson
 lm1 = lm(vars ~ means - 1)
 phi.fit = coef(lm1)
 # The -1 specifies a model with the intercept set to zero
 
 # Negative binomial
 lm2 = lm(vars ~ I(means^2) + offset(means) - 1)
 k.fit = 1/coef(lm2)
 # The offset() is used to specify that we want the group
 # means added as a term with its coefficient fixed to 1
 
 # Non-parametric loess fit
 Lfit = loess(vars ~ means)
 
 # The plot
 plot(vars ~ means, xlab = "Group means", ylab = "Group variances ")
 abline(a = 0, b = 1, lty = 2)
 text(-4, -5, "Poisson")
 curve(phi.fit * x, col = 2, add = TRUE)
 # bquote() is used to substitute numeric values in
 # equations with symbols
 text(4, 4, bquote(paste("QP: ", sigma^2 == .(round(phi.fit,
                                                    1)) * mu)), col = 2)
 curve(x * (1 + x/k.fit), col = 4, add = TRUE)
 text(5, 1, paste("NB: k = ", round(k.fit, 1), sep = ""),
      col = 4)
 mvec <- 0:120
 lines(mvec, predict(Lfit, mvec), col = 5)
 text(5, 6, "loess", col = 5)
 
 # logging
 summary_tibble = efficiency %>%
   dplyr::group_by(line,treatment,pool) %>%
   dplyr::summarise(log_means = log(mean(scaled_proportion)),
                    log_vars = log(var(scaled_proportion))) %>%
   dplyr::filter(!log_means %in% c(-Inf,Inf,NA,NaN))
 log_vars = summary_tibble$log_vars
 log_means = summary_tibble$log_means

 # normal
 lm0 = lm(log_vars ~ log_means)
 fit = coef(lm0)
 
 # Quasi-Poisson
 lm1 = lm(log_vars ~ log_means - 1)
 phi.fit = coef(lm1)
 # The -1 specifies a model with the intercept set to zero
 
 # Negative binomial
 lm2 = lm(log_vars ~ I(log_means^2) + offset(log_means) - 1)
 k.fit = 1/coef(lm2)
 # The offset() is used to specify that we want the group
 # means added as a term with its coefficient fixed to 1
 
 # Non-parametric loess fit
 Lfit = loess(log_vars ~ log_means)
 
 # The plot
 plot(log_vars ~ log_means, xlab = "Group means (log)", ylab = "Group variances (log)")
 abline(a = 0, b = 1, lty = 2)
 text(-4, -5, "Poisson")
 
 abline(fit) 
 text(-4, -11, "Normal")
 
 curve(phi.fit * x, col = 2, add = TRUE)
 # bquote() is used to substitute numeric values in
 # equations with symbols
 text(4, 4, bquote(paste("QP: ", sigma^2 == .(round(phi.fit,
                                                         1)) * mu)), col = 2)
 curve(x * (1 + x/k.fit), col = 4, add = TRUE)
 text(5, 1, paste("NB: k = ", round(k.fit, 1), sep = ""),
      col = 4)
 mvec <- 0:120
 lines(mvec, predict(Lfit, mvec), col = 5)
 text(5, 6, "loess", col = 5)
 
 # Quasi-Poisson, normal or loess?
 #https://r.qcbs.ca/workshop07/book-en/choose-an-error-distribution.html
 
 
 
 ####### association tests
 line_info=read.csv("../../../OTAR2065_phenotypic_QTLs/data/allinfo_hipsci_PD_AD_PRS_nondisease_feederfree_european.csv") %>%
   dplyr::mutate(sex = case_when(grepl(pattern = "Female",.$Sex) ~ "Female",
                                 grepl(pattern = "Male",.$Sex) ~ "Male")) %>%
   dplyr::select(Line,sex) %>%
   dplyr::rename(line=Line)
 
 efficiency = efficiency %>%
   dplyr::left_join(.,line_info)
 # load minor allele dosage file for all lines

 genotype = readr::read_csv("../../../OTAR2065_phenotypic_QTLs/data/genotypes/full_genotype/genotype_minor_allele_dosage_1.csv") %>%
   dplyr::relocate(rn) %>%
   dplyr::ungroup() %>%
   #dplyr::filter(rn %in% to_subset) %>%
   tidyr::pivot_longer(cols = !rn, values_to = "alt_allele_dosage", names_to = "line") %>%
   dplyr::mutate(alt_allele_dosage = dplyr::case_when(
     alt_allele_dosage == 0.5 ~ 1,
     alt_allele_dosage == 1   ~ 2,
     TRUE          ~ alt_allele_dosage  # Recode so var of allele dosage is 1, Keep 0 unchanged
   ))
 
 length(unique(genotype$rn)) # number of SNPs after filtering
 
 # remove genotypes under 10% MAF
 snps_to_keep = genotype %>%
   dplyr::mutate(alt_allele_dosage = as.factor(alt_allele_dosage)) %>%
   count(rn,alt_allele_dosage) %>%
   dplyr::group_by(rn) %>%
   dplyr::filter(min(n) >= (0.1 * length(unique(genotype$line)))) %>%
   ungroup() %>%
   dplyr::select(rn) %>%
   distinct() %>%
   unlist()
 
 
 genotype_subset = genotype %>%
   dplyr::filter(rn %in% snps_to_keep) %>%
   dplyr::mutate(rn = str_replace(rn, "chr", "")) 
 
 length(unique(genotype_subset$rn)) # number of SNPs after filtering
 
 if(nrow(genotype_subset)!=0){
   
 message("Adding genotype PCs")
 
 full_genotype_pcs = read.table("../../../OTAR2065_sc_eQTL/data/genotype/plink_genotypes/all_pools.genotype.MAF05.eigenvec")
 rownames(full_genotype_pcs) = ifelse(full_genotype_pcs$V1 ==full_genotype_pcs$V2,
                                      yes = full_genotype_pcs$V1,
                                      no = paste(full_genotype_pcs$V1,full_genotype_pcs$V2,sep = "_"))
 full_genotype_pcs = full_genotype_pcs[c(2:7)] %>%
   dplyr::rename(line=V2,genotypePC1 = V3,genotypePC2 = V4,genotypePC3 = V5,genotypePC4 = V6,genotypePC5 = V7)
 
 efficiency_w_info = efficiency %>%
   dplyr::inner_join(.,genotype_subset) %>%
   dplyr::distinct() %>%
   dplyr::inner_join(.,full_genotype_pcs) %>%
   dplyr::distinct() 
 
 rm(genotype,full_genotype_pcs)
 gc()
 
 # cleaning up Inf and NAs
 efficiency_w_info = efficiency_w_info %>%
   dplyr::mutate(log_scaled_proportion = log(scaled_proportion)) %>%
   dplyr::filter(!log_scaled_proportion %in% c(-Inf,Inf,NA,NaN)) 
 
 efficiency_w_info = efficiency_w_info %>%
   group_by(rn) %>%
   filter(n_distinct(alt_allele_dosage) >= 2) %>%
   ungroup() %>%
   distinct() %>%
   dplyr::mutate(treatment = factor(treatment)) %>% # setting untreated as baseline
   dplyr::mutate(treatment = relevel(treatment,ref = c("untreated"))) # setting untreated as baseline
 length(unique(efficiency_w_info$rn)) # number of SNPs after filtering
 
 # List to store the results
 res = list()
 result <- microbenchmark(
 for(i in 1:length(unique(efficiency_w_info$rn))){
   snp =  unique(efficiency_w_info$rn)[i]
   
   test = efficiency_w_info[efficiency_w_info$rn == snp, ] %>%
     # dplyr::filter(treatment == "untreated") %>%
     dplyr::filter(!is.na(sex)) # remove missing data causing problems
   
 
 my_fit =  lmerTest::lmer(log_scaled_proportion ~  genotypePC1 + genotypePC2 + genotypePC3 + genotypePC4 + genotypePC5 + 
                            sex + sequencing + 
                            treatment*alt_allele_dosage + # effect of each unit of change in allele dosage in treatment vs baseline
                            #(1|sex) + (1 | sequencing) +  
                            (1|line) + (1|pool),
              data = test, 
              control = lmerControl(optCtrl=list(maxfun=5000) ))
 
 sum_res = summary(my_fit)
 res[[i]] = list(
   coefficients = sum_res$coefficients,
   varcor =  sum_res$varcor,
   residuals =  sum_res$residuals,
   call =  sum_res$call,
   AICtab =  sum_res$AICtab)
 
 if(i%% 1000 ==0)  { print(paste0("i is ",i))}
 }
 
 ,times=1)
print(result)
 }
 
 # save these as rds
 saveRDS(object = res,paste0(output_dir,"/res_untreated.rds"))
# incorporate maximum proportion information in each pair, to then fit as covariate
# also associations with preMAC aging - timecourse model