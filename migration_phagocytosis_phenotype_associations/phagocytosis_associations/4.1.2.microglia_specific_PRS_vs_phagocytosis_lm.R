
library(patchwork)
library(tidyverse)
library(jtools)
source("functions.R")

output_path = "../../../data/results/phagocytosis/4.1.2.microglia_specific_PRS_vs_phagocytosis_lm/"
dir.create(output_path,recursive = T)

sample_info_path = "../../../data/all_pools_phagocytosis_sample_info.csv"
line_prop_changes_path = "../../../data/results/phagocytosis/1.check_line_proportions/line_prop_changes_averages_variances_1pct_filtered.csv"
#genotype_pc_path = "../../../../OTAR2065_sc_eQTL/data/genotype/plink_genotypes/all_pools.no_outliers.genotype.MAF05.eigenvec"
# no need to fit genotype PC because PRS already adjusts for this
other_info_path = "../../../../OTAR2065_differentiation_efficiency/data/donor_metadata_complete_with_imputed_sex.csv"

prs = list()
for(treat in c("untreated","IFN","LPS")){
  prs_path = paste0("../../../../hipsci_genotype_processing/data/prs_ad_bellenguez/PRSice/hipsci_polygenic_score_AD_Bellenguez_withAPOE_microglia_",treat,".tsv")
  
  prs[[treat]] = read_tsv(prs_path)  %>%
    dplyr::select(line, APOE_sum_scaled, prs_scaled, full_PRS) %>%
    dplyr::mutate(treatment = treat) %>%
    dplyr::rename(donor=line) %>%
    dplyr::mutate(line = case_when(donor== "Arlene" ~ "Arlene-003",
                                   donor=="Cindy" ~ "Cindy-005",
                                   donor=="Dexter" ~ "Dexter-006",
                                   donor=="Fiona" ~ "Fiona-010",
                                   donor=="Gilma" ~ "Gilma-009",
                                   donor=="Hector" ~ "Hector-011",
                                   donor=="Imani" ~ "Imani-012",
                                   donor=="Javier" ~ "Javier-013",
                                   donor=="Keoni" ~ "Keoni-014",
                                   donor=="Olaf" ~ "Olaf-018",
                                   donor=="Bertha" ~ "Bertha-004",
                                   donor=="Mindy" ~ "Mindy-016",
                                   donor=="Qiana" ~ "Qiana-022",
                                   donor=="Nestor" ~ "Nestor-017",
                                   .default = donor)) %>%
    dplyr::mutate(donor =  str_split_i(donor,"_",i=1))
}

prs = do.call("rbind",prs)


sample_info = read_csv(sample_info_path) 
other_info = read_csv(other_info_path)
line_prop_changes = read_csv(line_prop_changes_path) 

message("There are ", length(unique(line_prop_changes$line)), " lines in the analysis") # 143

line_prop_changes = line_prop_changes %>%
  # dplyr::filter(prop_unadjusted_max_value > 0.005) %>%
  dplyr::filter(!mean_outcome_pool %in% c(NA,Inf,-Inf)) %>%
  dplyr::select(mean_outcome_pool,line,treatment,n_pools,mean_min_prop_pool) %>%
  distinct() 

message("There are ", length(unique(line_prop_changes$line)), " lines in the analysis after Na and Inf filters") # 143

### adding PRS and other info
line_prop_changes = line_prop_changes %>%
  dplyr::left_join(prs) %>%
  dplyr::left_join(other_info, by = "donor") %>%
  dplyr::rename(sex=Sex)

### examining correlations


info = line_prop_changes %>%
  dplyr::distinct() %>%
  dplyr::select(treatment,mean_outcome_pool, line,n_pools,sex,mean_min_prop_pool,APOE_sum_scaled,prs_scaled,full_PRS) %>%
  as.data.frame()



form = ~ treatment + n_pools + mean_min_prop_pool + mean_outcome_pool  + APOE_sum_scaled + prs_scaled + full_PRS

# Compute Canonical Correlation Analysis (CCA)
# between all pairs of variables
# returns absolute correlation value
C = canCorPairs(form, info)

# Plot correlation matrix
# between all pairs of variables

pdf(paste0(output_path,"PRS_phagocytosis_canonical_correlation_analysis.pdf"),width = 6,height = 6)

plotCorrMatrix(C)

dev.off()
# I will include genotype PCs 1 and 3 because they are less correlated than the rest of the variables, and 
# won't include replicate because the categories are meaningless, the model will just see there are two replicates within each genotype PC and pool
# won't include age because high PRS IPMAR lines are quite old

# List to store the results
res = list()
res[["untreated"]] = list()
res[["IFN"]] = list()
res[["LPS"]] = list()
p_list = list()
message("Fitting")

for(treat in c("untreated","LPS","IFN")){
  totest = line_prop_changes %>%
    dplyr::filter(treatment == treat) %>%
    dplyr::distinct() 
  
  my_fit = lm(
    mean_outcome_pool ~ sex   + mean_min_prop_pool + full_PRS  ,
    data = totest
  )
  # checking regression plot with partial residuals (accounting for all covariates)
  p_list[[paste0(treat,"_full_PRS")]] = jtools::effect_plot(my_fit, pred = full_PRS, interval = FALSE, partial.residuals = TRUE) + 
    ggtitle("")
  
  sum_res = summary(my_fit)
  res[[treat]][["full_PRS"]] = list(
    coefficients = sum_res$coefficients,
    varcor =  sum_res$varcor,
    residuals =  sum_res$residuals,
    call =  sum_res$call,
    AICtab =  sum_res$AICtab,
    r2=sum_res$r.squared,
    adjr2=sum_res$adj.r.squared,
    treatment = treat,
    PRS = "full_PRS"
  )
  
  my_fit =lm(
    mean_outcome_pool ~ sex   + mean_min_prop_pool + prs_scaled  ,
    data = totest
  )
  
  p_list[[paste0(treat,"_polygenicHR")]] = jtools::effect_plot(my_fit, pred = prs_scaled, interval = FALSE, partial.residuals = TRUE)  + 
    ggtitle(treat)
  
  sum_res = summary(my_fit)
  res[[treat]][["polygenicHR"]] = list(
    coefficients = sum_res$coefficients,
    varcor =  sum_res$varcor,
    residuals =  sum_res$residuals,
    call =  sum_res$call,
    AICtab =  sum_res$AICtab,
    r2=sum_res$r.squared,
    adjr2=sum_res$adj.r.squared,
    treatment = treat,
    PRS = "polygenicHR"
  )
  
  my_fit = lm(
    mean_outcome_pool ~ sex   + mean_min_prop_pool + APOE_sum_scaled ,
    data = totest
  )
  p_list[[paste0(treat,"_APOE")]] = jtools::effect_plot(my_fit, pred = APOE_sum_scaled, interval = FALSE, partial.residuals = TRUE)  + 
    ggtitle("")
  
  sum_res = summary(my_fit)
  res[[treat]][["APOE"]] = list(
    coefficients = sum_res$coefficients,
    varcor =  sum_res$varcor,
    residuals =  sum_res$residuals,
    call =  sum_res$call,
    AICtab =  sum_res$AICtab,
    r2=sum_res$r.squared,
    adjr2=sum_res$adj.r.squared,
    treatment = treat,
    PRS = "APOE"
  )
  
  
}

results_full = list()
results_polygenic = list()
results_APOE = list()
for(treat in c("untreated","LPS","IFN")){
  results_full[[treat]] = data.frame(treatment = NA,p = NA,coef=NA,se=NA,r2=NA, adjr2 = NA,PRS = NA, call = NA)
  results_full[[treat]]$treatment = treat
  results_full[[treat]]$p= res[[treat]][["full_PRS"]]$coefficients["full_PRS", "Pr(>|t|)"] # pval
  results_full[[treat]]$coef = res[[treat]][["full_PRS"]]$coefficients["full_PRS", "Estimate"] # coefficient
  results_full[[treat]]$se = res[[treat]][["full_PRS"]]$coefficients["full_PRS", "Std. Error"] # std error
  results_full[[treat]]$r2 = res[[treat]][["full_PRS"]]$r2 # r squared
  results_full[[treat]]$adjr2 = res[[treat]][["full_PRS"]]$adjr2 # adjusted r squared for multiple coefficients
  results_full[[treat]]$PRS = "full_PRS"
  results_full[[treat]]$call = as.character(res[[treat]][["full_PRS"]]$call)[2]
  results_polygenic[[treat]] = data.frame(treatment = NA,p = NA,coef=NA,se=NA,r2=NA, adjr2 = NA, PRS = NA,  call = NA)
  results_polygenic[[treat]]$treatment = treat
  results_polygenic[[treat]]$p= res[[treat]][["polygenicHR"]]$coefficients["prs_scaled", "Pr(>|t|)"] # pval
  results_polygenic[[treat]]$coef = res[[treat]][["polygenicHR"]]$coefficients["prs_scaled", "Estimate"] # coefficient
  results_polygenic[[treat]]$se = res[[treat]][["polygenicHR"]]$coefficients["prs_scaled", "Std. Error"] # std error
  results_polygenic[[treat]]$r2 = res[[treat]][["polygenicHR"]]$r2 # r squared
  results_polygenic[[treat]]$adjr2 = res[[treat]][["polygenicHR"]]$adjr2 
  results_polygenic[[treat]]$PRS = "polygenicHR"
  results_polygenic[[treat]]$call = as.character(res[[treat]][["polygenicHR"]]$call)[2]
  results_APOE[[treat]] = data.frame(treatment = NA,p = NA,coef=NA,se=NA,r2=NA, adjr2 = NA, PRS = NA,  call = NA)
  results_APOE[[treat]]$treatment = treat
  results_APOE[[treat]]$p= res[[treat]][["APOE"]]$coefficients["APOE_sum_scaled", "Pr(>|t|)"] # pval
  results_APOE[[treat]]$coef = res[[treat]][["APOE"]]$coefficients["APOE_sum_scaled", "Estimate"] # coefficient
  results_APOE[[treat]]$se = res[[treat]][["APOE"]]$coefficients["APOE_sum_scaled", "Std. Error"] # std error
  results_APOE[[treat]]$r2 = res[[treat]][["APOE"]]$r2 # r squared
  results_APOE[[treat]]$adjr2 = res[[treat]][["APOE"]]$adjr2 
  results_APOE[[treat]]$PRS = "APOE"
  results_APOE[[treat]]$call = as.character(res[[treat]][["APOE"]]$call)[2]
}

results_full = do.call("rbind",results_full)
results_polygenic = do.call("rbind",results_polygenic)
results_APOE = do.call("rbind",results_APOE)

results = rbind(results_full,results_polygenic,results_APOE)

write_csv(results, paste0(output_path,"full_PRS_phagocytosis_LM.csv"))

# checking some assumptions
# https://vasishth.github.io/Freq_CogSci/the-essentials-of-linear-modeling-theory.html#some-further-important-topics-in-linear-modeling
acf(residuals(res$untreated$full_PRS))
car::qqPlot(residuals(res$untreated$full_PRS))
acf(residuals(res$untreated$polygenicHR))
car::qqPlot(residuals(res$untreated$polygenicHR))


pdf(paste0(output_path,"PRS_phagocytosis_LM_fit_partial_residuals.pdf"),width = 10,height = 10)
patchwork::wrap_plots(p_list,ncol = 3,nrow = 3)
dev.off()
