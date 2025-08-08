# examine proliferation SKATO tests
library(patchwork)
library(tidyverse)
library(pROC)
library(ggpubr)
library(jtools)
library(lmerTest)
source("./functions.R")
directory = "../../data/results/2.2.rare_vars_vs_prolif/"
dir.create(directory, recursive = TRUE)

options(future.globals.maxSize = 240000 * 1024^2) # 240Gb

comparisons = c("line_prop_changes_old_vs_young_premac","line_prop_changes_microglia_premac","line_prop_changes_premac_iPSC")

for(mut_type in c("del","ptv","syn")){
  res = list()
  
for(comp in comparisons){
  
    files = list.files(paste0("../../data/results/2.2.rare_vars_vs_prolif/",comp,"/", mut_type),
                       pattern = paste0("[a-zA-Z0-9]","_",mut_type,"_burden_lm_res.rds"))
    files = str_sort(files) 
    genes = str_split_i(files,pattern = "_",i = 1) #they may not be contiguous
    res[[comp]]  = lapply(paste0("../../data/results/2.2.rare_vars_vs_prolif/",
                                       comp,"/", mut_type,"/",files), function(x){
      read_rds(x)
    })
    names(res[[comp]] ) = genes

    
    res_pvals = purrr::map(res[[comp]] , ~extract_burden_results(.))
    
    res_pvals = do.call("rbind",res_pvals)
    
    res[[comp]]  = res_pvals %>%
      dplyr::relocate(gene_name)
    
  
    res[[comp]]$comparison = comp

    res[[comp]]$p_Bonf = p.adjust(res[[comp]]$p_val,method = "bonferroni")
    
}
  res = do.call("rbind",res)

  # bonferroni correction across tests for genes that were actually tested (enough vars present)
  
  
  write_csv(res,paste0(directory,mut_type,"_burden_scaled_centered_prop_pvals.csv"))
}



for (mut_type in c("del","ptv","syn")){
  res = read_csv(paste0(directory,mut_type,"_burden_scaled_centered_prop_pvals.csv"))
  
  
  res= res%>%
    dplyr::arrange(p_Bonf)
  
  table(res$p_Bonf<0.05,res$comparison)
  
  
  ### plot all significant genes
  sign = res %>%
    dplyr::filter(p_Bonf < 0.05)
  
  for(comp in unique(sign$comparison)){
    pdf(paste0(directory,comp,"/",mut_type,"_burden_vars_vs_proliferation.pdf"),
        width = 3, height = 3)
    for(gene in unique(sign$gene_name[sign$comparison == comp])){
      toplot = sign %>%
        dplyr::filter(gene_name== gene & comparison %in% comp) 
      plots = read_rds(paste0("../../data/results/2.2.rare_vars_vs_prolif/",
                              toplot$comparison,"/",mut_type,"/",toplot$gene_name,"_",mut_type,"_partial_residual_LMM_plots.rds"))
      
      
      plot(plots[[1]])
      
      
    }
    
    dev.off()
  }
  # any in Sam's list?
  sam_phago = read_tsv("../../../CRISPR/OTAR2065_phagocytosis_CRISPR/data/2024_07_sam_screen_results.tsv") %>%
    dplyr::arrange(FDR)
  
  table(sam_phago$id %in% res$gene_name) 
  
  sam_phago = sam_phago %>%
    dplyr::rename(gene_name = id) %>%
    dplyr::left_join(res[res$comparison=="untreated",]) %>%
    dplyr::filter(!is.na(p_Bonf)) %>%
    dplyr::arrange(p_Bonf)
  
  # compare with burden test directionality
  write_csv(sam_phago,paste0(directory,mut_type,"_untreated_sam_CRISPR_significance.csv") )
}


###### now per treatment for microglia vs premac


for(mut_type in c("del","ptv","syn")){
  res = list()
  
  for(treat in c("untreated","IFNg","LPS")){

    files = list.files(paste0("../../data/results/2.3.rare_vars_vs_prolif_per_treatment/line_prop_changes_microglia_premac/", mut_type,"/",treat),
                       pattern = paste0("[a-zA-Z0-9]","_",mut_type,"_burden_lm_res.rds"))
    files = str_sort(files) 
    genes = str_split_i(files,pattern = "_",i = 1) #they may not be contiguous
    res[[treat]]  = lapply(paste0("../../data/results/2.3.rare_vars_vs_prolif_per_treatment/line_prop_changes_microglia_premac/", 
                                 mut_type,"/",treat,"/",files), function(x){
                                   read_rds(x)
                                 })
    names(res[[treat]] ) = genes
    
    
    res_pvals = purrr::map(res[[treat]] , ~extract_burden_results(.))
    
    res_pvals = do.call("rbind",res_pvals)
    
    res[[treat]]  = res_pvals %>%
      dplyr::relocate(gene_name)
    
    
    res[[treat]]$comparison = "line_prop_changes_microglia_premac"
    res[[treat]]$treatment = treat
    
    res[[treat]]$p_Bonf = p.adjust(res[[treat]]$p_val,method = "bonferroni")
    
  }
  res = do.call("rbind",res)
  
  # bonferroni correction across tests for genes that were actually tested (enough vars present)
  write_csv(res,paste0(directory,mut_type,"_burden_scaled_centered_prop_pvals.csv"))
  
  }

### AUROC / log-odds aggregated deleterious burdens #####
# related to Daianna's plots
# deleterious burdens aggregated by the 
 load("../../../OTAR2065_Daianna/output_data/04_Burden_test_proliferation/01_Burden_tests/all_Del_Burdens_per_line.Rdata")

load("../../../OTAR2065_Daianna/output_data/04_Burden_test_proliferation/01_Burden_tests/line_prop_changes_premac_vs_iPSC.Rdata")
load("../../../OTAR2065_Daianna/output_data/04_Burden_test_proliferation/01_Burden_tests/line_prop_changes_old_vs_young_premac.Rdata")
load("../../../OTAR2065_Daianna/output_data/04_Burden_test_proliferation/01_Burden_tests/line_prop_changes_microglia_vs_premac.Rdata")

all_Del_Burdens_per_line = all_Del_Burdens_per_line %>%
  tibble() %>%
  dplyr::select(line,contains("pos"),contains("neg"))

# premac vs iPSC - negative effects
test = line_prop_changes_premac_vs_iPSC %>% 
  dplyr::right_join(all_Del_Burdens_per_line[c("line","Burden_del_premac_vs_iPSC_burden_neg")], 
                    by = "line") %>%
  drop_na(scaled_log_fraction) %>%
  dplyr::group_by(line, pool) %>%
  # summarising per line and gene across all pools
  dplyr::reframe(scaled_log_fraction_pool = mean(scaled_log_fraction, na.rm = TRUE),
                 prop_unadjusted_min_value_mean = mean(prop_unadjusted_min_value,na.rm = TRUE),
                 sex = sex, genotypePC1= PC1, genotypePC2 =PC2,
                 rare_burden = as.numeric(Burden_del_premac_vs_iPSC_burden_neg))
                 
my_fit=  lmerTest::lmer(scaled_log_fraction_pool ~ sex + prop_unadjusted_min_value_mean +
                          genotypePC1 + genotypePC2 + 
                          rare_burden + 
                          (1 | pool),
                        data = test, 
                        control = lmerControl(optCtrl=list(maxfun=5000) ))


sum_res = summary(my_fit)

# getting covariate corrected scaled data
corrected_data = jtools::make_predictions(my_fit, 
                                                                pred = "rare_burden", 
                                                                partial.residuals = TRUE)$data
corrected_data = corrected_data %>%
  dplyr::select(scaled_log_fraction_pool,prop_unadjusted_min_value_mean,rare_burden) %>%
  dplyr::mutate(threshold_low = mean(scaled_log_fraction_pool) - 2 * sd(scaled_log_fraction_pool),
                threshold_high = mean(scaled_log_fraction_pool) + 2 * sd(scaled_log_fraction_pool),
                rare_burden = as.numeric(rare_burden)) %>%
  dplyr::mutate(dropout = ifelse(scaled_log_fraction_pool <= threshold_low, 1, 0),
                takeover = ifelse(scaled_log_fraction_pool >= threshold_high, 1, 0)) %>%
  dplyr::mutate(recoded_predictor_neg = ifelse(rare_burden>0,1,0)) %>%
  drop_na(scaled_log_fraction_pool)

# normal
hist(corrected_data$scaled_log_fraction_pool)

roc_premac_ipsc_neg <- roc(corrected_data$dropout, 
               corrected_data$rare_burden)
auc(roc_premac_ipsc_neg)  # prints AUROC value
# 0.66
plot(roc_premac_ipsc_neg, main = "AUROC - Aggregate Deleterious Variant Score",
     col = "darkgreen", lwd = 3)

### microglia vs premac
# negative direction

test = line_prop_changes_microglia_vs_premac %>% 
  dplyr::right_join(all_Del_Burdens_per_line[c("line","Burden_del_microglia_vs_premac_burden_neg")], 
                    by = "line") %>%
  drop_na(scaled_log_fraction) %>%
  dplyr::group_by(line, pool, treatment) %>%
  # summarising per line and gene across all pools
  dplyr::reframe(scaled_log_fraction_pool = mean(scaled_log_fraction, na.rm = TRUE),
                 prop_unadjusted_min_value_mean = mean(prop_unadjusted_min_value,na.rm = TRUE),
                 sex = sex, genotypePC1= PC1, genotypePC2 =PC2,
                 rare_burden = as.numeric(Burden_del_microglia_vs_premac_burden_neg),
                 treatment = treatment)

my_fit=  lmerTest::lmer(scaled_log_fraction_pool ~ sex + prop_unadjusted_min_value_mean +
                          genotypePC1 + genotypePC2 + 
                          rare_burden + 
                          (1 | pool) +
                          (1|treatment),
                        data = test, 
                        control = lmerControl(optCtrl=list(maxfun=5000) ))


sum_res = summary(my_fit)

corrected_data = jtools::make_predictions(my_fit, 
                                          pred = "rare_burden", 
                                          partial.residuals = TRUE)$data
corrected_data = corrected_data %>%
  dplyr::select(scaled_log_fraction_pool,prop_unadjusted_min_value_mean,rare_burden) %>%
  dplyr::mutate(threshold_low = mean(scaled_log_fraction_pool) - 2 * sd(scaled_log_fraction_pool),
                threshold_high = mean(scaled_log_fraction_pool) + 2 * sd(scaled_log_fraction_pool),
                rare_burden = as.numeric(rare_burden)) %>%
  dplyr::mutate(dropout = ifelse(scaled_log_fraction_pool <= threshold_low, 1, 0),
                takeover = ifelse(scaled_log_fraction_pool >= threshold_high, 1, 0)) %>%
  dplyr::mutate(recoded_predictor_neg = ifelse(rare_burden>0,1,0)) %>%
  drop_na(scaled_log_fraction_pool)

# normal
hist(corrected_data$scaled_log_fraction_pool)

corrected_data %>%
  ggplot(aes(y = scaled_log_fraction_pool, x = rare_burden, group = rare_burden)) +
  geom_jitter(width = 0.1, alpha = 0.7, size = 0.9) +
  geom_violin(fill = NA, color = "black", linewidth = 0.4) +
  stat_smooth(geom="line", alpha = 1, size = 1, method = lm, aes(group=1), color='black') +
  theme_classic()


roc_microglia_premac_neg <- roc(response = corrected_data$dropout, 
               predictor = corrected_data$rare_burden)
auc(roc_microglia_premac_neg)  # prints AUROC value
# 0.60
plot(roc_microglia_premac_neg, main = "AUROC - Aggregate Deleterious Variant Score",
     col = "darkgreen", lwd = 3)

## positive
test = line_prop_changes_microglia_vs_premac %>% 
  dplyr::right_join(all_Del_Burdens_per_line[c("line","Burden_del_microglia_vs_premac_burden_pos")], 
                    by = "line") %>%
  drop_na(scaled_log_fraction) %>%
  dplyr::group_by(line, pool, treatment) %>%
  # summarising per line and gene across all pools
  dplyr::reframe(scaled_log_fraction_pool = mean(scaled_log_fraction, na.rm = TRUE),
                 prop_unadjusted_min_value_mean = mean(prop_unadjusted_min_value,na.rm = TRUE),
                 sex = sex, genotypePC1= PC1, genotypePC2 =PC2,
                 rare_burden = as.numeric(Burden_del_microglia_vs_premac_burden_pos),
                 treatment = treatment)

my_fit=  lmerTest::lmer(scaled_log_fraction_pool ~ sex + prop_unadjusted_min_value_mean +
                          genotypePC1 + genotypePC2 + 
                          rare_burden + 
                          (1 | pool) +
                          (1|treatment),
                        data = test, 
                        control = lmerControl(optCtrl=list(maxfun=5000) ))


sum_res = summary(my_fit)

corrected_data = jtools::make_predictions(my_fit, 
                                          pred = "rare_burden", 
                                          partial.residuals = TRUE)$data
corrected_data = corrected_data %>%
  dplyr::select(scaled_log_fraction_pool,prop_unadjusted_min_value_mean,rare_burden) %>%
  dplyr::mutate(threshold_low = mean(scaled_log_fraction_pool) - 2 * sd(scaled_log_fraction_pool),
                threshold_high = mean(scaled_log_fraction_pool) + 2 * sd(scaled_log_fraction_pool),
                rare_burden = as.numeric(rare_burden)) %>%
  dplyr::mutate(dropout = ifelse(scaled_log_fraction_pool <= threshold_low, 1, 0),
                takeover = ifelse(scaled_log_fraction_pool >= threshold_high, 1, 0)) %>%
  dplyr::mutate(recoded_predictor_pos = ifelse(rare_burden>10,1,0)) %>%
  drop_na(scaled_log_fraction_pool)

# normal
hist(corrected_data$scaled_log_fraction_pool)

corrected_data %>%
  ggplot(aes(y = scaled_log_fraction_pool, 
             x = rare_burden, group = rare_burden)) +
  geom_jitter(width = 0.1, alpha = 0.7, size = 0.9) +
  geom_violin(fill = NA, color = "black", linewidth = 0.4) +
  stat_smooth(geom="line", alpha = 1, size = 1, method = lm, aes(group=1), color='black') +
  theme_classic()


roc_microglia_premac_pos <- roc(response = corrected_data$takeover, 
               predictor = corrected_data$rare_burden)
auc(roc_microglia_premac_pos)  # prints AUROC value
# 0.79
plot(roc_microglia_premac_pos, main = "AUROC - Aggregate Deleterious Variant Score",
     col = "darkgreen", lwd = 3)
roc_microglia_premac_pos_threshold <- roc(response = corrected_data$takeover, 
               predictor = corrected_data$recoded_predictor_pos)
auc(roc_microglia_premac_pos_threshold)  # prints AUROC value
# 0.78 if chosen 10 aggregate burden as threshold
ggroc(roc_microglia_premac_pos,
     col = "darkgreen", lwd = 1.5) +
  geom_abline(intercept = 1) + 
  ggtitle("AUROC - Aggregate Deleterious Variant Score") +
  ggpubr::theme_pubr()


roc_combined = ggroc(list(roc_premac_ipsc_neg, 
                          roc_microglia_premac_neg,
                          roc_microglia_premac_pos)) +
  geom_abline(intercept = 1) + 
  ggpubr::theme_pubr() +
  ggtitle("") + 
  scale_color_manual(labels =  c('1' = 'iPSC->preMAC: negative', 
                                  '2' = 'preMAC->microglia: negative', 
                                  '3' = 'preMAC->microglia: positive'),
                     values = c('1'= "#44CFCB",
                                '2' = "#2A4494",
                                '3' = "#BA274A")) +
  geom_text(aes(x = 0.75, y = 0.60, 
                label = as.character(round(auc(roc_premac_ipsc_neg),2)), color = "#44CFCB")) +
  geom_text(aes(x = 0.75, y = 0.42, 
                label =as.character(round(auc(roc_microglia_premac_neg),2))), color = "#2A4494") +
  geom_text(aes(x = 0.90, y = 0.74, 
                label =as.character(round(auc(roc_microglia_premac_pos),2))), color = "#BA274A") + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 8),
        # legend.key.width = unit(0.4, "cm"),
        legend.position = "inside", legend.position.inside =  c(0.65, 0.20),   # Place legend inside the plot at (x, y) coordinates
        legend.direction = "vertical")   # Set legend items to be arranged vertically)

pdf(paste0(directory,"/aggregate_del_burden_vars_vs_proliferation_AUROC.pdf"),
    width = 4, height = 4.1)

roc_combined

dev.off()

sessionInfo()

# R version 4.4.1 (2024-06-14)
# Platform: x86_64-pc-linux-gnu
# Running under: Ubuntu 22.04.4 LTS
# 
# Matrix products: default
# BLAS:   /software/hgi/installs/softpack/rstudio/.spack/linux-ubuntu22.04-x86_64_v3/gcc-11.4.0/r-4.4.1-zf3d5qbxgbiqsk4ddke3fl6uluwcbqcu/rlib/R/lib/libRblas.so 
# LAPACK: /software/hgi/installs/softpack/rstudio/.spack/linux-ubuntu22.04-x86_64_v3/gcc-11.4.0/r-4.4.1-zf3d5qbxgbiqsk4ddke3fl6uluwcbqcu/rlib/R/lib/libRlapack.so;  LAPACK version 3.12.0
# 
# locale:
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8   
# [6] LC_MESSAGES=en_GB.UTF-8    LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Belfast
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggpubr_0.6.0    pROC_1.18.5     lmerTest_3.1-3  lme4_1.1-35.2   Matrix_1.6-5    jtools_2.2.2    lubridate_1.9.3 forcats_1.0.0  
# [9] stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2     readr_2.1.5     tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0
# [17] patchwork_1.2.0
# 
# loaded via a namespace (and not attached):
#   [1] gtable_0.3.4        rstatix_0.7.2       lattice_0.22-6      tzdb_0.4.0          numDeriv_2016.8-1.1 vctrs_0.6.5        
# [7] tools_4.4.1         generics_0.1.3      fansi_1.0.6         pkgconfig_2.0.3     lifecycle_1.0.4     compiler_4.4.1     
# [13] farver_2.1.1        munsell_0.5.1       ggforce_0.4.2       carData_3.0-5       pillar_1.9.0        car_3.1-2          
# [19] nloptr_2.0.3        crayon_1.5.2        MASS_7.3-65         boot_1.3-30         abind_1.4-5         nlme_3.1-164       
# [25] tidyselect_1.2.1    digest_0.6.35       stringi_1.8.3       pander_0.6.5        labeling_0.4.3      splines_4.4.1      
# [31] polyclip_1.10-6     grid_4.4.1          colorspace_2.1-0    cli_3.6.2           magrittr_2.0.3      utf8_1.2.4         
# [37] broom_1.0.5         withr_3.0.0         scales_1.3.0        backports_1.4.1     timechange_0.3.0    ggsignif_0.6.4     
# [43] hms_1.1.3           mgcv_1.9-1          rlang_1.1.3         Rcpp_1.0.12         glue_1.7.0          tweenr_2.0.3       
# [49] pkgload_1.3.4       rstudioapi_0.16.0   minqa_1.2.6         R6_2.5.1            plyr_1.8.9    