# rare variants vs proliferation

library(patchwork)
library(tidyverse)
library(future)
library(lmerTest)
library(variancePartition)
library(jtools)
library(ggforce)
library(ggpubr)

directory = "../../data/results/2.2.rare_vars_vs_prolif/"
dir.create(directory, recursive = T)

options(future.globals.maxSize = 240000 * 1024^2) # 240Gb


skat_results_path = "../../data/results/2.1.check_SKAT_WES/deleterious_SKATO_scaled_centered_prop_pvals.csv"
genotype_pc_path = "../../../OTAR2065_sc_eQTL/data/genotype/plink_genotypes/all_pools.all_donors.genotype.MAF05.eigenvec"

for(contrast in c("line_prop_changes_old_vs_young_premac","line_prop_changes_microglia_premac","line_prop_changes_premac_iPSC")){
line_prop_changes_path = paste0("../../data/results/1.2.scale_proliferation/",contrast,".csv")

proliferation_w_info = read_rds("../../../OTAR2065_differentiation_efficiency/data/results/2.efficiency/proliferation_w_info.rds") %>%
  # take non-proliferating
  dplyr::filter(proliferation_status == "Not_proliferating") %>%
  # dicotomising first burden: 0 == no, else == yes
  dplyr::mutate(rare_burden_dichotomy = factor(case_when(rare_burden == 0 ~ "no",
                                                         .default = "yes"))) %>%
  dplyr::mutate(rare_burden_dichotomy = relevel(rare_burden_dichotomy,ref = c("no"))) %>%  # setting no as baseline
  dplyr::select(line, sex, gene, rare_burden, rare_mutation_type, rare_burden_dichotomy) %>%
    dplyr::distinct()

### move the creation of this file from the "efficiency" (old code) to here
write_csv(proliferation_w_info,paste0(directory,"burden_info.csv"))

gc()

line_prop_changes = read_csv(line_prop_changes_path)

message("There are ", length(unique(line_prop_changes$line)), " lines in the analysis") # 261

skat_res = read_csv(skat_results_path) %>%
  dplyr::filter(comparison %in% str_remove(basename(line_prop_changes_path),".csv") & !is.na(resampling_pval)) %>%
  # only retain significant for testing
  dplyr::filter(resampling_pval < 0.05) 


line_prop_changes = line_prop_changes %>%
  # dplyr::filter(prop_unadjusted_max_value > 0.005) %>%
  dplyr::filter(!scaled_log_fraction %in% c(NA,Inf,-Inf)) %>%
  dplyr::select(any_of(c("scaled_log_fraction","line","donor","differentiation","treatment","pool","prop_unadjusted_min_value","type"))) %>%
  distinct() 

message("There are ", length(unique(line_prop_changes$line)), " lines in the analysis after Na and Inf filters") # 220 for old vs young premac
# 241 for microglia vs premac

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

line_prop_changes = line_prop_changes %>%
  dplyr::left_join(full_genotype_pcs)

length(unique(line_prop_changes$line))

# test only remaining lines after filters
# keep in mind I don't have info for the IPMAR donors

proliferation_w_info = proliferation_w_info %>%
  dplyr::filter(line %in% unique(line_prop_changes$line)) %>%
  # filter to significantresults in SKAt
  dplyr::filter(gene %in% skat_res$gene_name)

gc()


message("There are ", length(unique(proliferation_w_info$line)), 
        " lines in the analysis after missing line filters") # 178 old vs young premac
# 197 microglia vs premac
# 194 preMac vs iPSC
### adding other info
proliferation_w_info = proliferation_w_info %>%
  dplyr::left_join(line_prop_changes) %>%
  distinct()

# sort genes by number of samples in the burdened category
genes = proliferation_w_info %>% 
  dplyr::filter(rare_burden_dichotomy == "yes") %>%
  dplyr::group_by(rare_mutation_type,gene) %>% 
  dplyr::summarise(number=n()) %>% 
  dplyr::arrange(desc(number))  %>% 
  .$gene 


### fix this
# info = proliferation_w_info %>%
#   dplyr::select(treatment,scaled_log_fraction,replicate, line,pool,sex,prop_unadjusted_min_value, genotypePC1, genotypePC2,genotypePC3, genotypePC4, genotypePC5,
#                 rare_burden, rare_burden_dichotomy) %>%
#   as.data.frame()
# 
# 
# 
# form = ~ treatment + pool + replicate  + 
#   prop_unadjusted_min_value + scaled_log_fraction + 
#   genotypePC1 + genotypePC2 + genotypePC3 + genotypePC4 + genotypePC5 + 
#   rare_burden + rare_burden_dichotomy
# 
# # Compute Canonical Correlation Analysis (CCA)
# # between all pairs of variables
# # returns absolute correlation value
# C = canCorPairs(form, info)

# Plot correlation matrix
# between all pairs of variables


### regression
mut_types = unique(proliferation_w_info$rare_mutation_type)
for(mut_type in mut_types){
  message("Working on ", mut_type, " variants.")
  dir.create(paste0(directory,unique(skat_res$comparison),"/", mut_type), recursive = TRUE)
  
  if("treatment" %in% colnames(proliferation_w_info)){
    sub_table = proliferation_w_info %>%
      dplyr::filter(rare_mutation_type == mut_type) %>%
      dplyr::group_by(line, gene, pool, treatment) %>%
      # summarising per line and gene across all pools
      dplyr::reframe(scaled_log_fraction_pool = mean(scaled_log_fraction, na.rm = TRUE),
                     prop_unadjusted_min_value_mean = mean(prop_unadjusted_min_value,na.rm = TRUE),
                     sex = sex, genotypePC1= genotypePC1, genotypePC2 =genotypePC2,
                     rare_burden = rare_burden, rare_burden_dichotomy= rare_burden_dichotomy) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(gene_cat = paste(gene,rare_burden,sep = "_")) %>%
      dplyr::distinct() %>%
      dplyr::filter(!is.na(scaled_log_fraction_pool))
    
    length(unique(sub_table$gene))
    # remove categories of rare_burden within gene that have fewer than 3 donors each
    
    genes_cat_to_keep = sub_table %>%
      dplyr::select(gene,line,rare_burden) %>%
      dplyr::distinct() %>%
      dplyr::group_by(gene, rare_burden) %>%
      dplyr::summarise(n_lines = n()) %>%
      dplyr::ungroup() %>%
      dplyr::filter(n_lines >= 3) %>%
      dplyr::mutate(gene_cat = paste(gene,rare_burden,sep = "_")) %>%
      dplyr::select(gene_cat)
    
    sub_table = sub_table %>%
      dplyr::filter(gene_cat %in% genes_cat_to_keep$gene_cat )
    
    # and remove genes that have only one level of rare burden
    genes_to_keep = sub_table %>%
      dplyr::select(gene,rare_burden) %>%
      dplyr::distinct() %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(n_burden = n()) %>%
      dplyr::ungroup() %>%
      dplyr::filter(n_burden > 1) 
    
    sub_table = sub_table %>%
      dplyr::filter(gene %in% genes_to_keep$gene )
    
    length(unique(sub_table$gene))
    
    
    for(gene in unique(sub_table$gene)){
      res_dichotomy = list()
      res = list()
      p_list = list()
      
      
      test = sub_table[sub_table$gene == gene, ]
      
      
      
      # rare burden as continuous variable
      
      my_fit=  lmerTest::lmer(scaled_log_fraction_pool ~ sex + prop_unadjusted_min_value_mean +
                                genotypePC1 + genotypePC2 + 
                                rare_burden + 
                                 (1| treatment) + 
                                (1 | pool),
                              data = test, 
                              control = lmerControl(optCtrl=list(maxfun=5000) ))
      
      
      sum_res = summary(my_fit)
      res = list(
        coefficients = sum_res$coefficients,
        varcor =  sum_res$varcor,
        residuals =  sum_res$residuals,
        call =  sum_res$call,
        AICtab =  sum_res$AICtab,
        gene = gene)
      
      p_list[[paste0(gene,"_",mut_type)]]  = jtools::make_predictions(my_fit, 
                                                                      pred = "rare_burden", 
                                                                      partial.residuals = TRUE)$data %>%
        ggplot(., aes(x = rare_burden, y = scaled_log_fraction_pool, group = rare_burden)) + 
        geom_violin(fill = NA, color = "grey20", linewidth = 0.4) +
        ggforce::geom_sina(col = "grey40", alpha = 0.30, size = 0.5) +
        stat_smooth(geom="line", alpha = 1, size = 1, method = lm, aes(group=1), color='grey20') +
        scale_x_continuous(
          breaks = unique(test$rare_burden)) +
        ggpubr::theme_pubr() +
        labs(title = gene, 
             subtitle = paste("Beta:", signif(res$coefficients["rare_burden", "Estimate"], digits = 2), '    ', 
                              'p:', signif(res$coefficients["rare_burden", "Pr(>|t|)"], digits = 2) ), 
             x = "Burden", y = paste0("Differentiation efficiency")) 
      
      # save these as rds
      saveRDS(object = res,paste0(directory,unique(skat_res$comparison),"/",
                                  mut_type,"/",gene,"_",mut_type,"_burden_lm_res.rds"))
      saveRDS(object = p_list,paste0(directory,unique(skat_res$comparison),"/",
                                     mut_type,"/",gene,"_",mut_type,"_partial_residual_LMM_plots.rds"))
      
      
    }
    
    
  } else{
    sub_table = proliferation_w_info %>%
      dplyr::filter(rare_mutation_type == mut_type) %>%
      dplyr::group_by(line, gene, pool) %>%
      # summarising per line and gene across all pools
      dplyr::reframe(scaled_log_fraction_pool = mean(scaled_log_fraction),
                     prop_unadjusted_min_value_mean = mean(prop_unadjusted_min_value),
                     sex = sex, genotypePC1= genotypePC1, genotypePC2 =genotypePC2,
                     rare_burden = rare_burden, rare_burden_dichotomy= rare_burden_dichotomy) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(gene_cat = paste(gene,rare_burden,sep = "_")) %>%
      dplyr::distinct()
    
      length(unique(sub_table$gene))
    # remove categories of rare_burden within gene that have fewer than 3 donors each
    
    genes_cat_to_keep = sub_table %>%
      dplyr::select(gene,line,rare_burden) %>%
      dplyr::distinct() %>%
      dplyr::group_by(gene, rare_burden) %>%
      dplyr::summarise(n_lines = n()) %>%
      dplyr::ungroup() %>%
      dplyr::filter(n_lines >= 3) %>%
      dplyr::mutate(gene_cat = paste(gene,rare_burden,sep = "_")) %>%
      dplyr::select(gene_cat)
    
    sub_table = sub_table %>%
      dplyr::filter(gene_cat %in% genes_cat_to_keep$gene_cat )
    
    # and remove genes that have only one level of rare burden
    genes_to_keep = sub_table %>%
      dplyr::select(gene,rare_burden) %>%
      dplyr::distinct() %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(n_burden = n()) %>%
      dplyr::ungroup() %>%
      dplyr::filter(n_burden > 1) 
    
    sub_table = sub_table %>%
      dplyr::filter(gene %in% genes_to_keep$gene )
    
    length(unique(sub_table$gene))
    
    
    for(gene in unique(sub_table$gene)){
      res_dichotomy = list()
      res = list()
      p_list = list()
        
        
        test = sub_table[sub_table$gene == gene, ]
        
        
        
        # rare burden as continuous variable
        
        my_fit=  lmerTest::lmer(scaled_log_fraction_pool ~ sex + prop_unadjusted_min_value_mean +
                                  genotypePC1 + genotypePC2 + 
                                  rare_burden + 
                                  # (1|line) + 
                                  (1 | pool),
                                data = test, 
                                control = lmerControl(optCtrl=list(maxfun=5000) ))
        
        
        sum_res = summary(my_fit)
        res = list(
          coefficients = sum_res$coefficients,
          varcor =  sum_res$varcor,
          residuals =  sum_res$residuals,
          call =  sum_res$call,
          AICtab =  sum_res$AICtab,
          gene = gene)
        
        p_list[[paste0(gene,"_",mut_type)]]  = jtools::make_predictions(my_fit, 
                                 pred = "rare_burden", 
                                 partial.residuals = TRUE)$data %>%
          ggplot(., aes(x = rare_burden, y = scaled_log_fraction_pool, group = rare_burden)) + 
          geom_violin(fill = NA, color = "grey20", linewidth = 0.4) +
          ggforce::geom_sina(col = "grey40", alpha = 0.30, size = 0.5) +
          stat_smooth(geom="line", alpha = 1, size = 1, method = lm, aes(group=1), color='grey20') +
          scale_x_continuous(
                         breaks = unique(test$rare_burden)) +
          ggpubr::theme_pubr() +
          labs(title = gene, 
               subtitle = paste("Beta:", signif(res$coefficients["rare_burden", "Estimate"], digits = 2), '    ', 
                                'p:', signif(res$coefficients["rare_burden", "Pr(>|t|)"], digits = 2) ), 
               x = "Burden", y = paste0("Differentiation efficiency")) 
  
        
        
        # save these as rds
        saveRDS(object = res,paste0(directory,unique(skat_res$comparison),"/",
                                    mut_type,"/",gene,"_",mut_type,"_burden_lm_res.rds"))
          saveRDS(object = p_list,paste0(directory,unique(skat_res$comparison),"/",
                                         mut_type,"/",gene,"_",mut_type,"_partial_residual_LMM_plots.rds"))
        
        
  }
    
      
    }
  }
  
}

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
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
# [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8    LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Belfast
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggpubr_0.6.0             ggforce_0.4.2            jtools_2.2.2             variancePartition_1.32.5
# [5] BiocParallel_1.36.0      limma_3.58.1             lmerTest_3.1-3           lme4_1.1-35.2           
# [9] Matrix_1.6-5             future_1.33.2            lubridate_1.9.3          forcats_1.0.0           
# [13] stringr_1.5.1            dplyr_1.1.4              purrr_1.0.2              readr_2.1.5             
# [17] tidyr_1.3.1              tibble_3.2.1             ggplot2_3.5.1            tidyverse_2.0.0         
# [21] patchwork_1.2.0         
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.2.1    farver_2.1.1        bitops_1.0-7        tweenr_2.0.3        digest_0.6.35      
# [6] timechange_0.3.0    lifecycle_1.0.4     statmod_1.5.0       magrittr_2.0.3      compiler_4.4.1     
# [11] rlang_1.1.3         tools_4.4.1         utf8_1.2.4          ggsignif_0.6.4      plyr_1.8.9         
# [16] abind_1.4-5         KernSmooth_2.23-22  withr_3.0.0         numDeriv_2016.8-1.1 BiocGenerics_0.48.1
# [21] grid_4.4.1          polyclip_1.10-6     aod_1.3.3           fansi_1.0.6         caTools_1.18.2     
# [26] colorspace_2.1-0    globals_0.16.3      scales_1.3.0        gtools_3.9.5        iterators_1.0.14   
# [31] MASS_7.3-65         mvtnorm_1.2-4       cli_3.6.2           crayon_1.5.2        generics_0.1.3     
# [36] rstudioapi_0.16.0   reshape2_1.4.4      tzdb_0.4.0          minqa_1.2.6         pander_0.6.5       
# [41] splines_4.4.1       parallel_4.4.1      matrixStats_1.2.0   vctrs_0.6.5         boot_1.3-30        
# [46] carData_3.0-5       car_3.1-2           hms_1.1.3           pbkrtest_0.5.2      rstatix_0.7.2      
# [51] listenv_0.9.1       glue_1.7.0          parallelly_1.37.1   nloptr_2.0.3        codetools_0.2-20   
# [56] stringi_1.8.3       gtable_0.3.4        EnvStats_2.8.1      munsell_0.5.1       remaCor_0.0.18     
# [61] pillar_1.9.0        gplots_3.1.3.1      R6_2.5.1            Rdpack_2.6          lattice_0.22-6     
# [66] Biobase_2.62.0      rbibutils_2.2.16    backports_1.4.1     RhpcBLASctl_0.23-42 broom_1.0.5        
# [71] fANCOVA_0.6-1       corpcor_1.6.10      Rcpp_1.0.12         nlme_3.1-164        pkgconfig_2.0.3  
# 
# 
# 
