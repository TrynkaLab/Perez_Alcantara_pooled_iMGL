# rare variants vs migration

library(patchwork)
library(tidyverse)
library(future)
library(lmerTest)
library(variancePartition)
library(jtools)

skat_results_path = "../../../data/results/migration/5.1.check_SKAT_WES_vs_migration/deleterious_SKATO_scaled_centered_prop_pvals.csv"
directory = "../../../data/results/migration/5.2.rare_vars_vs_migration/"
dir.create(directory, recursive = T)

options(future.globals.maxSize = 240000 * 1024^2) # 240Gb

line_prop_changes_path = "../../../data/results/migration/1.check_line_proportions/line_prop_changes_per_well.csv"
genotype_pc_path = "../../../../OTAR2065_sc_eQTL/data/genotype/plink_genotypes/all_pools.all_donors.genotype.MAF05.eigenvec"
scaled_proliferation_path = "../../../../OTAR2065_WGS_iPSC_premac_micro/data/results/1.2.scale_proliferation/line_prop_changes_microglia_premac.csv"


proliferation_w_info = read_rds("../../../../OTAR2065_differentiation_efficiency/data/results/2.efficiency/proliferation_w_info.rds") %>%
  # take non-proliferating
  dplyr::filter(proliferation_status == "Not_proliferating") %>%
  # dicotomising first burden: 0 == no, else == yes
  dplyr::mutate(rare_burden_dichotomy = factor(case_when(rare_burden == 0 ~ "no",
                                                         .default = "yes"))) %>%
  dplyr::mutate(rare_burden_dichotomy = relevel(rare_burden_dichotomy,ref = c("no"))) %>%  # setting no as baseline
  dplyr::select(line, treatment, sex, gene, rare_burden, rare_mutation_type, rare_burden_dichotomy) %>%
    dplyr::distinct()


gc()


scaled_proliferation = readr::read_csv(scaled_proliferation_path) %>%
  dplyr::mutate(treatment = case_when(treatment == "IFNg" ~ "IFN",
                                      .default = treatment)) %>%
  dplyr::filter(!scaled_log_fraction %in% c(NA,Inf,-Inf)) %>%
  dplyr::group_by(line,pool,treatment) %>%
  dplyr::summarise(scaled_proliferation = mean(scaled_log_fraction,na.rm = TRUE)) %>%
  dplyr::ungroup()

line_prop_changes = read.csv(line_prop_changes_path) 

message("There are ", length(unique(line_prop_changes$line)), " lines in the analysis")

line_prop_changes = line_prop_changes %>%
  # dplyr::filter(prop_unadjusted_max_value > 0.005) %>%
  dplyr::filter(!scaled_log_fraction %in% c(NA,Inf,-Inf)) %>%
  dplyr::select(scaled_log_fraction,replicate,line,condition,treatment,pool,sex,prop_unadjusted_min_value) %>%
  distinct() 

message("There are ", length(unique(line_prop_changes$line)), " lines in the analysis after Na and Inf filters") # 216

#### remove genes we're not interested in (not significant in SKAT)
skat_res = readr::read_csv(skat_results_path) %>%
  dplyr::filter(!is.na(resampling_pval)) %>%
  dplyr::filter(resampling_pval < 0.05) 

table(skat_res$treatment)
# IFN       LPS untreated 
# 161       245       201 

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
  ) %>%
  dplyr::mutate(line = case_when(line=="Arlene-003" ~ "Arlene",
                                 line=="Cindy-005" ~ "Cindy",
                                 line=="Dexter-006" ~ "Dexter",
                                 line=="Fiona-010" ~ "Fiona",
                                 line=="Gilma-009" ~ "Gilma",
                                 line=="Hector-011" ~ "Hector",
                                 line=="Imani-012" ~ "Imani",
                                 line=="Javier-013" ~ "Javier",
                                 line=="Keoni-014" ~ "Keoni",
                                 line=="Olaf-018" ~ "Olaf",
                                 line=="Bertha-004" ~ "Bertha",
                                 line=="Mindy-016" ~ "Mindy",
                                 line=="Qiana-022" ~ "Qiana",
                                 line=="Nestor-017" ~ "Nestor",
                                 .default = line))

line_prop_changes = line_prop_changes %>%
  dplyr::left_join(full_genotype_pcs)

# test only remaining lines after filters
# keep in mind I don't have info for the IPMAR donors

proliferation_w_info = proliferation_w_info %>%
  dplyr::filter(line %in% unique(line_prop_changes$line)) %>%
  # filter to significantresults in SKAt
  dplyr::filter(gene %in% skat_res$gene_name)

gc()

message("There are ", length(unique(proliferation_w_info$line)), " lines in the analysis after missing line filters") # 176

### adding other info
proliferation_w_info = proliferation_w_info %>%
  dplyr::left_join(line_prop_changes) %>%
  dplyr::left_join(scaled_proliferation) %>%
  distinct()

# sort genes by number of samples in the burdened category
genes = proliferation_w_info %>% 
  dplyr::filter(rare_burden_dichotomy == "yes") %>%
  dplyr::group_by(rare_mutation_type,gene,treatment) %>% 
  dplyr::summarise(number=n()) %>% 
  dplyr::arrange(desc(number))  %>% 
  .$gene 

### regression
mut_types = unique(proliferation_w_info$rare_mutation_type)
for(mut_type in mut_types){
  message("Working on ", mut_type, " variants.")
  dir.create(paste0(directory, mut_type), recursive = TRUE)
  
  for(treat in c("untreated","LPS","IFN")){
    
    sub_table = proliferation_w_info %>%
      dplyr::filter(rare_mutation_type == mut_type  & treatment == treat) %>%
      dplyr::group_by(line, gene, pool) %>%
      # summarising per line and gene across all pools
      dplyr::reframe(scaled_log_fraction_pool = mean(scaled_log_fraction, na.rm = TRUE),
                     prop_unadjusted_min_value_mean = mean(prop_unadjusted_min_value,na.rm = TRUE),
                     sex = sex, genotypePC1= genotypePC1, genotypePC2 =genotypePC2,
                     rare_burden = rare_burden, rare_burden_dichotomy= rare_burden_dichotomy,
                     scaled_proliferation = scaled_proliferation) %>%
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
    
    # and remove genes that are not significant in specific skat treatment and
    # that have only one level of rare burden
    genes_to_keep = sub_table %>%
      dplyr::filter(gene %in% skat_res[skat_res$treatment == treat,"gene_name"]$gene_name) %>%
      dplyr::select(gene,rare_burden) %>%
      dplyr::distinct() %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(n_burden = n()) %>%
      dplyr::ungroup() %>%
      dplyr::filter(n_burden > 1)
    
    sub_table = sub_table %>%
      dplyr::filter(gene %in% genes_to_keep$gene )
    
    length(unique(sub_table$gene)) 
    
    p_list = list()
    res = list()
    for(gene in unique(sub_table$gene)){
      
      
      
      
      if(!gene %in% sub_table$gene){
        next()
      }else{
        
        
        test = sub_table[sub_table$gene == gene, ]
        
        
        
        # rare burden as continuous variable
        
        my_fit=  lmerTest::lmer(scaled_log_fraction_pool ~ rare_burden + # effect of each unit of change in allele dosage in treatment vs baseline
                                  
                                  sex + prop_unadjusted_min_value_mean + scaled_proliferation +
                                  genotypePC1 + genotypePC2 + 
                                  
                                  (1 | pool),
                                data = test, 
                                control = lmerControl(optCtrl=list(maxfun=5000) ))
        
        
        sum_res = summary(my_fit)
        res[[paste0(treat,"_",gene,"_",mut_type)]] = list(
          coefficients = sum_res$coefficients,
          varcor =  sum_res$varcor,
          residuals =  sum_res$residuals,
          call =  sum_res$call,
          AICtab =  sum_res$AICtab,
          gene = gene)
        
        p_list[[paste0(treat,"_",gene,"_",mut_type)]] = jtools::effect_plot(my_fit, 
                                                                            pred = rare_burden, 
                                                                            interval = FALSE, partial.residuals = TRUE,
                                                                            jitter=0.1,plot.points = TRUE)  + 
          ggtitle(paste0(treat," - ",gene),subtitle = paste0(mut_type," beta: ",
                                                             round(res[[paste0(treat,"_",gene,"_",mut_type)]]$coefficients["rare_burden","Estimate"],2))) + 
          xlab("Burden") + ylab("Scaled log fraction") +
          scale_x_continuous(breaks = sort(unique(test$rare_burden))) 
        
        
        
        
      }
      
    }
    # save these as rds
    saveRDS(object = res,paste0(directory,mut_type,"/",mut_type,"_",treat,"_burden_lm_res.rds"))
    if(mut_type != "syn"){
      saveRDS(object = p_list,paste0(directory,mut_type,"/",treat,"partial_residual_LMM_plots.rds"))
    }
    
    rm(res)
    gc()
  }
  
  gc()
}

