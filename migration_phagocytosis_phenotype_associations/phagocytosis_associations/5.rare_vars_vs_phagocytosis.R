# rare variants vs phagocytosis

library(patchwork)
library(tidyverse)
library(future)
library(lmerTest)
library(variancePartition)
library(jtools)

directory = "../../../data/results/phagocytosis/5.rare_vars_vs_phagocytosis/"
dir.create(directory, recursive = T)

options(future.globals.maxSize = 200000 * 1024^2) # 200Gb


line_prop_changes_path = "../../../data/results/phagocytosis/1.check_line_proportions/line_prop_changes_per_well.csv"
genotype_pc_path = "../../../../OTAR2065_sc_eQTL/data/genotype/plink_genotypes/all_pools.no_outliers.genotype.MAF05.eigenvec"
prs_path = "../../../../hipsci_genotype_processing/data/prs_ad_bellenguez/PRSice/hipsci_raw_polygenic_score_AD_Bellenguez_withAPOE.tsv"


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

prs = read_tsv(prs_path)  %>%
  dplyr::select(line, APOE_sum_scaled, prs_scaled, full_PRS)
line_prop_changes = read.csv(line_prop_changes_path) %>%
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

message("There are ", length(unique(line_prop_changes$line)), " lines in the analysis")

line_prop_changes = line_prop_changes %>%
  # dplyr::filter(prop_unadjusted_max_value > 0.005) %>%
  dplyr::filter(!scaled_log_fraction %in% c(NA,Inf,-Inf)) %>%
  dplyr::select(scaled_log_fraction,replicate,line,condition,treatment,pool,sex,prop_unadjusted_min_value) %>%
  distinct() 

message("There are ", length(unique(line_prop_changes$line)), " lines in the analysis after Na and Inf filters") # 202

### adding PRS info
line_prop_changes = line_prop_changes %>%
  dplyr::left_join(prs)

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
  dplyr::filter(line %in% unique(line_prop_changes$line))

gc()

message("There are ", length(unique(proliferation_w_info$line)), " lines in the analysis after missing line filters") # 165

### adding other info
proliferation_w_info = proliferation_w_info %>%
  dplyr::left_join(line_prop_changes) %>%
  distinct()

# sort genes by number of samples in the burdened category
genes = proliferation_w_info %>% 
  dplyr::filter(rare_burden_dichotomy == "yes") %>%
  dplyr::group_by(rare_mutation_type,gene,treatment) %>% 
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

# pdf(paste0(directory,"PRS_phagocytosis_canonical_correlation_analysis.pdf"),
#     width = 6,height = 6)
# 
# plotCorrMatrix(C)
# 
# dev.off()

### regression
mut_types = unique(proliferation_w_info$rare_mutation_type)
for(mut_type in mut_types){
  for(treat in c("untreated","LPS","IFN")){
  sub_table = proliferation_w_info %>%
    dplyr::filter(rare_mutation_type == mut_type & treatment == treat) %>%
    dplyr::mutate(gene_cat = paste(gene,rare_burden,sep = "_"))
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
  for(i in 1:length(unique(c("BCOR",genes)))){
    res_dichotomy = list()
    res = list()
    p_list = list()
    gene =  unique(c("BCOR",genes))[i]
    if(!gene %in% sub_table$gene){
      next()
    }else{

        
        test = sub_table[sub_table$gene == gene, ]
        
        
        
        # rare burden as continuous variable
        
        my_fit=  lmerTest::lmer(scaled_log_fraction ~ sex + prop_unadjusted_min_value +
                               genotypePC1 + genotypePC2 + genotypePC3 + genotypePC4 + genotypePC5 + 
                                         rare_burden + # effect of each unit of change in allele dosage in treatment vs baseline
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
        
        p_list[[paste0(treat,"_",gene,"_",mut_type)]] = jtools::effect_plot(my_fit, 
                                                               pred = rare_burden, 
                                                               interval = FALSE, partial.residuals = TRUE)  + 
          ggtitle(paste0(treat," - ",gene),subtitle = mut_type)
        
        if(i%% 500 ==0)  { print(paste0("i is ",i))}
      
        
        
      
      # save these as rds
      saveRDS(object = res,paste0(directory,mut_type,"/",i,"_",mut_type,"_",treat,"_burden_lm_res.rds"))
      saveRDS(object = p_list,paste0(directory,mut_type,"/",i,"_",treat,"partial_residual_LMM_plots.rds"))
      
      
    }
    rm(res)
    gc()
  }
  }
  
  
  # remove that mut type and clean to save memory
  
  proliferation_w_info = proliferation_w_info %>%
    dplyr::filter(rare_mutation_type != mut_type)
  gc()
}

