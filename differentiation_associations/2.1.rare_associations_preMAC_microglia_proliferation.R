# rare variant burden associations
.libPaths(c('/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/',"/software/teamtrynka/ma23/R4.1/libs",
            "/software/R-4.1.0/lib/R/library"))
library(patchwork)
library(tidyverse)
library(microbenchmark)
library(lmerTest)
source("./functions.R")

input_dir = "../../data/results/1.alluvial_plots"
output_dir = "../../data/results/2.1.rare_associations_preMAC_microglia_proliferation"
dir.create(output_dir, recursive = TRUE)
data = read.table(paste0(input_dir,
                         "/pools2-11_13-15_changing_props_iPSC_preMacs_microglia_WGS_sc.txt"),
                  header = TRUE)
# read in rare variant burden
#  3 types of variant consequence classes,  annotated by VEP program from Sanger (release 99) and Ensembl (v 75):
# 
#   - Deleterious (union of loss-of-funcion and missense pathogenic, missense variants with CADD score > 15)
# - Loss-of-function ( named as protein-truncating variants): 
# - Synonymous
rare_burden_allgenes = read.csv("../../../resources/Puigdevall_Neuroseq_efficiency_2023/TableS1.csv")
rare_burden_pergene = list()
rare_burden_pergene[["del"]] = readRDS("../../../resources/Puigdevall_Neuroseq_efficiency_2023/mutBurdenTabs/mutBurden_del.RDS")
rare_burden_pergene[["ptv"]] = readRDS("../../../resources/Puigdevall_Neuroseq_efficiency_2023/mutBurdenTabs/mutBurden_ptv.RDS")
rare_burden_pergene[["syn"]] = readRDS("../../../resources/Puigdevall_Neuroseq_efficiency_2023/mutBurdenTabs/mutBurden_synonymous.RDS") 

# change rownames to gene column, and rename other columns to simple line names
rare_burden_pergene <- purrr::map(rare_burden_pergene, ~fix_burden_matrices(.))


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

line_info=read.csv("../../../OTAR2065_phenotypic_QTLs/data/allinfo_hipsci_PD_AD_PRS_nondisease_feederfree_european.csv") %>%
  dplyr::mutate(sex = case_when(grepl(pattern = "Female",.$Sex) ~ "Female",
                                grepl(pattern = "Male",.$Sex) ~ "Male")) %>%
  dplyr::select(Line,sex) %>%
  dplyr::rename(line=Line)

efficiency = efficiency %>%
  dplyr::left_join(.,line_info)


message("Adding genotype PCs")

full_genotype_pcs = read.table("../../../OTAR2065_sc_eQTL/data/genotype/plink_genotypes/all_pools.genotype.MAF05.eigenvec")
rownames(full_genotype_pcs) = ifelse(full_genotype_pcs$V1 ==full_genotype_pcs$V2,
                                     yes = full_genotype_pcs$V1,
                                     no = paste(full_genotype_pcs$V1,full_genotype_pcs$V2,sep = "_"))
full_genotype_pcs = full_genotype_pcs[c(2:7)] %>%
  dplyr::rename(line=V2,genotypePC1 = V3,genotypePC2 = V4,genotypePC3 = V5,genotypePC4 = V6,genotypePC5 = V7)

efficiency_w_info = efficiency %>%
  dplyr::inner_join(.,full_genotype_pcs) %>%
  dplyr::distinct() 

rm(full_genotype_pcs)
gc()

# cleaning up Inf and NAs
efficiency_w_info = efficiency_w_info %>%
  dplyr::mutate(log_scaled_proportion = log(scaled_proportion)) %>%
  dplyr::filter(!log_scaled_proportion %in% c(-Inf,Inf,NA,NaN)) %>%
  dplyr::filter(!is.na(sex))# remove missing data causing problems - 4 donors


# clean up rare variant burden table to contain only donors remaining after QC filters
rare_burden_pergene = purrr::map(rare_burden_pergene, ~subset_burden_shared_lines(x=.,y = unique(efficiency_w_info$line)))
length(colnames(rare_burden_pergene$del))-1 # 170 donors in both rare burden matrix and my post QC proliferation data
nrow(rare_burden_pergene$del) # 19,653 genes

# retain for each type of rare variant and gene that is present (burden>=1) in over 5% of donors (8 donors) AND absent in at least 5% donors (otherwise dichotomised won't work)
rare_burden_pergene = purrr::map(rare_burden_pergene, ~filter_rare_variant_genes(.,fraction=0.05))
nrow(rare_burden_pergene$del) # 7,757 genes
nrow(rare_burden_pergene$ptv) # 1,012 genes
nrow(rare_burden_pergene$syn) # 9,080 genes

# check that the gene is expressed in microglia to at least 1CPM in pseudobulk
# expressed_genes = list()
# for(condition in c("untreated","IFN","LPS")){
# expressed_genes[[condition]] = read_table(paste0("../../../OTAR2065_differentiation_efficiency/data/results/1.2.Inspect_integration/",
#                                                  condition,
#                                                  "_pseudobulk_aggregated_raw_counts_CPMs.txt"))
# expressed_genes[[condition]] = expressed_genes[[condition]] %>%
#   dplyr::filter(CPM > 1)
# }
# expressed_genes = unique(c(expressed_genes$untreated$gene,
#                          expressed_genes$IFN$gene,
#                          expressed_genes$LPS$gene))
# 
# rare_burden_pergene$del = rare_burden_pergene$del %>%
#   dplyr::filter(gene %in% expressed_genes)
# rare_burden_pergene$ptv = rare_burden_pergene$ptv %>%
#   dplyr::filter(gene %in% expressed_genes)
# rare_burden_pergene$syn = rare_burden_pergene$syn %>%
#   dplyr::filter(gene %in% expressed_genes)
# 
# nrow(rare_burden_pergene$del) # 4,324 genes
# nrow(rare_burden_pergene$ptv) # 520 genes
# nrow(rare_burden_pergene$syn) # 5,370 genes
# convert rare variant burden to long format
rare_burden_pergene = rare_variant_to_long(rare_burden_pergene)

#add to efficiency table
efficiency_w_info = efficiency_w_info %>%
  dplyr::inner_join(.,rare_burden_pergene) %>%
  dplyr::distinct()

rm(rare_burden_pergene)
gc()
######## regression

res_dichotomy = list()
res = list()
result <- microbenchmark(
  for(mut_type in unique(efficiency_w_info$rare_mutation_type)){
    sub_table = efficiency_w_info %>%
      dplyr::filter(rare_mutation_type == mut_type)
    res_dichotomy[[mut_type]] = list()
    res[[mut_type]] = list()
    
  for(i in 1:length(unique(sub_table$gene))){
    gene =  unique(sub_table$gene)[i]
    
    test = sub_table[sub_table$gene == gene, ] %>%
    # dicotomising first burden: 0 == no, else == yes
    dplyr::mutate(rare_burden_dichotomy = factor(case_when(rare_burden == 0 ~ "no",
                                                    .default = "yes")),
                  treatment = factor(treatment)) %>%
      dplyr::mutate(rare_burden_dichotomy = relevel(rare_burden_dichotomy,ref = c("no"), # setting no as baseline
                                                    ),
                    treatment = relevel(treatment,ref = c("untreated"))) # setting untreated as baseline

    
    res_dichotomy[[mut_type]][[i]]=  lmerTest::lmer(log_scaled_proportion ~ sex + sequencing + 
                              #genotypePC1 + genotypePC2 + genotypePC3 + genotypePC4 + genotypePC5 + 
                               treatment*rare_burden_dichotomy + # effect of each unit of change in allele dosage in treatment vs baseline
                               #(1|sex) + (1 | sequencing) +  
                               (1|line) + (1|pool),
                             data = test, 
                             control = lmerControl(optCtrl=list(maxfun=5000) ))

    
    # sum_res = summary(my_fit)
    # res_dichotomy[[mut_type]][[i]] = list(
    #   coefficients = sum_res$coefficients,
    #   varcor =  sum_res$varcor,
    #   residuals =  sum_res$residuals,
    #   call =  sum_res$call,
    #   AICtab =  sum_res$AICtab,
    #   gene = gene)
    
  # rare burden as continuous variable
    
    res[[mut_type]][[i]]=  lmerTest::lmer(log_scaled_proportion ~  sex + sequencing + 
                              #genotypePC1 + genotypePC2 + genotypePC3 + genotypePC4 + genotypePC5 + 
                              treatment*rare_burden + # effect of each unit of change in allele dosage in treatment vs baseline
                              #(1|sex) + (1 | sequencing) +  
                              (1|line) + (1|pool),
                            data = test, 
                            control = lmerControl(optCtrl=list(maxfun=5000) ))
    
    
    # sum_res = summary(my_fit)
    # res[[mut_type]][[i]] = list(
    #   coefficients = sum_res$coefficients,
    #   varcor =  sum_res$varcor,
    #   residuals =  sum_res$residuals,
    #   call =  sum_res$call,
    #   AICtab =  sum_res$AICtab,
    #   gene = gene)
    
    if(i%% 1000 ==0)  { print(paste0("i is ",i))}
  }
  
  
  }
  ,times=1)

# save these as rds
saveRDS(object = res,paste0(output_dir,"/res.rds"))
saveRDS(object = res_dichotomy,paste0(output_dir,"/res_dichotomy.rds"))




