# examine proliferation SKATO tests
library(patchwork)
library(tidyverse)
library(qvalue)
source("./functions.R")
directory = "../../../data/results/phagocytosis/5.1.2.check_SKAT_WES_vs_phagocytosis/"
dir.create(directory, recursive = TRUE)

options(future.globals.maxSize = 20000 * 1024^2) # 20Gb


for(mut_type in c("deleterious")){
  res_per_treat = list()
  for(treat in c("untreated","IFN","LPS")){
  res = list()
  files = list.files("../../../data/results/phagocytosis/5.2.rare_var_SKAT_WES/",
                            pattern = paste0("*","_",treat,"_",mut_type,"_SKATO_result.RDS"))
  files = str_sort(files, numeric = TRUE) # to later recover right plots by number
  numbers = str_split_i(files,pattern = "_",i = 1) #they may not be contiguous
  res = lapply(paste0("../../../data/results/phagocytosis/5.2.rare_var_SKAT_WES/",files), function(x){
    read_rds(x) 
  })
  
  res_pvals <- purrr::map(res, ~extract_SKATO_results(.))

  res_pvals = do.call("rbind",res_pvals)

  res_per_treat[[treat]] = res_pvals %>%
    dplyr::relocate(gene_name,gene_id)
    
  res_per_treat[[treat]]$p_Bonf = p.adjust(res_per_treat[[treat]]$p_val,method = "bonferroni")
  res_per_treat[[treat]]$treatment = treat
  res_per_treat[[treat]]$number = numbers

  res_per_treat[[treat]] = res_per_treat[[treat]] %>%
    dplyr::arrange(p_Bonf)

  }
  res = do.call("rbind",res_per_treat)
  hist(res$p_val) # why are p-values so suspiciously significant all over?
  hist(res$p_Bonf) 
  write_csv(res,paste0(directory,mut_type,"_SKATO_scaled_centered_prop_pvals.csv"))
}


table(res$p_Bonf<0.05,res$treatment)
# IFN  LPS untreated
# FALSE 5256 7703      7656
# TRUE   143   20        67

### plot all significant genes
# check variants individually
WES_missense_del_gt_path ="../../../../hipsci_genotype_processing/data/WES/VEP/missense_clean_GRCh38.vcf.gz"
WES_missense_del_cleaninfo_path ="../../../../hipsci_genotype_processing/data/WES/VEP/missense_clean_GRCh38.txt"
line_prop_changes_path = "../../../data/results/phagocytosis/1.check_line_proportions/line_prop_changes_per_well.csv"
genotype_pc_path = "../../../../OTAR2065_sc_eQTL/data/genotype/plink_genotypes/all_pools.no_outliers.genotype.MAF05.eigenvec"
 line_info_path="../../../../OTAR2065_differentiation_efficiency/data/donor_metadata_complete_with_imputed_sex.csv"
# general metadata
line_info = read_csv(line_info_path) %>%
  dplyr::select(donor,Sex)
# PCs
geno_PCs = read_tsv(genotype_pc_path) %>%
  dplyr::select(-`#FID`) %>%
  dplyr::rename(
    line = IID,
    genotypePC1 = PC1,
    genotypePC2 = PC2,
    genotypePC3 = PC3,
    genotypePC4 = PC4,
    genotypePC5 = PC5
  ) %>%
  dplyr::mutate(donor = if_else(str_detect(line,"-"),
                                str_split_i(line,pattern="-",1),
                                str_split_i(line,pattern="_",1))) %>%
  dplyr:::select(donor,line:genotypePC5)

# add to line info

line_info = line_info %>%
  dplyr::left_join(geno_PCs)

# line proportion changes (phagocytosis)
line_prop_changes = read_csv(line_prop_changes_path) %>%
  dplyr::mutate(donor = if_else(str_detect(line,"-"),
                                str_split_i(line,pattern="-",1),
                                str_split_i(line,pattern="_",1))) %>%
  dplyr::filter(!log_fraction_mean %in% c(NA,Inf,-Inf)) %>%
  dplyr::select(log_fraction_mean,replicate,line,donor,condition,treatment,pool,prop_unadjusted_min_value) %>%
  distinct()

message("There are ", length(unique(line_prop_changes$line)), " lines in the analysis after Na and Inf filters") # 202

## agregate all phenotype and covariate information

full_pheno_info = line_prop_changes %>%
  dplyr::left_join(line_info) %>%
  # omit donors that are outliers (no genotype PC info) and sh5y5y
  dplyr::filter(!is.na(genotypePC1))
colnames(full_pheno_info)

# clean VEP info ###
WES_missense_del_info = read_tsv(WES_missense_del_cleaninfo_path) %>%
  dplyr::mutate(CADD_PHRED = case_when(CADD_PHRED == "." ~ NA,
                                       .default = as.double(CADD_PHRED))) %>%
  # defining as deleterious - union of LoF and missense pathogenic mutations
  # LoF -  frameshift, stop-gain, transcript_ablation, splice acceptor or splice donor variants
  # see https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
  # missense pathogenic - annotated as missense or start loss with a CADD Phred score cutoff > 15 (as Pau did)
  dplyr::mutate(category = case_when((str_detect(Consequence,pattern="missense") | 
                                        str_detect(Consequence,pattern="start_lost") | 
                                        str_detect(Consequence,pattern="protein_altering")) & CADD_PHRED >=15 ~ "deleterious",
                                     (str_detect(Consequence,pattern="frameshift")|
                                        str_detect(Consequence,pattern="stop_gained")|
                                        str_detect(Consequence,pattern="transcript_ablation")|
                                        str_detect(Consequence,pattern="splice_acceptor")|
                                        str_detect(Consequence,pattern="splice_donor")
                                     ) ~ "deleterious",
                                     .default = "missense_non_deleterious"),
                ID = paste(CHROM,POS,REF,ALT,sep = "_")) %>%
  dplyr::filter(!duplicated(ID))

length(unique(WES_missense_del_info$Gene)) # 18,588 genes, 18,581 gene IDs (SYMBOL)
table(WES_missense_del_info$category) # 69,223 deleterious, 62,022 missense non-deleterious variants

# clean VCF ###
WES_missense_del_genotype = vcfR::read.vcfR(WES_missense_del_gt_path)

# INFO section
info = WES_missense_del_genotype@fix[,1:5] %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(ID = paste(CHROM,POS,REF,ALT,sep = "_"))
# GT section
WES_missense_del_genotype = vcfR::extract_gt_tidy(WES_missense_del_genotype)

WES_missense_del_genotype = WES_missense_del_genotype %>%
  dplyr::mutate(genotype = case_when(gt_GT == "1/1" ~ 2,
                                     gt_GT %in% c("0/1","1/0") ~ 1,
                                     gt_GT == "0/0" ~ 0,
                                     .default = NA),
                line = if_else(str_detect(Indiv,"HPSI"), 
                               str_split_i(Indiv,"-",i=2),
                               NA)) %>%
  dplyr::mutate(genotype = as.numeric(genotype),
                donor = str_split_i(line,"_",i=1)) 



table(info$ID %in% WES_missense_del_info$ID) # lost a handful after filters
table(WES_missense_del_info$ID %in% info$ID) 

# join genotype info
WES_missense_del_genotype = cbind(WES_missense_del_genotype,info)

### subset genes


for(treat in c("untreated","IFN","LPS")){
  message("Working on ", treat, " treatment")
  message("\n")
  
  
  gene_df= res %>%
    dplyr::filter(treatment == treat & p_Bonf < 0.05) %>%
    distinct()
  
  vars_within_mut_type = WES_missense_del_info %>%
    dplyr::filter(category == mut_type & Gene %in% gene_df$gene_id ) %>%
    dplyr::select(ID,Gene,SYMBOL) 
  message("subset right metadata") 
  
  meta_subset = full_pheno_info %>%
    # filter objects to shared lines
    dplyr::filter((treatment == treat) & (line %in% WES_missense_del_genotype$line)) %>%
    dplyr::mutate(sex_numeric = case_when(Sex == "Male" ~ 1, # female as reference (0)
                                          .default = 0)) %>%
    dplyr::group_by(line) %>%
    # summarising per line, pool, across all reps
    dplyr::reframe(log_fraction_mean_pool = mean(log_fraction_mean),
                   prop_unadjusted_min_value_mean = mean(prop_unadjusted_min_value),
                   sex_numeric = sex_numeric, genotypePC1= genotypePC1, genotypePC2 =genotypePC2) %>%
    dplyr::select(line,log_fraction_mean_pool,prop_unadjusted_min_value_mean,sex_numeric, genotypePC1, genotypePC2) %>%
    distinct()

  
  pdf(paste0(directory,mut_type,"_",treat,"_SKATO_vars_vs_phago.pdf"),
      width = 4, height = 5)
for(gene in gene_df$gene_id){
  vars = vars_within_mut_type %>%
    dplyr::filter(Gene == gene) %>%
    dplyr::left_join(WES_missense_del_genotype) %>%
    dplyr::select(ID,Gene, SYMBOL, line,genotype,gt_GT_alleles) %>%
    dplyr::left_join(meta_subset)
  plist = list()
  for(var in unique(vars$ID)){
    
  plist[[var]] =  vars %>%
    dplyr::filter(ID == var) %>%
    dplyr::arrange(genotype) %>%
    dplyr::mutate(gt_GT_alleles = factor(gt_GT_alleles, levels = unique(gt_GT_alleles), ordered = TRUE)) %>%
    ggplot(aes(x=gt_GT_alleles,y=log_fraction_mean_pool)) +
    geom_boxplot() +
    geom_point(alpha = 0.1) +
    theme_bw()
  }
  p = patchwork::wrap_plots(plist)
 plot(p)
}
  
  dev.off()
}


################################
# any in Sam's list?
sam_phago = read_tsv("../../../../CRISPR/OTAR2065_phagocytosis_CRISPR/data/2024_07_sam_screen_results.tsv") %>%
  dplyr::arrange(FDR) %>%
  dplyr::rename(gene = id)

directionality = sam_phago %>%
  dplyr::left_join(res[res$treatment=="untreated",]) %>%
  dplyr::mutate(direction = case_when((estimate < 0 & Score > 0 & p_Bonf < 0.05) | (estimate > 0 & Score < 0 & p_Bonf < 0.05) ~ "same direction",
                                      (estimate < 0 & Score < 0 & p_Bonf < 0.05) | (estimate > 0 & Score > 0 & p_Bonf < 0.05)  ~ "opposite direction",
                                      .default = "NS burden"),
                direction_no_significance = case_when((estimate < 0 & Score > 0 ) | (estimate > 0 & Score < 0 ) ~ "same direction",
                                                      (estimate < 0 & Score < 0 ) | (estimate > 0 & Score > 0 )  ~ "opposite direction",
                                                      .default = "Not present"))


table(directionality$direction)
table(directionality$direction_no_significance)
# coin toss for chance of same direction
write_csv(directionality,paste0("../../../data/results/phagocytosis/5.rare_vars_vs_phagocytosis/del/", number,"_untreated_sam_CRISPR_directionality.csv") )
