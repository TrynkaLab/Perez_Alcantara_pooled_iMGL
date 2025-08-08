# checking PTK2B region

library(coloc)
library(patchwork)
library(tidyverse)
source("./functions.R")

input_path = "../../data/results/5.proportion_QTL/"
output_path="../../data/results/8.check_results/"
dir.create(file.path(output_path),recursive = TRUE)

pqtl_res_files = list()
pqtl_res_files[["migration"]] = read.csv("/lustre/scratch123/hgi/projects/otar2065/OTAR2065_phenotypic_QTLs/data/results/1.2.proportion_QTL_mixed_models/all_res_253.csv")
pqtl_res_files[["phagocytosis"]] = read.csv("../../data/results/phagocytosis/1.2.proportion_QTL_mixed_models/prop_QTL_253.csv")


# preparing harmonised GWAS for coloc

GWAS_folder = as.character("/lustre/scratch123/hgi/teams/trynka/resources/summary_statistics/public/GCST90027158") 

# GWAS data
GWAS_path = list.files(paste0(GWAS_folder,"/harmonised"),pattern = "*.h.tsv.gz",full.names = TRUE) #there should be just one file with that pattern

gwas = readr::read_tsv(GWAS_path) %>%
  dplyr::rename(standard_error = if_else(condition = "se" %in% colnames(.),true = "se", false = "standard_error")) %>%
  dplyr::select(starts_with("hm_var"), hm_beta,p_value,hm_effect_allele, hm_chrom,hm_pos,
                standard_error,standard_error,effect_allele_frequency) %>%
  dplyr::rename(variant_id = starts_with("hm_var"),
                beta = hm_beta,
                pval_nominal = p_value,
                effect_allele = hm_effect_allele,
                maf=effect_allele_frequency) %>%
  dplyr::mutate(chr = as.numeric(hm_chrom),
                pos=as.numeric(hm_pos), # ensuring types are correct
                varbeta = standard_error^2) %>%
  # no need to adjust for N, see https://github.com/chr1swallace/coloc/issues/14
  
  dplyr::select(chr, pos,variant_id,effect_allele, beta,varbeta,pval_nominal)


message("Printing colnames of GWAS file:")
print( colnames(gwas))

gwas_list = list()
for(assay in c("migration","phagocytosis")){
  gwas_list[[assay]] = gwas %>%
    dplyr::filter(variant_id %in% pqtl_res_files[[assay]]$snp)
}

rm(gwas)

message("Formatting phenotype GWAS data")

pqtl_res_files_format = list()
for(assay in c("migration","phagocytosis")){
  for(treatment in c("untreated","IFN","LPS")){
    # split by treatment
    pqtl_res_files_format[[paste(assay,treatment,sep = "_")]] =  pqtl_res_files[[assay]] %>%
      dplyr::select(contains(treatment),snp) %>%
      dplyr::rename(beta = contains("coef"), # the coef is the beta
                    pval_nominal = contains("p_"),
                    se = contains("se_"),
                    variant_id=snp)  %>%
      tidyr::separate_wider_delim(cols =  variant_id, names = c("chr", "pos","major_allele","minor_allele"), 
                                  delim = "_",cols_remove = FALSE) %>%
      dplyr::mutate(varbeta = se^2,
                    chr = as.numeric(chr),
                    pos=as.numeric(pos)) %>%
      dplyr::select(chr, pos,minor_allele,major_allele, variant_id , beta,varbeta,pval_nominal) %>%
      dplyr::arrange(chr,pos)
  }
  
}


write_rds(pqtl_res_files_format,"../../data/results/PTK2B/phenotype_associations_formatted.rds")
write_rds(gwas_list,"../../data/results/PTK2B/bellenguez_AS_GWAS_formatted.rds")

### colocalisation #####

pqtl_res_files_format = read_rds("../../data/results/PTK2B/phenotype_associations_formatted.rds")
gwas_list = read_rds("../../data/results/PTK2B/bellenguez_AS_GWAS_formatted.rds")

message("Performing colocalisation with coloc")

coloc_PTK2B= list()
for(pheno_treat in names(pqtl_res_files_format)){
  assay = str_split_1(pheno_treat,pattern = "_")[1]
  
  # remove all NA values from beta in phenotype file and gwas file
  pqtl_res_files_format[[pheno_treat]] = pqtl_res_files_format[[pheno_treat]] %>%
  dplyr::filter(!is.na(beta))
  gwas_list[[assay]] = gwas_list[[assay]] %>%
    dplyr::filter(!is.na(beta))
  
  
  subset_coloc = list(beta = pqtl_res_files_format[[pheno_treat]]$beta,
                           varbeta = pqtl_res_files_format[[pheno_treat]]$varbeta,
                           snp = pqtl_res_files_format[[pheno_treat]]$variant_id,
                           position=pqtl_res_files_format[[pheno_treat]]$pos,
                           type="quant",
                           sdY= 1 
  )
  if(is.null(check_dataset(subset_coloc)) & length(pqtl_res_files_format[[pheno_treat]]$variant_id)>50){
    message("phenotype GWAS list for coloc is valid")
  } else{
    message("WARNING: phenotype GWAS list for coloc is NOT valid")
    
  }
  
  
  # Convert each column to a named list element
  subset_gwas_variants = gwas_list[[assay]] %>%
    dplyr::filter(variant_id %in% pqtl_res_files_format[[pheno_treat]]$variant_id)
  
  if(nrow(subset_gwas_variants) == nrow(pqtl_res_files_format[[pheno_treat]])){
    message("phenotype GWAS and AD GWAS subsets have the same number of variants")
  } else{
    message("WARNING: phenotype GWAS and AD GWAS subsets DO NOT have the same number of variants")
    message("fixing...")
    pqtl_res_files_format[[pheno_treat]] = pqtl_res_files_format[[pheno_treat]] %>%
      dplyr::filter(variant_id %in% gwas_list[[assay]]$variant_id)
    
    if(nrow(subset_gwas_variants) == nrow(pqtl_res_files_format[[pheno_treat]])){
      message("phenotype GWAS and AD GWAS subsets now have the same number of variants")
    } else{
      message("WARNING: phenotype GWAS and AD GWAS subsets STILL DO NOT have the same number of variants")
    }
  }
  
  subset_gwas_coloc = list(beta = subset_gwas_variants$beta,
                           varbeta = subset_gwas_variants$varbeta,
                           snp = subset_gwas_variants$variant_id,
                           position=subset_gwas_variants$pos,
                           type="cc"
  )
  
  if(is.null(check_dataset(subset_gwas_coloc)) & length(gwas_list[[assay]]$variant_id)>50){
    message("GWAS list for coloc is valid")
  } else{
    message("WARNING: GWAS list for coloc is NOT valid")
    
  }
  
  message("Running coloc")
  
  if(is.null(check_dataset(subset_gwas_coloc)) & length(gwas_list[[assay]]$variant_id)>50 & is.null(check_dataset(subset_coloc)) & length(pqtl_res_files_format[[pheno_treat]]$variant_id)>50 ){
    coloc_PTK2B[[pheno_treat]]= coloc.abf(dataset1=subset_gwas_coloc, # AD GWAS
                                                   dataset2=subset_coloc, # out phenotype association
                                                   p12=1e-6)
  } else{
    coloc_PTK2B[[pheno_treat]]= NULL
  }
  # gives posterior probability (Bayesian: not p-vals!!) for each of the following hypotheses:
  # H0: neither trait has a genetic association in the region
  # H1: only trait 1 has a genetic association in the region
  # H2: only trait 2 has a genetic association in the region
  # H3: both traits are associated, but with different causal variants
  # H4: both traits are associated and share a single causal variant
}


for(pheno_treat in names(coloc_PTK2B)){
  pdf(paste0("../../data/results/PTK2B/coloc_",pheno_treat,".pdf"),
      width = 6, height = 5)
  par(mfcol = c(2, 2))
  sensitivity(coloc_PTK2B[[pheno_treat]],rule="H4 > 0.5",
              plot.manhattans = TRUE,preserve.par = TRUE)
  dev.off()
  }


# nothing looks particularly good 

#  plot association lead pQTL variant -boxplot


# and the one associated by AD GWAS
