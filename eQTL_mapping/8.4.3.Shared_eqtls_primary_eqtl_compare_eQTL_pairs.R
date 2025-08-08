# comparing eGene-eQTL pairs between our eQTL and primary
library(tidyverse)
library(qvalue)
source("./functions.R")
set.seed(123)
outDir = "../../data/results/8.colocalisation_analysis/coloc_results/primary_eQTL_Fujita_2024_ROSMAP/"
name = basename(outDir)
external_eQTL_folder = as.character("../../../resources/primary_eQTL_Fujita_2024_ROSMAP/") 


### load lead variants per eGene, for proportion analysis
external_eQTL_lead_path = list.files(external_eQTL_folder,pattern = "*permuted.tsv.gz",full.names = TRUE) 

if(!is_empty(external_eQTL_lead_path)){
  external_eQTL_variants_path = list.files(external_eQTL_folder,pattern = "*permuted.tsv.gz",full.names = TRUE) 
  
  
  external_eQTL_lead =  readr::read_tsv(external_eQTL_lead_path) %>%
    dplyr::rename(variant_id = variant,chr = chromosome,snp_pos = position) %>%
    # omit genes without significant eQTLs
    dplyr::filter(!is.na(variant_id)) %>%
    dplyr::filter(!is.na(p_beta)) %>%
    # calculating q values
    dplyr::mutate(qval =  qvalue(p_beta)$qvalues) %>%
    dplyr::filter(qval <0.05) %>%
    dplyr::rename(gene_id = molecular_trait_id) %>%
    dplyr::mutate(gene_variant_pair = paste0(gene_id,"-",str_remove(variant_id,"chr")))
  
} else{
  external_eQTL_path = list.files(external_eQTL_folder,pattern = "*.all.tsv.gz",full.names = TRUE) 
  # exclude index, if present
  external_eQTL_path = external_eQTL_path[!grepl("tbi",external_eQTL_path)]
  
  external_eQTL_lead = readr::read_tsv(external_eQTL_path) %>%
    # check I'm taking the correct pval - should be p_beta but some (Fujita) don't have that
    dplyr::rename(standard_error = if_else(condition = "se" %in% colnames(.),true = "se", false = "standard_error")) %>%
    dplyr::rename(ref = REF) %>%
    dplyr::mutate(variant = paste(chromosome,position,ref,alt,sep = "_")) %>%
    dplyr::select(molecular_trait_id,variant, beta,pvalue,alt, chromosome,position,
                  standard_error,maf,significant_by_2step_FDR) %>%
    dplyr::rename(variant_id = variant,
                  pval_nominal = pvalue,
                  effect_allele = alt,
                  gene_id = molecular_trait_id) %>%
    dplyr::mutate(chr = chromosome,
                  pos=position, # ensuring types are correct
                  varbeta = standard_error^2) %>%
    # no need to adjust for N, see https://github.com/chr1swallace/coloc/issues/14
    
    dplyr::select(gene_id,chr, pos,variant_id,effect_allele, beta,varbeta,pval_nominal,significant_by_2step_FDR) %>%
    dplyr::filter(significant_by_2step_FDR == "Yes")
  
  external_eQTL_lead = external_eQTL_lead %>%
    dplyr::group_by(gene_id) %>%
    dplyr::arrange(pval_nominal) %>%
    dplyr::slice_head(n=1)  %>% # take the lead variant per eGene as in tensorQTL
    dplyr::ungroup() %>%
    dplyr::mutate(gene_variant_pair = paste0(gene_id,"-",str_remove(variant_id,"chr")))
  
  # check I'm taking the correct pval - should be pval_beta but some (Fujita) don't have that
  
  external_eQTL_lead = external_eQTL_lead %>%
    dplyr::rename(pval_beta = pval_nominal)
  
}


# load our significant results
eqtl = readr::read_csv("../../data/results/4.Inspect_eQTL_results/tensorQTL_variant_gene_60PCs.csv")
# extract significant results
signif = eqtl %>%
  dplyr::filter(group %in% c(paste0("60_",c("untreated","IFN","LPS"),"_Not_proliferating"))) %>%
  dplyr::filter(qval<0.05) %>%
  dplyr::mutate(id_variant_pair = paste0(gene_id,"-",variant_id))

# from external eQTL dataset, load all the variants and subset to those shared with ours
external_eQTL_path = list.files(external_eQTL_folder,pattern = "*.all.tsv.gz",full.names = TRUE) 
# exclude index, if present
external_eQTL_path = external_eQTL_path[!grepl("tbi",external_eQTL_path)]

if(is_empty(external_eQTL_lead_path)){
  
external_eQTL_shared_nominal = readr::read_tsv(external_eQTL_path) %>%
  dplyr::rename(standard_error = if_else(condition = "se" %in% colnames(.),true = "se", false = "standard_error")) %>%
  dplyr::rename(ref = REF) %>%
  dplyr::mutate(variant = paste(chromosome,position,ref,alt,sep = "_")) %>%
  dplyr::select(molecular_trait_id,variant, beta,pvalue,alt, chromosome,position,
                standard_error,maf) %>%
  dplyr::rename(variant_id = variant,
                pval_nominal = pvalue,
                effect_allele = alt,
                gene_id = molecular_trait_id) %>%
  dplyr::mutate(chr = chromosome,
                pos=position, # ensuring types are correct
                varbeta = standard_error^2) %>%
  # no need to adjust for N, see https://github.com/chr1swallace/coloc/issues/14
  
  dplyr::select(gene_id,chr, pos,variant_id,effect_allele, beta,varbeta,pval_nominal) %>%
  dplyr::mutate(id_variant_pair = paste0(gene_id,"-",str_remove(variant_id,"chr"))) %>%
  dplyr::filter(id_variant_pair %in% signif$id_variant_pair) %>%
  dplyr::rename(p_beta = pval_nominal)


} else {
  # take p_beta directly
}
### getting the same eGene-eQTL pairs
# then testing the proportion of significant eGenes - eQTL pairs (non-null p-vals)

prop_null_p_beta = qvalue(external_eQTL_shared_nominal$p_beta)$pi0
shared_eQTL_prop = 1-prop_null_p_beta

shared_eQTL_list = list(shared_eQTL = external_eQTL_shared_nominal,
                        shared_eQTL_prop = shared_eQTL_prop)
write_rds(shared_eQTL_list,paste0(outDir,"shared_eQTL_prop_list_pi0.rds"))

