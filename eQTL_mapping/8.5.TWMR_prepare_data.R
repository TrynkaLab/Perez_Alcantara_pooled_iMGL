# adapted from Daianna's edits based on my initial script
library(tidyverse)
library(dplyr)
# library(biomaRt)
library(sessioninfo)
source("./functions.R")

################################################################################
##            Transcriptome Wide Mendelian Randomization (TWMR)
################################################################################

#-------------------------------------------------------------------------------
#     1 Selection of instrumental variables and exposures per focal gene*
#-------------------------------------------------------------------------------
#  Code to define the final set of genes (exposures) and their SNP instrumental 
#  variables (IVs) for multivariable MR to test the association between the 
#  microglia expression of each gene with proliferation from preMac -> microglia  
#  in IFN/LPS/untreated.
#  * Note: code based on Marta's code and analysis.
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

## Job index


args = commandArgs(trailingOnly = TRUE)

# Checking arguments
if (length(args) < 4) {
  stop("4 arguments must be supplied", call. = FALSE)
}

focus_gene = as.character(args[1]) # gene to analyse as exposure for TWMR
treatment = as.character(args[2]) # treatment, for coloc, eQTL and GWAS files
phenotype = as.character(args[3]) # phenotype for GWAS
output_dir = as.character(args[4]) # Output directory. Will create a treatment subdirectory

#### for testing
# focus_gene = "SCAPER"
# treatment = "untreated"
# phenotype = "phagocytosis"
# output_dir = "../../data/results/8.5.eQTL_MR/TWMR/input/"

###
dir.create(output_dir, recursive = T, showWarnings = F)


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#             1.1 Multivariable model definition per focal gene 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


## Function to find lead coloc variant per focal gene 
find_lead_coloc_variant <- function(focus_gene, treatment){
  
  coloc_res = coloc_res[[treatment]]
  
  ## Find top coloc variant for focal gene, if available
  if (sum(str_split_i(names(coloc_res), "_", i = 2) %in% focus_gene) == 0) {
    lead_coloc_var = NA # gene not tested
  }
  else{
    coloc_res = coloc_res[[names(coloc_res)[str_split_i(names(coloc_res), "_", i = 2) %in% focus_gene]]]
    
    if (is_null(coloc_res)) {
      lead_coloc_var = NA # no coloc results as gene didn't fill the min num of SNPs needed for coloc
    }
    else{
      
      if (coloc_res$summary["PP.H4.abf"] < 0.70) {
        lead_coloc_var = NA # no significant colocs 
        
      } else{
        # select top SNP with largest PP for H4
        lead_coloc_var = coloc_res$results %>%
          dplyr::select("snp", "SNP.PP.H4") %>%
          dplyr::arrange(desc(SNP.PP.H4)) %>%
                                        # snp                               SNP.PP.H4
        #   1                           12_40177839_TA_T 0.3184834028059976840374645234987838194
        # 2                            12_40178345_T_C 0.1501720073409992806112711605237564072
        # 3                            12_40189297_T_G 0.1501720073409992806112711605237564072
        # 4                            12_40179672_C_A 0.1182284369180451738534642913691641297
        # 5                            12_40183894_T_A 0.1182284369180451738534642913691641297
          slice_max(SNP.PP.H4, with_ties = FALSE, n = 1) %>% 
          dplyr::pull(snp)
      }
    }
  }
  
  return(lead_coloc_var)
}


## Function to find pairs of variants in high LD (R^2 >0.2)
## (this function expects a list of variants in the same chr)
variants_in_LD <- function(variants){
  
  ## LD information
  ## All cis-eQTLs for eGenes included in a model are in the same chr
  ld = readr::read_delim(
    paste0(
      "/lustre/scratch125/humgen/resources/TopLD/release_28032022/EUR/SNV/EUR_chr",
      unique(stringr::str_split_i(variants, pattern = "_", i = 1)),
      "_no_filter_0.2_1000000_LD.csv.gz")) %>% 
    dplyr::filter(
      SNP1 %in% stringr::str_split_i(variants, pattern = "_", i = 2) &
        SNP2 %in% stringr::str_split_i(variants, pattern = "_", i = 2)) %>%
    dplyr::distinct()
  
  ## Subset to our variants 
  ld_filt = ld %>%
    dplyr::mutate(
      
      SNP1 = paste(
        stringr::str_split_i(variants[1], pattern = "_", i = 1),
        str_split_i(Uniq_ID_1, pattern = ":", i = 1),
        str_split_i(Uniq_ID_1, pattern = ":", i = 2),
        str_split_i(Uniq_ID_1, pattern = ":", i = 3),
        sep = "_"),
      
      SNP2 = paste(
        stringr::str_split_i(variants[1], pattern = "_", i = 1),
        str_split_i(Uniq_ID_2, pattern = ":", i = 1),
        str_split_i(Uniq_ID_2, pattern = ":", i = 2),
        str_split_i(Uniq_ID_2, pattern = ":", i = 3),
        sep = "_")
    ) %>%
    dplyr::mutate(SNP1 = pmin(SNP1, SNP2), 
                  SNP2 = pmax(SNP1, SNP2)) %>%
    dplyr::select(SNP1, SNP2, R2) %>%
    dplyr::filter(SNP1 %in% variants & SNP2 %in% variants) %>% 
    dplyr::filter(R2 > 0.2) %>% 
    dplyr::distinct()  
  
  return(ld_filt)
  
}


## Function to compute corr between eQTL effects for all gene pairs in a gene list 
gene_pair_cor <- function(initial_genes_in_model){
  
  cors <- matrix(NA, nrow = length(initial_genes_in_model), ncol = length(initial_genes_in_model))
  colnames(cors) <- rownames(cors) <- initial_genes_in_model
  
  for(i in 1:length(initial_genes_in_model)){
    for(j in 1:length(initial_genes_in_model)){
      
      gene1 = initial_genes_in_model[i]
      gene2 = initial_genes_in_model[j]
      
      ## all cis-eQTLs for gene1
      gene1_eQTLs <- nominal_eqtl_results %>% 
        dplyr::select(gene, variant_id, pval_nominal, slope) %>% 
        dplyr::filter(gene == gene1) %>% 
        unique
      ## all cis-eQTLs for gene2
      gene2_eQTLs <- nominal_eqtl_results %>% 
        dplyr::select(gene, variant_id, pval_nominal, slope) %>% 
        dplyr::filter(gene == gene2) %>% 
        unique
      
      ## Subset to shared variants
      shared_variants <- inner_join(gene1_eQTLs, gene2_eQTLs, by = "variant_id", suffix = c("1", "2"))
      
      ## Eliminate duplicates (likely to be variants in high LD with same effect sizes)
      shared_indep_variants <- shared_variants[-which(duplicated(shared_variants[, c("slope1", "slope2")])), ]
      
      if(gene1 == gene2){
        cor = 1
      }
      else{
        # If less than 3 shared assumed indep. variants, set cor = 0 (very few common eQTLs between genes)
        if(length(shared_indep_variants$variant_id)<3){
          cor = 0
        }
        else{
          cor = cor(shared_indep_variants$slope1, shared_indep_variants$slope2, method = "pearson")
        }
        
      }
      
      cors[gene1, gene2] <- cor
      
      # Previous code to obtain an independent set of shared eQTLs:
      # if(!gene1 == gene2){
      #   ## Prun variants in high LD
      #   shared_variants_in_LD <- variants_in_LD(shared_variants$variant_id)
      #   
      #   remaining_shared_variants <- shared_variants$variant_id
      #   independent_shared_variants <- shared_variants$variant_id
      #   
      #   while(length(remaining_shared_variants) >0){
      #     
      #     ## Pick a random variant
      #     random_var <- remaining_shared_variants[1]
      #     
      #     ## Exclude cor variants 
      #     random_and_cor_variants <- shared_variants_in_LD %>% filter(SNP1 %in% random_var | SNP2 %in% random_var) %>% .[, c("SNP1", "SNP2")] %>% unlist %>% 
      #       as.vector() %>% append(random_var) %>%  unique()
      #     cor_variants <- random_and_cor_variants[-which(random_and_cor_variants == random_var)]
      #     
      #     if(any(cor_variants %in% independent_shared_variants)){
      #       independent_shared_variants <- independent_shared_variants[-which(independent_shared_variants %in% cor_variants)]
      #     }
      #     
      #     ## Repeat on remaining variants
      #     remaining_shared_variants <-  remaining_shared_variants[-which(remaining_shared_variants %in% random_and_cor_variants)]
      #   } 
      #   
      #   if(length(independent_shared_variants)>=3){
      #     
      #     independent_shared_variants_effects = shared_variants %>% filter(variant_id %in% independent_shared_variants)
      #     
      #     ## Correlation between shared indep. cis-eQTL effects (not subsetting to significant ones!!)
      #     cor = cor(independent_shared_variants_effects$slope1, independent_shared_variants_effects$slope2, method = "pearson")
      #     
      #   }
      #   ## 
      #   else{
      #     cor = 0
      #   }
      #   
      # }
      # else{
      #   cor = 1
      # }
      
    }
  }
  
  return(cors)
}


## Function to define exposures and IVs x focal gene in multivariable model
select_IVs_and_exposures <- function(focus_gene, treatment){
  
  # ****************************************************************************
  #    Step 1:  Select causal + strongest independent cis-eQTLs of focal gene
  # ****************************************************************************
  
  # Select strongest cis-eQTL(s) x gene (with min pval)
  gene_cis_eQTLs <- nominal_eqtl_results %>% 
    dplyr::select(variant_id, gene, pval_nominal, slope) %>%
    dplyr::filter(gene == focus_gene) %>% 
    dplyr::slice_min(pval_nominal, na_rm = TRUE)  
  
  ## Remove variants with duplicated pval and effect size (SNPs in LD)
  if(length(which(duplicated(gene_cis_eQTLs$slope)))>0){
    gene_cis_eQTLs <- gene_cis_eQTLs[-which(duplicated(gene_cis_eQTLs$slope)),]
  }
  
  ## Select lead coloc cis-eQTL (if any)
  lead_coloc_var <- find_lead_coloc_variant(focus_gene, treatment)
  if(!is.na(lead_coloc_var)){
    gene_lead_coloc_cis_eQTL <- nominal_eqtl_results %>% 
      dplyr::select(variant_id, gene, pval_nominal, slope) %>%
      dplyr::filter(gene == focus_gene) %>% 
      dplyr::filter(variant_id == lead_coloc_var) 
    
    ## If lead coloc is in LD with any cis-eQTL already included, in step 1.1 the latter will be removed
    gene_cis_eQTLs <- rbind(gene_cis_eQTLs, gene_lead_coloc_cis_eQTL) %>% 
      dplyr::distinct() 
  }
  
  # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  #     Step 1.1: if >1 cis-eQTL keep the independent strongest ones 
  # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  if(length(gene_cis_eQTLs$variant_id) > 1){
    
    ####### IMPORTANT ##########
    #
    # Variant pairs in high LD if R^2 >0.2
    #
    ##########
    variant_pairs_in_LD = variants_in_LD(gene_cis_eQTLs$variant_id) 
    
    # 1) first discard all variants in LD with the lead coloc variant 
    # 2) then discard all variants in LD with the variant with the strongest effect size
    # 3) then repeat with the second strongest variant among the remaining and so on ...
    
    if(is.na(lead_coloc_var)){
      ordered_eQTLs <- gene_cis_eQTLs %>% arrange(-abs(slope)) %>% dplyr::select(variant_id) %>% unlist %>% as.vector()
    }
    else{
      ordered_eQTLs <- c(lead_coloc_var,
                         gene_cis_eQTLs %>% 
                           dplyr::filter(!variant_id == lead_coloc_var) %>% 
                           dplyr::arrange(-abs(slope)) %>% 
                           dplyr::select(variant_id) %>% 
                           unlist %>% as.vector())
    }
    
    independent_variants <- gene_cis_eQTLs$variant_id
    while(length(ordered_eQTLs)>0){
      
      ## Causal/strongest variant
      strong_eQTL = ordered_eQTLs[1]
      
      ## Correlated variants + strongest
      strong_and_cor_variants = variant_pairs_in_LD %>% 
        dplyr::filter(SNP1 %in% strong_eQTL | SNP2 %in% strong_eQTL) %>% 
        .[, c("SNP1", "SNP2")] %>% 
        unlist() %>% 
        append(strong_eQTL) %>% unique
      
      ## Keep strongest variant and discard the ones correlated with it
      cor_variants = strong_and_cor_variants[-which(strong_and_cor_variants == strong_eQTL)]
      if(any(cor_variants %in% independent_variants)){
        independent_variants <- independent_variants[-which(independent_variants %in% cor_variants)]
      }
      ## Repeat with remaining variants 
      ordered_eQTLs <- ordered_eQTLs[-which(ordered_eQTLs %in% strong_and_cor_variants)]
    }
    
    ## Confirm lead coloc variant included
    if((!is.na(lead_coloc_var)) & (!lead_coloc_var %in% independent_variants)){
      message("Lead coloc variant not included in the set of independent cis-eQTLs for focal gene.")
      stop()
    }
  }
  
  else{
    independent_variants = gene_cis_eQTLs$variant_id
  }
  gene_indep_cis_eQTLs <- gene_cis_eQTLs %>% 
    dplyr::filter(variant_id %in% independent_variants)
  
  
  # ****************************************************************************
  #      Step 2: Add all eGenes of included independent stronger cis-eQTL(s)
  # ****************************************************************************
  indep_cis_eQTLs_eGenes <- nominal_eqtl_results %>% 
    dplyr::select(variant_id, gene, pval_nominal, slope) %>%
  dplyr::filter(variant_id %in% independent_variants) 
  
  gene_indep_cis_eQTLs_with_eGenes = rbind(gene_indep_cis_eQTLs, indep_cis_eQTLs_eGenes) %>% 
    unique
  
  # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  #          Step 2.1: exclude highly-correlated genes (r^2 >=0.4)
  # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  initial_genes_in_model <- gene_indep_cis_eQTLs_with_eGenes$gene %>% unique()
  
  if(length(initial_genes_in_model)>1){
    ## Compute cor between independent cis-eQTL effects for each pair of genes
    cors <- gene_pair_cor(initial_genes_in_model)
    # heatmap(cors, Rowv = NA, Colv = NA, scale = "none")
    
    ## First discard eGenes correlated with focal gene
    uncorr_eGenes <- initial_genes_in_model
    focus_gene_cor <- cors[focus_gene, ] %>% .[abs(.)>=0.4] %>% names %>% .[!.== focus_gene]
    if(any(focus_gene_cor %in% uncorr_eGenes)){
      uncorr_eGenes <- uncorr_eGenes[-which(uncorr_eGenes %in% focus_gene_cor)]
    }
    remaining_eGenes <- initial_genes_in_model[-which(initial_genes_in_model %in% c(focus_gene_cor, focus_gene))]
    
    ## Second, take less corr gene with focal gene and discard its correlated genes
    remaining_eGenes <- cors[focus_gene, ] %>% .[remaining_eGenes] %>% abs() %>% sort() %>% names
    
    while(length(remaining_eGenes)>0){
      
      ## Take less corr gene
      gene = remaining_eGenes[1]
      
      corr_genes <- cors[gene, ] %>% .[abs(.)>=0.4] %>% names %>% .[!. == gene] 
      if(any(corr_genes %in% uncorr_eGenes)){
        uncorr_eGenes <- uncorr_eGenes[-which(uncorr_eGenes %in% corr_genes)]
      }
      remaining_eGenes <- remaining_eGenes[-which(remaining_eGenes %in% c(gene, corr_genes))]
    }
    # heatmap(cors[uncorr_eGenes, uncorr_eGenes], Rowv = NA, Colv = NA, scale = "none")
  }
  
  else{
    uncorr_eGenes <- focus_gene
  }
  
  gene_indep_cis_eQTLs_with_indep_eGenes <- gene_indep_cis_eQTLs_with_eGenes %>% 
    dplyr::filter(gene %in% uncorr_eGenes)
  
  
  # ****************************************************************************
  #          Step 3: Add all signif (p < 10-6) cis-eQTLs x added eGene
  # ****************************************************************************
  indep_eGenes_signif_cis_eQTLs <- nominal_eqtl_results %>% 
    dplyr::select(variant_id, gene, pval_nominal, slope) %>%
  dplyr::filter(gene %in% gene_indep_cis_eQTLs_with_indep_eGenes$gene & pval_nominal < 0.000001) %>%
  dplyr::filter(! gene == focus_gene) # don't add cis-eQTLs for focal gene that were not added in step 1 !!!
  
  duplicates <- which(duplicated(indep_eGenes_signif_cis_eQTLs[, c("gene", "pval_nominal", "slope")]))
  
  indep_eGenes_signif_cis_eQTLs <- indep_eGenes_signif_cis_eQTLs[-duplicates, ]
  
  # ****************************************************************************
  #                  Step 4: Prun variants in LD (R^2 >= 0.1) 
  # ****************************************************************************
  if(dim(indep_eGenes_signif_cis_eQTLs)[1] > 0){
    
    all_variants <- unique(c(indep_eGenes_signif_cis_eQTLs$variant_id, gene_indep_cis_eQTLs$variant_id))
    all_variants_in_LD <- variants_in_LD(all_variants)
    
    ## Keep focal gene indep. cis-eQTLs (from Step 1) and discard its corr variants
    ordered_final_variants <- unique(c(gene_indep_cis_eQTLs$variant_id, indep_eGenes_signif_cis_eQTLs$variant_id)) 
    independent_final_variants <- ordered_final_variants
    
    while(length(ordered_final_variants)>0){
      
      eQTL = ordered_final_variants[1]
      
      ## Variants in LD with cis-eQTL
      eQTL_and_cor_variants = all_variants_in_LD %>% 
      dplyr::filter(SNP1 %in% eQTL | SNP2 %in% eQTL) %>% .[, c("SNP1", "SNP2")] %>% unlist() %>% append(eQTL) %>% unique
      
      cor_variants = eQTL_and_cor_variants[-which(eQTL_and_cor_variants == eQTL)]
      if(any(cor_variants %in% independent_final_variants)){
        independent_final_variants <- independent_final_variants[-which(independent_final_variants %in% cor_variants)]
      }
      ## Repeat with remaining variants 
      ordered_final_variants <- ordered_final_variants[-which(ordered_final_variants %in% eQTL_and_cor_variants)]
    }
    
    ## Confirm all focal gene causal and strong cis-eQTLs are included 
    if(length(which(!gene_indep_cis_eQTLs$variant_id %in% independent_final_variants)) >0){
      message(paste0("Not all independent cis-eQTLs of focal gene ", focus_gene, " are included in the final model."))
      stop()
    }
    if((!is.na(lead_coloc_var)) & !(lead_coloc_var %in% independent_final_variants)){
      message(paste0("Lead coloc cis-eQTL of focal gene ", focus_gene, " not included in the final model."))
      stop()
    }
    
  }
  else{
    independent_final_variants <- gene_indep_cis_eQTLs$variant_id
  }
  
  
  ## Subset to final set of independent eQTLs
  gene_indep_cis_eQTLs_with_indep_eGenes_signif_cis_eQTLs <- rbind(gene_indep_cis_eQTLs_with_indep_eGenes, indep_eGenes_signif_cis_eQTLs) %>% unique
  gene_indep_cis_eQTLs_with_indep_eGenes_signif_cis_eQTLs = gene_indep_cis_eQTLs_with_indep_eGenes_signif_cis_eQTLs %>% 
  dplyr::filter(variant_id %in% independent_final_variants)
  
  
  ## Final set of IVs and exposures
  final_model = gene_indep_cis_eQTLs_with_indep_eGenes_signif_cis_eQTLs %>% 
    dplyr::select(variant_id, gene, slope) %>%
    tidyr::pivot_wider(
      id_cols = variant_id,
      names_from = gene,
      values_from = slope,
      values_fill = 0) # 0 effect for eQTLs ns in other genes
  
  ## Confirm final set of uncorr genes (with focal gene) and independent cis-eQTLs
  if(length(which(!gene_indep_cis_eQTLs$variant_id %in% final_model$variant_id)) >0){
    message(paste0("Not all independent cis-eQTLs of focal gene ",
                   focus_gene,
                   " are included in the final model in ", treatment, "."))
    stop()
  }
  
  if(!focus_gene %in% colnames(final_model)){
    message(paste0("Focal gene ", focus_gene, " not included in the final model in ", treatment, "."))
    stop()
  }
  
  if(!setequal(colnames(final_model)[-1], uncorr_eGenes)){
    message(paste0("Not all uncorrelated genes for ", focus_gene, 
                   " included in the final model in ", treatment, "."))
    stop()
  }
  
  print(paste0("Final model for ", focus_gene, " in ", treatment, ": ", dim(final_model)[2]-2, " uncorrelated eGenes as exposures and ", 
               dim(final_model)[1], " independent cis-eQTLs as IVs"))
  
  final_model = final_model %>% 
    column_to_rownames("variant_id")
  return(final_model)
  
}

load_nominal_results = function(treatment) {
  nominal_res = readr::read_delim(
    paste0(
      "../../data/results/tensorqtl/best_results/sum_sizefactorsNorm_log2_scaled_centered_",
      treatment,
      "_Not_proliferating_common_500kb_window_tensorQTL_nominal.txt"
    )
  )
  
  
  # mapping ensembl to symbol
  ensembl_to_symbol = ensembl_to_gene_name(
    IDs = unique(nominal_res$phenotype_id),
    IDFrom = "GENEID",
    IDTo = "SYMBOL"
  )$map %>%
    dplyr::rename(phenotype_id = From, gene = To) %>%
    as_tibble() %>%
    dplyr::distinct(phenotype_id, .keep_all = TRUE)
  
  nominal_res = nominal_res %>%
    dplyr::left_join(ensembl_to_symbol)
  # adding info and filtering
  return(nominal_res)
  
}


## Test function with example focal genes to avoid headaches later :)
## (examples in untreated)

# * genes with 1 independent cis-eQTL that has only 1 eGene (the focal gene):
# E <- select_IVs_and_exposures(focus_gene = "VANGL1", treatment = "untreated")
# [1] "Final model for VANGL1 in untreated: 0 uncorrelated eGenes as exposures and 1 independent cis-eQTLs as IVs" 
# E <- select_IVs_and_exposures(focus_gene = "ROCK1", treatment = "untreated")
# [1] "Final model for ROCK1 in untreated: 0 uncorrelated eGenes as exposures and 1 independent cis-eQTLs as IVs"

# * genes with lead coloc variant as unique independent cis-eQTL: 
# E <- select_IVs_and_exposures(focus_gene = "OTOA", treatment = "untreated")
# [1] "Final model for OTOA in untreated: 1 uncorrelated eGenes as exposures and 1 independent cis-eQTLs as IVs"
# E <- select_IVs_and_exposures(focus_gene = "LY86", treatment = "untreated")
# [1] "Final model for LY86 in untreated: 0 uncorrelated eGenes as exposures and 1 independent cis-eQTLs as IVs"

# * gene with 1 independent cis-eQTL with multiple eGenes with signif eQTLs
# E <- select_IVs_and_exposures(focus_gene = "A2M", treatment = "untreated")
# [1] "Final model for A2M in untreated: 2 uncorrelated eGenes as exposures and 24 independent cis-eQTLs as IVs"

# * gene with 1 independent cis-eQTL with multiple eGenes with 0 signif eQTLs
# E <- select_IVs_and_exposures(focus_gene = "SAMD11", treatment = "untreated")
# [1] "Final model for SAMD11 in untreated: 4 uncorrelated eGenes as exposures and 1 independent cis-eQTLs as IVs"

# * gene with 2 independent cis-eQTLs with 1 eGene only (the focal gene)
# E <- select_IVs_and_exposures(focus_gene = "APBB2", treatment = "untreated")
# [1] "Final model for APBB2 in untreated: 0 uncorrelated eGenes as exposures and 2 independent cis-eQTLs as IVs"

# * gene with 2 (initial) independent cis-eQTLs with multiple eGenes with signif eQTLs
# E <- select_IVs_and_exposures(focus_gene = "BLOC1S3", treatment = "untreated")
# [1] "Final model for BLOC1S3 in untreated: 8 uncorrelated eGenes as exposures and 21 independent cis-eQTLs as IVs"

# * gene with 2 independent cis-eQTLs with multiple eGenes with 0 signif eQTLs
# E <- select_IVs_and_exposures(focus_gene = "BCKDK", treatment = "untreated")
# [1] "Final model for BCKDK in untreated: 2 uncorrelated eGenes as exposures and 2 independent cis-eQTLs as IVs"



## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#                    1.2 Prepare input matrices for TWMR
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# for(treatment in c("IFN", "LPS", "untreated")){
  
  ## Define output dir
  dir.create(paste(output_dir,treatment,phenotype,sep="/"), recursive = T, showWarnings = F)
  
  ## Read cis-eQTL summary stats x treatment 
  nominal_eqtl_results = load_nominal_results(treatment)
  
  ## Read GWAS summary stats x treatment 
  pheno_gwas = readr::read_csv(
    file = paste0(
      "../../../OTAR2065_phenotypic_QTLs/data/results/",
      phenotype,
      "/2.check_association_results/lm_1pct_filtered_deflated/all_res_pqtl_",
      phenotype,
      ".csv"
    )
  )

  ## Coloc results
  coloc_res = readr::read_rds(
    file = paste0(
      "../../../OTAR2065_phenotypic_QTLs/data/results/",
      phenotype,
      "/3.coloc/coloc_results/my_GWAS_my_eQTL/coloc_results_single_causal_variant_500kb.rds"
    )
  )

  
  ## Find model x focal gene 
    ## E: IVs x eGenes
    E = select_IVs_and_exposures(focus_gene, treatment) 
    
    ## * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    #  Warning 1:                   !!!                            *
    #  Effect of variants with ns/non-tested effects on any of the *
    #  eGenes included in E are assumed to be 0.                   *                              
    ## * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    
    ## Make sure focal gene is the 1st gene column in E
    if(which(colnames(E) == focus_gene) != 1){
      eGenes <- c(colnames(E), colnames(E)[1])
      eGenes[1] <- focus_gene
      eGenes <- eGenes[-which(duplicated(eGenes))]
      E <- E[, eGenes]
    }
    
    ## Find pair-wise LD matrix for IVs (C: IV1 x IV2)
    C = matrix(0, ncol = nrow(E), nrow = nrow(E))
    rownames(C) <- colnames(C) <- rownames(E)
    
    ## * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    #  Warning 2:                   !!!                            *
    #  Variant pairs not present in LD file are assumed to be 0.   *
    ## * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    
    if(nrow(E) >1){
      ld = readr::read_delim(
        paste0(
          "/lustre/scratch125/humgen/resources/TopLD/release_28032022/EUR/SNV/EUR_chr",
          unique(stringr::str_split_i(rownames(E), pattern = "_", i = 1)),
          "_no_filter_0.2_1000000_LD.csv.gz")) %>% 
      dplyr::filter(
          SNP1 %in% stringr::str_split_i(rownames(E), pattern = "_", i = 2) &
            SNP2 %in% stringr::str_split_i(rownames(E), pattern = "_", i = 2)) %>% 
        dplyr::mutate(
          SNP1 = paste(
            stringr::str_split_i(rownames(E)[1], pattern = "_", i = 1),
            str_split_i(Uniq_ID_1, pattern = ":", i = 1),
            str_split_i(Uniq_ID_1, pattern = ":", i = 2),
            str_split_i(Uniq_ID_1, pattern = ":", i = 3),
            sep = "_"),
          SNP2 = paste(
            stringr::str_split_i(rownames(E)[1], pattern = "_", i = 1),
            str_split_i(Uniq_ID_2, pattern = ":", i = 1),
            str_split_i(Uniq_ID_2, pattern = ":", i = 2),
            str_split_i(Uniq_ID_2, pattern = ":", i = 3),
            sep = "_")
        ) %>%
        dplyr::mutate(SNP1 = pmin(SNP1, SNP2), 
                      SNP2 = pmax(SNP1, SNP2)) %>%
        ## Subset to our IVs in E 
        dplyr::filter(SNP1 %in% rownames(E) & SNP2 %in% rownames(E)) %>% 
        dplyr::select(SNP1, SNP2, R2) %>%
        dplyr::distinct() 
      
      ## Fill matrix
      if(dim(ld)[1] > 0){
        for(i in 1:dim(ld)[1]){
          pair = ld[i, ]
          SNP1 = pair[, "SNP1"] %>% as.character()
          SNP2 = pair[, "SNP2"] %>% as.character()
          C[SNP1, SNP2] <-  C[SNP2, SNP1] <- pair[, "R2"] %>% as.numeric()
        }
      }
    }
    
    ## SNP LD with itself = 1
    diag(C) <- 1
    
    ## Confirm C contains all IVs, is symmetric, and has 1s on the diagonal
    if(any(rownames(C) != rownames(E)) | any(colnames(C) != rownames(E))){
      message(paste0("Not all IVs for ", focus_gene, " are present in C in ", treatment, "."))
      stop()
    }
    if(!isSymmetric(C)){
      message(paste0("C for ", focus_gene, " is not symmetric in ", treatment, "."))
      stop()
    }
    if(any(! diag(C) == 1)){
      message(paste0("C for ", focus_gene, " has not all 1s on the diagonal in ", treatment, "."))
      stop()
    } 
    
    ## IVs effects on proliferation from preMac -> microglia in treatment (G: IVs)
    pheno_gwas_subset = pheno_gwas %>%
      dplyr::filter(snp %in% rownames(E)) %>%
      dplyr::select(contains(treatment), snp) %>%
      dplyr::select(contains("coef"), snp) %>%
      dplyr::rename(variant_id = snp,
                    BETA_GWAS = paste0("coef_", treatment)) %>%
      dplyr::mutate(BETA_GWAS = tidyr::replace_na(BETA_GWAS, 0)) %>% #  NA variant effects are 0
      column_to_rownames("variant_id")
      # reorder like E
      G = pheno_gwas_subset[rownames(E),]
      names(G) = rownames(E)
    
    ## * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    #  Warning 3:                   !!!                            *
    #  Effect of variants with no GWAS effect available are set to *
    #  0.                                                          *
    ## * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    
    ## Sample sizes for eQTL, GWAS, and LD  
      Ngwas <- case_when(
        treatment == "untreated" & phenotype == "phagocytosis"  ~ 133,
        treatment == "IFN" & phenotype == "phagocytosis" ~ 139,
        treatment == "LPS"& phenotype == "phagocytosis"  ~ 140,
        treatment == "untreated" & phenotype == "migration"  ~ 147,
        treatment == "IFN" & phenotype == "migration" ~ 92,
        treatment == "LPS"& phenotype == "migration"  ~ 147,
        )
      
      N_eQTLs = case_when(
        treatment == "untreated" ~ 189,
        treatment == "IFN" ~ 188,
        treatment == "LPS" ~ 188)
      
      n_for_LD = 13160 # info from TopMED readme
   
    
    inputs <- list("E" = E, "G" = G, "C" = C, 
                   "N_eQTLs" = N_eQTLs, "Ngwas" = Ngwas, "n_for_LD" = n_for_LD)
    saveRDS(inputs, file = paste0(output_dir, "/",treatment,"/",phenotype,
                                  "/", focus_gene, "_input_matrices.rds"))
    
  


## Try example genes (in untreated):
# * focal gene with 1 IV only: CLPTM1
# * focal gene with IVs not in LD file: A2M







# Reproducibility info
# session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
    # R version 4.4.1 (2024-06-14)
    # Platform: x86_64-pc-linux-gnu
    # Running under: Ubuntu 22.04.4 LTS
    # 
    # Matrix products: default
    # BLAS:   /software/hgi/installs/softpack/rstudio/.spack/linux-ubuntu22.04-x86_64_v3/gcc-11.4.0/r-4.4.1-zf3d5qbxgbiqsk4ddke3fl6uluwcbqcu/rlib/R/lib/libRblas.so 
    # LAPACK: /software/hgi/installs/softpack/rstudio/.spack/linux-ubuntu22.04-x86_64_v3/gcc-11.4.0/r-4.4.1-zf3d5qbxgbiqsk4ddke3fl6uluwcbqcu/rlib/R/lib/libRlapack.so;  LAPACK version 3.12.0
    # 
    # locale:
    #   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8       
    # [4] LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
    # [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
    # [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
    # 
    # time zone: Europe/Belfast
    # tzcode source: internal
    # 
    # attached base packages:
    #   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
    # 
    # other attached packages:
    #   [1] sessioninfo_1.2.2       TwoSampleMR_0.5.6       ensembldb_2.26.0        AnnotationFilter_1.26.0
    # [5] GenomicFeatures_1.54.4  AnnotationDbi_1.64.1    Biobase_2.62.0          GenomicRanges_1.54.1   
    # [9] GenomeInfoDb_1.38.8     IRanges_2.36.0          S4Vectors_0.40.2        AnnotationHub_3.10.0   
    # [13] BiocFileCache_2.10.2    dbplyr_2.5.0            BiocGenerics_0.48.1     ldscr_0.1.0            
    # [17] mrSampleOverlap_0.1.1   lubridate_1.9.3         forcats_1.0.0           stringr_1.5.1          
    # [21] dplyr_1.1.4             purrr_1.0.2             readr_2.1.5             tidyr_1.3.1            
    # [25] tibble_3.2.1            ggplot2_3.5.1           tidyverse_2.0.0        
    # 
    # loaded via a namespace (and not attached):
    #   [1] splines_4.4.1                 later_1.3.2                   BiocIO_1.12.0                
    # [4] bitops_1.0-7                  filelock_1.0.3                preprocessCore_1.64.0        
    # [7] XML_3.99-0.16.1               lifecycle_1.0.4               mixsqp_0.3-54                
    # [10] rstatix_0.7.2                 lattice_0.22-6                vroom_1.6.5                  
    # [13] MASS_7.3-65                   backports_1.4.1               magrittr_2.0.3               
    # [16] rmarkdown_2.26                plotly_4.10.4                 locuszoomr_0.2.1             
    # [19] limma_3.58.1                  vcfR_1.15.0                   yaml_2.3.8                   
    # [22] httpuv_1.6.15                 cowplot_1.1.3                 DBI_1.2.2                    
    # [25] abind_1.4-5                   zlibbioc_1.48.2               RCurl_1.98-1.14              
    # [28] rappdirs_0.3.3                GenomeInfoDbData_1.2.11       ggrepel_0.9.5                
    # [31] irlba_2.3.5.1                 gdata_3.0.0                   nortest_1.0-4                
    # [34] vegan_2.6-4                   mr.raps_0.2                   permute_0.9-7                
    # [37] codetools_0.2-20              DelayedArray_0.28.0           xml2_1.3.6                   
    # [40] tidyselect_1.2.1              shape_1.4.6.1                 farver_2.1.1                 
    # [43] viridis_0.6.5                 matrixStats_1.2.0             jsonlite_1.8.8               
    # [46] GenomicAlignments_1.38.2      pinfsc50_1.3.0                survival_3.5-8               
    # [49] iterators_1.0.14              foreach_1.5.2                 tools_4.4.1                  
    # [52] progress_1.2.3                Rcpp_1.0.12                   glue_1.7.0                   
    # [55] gridExtra_2.3                 SparseArray_1.2.4             xfun_0.43                    
    # [58] mgcv_1.9-1                    MatrixGenerics_1.14.0         withr_3.0.0                  
    # [61] BiocManager_1.30.22           fastmap_1.1.1                 LDlinkR_1.4.0                
    # [64] fansi_1.0.6                   digest_0.6.35                 timechange_0.3.0             
    # [67] R6_2.5.1                      mime_0.12                     colorspace_2.1-0             
    # [70] gtools_3.9.5                  biomaRt_2.58.2                RSQLite_2.3.6                
    # [73] UpSetR_1.4.0                  utf8_1.2.4                    generics_0.1.3               
    # [76] MRlap_0.0.3.3                 data.table_1.15.4             rtracklayer_1.62.0           
    # [79] htmlwidgets_1.6.4             prettyunits_1.2.0             httr_1.4.7                   
    # [82] S4Arrays_1.2.1                pkgconfig_2.0.3               gtable_0.3.4                 
    # [85] blob_1.2.4                    XVector_0.42.0                htmltools_0.5.8              
    # [88] susieR_0.12.35                carData_3.0-5                 gggrid_0.2-0                 
    # [91] ProtGenerics_1.34.0           scales_1.3.0                  png_0.1-8                    
    # [94] knitr_1.45                    rstudioapi_0.16.0             tzdb_0.4.0                   
    # [97] rjson_0.2.21                  checkmate_2.3.1               nlme_3.1-164                 
    # [100] curl_5.2.1                    zoo_1.8-12                    cachem_1.0.8                 
    # [103] BiocVersion_3.18.1            parallel_4.4.1                vsn_3.70.0                   
    # [106] restfulr_0.0.15               reshape_0.8.9                 pillar_1.9.0                 
    # [109] grid_4.4.1                    vctrs_0.6.5                   promises_1.2.1               
    # [112] ggpubr_0.6.0                  car_3.1-2                     xtable_1.8-4                 
    # [115] cluster_2.1.6                 dtplyr_1.3.1                  evaluate_0.23                
    # [118] cli_3.6.2                     compiler_4.4.1                Rsamtools_2.18.0             
    # [121] rlang_1.1.3                   crayon_1.5.2                  ggsignif_0.6.4               
    # [124] labeling_0.4.3                affy_1.80.0                   plyr_1.8.9                   
    # [127] fs_1.6.3                      stringi_1.8.3                 SCpubr_2.0.2                 
    # [130] coloc_5.2.3                   viridisLite_0.4.2             BiocParallel_1.36.0          
    # [133] munsell_0.5.1                 Biostrings_2.70.3             lazyeval_0.2.2               
    # [136] glmnet_4.1-8                  Matrix_1.6-5                  hms_1.1.3                    
    # [139] patchwork_1.2.0               bit64_4.0.5                   KEGGREST_1.42.0              
    # [142] statmod_1.5.0                 shiny_1.8.1.1                 SummarizedExperiment_1.32.0  
    # [145] interactiveDisplayBase_1.40.0 broom_1.0.5                   memoise_2.0.1                
    # [148] affyio_1.72.0                 bit_4.0.5                     ape_5.7-1    