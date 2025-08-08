# for coloc - running on my own GWAS and my eQTLs
library(tidyverse)
library(tidyr)
library(stringr)
library(purrr)

options(scipen = 999) # prevents from showing very small numbers as 0

args = commandArgs(trailingOnly=TRUE)

# Checking arguments
if (length(args)<4) {
  stop("4 arguments must be supplied", call. = FALSE)
}

my_GWAS_file = as.character(args[1]) # my phenotype GWAS file
eQTL_path = as.character(args[2]) # file path of the eQTL results of interest (created in 4.Inspect_eQTL_results.R)
eQTL_nominal_path=as.character(args[3]) # the nominal eQTL results for same conditions (see Snakefile)
distance_threshold=as.numeric(args[4]) # Define the distance threshold around lead variant in bp (e.g., 500,000 bp)


# 
# # test data
# my_GWAS_file = as.character("../../../data/results/migration/2.check_association_results/lm_1pct_filtered_deflated/untreated/untreated_GWAS_GRCh38.tsv.gz")
 # eQTL_path = as.character("../../../../OTAR2065_sc_eQTL/data/results/4.Inspect_eQTL_results/tensorQTL_variant_gene.csv")
# eQTL_nominal_path=as.character("../../../../OTAR2065_sc_eQTL/data/results/tensorqtl/best_results/")
# distance_threshold = 500000

outDir = "../../../data/results/migration/3.coloc/myeQTL_myGWAS_subsets_for_coloc/"
treat = stringr::str_split(my_GWAS_file,"/")[[1]][9]

# GWAS data
gwas = readr::read_tsv(my_GWAS_file) %>%
  dplyr::select(variant_id ,coef ,p_value,    se,   chr, snp_pos,   ALT) %>%
  dplyr::rename(
                beta = coef, # otherwise beta = log(OR) https://github.com/chr1swallace/coloc/issues/149
                pval_nominal = p_value,
                standard_error = se,
                effect_allele = ALT) %>%
  dplyr::mutate(chr = as.numeric(chr),
                pos=as.numeric(snp_pos), # ensuring types are correct
                varbeta = standard_error^2) %>%
  # no need to adjust for N, see https://github.com/chr1swallace/coloc/issues/14
  
  dplyr::select(chr, pos,variant_id,effect_allele, beta,varbeta,pval_nominal)

message("Printing colnames of GWAS file:")
print( colnames(gwas))



## eQTL data
message("Reading in eQTL results from:")
message(eQTL_path)
eqtl_full = readr::read_csv(eQTL_path)
message("Head of the eQTL result file:")
print(head(eqtl_full))

# Subset to Non proliferating and the right treatment, and significant q value
eqtl_full = eqtl_full %>%
  dplyr::filter(cluster == "Not_proliferating" & treatment == treat & qval<0.05)

# run

  eqtl = eqtl_full %>%
    dplyr::rename(beta = slope) %>% # the slope is the beta
    tidyr::separate_wider_delim(cols =  position, names = c("chr", "pos"), delim = "_") %>%
    dplyr::mutate(varbeta = slope_se^2,
                  chr = as.numeric(chr),
                  pos=as.numeric(pos)) %>%
    dplyr::select(chr, pos,gene_id, gene_name, variant_id , beta,varbeta,pval_nominal, qval) %>%
    dplyr::filter(complete.cases(.)) %>% # will remove rows with NA, which are those from the X chr
    #group_by(variant_id) %>% # select lead associations per variant
    #dplyr::filter(pval_nominal == min(pval_nominal)) %>%
    #ungroup() %>%
    dplyr::arrange(chr,pos) 
  
  message("Printing colnames of eqtl file:")
  print( colnames(eqtl))
  
  # nominal results
  # read in from each condition (treatment x cluster)
  eqtl_nominal=read_tsv(paste0(eQTL_nominal_path,
                               "sum_sizefactorsNorm_log2_scaled_centered_",
                               str_split_fixed(unique(eqtl_full$group),"_",2)[2],
                               "_common_500kb_window_tensorQTL_nominal.txt")) %>%
    dplyr::rename(beta = slope, # the slope is the beta
                  gene_id = phenotype_id
    ) %>%
    tidyr::separate_wider_delim(cols =  variant_id, names = c("chr", "pos","minor_allele","major_allele"), delim = "_",cols_remove = FALSE) %>%
    dplyr::mutate(varbeta = slope_se^2,
                  chr = as.numeric(chr),
                  pos=as.numeric(pos)) %>%
    dplyr::select(chr, pos,minor_allele,major_allele, gene_id, variant_id , beta,varbeta,pval_nominal) %>%
    dplyr::filter(complete.cases(.)) %>% # will remove rows with NA, which are those from the X chr
    #group_by(variant_id) %>% # select lead associations per variant
    #plyr::filter(pval_nominal == min(pval_nominal)) %>%
    #ungroup() %>%
    dplyr::arrange(chr,pos)
  
  ### gathering eqtl and lead GWAS variants that are close enough #####
  
  # from previous paper: https://www.nature.com/articles/s41588-022-01066-3#Sec9
  # For each trait–cell type pair, we applied colocalization to any locus where a 
  # lead variant for a significant eQTL (q value < 0.1) was located within 100kb (100,000 bp) 
  # and in high LD (r2 > 0.5) with a significant GWAS variant (i.e., any GWAS variant
  # with nominal P value < 1 ×10−5, which enabled us to capture suggestive association 
  # signals). In addition, we required at least 50 variants to be available for 
  # testing at each candidate locus. 
  
  # from jeremy's paper:
  # colocalisation tests between GWAS and eQTL signals where the lead variants were within 500 kb of each other, 
  # and passed to coloc all variants within 200 kb of each lead variant.
  
  # subset GWAS summary at same location, and for shared variants with eqtl
  
  # subset in this case by signficant eGenes
  subset_gwas = list()
  subset_eqtl = list()
  n = 1
  for(eGene in paste(eqtl$gene_id, eqtl$gene_name,sep = "_")){
    
    subset_gwas[[eGene]] = gwas %>%
      dplyr::filter(chr == eqtl[n,]$chr & pos <  (eqtl[n,]$pos + distance_threshold) & pos >  (eqtl[n,]$pos - distance_threshold))
    subset_gwas[[eGene]] = subset_gwas[[eGene]]  %>%
      dplyr::filter(!is.na(beta )) 
    
    
    subset_eqtl[[eGene]] =  eqtl_nominal %>%
      dplyr::filter(gene_id %in% stringr::str_split_i(eGene,"_",i=1)) %>%
      dplyr::filter(chr == eqtl[n,]$chr & pos <  (eqtl[n,]$pos + distance_threshold) & pos >  (eqtl[n,]$pos - distance_threshold)) %>%
      dplyr::filter(!is.na(beta ))
    
    if (nrow(subset_eqtl[[eGene]]) == 0 | nrow(subset_gwas[[eGene]]) == 0 ) {
      message(
        "There are no significant eQTL or GWAS variants 500kb around the lead eQTL variant in locus ",
        eGene
      )
      
    } else{
      # check there are no swaps
      
      # if there are, flip the beta
      subset_eqtl[[eGene]] =  subset_eqtl[[eGene]] %>%
        dplyr::mutate(swapped_variant_id = paste(chr,pos,major_allele,minor_allele,sep = "_")) 
      if(any( subset_eqtl[[eGene]]$swapped_variant_id %in% subset_gwas[[eGene]]$variant_id)){
        message(
          "There are swaps in the variant ID in eQTL vs phenotype GWAS for ", eGene, ". Swapping variant_id and beta's sign... "
          
        )
        subset_eqtl[[eGene]] =  subset_eqtl[[eGene]] %>%
          dplyr::mutate(variant_id = case_when(swapped_variant_id %in% subset_gwas[[eGene]]$variant_id == TRUE ~ swapped_variant_id,
                                               .default = variant_id),
                        beta = case_when(swapped_variant_id %in% subset_gwas[[eGene]]$variant_id == TRUE ~ -beta,
                                         .default = beta)) %>%
        
        # if there are duplicated variants now, keep the one with the smallest p-val
        dplyr::group_by(variant_id) %>%
          dplyr::slice_min(pval_nominal) %>%
          dplyr::ungroup()
          
      }
      
      # finally subset to shared variants
      subset_eqtl[[eGene]] = subset_eqtl[[eGene]] %>%
        dplyr::filter(variant_id %in% subset_gwas[[eGene]]$variant_id) %>%
        dplyr::arrange(variant_id)
      subset_gwas[[eGene]] = subset_gwas[[eGene]] %>%
        dplyr::filter(variant_id %in% subset_eqtl[[eGene]]$variant_id) %>%
        dplyr::arrange(variant_id)
      
      if (!identical(subset_eqtl[[eGene]]$variant_id, subset_gwas[[eGene]]$variant_id) ) {
        message(
          "Warning: eQTL and GWAS subset are not identical"
        )
      }else{
        if (nrow(subset_eqtl[[eGene]] ) < 50) {
        message(
          "Warning: There are not enough variants for coloc at locus ", eGene
        )
        }
      }
      
    }
    n = n + 1
  }
  
  
  dir.create(paste0(outDir,treat), recursive = TRUE)
  
  readr::write_rds(x = subset_gwas,file =  paste0(outDir,treat,"_my_GWAS_subsets.rds"))
  readr::write_rds(x = subset_eqtl,file = paste0(outDir,treat,"_eQTL_subsets.rds"))

  

