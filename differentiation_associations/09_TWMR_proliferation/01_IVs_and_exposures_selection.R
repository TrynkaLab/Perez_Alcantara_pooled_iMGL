
library(tidyverse)
library(dplyr)
# library(biomaRt)
library(sessioninfo)

################################################################################
##           9. Transcriptome Wide Mendelian Randomization (TWMR)
################################################################################

#-------------------------------------------------------------------------------
#     9.1 Selection of instrumental variables and exposures per focal gene*
#-------------------------------------------------------------------------------
#  Code to define the final set of genes (exposures) and their SNP instrumental 
#  variables (IVs) for multivariable MR to test the association between the 
#  microglia expression of each gene with proliferation from preMac -> microglia  
#  in IFN/LPS/untreated.
#  * Note: code based on Marta's code and analysis.
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

## Job index
job_i = commandArgs(trailingOnly = TRUE)[1]

## Define dirs
input_dir = paste("..", "..", "input_data", "09_TWMR_proliferation", "01_IVs_and_exposures_selection", sep = "/")
plotdir = paste("..", "..", "plots", "09_TWMR_proliferation", "01_IVs_and_exposures_selection", sep = "/")

dir.create(input_dir, recursive = T, showWarnings = F)
dir.create(plotdir, recursive = T, showWarnings = F)


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#                      9.1.0 Explore input data sets
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

##########################  eQTL summary statistics  ########################### 

## Unique gene-variant pairs
# which(duplicated(nominal_eqtl_results[, c("variant_id", "gene")]))
# integer(0)

# dim(nominal_eqtl_results)                  
# length(unique(nominal_eqtl_results$gene))  

## IFN: 11,664,422 eQTLs across ~10k genes
## LPS: 11,730,468 eQTLs across ~10k genes
## untreated: 12,044,074 eQTLs across ~10k genes


## Add gene symbols
# mart <- useEnsembl(biomart="ENSEMBL_MART_ENSEMBL",
#                    dataset="hsapiens_gene_ensembl",
#                    GRCh = "GRCh38")
# 
# for(treatment in c("IFN", "LPS", "untreated")){
#   nominal_eqtl_results = readr::read_delim(paste0("../../../OTAR2065_sc_eQTL/data/results/tensorqtl/best_results/sum_sizefactorsNorm_log2_scaled_centered_", 
#                                                   treatment, "_Not_proliferating_common_500kb_window_tensorQTL_nominal.txt")) %>% rename(gene = phenotype_id)
#   
#   symbols <- getBM(values = unique(nominal_eqtl_results$gene),
#                    mart = mart,
#                    attributes = c("external_gene_name", "ensembl_gene_id"),
#                    filters = "ensembl_gene_id") %>% rename(gene_symbol = external_gene_name,
#                                                            gene = ensembl_gene_id)
#   
#   nominal_eqtl_results = nominal_eqtl_results %>% left_join(symbols) 
#   write_csv(nominal_eqtl_results, 
#        paste0(input_dir, "/nominal_eqtl_results_", treatment, ".csv"))
# }


## Num eQTLs x gene:
# nominal_eqtl_results %>% .[,"gene"] %>% table %>% as.vector() %>% summary()
# ********** IFN# ********** 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   1     821    1084    1105    1350    9874 
# ********** LPS# ********** 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   1     823    1086    1106    1353    9874 
# ********** untreated# ********** 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   1     822    1085    1105    1352    9874 

## Num of eGenes x eQTL:
# nominal_eqtl_results %>% .[,"variant_id"] %>% table %>% as.vector() %>% summary()
# ********** IFN ********** 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   2.000   3.189   4.000  31.000 
# ********** LPS ********** 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   2.000   3.182   4.000  31.000 
# ********** untreated ********** 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   2.000   3.236   4.000  31.000 


## Num signif eQTLs x gene:
# nominal_eqtl_results %>% filter(pval_nominal < 0.05) %>% .[,"gene"] %>% table %>% as.vector() %>% summary()
# ********** IFN ********** 
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    1      41     120     188     265    5725 
# ********** LPS ********** 
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1.0    42.0   122.0   188.1   267.0  4976.0 
# ********** untreated ********** 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1.0    45.0   128.0   196.1   278.0  5032.0 


## Num of eGenes x signif eQTL:
# nominal_eqtl_results %>% filter(pval_nominal < 0.05) %>% .[,"variant_id"] %>% table %>% as.vector() %>% summary()
# ********** IFN **********
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   1.000   1.485   2.000  11.000 
# ********** LPS **********
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   1.000   1.455   2.000  11.000 
# ********** untreated **********
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1.0     1.0     1.0     1.5     2.0    11.0 


## Num of cis-eQTLs with min pval x gene
# nominal_eqtl_results %>% group_by(gene_symbol) %>% slice_min(pval_nominal) %>% .[, "gene_symbol"] %>% table %>% table
# Up: num of cis-eQTLs with min pval x gene
# Down: num of genes with such num of cis-eQTLs

# ********** IFN **********
#    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22   23   24   25   26 
# 7710 1062  478  290  202  133   94   77   58   35   46   40   33   36   26   20   17   15   16   10   10    7   13    6    8    6 
# 27   28   29   30   31   32   33   34   35   36   37   38   39   40   41   42   43   44   45   47   48   52   54   55   56   60 
#  3    4    7    4    8    4    4    4    2    2    3    1    3    2    2    1    2    2    2    2    3    1    1    2    2    1 
# 62   63   66   69   70   72   73   74   75   77   81   82   86   90   91   98   99  110  112  120  134  162  165 
#  1    1    1    2    2    2    1    1    2    1    1    1    1    1    1    2    1    2    1    1    1    1    1 
# ********** LPS **********
#    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22   23   24   25   26 
# 7755 1071  485  289  218  142   92   91   55   43   39   37   16   21   23   22   15   15   12    7   12    9    8    8   10    7 
# 27   28   29   30   31   32   33   34   35   36   37   38   39   42   43   44   45   48   50   51   53   56   58   60   62   63 
# 10    5    5    5    3    3    9    2    2    1    4    3    2    2    3    3    1    4    2    2    2    1    1    2    1    1 
# 66   67   69   70   72   74   75   76   77   86  112  134  181  218 
#  2    1    2    1    1    3    1    1    1    1    3    1    1    1 
# ********** untreated **********
#    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32   33   34   37 
# 8006 1072  490  279  210  143   99  101   57   35   42   42   30   27   29   17   19   22   15   10   13   12   10   10    4   10    6    5    4    2    5    1    7    5    4 
# 38   39   41   42   43   44   45   46   48   50   52   53   55   56   59   60   62   64   65   66   67   69   71   72   73   78   79   81   84   86   90   97  104  112  120 
#  2    1    1    2    3    2    2    2    4    1    1    1    1    1    2    1    2    1    1    1    2    1    2    1    1    1    2    1    2    1    1    1    1    2    1 
# 195 
#   1 

## Num of cis-eQTLs with min pval and diff slope x gene (i.e. independent cis-eQTL)
# min_eqtls <- nominal_eqtl_results %>% group_by(gene_symbol) %>% slice_min(pval_nominal)
# min_eqtls %>% group_by(gene_symbol) %>% dplyr::select(slope) %>% 
#   ungroup %>% table() %>% as.data.frame() %>% 
#   filter(Freq > 0 ) %>% 
#   dplyr::select(gene_symbol) %>% 
#   table %>% table
# ********** IFN **********
#     1     2 
# 10455    93 
# ********** LPS **********
#     1     2 
# 10541    59 
# ********** untreated **********
#     1     2 
# 10816    79 


###########################  Colocalization results  ########################### 

# Colocalization results
# coloc_res = readr::read_rds(
#   file = paste0(
#     "../../../OTAR2065_WGS_iPSC_premac_micro/data/results/",
#     "/3.coloc/coloc_results/my_GWAS_my_eQTL/coloc_results_single_causal_variant_500kb.rds"))

## Num of genes with coloc results in each treatment 
# lapply(coloc_res, length)
# ********** IFN **********
# [1] 4402
# ********** LPS **********
# [1] 4414
# ******* untreated *******
# [1] 5082


## Genes with signif coloc variants
# lapply(coloc_res, function(treatment){table(lapply(treatment, function(gene){gene$summary["PP.H4.abf"]}) > 0.7)})
# ********** IFN **********
# FALSE 
# 4402 
# ********** LPS **********
# FALSE  TRUE 
# 4411     3 (FRYL, ELP3 & SUSD3)
# ******* untreated *******
# FALSE  TRUE 
# 5080     2 (LY86 & OTOA)


##########################  GWAS summary statistics  ########################### 

# load(paste0("../", "../", "output_data/", "08_GWAS_proliferation/", "02_Check_association_results/", 
#             "summary_stats_microglia_vs_premac_", treatment, ".Rdata"), verbose = T)
# [1] "summary_stats"



## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#             9.1.1 Multivariable model definition per focal gene 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

## Function to find lead coloc variant per focal gene 
find_lead_coloc_variant <- function(focal_gene, treatment){
  
  coloc_res = coloc_res[[treatment]]
  
  ## Find top coloc variant for focal gene, if available
  if (sum(str_split_i(names(coloc_res), "_", i = 2) %in% focal_gene) == 0) {
    lead_coloc_var = NA # gene not tested
  }
  else{
    coloc_res = coloc_res[[names(coloc_res)[str_split_i(names(coloc_res), "_", i = 2) %in% focal_gene]]]
    
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
      gene1_eQTLs <- nominal_eqtl_results %>% dplyr::select(gene_symbol, variant_id, pval_nominal, slope) %>% filter(gene_symbol == gene1) %>% unique
      ## all cis-eQTLs for gene2
      gene2_eQTLs <- nominal_eqtl_results %>% dplyr::select(gene_symbol, variant_id, pval_nominal, slope) %>% filter(gene_symbol == gene2) %>% unique
      
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
select_IVs_and_exposures <- function(focal_gene, treatment){
  
  # ****************************************************************************
  #    Step 1:  Select causal + strongest independent cis-eQTLs of focal gene
  # ****************************************************************************
  
  # Select strongest cis-eQTL(s) x gene (with min pval)
  gene_cis_eQTLs <- nominal_eqtl_results %>% 
    dplyr::select(variant_id, gene_symbol, pval_nominal, slope) %>%
    dplyr::filter(gene_symbol == focal_gene) %>% 
    dplyr::slice_min(pval_nominal, na_rm = TRUE)  
  
  ## Remove variants with duplicated pval and effect size (SNPs in LD)
  if(length(which(duplicated(gene_cis_eQTLs$slope)))>0){
    gene_cis_eQTLs <- gene_cis_eQTLs[-which(duplicated(gene_cis_eQTLs$slope)),]
  }
  
  ## Select lead coloc cis-eQTL (if any)
  lead_coloc_var <- find_lead_coloc_variant(focal_gene, treatment)
  if(!is.na(lead_coloc_var)){
    gene_lead_coloc_cis_eQTL <- nominal_eqtl_results %>% 
      dplyr::select(variant_id, gene_symbol, pval_nominal, slope) %>%
      dplyr::filter(gene_symbol == focal_gene) %>% 
      dplyr::filter(variant_id == lead_coloc_var) 
    
    ## If lead coloc is in LD with any cis-eQTL already included, in step 1.1 the latter will be removed
    gene_cis_eQTLs <- rbind(gene_cis_eQTLs, gene_lead_coloc_cis_eQTL) %>% unique 
  }
  
  # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  #     Step 1.1: if >1 cis-eQTL keep the independent strongest ones 
  # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  if(length(gene_cis_eQTLs$variant_id) > 1){
    
    ## Variant pairs in high LD (R^2 >0.2)
    variant_pairs_in_LD = variants_in_LD(gene_cis_eQTLs$variant_id) 
    
    # 1) first discard all variants in LD with the lead coloc variant 
    # 2) then discard all variants in LD with the variant with the strongest effect size
    # 3) then repeat with the second strongest variant among the remaining and so on ...
    
    if(is.na(lead_coloc_var)){
      ordered_eQTLs <- gene_cis_eQTLs %>% arrange(-abs(slope)) %>% dplyr::select(variant_id) %>% unlist %>% as.vector()
    }
    else{
      ordered_eQTLs <- c(lead_coloc_var,
                         gene_cis_eQTLs %>% filter(!variant_id == lead_coloc_var) %>% arrange(-abs(slope)) %>% dplyr::select(variant_id) %>% unlist %>% as.vector())
    }
    
    independent_variants <- gene_cis_eQTLs$variant_id
    while(length(ordered_eQTLs)>0){
      
      ## Causal/strongest variant
      strong_eQTL = ordered_eQTLs[1]
      
      ## Correlated variants + strongest
      strong_and_cor_variants = variant_pairs_in_LD %>% 
        filter(SNP1 %in% strong_eQTL | SNP2 %in% strong_eQTL) %>% .[, c("SNP1", "SNP2")] %>% unlist() %>% append(strong_eQTL) %>% unique
      
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
  gene_indep_cis_eQTLs <- gene_cis_eQTLs %>% filter(variant_id %in% independent_variants)
  
  
  # ****************************************************************************
  #      Step 2: Add all eGenes of included independent stronger cis-eQTL(s)
  # ****************************************************************************
  indep_cis_eQTLs_eGenes <- nominal_eqtl_results %>% 
    dplyr::select(variant_id, gene_symbol, pval_nominal, slope) %>%
    filter(variant_id %in% independent_variants) 
  
  gene_indep_cis_eQTLs_with_eGenes = rbind(gene_indep_cis_eQTLs, indep_cis_eQTLs_eGenes) %>% unique
  
  # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  #          Step 2.1: exclude highly-correlated genes (r^2 >=0.4)
  # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  initial_genes_in_model <- gene_indep_cis_eQTLs_with_eGenes$gene_symbol %>% unique()
  
  if(length(initial_genes_in_model)>1){
    ## Compute cor between independent cis-eQTL effects for each pair of genes
    cors <- gene_pair_cor(initial_genes_in_model)
    # heatmap(cors, Rowv = NA, Colv = NA, scale = "none")
    
    ## First discard eGenes correlated with focal gene
    uncorr_eGenes <- initial_genes_in_model
    focal_gene_cor <- cors[focal_gene, ] %>% .[abs(.)>=0.4] %>% names %>% .[!.== focal_gene]
    if(any(focal_gene_cor %in% uncorr_eGenes)){
      uncorr_eGenes <- uncorr_eGenes[-which(uncorr_eGenes %in% focal_gene_cor)]
    }
    remaining_eGenes <- initial_genes_in_model[-which(initial_genes_in_model %in% c(focal_gene_cor, focal_gene))]
    
    ## Second, take less corr gene with focal gene and discard its correlated genes
    remaining_eGenes <- cors[focal_gene, ] %>% .[remaining_eGenes] %>% abs() %>% sort() %>% names
    
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
    uncorr_eGenes <- focal_gene
  }
  
  gene_indep_cis_eQTLs_with_indep_eGenes <- gene_indep_cis_eQTLs_with_eGenes %>% filter(gene_symbol %in% uncorr_eGenes)
  
  
  # ****************************************************************************
  #          Step 3: Add all signif (p < 10-6) cis-eQTLs x added eGene
  # ****************************************************************************
  indep_eGenes_signif_cis_eQTLs <- nominal_eqtl_results %>% 
    dplyr::select(variant_id, gene_symbol, pval_nominal, slope) %>%
    filter(gene_symbol %in% gene_indep_cis_eQTLs_with_indep_eGenes$gene_symbol & pval_nominal < 0.000001) %>%
    filter(! gene_symbol == focal_gene) # don't add cis-eQTLs for focal gene that were not added in step 1 !!!
  
  duplicates <- which(duplicated(indep_eGenes_signif_cis_eQTLs[, c("gene_symbol", "pval_nominal", "slope")]))
  
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
        filter(SNP1 %in% eQTL | SNP2 %in% eQTL) %>% .[, c("SNP1", "SNP2")] %>% unlist() %>% append(eQTL) %>% unique
      
      cor_variants = eQTL_and_cor_variants[-which(eQTL_and_cor_variants == eQTL)]
      if(any(cor_variants %in% independent_final_variants)){
        independent_final_variants <- independent_final_variants[-which(independent_final_variants %in% cor_variants)]
      }
      ## Repeat with remaining variants 
      ordered_final_variants <- ordered_final_variants[-which(ordered_final_variants %in% eQTL_and_cor_variants)]
    }
    
    ## Confirm all focal gene causal and strong cis-eQTLs are included 
    if(length(which(!gene_indep_cis_eQTLs$variant_id %in% independent_final_variants)) >0){
      message(paste0("Not all independent cis-eQTLs of focal gene ", focal_gene, " are included in the final model."))
      stop()
    }
    if((!is.na(lead_coloc_var)) & !(lead_coloc_var %in% independent_final_variants)){
      message(paste0("Lead coloc cis-eQTL of focal gene ", focal_gene, " not included in the final model."))
      stop()
    }
    
  }
  else{
    independent_final_variants <- gene_indep_cis_eQTLs$variant_id
  }
  
  
  ## Subset to final set of independent eQTLs
  gene_indep_cis_eQTLs_with_indep_eGenes_signif_cis_eQTLs <- rbind(gene_indep_cis_eQTLs_with_indep_eGenes, indep_eGenes_signif_cis_eQTLs) %>% unique
  gene_indep_cis_eQTLs_with_indep_eGenes_signif_cis_eQTLs = gene_indep_cis_eQTLs_with_indep_eGenes_signif_cis_eQTLs %>% 
    filter(variant_id %in% independent_final_variants)
  
  
  ## Final set of IVs and exposures
  final_model = gene_indep_cis_eQTLs_with_indep_eGenes_signif_cis_eQTLs %>% 
    dplyr::select(variant_id, gene_symbol, slope) %>%
    tidyr::pivot_wider(
      id_cols = variant_id,
      names_from = gene_symbol,
      values_from = slope,
      values_fill = 0) # 0 effect for eQTLs ns in other genes
  
  ## Confirm final set of uncorr genes (with focal gene) and independent cis-eQTLs
  if(length(which(!gene_indep_cis_eQTLs$variant_id %in% final_model$variant_id)) >0){
    message(paste0("Not all independent cis-eQTLs of focal gene ", focal_gene, " are included in the final model in ", treatment, "."))
    stop()
  }
  
  if(!focal_gene %in% colnames(final_model)){
    message(paste0("Focal gene ", focal_gene, " not included in the final model in ", treatment, "."))
    stop()
  }
  
  if(!setequal(colnames(final_model)[-1], uncorr_eGenes)){
    message(paste0("Not all uncorrelated genes for ", focal_gene, " included in the final model in ", treatment, "."))
    stop()
  }
  
  print(paste0("Final model for ", focal_gene, " in ", treatment, ": ", dim(final_model)[2]-2, " uncorrelated eGenes as exposures and ", 
               dim(final_model)[1], " independent cis-eQTLs as IVs"))
  
  final_model = final_model %>% column_to_rownames("variant_id")
  return(final_model)
  
}


## Test function with example focal genes to avoid headaches later :)
## (examples in untreated)

# * genes with 1 independent cis-eQTL that has only 1 eGene (the focal gene):
# E <- select_IVs_and_exposures(focal_gene = "VANGL1", treatment = "untreated")
# [1] "Final model for VANGL1 in untreated: 0 uncorrelated eGenes as exposures and 1 independent cis-eQTLs as IVs" 
# E <- select_IVs_and_exposures(focal_gene = "ROCK1", treatment = "untreated")
# [1] "Final model for ROCK1 in untreated: 0 uncorrelated eGenes as exposures and 1 independent cis-eQTLs as IVs"

# * genes with lead coloc variant as unique independent cis-eQTL: 
# E <- select_IVs_and_exposures(focal_gene = "OTOA", treatment = "untreated")
# [1] "Final model for OTOA in untreated: 1 uncorrelated eGenes as exposures and 1 independent cis-eQTLs as IVs"
# E <- select_IVs_and_exposures(focal_gene = "LY86", treatment = "untreated")
# [1] "Final model for LY86 in untreated: 0 uncorrelated eGenes as exposures and 1 independent cis-eQTLs as IVs"

# * gene with 1 independent cis-eQTL with multiple eGenes with signif eQTLs
# E <- select_IVs_and_exposures(focal_gene = "A2M", treatment = "untreated")
# [1] "Final model for A2M in untreated: 2 uncorrelated eGenes as exposures and 24 independent cis-eQTLs as IVs"

# * gene with 1 independent cis-eQTL with multiple eGenes with 0 signif eQTLs
# E <- select_IVs_and_exposures(focal_gene = "SAMD11", treatment = "untreated")
# [1] "Final model for SAMD11 in untreated: 4 uncorrelated eGenes as exposures and 1 independent cis-eQTLs as IVs"

# * gene with 2 independent cis-eQTLs with 1 eGene only (the focal gene)
# E <- select_IVs_and_exposures(focal_gene = "APBB2", treatment = "untreated")
# [1] "Final model for APBB2 in untreated: 0 uncorrelated eGenes as exposures and 2 independent cis-eQTLs as IVs"

# * gene with 2 (initial) independent cis-eQTLs with multiple eGenes with signif eQTLs
# E <- select_IVs_and_exposures(focal_gene = "BLOC1S3", treatment = "untreated")
# [1] "Final model for BLOC1S3 in untreated: 8 uncorrelated eGenes as exposures and 21 independent cis-eQTLs as IVs"

# * gene with 2 independent cis-eQTLs with multiple eGenes with 0 signif eQTLs
# E <- select_IVs_and_exposures(focal_gene = "BCKDK", treatment = "untreated")
# [1] "Final model for BCKDK in untreated: 2 uncorrelated eGenes as exposures and 2 independent cis-eQTLs as IVs"



## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#                    9.1.2 Prepare input matrices for TWMR
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

for(treatment in c("IFN", "LPS", "untreated")){
  
  ## Define output dir
  outdir = paste("..", "..", "output_data", "09_TWMR_proliferation", "01_IVs_and_exposures_selection", treatment, sep = "/")
  dir.create(outdir, recursive = T, showWarnings = F)

  ## Read cis-eQTL summary stats x treatment 
  nominal_eqtl_results = readr::read_delim(paste0(input_dir, "/nominal_eqtl_results_", treatment, ".csv")) %>% 
    filter(! is.na(gene_symbol))
  
  ## Read GWAS summary stats x treatment 
  load(paste0("../", "../", "output_data/", "08_GWAS_proliferation/", "02_Check_association_results/", 
              "summary_stats_microglia_vs_premac_", treatment, ".Rdata"), verbose = T)
  # [1] "summary_stats"
  
  ## Coloc results
  coloc_res = readr::read_rds(
    file = paste0(
      "../../../OTAR2065_WGS_iPSC_premac_micro/data/results/",
      "/3.coloc/coloc_results/my_GWAS_my_eQTL/coloc_results_single_causal_variant_500kb.rds"))
  
  ## Split genes in 1000 chunks (10-11 genes each)
  all_genes <- unique(nominal_eqtl_results$gene_symbol) 
  all_genes_chunks <- split(all_genes, cut(seq_along(all_genes), 1000, labels = FALSE))
  focal_genes <- all_genes_chunks[job_i] %>% unlist
  
  ## Find model x focal gene 
  for(focal_gene in focal_genes){
    
    ## Ignore genes that have been read already
    if(focal_gene %in% str_split_fixed(list.files(outdir), "_", 2)[,1]){
      message("Gene already read.")
      next()
    }
    
    ## E: IVs x eGenes
    E = select_IVs_and_exposures(focal_gene, treatment) 
    
    ## * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    #  Warning 1:                   !!!                            *
    #  Effect of variants with ns/non-tested effects on any of the *
    #  eGenes included in E are assumed to be 0.                   *                              
    ## * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    
    ## Make sure focal gene is the 1st gene column in E
    if(which(colnames(E) == focal_gene) != 1){
      eGenes <- c(colnames(E), colnames(E)[1])
      eGenes[1] <- focal_gene
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
        filter(
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
      message(paste0("Not all IVs for ", focal_gene, " are present in C in ", treatment, "."))
      stop()
    }
    if(!isSymmetric(C)){
      message(paste0("C for ", focal_gene, " is not symmetric in ", treatment, "."))
      stop()
    }
    if(any(! diag(C) == 1)){
      message(paste0("C for ", focal_gene, " has not all 1s on the diagonal in ", treatment, "."))
      stop()
    } 
   
    ## IVs effects on proliferation from preMac -> microglia in treatment (G: IVs)
    G = summary_stats[rownames(E), "Estimate"] %>% 
      tidyr::replace_na(0) # NA variant effects are 0
    names(G) = rownames(E)
    
    ## * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    #  Warning 3:                   !!!                            *
    #  Effect of variants with no GWAS effect available are set to *
    #  0.                                                          *
    ## * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    
    ## Sample sizes for eQTL, GWAS, and LD  
    N_eQTLs = 261
    Ngwas <- case_when(
      treatment == "IFN" ~ 230,
      treatment == "LPS" ~ 229,
      treatment == "untreated" ~ 209)
    n_for_LD = 13160 # from TopMED README
    
    inputs <- list("E" = E, "G" = G, "C" = C, 
                   "N_eQTLs" = N_eQTLs, "Ngwas" = Ngwas, "n_for_LD" = n_for_LD)
    save(inputs, file = paste0(outdir, "/", focal_gene, "_input_matrices.Rdata"))
     
   }
}

## Try example genes (in untreated):
# * focal gene with 1 IV only: CLPTM1
# * focal gene with IVs not in LD file: A2M







# Reproducibility info
# session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.1 (2023-06-16)
# os       Ubuntu 22.04.4 LTS
# system   x86_64, linux-gnu
# ui       RStudio
# language (EN)
# collate  en_GB.UTF-8
# ctype    en_GB.UTF-8
# tz       Europe/Belfast
# date     2025-03-18
# rstudio  2024.04.0+735 Chocolate Cosmos (server)
# pandoc   3.1.12.3 @ /opt/view/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version     date (UTC) lib source
# abind                  1.4-5       2016-07-21 [1] CRAN (R 4.3.1)
# AnnotationDbi          1.64.1      2023-11-03 [1] Bioconductor
# AnnotationFilter       1.26.0      2023-10-24 [1] Bioconductor
# ape                    5.7-1       2023-03-13 [1] CRAN (R 4.3.1)
# backports              1.4.1       2021-12-13 [1] CRAN (R 4.3.1)
# beachmat               2.18.1      2024-02-14 [1] Bioconductor 3.18 (R 4.3.1)
# Biobase                2.62.0      2023-10-24 [1] Bioconductor
# BiocFileCache          2.10.2      2024-03-27 [1] Bioconductor 3.18 (R 4.3.1)
# BiocGenerics           0.48.1      2023-11-01 [1] Bioconductor
# BiocIO                 1.12.0      2023-10-24 [1] Bioconductor
# BiocParallel           1.36.0      2023-10-24 [1] Bioconductor
# BiocSingular           1.18.0      2023-10-24 [1] Bioconductor
# biomaRt              * 2.58.2      2024-01-30 [1] Bioconductor 3.18 (R 4.3.1)
# Biostrings             2.70.3      2024-03-13 [1] Bioconductor 3.18 (R 4.3.1)
# bit                    4.0.5       2022-11-15 [1] CRAN (R 4.3.1)
# bit64                  4.0.5       2020-08-30 [1] CRAN (R 4.3.1)
# bitops                 1.0-7       2021-04-24 [1] CRAN (R 4.3.1)
# blob                   1.2.4       2023-03-17 [1] CRAN (R 4.3.1)
# broom                  1.0.5       2023-06-09 [1] CRAN (R 4.3.1)
# cachem                 1.0.8       2023-05-01 [1] CRAN (R 4.3.1)
# car                    3.1-2       2023-03-30 [1] CRAN (R 4.3.1)
# carData                3.0-5       2022-01-06 [1] CRAN (R 4.3.1)
# circlize               0.4.16      2024-02-20 [1] CRAN (R 4.3.1)
# cli                    3.6.2       2023-12-11 [1] CRAN (R 4.3.1)
# cluster                2.1.6       2023-12-01 [1] CRAN (R 4.3.1)
# codetools              0.2-20      2024-03-31 [1] CRAN (R 4.3.1)
# coloc                  5.2.3       2023-10-03 [1] CRAN (R 4.3.1)
# colorspace             2.1-0       2023-01-23 [1] CRAN (R 4.3.1)
# cowplot              * 1.1.3       2024-01-22 [1] CRAN (R 4.3.1)
# crayon                 1.5.2       2022-09-29 [1] CRAN (R 4.3.1)
# curl                   5.2.1       2024-03-01 [1] CRAN (R 4.3.1)
# data.table             1.15.4      2024-03-30 [1] CRAN (R 4.3.1)
# DBI                    1.2.2       2024-02-16 [1] CRAN (R 4.3.1)
# dbplyr                 2.5.0       2024-03-19 [1] CRAN (R 4.3.1)
# DelayedArray           0.28.0      2023-10-24 [1] Bioconductor
# digest                 0.6.35      2024-03-11 [1] CRAN (R 4.3.1)
# dplyr                * 1.1.4       2023-11-17 [1] CRAN (R 4.3.1)
# edgeR                  4.0.16      2024-02-18 [1] Bioconductor 3.18 (R 4.3.1)
# ensembldb              2.26.0      2023-10-24 [1] Bioconductor
# fansi                  1.0.6       2023-12-08 [1] CRAN (R 4.3.1)
# fastmap                1.1.1       2023-02-24 [1] CRAN (R 4.3.1)
# filelock               1.0.3       2023-12-11 [1] CRAN (R 4.3.1)
# forcats              * 1.0.0       2023-01-29 [1] CRAN (R 4.3.1)
# generics               0.1.3       2022-07-05 [1] CRAN (R 4.3.1)
# GenomeInfoDb           1.38.8      2024-03-15 [1] Bioconductor 3.18 (R 4.3.1)
# GenomeInfoDbData       1.2.11      2025-02-07 [1] Bioconductor
# GenomicAlignments      1.38.2      2024-01-16 [1] Bioconductor 3.18 (R 4.3.1)
# GenomicFeatures        1.54.4      2024-03-13 [1] Bioconductor 3.18 (R 4.3.1)
# GenomicRanges          1.54.1      2023-10-29 [1] Bioconductor
# gggrid                 0.2-0       2022-01-11 [1] CRAN (R 4.3.1)
# ggplot2              * 3.5.1       2024-04-23 [1] CRAN (R 4.3.1)
# ggpubr                 0.6.0       2023-02-10 [1] CRAN (R 4.3.1)
# ggrepel                0.9.5       2024-01-10 [1] CRAN (R 4.3.1)
# ggsignif               0.6.4       2022-10-13 [1] CRAN (R 4.3.1)
# GlobalOptions          0.1.2       2020-06-10 [1] CRAN (R 4.3.1)
# glue                   1.7.0       2024-01-09 [1] CRAN (R 4.3.1)
# gridExtra              2.3         2017-09-09 [1] CRAN (R 4.3.1)
# gtable                 0.3.4       2023-08-21 [1] CRAN (R 4.3.1)
# here                   1.0.1       2020-12-13 [1] CRAN (R 4.3.1)
# hms                    1.1.3       2023-03-21 [1] CRAN (R 4.3.1)
# htmltools              0.5.8       2024-03-25 [1] CRAN (R 4.3.1)
# htmlwidgets            1.6.4       2023-12-06 [1] CRAN (R 4.3.1)
# httr                 * 1.4.7       2023-08-15 [1] CRAN (R 4.3.1)
# IRanges                2.36.0      2023-10-24 [1] Bioconductor
# irlba                  2.3.5.1     2022-10-03 [1] CRAN (R 4.3.1)
# jsonlite               1.8.8       2023-12-04 [1] CRAN (R 4.3.1)
# KEGGREST               1.42.0      2023-10-24 [1] Bioconductor
# lattice                0.22-6      2024-03-20 [1] CRAN (R 4.3.1)
# lazyeval               0.2.2       2019-03-15 [1] CRAN (R 4.3.1)
# LDlinkR                1.4.0       2024-04-10 [1] CRAN (R 4.3.1)
# lifecycle              1.0.4       2023-11-07 [1] CRAN (R 4.3.1)
# limma                  3.58.1      2023-10-31 [1] Bioconductor
# locfit                 1.5-9.9     2024-03-01 [1] CRAN (R 4.3.1)
# locuszoomr             0.2.1       2024-02-17 [1] CRAN (R 4.3.1)
# lubridate            * 1.9.3       2023-09-27 [1] CRAN (R 4.3.1)
# magrittr               2.0.3       2022-03-30 [1] CRAN (R 4.3.1)
# MASS                   7.3-60.0.1  2024-01-13 [1] CRAN (R 4.3.1)
# Matrix                 1.6-5       2024-01-11 [1] CRAN (R 4.3.1)
# MatrixGenerics         1.14.0      2023-10-24 [1] Bioconductor
# matrixStats            1.2.0       2023-12-11 [1] CRAN (R 4.3.1)
# memoise                2.0.1       2021-11-26 [1] CRAN (R 4.3.1)
# mgcv                   1.9-1       2023-12-21 [1] CRAN (R 4.3.1)
# mixsqp                 0.3-54      2023-12-20 [1] CRAN (R 4.3.1)
# munsell                0.5.1       2024-04-01 [1] CRAN (R 4.3.1)
# nlme                   3.1-164     2023-11-27 [1] CRAN (R 4.3.1)
# patchwork              1.2.0       2024-01-08 [1] CRAN (R 4.3.1)
# permute                0.9-7       2022-01-27 [1] CRAN (R 4.3.1)
# pillar                 1.9.0       2023-03-22 [1] CRAN (R 4.3.1)
# pinfsc50               1.3.0       2023-12-05 [1] CRAN (R 4.3.1)
# pkgconfig              2.0.3       2019-09-22 [1] CRAN (R 4.3.1)
# pkgload                1.3.4       2024-01-16 [1] CRAN (R 4.3.1)
# plotly                 4.10.4      2024-01-13 [1] CRAN (R 4.3.1)
# plyr                   1.8.9       2023-10-02 [1] CRAN (R 4.3.1)
# png                    0.1-8       2022-11-29 [1] CRAN (R 4.3.1)
# prettyunits            1.2.0       2023-09-24 [1] CRAN (R 4.3.1)
# progress               1.2.3       2023-12-06 [1] CRAN (R 4.3.1)
# ProtGenerics           1.34.0      2023-10-24 [1] Bioconductor
# purrr                * 1.0.2       2023-08-10 [1] CRAN (R 4.3.1)
# R6                     2.5.1       2021-08-19 [1] CRAN (R 4.3.1)
# rappdirs               0.3.3       2021-01-31 [1] CRAN (R 4.3.1)
# Rcpp                   1.0.12      2024-01-09 [1] CRAN (R 4.3.1)
# RCurl                  1.98-1.14   2024-01-09 [1] CRAN (R 4.3.1)
# readr                * 2.1.5       2024-01-10 [1] CRAN (R 4.3.1)
# reshape                0.8.9       2022-04-12 [1] CRAN (R 4.3.1)
# restfulr               0.0.15      2022-06-16 [1] CRAN (R 4.3.1)
# reticulate             1.35.0      2024-01-31 [1] CRAN (R 4.3.1)
# rjson                  0.2.21      2022-01-09 [1] CRAN (R 4.3.1)
# rlang                  1.1.3       2024-01-10 [1] CRAN (R 4.3.1)
# rprojroot              2.0.4       2023-11-05 [1] CRAN (R 4.3.1)
# Rsamtools              2.18.0      2023-10-24 [1] Bioconductor
# RSQLite                2.3.6       2024-03-31 [1] CRAN (R 4.3.1)
# rstatix                0.7.2       2023-02-01 [1] CRAN (R 4.3.1)
# rstudioapi             0.16.0      2024-03-24 [1] CRAN (R 4.3.1)
# rsvd                   1.0.5       2021-04-16 [1] CRAN (R 4.3.1)
# rtracklayer            1.62.0      2023-10-24 [1] Bioconductor
# S4Arrays               1.2.1       2024-03-04 [1] Bioconductor 3.18 (R 4.3.1)
# S4Vectors              0.40.2      2023-11-23 [1] Bioconductor 3.18 (R 4.3.1)
# ScaledMatrix           1.10.0      2023-10-24 [1] Bioconductor
# scales                 1.3.0       2023-11-28 [1] CRAN (R 4.3.1)
# sessioninfo          * 1.2.2       2021-12-06 [1] CRAN (R 4.3.1)
# shape                  1.4.6.1     2024-02-23 [1] CRAN (R 4.3.1)
# SparseArray            1.2.4       2024-02-11 [1] Bioconductor 3.18 (R 4.3.1)
# statmod                1.5.0       2023-01-06 [1] CRAN (R 4.3.1)
# stringi                1.8.3       2023-12-11 [1] CRAN (R 4.3.1)
# stringr              * 1.5.1       2023-11-14 [1] CRAN (R 4.3.1)
# SummarizedExperiment   1.32.0      2023-10-24 [1] Bioconductor
# susieR                 0.12.35     2023-02-17 [1] CRAN (R 4.3.1)
# tibble               * 3.2.1       2023-03-20 [1] CRAN (R 4.3.1)
# tidyr                * 1.3.1       2024-01-24 [1] CRAN (R 4.3.1)
# tidyselect             1.2.1       2024-03-11 [1] CRAN (R 4.3.1)
# tidyverse            * 2.0.0       2023-02-22 [1] CRAN (R 4.3.1)
# timechange             0.3.0       2024-01-18 [1] CRAN (R 4.3.1)
# tzdb                   0.4.0       2023-05-12 [1] CRAN (R 4.3.1)
# utf8                   1.2.4       2023-10-22 [1] CRAN (R 4.3.1)
# vcfR                   1.15.0      2023-12-08 [1] CRAN (R 4.3.1)
# vctrs                  0.6.5       2023-12-01 [1] CRAN (R 4.3.1)
# vegan                  2.6-4       2022-10-11 [1] CRAN (R 4.3.1)
# viridis                0.6.5       2024-01-29 [1] CRAN (R 4.3.1)
# viridisLite            0.4.2       2023-05-02 [1] CRAN (R 4.3.1)
# vroom                  1.6.5       2023-12-05 [1] CRAN (R 4.3.1)
# withr                  3.0.0       2024-01-16 [1] CRAN (R 4.3.1)
# XML                    3.99-0.16.1 2024-01-22 [1] CRAN (R 4.3.1)
# xml2                   1.3.6       2023-12-04 [1] CRAN (R 4.3.1)
# XVector                0.42.0      2023-10-24 [1] Bioconductor
# yaml                   2.3.8       2023-12-11 [1] CRAN (R 4.3.1)
# zlibbioc               1.48.2      2024-03-13 [1] Bioconductor 3.18 (R 4.3.1)
# zoo                    1.8-12      2023-04-13 [1] CRAN (R 4.3.1)
# 
# [1] /opt/view/rlib/R/library
# [2] /software/hgi/installs/softpack/rstudio/.spack/linux-ubuntu22.04-x86_64_v3/gcc-11.4.0/r-4.3.1-bfwldrk76z6f52upk47zepliekn7ayqz/rlib/R/library
# 
# ─ Python configuration ───────────────────────────────────────────────────────────────────────────────────────────────
# python:         /opt/view/bin/python3
# libpython:      /opt/software/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_placeholder__/__spac/linux-ubuntu22.04-x86_64_v3/gcc-11.4.0/python-3.10.10-xe7qdrd7kcyaemvzf72gss3fwgtnn724/lib/libpython3.10.so
# pythonhome:     /opt/._view/ytxdco23iou6p3zeecg3ovplyjmzzkur:/opt/._view/ytxdco23iou6p3zeecg3ovplyjmzzkur
# version:        3.10.10 (main, Nov  4 2024, 10:12:22) [GCC 11.4.0]
# numpy:          /opt/._view/ytxdco23iou6p3zeecg3ovplyjmzzkur/lib/python3.10/site-packages/numpy
# numpy_version:  1.26.4
# keyword:        /opt/._view/ytxdco23iou6p3zeecg3ovplyjmzzkur/lib/python3.10
# 
# NOTE: Python version was forced by PATH
# 
# python versions found: 
#   /opt/view/bin/python3
# /opt/view/bin/python
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

