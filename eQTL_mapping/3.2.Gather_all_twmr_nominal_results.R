# gathering the nominal twmr trans results in rds objects that contain the full information per gene
# table number 3 (af_s) has been transposed in bash first
library(tidyverse)
library(qs)
options(future.globals.maxSize = 40000 * 1024^2) # 40Gb

treatments = c("LPS","IFN","untreated")

# outputs are
# 0:pval_df
# 1: b_df
# 2:b_se_df
# 3:af_s

for(treat in treatments){
  trans_nominal_res = list()
  
  for(n in c(3,0,1,2)){
    # starting with 3 to get the names of the variants
    if(n==3){
      # this has been transposed in bash
      trans_nominal_res[[paste0(treat,"_element_",n)]] = read_csv(paste0("../../data/results/tensorqtl/best_results/trans_twmr/sum_sizefactorsNorm_log2_scaled_centered_",
                                                                         treat,
                                                                         "_Not_proliferating_common_tensorQTL_nominal_trans_qtl_pairs_element_",n,".csv"), col_names = FALSE)
      colnames( trans_nominal_res[[paste0(treat,"_element_",n)]]) = c("variant_id","effect_AF")
     trans_nominal_res[[paste0(treat,"_element_",n)]] = trans_nominal_res[[paste0(treat,"_element_",n)]] %>%
       tidyr::separate(variant_id, into = c("CHR", "POS", "REF", "ALT"), sep = "_",remove = FALSE)
      
    }
    else{
      trans_nominal_res[[paste0(treat,"_element_",n)]] = readr::read_csv(paste0("../../data/results/tensorqtl/best_results/trans_twmr/sum_sizefactorsNorm_log2_scaled_centered_",
                                                                                treat,
                                                                                "_Not_proliferating_common_tensorQTL_nominal_trans_qtl_pairs_element_",n,".csv"),
                                                                         col_names = TRUE)
     if(n==2){
       # get gene names
       all_genes = names(lapply(trans_nominal_res[[paste0(treat,"_element_",n)]], colnames))
       
     }
      }
  }
  # reshaping
  # For each gene, extract the corresponding columns across dataframes
  gene_dfs = lapply(all_genes, function(gene) {
    # For each dataframe, extract the gene column
    gene_columns <- lapply(trans_nominal_res, function(df) {
      if (gene %in% colnames(df)) {
        df[[gene]]
      } else {
        NA  # If the gene is not found in a particular df, fill with NAs
      }
    })
    # Combine the columns into a new dataframe
    gene_df = as.data.frame(gene_columns)
    colnames(gene_df) = paste0("metric_", names(trans_nominal_res))  
    gene_df = gene_df %>%
      dplyr::select(!ends_with("_element_3")) %>%
      dplyr::rename(nominal_pvalue = ends_with("_element_0"),
                    beta = ends_with("_element_1"),
                    beta_se = ends_with("_element_2"))
    gene_df = cbind(gene_df,trans_nominal_res[[paste0(treat,"_element_",3)]])
    return(gene_df)
  })
  names(gene_dfs) = all_genes

  qsave(gene_dfs,paste0("../../data/results/tensorqtl/best_results/trans_twmr/formated_",
                          treat,
                          "_Not_proliferating_common_tensorQTL_nominal_trans_qtl_pairs_twmr.qs"))
  
  rm(trans_nominal_res)
}
