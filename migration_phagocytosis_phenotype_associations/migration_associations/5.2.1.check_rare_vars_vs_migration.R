# examine proliferation burden tests
library(patchwork)
library(tidyverse)
library(qvalue)
source("../functions.R")
directory = "../../../data/results/migration/5.2.1.check_rare_vars_vs_migration/"
dir.create(directory, recursive = TRUE)

options(future.globals.maxSize = 240000 * 1024^2) # 240Gb

skat_res = read_csv(paste0("../../../data/results/migration/5.1.check_SKAT_WES_vs_migration/deleterious_SKATO_scaled_centered_prop_pvals.csv"))%>%
  tidyr::drop_na(resampling_pval)


for(mut_type in c("ptv","syn","del")){
  res_per_treat = list()
  for(treat in c("untreated","IFN","LPS")){
  res = list()
  files = list.files(paste0("../../../data/results/migration/5.2.rare_vars_vs_migration/",mut_type,"/"), pattern = paste0("*",treat,"_burden_lm_res.rds"))
  files = str_sort(files, numeric = TRUE) # to later recover right plots by number
  numbers = str_split_i(files,pattern = "_",i = 1) #they may not be contiguous
  res = lapply(paste0("../../../data/results/migration/5.2.rare_vars_vs_migration/",mut_type,"/",files), function(x){
    read_rds(x) 
  })
  
  res_pvals <- purrr::map(res[[1]], ~extract_burden_results(.))
  
  res_pvals = do.call("rbind",res_pvals)
  
  res_per_treat[[treat]] = res_pvals %>%
    tidyr::pivot_wider(id_cols = gene, names_from = "coefficients", values_from = "value")
  
  res_per_treat[[treat]]$p_Bonf = p.adjust(res_per_treat[[treat]]$p_val,method = "bonferroni")
  res_per_treat[[treat]]$p_BH = p.adjust(res_per_treat[[treat]]$p_val,method = "BH")
  res_per_treat[[treat]]$treatment = treat
  res_per_treat[[treat]]$type = mut_type
  
  res_per_treat[[treat]] = res_per_treat[[treat]] %>%
    dplyr::left_join(skat_res[skat_res$treatment == treat,c("gene_name","resampling_pval")],
                     by = join_by(gene == gene_name)) %>%
    dplyr::rename(resampling_pval_skato = resampling_pval  )
  
  res_per_treat[[treat]] = res_per_treat[[treat]] %>%
    dplyr::arrange(p_Bonf)
  
  }
  res = do.call("rbind",res_per_treat)
  

  write_csv(res,paste0(directory,mut_type,"_burden_scaled_centered_prop_pvals.csv"))
}

table(res$p_Bonf<0.05,res$treatment)
# IFN LPS untreated
# FALSE 100 169       126
# TRUE    2   0         1

table(res$p_BH<0.05,res$treatment)
# IFN LPS untreated
# FALSE  99 169       124
# TRUE    3   0         3

table(res$resampling_pval_skato<0.05,res$treatment)
# IFN LPS untreated
# TRUE 102 169       127     # all true as expected

hist(res$p_val) # slightly left-leaning because they are pre-selected from significant SKAT
hist(res$p_BH)
hist(res$p_Bonf) # maybe better to use p-value from SKAT-O because it does permutations?

for(mut_type in c("ptv","del")){
  for(treat in c("untreated","IFN","LPS")){
    
    res = read_csv(paste0(directory,mut_type,"_burden_scaled_centered_prop_pvals.csv")) %>%
      dplyr::filter(treatment == treat & type == mut_type ) %>%
      dplyr::arrange(resampling_pval_skato)
    
    plist = read_rds(file = paste0("../../../data/results/migration/5.2.rare_vars_vs_migration/",
                                   mut_type,"/",treat,"partial_residual_LMM_plots.rds"))
    
    pdf(paste0(directory,mut_type,"_",treat, "_sign_burden_scaled_centered.pdf"),
        width = 4,height = 5)
    
    for(gene in res$gene){
      plot(plist[[paste(treat,gene,mut_type,sep = "_")]])
    }
    dev.off()
  }
  
  
}
