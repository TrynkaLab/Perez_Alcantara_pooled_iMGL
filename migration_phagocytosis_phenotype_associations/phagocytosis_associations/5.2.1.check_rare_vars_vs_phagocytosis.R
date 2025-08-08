# examine proliferation burden tests
library(patchwork)
library(tidyverse)
library(qvalue)
source("../functions.R")
directory = "../../../data/results/phagocytosis/5.2.check_rare_vars_vs_phagocytosis/"
dir.create(directory, recursive = TRUE)

skat_res = read_csv(paste0("../../../data/results/phagocytosis/5.1.check_SKAT_WES_vs_phagocytosis/deleterious_SKATO_scaled_centered_prop_pvals.csv")) %>%
  tidyr::drop_na(resampling_pval)

# options(future.globals.maxSize = 240000 * 1024^2) # 240Gb
options(future.globals.maxSize = 35000 * 1024^2) # 35Gb


for(mut_type in c("ptv","syn","del")){
  res_per_treat = list()
  for(treat in c("untreated","IFN","LPS")){
  res = list()
  files = list.files(paste0("../../../data/results/phagocytosis/5.2.rare_vars_vs_phagocytosis/",mut_type,"/"), pattern = paste0("*",treat,"_burden_lm_res.rds"))
  files = str_sort(files, numeric = TRUE) # to later recover right plots by number
  res = read_rds(paste0("../../../data/results/phagocytosis/5.2.rare_vars_vs_phagocytosis/",mut_type,"/",files))
  
  res_pvals <- purrr::map(res, ~extract_burden_results(.))

  res_pvals = do.call("rbind",res_pvals)

  res_per_treat[[treat]] = res_pvals %>%
    tidyr::pivot_wider(id_cols = gene, names_from = "coefficients", values_from = "value")
    
  res_per_treat[[treat]]$p_Bonf = p.adjust(res_per_treat[[treat]]$p_val,method = "bonferroni")
  res_per_treat[[treat]]$p_BH = p.adjust(res_per_treat[[treat]]$p_val,method = "BH")
  res_per_treat[[treat]]$treatment = treat
  res_per_treat[[treat]]$type = mut_type

  res_per_treat[[treat]] = res_per_treat[[treat]] %>%
    dplyr::inner_join(skat_res[skat_res$treatment == treat,c("gene_name","resampling_pval")],
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
# FALSE 177 120       129
# TRUE    3   2         7
table(res$p_BH<0.05,res$treatment)
# IFN LPS untreated
# FALSE 161 110       121
# TRUE   19  12        15
table(res$resampling_pval_skato<0.05,res$treatment)
# IFN LPS untreated
# TRUE 180 122       136

### plot all significant genes
# from the SKAT analysis
hist(res$p_val) # slightly left-leaning because they are pre-selected from significant SKAT
hist(res$p_BH)
hist(res$p_Bonf) # maybe better to use p-value from SKAT-O because it does permutations?

for(mut_type in c("ptv","del")){
  for(treat in c("untreated","IFN","LPS")){

    res = read_csv(paste0(directory,mut_type,"_burden_scaled_centered_prop_pvals.csv")) %>%
      dplyr::filter(treatment == treat & type == mut_type ) %>%
      dplyr::arrange(resampling_pval_skato)
    
    plist = read_rds(file = paste0("../../../data/results/phagocytosis/5.2.rare_vars_vs_phagocytosis/",
                                   mut_type,"/",treat,"partial_residual_LMM_plots.rds"))
      
    pdf(paste0(directory,mut_type,"_",treat, "_sign_burden_scaled_centered.pdf"),
        width = 4,height = 5)
    
    for(gene in res$gene){
    plot(plist[[paste(treat,gene,mut_type,sep = "_")]])
    }
    dev.off()
    }
  
  
  }


# skat results
skat_res = read_csv("../../../data/results/phagocytosis/5.1.check_SKAT_WES_vs_phagocytosis/deleterious_SKATO_scaled_centered_prop_pvals.csv") %>%
  dplyr::filter(treatment == "untreated") %>%
  tidyr::drop_na(resampling_pval)
# any in Sam's list?
sam_phago = read_tsv("../../../../CRISPR/OTAR2065_phagocytosis_CRISPR/data/2024_07_sam_screen_results.tsv") %>%
  dplyr::arrange(FDR) %>%
  dplyr::rename(gene = id)
directionality = sam_phago %>%
    dplyr::left_join(res[res$treatment=="untreated",c("gene","estimate")]) %>%
  dplyr::left_join(skat_res[skat_res$treatment=="untreated",c("gene_name","resampling_pval")], 
                    by = join_by(gene == gene_name)) %>%
  dplyr::mutate(direction = case_when((estimate < 0 & Score > 0 & resampling_pval < 0.05) | (estimate > 0 & Score < 0 & resampling_pval < 0.05) ~ "same direction",
                                      (estimate < 0 & Score < 0 & resampling_pval < 0.05) | (estimate > 0 & Score > 0 & resampling_pval < 0.05)  ~ "opposite direction",
                                      .default = "NS SKAT")) 


# INPP5D and TRAPPC5 significant in SKAT, burden shows same direction
# ANK3 not signficant, but burden in same direction
# UBA1 not signfiicant, burden in opposite direction
# the rest are not tested
write_csv(directionality,paste0("../../../data/results/phagocytosis/5.2.check_rare_vars_vs_phagocytosis/untreated_sam_CRISPR_directionality_SKAT.csv") )
