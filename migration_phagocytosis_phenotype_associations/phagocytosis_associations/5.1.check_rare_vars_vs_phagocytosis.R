# examine proliferation burden tests
library(patchwork)
library(tidyverse)
library(qvalue)
source("./functions.R")
directory = "../../../data/results/phagocytosis/5.1.check_rare_vars_vs_phagocytosis/"
dir.create(directory, recursive = TRUE)

options(future.globals.maxSize = 100000 * 1024^2) # 100Gb


for(mut_type in c("del")){
  res_per_treat = list()
  for(treat in c("untreated","IFN","LPS")){
  res = list()
  files = list.files(paste0("../../../data/results/phagocytosis/5.rare_vars_vs_phagocytosis/",mut_type,"/"), pattern = paste0("*",treat,"_burden_lm_res.rds"))
  files = str_sort(files, numeric = TRUE) # to later recover right plots by number
  numbers = str_split_i(files,pattern = "_",i = 1) #they may not be contiguous
  res = lapply(paste0("../../../data/results/phagocytosis/5.rare_vars_vs_phagocytosis/",mut_type,"/",files), function(x){
    read_rds(x) 
  })
  
  res_pvals <- purrr::map(res, ~extract_burden_results(.))

  res_pvals = do.call("rbind",res_pvals)

  res_per_treat[[treat]] = res_pvals %>%
    tidyr::pivot_wider(id_cols = gene, names_from = "coefficients", values_from = "value")
    
  res_per_treat[[treat]]$p_Bonf = p.adjust(res_per_treat[[treat]]$p_val,method = "bonferroni")
  res_per_treat[[treat]]$treatment = treat
  res_per_treat[[treat]]$number = numbers

  res_per_treat[[treat]] = res_per_treat[[treat]] %>%
    dplyr::arrange(p_Bonf)

  }
  res = do.call("rbind",res_per_treat)
  hist(res$p_val) # why are p-values so suspiciously significant all over?
  hist(res$p_Bonf) 
  write_csv(res,paste0(directory,mut_type,"_burden_scaled_centered_prop_pvals.csv"))
}


table(res$p_Bonf<0.05,res$treatment)
# IFN  LPS untreated
# FALSE 5256 7703      7656
# TRUE   143   20        67

### plot all significant genes

for(treat in c("untreated","IFN","LPS")){
  gene_df= res %>%
    dplyr::filter(treatment == treat & p_Bonf < 0.05) %>%
    dplyr::select(gene,number) %>%
    distinct()
  
  pdf(paste0(directory,"del_",treat,"_burden_scaled_centered_prop_pvals.pdf"),
      width = 4, height = 5)
for(number in gene_df$number){
  plots =   read_rds( paste0("../../../data/results/phagocytosis/5.rare_vars_vs_phagocytosis/del/", number,"_",treat,"partial_residual_LMM_plots.rds"))
  
  plot(plots[[1]])
 
}
  
  dev.off()
}

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
