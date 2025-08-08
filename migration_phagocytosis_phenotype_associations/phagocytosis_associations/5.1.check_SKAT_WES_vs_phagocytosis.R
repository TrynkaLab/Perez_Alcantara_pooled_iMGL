# examine phagocytosis SKATO tests
library(patchwork)
library(tidyverse)
treatment_cols =  c(untreated = "#8D918B", IFN = "#3A5683", LPS = "#F8766D")
source("./functions.R")
directory = "../../../data/results/phagocytosis/5.1.check_SKAT_WES_vs_phagocytosis/"
dir.create(directory, recursive = TRUE)

options(future.globals.maxSize = 240000 * 1024^2) # 240Gb


for(mut_type in c("missense_non_deleterious","deleterious")){
  res_per_treat = list()
  for(treat in c("untreated","IFN","LPS")){
  res = list()
  files = list.files("../../../data/results/phagocytosis/5.rare_var_SKAT_WES/",
                            pattern = paste0("[0-9]","_",treat,"_",mut_type,"_SKATO_result.RDS"))
  files = str_sort(files, numeric = TRUE) # to later recover right plots by number
  numbers = str_split_i(files,pattern = "_",i = 1) #they may not be contiguous
  res = lapply(paste0("../../../data/results/phagocytosis/5.rare_var_SKAT_WES/",files), function(x){
    read_rds(x) 
  })
  names(res) = numbers
  
  res_pvals <- purrr::map(res, ~extract_SKATO_results(.))

  res_pvals = do.call("rbind",res_pvals)

  res_per_treat[[treat]] = res_pvals %>%
    dplyr::relocate(gene_name,gene_id)
    
  res_per_treat[[treat]]$treatment = treat
  res_per_treat[[treat]]$number = numbers

  }
  res_per_treat = do.call("rbind",res_per_treat)

  write_csv(res_per_treat,paste0(directory,mut_type,"_SKATO_scaled_centered_prop_pvals.csv"))
}

# remove genes with NAs (genes that were not tested)
res_per_treat = res_per_treat %>%
  dplyr::filter(!is.na(p_val))
# bonferroni correction across tests for genes that were actually tested (enough vars present)

res_per_treat$p_Bonf = p.adjust(res_per_treat$p_val,method = "bonferroni")

res_per_treat= res_per_treat%>%
  dplyr::arrange(p_Bonf)

table(res_per_treat$p_Bonf<0.05,res_per_treat$treatment)
#       IFN  LPS untreated
# FALSE 3681 3667      3654
# TRUE   0   0        0
table(res_per_treat$resampling_pval<0.05,res_per_treat$treatment)
# IFN  LPS untreated
# FALSE 3431 3503      3455
# TRUE   250  164       199

hist(res_per_treat$p_val) # looks uniform, even right leaning
hist(res_per_treat$p_Bonf,breaks = 100) # all 1s
hist(res_per_treat$resampling_pval) # looks uniform, even right leaning

# qqplot
plist = list()
for(treat in c("untreated","IFN","LPS")){
  
  qq_res = res_per_treat %>%
    dplyr::filter(treatment == treat & !is.na(p_val))
  
  expect.stats = qchisq(ppoints(nrow(qq_res)), df = 1, lower = F)
  obs.stats = qchisq(qq_res$p_val, df = 1, lower = F)
  lambda = median(na.omit(obs.stats)) / median(expect.stats) #GC lambda = ratio at medians
  

  plist[[treat]] = gg_qqplot(qq_res$p_val) +
    ggtitle(treat,subtitle = substitute(paste(lambda, "=", lam), list(lam = signif(lambda,3))))
  
  rm(qq_res)
}

pdf(paste0(directory,"QQ_plot_phagocytosis_SKAT_pvals.pdf"),
    width = 8, height = 4 )
plot(patchwork::wrap_plots(plist) + patchwork::plot_annotation(title = "Q-Q plots from SKAT-O results"))
dev.off()
### plot all significant genes

sign = res_per_treat %>%
  dplyr::filter(resampling_pval < 0.05) %>%
  dplyr::select(gene_name,gene_id,number, treatment,resampling_pval, p_Bonf )

  pdf(paste0(directory,mut_type,"_",treat,"_SKATO_vars_vs_phago.pdf"),
      width = 10, height = 4)
for(gene in unique(sign$gene_id)){
  num = sign %>%
    dplyr::filter(gene_id == gene) %>%
    dplyr::select(number) %>%
    distinct() %>%
    as.character(.$number)
  
  metadata = res[[num]][[1]]$metadata
  
  vars = res[[num]][[1]]$genotype %>%
    tidyr::pivot_longer(cols = -line,names_to = "ID",values_to = "genotype") %>%
    # subset to SNPs that were actually used in the test
    dplyr::filter(ID %in% names(res[[num]][[1]]$test.snp.mac)) %>%
    dplyr::left_join(metadata) 

  
  
  plist = list()
  for(var in unique(vars$ID)){
    
  plist[[var]] =  vars %>%
    dplyr::filter(ID == var) %>%
    dplyr::arrange(genotype) %>%
    ggplot(aes(x=genotype,group = genotype,y=scaled_log_fraction_pool)) +
    geom_boxplot() +
    geom_point(alpha = 0.1) +
    theme_bw() + 
    scale_x_continuous(breaks = c(0,1,2)) + 
    ggtitle(var)
  }
  p = patchwork::wrap_plots(plist,ncol = 3) + patchwork::plot_annotation(title = unique(sign[sign$gene_id == gene,"gene_name"]),
                                                                subtitle = gene)
 plot(p)
}
  
  dev.off()

### calculate aggregate effect for lines that are heterozygous or homozygous for the deletereous variants?
# SKAT developers don't provide any directionality, beta or anything similar, even for burden test

  # overlaps across treatments? directionality in burden test?
  
  
################################
# any in Sam's list?
sam_phago = read_tsv("../../../../CRISPR/OTAR2065_phagocytosis_CRISPR/data/2024_07_sam_screen_results.tsv") %>%
  dplyr::arrange(FDR)
  
  table(sam_phago$id %in% res_per_treat$gene_name) 
  
  sam_phago = sam_phago %>%
  dplyr::rename(gene_name = id) %>%
    dplyr::left_join(res_per_treat[res_per_treat$treatment=="untreated",]) %>%
    dplyr::filter(!is.na(resampling_pval)) %>%
    dplyr::arrange(resampling_pval)

# compare with burden test directionality
write_csv(sam_phago,paste0(directory,mut_type,"_untreated_sam_CRISPR_significance.csv") )
