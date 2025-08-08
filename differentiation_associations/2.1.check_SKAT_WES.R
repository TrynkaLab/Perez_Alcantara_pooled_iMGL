# examine proliferation SKATO tests
library(patchwork)
library(tidyverse)
source("./functions.R")
directory = "../../data/results/2.1.check_SKAT_WES/"
dir.create(directory, recursive = TRUE)

options(future.globals.maxSize = 240000 * 1024^2) # 240Gb

comparisons = c("line_prop_changes_old_vs_young_premac","line_prop_changes_microglia_premac","line_prop_changes_premac_iPSC")
for(mut_type in c("missense_non_deleterious","deleterious")){
  res = list()
  
for(comparison in comparisons){
  
    files = list.files(paste0("../../data/results/2.rare_var_SKAT_WES/",comparison,"/"),
                       pattern = paste0("[0-9]","_",mut_type,"_SKATO_result.RDS"))
    files = str_sort(files, numeric = TRUE) # to later recover right plots by number
    numbers = str_split_i(files,pattern = "_",i = 1) #they may not be contiguous
    res[[comparison]]  = lapply(paste0("../../data/results/2.rare_var_SKAT_WES/",comparison,"/",files), function(x){
      read_rds(x)
    })
    names(res[[comparison]] ) = numbers

    
    res_pvals = purrr::map(res[[comparison]] , ~extract_SKATO_results(.))
    
    res_pvals = do.call("rbind",res_pvals)
    
    res[[comparison]]  = res_pvals %>%
      dplyr::relocate(gene_name,gene_id)
    
  
    res[[comparison]] $number = numbers
    res[[comparison]] $comparison = comparison
    
  

}
  res = do.call("rbind",res)
  
  write_csv(res,paste0(directory,mut_type,"_SKATO_scaled_centered_prop_pvals.csv"))
}
# remove genes with NAs (genes that were not tested)
res = res %>%
  dplyr::filter(!is.na(p_val))
# bonferroni correction across tests for genes that were actually tested (enough vars present)

res$p_Bonf = p.adjust(res$p_val,method = "bonferroni")

res= res%>%
  dplyr::arrange(p_Bonf)

table(res$p_Bonf<0.05,res$comparison)
#       line_prop_changes_microglia_premac line_prop_changes_old_vs_young_premac line_prop_changes_premac_iPSC
# FALSE 3840 3746      3824

table(res$resampling_pval<0.05,res$comparison)
# line_prop_changes_microglia_premac line_prop_changes_old_vs_young_premac line_prop_changes_premac_iPSC
# FALSE                               3575                                  3511                          3641
# TRUE                                 265                                   235                           183

hist(res$p_val) # looks uniform, even right leaning
hist(res$p_Bonf) # right tailed as expected
hist(res$resampling_pval) # looks uniform, even right leaning

# qqplot
plist = list()
for(comp in comparisons){
  
  qq_res = res %>%
    dplyr::filter(comparison == comp & !is.na(p_val))
  
  expect.stats = qchisq(ppoints(nrow(qq_res)), df = 1, lower = F)
  obs.stats = qchisq(qq_res$p_val, df = 1, lower = F)
  lambda = median(na.omit(obs.stats)) / median(expect.stats) #GC lambda = ratio at medians
  
  
  plist[[comp]] = gg_qqplot(qq_res$p_val) +
    ggtitle(comp,subtitle = substitute(paste(lambda, "=", lam), list(lam = signif(lambda,3))))
  
  rm(qq_res)
}




#### below we're missing metadata genotype information

pdf(paste0(directory,"QQ_plot_proliferation_GWAS_pvals.pdf"),
    width = 8, height = 4 )
plot(patchwork::wrap_plots(plist) + patchwork::plot_annotation(title = "Q-Q plots from SKAT-O results"))
dev.off()
### plot all significant genes

sign = res %>%
  dplyr::filter(resampling_pval < 0.05) %>%
  dplyr::select(gene_name,gene_id,number, comparison,resampling_pval, p_Bonf )

pdf(paste0(directory,mut_type,"_",treat,"_SKATO_vars_vs_proliferation.pdf"),
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

# overlaps across comparisons? directionality in burden test?
# compare with prolif sign genes - GSEA to check the sign ones here are not overrepresented at the top of phago/migration genes

################################
# any in Sam's list?
sam_phago = read_tsv("../../../../CRISPR/OTAR2065_phagocytosis_CRISPR/data/2024_07_sam_screen_results.tsv") %>%
  dplyr::arrange(FDR)

table(sam_phago$id %in% res$gene_name) 

sam_phago = sam_phago %>%
  dplyr::rename(gene_name = id) %>%
  dplyr::left_join(res[res$comparison=="untreated",]) %>%
  dplyr::filter(!is.na(resampling_pval)) %>%
  dplyr::arrange(resampling_pval)

# compare with burden test directionality
write_csv(sam_phago,paste0(directory,mut_type,"_untreated_sam_CRISPR_significance.csv") )
