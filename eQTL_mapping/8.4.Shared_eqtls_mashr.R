# checking shared/specific eQTLs with mashr
# after loading module containing R libraries with softpack
# module load HGI/softpack/groups/otar2065/otar2065_test2/1
# based on https://cran.r-project.org/web/packages/mashr/vignettes/intro_mash.html
# https://cran.r-project.org/web/packages/mashr/vignettes/eQTL_outline.html
# and https://cran.r-project.org/web/packages/mashr/vignettes/mash_sampling.html
library(tidyverse)
library(mashr)
library(UpSetR)
source("./functions.R")

set.seed(123)
outdir = "../../data/results/8.4.Shared_eqtls_mashr"
treatment_cols =  c(untreated = "#8D918B", IFN = "#3A5683", LPS = "#F8766D")

dir.create(outdir,recursive = TRUE)
#mashr needs a strong subset of eQTLs (i.e. significant). Reasonable size would be 16k tests for 44 conditions (the top eQTL per gene for all conditions)
# and a random subset ~ 20k randomly-selected tests for the 16k strong results from above.
eqtl = readr::read_csv("../../data/results/4.Inspect_eQTL_results/tensorQTL_variant_gene.csv")
# extract significant results
signif = eqtl %>%
  dplyr::filter(cluster %in% "Not_proliferating") %>%
  dplyr::filter(qval<0.05)

# extract the top eQTL per gene for all conditions
strong = eqtl %>%
  dplyr::filter(cluster %in% "Not_proliferating") %>%
  dplyr::filter(gene_variant_pair %in% unique(signif$gene_variant_pair)) %>%
  dplyr::mutate(ensembl_id_var = paste(gene_id,variant_id,sep = "-"))
# extract matrices of observed estimates (slope x condition) and their standard errors (s.e.of slope x condition)
eqtl_nominal = list()
# Not loading Prolif subsets
for(condition in c(paste0(c(63,73,84),"/sum_sizefactorsNorm_log2_scaled_centered_",c("untreated","IFN","LPS"),"_Not_proliferating"))){
  eqtl_nominal[[condition]] = readr::read_delim(paste0("../../data/results/tensorqtl/",
                                                       condition,"_common_500kb_window_tensorQTL_nominal.txt")) %>%
    dplyr::mutate(group = paste(str_split_i(condition, pattern = "/",1),str_split_i(condition, pattern = "_",6),"Not_proliferating",sep="_"),
                  ensembl_id_var = paste(phenotype_id,variant_id,sep = "-")) %>%
    dplyr::filter(ensembl_id_var %in% strong$ensembl_id_var) %>% # extracting same gene-variant pair for all
    
    dplyr::select(slope,slope_se,group,ensembl_id_var)
  
  
}

eqtl_nominal = do.call("rbind",eqtl_nominal)

strong = strong %>%
  dplyr::select(slope,slope_se,group,ensembl_id_var) %>%
  dplyr::bind_rows(eqtl_nominal) %>%
  # decimal rounding in nominal and cis_map is slightly different after the 5th decimal!
  # selecting the first value that appears after binding (the cis_map one)
  dplyr::mutate(group_var = paste(group,ensembl_id_var,sep = "-")) %>%
  dplyr::distinct(group_var, .keep_all = TRUE)
# pivot to wider
strong_Bhat = strong %>%
  tidyr::pivot_wider(names_from = group, values_from = slope,id_cols = ensembl_id_var) %>%
  tibble::column_to_rownames(var = "ensembl_id_var") %>%
  as.matrix()
nrow(strong_Bhat) # 11238

strong_Shat = strong %>%
  tidyr::pivot_wider(names_from = group, values_from = slope_se,id_cols = ensembl_id_var) %>%
  tibble::column_to_rownames(var = "ensembl_id_var") %>%
  as.matrix()

nrow(strong_Shat) == nrow(strong_Bhat)
# extract random sample from the nominal
# loading again nominal without subseting to strong variant pair positions

# select genes that are present in the strong signals to then extract 5 random snps from them
genes = unique(stringr::str_split_fixed(strong$ensembl_id_var,pattern = "-",n=2)[,1])

eqtl_nominal = list()
for(condition in c(paste0(c(63,73,84),"/sum_sizefactorsNorm_log2_scaled_centered_",c("untreated","IFN","LPS"),"_Not_proliferating"))){
  eqtl_nominal[[condition]] = readr::read_delim(paste0("../../data/results/tensorqtl/",
                                                       condition,"_common_500kb_window_tensorQTL_nominal.txt")) %>%
    dplyr::filter(phenotype_id %in% genes) %>%
    dplyr::mutate(group = paste(str_split_i(condition, pattern = "/",1),str_split_i(condition, pattern = "_",6),"Not_proliferating",sep="_"),
                  ensembl_id_var = paste(phenotype_id,variant_id,sep = "-")) %>%
    dplyr::select(slope,slope_se,group,ensembl_id_var,phenotype_id,variant_id)
  
  
}

eqtl_nominal = do.call("rbind",eqtl_nominal)

random_subset_gene_variant = eqtl_nominal %>%
  dplyr::group_by(phenotype_id) %>%
  dplyr::slice_sample(n= 5) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(ensembl_id_var) %>%
  .$ensembl_id_var

length(unique(random_subset_gene_variant)) # 35308

# and transform to mashr format

random_Bhat = eqtl_nominal %>%
  dplyr::filter(ensembl_id_var %in% random_subset_gene_variant) %>%
  dplyr::select(slope,group,ensembl_id_var) %>%
  tidyr::pivot_wider(names_from = group, values_from = slope, id_cols = ensembl_id_var) %>%
  tibble::column_to_rownames(var = "ensembl_id_var") %>%
  as.matrix()

random_Shat = eqtl_nominal %>%
  dplyr::filter(ensembl_id_var %in% random_subset_gene_variant) %>%
  dplyr::select(slope_se,group,ensembl_id_var) %>%
  tidyr::pivot_wider(names_from = group, values_from = slope_se, id_cols = ensembl_id_var) %>%
  tibble::column_to_rownames(var = "ensembl_id_var") %>%
  as.matrix()

nrow(random_Bhat)
nrow(random_Shat) == nrow(random_Bhat) 

# 1. Learn correlation structure among null tests using random test.
# very important that random test is not enriched for strong effects
summary(random_Bhat)
summary(strong_Bhat) # the quartiles are much more extreme, the median is strongly shifted

data.temp = mashr::mash_set_data(Bhat = random_Bhat,Shat = random_Shat)
Vhat = mashr::estimate_null_correlation_simple(data.temp)
rm(data.temp)


# set up main objects
data.random = mashr::mash_set_data(random_Bhat,random_Shat,V=Vhat)
data.strong = mashr::mash_set_data(strong_Bhat,strong_Shat, V=Vhat)


# 2. Learn data-driven covariance matrices using strong tests.
U.pca = mashr::cov_pca(data.strong,npc = length(unique(strong$group)))
U.ed = mashr::cov_ed(data.strong, U.pca)

#3. Fit the mashr model to the random tests, to learn the mixture weights on all 
#the different covariance matrices and scaling coefficients.
U.c = mashr::cov_canonical(data.random)
m = mashr::mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1)

# 4. Compute posterior summaries on the strong tests, using the model fit from step 2. 
m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE)

readr::write_rds(m2,paste0(outdir,"/mashr_results.rds"))

# m2 = readr::read_rds(paste0(outdir,"/mashr_results.rds"))

# barplot(mashr::get_estimated_pi(m),las=2)
# # get posterior summaries
# head(ashr::get_pm(m2))
# mash_plot_meta(m2,get_significant_results(m2)[1])

# what proportion of significant effects have the same sign and similar magnitude for each pair of conditions?
print(mashr::get_pairwise_sharing(m2, factor=0.5, lfsr_thresh = 0.05)) # same sign and within a factor of 0.5 from each other (default)
# For each pair of tissues, first identify the effects that are significant (by lfsr<lfsr_thresh) IN AT LEAST ONE of the two tissues.
# check if what is considered significant here can be insignificant in my eQTL data

# Then compute what fraction of these have an estimated (posterior mean) effect size within a factor 'factor' of one another. 
corrplot::corrplot(mashr::get_pairwise_sharing(m2, factor=0.5, lfsr_thresh = 0.05),
                   method='color', col.lim=c(0,1), type='upper', addCoef.col = "black", tl.col="black", tl.srt=45, title = 'Pairwise Sharing by Magnitude', mar = c(4,0,4,0))
# this is very high sharing compared to blagoje's paper
# it may be because in his paper he's keeping the top variant per gene, and not doing
# the additional qvalue filter. This may give more room for comparison effects with different magnitudes and directions per cell type/condition

ann_color = list("prolif_x_treatment" = c("Not_proliferating_LPS" = "yellow","Proliferating_LPS" = "gold", 
                                          "Not_proliferating_IFN" = "coral","Proliferating_IFN" = "red",
                                          "Not_proliferating_untreated" = "grey","Proliferating_untreated" = "black"),
                 "treatment"=c(untreated = "#8D918B", IFN = "#3A5683", LPS = "#F8766D"),
                 "proliferation_status"=c("Not_proliferating"="lightblue","Proliferating"="darkblue"))
metadata = data.frame("sample" = colnames(mashr::get_pairwise_sharing(m2, factor=0.5, lfsr_thresh = 0.05)),
                      "proliferation_status" = rep("Not_proliferating",3))
rownames(metadata) = metadata$sample
metadata$treatment = str_split_i(metadata$sample,pattern="_",i=2)

png(paste0(outdir,"/mashr_heatmap_treatments.png"),
    width = 4.5, height = 3.2, res = 400,units = "in", type = "cairo")
mashr::get_pairwise_sharing(m2, factor=0.5, lfsr_thresh = 0.05) %>%
  magrittr::set_rownames(str_split_i(rownames(.),pattern="_",i=2)) %>% # careful here
  pheatmap::pheatmap(., annotation = metadata[, c("treatment"), drop=F],
                     show_colnames = FALSE,
                     annotation_colors = ann_color,
                     display_numbers = TRUE,
                     fontsize_number = 14,
                     color = colorRampPalette(c( "#dddddd","#aaaaaa"))( 10 ))

dev.off()

# just same sign
print(mashr::get_pairwise_sharing(m2, factor=0, lfsr_thresh = 0.05)) 

fsr = ashr::get_lfsr(m2) # get local false sign rate (significance)
gene_variants = rownames(fsr)
colSums(fsr<0.05) # variants with significant effects per condition
pm = ashr::get_pm(m2) # get posterior mean

table(get_n_significant_conditions(m2))

full_share_significant = get_n_significant_conditions(m2)[get_n_significant_conditions(m2)==3] # # get genes that are significant across all three conditions
two_share_significant = get_n_significant_conditions(m2)[get_n_significant_conditions(m2)==2] # two conditions
unique_significant = get_n_significant_conditions(m2)[get_n_significant_conditions(m2)==1] # only one

fsr_two_share_significant = fsr[names(two_share_significant),] %>%
  as.data.frame() %>%
  rownames_to_column("rows") %>%
  dplyr::rename_with(.,~ paste0("fsr_", .x, recycle0 = TRUE),ends_with("Not_proliferating"))
# reshape effect table 
pm_rownames = rownames(pm)
effect_reshaped = pm %>%
  as_tibble() %>%
  dplyr::mutate(rows = pm_rownames) %>%
  dplyr::mutate(ensembl = str_split_i(rows,pattern = "-",i = 1),
                variant = str_split_i(rows,pattern = "-",i = 2))

ensembl =read.csv("/lustre/scratch123/hgi/teams/trynka/resources/biomart/Homo_sapiens.GRCh38.111.genes.csv") %>%
  dplyr::rename("ensembl"="gene_id") %>%
  dplyr::select("ensembl","gene_name")

effect_reshaped = effect_reshaped %>%
  dplyr::left_join(ensembl)

# get shared genes within effect <0.5 
# get shared genes regardless of effect (same sign)
# get shared genes different sign

shared_genes_all = effect_reshaped %>%
  dplyr::filter(rows %in% names(full_share_significant)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    same_sign = all(sign(c_across(`63_untreated_Not_proliferating` :`84_LPS_Not_proliferating`)) == sign(`63_untreated_Not_proliferating`)), # Check if all values have the same sign
    within_05 = max(c_across(`63_untreated_Not_proliferating` :`84_LPS_Not_proliferating`)) - min(c_across(`63_untreated_Not_proliferating` :`84_LPS_Not_proliferating`)) <= 0.5 # Check if the values are within 0.5 of each other
  ) %>%
  dplyr::ungroup() %>%
  mutate(check_shared_genes_05_effect = same_sign & within_05,
         check_shared_genes_distinct_effect =  same_sign & !(within_05),
         check_shared_genes_different_sign = !same_sign)
# numbers
shared_genes_all %>% dplyr::filter(check_shared_genes_05_effect == TRUE) %>% distinct(gene_name) %>% nrow() # 5057
shared_genes_all %>% dplyr::filter(check_shared_genes_distinct_effect == TRUE) %>% distinct(gene_name) %>% nrow() # 44
shared_genes_all %>% dplyr::filter(check_shared_genes_different_sign == TRUE) %>% distinct(gene_name) %>% nrow() # 26


# genes for upsetr
toplot = shared_genes_all %>% dplyr::filter(check_shared_genes_05_effect == TRUE) %>% distinct(gene_name) %>%
  dplyr::mutate(category = rep("Significant across all with same effect",
                               shared_genes_all %>% dplyr::filter(check_shared_genes_05_effect == TRUE) %>% distinct(gene_name) %>% nrow() ))
toplot = rbind(toplot,shared_genes_all %>% dplyr::filter(check_shared_genes_distinct_effect == TRUE) %>% distinct(gene_name) %>%
                 dplyr::mutate(category = rep("Significant across all with varying effects",shared_genes_all %>% dplyr::filter(check_shared_genes_distinct_effect == TRUE) %>% distinct(gene_name) %>% nrow() )))

toplot = rbind(toplot,shared_genes_all %>% dplyr::filter(check_shared_genes_different_sign == TRUE) %>% distinct(gene_name)  %>%
                 dplyr::mutate(category = rep("Significant across all with opposite signs", shared_genes_all %>% dplyr::filter(check_shared_genes_different_sign == TRUE) %>% distinct(gene_name) %>% nrow())))

# get unique significant genes - remove also those that present a different variant that is significant and shared across all 
p_vals_unique = fsr[rownames(fsr) %in% names(unique_significant),]
p_vals_unique_untreated = names((p_vals_unique[,"63_untreated_Not_proliferating"]<0.05)[p_vals_unique[,"63_untreated_Not_proliferating"]<0.05])
p_vals_unique_IFN= names((p_vals_unique[,"73_IFN_Not_proliferating"]<0.05)[p_vals_unique[,"73_IFN_Not_proliferating"]<0.05])
p_vals_unique_LPS= names((p_vals_unique[,"84_LPS_Not_proliferating"]<0.05)[p_vals_unique[,"84_LPS_Not_proliferating"]<0.05])

unique_genes = effect_reshaped %>%
  dplyr::filter(rows %in% names(unique_significant)) %>%
  dplyr::mutate(
    unique = case_when( rows %in% p_vals_unique_untreated ~ "untreated",
                        rows %in% p_vals_unique_IFN ~ "IFN",
                        rows %in% p_vals_unique_LPS ~ "LPS"),
    also_shared_gene = gene_name %in% shared_genes_all$gene_name
  ) 


unique_genes = unique_genes %>%
  dplyr::filter(also_shared_gene == FALSE) # removing those shared eGenes


toplot = rbind(toplot,unique_genes  %>% dplyr::filter(unique == "untreated") %>% distinct(gene_name)  %>%
                 dplyr::mutate(category = rep("Significant only in untreated", unique_genes  %>% dplyr::filter(unique == "untreated") %>% distinct(gene_name)  %>% nrow())))
toplot = rbind(toplot,unique_genes  %>% dplyr::filter(unique == "IFN") %>% distinct(gene_name)  %>%
                 dplyr::mutate(category = rep("Significant only in IFN", unique_genes  %>% dplyr::filter(unique == "IFN") %>% distinct(gene_name)  %>% nrow())))
toplot = rbind(toplot,unique_genes  %>% dplyr::filter(unique == "LPS") %>% distinct(gene_name)  %>%
                 dplyr::mutate(category = rep("Significant only in LPS", unique_genes  %>% dplyr::filter(unique == "LPS") %>% distinct(gene_name)  %>% nrow())))

# get pairwise shared IFN-LPS and IFN-untreated and LPS-untreated ONLY

shared_pairwise = effect_reshaped %>%
  left_join(fsr_two_share_significant ) %>%
  dplyr::filter(rows %in% names(two_share_significant)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    same_sign_IFN_untreated = all(sign(c_across(`63_untreated_Not_proliferating` :`73_IFN_Not_proliferating`)) == sign(`63_untreated_Not_proliferating`)), # Check if all values have the same sign
    within_05_IFN_untreated = max(c_across(`63_untreated_Not_proliferating` :`73_IFN_Not_proliferating`)) - min(c_across(`63_untreated_Not_proliferating` :`73_IFN_Not_proliferating`)) <= 0.5, # Check if the values are within 0.5 of each other
    sign_IFN_untreated = fsr_63_untreated_Not_proliferating < 0.05 & fsr_73_IFN_Not_proliferating < 0.05,
    same_sign_LPS_untreated = all(sign(c_across(c(`63_untreated_Not_proliferating` ,`84_LPS_Not_proliferating`))) == sign(`63_untreated_Not_proliferating`)), 
    within_05_LPS_untreated = max(c_across(c(`63_untreated_Not_proliferating` ,`84_LPS_Not_proliferating`))) - min(c_across(c(`63_untreated_Not_proliferating` ,`84_LPS_Not_proliferating`))) <= 0.5, 
    sign_LPS_untreated = fsr_63_untreated_Not_proliferating < 0.05 & fsr_84_LPS_Not_proliferating < 0.05,
    same_sign_IFN_LPS = all(sign(c_across(c(`73_IFN_Not_proliferating` ,`84_LPS_Not_proliferating`))) == sign(`73_IFN_Not_proliferating`)), 
    within_05_IFN_LPS = max(c_across(c(`73_IFN_Not_proliferating` ,`84_LPS_Not_proliferating`))) - min(c_across(c(`73_IFN_Not_proliferating` ,`84_LPS_Not_proliferating`))) <= 0.5 ,
    sign_IFN_LPS = fsr_73_IFN_Not_proliferating < 0.05 & fsr_84_LPS_Not_proliferating < 0.05
    
  ) %>%
  dplyr::ungroup() %>%
  mutate(check_shared_genes_05_effect_IFN_untreated = same_sign_IFN_untreated & within_05_IFN_untreated & sign_IFN_untreated,
         check_shared_genes_05_effect_LPS_untreated = same_sign_LPS_untreated & within_05_LPS_untreated & sign_LPS_untreated,
         check_shared_genes_05_effect_IFN_LPS = same_sign_IFN_LPS & within_05_IFN_LPS & sign_IFN_LPS
  )




toplot = rbind(toplot,shared_pairwise  %>% dplyr::filter(check_shared_genes_05_effect_IFN_untreated == TRUE) %>% distinct(gene_name)  %>%
                 dplyr::mutate(category = rep("Significant shared effects IFN-untreated",shared_pairwise  %>% dplyr::filter(check_shared_genes_05_effect_IFN_untreated == TRUE) %>% distinct(gene_name)  %>% nrow())))

toplot = rbind(toplot,shared_pairwise  %>% dplyr::filter(check_shared_genes_05_effect_LPS_untreated == TRUE) %>% distinct(gene_name)  %>%
                 dplyr::mutate(category = rep("Significant shared effects LPS-untreated",shared_pairwise  %>% dplyr::filter(check_shared_genes_05_effect_LPS_untreated == TRUE) %>% distinct(gene_name)  %>% nrow())))


toplot = rbind(toplot,shared_pairwise  %>% dplyr::filter(check_shared_genes_05_effect_IFN_LPS == TRUE) %>% distinct(gene_name)  %>%
                 dplyr::mutate(category = rep("Significant shared effects IFN-LPS",shared_pairwise  %>% dplyr::filter(check_shared_genes_05_effect_IFN_LPS == TRUE) %>% distinct(gene_name)  %>% nrow())))

write_csv(shared_pairwise,paste0(outdir,"/mashr_shared_genes_pairwise_effects_2_conditions_within_05.csv"))
write_csv(unique_genes,paste0(outdir,"/mashr_unique_genes.csv"))
write_csv(shared_genes_all,paste0(outdir,"/mashr_shared_genes_3_conditions.csv"))
write_csv(toplot,paste0(outdir,"/mashr_summary_shared_effects_gene_lists.csv"))

#### barplots with results
pct_format = scales::percent_format(accuracy = .1)


p = toplot %>%
  dplyr::mutate(sharing = case_when(category %in% c("Significant across all with same effect" ,   
                                                    "Significant across all with varying effects",
                                                    "Significant across all with opposite signs") ~ "Shared across all",
                                    category %in% c("Significant only in untreated"   ,           
                                                    "Significant only in IFN" , "Significant only in LPS"    ) ~ "Unique",
                                    category %in% c("Significant shared effects IFN-untreated",
                                                    "Significant shared effects LPS-untreated"  ,  
                                                    "Significant shared effects IFN-LPS"    ) ~ "Shared across two")) %>%
  dplyr::mutate(category = factor(category,levels = c("Significant across all with same effect",
                                                      "Significant across all with varying effects",
                                                      "Significant across all with opposite signs",
                                                      "Significant shared effects IFN-untreated",
                                                      "Significant shared effects LPS-untreated",
                                                      "Significant shared effects IFN-LPS" ,
                                                      "Significant only in untreated" ,
                                                      "Significant only in LPS" ,
                                                      "Significant only in IFN" ), ordered = TRUE)) %>%
  ggplot(aes(x =   category)) + 
  geom_bar(color = "grey40",fill = "grey40") +
  theme_minimal()  +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_y_continuous(transform='pseudo_log',breaks = c(0,100,200,500,2000),limits = c(0,100000)) + 
  xlab("") + ylab("") +
  geom_text(
    aes(
      label = sprintf(
        '%d (%s)',
        after_stat(count),
        pct_format(after_stat(count) / sum(after_stat(count)))
      )
    ),
    stat = 'count',
    nudge_y = 0.1,hjust = 0,
    colour = 'grey30',
    size = 4
  )

p

pdf(paste0(outdir,"/mashr_shared_effects_barplot.pdf"),width = 7,height = 3.4)
plot(p)

dev.off()
