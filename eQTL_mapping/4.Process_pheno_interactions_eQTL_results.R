# processing phenotype interactions in eqtl results
library(tidyverse)
source("./functions.R")
outDir = "../../data/results/4.Process_pheno_interactions_eQTL_results"
dir.create(outDir)

# getting gene names
message("Annotating gene names from Ensembl Ids")
# ensembl v 111, GRCh38.p14 Downloaded from https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz
# and filtered to genes
ensembl = read.csv("../../../resources/ENSEMBL_human_v111_2023/Homo_sapiens.GRCh38.111.genes.csv") %>%
  dplyr::select(gene_id,gene_name)


int = list()
pheno = list()
expr = list()
for(cond in paste(rep(c("phagocytosis","migration"), each = 3),c("untreated","IFN","LPS"),c("Not_proliferating"),sep = "_")){
  print(cond)
  int[[cond]] = read.table(paste0("../../data/results/tensorqtl/15/interactions_sum_sizefactorsNorm_log2_scaled_centered_",
                                  cond,
                                  "_common_500kb_window_tensorQT.cis_qtl_top_assoc.txt.gz"), 
                           header = TRUE)
  pheno[[cond]] = read.table(paste0("../../data/for_tensorQTL/15/interaction_df_mean_",
                                    cond,
                                    ".txt"), header=TRUE)
  expr[[cond]] = read.delim(paste0("../../data/for_tensorQTL/15/expr_sum_sizefactorsNorm_log2_scaled_centered_",
                                   cond,
                                   "_for_interaction.bed.gz")) %>%
    dplyr::rename(`#chr` = X.chr)
  
  int[[cond]] = int[[cond]] %>%
    dplyr::rename(gene_id = phenotype_id) %>%
    dplyr::left_join(ensembl)
  
  
}


# the columns b_g, b_g_se, pval_g are the effect size, 
# standard error, and p-value of g in the model, with matching columns for i and gi.
lapply(int,summary)

lapply(int,function(x) min(x$pval_adj_bh))

hist(int$migration_untreated_Not_proliferating$pval_adj_bh)
hist(int$phagocytosis_untreated_Not_proliferating$pval_adj_bh)
hist(int$migration_LPS_Not_proliferating$pval_adj_bh)

lapply(int, function(x) table(x$pval_adj_bh<0.05))

sign_res = list()
for(nam in names(int)){
  sign_res[[nam]] = int[[nam]] %>%
    dplyr::filter(pval_adj_bh < 0.05) %>%
    dplyr::arrange(pval_adj_bh) %>%
    dplyr::mutate(condition = nam)
  
}

sign_res = do.call("rbind",sign_res)
length(unique(sign_res$gene_name))

## plots top gi
# y - expr
# x - phago

# slopes - genotypes
# phagocytosis_unt = int$phagocytosis_untreated_Not_proliferating %>%
#   dplyr::filter(pval_adj_bh < 0.05) %>%
#   dplyr::arrange(pval_adj_bh)
# 
# migration_unt = int$migration_untreated_Not_proliferating %>%
#   dplyr::filter(pval_adj_bh < 0.05) %>%
#   dplyr::arrange(pval_adj_bh)
# 
# # phagocytosis test
# expr_phago_unt = expr$phagocytosis_untreated_Not_proliferating
# pheno_phago_unt = pheno$phagocytosis_untreated_Not_proliferating

# fix this to make the table come from the interaction results themselves (this is from results without interactions)
ref_alt_donor = read_rds(paste0("../../data/results/4.Inspect_eQTL_results","/vcf_tibble_sign_tensorQTL_REF_ALT.rds")) %>%
  dplyr::filter(names %in% sign_res$variant_id)

int_subset = sign_res %>%
  dplyr::filter(variant_id %in% ref_alt_donor$names)


# 
# var_to_plot = unlist(int_subset[int_subset$pval_adj_bh == min(int_subset$pval_adj_bh),"variant_id"])
# 
# # subset expression data
# expr_long = expr_phago_unt %>%
#   dplyr::filter(gene_id == unlist(int_subset[int_subset$variant_id == var_to_plot,"gene_id"])) %>%
#   dplyr::select(-`#chr`,-start,-end) %>%
#   tidyr::pivot_longer(cols = !gene_id,
#                       names_to = "line",values_to = "scaled_expression") %>%
#   dplyr::left_join(ensembl) %>%
#   dplyr::select(gene_name, line, scaled_expression)
# 
# 
# test = ref_alt_donor %>%
#   dplyr::rename(variant_id = names) %>%
#   dplyr::filter(variant_id == var_to_plot) %>%
#   dplyr::select(-position) %>%
#   dplyr::left_join(int_subset,by = "variant_id") %>%
#   tidyr::pivot_longer(cols = !c(colnames(int_subset)),
#                                 names_to = "line",values_to = "genotype") %>%
#   dplyr::left_join(pheno_phago_unt,by= "line") %>%
#   dplyr::left_join(expr_long) 
# 
# # counts phased (eg C|T and T|C) as two different genotypes - fixing
# 
# 
# test = test %>%
#   dplyr::mutate(genotype = case_when(genotype == paste(str_split_1(var_to_plot,pattern = "_")[4],
#                                                        str_split_1(var_to_plot,pattern = "_")[3],
#                                                        sep = "|") ~ paste(str_split_1(var_to_plot,pattern = "_")[3],
#                                                   str_split_1(var_to_plot,pattern = "_")[4],sep = "|"), # if genotype is opposite of REF -> ALT het, swap
#                                      .default = genotype)) # if any other genotype combination, leave as is
# 
# ### plotting different regression lines per group, for visualization of slopes
# lm_results = test %>%
#   group_by(genotype) %>%
#   do(model =  broom::tidy(lm(scaled_expression ~ log_fraction_mean, data = .))) %>%
#   tidyr::unnest(cols = c(model)) %>%
#   select(genotype, term, estimate)
# 
# pdf(paste0(outDir,"/",unique(test$variant_id), unique(test$gene_name),"_phagocytosis_untreated_Not_proliferating_60PCs.pdf"), 
#     width = 7, height = 7)
#   p = ggplot(test, aes(x = log_fraction_mean,y=scaled_expression, col = genotype)) +
#   geom_point() + 
#   theme_minimal() +
#     geom_smooth(method = "lm", se = FALSE, aes(group = genotype, color =genotype)) +
#     
#   xlab("Average of log(mCherry+ prop. / mCherry- prop.)") +
#   ylab("Scaled gene expression") + 
#     # geom_text(data = lm_results, aes(label = paste("Estimate (", term, ") =", round(estimate, 2))),
#     #           x = Inf, y = -Inf, hjust = 1, vjust = 0) +
#   ggtitle(subtitle = paste(unique(test$variant_id), unique(test$gene_name), sep = " - "),
#           label = "Interaction of eQTL with phagocytosis phenotype")
#   plot(p)
# 
# dev.off()

# migration test
migration_unt = int$migration_untreated_Not_proliferating %>%
  dplyr::filter(pval_adj_bh < 0.05) %>%
  dplyr::arrange(pval_adj_bh)

# phagocytosis test
expr_migr_unt = expr$migration_untreated_Not_proliferating
pheno_migr_unt = pheno$migration_untreated_Not_proliferating

# fix this to make the table come from the interaction results themselves (this is from results without interactions)
ref_alt_donor = read_rds(paste0("../../data/results/4.Inspect_eQTL_results","/vcf_tibble_sign_tensorQTL_REF_ALT.rds")) %>%
  dplyr::filter(names %in% migration_unt$variant_id)

int_subset = migration_unt %>%
  dplyr::filter(variant_id %in% ref_alt_donor$names)

var_to_plot = unlist(int_subset[int_subset$pval_adj_bh == min(int_subset$pval_adj_bh),"variant_id"])

# the top hit seems to be driven by outliers

var_to_plot = unlist(int_subset[int_subset$gene_name == "CDIPT","variant_id"])
# subset expression data
expr_long = expr_migr_unt %>%
  dplyr::filter(gene_id == unlist(int_subset[int_subset$variant_id == var_to_plot,"gene_id"])) %>%
  dplyr::select(-`#chr`,-start,-end) %>%
  tidyr::pivot_longer(cols = !gene_id,
                      names_to = "line",values_to = "scaled_expression") %>%
  dplyr::left_join(ensembl) %>%
  dplyr::select(gene_name, line, scaled_expression)


test = ref_alt_donor %>%
  dplyr::rename(variant_id = names) %>%
  dplyr::filter(variant_id == var_to_plot) %>%
  dplyr::select(-position) %>%
  dplyr::left_join(int_subset,by = "variant_id") %>%
  tidyr::pivot_longer(cols = !c(colnames(int_subset)),
                      names_to = "line",values_to = "genotype") %>%
  dplyr::left_join(pheno_migr_unt,by= "line") %>%
  dplyr::left_join(expr_long) 

# counts phased (eg C|T and T|C) as two different genotypes - fixing


test = test %>%
  dplyr::mutate(genotype = case_when(genotype == paste(str_split_1(var_to_plot,pattern = "_")[4],
                                                       str_split_1(var_to_plot,pattern = "_")[3],
                                                       sep = "|") ~ paste(str_split_1(var_to_plot,pattern = "_")[3],
                                                                          str_split_1(var_to_plot,pattern = "_")[4],sep = "|"), # if genotype is opposite of REF -> ALT het, swap
                                     .default = genotype)) # if any other genotype combination, leave as is

### plotting different regression lines per group, for visualization of slopes
lm_results = test %>%
  group_by(genotype) %>%
  do(model =  broom::tidy(lm(scaled_expression ~ log_fraction_mean, data = .))) %>%
  tidyr::unnest(cols = c(model)) %>%
  dplyr::select(genotype, term, estimate)

pdf(paste0(outDir,"/",unique(test$variant_id), "_",unique(test$gene_name),"_migration_untreated_Not_proliferating_60PCs.pdf"), 
    width = 7, height = 7)
p = ggplot(test, aes(x = log_fraction_mean,y=scaled_expression, col = genotype)) +
  geom_point() + 
  theme_minimal() +
  geom_smooth(method = "lm", se = FALSE, aes(group = genotype, color =genotype)) +
  
  xlab("Average of log(bottom prop. / top prop.)") +
  ylab("Scaled gene expression") + 
  # geom_text(data = lm_results, aes(label = paste("Estimate (", term, ") =", round(estimate, 2))),
  #           x = Inf, y = -Inf, hjust = 1, vjust = 0) +
  ggtitle(subtitle = paste(unique(test$variant_id), unique(test$gene_name), sep = " - "),
          label = "Interaction of eQTL with migration phenotype")
plot(p)

dev.off()

