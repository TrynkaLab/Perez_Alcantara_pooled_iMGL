# gather proportion QTL results and visualize

.libPaths(c("/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/",
            "/software/teamtrynka/conda/otar2065/lib/R/library",
            "/software/teamtrynka/ma23/R4.1/libs"))
library(patchwork)
library(tidyverse)
library(qvalue)
source("./functions.R")

input_path = "../../data/results/5.proportion_QTL/"
output_path="../../data/results/8.check_results/"
dir.create(file.path(output_path),recursive = TRUE)

pqtl_res_files = list.files(path =input_path, pattern = "prop_QTL_\\d+\\.csv")

# read in csv files from lm() in order
# Extract the numerical order from each file name and sort
order_list = as.integer(sub("prop_QTL_(\\d+)\\.csv", "\\1", pqtl_res_files))
sorted_file_names = pqtl_res_files[order(order_list)]

# Read and append the files
pqtl_res = NULL
for (file in sorted_file_names) {
   data = read.csv(paste0(input_path,file))
   pqtl_res = rbind(pqtl_res, data)
}

# omit NAs - regions without eqtl overlaps
pqtl_res = pqtl_res[!is.na(pqtl_res$snp),]
anyDuplicated(pqtl_res$snp)

# double check variants
genotype = readr::read_csv("../../data/genotypes/full_genotype/genotype_minor_allele_dosage.csv") %>%
      dplyr::relocate(rn)

table(genotype$rn %in% paste0("chr",pqtl_res$snp)) # it's normal that many snps from genotype are not in results
table(paste0("chr",pqtl_res$snp) %in% genotype$rn ) # all snps from results should be in genotype (500k)

# multiple correction adjustment ####
# q-values, local FDR:
pqtl_res$q_untreated= qvalue(pqtl_res$p_untreated)$qvalues
pqtl_res$q_IFN = qvalue(pqtl_res$p_IFN)$qvalues
pqtl_res$q_LPS = qvalue(pqtl_res$p_LPS)$qvalues

# Bonferroni, the usual correction in GWAS
pqtl_res$p_bonf_untreated = p.adjust(pqtl_res$p_untreated, method="bonferroni")
pqtl_res$p_bonf_IFN = p.adjust(pqtl_res$p_IFN, method="bonferroni")
pqtl_res$p_bonf_LPS = p.adjust(pqtl_res$p_LPS, method="bonferroni")

pqtl_res$p_BH_untreated = p.adjust(pqtl_res$p_untreated, method="BH")
pqtl_res$p_BH_IFN = p.adjust(pqtl_res$p_IFN, method="BH")
pqtl_res$p_BH_LPS = p.adjust(pqtl_res$p_LPS, method="BH")

nrow(pqtl_res) # 529,669 
length(unique(pqtl_res$snp)) # 529,669

# 
# # raw p-values distribution - untreated
# p1 = ggplot(pqtl_res,aes(x = p_untreated)) + 
#    geom_histogram() +
#    theme_bw() 
# p1
# 
# summary(pqtl_res$p_untreated) #median of P-values should be 0.5 under the null
# # median = 0.32 - this is not the null so we'd expect lower p-values
# ggplot(pqtl_res[1:10000,],aes(x = 1:10000, y = -log10(p_untreated))) + 
#    geom_point() + 
#    theme_bw() + 
#    ylab("-log10 P-value") + 
#    xlab("variant")

## q-q plot using chi squared statistics under null and in our data
#Under NULL p-values are Uniformly distributed between 0 and 1,
#hence chisq-stats are expected to be:
expect.stats = qchisq(ppoints(nrow(pqtl_res)), df = 1, lower = F)
obs.stats = qchisq(pqtl_res$p_untreated, df = 1, lower = F)
lambda = median(obs.stats) / median(expect.stats) #GC lambda = ratio at medians

pdf(paste0(output_path,"QQ_plot_migration_GWAS_untreated_pvals.pdf "),
    width = 13, height = 10 )
qqplot(expect.stats, obs.stats, xlab = "chisq expected", ylab = "chisq observed",
       sub = substitute(paste(lambda, "=", lam), list(lam = signif(lambda,3))), 
       cex = 0.8, col = "red")
abline(0,1)
dev.off()

# Deflated QQ plot, see here for an explanation: 
# https://stats.stackexchange.com/questions/185284/deflated-qq-plots-in-genome-wide-association-studies
# number of significant results
# local FDR
sum(pqtl_res$q_untreated < 0.05) # 98033 -> 0
sum(pqtl_res$q_IFN < 0.05) # 127481 -> 0
sum(pqtl_res$q_LPS < 0.05) # 85566 -> 0

min(pqtl_res$p_untreated) # 1.041247e-05 (800k, -genoPCs) -> 3.021144e-06 (500k, + genoPCs)
min(pqtl_res$p_IFN) #8.189842e-06 (800k, -genoPCs) -> 1.020174e-05 (500k, + genoPCs)
min(pqtl_res$p_LPS)

# BH
sum(pqtl_res$p_BH_untreated < 0.05) # 63191 -> 113781 --> 0
sum(pqtl_res$p_BH_IFN < 0.05) # 82629  -> 18972 --> 0
sum(pqtl_res$p_BH_LPS < 0.05) # 52582 -> 7241 --> 0
min(pqtl_res$p_BH_untreated )

# Bonferroni
sum(pqtl_res$p_bonf_untreated < 0.05) # 241 -> 500 --> 0
sum(pqtl_res$p_bonf_IFN < 0.05) # 149 -> 42 --> 0
sum(pqtl_res$p_bonf_LPS < 0.05) # 92 -> 0 --> 0
min(pqtl_res$p_bonf_untreated ) # 1

# bonferroni-adjusted p-value distribution
p2 = ggplot(pqtl_res,aes(x = p_bonf_untreated)) + 
   geom_histogram(bins = 100) +
   theme_bw() +
   scale_y_log10() + 
   ylab("log10 counts")

p2

ggplot(pqtl_res[pqtl_res$p_bonf_untreated < 0.10,],aes(x = p_BH_untreated, y = p_bonf_untreated)) + 
   geom_point() + 
   geom_abline(slope = 1) + 
   theme_bw() 


# write all results for coloc
pqtl_res %>%
   readr::write_csv(.,"../../data/results/5.proportion_QTL/all_res_pqtl_migration.csv")

# filter to significant

sign_res_untreated = pqtl_res[pqtl_res$p_bonf_untreated < 0.05,c(1,grep(pattern = "untreated",colnames(pqtl_res)))]
sign_res_IFN = pqtl_res[pqtl_res$p_bonf_IFN < 0.05,c(1,grep(pattern = "IFN",colnames(pqtl_res)))]
sign_res_LPS= pqtl_res[pqtl_res$p_bonf_LPS < 0.05,c(1,grep(pattern = "LPS",colnames(pqtl_res)))]

# needed for plots
# # load minor allele dosage file for all donors
# sample_info = read.csv("../../data/all_pools_migration_sample_info.csv")
donor_prop_changes = read.csv("../../data/results/4.check_donor_proportions/line_prop_changes_averages_variances_1pct_filtered.csv")  %>%
   dplyr::filter(condition %in% c("3nM_C5a")) %>%
   dplyr::filter(!log_fraction_mean %in% c(NA,Inf,-Inf))

genotype = readr::read_csv("../../data/genotypes/full_genotype/genotype_minor_allele_dosage.csv") %>%
   dplyr::relocate(rn)
genotype = genotype %>% 
   dplyr::select(colnames(genotype)[colnames(genotype) %in% c("rn",unique(donor_prop_changes$donor))]) 


genotype_subset = genotype %>%
   dplyr::ungroup() %>%
   #dplyr::filter(rn %in% to_subset) %>%
   tidyr::pivot_longer(cols = !rn, values_to = "alt_allele_dosage", names_to = "donor")

length(unique(genotype_subset$rn)) # number of SNPs after filtering

genotype_subset = genotype_subset %>%
   dplyr::mutate(rn = stringr::str_replace(rn, "chr", "")) 

sign_eqtl = read.csv("../../../OTAR2065_sc_eQTL/data/results/4.Inspect_eQTL_results/tensorQTL_variant_gene_60PCs.csv")

sign_eqtl = sign_eqtl %>% 
   dplyr::filter(variant_id %in% unique(genotype_subset$rn) | swapped_tensor_alleles %in% unique(genotype_subset$rn))

genotype_subset_eqtl = genotype_subset %>%
   dplyr::filter(rn %in% sign_eqtl$variant_id | rn %in% sign_eqtl$swapped_tensor_alleles)

rm(genotype,genotype_subset)
## filter to eQTL

donor_prop_changes_IFN = donor_prop_changes %>%
   ungroup() %>%
   inner_join(.,genotype_subset_eqtl) %>%
   dplyr::filter(treatment == "IFN") %>%
   distinct()
donor_prop_changes_untreated = donor_prop_changes %>%
   ungroup() %>%
   left_join(.,genotype_subset_eqtl) %>%
   dplyr::filter(treatment == "Untreated")  %>%
   distinct()
donor_prop_changes_LPS = donor_prop_changes %>%
   ungroup() %>%
   left_join(.,genotype_subset_eqtl) %>%
   dplyr::filter(treatment == "LPS") %>%
   distinct()


pdf(file = "../../data/results/5.proportion_QTL/untreated_alele_dosage_testsnps_topSNP.pdf",
    width = 5, height = 5,
)
which.min(sign_res_untreated$p_bonf_untreated)
snp = sign_res_untreated[which.min(sign_res_untreated$p_bonf_untreated),"snp"]

# for(snp in names(sign_res_untreated)){
p = boxplot_pQTL(variant = snp,phenotype_with_genotype_info_df = donor_prop_changes_untreated,
                 qval = sign_res_untreated[sign_res_untreated$snp==snp,"p_bonf_untreated"],
                 coef = sign_res_untreated[sign_res_untreated$snp==snp,"coef_untreated"],color_by="pool")
plot(p)


dev.off()

pdf(file = "../../data/results/5.proportion_QTL/IFN_alele_dosage_testsnps_topSNP.pdf",
    width = 5, height = 5,
)
which.min(sign_res_IFN$p_bonf_IFN)
snp = sign_res_IFN[which.min(sign_res_IFN$p_bonf_IFN),"snp"]

# for(snp in names(sign_res_IFN)){
p = boxplot_pQTL(variant = snp,phenotype_with_genotype_info_df = donor_prop_changes_IFN,
                 qval = sign_res_IFN[sign_res_IFN$snp==snp,"p_bonf_IFN"],
                 coef = sign_res_IFN[sign_res_IFN$snp==snp,"coef_IFN"],color_by="pool")
plot(p)


dev.off()

pdf(file = "../../data/results/5.proportion_QTL/LPS_alele_dosage_testsnps_topSNP.pdf",
    width = 5, height = 5,
)
which.min(sign_res_LPS$p_bonf_LPS)
snp = sign_res_LPS[which.min(sign_res_LPS$p_bonf_LPS),"snp"]

p = boxplot_pQTL(variant = snp,phenotype_with_genotype_info_df = donor_prop_changes_LPS,
                 qval = sign_res_LPS[sign_res_LPS$snp==snp,"p_bonf_LPS"],
                 coef = sign_res_LPS[sign_res_LPS$snp==snp,"coef_LPS"],color_by="pool")
plot(p)


dev.off()

# save sign res
sign_res_IFN = sign_res_IFN %>%
   dplyr::rename(coefficient = coef_IFN, 
                 p = p_IFN, q = q_IFN, p_BH = p_BH_IFN,p_bonferroni=p_bonf_IFN ) %>%
   dplyr::mutate(treatment = "IFN")

sign_res_LPS = sign_res_LPS %>%
   dplyr::rename(coefficient = coef_LPS, p = p_LPS, q = q_LPS, p_BH = p_BH_LPS,
                 p_bonferroni=p_bonf_LPS) %>%
   dplyr::mutate(treatment = "LPS")

sign_res_untreated = sign_res_untreated %>%
   dplyr::rename(coefficient = coef_untreated, 
                 p = p_untreated, q = q_untreated, p_BH = p_BH_untreated,
                 p_bonferroni=p_bonf_untreated) %>%
   dplyr::mutate(treatment = "untreated")

sign_res_untreated %>%
   dplyr::rows_append(sign_res_IFN) %>%
   dplyr::rows_append(sign_res_LPS) %>%
   readr::write_csv(.,"../../data/results/5.proportion_QTL/sign_res_pqtl_migration.csv")

donor_prop_changes_IFN %>%
   dplyr::filter(rn == snp) %>%
   group_by(alt_allele_dosage) %>%
   dplyr::summarise(n = n())

donor_prop_changes_untreated %>%
   dplyr::filter(rn == snp) %>%
   group_by(alt_allele_dosage) %>%
   dplyr::summarise(n = n())

donor_prop_changes_LPS %>%
   dplyr::filter(rn == snp) %>%
   group_by(alt_allele_dosage) %>%
   dplyr::summarise(n = n())

donor_prop_changes_IFN %>%
   dplyr::filter(rn == snp) %>%
   summarise(Unique_Elements = n_distinct(donor))  # 23 unique donors

donor_prop_changes_untreated %>%
   dplyr::filter(rn == snp) %>%
   summarise(Unique_Elements = n_distinct(donor))  # 20 unique donors

donor_prop_changes_LPS %>%
   dplyr::filter(rn == snp) %>%
   summarise(Unique_Elements = n_distinct(donor))  # 20 unique donors

# preliminary look - proper way is coloc
sign_eqtl_untreated_pqtl = sign_eqtl[sign_eqtl$variant_id %in% sign_res_untreated$snp | sign_eqtl$swapped_tensor_alleles %in% sign_res_untreated$snp,] %>%
   dplyr::filter(treatment == "untreated")
sign_eqtl_IFN_pqtl = sign_eqtl[sign_eqtl$variant_id %in% sign_res_IFN$snp | sign_eqtl$swapped_tensor_alleles %in% sign_res_IFN$snp,] %>%
   dplyr::filter(treatment == "IFN")
sign_eqtl_LPS_pqtl = sign_eqtl[sign_eqtl$variant_id %in% sign_res_LPS$snp | sign_eqtl$swapped_tensor_alleles %in% sign_res_LPS$snp,] %>%
   dplyr::filter(treatment == "LPS")

# SLC2A13
pdf(file = "../../data/results/5.proportion_QTL/untreated_alele_dosage_testsnps_topSNP_SLC2A13-12_40147825_A_G.pdf",
    width = 5, height = 5,
)
snp = "12_40147825_A_G"

# for(snp in names(sign_res_untreated)){
p = boxplot_pQTL(variant = snp,phenotype_with_genotype_info_df = donor_prop_changes_untreated,
                 qval = sign_res_untreated[sign_res_untreated$snp==snp,"p_BH_untreated"],
                 coef = sign_res_untreated[sign_res_untreated$snp==snp,"coef_untreated"],color_by="pool")
plot(p)


dev.off()
## legacy code #########################
# 
# # replicates per pool and treatment
# sample_info %>%
#   filter(top_or_bottom  == "All") %>%
#   count(treatment, pool, sort = TRUE)
# 
# # add DNA input info for all
# 
# dna_info_grouped = sample_info %>% 
#   group_by(well) %>% 
#   summarise(DNA_input_ng = sum(DNA_input_ng, na.rm = T), top_or_bottom = "All") %>%
#   ungroup() %>%
#   unique() 
# 
# sample_info = left_join(sample_info, dna_info_grouped, by = c("well","top_or_bottom")) %>%
#   dplyr::mutate(DNA_input_ng = coalesce(DNA_input_ng.x, DNA_input_ng.y)) %>%
#   dplyr::select(-c(DNA_input_ng.x, DNA_input_ng.y))
# 
# w = list()
# for(name in sample_info[["ID"]]){
#   pool = sample_info %>%
#     filter(ID == name) %>%
#     .[["pool"]]
#   w[[name]] = read.table(paste0("../../data/",pool,"/w/",name,"_w_estimate.txt"))
#   colnames(w[[name]]) = c("donor","proportion")
#   w[[name]]$ID = name
# }
# 
# w = do.call("rbind",w)
# # convert slightly negative values to 0
# w[w$proportion<0, "proportion"] = 0
# 
# samples = merge(w, sample_info, all.y = TRUE) %>%
#   filter(!is.na(donor)) # sample P4_170 is misbehaving for some reason
# 
# samples$well = factor(samples$well, ordered = T, levels = unique(sort(samples$well)))
# samples$top_or_bottom = factor(samples$top_or_bottom, ordered = T, levels = c("All","Top","Bottom"))
# 
# 
# # Compare with errors from simulations
# # Assuming coverage of 2 (min mean depth for these samples = 2.5, shouldn't matter much in terms of max errors)
# w_estimate_sim_C2 =  read.table(
#   paste0(
#     "../../../simulation_WGS_exome_data/data/pool_19d_another/WGS/w/w_estimate_donor_proportions_fixed_100_trials_WGS_C2.txt.gz"))
# 
# colnames(w_estimate_sim_C2) = c("donor","w_est","trial","coverage")
# w_estimate_sim_C2$coverage = gsub("C","",w_estimate_sim_C2$coverage)
# w_real = read.table(
#   paste0(
#     "../../../simulation_WGS_exome_data/data/pool_19d_another/w/donor_proportions_fixed_10_sampling_reps_100_trials.txt.gz"), header = T)
# # retain only the first sampling replicate (only one taken to simulate WGS)
# w_real = w_real[w_real$sampling_rep==1,1:ncol(w_real)-1]
# w_real = w_real %>% tidyr::pivot_longer(!trial,values_to = "w_real",names_to = "donor")
# w_errors = inner_join(w_estimate_sim_C2, w_real, by = c("trial", "donor"))
# w_errors$w_diff = w_errors$w_real - w_errors$w_est
# 
# w_errors = w_errors %>% 
#   group_by(w_real) %>%
#   summarise(max_error = max(abs(w_diff)))
# 
# # add errors to our sample values depending on the proportion magnitude
# matched = sapply(samples$proportion, function(x) which.min(abs(w_errors$w_real - x)))
# samples[["proportion_error"]] = w_errors[matched,]$max_error
# 
# # adjust minimum and maximum proportions per well and donor (only top and bottom for now)
# min_adjusted_props = samples %>% 
#   filter(top_or_bottom != "All") %>%
#   group_by(well, donor) %>%
#   slice_min(proportion) %>%
#   mutate(adjusted_proportion = proportion + proportion_error) %>%
#   ungroup() %>%
#  dplyr::select(ID, donor, adjusted_proportion)
# 
# max_adjusted_props = samples %>% 
#   filter(top_or_bottom != "All") %>%
#   group_by(well, donor) %>%
#   slice_max(proportion) %>%
#   mutate(adjusted_proportion.y = proportion - proportion_error) %>%
#   mutate(adjusted_proportion.y = case_when(adjusted_proportion.y >=0 ~ adjusted_proportion.y,
#                                            adjusted_proportion.y < 0 ~ 0.001)) %>%
#   ungroup() %>%
#  dplyr::select(ID, donor, adjusted_proportion.y)
# 
# samples = samples %>%
#   left_join(min_adjusted_props) %>%
#   left_join(max_adjusted_props) %>%
#   mutate(adjusted_proportion = coalesce(adjusted_proportion, adjusted_proportion.y)) %>%
#  dplyr::select(-c(adjusted_proportion.y)) 
# 
# samples$donor_well = paste0(samples$donor,"_",samples$well)
# 
# 
# p1 = samples %>%
#   group_by(well,donor) %>%
#   summarise(bottom_top_fraction = adjusted_proportion[top_or_bottom == "Bottom"] / adjusted_proportion[top_or_bottom == "Top"],
#             treatment = treatment, condition = condition, pool = pool) %>%
#   dplyr::select(well:pool) %>%
#   unique() %>%
#   ggplot(., aes(bottom_top_fraction)) +
#   geom_histogram(bins = 100) +
#   geom_vline(xintercept = 1,col = "red") +
#   facet_wrap(vars(treatment), scales = "free_x") + 
#   theme_bw()
# p2 = samples %>%
#   group_by(well,donor) %>%
#   summarise(bottom_top_fraction = adjusted_proportion[top_or_bottom == "Bottom"] / adjusted_proportion[top_or_bottom == "Top"],
#             treatment = treatment, condition = condition, pool = pool) %>%
#   dplyr::select(well:pool) %>%
#   unique() %>%
#   dplyr::select(well:pool) %>%
#   unique() %>%
#   ggplot(., aes(log(bottom_top_fraction))) +
#   geom_histogram(bins = 10) +
#   geom_vline(xintercept = 0,col = "red") +
#   facet_wrap(vars(treatment), scales = "free_x") + 
#   theme_bw()
# 
# png("../../data/QTL_results/1.proportion_QTL/histogram_fractions_C5a_prefilter.png",
#     width = 9, height = 4, units = "in", res = 400)
# p1 / p2 + patchwork::plot_annotation(title = "Distribution of top/bottom fraction",
#                                      subtitle = "With C5a. As-is (top) or logged (bottom). Red line = no change.")
# 
# dev.off()
# 
# png("../../data/QTL_results/1.proportion_QTL/change_in_fractions_per_donor.png",
#     width = 13, height = 5, units = "in", res = 400)
# p = samples %>%
#   group_by(well,donor) %>%
#   summarise(bottom_top_fraction = adjusted_proportion[top_or_bottom == "Bottom"] / adjusted_proportion[top_or_bottom == "Top"],
#             treatment = treatment, condition = condition, pool = pool, proportion = proportion,
#             treatment_condition = paste(treatment, condition, sep = "_")) %>%
#   mutate(proportion = median(proportion) + 0.0000001) %>% # keeping proportion from "all" (should be the median of top and bottom)
#   ungroup() %>%
#   unique() %>%
#   ggplot( aes(x=donor,y=log(bottom_top_fraction), color=treatment_condition)) + 
#   geom_point(aes(size = proportion), binaxis='y', stackdir='center',position = position_dodge(), alpha = 0.6) + 
#   geom_hline(yintercept = 0) +
#   theme_bw() + 
#   theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust =1)) +
#   scale_size_continuous(range  = c(0.1, 8), 
#                         limits = c(0, 0.8), 
#                         breaks = c(0, 0.1, 0.2, 0.3,0.6)) +
#   ylab("log([bottom/top] proportions)") + 
#   ggtitle("Fractions per treatment + condition per donor") 
#   plot(p)
# dev.off()
# 
# png("../../data/QTL_results/1.proportion_QTL/change_in_fractions_per_donor_blackandwhite.png",
#     width = 13, height = 5, units = "in", res = 400)
# p = samples %>%
#   group_by(well,donor) %>%
#   summarise(bottom_top_fraction = adjusted_proportion[top_or_bottom == "Bottom"] / adjusted_proportion[top_or_bottom == "Top"],
#             treatment = treatment, condition = condition, pool = pool, proportion = proportion,
#             treatment_condition = paste(treatment, condition, sep = "_")) %>%
#   mutate(proportion = median(proportion) + 0.0000001) %>% # keeping proportion from "all" (should be the median of top and bottom)
#   ungroup() %>%
#   unique() %>%
#   ggplot( aes(x=donor,y=log(bottom_top_fraction), color=ifelse(proportion < 0.01,"grey","black"))) + 
#   geom_point(aes(size = proportion), binaxis='y', stackdir='center',position = position_dodge(), alpha = 0.9) + 
#   geom_hline(yintercept = 0) +
#   scale_color_identity() +
#   theme_bw() + 
#   theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust =1)) +
#   scale_size_continuous(range  = c(0.1, 8), 
#                         limits = c(0, 0.8), 
#                         breaks = c(0, 0.1, 0.2, 0.3,0.6)) +
#   ylab("log([bottom/top] proportions)") + 
#   ggtitle("Fractions where the total donor proportion is over (black) or under 1% per donor") 
# plot(p)
# dev.off()
# 
# 
# # correlations in within-pool replicates
# 
# to_plot = samples %>%
#   dplyr::filter(well %in% c(100,87))%>%
#   dplyr::group_by(well,donor) %>%
#   dplyr::summarise(bottom_top_fraction = adjusted_proportion[top_or_bottom == "Bottom"] / adjusted_proportion[top_or_bottom == "Top"]) %>%
#   ungroup() %>%
#   unique() %>%  
#   tidyr::pivot_wider(values_from = "bottom_top_fraction", names_from = c("well"), id_cols = c("donor"))
# 
# p1 = ggplot(to_plot,aes(x = `87`, y = `100`)) + 
#   geom_point() +
#   ggrepel::geom_text_repel(data = to_plot, aes(label = donor), 
#                   box.padding = unit(0.45, "lines")) +
#   geom_abline() +
#   ggpubr::stat_cor(p.accuracy = 0.001, r.accuracy = 0.01) +
#   theme_bw() + theme(title = element_text(size = 9)) +
#   ggtitle("bottom/top fraction correlations in replicates")
# 
# to_plot = samples %>%
#   dplyr::filter(well %in% c(100,87) &   top_or_bottom =="Bottom") %>%
#   ungroup() %>%
#   unique() %>%  
#   tidyr::pivot_wider(values_from = "proportion", names_from = c("well"), id_cols = c("donor"))
# 
# p2 = ggplot(to_plot,aes(x = `87`, y = `100`)) + 
#   geom_point() +
#   ggrepel::geom_text_repel(data = to_plot, aes(label = donor), 
#                            box.padding = unit(0.45, "lines")) +
#   geom_abline() +
#   ggpubr::stat_cor(p.accuracy = 0.001, r.accuracy = 0.01) +
#   theme_bw() + theme(title = element_text(size = 9)) +
#   ggtitle("Donor proportion (bottom well) correlations in replicates")
# 
# p3 = ggplot(to_plot,aes(x = `87`, y = `100`)) + 
#   geom_point() +
#   ggrepel::geom_text_repel(data = to_plot, aes(label = donor), 
#                            box.padding = unit(0.45, "lines")) +
#   scale_x_continuous(limits = c(0,0.05)) +
#   scale_y_continuous(limits = c(0,0.05)) +
#   geom_abline() +
#   ggpubr::stat_cor(p.accuracy = 0.001, r.accuracy = 0.01) +
#   theme_bw() + theme(title = element_text(size = 9)) +
#   ggtitle("Donor proportion (bottom well) correlations in replicates (prop<0.05)")
# 
# png("../../data/QTL_results/1.proportion_QTL/scatterplot_well_replicates_87_100_bottomprop.png", 
#     width = 13, height = 6, units = "in", res = 400)
# (p1 / p2) | p3
# 
# dev.off()
# 
# # previous p1 but filtering out donors with fewer than 1% proportion
# to_plot = samples %>%
#   dplyr::filter(well %in% c(100,87))%>%
#   dplyr::group_by(well,donor) %>%
#   dplyr::mutate(proportion_all = median(proportion)) %>%
#   dplyr::filter(proportion_all > 0.01) %>%
#   dplyr::summarise(bottom_top_fraction = adjusted_proportion[top_or_bottom == "Bottom"] / adjusted_proportion[top_or_bottom == "Top"]) %>%
#   ungroup() %>%
#   unique() %>%  
#   tidyr::pivot_wider(values_from = "bottom_top_fraction", names_from = c("well"), id_cols = c("donor"))
# 
# p1 = ggplot(to_plot,aes(x = `87`, y = `100`)) + 
#   geom_point() +
#   ggrepel::geom_text_repel(data = to_plot, aes(label = donor), 
#                            box.padding = unit(0.45, "lines")) +
#   geom_abline() +
#   ggpubr::stat_cor(p.accuracy = 0.001, r.accuracy = 0.01) +
#   theme_bw() + theme(title = element_text(size = 9)) +
#   ggtitle("bottom/top fraction correlations in replicates (samples < 1% prop.)")
# 
# png("../../data/QTL_results/1.proportion_QTL/scatterplot_well_replicates_87_100_bottomprop_fractions_under_onepercent.png", 
#     width = 6, height = 6, units = "in", res = 400)
# p1
# 
# dev.off()
# 
# 
# to_plot = samples %>%
#   dplyr::filter(well %in% c(138,139))%>%
#   dplyr::group_by(well,donor) %>%
#   dplyr::summarise(bottom_top_fraction = adjusted_proportion[top_or_bottom == "Bottom"] / adjusted_proportion[top_or_bottom == "Top"]) %>%
#   ungroup() %>%
#   unique() %>%  
#   tidyr::pivot_wider(values_from = "bottom_top_fraction", names_from = c("well"), id_cols = c("donor"))
# 
# p1 = ggplot(to_plot,aes(x = `138`, y = `139`)) + 
#   geom_point() +
#   ggrepel::geom_text_repel(data = to_plot, aes(label = donor), 
#                            box.padding = unit(0.45, "lines")) +
#   geom_abline() +
#   ggpubr::stat_cor(p.accuracy = 0.001, r.accuracy = 0.01) +
#   theme_bw() + theme(title = element_text(size = 9)) +
#   ggtitle("bottom/top fraction correlations in replicates")
# 
# to_plot = samples %>%
#   dplyr::filter(well %in% c(138,139) &   top_or_bottom =="Bottom") %>%
#   ungroup() %>%
#   unique() %>%  
#   tidyr::pivot_wider(values_from = "proportion", names_from = c("well"), id_cols = c("donor"))
# 
# p2 = ggplot(to_plot,aes(x = `138`, y = `139`)) + 
#   geom_point() +
#   ggrepel::geom_text_repel(data = to_plot, aes(label = donor), 
#                            box.padding = unit(0.45, "lines")) +
#   geom_abline() +
#   ggpubr::stat_cor(p.accuracy = 0.001, r.accuracy = 0.01) +
#   theme_bw() + theme(title = element_text(size = 9)) +
#   ggtitle("Donor proportion (bottom well) correlations in replicates")
# 
# p3 = ggplot(to_plot,aes(x = `138`, y = `139`)) + 
#   geom_point() +
#   ggrepel::geom_text_repel(data = to_plot, aes(label = donor), 
#                            box.padding = unit(0.45, "lines")) +
#   scale_x_continuous(limits = c(0,0.05)) +
#   scale_y_continuous(limits = c(0,0.05)) +
#   geom_abline() +
#   ggpubr::stat_cor(p.accuracy = 0.001, r.accuracy = 0.01) +
#   theme_bw() + theme(title = element_text(size = 9)) +
#   ggtitle("Donor proportion (bottom well) correlations in replicates (prop<0.05)")
# 
# png("../../data/QTL_results/1.proportion_QTL/scatterplot_well_replicates_138-9_bottomprop.png", 
#     width = 13, height = 6, units = "in", res = 400)
# (p1 / p2) | p3
# 
# dev.off()
# 
# # previous p1 but filtering out donors with fewer than 1% proportion
# to_plot = samples %>%
#   dplyr::filter(well %in% c(138,139))%>%
#   dplyr::group_by(well,donor) %>%
#   dplyr::mutate(proportion_all = median(proportion)) %>%
#   dplyr::filter(proportion_all > 0.01) %>%
#   dplyr::summarise(bottom_top_fraction = adjusted_proportion[top_or_bottom == "Bottom"] / adjusted_proportion[top_or_bottom == "Top"]) %>%
#   ungroup() %>%
#   unique() %>%  
#   tidyr::pivot_wider(values_from = "bottom_top_fraction", names_from = c("well"), id_cols = c("donor"))
# 
# p1 = ggplot(to_plot,aes(x = `138`, y = `139`)) + 
#   geom_point() +
#   ggrepel::geom_text_repel(data = to_plot, aes(label = donor), 
#                            box.padding = unit(0.45, "lines")) +
#   geom_abline() +
#   ggpubr::stat_cor(p.accuracy = 0.001, r.accuracy = 0.01) +
#   theme_bw() + theme(title = element_text(size = 9)) +
#   ggtitle("bottom/top fraction correlations in replicates (samples < 1% prop.)")
# 
# png("../../data/QTL_results/1.proportion_QTL/scatterplot_well_replicates_138-9_bottomprop_fractions_under_onepercent.png", 
#     width = 6, height = 6, units = "in", res = 400)
# p1
# 
# dev.off()
# 
# # correlations all wells for shared donors
# cor_res = samples %>%
#   dplyr::filter(top_or_bottom =="Bottom") %>%
#   dplyr::mutate(well_treatment_condition_pool = paste(well, treatment, condition, pool,sep = "_")) %>%
#   dplyr::select(donor,well_treatment_condition_pool,proportion) %>%
#   tidyr::pivot_wider(names_from  = well_treatment_condition_pool,values_from = proportion,id_cols = donor) %>%
#   dplyr::arrange(donor) %>%
#   dplyr::select(!donor) %>%
#   tidyr::drop_na() %>%
#   corr.test()
# 
# png("../../data/QTL_results/1.proportion_QTL/correlation_donor_proportions_shared_donors_bottomprop.png", 
#     width = 10, height = 10, units = "in", res = 400)
# p1 = nicer_heatmap(cor_res)
# cowplot::plot_grid(p1) 
# 
# dev.off()
# 
# # correlations all wells for same pool
# cor_res_pool2 = samples %>%
#   dplyr::filter(top_or_bottom =="Bottom" & pool == "pool2") %>%
#   dplyr::mutate(well_treatment_condition_pool = paste(well, treatment, condition, pool,sep = "_")) %>%
#   dplyr::select(donor,well_treatment_condition_pool,proportion) %>%
#   tidyr::pivot_wider(names_from  = well_treatment_condition_pool,values_from = proportion,id_cols = donor) %>%
#   dplyr::arrange(donor) %>%
#   dplyr::select(!donor) %>%
#   corr.test()
# 
# png("../../data/QTL_results/1.proportion_QTL/correlation_donor_proportions_pool2_alldonors_bottomprop.png", 
#     width = 6, height = 6, units = "in", res = 400)
# p1 = nicer_heatmap(cor_res_pool2)
# cowplot::plot_grid(p1) 
# dev.off()
# 
# cor_res_pool3 = samples %>%
#   dplyr::filter(top_or_bottom =="Bottom" & pool == "pool3") %>%
#   dplyr::mutate(well_treatment_condition_pool = paste(well, treatment, condition, pool,sep = "_")) %>%
#   dplyr::select(donor,well_treatment_condition_pool,proportion) %>%
#   tidyr::pivot_wider(names_from  = well_treatment_condition_pool,values_from = proportion,id_cols = donor) %>%
#   dplyr::arrange(donor) %>%
#   dplyr::select(!donor) %>%
#   corr.test()
# 
# png("../../data/QTL_results/1.proportion_QTL/correlation_donor_proportions_pool3_alldonors_bottomprop.png", 
#     width = 6, height = 6, units = "in", res = 400)
# p1 = nicer_heatmap(cor_res_pool3)
# cowplot::plot_grid(p1) 
# dev.off()
# 
# cor_res_pool4 = samples %>%
#   dplyr::filter(top_or_bottom =="Bottom" & pool == "pool4") %>%
#   dplyr::mutate(well_treatment_condition_pool = paste(well, treatment, condition, pool,sep = "_")) %>%
#   dplyr::select(donor,well_treatment_condition_pool,proportion) %>%
#   tidyr::pivot_wider(names_from  = well_treatment_condition_pool,values_from = proportion,id_cols = donor) %>%
#   dplyr::arrange(donor) %>%
#   dplyr::select(!donor) %>%
#   corr.test()
# 
# png("../../data/QTL_results/1.proportion_QTL/correlation_donor_proportions_pool4_alldonors_bottom_prop.png", 
#     width = 6, height = 6, units = "in", res = 400)
# p1 = nicer_heatmap(cor_res_pool4)
# cowplot::plot_grid(p1) 
# dev.off()
# 
# 
# ### same for unfiltered bottom/top fractions
# 
#   cor_res = samples %>%
#     dplyr::group_by(well,donor) %>%
#     dplyr::summarise(bottom_top_fraction = adjusted_proportion[top_or_bottom == "Bottom"] / adjusted_proportion[top_or_bottom == "Top"],
#                      treatment = treatment, condition = condition, pool = pool) %>%
#     ungroup() %>%
#     unique() %>%
#     dplyr::mutate(well_treatment_condition_pool = paste(well, treatment, condition, pool,sep = "_")) %>%
#   dplyr::select(donor,well_treatment_condition_pool,bottom_top_fraction) %>%
#   tidyr::pivot_wider(names_from  = well_treatment_condition_pool,values_from = bottom_top_fraction,id_cols = donor) %>%
#   dplyr::arrange(donor) %>%
#   dplyr::select(!donor) %>%
#   tidyr::drop_na() %>%
#   corr.test()
# 
# png("../../data/QTL_results/1.proportion_QTL/correlation_donor_bottom_top_fractions_shared_donors.png", 
#     width = 10, height = 10, units = "in", res = 400)
# p1 = nicer_heatmap(cor_res)
# cowplot::plot_grid(p1) 
# dev.off()
# 
# # correlations all wells for same pool
# cor_res_pool2 = samples %>%
#   dplyr::filter( pool == "pool2") %>%
#   dplyr::group_by(well,donor) %>%
#   dplyr::summarise(bottom_top_fraction = adjusted_proportion[top_or_bottom == "Bottom"] / adjusted_proportion[top_or_bottom == "Top"],
#                    treatment = treatment, condition = condition, pool = pool) %>%
#   ungroup() %>%
#   unique() %>%
#   dplyr::mutate(well_treatment_condition_pool = paste(well, treatment, condition, pool,sep = "_")) %>%
#   dplyr::select(donor,well_treatment_condition_pool,bottom_top_fraction) %>%
#   tidyr::pivot_wider(names_from  = well_treatment_condition_pool,values_from = bottom_top_fraction,id_cols = donor) %>%
#   dplyr::arrange(donor) %>%
#   dplyr::select(!donor) %>%
#   tidyr::drop_na() %>%
#   corr.test()
# 
# png("../../data/QTL_results/1.proportion_QTL/correlation_donor_pool2_bottom_top_fractions.png", 
#     width = 6, height = 6, units = "in", res = 400)
# p1 = nicer_heatmap(cor_res_pool2)
# cowplot::plot_grid(p1) 
# dev.off()
# 
# cor_res_pool3 = samples %>%
#   dplyr::filter( pool == "pool3") %>%
#   dplyr::group_by(well,donor) %>%
#   dplyr::summarise(bottom_top_fraction = adjusted_proportion[top_or_bottom == "Bottom"] / adjusted_proportion[top_or_bottom == "Top"],
#                    treatment = treatment, condition = condition, pool = pool) %>%
#   ungroup() %>%
#   unique() %>%
#   dplyr::mutate(well_treatment_condition_pool = paste(well, treatment, condition, pool,sep = "_")) %>%
#   dplyr::select(donor,well_treatment_condition_pool,bottom_top_fraction) %>%
#   tidyr::pivot_wider(names_from  = well_treatment_condition_pool,values_from = bottom_top_fraction,id_cols = donor) %>%
#   dplyr::arrange(donor) %>%
#   dplyr::select(!donor) %>%
#   tidyr::drop_na() %>%
#   corr.test()
# 
# png("../../data/QTL_results/1.proportion_QTL/correlation_donor_pool3_bottom_top_fractions.png", 
#     width = 6, height = 6, units = "in", res = 400)
# p1 = nicer_heatmap(cor_res_pool3)
# cowplot::plot_grid(p1) 
# dev.off()
# 
# cor_res_pool4 = samples %>%
#   dplyr::filter(pool == "pool4") %>%
#   dplyr::group_by(well,donor) %>%
#   dplyr::summarise(bottom_top_fraction = adjusted_proportion[top_or_bottom == "Bottom"] / adjusted_proportion[top_or_bottom == "Top"],
#                    treatment = treatment, condition = condition, pool = pool) %>%
#   ungroup() %>%
#   unique() %>%
#   dplyr::mutate(well_treatment_condition_pool = paste(well, treatment, condition, pool,sep = "_")) %>%
#   dplyr::select(donor,well_treatment_condition_pool,bottom_top_fraction) %>%
#   tidyr::pivot_wider(names_from  = well_treatment_condition_pool,values_from = bottom_top_fraction,id_cols = donor) %>%
#   dplyr::arrange(donor) %>%
#   dplyr::select(!donor) %>%
#   tidyr::drop_na() %>%
#   corr.test()
# 
# png("../../data/QTL_results/1.proportion_QTL/correlation_donor_pool4_bottom_top_fractions.png", 
#     width = 6, height = 6, units = "in", res = 400)
# p1 = nicer_heatmap(cor_res_pool4)
# cowplot::plot_grid(p1) 
# dev.off()
# 
# # bottom / top proportions - do not include samples that have < 1% (0.01) difference between bottom and top
# # for adjusted proportions
# # do not include donors with 0% - 1% in "All" (samples with such low proportions may cause trouble in the fractions)
# to_keep = samples[samples$top_or_bottom == "All" & samples$proportion >= 0.01,"donor_well"]
# 
# bottom_top_samples = samples %>%
#   mutate(donor_well = paste0(donor,"_",well)) %>%
#   group_by(well,donor) %>%
#   filter(donor_well %in% to_keep) %>%
#   mutate(difference_max_min = max(adjusted_proportion,na.rm = T) - min(adjusted_proportion, na.rm = T)) %>%
#   summarise(bottom_top_fraction = adjusted_proportion[top_or_bottom == "Bottom"] / adjusted_proportion[top_or_bottom == "Top"],
#             treatment = treatment, condition = condition, pool = pool) %>%
#   unique() %>%
#   ungroup()
# 
# 
# 
# #### load genotype data
# # right now I'm only loading the dosage file with the largest number of shared variants covered by sequencing
# # I should include all non-identical SNP positions across all donors 
# pool2 = read_csv("../../data/pool2/genotypes/99_genotype_minor_allele_dosage.csv")
# pool3 = read_csv("../../data/pool3/genotypes/P3_138_genotype_minor_allele_dosage.csv")
# pool4 = read_csv("../../data/pool4/genotypes/P4_161_genotype_minor_allele_dosage.csv")
# 
# # joining 3 pools
# genotype = inner_join(pool2,pool3)
# # check genotypes at same positions for same donors are identical
# identical(genotype[genotype$rn, ]$aowh_2,pool2[genotype$rn, ]$aowh_2)
# identical(pool2[genotype$rn , ]$aowh_2,pool3[genotype$rn, ]$aowh_2)
# 
# genotype = inner_join(genotype,pool4) %>%
#   relocate(rn)
# # check genotypes at same positions for same donors are identical
# identical(genotype[genotype$rn, ]$hegp_3,pool2[genotype$rn, ]$hegp_3)
# identical(pool2[genotype$rn , ]$hegp_3,pool3[genotype$rn, ]$hegp_3)
# identical(genotype[genotype$rn , ]$hegp_3,pool4[genotype$rn, ]$hegp_3)
# 
# # remove genotypes that don't pass filters
# # more than 20% of samples in different bins ( for 53 donors that's 11)
# to_subset = genotype %>% 
#   ungroup() %>%
#   dplyr::select(!rn) %>% 
#   dplyr::mutate(hom_0 = rowSums(. == 0), hom_1 = rowSums(. == 1), het = rowSums(. == 0.5), 
#          threshold = ceiling(0.2 * ncol(.)),rn = genotype$rn) %>%
#   dplyr::filter(hom_0 > threshold & hom_1 > threshold & het > threshold) %>%
#   .$rn
# 
# genotype_subset = genotype %>%
#   ungroup() %>%
#   filter(rn %in% to_subset) %>%
#   pivot_longer(cols = !rn, values_to = "alt_allele_dosage", names_to = "donor")
# 
# length(unique(genotype_subset$rn)) # number of SNPs after filtering
# 
# subset_IFN = bottom_top_samples %>%
#   ungroup() %>%
#   inner_join(.,genotype_subset) %>%
#   filter(treatment == "IFN")
# 
# subset_untreated = bottom_top_samples %>%
#   ungroup() %>%
#   left_join(.,genotype_subset) %>%
#   filter(treatment == "Untreated")
# 
# # test mixed model
# # checking distribution again
# 
# p1 = rbind(subset_IFN,subset_untreated) %>%
#   dplyr::select(well:pool) %>%
#   unique() %>%
#   ggplot(., aes(bottom_top_fraction)) +
#   geom_histogram(bins = 10) +
#   geom_vline(xintercept = 1,col = "red") +
#   facet_wrap(vars(treatment), scales = "free_x") + 
#   theme_bw()
# p2 = rbind(subset_IFN,subset_untreated) %>%
#   dplyr::select(well:pool) %>%
#   unique() %>%
#   ggplot(., aes(log(bottom_top_fraction))) +
#   geom_histogram(bins = 10) +
#   geom_vline(xintercept = 0,col = "red") +
#   facet_wrap(vars(treatment), scales = "free_x") + 
#   theme_bw()
# 
# png("../../data/QTL_results/1.proportion_QTL/histogram_fractions_C5a_postfilter.png",
#     width = 9, height = 4, units = "in", res = 400)
# p1 / p2 + patchwork::plot_annotation(title = "Distribution of top/bottom fraction",
#                                      subtitle = "With C5a. As-is (top) or logged (bottom). Red line = no change.")
# 
# dev.off()
# 
# # subset of genotype first
# res_IFN = list()
# res_IFN_afex = list()
# res_untreated = list()
# # 2 secs per SNP
# for(snp in genotype_subset[1:3000,]$rn){
# 
#   res_IFN_afex[[snp]] = afex::mixed(bottom_top_fraction ~ alt_allele_dosage + (1 | donor) +( 1 | pool), 
#                                     subset_IFN[subset_IFN$rn == snp,], 
#                                     control = lmerControl(optCtrl = list(maxfun = 1e6)))
#     
#   ### classic lm
#   res_IFN[[snp]] = lm(formula = bottom_top_fraction ~ alt_allele_dosage,
#                data = subset_IFN[subset_IFN$rn == snp,] )
#   res_untreated[[snp]] = lm(formula = bottom_top_fraction ~ alt_allele_dosage,
#                       data = subset_untreated[subset_untreated$rn == snp,] )
# 
# }
# 
# save.image( file = "~/proportion_allele_effect.RData")
# load("~/proportion_allele_effect.RData")
# #unlink("~/proportion_allele_effect.RData")
# # check significant SNPs and plot them
# # how to use the media only vs C5a effect?
# # should alt allele dosage be a factor?
# res_IFN_afex[[snp]]$anova_table # effect of allele dosage vs reduced model w. intercept
# res_IFN_afex[[snp]]$full_model #fit of full model
# summary(res_IFN_afex[[snp]])$varcor 
# summary(res_IFN[[snp]])
# anova(res_IFN[[snp]])
# # extract p-value for allelic dosage coefficient
# p_IfN_afex  = sapply(res_IFN_afex, FUN = function(x) x$anova_table$`Pr(>F)`[1])
# p_IFN = sapply(res_IFN,FUN = function(x) anova(x)$`Pr(>F)`[1])
# p_untreated = sapply(res_untreated,FUN = function(x) anova(x)$`Pr(>F)`[1])
# 
# sum(p_IfN_afex < 0.05)
# sum(p_IFN < 0.05)
# sum(p_untreated < 0.05)
# 
# res_untreated = res_untreated[p_untreated < 0.05]
# res_IFN = res_IFN[p_IFN < 0.05]
# # car::Anova(res_IFN)# need to provide design matrix
# # plot(res_IFN)
# pdf(file = "../../data/QTL_results/1.proportion_QTL/untreated_alele_dosage_3000snps.pdf",
#     width = 9, height = 5,
#     )
# 
# for(snp in names(res_untreated)){
#   p1 = subset_untreated %>% 
#     filter(rn == snp) %>%
#     ggplot(aes(x=alt_allele_dosage, y=bottom_top_fraction,col = well, shape=pool)) + 
#     scale_x_continuous(breaks = c(0,0.5,1), limits = c(-0.1,1.1)) + 
#     geom_beeswarm(priority='random',cex=3.5, groupOnX=T) +
#     guides(col=guide_legend(ncol=2)) +
#     theme_bw()
#   p2 = subset_untreated %>% 
#     filter(rn == snp) %>%
#     ggplot(aes(x=alt_allele_dosage, y=bottom_top_fraction,col = donor, shape=pool)) + 
#     scale_x_continuous(breaks = c(0,0.5,1), limits = c(-0.1,1.1)) + 
#     geom_beeswarm(priority='random',cex=3.5, groupOnX=T) +
#     guides(col=guide_legend(ncol=2)) +
#     theme_bw()
#   
#   p = p1 + p2 + patchwork::plot_annotation(title = paste0("Allelic effects at SNP ", snp),
#                                        subtitle = paste0("In ", unique(subset_untreated$treatment), 
#                                                          " + ",unique(subset_untreated$condition),
#                                                          ". Unadjusted p-val = ",
#                                                          p_untreated[[snp]] ))
#   plot(p)
# }
# 
# dev.off()
# 
# pdf(file = "../../data/QTL_results/1.proportion_QTL/IFN_alele_dosage_3000snps.pdf",
#     width = 9, height = 5,
# )
# 
# for(snp in names(res_IFN)){
#   
#   p1 = subset_IFN %>% 
#     filter(rn == snp) %>%
#     ggplot(.,aes(x=alt_allele_dosage, y=bottom_top_fraction,col = well, shape=pool)) + 
#     scale_x_continuous(breaks = c(0,0.5,1), limits = c(-0.1,1.1)) + 
#     geom_beeswarm(priority='random',cex=3.5, groupOnX=T) +
#     guides(col=guide_legend(ncol=2)) +
#     theme_bw()
#   p2 = subset_IFN %>% 
#     filter(rn == snp) %>%
#     ggplot(.,aes(x=alt_allele_dosage, y=bottom_top_fraction,col = donor, shape=pool)) + 
#     scale_x_continuous(breaks = c(0,0.5,1), limits = c(-0.1,1.1)) + 
#     geom_beeswarm(priority='random',cex=3.5, groupOnX=T) +
#     guides(col=guide_legend(ncol=2)) +
#     theme_bw()
#   
#   p = p1 + p2 + patchwork::plot_annotation(title = paste0("Allelic effects at SNP ", snp),
#                                            subtitle = paste0("In ", unique(subset_IFN$treatment), 
#                                                              " + ",unique(subset_IFN$condition),
#                                                              ". Unadjusted p-val = ",
#                                                              p_IFN[[snp]] ))
#   plot(p)
# }
# 
# dev.off()
# 
# subset_IFN %>% 
#   dplyr::filter(rn == snp) %>%
#   group_by(alt_allele_dosage) %>%
#   dplyr::summarise(n = n())
# 
# subset_untreated %>% 
#   dplyr::filter(rn == snp) %>%
#   group_by(alt_allele_dosage) %>%
#   dplyr::summarise(n = n())
# 
# subset_IFN %>% 
#   dplyr::filter(rn == snp) %>%
#   summarise(Unique_Elements = n_distinct(donor))  # 23 unique donors
# 
# subset_untreated %>% 
#   dplyr::filter(rn == snp) %>%
#   summarise(Unique_Elements = n_distinct(donor))  # 20 unique donors
# ##
# ## mean/median top/bottom condition x treatment
# # color by treatment
