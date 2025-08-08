# 8.5. Mendelian randomisation with eQTL results and GWAS

library(tidyverse)
library(mrSampleOverlap)
library(TwoSampleMR)
source("./functions.R")
options(scipen = 999) # prevents from showing very small numbers as 0


output_path = "../../data/results/8.6.MR_eQTL_pheno_bias_correction"
dir.create(output_path,recursive = TRUE)

TWMR_p = readRDS("../../data/results/8.5.eQTL_MR/TWMR/output/phagocytosis_twmr_p.rds")
TWMR_res = readRDS("../../data/results/8.5.eQTL_MR/TWMR/output/phagocytosis_twmr_res.rds")

TWMR_res[["LPS_LRRK2"]]
# only one gene: itself
treatment = "LPS"
focus_gene = "LRRK2"
phenotype = "phagocytosis"

# get variants used for TWMR calculations

snps = data.frame("SNP"=names(TWMR_res[[paste0(treatment,"_",focus_gene)]]$LD_matrix),
                  "A1"="A1","A2"="A2")
snps$A2 = stringr::str_split_i(snps$SNP,pattern = "_",3)
snps$A1 = stringr::str_split_i(snps$SNP,pattern = "_",4)

# reading input and output files


load_nominal_results = function(treatment) {
  nominal_res = readr::read_delim(
    paste0(
      "../../data/results/tensorqtl/best_results/sum_sizefactorsNorm_log2_scaled_centered_",
      treatment,
      "_Not_proliferating_common_500kb_window_tensorQTL_nominal.txt"
    )
  )
  
  
  # mapping ensembl to symbol
  ensembl_to_symbol = ensembl_to_gene_name(
    IDs = unique(nominal_res$phenotype_id),
    IDFrom = "GENEID",
    IDTo = "SYMBOL"
  )$map %>%
    dplyr::rename(phenotype_id = From, gene = To) %>%
    as_tibble() %>%
    dplyr::distinct(phenotype_id, .keep_all = TRUE)
  
  nominal_res = nominal_res %>%
    dplyr::left_join(ensembl_to_symbol)
  # adding info and filtering
  return(nominal_res)
  
}

# for(treatment in c("untreated","IFN","LPS")){
#   nominal_res = load_nominal_results(treatment)
# 
#   # saving chr pos ref alt to query dbsnp file later
#   nominal_res %>%
#     dplyr::mutate(chr = stringr::str_split_i(variant_id,pattern = "_",i=1),
#                   pos = stringr::str_split_i(variant_id,pattern = "_",i=2),
#                   ref = stringr::str_split_i(variant_id,pattern = "_",i=3),
#                   alt = stringr::str_split_i(variant_id,pattern = "_",i=4)) %>%
#     dplyr::select(chr,pos) %>%
#     dplyr::distinct() %>%
#     readr::write_tsv(paste0(output_path,"/",treatment,"_pos_to_subset_dbsnp151.txt"),col_names = FALSE)
#   rm(nominal_res)
#   gc()
# }

## load dbsnp files subset from nominal results - GRCh38
# dbsnp 151 because it's the one used for the 1000G weight files
# dbsnp = readr::read_delim(paste0("../../data/results/8.6.MR_eQTL_pheno_bias_correction/for_rsid_matching/",
# treatment,"_nominal.tsv"), col_names = FALSE,delim = " ") %>%
#   dplyr::rename(chr = X1,
#                 pos = X2,
#                 rsid = X3,
#                 ref = X4,
#                 ALTS = X5)
# exposures
nominal_res = load_nominal_results(treatment)
  
# nominal results already tests for all variants 500kb around the gene
exposure = nominal_res %>%
   dplyr::filter(gene == focus_gene) %>%
  dplyr::filter(variant_id %in% snps$SNP) %>%
  dplyr::mutate(
    chr=stringr::str_split_i(variant_id,pattern = "_",1),
    pos=stringr::str_split_i(variant_id,pattern = "_",2),
    ref=stringr::str_split_i(variant_id,pattern = "_",3),
    alt=stringr::str_split_i(variant_id,pattern = "_",4),
    info = 1, # dummy (is this used?)
    z = slope / slope_se,
    minuslog10p = -log10(pval_nominal),
    samplesize.exposure = ifelse(treatment == "untreated",189,188)
  ) %>%
  dplyr::rename(beta = slope,
                se = slope_se,
                pval.exposure = pval_nominal) %>%
  dplyr::mutate(pos = as.double(pos),
                chr = as.double(chr)) %>%
  dplyr::rename(SNP = variant_id,
                beta.exposure = beta,
                se.exposure = se) %>%
  dplyr::mutate(id.exposure =  paste0("eQTL_", treatment), 
                exposure = paste0("eQTL_", treatment),
                effect_allele.exposure = alt, 
                other_allele.exposure = ref, 
                eaf.exposure = af,
                units.exposure = "SD units") 

# outcomes
# extract GWAS beta for all variants

pheno_gwas = readr::read_csv(
  file = paste0(
    "../../../OTAR2065_phenotypic_QTLs/data/results/",
    phenotype,
    "/2.check_association_results/lm_1pct_filtered_deflated/all_res_pqtl_",
    phenotype,
    ".csv"
  )
)

outcome = pheno_gwas %>%
  dplyr::filter(snp %in% exposure$SNP) %>%
  dplyr::select(contains(treatment), snp) %>%
  dplyr::rename(beta = paste0("coef_", treatment),
                se = paste0("se_", treatment),
                pval.outcome = paste0("p_", treatment),
                SNP = snp
                ) %>%
  dplyr::mutate(beta = tidyr::replace_na(beta, 0), #  NA variant effects are 0
                chr=stringr::str_split_i(SNP,pattern = "_",1),
                pos=stringr::str_split_i(SNP,pattern = "_",2),
                ref=stringr::str_split_i(SNP,pattern = "_",3),
                alt=stringr::str_split_i(SNP,pattern = "_",4),
                info = 1, # dummy (is this used?)
                z = beta / se,
                samplesize.outcome = case_when(treatment == "untreated"~ 133,
                              treatment == "IFN" ~ 137,
                              treatment == "LPS" ~ 140))  %>%
  dplyr::mutate(pos = as.double(pos),
                chr = as.double(chr)) %>%
  dplyr::left_join(exposure[,c("SNP","alt","af")], by = c("SNP","alt")) %>% # taking af from larger eQTL analysis (same pop so shouldn't matter)
  dplyr::rename(beta.outcome = beta,
                se.outcome = se,
                effect_allele.outcome = alt, 
                other_allele.outcome = ref, 
                eaf.outcome = af) %>%
  dplyr::mutate(id.outcome = "phagocytosis_LPS", ,
                outcome = "phagocytosis_LPS",
                units.outcome = "SD units") %>%
  dplyr::distinct()

# inputs must be dataframes with the following columns:
# outcome_dat: SNP, beta.outcome, se.outcome, effect_allele.outcome, other_allele.outcome, eaf.outcome, outcome
# exposure_dat: SNP, beta.exposure, se.exposure, effect_allele.exposure, other_allele.exposure, eaf.exposure
#  

dat_harmonized = TwoSampleMR::harmonise_data(exposure, outcome) %>%
  dplyr::filter(mr_keep == TRUE)

mr_results = TwoSampleMR::mr(dat_harmonized)

dat_summarized <- dat_harmonized %>%
  TwoSampleMR::add_rsq() %>%
  dplyr::group_by(exposure) %>%
  dplyr::summarize(rsq_exposure = sum(rsq.exposure),
                   n_variants = n(), 
                   samplesize_exposure = max(samplesize.exposure),
                   samplesize_outcome = max(samplesize.outcome))

grid <- tidyr::crossing(overlap_prop = 1, # overlap is always 100%
                        ols_bias = seq(0, 1, 0.2))


bias_res <- dat_summarized %>%
  tidyr::crossing(grid) %>%
  dplyr::mutate(res = mrSampleOverlap::estimate_overlap_bias(samplesize_exposure = samplesize_exposure, 
                                            samplesize_outcome = samplesize_outcome, 
                                            n_variants = n_variants, 
                                            rsq_exposure = rsq_exposure, 
                                            overlap_prop = overlap_prop, 
                                            ols_bias = ols_bias)) %>%
  unnest(res)

# calculating bias under the null as in https://onlinelibrary.wiley.com/doi/10.1002/gepi.21998
# formula (3): bias under the null = OLS estimate bias x (Percentage overlap/100) x relative bias
# where relative bias. = 1/E(F)

# The following function is unhelpfully called estimate_F but looking inside is the expectation
# E_F = mrSampleOverlap::estimate_f(samplesize_exposure = dat_summarized$samplesize_exposure,
#                                   n_variants = dat_summarized$n_variants,
#                                   rsq_exposure = dat_summarized$rsq_exposure,
#                                   lci_95 = TRUE) # returning a more conservative estimate 
# # A more liberal estimation per IV:
# #calculate_F_statistic_MR(pval = exposure$pval.exposure[1],samplesize = exposure$samplesize.exposure[1])
# # over 10 is considered strong
# 
# # We know our OLS:the ordinary least squares (OLS, also known as standard least squares regression) 
# # estimate is obtained by regressing the outcome on the risk factor
# fit = lm(beta.outcome ~ beta.exposure,data = dat_harmonized)
# OLS_estimate = fit$coefficients["beta.exposure"]
# # how is the OLS bias calculated?
# bias = unname(OLS_estimate_bias * 1 * (1/E_F))

# I'm assuming the worst case scenario of OLS bias
bias_res %>%
  split_exposure() %>%
  pivot_longer(cols = c(bias, type1_error)) %>%
  ggplot(aes(x = ols_bias, y = value, color = as.character(ols_bias))) +
  geom_point() +
  facet_grid(rows = vars(exposure), 
             cols = vars(name),
             scales = "free_y") +
  labs(x = "Bias in observational estimate",
       y = "Value") +
  scale_color_discrete(name = "Bias in \nObservational \nEstimate") +
  ggpubr::theme_pubr() +
  ggtitle(focus_gene)

### example data of this bias correction method

# extract genetic instruments for BMI
bmi_file <- system.file("extdata", "bmi.txt", package = "TwoSampleMR")
bmi_exp_dat <- read_exposure_data(bmi_file)

cardio_file = system.file("extdata", "bmi.txt", package = "TwoSampleMR")
cardio_out_data = read_outcome_data(cardio_file)

# harmonize effect alleles, and keep only alleles present in both exposure and outcome data
dat_harmonized = TwoSampleMR::harmonise_data(bmi_exp_dat, cardio_out_data) %>%
  dplyr::filter(mr_keep == TRUE)

mr_results = TwoSampleMR::mr(dat_harmonized)

dat_summarized <- dat_harmonized %>%
  add_rsq() %>%
  dplyr::group_by(exposure) %>%
  dplyr::summarize(rsq_exposure = sum(rsq.exposure),
                   n_variants = n(), 
                   samplesize_exposure = max(samplesize.exposure),
                   samplesize_outcome = max(samplesize.outcome))

grid <- tidyr::crossing(overlap_prop = seq(0, 1, 0.1),
                        ols_bias = seq(0, 1, 0.2))

bias_res <- dat_summarized %>%
  tidyr::crossing(grid) %>%
  dplyr::mutate(res = estimate_overlap_bias(samplesize_exposure = samplesize_exposure, 
                                            samplesize_outcome = samplesize_outcome, 
                                            n_variants = n_variants, 
                                            rsq_exposure = rsq_exposure, 
                                            overlap_prop = overlap_prop, 
                                            ols_bias = ols_bias, case_prop = 122733/547261)) %>%
  unnest(res)

bias_res %>%
  split_exposure() %>%
  pivot_longer(cols = c(bias, type1_error)) %>%
  ggplot(aes(overlap_prop, value, group = ols_bias, color = as.character(ols_bias))) +
  geom_point() +
  geom_line() +
  facet_grid(rows = vars(exposure), 
             cols = vars(name),
             scales = "free_y") +
  labs(x = "Proportion of Overlapping Participants",
       y = "Value") +
  scale_color_discrete(name = "Bias in \nObservational \nEstimate") +
  theme_bw(base_size = 14) 


