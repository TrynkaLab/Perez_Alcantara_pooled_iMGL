# 8.5. Mendelian randomisation with eQTL results and GWAS

library(tidyverse)
library(ggbreak)
library(ggrepel)
library(mrSampleOverlap)
library(TwoSampleMR)
source("./functions.R")
options(scipen = 999) # prevents from showing very small numbers as 0


output_path = "../../data/results/8.6.MR_eQTL_pheno_bias_correction"
dir.create(output_path,recursive = TRUE)

#### functions #####

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

###

### simple version after farmaggedon ####
TWMR_p = readr::read_csv("../../data/results/8.5.eQTL_MR/TWMR/output/phagocytosis_migration_twmr_p.csv") 

lead_variants_tensorqtl = read_csv("../../data/for_tensorQTL/tensorQTL_eGene_topSNP.csv") %>%
  filter(gene_name %in% TWMR_p$gene & group %in% c("73_IFN_Not_proliferating" , "63_untreated_Not_proliferating" ,"84_LPS_Not_proliferating")) %>%
  rename(gene = gene_name, pval.exposure = pval_nominal) %>%
  mutate(samplesize.exposure = if_else(treatment == "untreated",189,188))
# some twmr results use more than 1 variant - but if the lead one is not strong enough, the second variant used won't be any good

TWMR_p = TWMR_p %>%
  left_join(lead_variants_tensorqtl %>% select(gene,treatment,pval.exposure,samplesize.exposure))

# # A more liberal estimation per IV than the one used in mrSampleOverlap:
TWMR_p = TWMR_p %>%
mutate(F_statistic = pmap_dbl(
  list(pval = pval.exposure, samplesize = samplesize.exposure),
  calculate_F_statistic_MR
))
write_csv(TWMR_p,"../../data/results/8.5.eQTL_MR/TWMR/output/phagocytosis_migration_twmr_p_with_F.csv")

# retain those with F stat over 10
TWMR_p = TWMR_p %>%
  filter(F_statistic > 10)
# 32 unique genes with F stat over 10 & p_Bonf < 0.05

### plot manhattan plot of TWMR gene associations
genes_to_location = readr::read_csv("../../../resources/biomart/Homo_sapiens.GRCh38.111.genes.csv") %>%
  dplyr::rename(CHR = seqname,BP = start, gene = gene_name) %>%
  dplyr::right_join(TWMR_p) %>%
  tidyr::drop_na() %>%
  dplyr::rename(P = p) %>% # selecting nominal p-value to plot
  dplyr::select(gene,CHR,BP,P, treatment, phenotype,p_Bonf) %>%
  dplyr::mutate(CHR = as.numeric(CHR),
                BP = as.numeric(BP)) %>%
  dplyr::mutate(logP_TWMR = -log10(P)) %>%
  dplyr::select(CHR,BP,logP_TWMR,gene,treatment, phenotype,p_Bonf) 



toplot_unt = genes_to_location %>% 
  filter(treatment == "untreated") 

p_Bonf_threshold = toplot_unt %>%
  filter(p_Bonf < 0.05) %>%
  slice_min(logP_TWMR) %>%
  pull(logP_TWMR)

toplot_unt = toplot_unt %>%
  dplyr::arrange(CHR,BP) %>%
  dplyr::mutate(new_pos = 1:nrow(.),
                
                col_TWMR = case_when(CHR %% 2 == 1 ~ "grey70",
                                     CHR %% 2 != 1 ~ "grey20"
                )) %>%
  dplyr::mutate(
    col_TWMR = case_when(logP_TWMR >p_Bonf_threshold~ "darkred",
                         .default = col_TWMR)) %>%
  dplyr::group_by(CHR) %>%
  dplyr::mutate(midpoints= mean(new_pos)) %>% # creating midpoints for label placement later
  dplyr::ungroup()

# untreated
p = toplot_unt%>%
  ggplot() +
  geom_hline(yintercept = p_Bonf_threshold, linetype = "dotted", col = "darkred") +
  geom_point(aes(x = interaction(new_pos, CHR), y = logP_TWMR, 
                 color = col_TWMR,
                 shape = phenotype),
             ) +
  geom_text_repel(
    aes(x = new_pos, y = logP_TWMR, label = ifelse(col_TWMR == "darkred", gene, NA)),
    size = 3, color = "grey10",
    vjust = -0.5, hjust = 0.5,
  ) +
  scale_color_identity() +
  
  ## Add y-axis breaks with clear tick marks
  scale_y_break(c(15, 63), scales = 0.1, ticklabels = c(63,65) ) +
  scale_y_break(c(65, 151), scales = 0.1, ticklabels = c(151,153)) +
  
  labs(
    x = "Chromosome",
    y = expression(-log[10]~p-value),
  ) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_blank(),
    plot.margin = unit(c(1, 1, 3, 1), units = "lines")
  ) +
  coord_cartesian(clip = "off") + 
  annotate(geom = "text", x = unique(toplot_unt$midpoints),
           y = -0.5, angle = 45,
           label = unique(toplot_unt$CHR), size = 3) +
  geom_hline(yintercept = 0.0001, linetype = "solid", col = "grey10", linewidth = 1) +
  geom_segment(
    aes(x = 0.4, xend = 0.4, y = 0, yend = 153),
    linetype = "solid", color = "black", linewidth = 1
  ) 

pdf(file = "../../data/results/8.5.eQTL_MR/TWMR/output/TWMR_plot_untreated.pdf", width = 10, height=5)
p
dev.off()

#IFN
toplot_ifn = genes_to_location %>% 
  filter(treatment == "IFN") 

p_Bonf_threshold = toplot_ifn %>%
  filter(p_Bonf < 0.05) %>%
  slice_min(logP_TWMR) %>%
  pull(logP_TWMR)

toplot_ifn = toplot_ifn %>%
  dplyr::arrange(CHR,BP) %>%
  dplyr::mutate(new_pos = 1:nrow(.),
                
                col_TWMR = case_when(CHR %% 2 == 1 ~ "grey70",
                                     CHR %% 2 != 1 ~ "grey20"
                )) %>%
  dplyr::mutate(
    col_TWMR = case_when(logP_TWMR > p_Bonf_threshold~ "darkred",
                         .default = col_TWMR)) %>%
  dplyr::group_by(CHR) %>%
  dplyr::mutate(midpoints= mean(new_pos)) %>% # creating midpoints for label placement later
  dplyr::ungroup()

p = toplot_ifn%>%
  ggplot() +
  geom_hline(yintercept = p_Bonf_threshold, linetype = "dotted", col = "darkred") +
  geom_point(aes(x = interaction(new_pos, CHR), y = logP_TWMR, 
                 color = col_TWMR,
                 shape = phenotype),
  ) +
  geom_text_repel(
    aes(x = new_pos, y = logP_TWMR, label = ifelse(col_TWMR == "darkred", gene, NA)),
    size = 3, color = "grey10",
    vjust = -0.5, hjust = 0.5,
  ) +
  scale_color_identity() +
  
  labs(
    x = "Chromosome",
    y = expression(-log[10]~p-value),
  ) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_blank(),
    plot.margin = unit(c(1, 1, 3, 1), units = "lines")
  ) +
  coord_cartesian(clip = "off") + 
  annotate(geom = "text", x = unique(toplot_ifn$midpoints),
           y = -0.5, angle = 45,
           label = unique(toplot_ifn$CHR), size = 3) +
  geom_hline(yintercept = 0.0001, linetype = "solid", col = "grey10", linewidth = 1) +
  geom_segment(
    aes(x = 0.4, xend = 0.4, y = 0, yend = 18),
    linetype = "solid", color = "black", linewidth = 1
  ) 

pdf(file = "../../data/results/8.5.eQTL_MR/TWMR/output/TWMR_plot_IFN.pdf", width = 10, height=5)
p
dev.off()

## LPS
toplot_lps = genes_to_location %>% 
  filter(treatment == "LPS") 

p_Bonf_threshold = toplot_lps %>%
  filter(p_Bonf < 0.05) %>%
  slice_min(logP_TWMR) %>%
  pull(logP_TWMR)

toplot_lps = toplot_lps %>%
  dplyr::arrange(CHR,BP) %>%
  dplyr::mutate(new_pos = 1:nrow(.),
                
                col_TWMR = case_when(CHR %% 2 == 1 ~ "grey70",
                                     CHR %% 2 != 1 ~ "grey20"
                )) %>%
  dplyr::mutate(
    col_TWMR = case_when(logP_TWMR > p_Bonf_threshold~ "darkred",
                         .default = col_TWMR)) %>%
  dplyr::group_by(CHR) %>%
  dplyr::mutate(midpoints= mean(new_pos)) %>% # creating midpoints for label placement later
  dplyr::ungroup()


p = toplot_lps%>%
  ggplot() +
  geom_hline(yintercept = p_Bonf_threshold, linetype = "dotted", col = "darkred") +
  geom_point(aes(x = interaction(new_pos, CHR), y = logP_TWMR, 
                 color = col_TWMR,
                 shape = phenotype),
  ) +
  geom_text_repel(
    aes(x = new_pos, y = logP_TWMR, label = ifelse(col_TWMR == "darkred", gene, NA)),
    size = 3, color = "grey10",
    vjust = -0.5, hjust = 0.5,
  ) +
  scale_color_identity() +
  
  labs(
    x = "Chromosome",
    y = expression(-log[10]~p-value),
  ) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_blank(),
    plot.margin = unit(c(1, 1, 3, 1), units = "lines")
  ) +
  coord_cartesian(clip = "off") + 
  annotate(geom = "text", x = unique(toplot_lps$midpoints),
           y = -0.5, angle = 45,
           label = unique(toplot_lps$CHR), size = 3) +
  geom_hline(yintercept = 0.0001, linetype = "solid", col = "grey10", linewidth = 1) +
  geom_segment(
    aes(x = 0.4, xend = 0.4, y = 0, yend = 18),
    linetype = "solid", color = "black", linewidth = 1
  ) 

pdf(file = "../../data/results/8.5.eQTL_MR/TWMR/output/TWMR_plot_LPS.pdf", width = 10, height=5)
p
dev.off()
#### previous version with more information

# exposures
nominal_res = list()
nominal_res = lapply(c("untreated","IFN","LPS"), load_nominal_results)
names(nominal_res) = c("untreated","IFN","LPS")

# outcomes

pheno_gwas = lapply(c("phagocytosis","migration"),function(x) {
  readr::read_csv(
  file = paste0(
    "../../../OTAR2065_phenotypic_QTLs/data/results/",
    x,
    "/2.check_association_results/lm_1pct_filtered_deflated/all_res_pqtl_",
    x,
    ".csv"
  )
)
  }
)
names(pheno_gwas) = c("phagocytosis","migration")


# TWMR tables
TWMR_p = readr::read_csv("../../data/results/8.5.eQTL_MR/TWMR/output/phagocytosis_migration_twmr_p.csv")
TWMR_res = readRDS("../../data/results/8.5.eQTL_MR/TWMR/output/phagocytosis_migration_twmr_res.rds")

# for all significant results
sign_res = TWMR_p %>%
  dplyr::filter(p_Bonf < 0.05)

E_F = list() # F parameter
for(n in 1:nrow(sign_res)){
  

  focus_gene = sign_res[n,] %>%
    pull(gene)
  treatment = sign_res[n,] %>%
    pull(treatment)
  phenotype = sign_res[n,] %>%
    pull(phenotype)

# TWMR_res[["LPS_LRRK2"]]
# # only one gene: itself
# treatment = "LPS"
# focus_gene = "LRRK2"
# phenotype = "phagocytosis"

# get variants used for TWMR calculations

  if(nrow(TWMR_res[[paste0(phenotype,"_",treatment,"_",focus_gene)]]$LD_matrix) == 1){
    snps = list()
    snps$SNP = rownames(TWMR_res[[paste0(phenotype,"_",treatment,"_",focus_gene)]]$LD_matrix)
    snps$A2 = stringr::str_split_i(snps$SNP,pattern = "_",3)
    snps$A1 = stringr::str_split_i(snps$SNP,pattern = "_",4)
  }else{
    
  
snps = data.frame("SNP"=names(TWMR_res[[paste0(phenotype,"_",treatment,"_",focus_gene)]]$LD_matrix),
                  "A1"="A1","A2"="A2")
snps$A2 = stringr::str_split_i(snps$SNP,pattern = "_",3)
snps$A1 = stringr::str_split_i(snps$SNP,pattern = "_",4)
}
# nominal results already tests for all variants 250kb around the gene
exposure = nominal_res[[treatment]] %>%
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


outcome = pheno_gwas[[phenotype]] %>%
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
                        ols_bias = seq(-7, 7, 0.5))


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
E_F[[paste0(focus_gene,"_",treatment,"_",phenotype)]]$F_stat_estimate = mrSampleOverlap::estimate_f(samplesize_exposure = dat_summarized$samplesize_exposure,
                                  n_variants = dat_summarized$n_variants,
                                  rsq_exposure = dat_summarized$rsq_exposure,
                                  lci_95 = FALSE)
E_F[[paste0(focus_gene,"_",treatment,"_",phenotype)]]$F_lower_limit_95CI = mrSampleOverlap::estimate_f(samplesize_exposure = dat_summarized$samplesize_exposure,
                                                                                                    n_variants = dat_summarized$n_variants,
                                                                                                    rsq_exposure = dat_summarized$rsq_exposure,
                                                                                                    lci_95 = TRUE)# returning a more conservative estimate
# # A more liberal estimation per IV:
# calculate_F_statistic_MR(pval = exposure$pval.exposure[1],samplesize = exposure$samplesize.exposure[1])
# # over 10 is considered strong
# 
# # We know our OLS:the ordinary least squares (OLS, also known as standard least squares regression) 
# # estimate is obtained by regressing the outcome on the risk factor
 fit = lm(beta.outcome ~ beta.exposure,data = dat_harmonized)
 OLS_estimate = fit$coefficients["beta.exposure"]
# # is the OLS esimate the OLS bias? 
# bias = unname(OLS_estimate_bias * 1 * (1/E_F))

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

### this only accounts for weak instrument bias and not for winner's curse
# while MRLap accounts for both

}
### F parameter
F_parameter = tibble(
  name = names(E_F),
  F_stat_estimate = map_dbl(E_F, 1),
  F_lower_limit_95CI = map_dbl(E_F, 2)  # use map_chr or map_dbl depending on type
) %>%
  tidyr::separate(name,into = c("gene","treatment","phenotype"),sep = "_")

readr::write_csv(F_parameter,paste0(outputDir,"F_parameter_sign_twmr_genes.csv"))

TWMR_p %>%
  dplyr::filter(p_Bonf < 0.05) %>%
  dplyr::left_join(F_parameter) %>%
  readr::write_csv(paste0(outputDir,"TWMR_p_sign_F_parameter.csv"))


####
# ### example data of this bias correction method

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