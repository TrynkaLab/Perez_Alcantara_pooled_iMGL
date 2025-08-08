# 8.5. Mendelian randomisation with eQTL results and GWAS

library(tidyverse)
library(TwoSampleMR)
library(MRInstruments)


args = commandArgs(trailingOnly=TRUE)


if (length(args)<6) {
  stop("You need to detail 5 input and 1 output file paths.n",
       call. = FALSE)
} else if (length(args) == 6) {
  sample_info_path = args[1]
  line_prop_changes_path = args[2]
  genotype_path = args[3]
  eqtl_path=args[4]
  output_path=args[5]
  output_csv=args[6]
}


output_path = "../../data/results/8.5.eQTL_MR"
dir.create(output_path,recursive = TRUE)
######### prepare exposure data (eQTL) ########
#https://mrcieu.github.io/TwoSampleMR/articles/exposure.html


exposure_path =  "../../data/results/tensorqtl/35/sum_sizefactorsNorm_log2_scaled_centered_untreated_Not_proliferating_common_500kb_window_tensorQTL_nominal.txt"
exposure = read_tsv(exposure_path) %>%
  dplyr::mutate(effect_allele.exposure = str_split_i(variant_id,pattern = "_",i = 4),
                other_allele.exposure = str_split_i(variant_id,pattern = "_",i = 3))

# convert to MR format
exposure =  format_data(exposure, type = "exposure",
                        beta_col = "slope",
                        se_col = "slope_se",
                        snp_col = "variant_id",
                        eaf_col = "af",
                        effect_allele_col = "effect_allele.exposure",
                        other_allele_col = "other_allele.exposure",
                        gene_col = "phenotype_id",
                        pval_col = "pval_nominal")

head(exposure)

table(exposure$mr_keep.exposure)
exposure$exposure = "eQTL_untreated_Not_proliferating"

exposure$snp_upper = str_to_upper(exposure$SNP)

## need independent SNPs - loading LD clumped from PLINK2
# query gives timeout
indep_snps = read_tsv("../../data/results/?",col_names = FALSE) 

# 
# options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
# 
# clumped = ieugwasr::ld_clump(dplyr::tibble(rsid=exposure$snp_upper[1:100], pval=exposure$pval.exposure[1:100], id=exposure$exposure[1:100]))
# 

######### prepare outcome data (GWAS) ########
gwas_path="/lustre/scratch123/hgi/teams/trynka/resources/summary_statistics/public"
gwas_name = "GCST90012877"
filename_path = list.files(paste0(gwas_path,"/",gwas_name,"/harmonised"),pattern = "*.h.tsv.gz",full.names = TRUE)
if(gwas_name ==  "GCST90012877") { # AD Bellenguez
  outcome_dat = read_outcome_data(
    filename = filename_path,
    sep = "\t",
    snp_col = "hm_variant_id", # careful here
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency",
    pval_col = "p_value"
  )
  
  outcome_dat$outcome = "AD_Bellenguez"
  
}


#### harmonise data against one another #######

dat = harmonise_data(
  exposure_dat = exposure, 
  outcome_dat = outcome_dat
)

# double check SNP id is not losing alleles

##### perform MR #######

res = mr(dat)

res

write_csv(paste0(output_path,"/",unique(res$exposure),"_",unique(res$outcome),".csv"))
