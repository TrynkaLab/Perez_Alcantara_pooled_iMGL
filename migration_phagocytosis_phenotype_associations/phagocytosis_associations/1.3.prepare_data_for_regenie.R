# preparing data for regenie


library(patchwork)
library(tidyverse)
args = commandArgs(trailingOnly=TRUE)

message(length(args)," arguments provided")
if (length(args)<6) {
  stop("You need to detail 3 input and 2 output file paths, and a blacklist of lines.n",
       call. = FALSE)
} else if (length(args) == 6) {
  sample_info_path = args[1]
  donor_prop_changes_path = args[2]
  genotype_pc_path = args[3]
  output_pheno_path=args[4]
  output_cov_path=args[5]
  donor_blacklist=args[6]
}



# to test
# output_dir = "../../data/results/1.3.data_for_regenie/"
# dir.create(file.path(output_dir),recursive = T)
# output_pheno_path = paste0(output_dir,"phenotype_for_regenie.txt")
# output_cov_path = paste0(output_dir,"covariates_for_regenie.txt")
# sample_info = read.csv("../../../data/all_pools_phagocytosis_sample_info.csv")
# donor_prop_changes = read.csv("../../../data/results/phagocytosis/1.check_line_proportions/line_prop_changes_averages_variances_1pct_filtered.csv")
# full_genotype_pcs = read.table("../../../../OTAR2065_sc_eQTL/data/genotype/plink_genotypes/all_pools.no_outliers.genotype.MAF05.eigenvec")
# donor_blacklist="letw_5,lizq_3,zaie_1,romx_2,sebn_4,seru_7,qonc_2,boqx_2,garx_2,sojd_3,yoch_6,sh5y5y"

donor_blacklist = str_split(donor_blacklist,pattern = ",")[[1]]
full_genotype_pcs = read.table(genotype_pc_path)
rownames(full_genotype_pcs) = ifelse(full_genotype_pcs$V1 ==full_genotype_pcs$V2,
                                     yes = full_genotype_pcs$V1,
                                     no = full_genotype_pcs$V2)
full_genotype_pcs = full_genotype_pcs[c(2:22)]
colnames(full_genotype_pcs) = c("IID",paste0("genotypePC",1:20))


sample_info = read.csv(sample_info_path)
donor_prop_changes = read.csv(donor_prop_changes_path)
length(unique(donor_prop_changes$line))
# load minor allele dosage file for all donors
# 

# load significant eQTL results and subset to those there
donor_prop_changes = donor_prop_changes %>%
  dplyr::select(mean_outcome_pool,line,treatment,n_pools,sex)

length(unique(donor_prop_changes$line)) # 143 lines

# creating tab-separated phenotype file
# FID IID Untreated IFN LPS

pheno_file = donor_prop_changes %>%
  dplyr::mutate(FID=stringr::str_split_fixed(line,pattern = "_",n=2)[,1]) %>%
  dplyr::rename(IID = line) %>%
  dplyr::select(FID,IID,mean_outcome_pool, treatment) %>%
  dplyr::filter(!IID %in% donor_blacklist) %>%
  distinct() %>%
  tidyr::pivot_wider(names_from = treatment, values_from = mean_outcome_pool)

pheno_file %>%
  write.table(.,file = output_pheno_path, col.names = TRUE,row.names = FALSE,sep = "\t",quote = FALSE)
# covariate file
# FID IID C1 C2 C3

cov_file = donor_prop_changes %>%
  dplyr::mutate(FID=stringr::str_split_fixed(line,pattern = "_",n=2)[,1]) %>%
  dplyr::rename(IID = line) %>%
  dplyr::select(FID,IID, sex,n_pools) %>%
  distinct() %>%
  dplyr::left_join(.,full_genotype_pcs) %>%
  dplyr::filter(!IID %in% donor_blacklist) %>%
  dplyr::filter(IID %in% pheno_file$IID) %>%
  dplyr::group_by(IID) %>%
  dplyr::mutate(n_pools= max(n_pools)) %>% # fixing cases where number of pool replicates is different between treatments
  distinct() 
  

table(pheno_file$IID %in% cov_file$IID)
table(cov_file$IID %in% pheno_file$IID)
setdiff(pheno_file$IID, cov_file$IID)
setdiff(cov_file$IID,pheno_file$IID)
nrow(pheno_file)
nrow(cov_file)

cov_file %>%
  write.table(.,file = output_cov_path, col.names = TRUE,row.names = FALSE,sep = "\t",quote = FALSE)



