# gather proportion QTL results and visualize

library(patchwork)
library(tidyverse)
library(qvalue)
source("../functions.R")
options(future.globals.maxSize = 20000 * 1024^2) # 20Gb

input_path = "../../../data/results/migration/1.2.proportion_QTL_lm_mean_filtered/"
# input_path = "../../../data/results/migration/1.2.proportion_QTL_mixed_models/"
output_path="../../../data/results/migration/2.check_association_results/lm_1pct_filtered_deflated"
dir.create(file.path(output_path),recursive = TRUE)

pqtl_res_files = list.files(path =input_path, pattern = "prop_QTL_\\d+\\.csv")

# read in csv files from lm() in order
# Extract the numerical order from each file name and sort
order_list = as.integer(sub("prop_QTL_(\\d+)\\.csv", "\\1", pqtl_res_files))
sorted_file_names = pqtl_res_files[order(order_list)]

# Read all files at the same time, otherwise it takes forever

pqtl_res = paste0(input_path,sorted_file_names) %>% 
  readr::read_csv(
  id = "file_path"
)

# omit NAs - regions without eqtl overlaps
pqtl_res = pqtl_res[!is.na(pqtl_res$snp),]
anyDuplicated(pqtl_res$snp)

# remove duplicated snps

# remove rows in which 9 or more columns are NA
pqtl_res = pqtl_res[rowSums(is.na(pqtl_res)) < 9, ]
pqtl_res = pqtl_res %>%
  dplyr::select(-file_path) %>%
  distinct()

# # double check variants
# genotype = readr::read_csv("../../../data/genotypes/full_genotype/genotype_minor_allele_dosage.csv") %>%
#       dplyr::relocate(rn)
# 
# table(genotype$rn %in% pqtl_res$snp) # it's normal that many snps from genotype are not in results
# table(pqtl_res$snp %in% genotype$rn ) # all snps from results should be in genotype 


### checking nominal significance
summary(pqtl_res)
sum(pqtl_res$p_untreated < 5e-8)
sum(pqtl_res$p_LPS < 5e-8)
sum(pqtl_res$p_IFN < 5e-8)

min(pqtl_res$p_untreated) #  2.9e-06
min(pqtl_res$p_IFN) #7.9e-08
min(pqtl_res$p_LPS) # 5.7e-07
table(pqtl_res$p_untreated < 10e-4) #  3434
table(pqtl_res$p_untreated < 10e-5) #  259
table(pqtl_res$p_untreated < 10e-6) #  30

table(pqtl_res$p_IFN < 10e-4) #  4210
table(pqtl_res$p_IFN < 10e-5) #  444
table(pqtl_res$p_IFN < 10e-6) #  97
table(pqtl_res$p_IFN < 10e-7) #  10
# multiple correction adjustment ####

# Bonferroni, the usual correction in GWAS
pqtl_res$p_bonf_untreated = p.adjust(pqtl_res$p_untreated, method="bonferroni")
pqtl_res$p_bonf_IFN = p.adjust(pqtl_res$p_IFN, method="bonferroni")
pqtl_res$p_bonf_LPS = p.adjust(pqtl_res$p_LPS, method="bonferroni")

pqtl_res$p_BH_untreated = p.adjust(pqtl_res$p_untreated, method="BH")
pqtl_res$p_BH_IFN = p.adjust(pqtl_res$p_IFN, method="BH")
pqtl_res$p_BH_LPS = p.adjust(pqtl_res$p_LPS, method="BH")

nrow(pqtl_res) #  3,842,190
nrow(pqtl_res) == length(unique(pqtl_res$snp))

# 
# # raw p-values distribution - untreated
p1 = ggplot(pqtl_res,aes(x = p_untreated)) +
   geom_histogram(bins = 100) +
   theme_bw()
pdf(paste0(output_path,"GWAS_untreated_pvals_distrib.pdf"),
    width = 7, height = 4 )
plot(p1)
dev.off()
# 
 summary(pqtl_res$p_untreated) #median of P-values should be 0.5 under the null


## q-q plot using chi squared statistics under null and in our data
#Under NULL p-values are Uniformly distributed between 0 and 1,
#hence chisq-stats are expected to be:
expect.stats = qchisq(ppoints(nrow(pqtl_res)), df = 1, lower = F)
obs.stats = qchisq(pqtl_res$p_untreated, df = 1, lower = F)
lambda = median(na.omit(obs.stats)) / median(expect.stats) #GC lambda = ratio at medians

# subsample to variants not in LD

not_LD = read_tsv("../../../../OTAR2065_sc_eQTL/data/genotype/plink_genotypes/all_pools.all_donors.genotype.MAF05.LDpruned.prune.in",
                  col_names = FALSE) %>%
  dplyr::mutate(X1 = str_remove(X1,"chr"))

positions = which(pqtl_res$snp %in% not_LD$X1)

x = expect.stats[positions]
y = obs.stats[positions]

pdf(paste0(output_path,"QQ_plot_migration_GWAS_untreated_pvals.pdf"),
       width = 13, height = 10 )
qqplot(x, y, xlab = "chisq expected", ylab = "chisq observed",
       sub = substitute(paste(lambda, "=", lam), list(lam = signif(lambda,3))), 
       cex = 0.8, col = "red")
abline(0,1)
dev.off()

### Inflated including reps and pools, deflated when averaging across pools
obs.stats = qchisq(pqtl_res$p_IFN, df = 1, lower = F)
lambda = median(na.omit(obs.stats)) / median(expect.stats) #GC lambda = ratio at medians
x = expect.stats[positions]
y = obs.stats[positions]

pdf(paste0(output_path,"QQ_plot_migration_GWAS_IFN_pvals.pdf"),
    width = 13, height = 10 )
qqplot(x, y, xlab = "chisq expected", ylab = "chisq observed",
       sub = substitute(paste(lambda, "=", lam), list(lam = signif(lambda,3))), 
       cex = 0.8, col = "red")
abline(0,1)

dev.off()

obs.stats = qchisq(pqtl_res$p_LPS, df = 1, lower = F)
lambda = median(na.omit(obs.stats)) / median(expect.stats) #GC lambda = ratio at medians
x = expect.stats[positions]
y = obs.stats[positions]

pdf(paste0(output_path,"QQ_plot_migration_GWAS_LPS_pvals.pdf"),
    width = 13, height = 10 )
qqplot(x, y, xlab = "chisq expected", ylab = "chisq observed",
       sub = substitute(paste(lambda, "=", lam), list(lam = signif(lambda,3))), 
       cex = 0.8, col = "red")
abline(0,1)

dev.off()

# Deflated QQ plot, see here for an explanation: 
# https://stats.stackexchange.com/questions/185284/deflated-qq-plots-in-genome-wide-association-studies
# number of significant results

# BH

min(pqtl_res$p_BH_untreated ) # 0.95
min(pqtl_res$p_BH_IFN ) # 0.17

# Bonferroni
min(pqtl_res$p_bonf_untreated ) # 1
min(pqtl_res$p_bonf_IFN ) # 0.30

# write all results for coloc
pqtl_res %>%
   readr::write_csv(.,"../../../data/results/migration/2.check_association_results/all_res_pqtl_migration.csv")

# manhattan plot
test = pqtl_res %>%
   dplyr::select(snp:p_untreated) %>%
   dplyr::mutate(CHR=as.numeric(stringr::str_split_fixed(string = snp,pattern = "_",n = 4)[,1]),
                 BP=as.numeric(stringr::str_split_fixed(string = snp,pattern = "_",n = 4)[,2])) %>%
   dplyr::rename(P=p_untreated,SNP=snp) %>%
  dplyr::mutate(row_id = row_number()) %>% # Add a row id for tracking
  dplyr::group_by(P >0.05) %>%               # Group by the condition
  dplyr::mutate(
    keep_row = if_else(P >0.05, row_number() %% 20 == 1, TRUE) # Keep every 20th row if P >0.05
  ) %>%
  dplyr::filter(keep_row) %>%              # Filter rows to keep
  dplyr::ungroup() %>%
  dplyr::select(SNP,P,CHR,BP)                


myManhattan( test, graph.title = "migration GWAS: untreated", 
             suggestiveline = 1e-5,
             suggestivecolor="lightgrey",
             col = c("grey30","lightgrey"),
                 genomewideline = 5e-8, 
             even.facet = T,highlight = 1e-5,
             highlight.col = "orange2") 
ggsave(paste0(output_path,"migration_GWAS_manhattan_untreated.png"),
    width = 12, height = 7 )

ggsave(paste0(output_path,"migration_GWAS_manhattan_untreated.svg"),
       width = 12, height = 7 )
# IFN, LPS
test = pqtl_res %>%
   dplyr::select(snp:p_IFN) %>%
   dplyr::mutate(CHR=as.numeric(stringr::str_split_fixed(string = snp,pattern = "_",n = 4)[,1]),
                 BP=as.numeric(stringr::str_split_fixed(string = snp,pattern = "_",n = 4)[,2])) %>%
   dplyr::rename(P=p_IFN,SNP=snp)  %>%
  dplyr::mutate(row_id = row_number()) %>% # Add a row id for tracking
  dplyr::group_by(P >0.05) %>%               # Group by the condition
  dplyr::mutate(
    keep_row = if_else(P >0.05, row_number() %% 20 == 1, TRUE) # Keep every 20th row if P >0.05
  ) %>%
  dplyr::filter(keep_row) %>%              # Filter rows to keep
  dplyr::ungroup() %>%
  dplyr::select(SNP,P,CHR,BP)                


myManhattan( test, graph.title = "migration GWAS: IFN", 
             suggestiveline = 1e-5,
             suggestivecolor="lightgrey",
             col = c("grey30","lightgrey"),
             genomewideline = 5e-8, 
             even.facet = T,highlight = 1e-5,
             highlight.col = "orange2") 
ggsave(paste0(output_path,"migration_GWAS_manhattan_IFN.png"),
       width = 12, height = 7 )
ggsave(paste0(output_path,"migration_GWAS_manhattan_IFN.svg"),
       width = 12, height = 7 )

test = pqtl_res %>%
   dplyr::select(snp:p_LPS) %>%
   dplyr::mutate(CHR=as.numeric(stringr::str_split_fixed(string = snp,pattern = "_",n = 4)[,1]),
                 BP=as.numeric(stringr::str_split_fixed(string = snp,pattern = "_",n = 4)[,2])) %>%
   dplyr::rename(P=p_LPS,SNP=snp)  %>%
  dplyr::mutate(row_id = row_number()) %>% # Add a row id for tracking
  dplyr::group_by(P >0.05) %>%               # Group by the condition
  dplyr::mutate(
    keep_row = if_else(P >0.05, row_number() %% 20 == 1, TRUE) # Keep every 20th row if P >0.05
  ) %>%
  dplyr::filter(keep_row) %>%              # Filter rows to keep
  dplyr::ungroup() %>%
  dplyr::select(SNP,P,CHR,BP)                


myManhattan( test, graph.title = "migration GWAS: LPS", 
             suggestiveline = 1e-5,
             suggestivecolor="lightgrey",
             col = c("grey30","lightgrey"),
             genomewideline = 5e-8, 
             even.facet = T,highlight = 1e-5,
             highlight.col = "orange2")
ggsave(paste0(output_path,"migration_GWAS_manhattan_LPS.png"),
       width = 12, height = 7 )
ggsave(paste0(output_path,"migration_GWAS_manhattan_LPS.svg"),
       width = 12, height = 7 )


###########
# filter to (nominally) significant
pqtl_res %>%
  dplyr::filter(p_untreated < 10e-5 | p_IFN < 10e-5 | p_LPS < 10e-5 ) %>%
  write_csv(paste0(output_path,"nominally_significant_results.csv"))
