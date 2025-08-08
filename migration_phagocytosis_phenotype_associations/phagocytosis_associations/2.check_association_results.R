# gather proportion QTL results and visualize

library(patchwork)
library(tidyverse)
library(qvalue)
source("../functions.R")
options(future.globals.maxSize = 40000 * 1024^2) # 40Gb

input_path = "../../../data/results/phagocytosis/1.2.proportion_QTL_lm_mean_filtered/"
# input_path = "../../../data/results/phagocytosis/1.2.proportion_QTL_mixed_models/"
output_path="../../../data/results/phagocytosis/2.check_association_results/lm_1pct_filtered_deflated/"
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
sum(pqtl_res$p_untreated < 5e-8) # 0
sum(pqtl_res$p_LPS < 5e-8) # 0
sum(pqtl_res$p_IFN < 5e-8) # 0

min(pqtl_res$p_untreated) #  1.5e-06
min(pqtl_res$p_IFN) #1.1e-06
min(pqtl_res$p_LPS) # 1.9e-06
table(pqtl_res$p_untreated < 10e-4) #  3695
table(pqtl_res$p_untreated < 10e-5) #  514
table(pqtl_res$p_untreated < 10e-6) #  45
# multiple correction adjustment ####

# Bonferroni, the usual correction in GWAS
pqtl_res$p_bonf_untreated = p.adjust(pqtl_res$p_untreated, method="bonferroni")
pqtl_res$p_bonf_IFN = p.adjust(pqtl_res$p_IFN, method="bonferroni")
pqtl_res$p_bonf_LPS = p.adjust(pqtl_res$p_LPS, method="bonferroni")

pqtl_res$p_BH_untreated = p.adjust(pqtl_res$p_untreated, method="BH")
pqtl_res$p_BH_IFN = p.adjust(pqtl_res$p_IFN, method="BH")
pqtl_res$p_BH_LPS = p.adjust(pqtl_res$p_LPS, method="BH")

nrow(pqtl_res) # 3,972,431
nrow(pqtl_res) == length(unique(pqtl_res$snp)) # TRUE

# 
# # raw p-values distribution - untreated
p1 = ggplot(pqtl_res,aes(x = p_untreated)) +
   geom_histogram(bins = 100) +
   theme_bw()
p2 = ggplot(pqtl_res,aes(x = p_IFN)) +
  geom_histogram(bins = 100) +
  theme_bw()
p3 = ggplot(pqtl_res,aes(x = p_LPS)) +
  geom_histogram(bins = 100) +
  theme_bw()
pdf(paste0(output_path,"GWAS_pvals_distrib.pdf"),
    width = 7, height = 4 )
plot(p1 + p2 + p3)
dev.off()
# 
 summary(pqtl_res$p_untreated) #median of P-values should be 0.5 under the null, and it is


## q-q plot using chi squared statistics under null and in our data


# subsample to variants not in LD

not_LD = read_tsv("../../../../OTAR2065_sc_eQTL/data/genotype/plink_genotypes/all_pools.all_donors.genotype.MAF05.LDpruned.prune.in",
                  col_names = FALSE) %>%
  dplyr::mutate(X1 = str_remove(X1,"chr"))

positions = which(pqtl_res$snp %in% not_LD$X1)

build_qqplot_stats = function(pvals, positions){
  #Under NULL p-values are Uniformly distributed between 0 and 1,
  #hence chisq-stats are expected to be:
  expect.stats = qchisq(ppoints(length(pvals)), df = 1, lower = F)
  obs.stats = qchisq(pvals, df = 1, lower = F)
  lambda = median(na.omit(obs.stats)) / median(expect.stats) #GC lambda = ratio at medians
  
  df = data.frame(x = expect.stats[positions], y = obs.stats[positions], lambda = lambda)
  return(df)
}

untreated_qq = build_qqplot_stats(pqtl_res$p_untreated, positions)
IFN_qq = build_qqplot_stats(pqtl_res$p_IFN, positions)
LPS_qq = build_qqplot_stats(pqtl_res$p_LPS, positions)

pdf(paste0(output_path,"QQ_plot_phagocytosis_GWAS_pvals_untreated.pdf"),
       width = 7, height = 7 )
qqplot(untreated_qq$x, untreated_qq$y, xlab = "chisq expected", ylab = "chisq observed",
       sub = substitute(paste(lambda, "=", lam), list(lam = signif(untreated_qq$lambda,3))), 
       cex = 0.8, col = "red")
abline(0,1)
dev.off()

pdf(paste0(output_path,"QQ_plot_phagocytosis_GWAS_pvals_IFN.pdf"),
    width = 7, height = 7 )
qqplot(IFN_qq$x, IFN_qq$y, xlab = "chisq expected", ylab = "chisq observed",
       sub = substitute(paste(lambda, "=", lam), list(lam = signif(IFN_qq$lambda,3))), 
       cex = 0.8, col = "red")
abline(0,1)
dev.off()

pdf(paste0(output_path,"QQ_plot_phagocytosis_GWAS_pvals_LPS.pdf"),
    width = 7, height = 7 )
qqplot(LPS_qq$x, LPS_qq$y, xlab = "chisq expected", ylab = "chisq observed",
       sub = substitute(paste(lambda, "=", lam), list(lam = signif(LPS_qq$lambda,3))), 
       cex = 0.8, col = "red")
abline(0,1)
dev.off()

### Inflated including reps and pools, deflated when averaging across pools
# Deflated QQ plot, see here for an explanation: 
# https://stats.stackexchange.com/questions/185284/deflated-qq-plots-in-genome-wide-association-studies
# number of significant results

# BH

min(pqtl_res$p_BH_untreated ) # 0.36

# Bonferroni
min(pqtl_res$p_bonf_untreated ) # 1

# write all results for coloc
pqtl_res %>%
   readr::write_csv(.,"../../../data/results/phagocytosis/2.check_association_results/all_res_pqtl_phagocytosis.csv")

# manhattan plot
manhattan = pqtl_res %>%
   dplyr::select(snp:p_untreated) %>%
   dplyr::mutate(CHR=as.numeric(stringr::str_split_fixed(string = snp,pattern = "_",n = 4)[,1]),
                 BP=as.numeric(stringr::str_split_fixed(string = snp,pattern = "_",n = 4)[,2])) %>%
   dplyr::rename(P=p_untreated,SNP=snp) %>%
  dplyr::mutate(row_id = row_number()) %>% # Add a row id for tracking
  dplyr::group_by(P >0.05) %>%               # Group by the condition
  dplyr::mutate(
    keep_row = if_else(P >0.30, row_number() %% 30 == 1, # Keep every 30th row if P >0.30
                       if_else(P >0.05, row_number() %% 10 == 1,  # Keep every 10th row if P >0.05
                               if_else(P >0.01, row_number() %% 2 == 1, TRUE))) # Keep every other row if P >0.01
  ) %>%
  dplyr::filter(keep_row) %>%              # Filter rows to keep
  dplyr::ungroup() %>%
  dplyr::select(SNP,P,CHR,BP)                


myManhattan( manhattan, graph.title = "Phagocytosis GWAS: untreated", 
             suggestiveline = 1e-5,
             suggestivecolor="lightgrey",
             col = c("grey30","lightgrey"),
                 genomewideline = 5e-8, 
             even.facet = T,highlight = 1e-5,
             highlight.col = "orange2",
             genomewidecolor = "darkred") 
ggsave(paste0(output_path,"phagocytosis_GWAS_manhattan_untreated.png"),
    width = 12, height = 7 )

ggsave(paste0(output_path,"phagocytosis_GWAS_manhattan_untreated.svg"),
       width = 12, height = 7 )
# IFN, LPS
manhattan = pqtl_res %>%
   dplyr::select(snp:p_IFN) %>%
   dplyr::mutate(CHR=as.numeric(stringr::str_split_fixed(string = snp,pattern = "_",n = 4)[,1]),
                 BP=as.numeric(stringr::str_split_fixed(string = snp,pattern = "_",n = 4)[,2])) %>%
   dplyr::rename(P=p_IFN,SNP=snp)  %>%
  dplyr::mutate(row_id = row_number()) %>% # Add a row id for tracking
  dplyr::group_by(P >0.05) %>%               # Group by the condition
  dplyr::mutate(
    keep_row = if_else(P >0.30, row_number() %% 30 == 1, # Keep every 30th row if P >0.30
                       if_else(P >0.05, row_number() %% 10 == 1,  # Keep every 10th row if P >0.05
                               if_else(P >0.01, row_number() %% 2 == 1, TRUE))) # Keep every other row if P >0.01
  ) %>%
  dplyr::filter(keep_row) %>%              # Filter rows to keep
  dplyr::ungroup() %>%
  dplyr::select(SNP,P,CHR,BP)                


myManhattan( manhattan, graph.title = "Phagocytosis GWAS: IFN", 
             suggestiveline = 1e-5,
             suggestivecolor="lightgrey",
             col = c("grey30","lightgrey"),
             genomewideline = 5e-8, 
             even.facet = T,highlight = 1e-5,
             highlight.col = "orange2",
             genomewidecolor = "darkred") 
ggsave(paste0(output_path,"phagocytosis_GWAS_manhattan_IFN.png"),
       width = 12, height = 7 )
ggsave(paste0(output_path,"phagocytosis_GWAS_manhattan_IFN.svg"),
       width = 12, height = 7 )

manhattan = pqtl_res %>%
   dplyr::select(snp:p_LPS) %>%
   dplyr::mutate(CHR=as.numeric(stringr::str_split_fixed(string = snp,pattern = "_",n = 4)[,1]),
                 BP=as.numeric(stringr::str_split_fixed(string = snp,pattern = "_",n = 4)[,2])) %>%
   dplyr::rename(P=p_LPS,SNP=snp)  %>%
  dplyr::mutate(row_id = row_number()) %>% # Add a row id for tracking
  dplyr::group_by(P >0.05) %>%               # Group by the condition
  dplyr::mutate(
    keep_row = if_else(P >0.30, row_number() %% 30 == 1, # Keep every 30th row if P >0.30
                       if_else(P >0.05, row_number() %% 10 == 1,  # Keep every 10th row if P >0.05
                               if_else(P >0.01, row_number() %% 2 == 1, TRUE))) # Keep every other row if P >0.01
  ) %>%
  dplyr::filter(keep_row) %>%              # Filter rows to keep
  dplyr::ungroup() %>%
  dplyr::select(SNP,P,CHR,BP)                


myManhattan( manhattan, graph.title = "Phagocytosis GWAS: LPS", 
             suggestiveline = 1e-5,
             suggestivecolor="lightgrey",
             col = c("grey30","lightgrey"),
             genomewideline = 5e-8, 
             even.facet = T,highlight = 1e-5,
             highlight.col = "orange2",
             genomewidecolor = "darkred") 
ggsave(paste0(output_path,"phagocytosis_GWAS_manhattan_LPS.png"),
       width = 12, height = 7 )
ggsave(paste0(output_path,"phagocytosis_GWAS_manhattan_LPS.svg"),
       width = 12, height = 7 )


###########
# filter to (nominally) significant
pqtl_res %>%
  dplyr::filter(p_untreated < 10e-5 | p_IFN < 10e-5 | p_LPS < 10e-5 ) %>%
  write_csv(paste0(output_path,"nominally_significant_results.csv"))


#### splitting and saving data in preparation for for coloc

outdir = "../../../data/results/phagocytosis/2.check_association_results/lm_1pct_filtered_deflated/"
gwas = read_csv(paste0(outdir,"nominally_significant_results.csv"))
all_res = read_csv(paste0(outdir,"all_res_pqtl_phagocytosis.csv"))
selected = list()
for(treat in c("untreated","IFN","LPS")){
  
  selected[[treat]] = gwas %>% 
    dplyr::select(snp, ends_with(treat)) %>%
    dplyr::mutate(chr = as.numeric(str_split_i(snp,pattern = "_",i = 1)),
                  snp_pos = as.numeric(str_split_i(snp,pattern = "_",i = 2)),
                  REF = str_split_i(snp,pattern = "_",i = 3),
                  ALT = str_split_i(snp,pattern = "_",i = 4),
                  locus_name = snp) %>%
    dplyr::rename(coef = paste("coef",treat,sep = "_"),
                  p_value = paste("p",treat,sep = "_"),
                  se = paste("se",treat,sep = "_"),
                  variant_id = snp) %>%
    dplyr::select(variant_id,locus_name,coef,p_value,se,chr,snp_pos,REF,ALT) %>%
    dplyr::arrange(chr,snp_pos) 
  
  # save nominal results
  all_res %>%
    dplyr::select(snp, ends_with(treat)) %>%
    dplyr::mutate(chr = as.numeric(str_split_i(snp,pattern = "_",i = 1)),
                  snp_pos = as.numeric(str_split_i(snp,pattern = "_",i = 2)),
                  REF = str_split_i(snp,pattern = "_",i = 3),
                  ALT = str_split_i(snp,pattern = "_",i = 4),
                  locus_name = snp) %>%
    dplyr::rename(coef = paste("coef",treat,sep = "_"),
                  p_value = paste("p",treat,sep = "_"),
                  se = paste("se",treat,sep = "_"),
                  variant_id = snp) %>%
    dplyr::select(variant_id,locus_name,coef,p_value,se,chr,snp_pos,REF,ALT) %>%
    dplyr::arrange(chr,snp_pos)  %>%
    write_tsv(.,paste0(outdir,treat,"/",treat,"_GWAS_GRCh38.tsv.gz"))
  
}

# specific association results: directions


pqtl_res = readr::read_csv("../../../data/results/phagocytosis/2.check_association_results/lm_1pct_filtered_deflated/all_res_pqtl_phagocytosis.csv")
 

var = c("12_40189297_T_G" # LRRK2 LPS phago coloc var, rs1844922
        ) # 

pqtl_res %>%
  dplyr::filter(snp == var) %>%
  dplyr::select(contains("coef"),contains("p"),contains("snp"))


  ##### should calculate lead results by clumping based on LD to identify independent signals
  # could do with PLINK
  
  
