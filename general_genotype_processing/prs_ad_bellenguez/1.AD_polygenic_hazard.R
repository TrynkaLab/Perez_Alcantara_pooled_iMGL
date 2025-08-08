
library(tidyverse)
library(patchwork)
library(ggrepel)
options(stringsAsFactors=F)

args = commandArgs(trailingOnly = TRUE)
hipsciVCF.file = args[1]
phsSNPDetails.file = args[2]
outputPath = args[3]

# manually loaded for testing
hipsciVCF.file = "../../data/imputed_from_kaur/microglia_samples.GRCh38.filtered.vcf.gz"
phsSNPDetails.file = "../../../resources/PGS_Bellenguez_AD/PGS002280_hmPOS_GRCh38_noheader.txt"
outputPath = "../../data/prs_ad_bellenguez"

hipsci.df = readr::read_tsv(hipsciVCF.file,comment = "##")
phs.snp.df = readr::read_tsv(phsSNPDetails.file) %>%
  dplyr::rename(log_hr = effect_weight) %>%
  dplyr::mutate(ID = paste0("chr",chr_name,"_",chr_position,"_",effect_allele,"_",other_allele), # REF is effect, ALT is other
                ID_swapped = paste0("chr",chr_name,"_",chr_position,"_",other_allele,"_",effect_allele)) # ALT is effect
# effect weight in harmonised file == log_hr


getGTString = function(gtStr) {
  sapply(strsplit(gtStr, ":", fixed=T), function(l) l[[1]])
}

fixed = hipsci.df %>%
  dplyr::select(c(`#CHROM`,POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT)) %>%
  dplyr::filter(ID %in% c(phs.snp.df$ID,phs.snp.df$ID_swapped,"chr19_44908684_T_C","chr19_44908822_C_T"))

genotype = hipsci.df %>%
  dplyr::filter(ID %in% c(phs.snp.df$ID,phs.snp.df$ID_swapped,"chr19_44908684_T_C","chr19_44908822_C_T")) %>%
  tidyr::pivot_longer(cols =!c(`#CHROM`,POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT),
                      names_to = "donor",values_to = "genotype") %>%
  dplyr::mutate(genotype = str_split(genotype, ":", simplify = TRUE)[,1]) %>%
  dplyr::mutate(allele1 = as.numeric(str_split_i(genotype, "\\|", i=1)),
                allele2=as.numeric(str_split_i(genotype, "\\|", i=2))) %>%
  dplyr::mutate(gtCount = allele1 + allele2,
                line = ifelse(grepl("HPS",donor),
                              no = str_split_i(donor,pattern = "-",i = 1),
                              yes = str_split_i(donor,pattern = "-",i = 2)))
  
  
  hipsci.apoe.df = genotype %>%
  dplyr::group_by(donor) %>%
  dplyr::filter(ID %in% c("chr19_44908684_T_C","chr19_44908822_C_T")) %>%
  tidyr::pivot_wider(id_cols = line,names_from = c(ID),values_from = c(allele1,allele2)) %>%
  group_by(line) %>%
  dplyr::summarise(apoe_allele1 = case_when(allele1_chr19_44908684_T_C== "0" & allele1_chr19_44908822_C_T ==  "0" ~ "e3",
                                            allele1_chr19_44908684_T_C== "0" &  allele1_chr19_44908822_C_T ==  "1" ~ "e2",
                                            allele1_chr19_44908684_T_C == "1" &  allele1_chr19_44908822_C_T == "0" ~ "e4",
                                            allele1_chr19_44908684_T_C == "1" &  allele1_chr19_44908822_C_T ==  "1" ~ "E5"), # error
                   apoe_allele2 = case_when(allele2_chr19_44908684_T_C == "0" & allele2_chr19_44908822_C_T == "0" ~ "e3",
                                            allele2_chr19_44908684_T_C == "0" &  allele2_chr19_44908822_C_T == "1"~ "e2",
                                            allele2_chr19_44908684_T_C == "1" & allele2_chr19_44908822_C_T== "0" ~ "e4",
                                            allele2_chr19_44908684_T_C == "1" & allele2_chr19_44908822_C_T ==  "1" ~ "E5"))

table(hipsci.apoe.df$apoe_allele1)
table(hipsci.apoe.df$apoe_allele2)

# Check allele frequency of APOE SNPs for sanity
# checking allele frequency E2/E2 (1%), E2/E3 (11%), E2/E4 (2%), E3/E3 (61%), E3/E4 (23%),E4/E4 (2%)
(sum(hipsci.apoe.df$apoe_allele1 == "e2" & hipsci.apoe.df$apoe_allele2== "e2")/ nrow(hipsci.apoe.df))*100 # ~ 0.8%, OK
(sum(hipsci.apoe.df$apoe_allele1 == "e4" & hipsci.apoe.df$apoe_allele2== "e4")/ nrow(hipsci.apoe.df))*100  # ~ 4%, OK
(sum(hipsci.apoe.df$apoe_allele1 == "e3" & hipsci.apoe.df$apoe_allele2== "e3")/ nrow(hipsci.apoe.df))*100 # ~ 58%, OK


# Compute hazard ratios. First do it for APOE, then for other SNPs,
# and finally combine these together.
# betas from here: https://www.nature.com/articles/s41467-021-24082-z#Sec8

apoe.snp.df = data.frame(ID = c("e2", "e3", "e4"), hr = c(exp(-0.47), exp(0), exp(1.12)))
hipsci.apoe.hr = hipsci.apoe.df %>%
  dplyr::left_join(apoe.snp.df %>% dplyr::rename(allele1HR = hr), by=c("apoe_allele1" = "ID")) %>%
  dplyr::left_join(apoe.snp.df %>% dplyr::rename(allele2HR = hr), by=c("apoe_allele2" = "ID")) %>%
  dplyr::mutate(apoeHR = allele1HR * allele2HR)

# ensure the alleles are in right order so I extract the log_hr with the correct sign
phs.snp.df = phs.snp.df %>%
  dplyr::mutate(log_hr_fixed = case_when(ID %in% genotype$ID ~ -log_hr, # in the genotype tibble the 1/1 is actually other/other for these positions
                                         .default = log_hr))

hipsci.all.hr = genotype %>%
  dplyr::left_join(phs.snp.df %>% dplyr::select(ID, log_hr_fixed), by="ID") %>%
  na.omit() %>%
  dplyr::group_by(line) %>%
  dplyr::summarise(polygenicHR = exp(sum(log_hr_fixed * gtCount))) %>%
  dplyr::left_join(hipsci.apoe.hr %>% dplyr::select(line, apoe_allele1, apoe_allele2, apoeHR), by="line") %>%
  dplyr::mutate(overallHR = polygenicHR * apoeHR) %>%
  dplyr::arrange(-overallHR)


getmode <- function(v) {
  uniqv = unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

p1 = ggplot(hipsci.all.hr, aes(x=polygenicHR)) + geom_histogram() + theme_minimal() +ggtitle("polygenic risk score")
p2 = ggplot(hipsci.all.hr, aes(x=overallHR)) + geom_histogram() + theme_minimal() + ggtitle("polygenic*APOE risk score")

pdf(file=paste0(outputPath, "/hist_polygenic_hazard_hipsci.pdf"),width = 7,height = 3)
p1 + p2
dev.off()


# checking where these fall in the quantiles of hipsci PRS from Bellenguez

hipsci.all.hr$polygenicHR_quantile = ecdf(hipsci.all.hr$polygenicHR)(hipsci.all.hr$polygenicHR)
hipsci.all.hr$overallHR_quantile = ecdf(hipsci.all.hr$overallHR)(hipsci.all.hr$overallHR)

## Binning per quartile
hipsci.all.hr$overallHR_quartile = rep(NA, nrow(hipsci.all.hr))
hipsci.all.hr[hipsci.all.hr$overallHR_quantile <= 0.25,"overallHR_quartile"] = "Q1"
hipsci.all.hr[hipsci.all.hr$overallHR_quantile <= 0.50 & hipsci.all.hr$overallHR_quantile > 0.25,"overallHR_quartile"] = "Q2"
hipsci.all.hr[hipsci.all.hr$overallHR_quantile <= 0.75 & hipsci.all.hr$overallHR_quantile > 0.50,"overallHR_quartile"] = "Q3"
hipsci.all.hr[hipsci.all.hr$overallHR_quantile > 0.75,"overallHR_quartile"] = "Q4"

hipsci.all.hr$polygenicHR_quartile = rep(NA, nrow(hipsci.all.hr))
hipsci.all.hr[hipsci.all.hr$polygenicHR_quantile <= 0.25,"polygenicHR_quartile"] = "Q1"
hipsci.all.hr[hipsci.all.hr$polygenicHR_quantile <= 0.50 & hipsci.all.hr$polygenicHR_quantile > 0.25,"polygenicHR_quartile"] = "Q2"
hipsci.all.hr[hipsci.all.hr$polygenicHR_quantile <= 0.75 & hipsci.all.hr$polygenicHR_quantile > 0.50,"polygenicHR_quartile"] = "Q3"
hipsci.all.hr[hipsci.all.hr$polygenicHR_quantile > 0.75,"polygenicHR_quartile"] = "Q4"

hipsci.all.hr =  hipsci.all.hr %>%
  dplyr::arrange(overallHR) %>%
  dplyr::mutate(overallHR_rank = row_number()) 

write.table(hipsci.all.hr, file=paste0(outputPath,"/Jeremy_AD_polygenic_hazard.hipsci_ipmar_donors.txt"),
            sep="\t", quote=F, col.names=T, row.names=F)

# compare PRSice and this

prsice = read_tsv("../../data/prs_ad_bellenguez/PRSice/hipsci_polygenic_score_AD_Bellenguez_withAPOE.tsv")

full = prsice %>%
  dplyr::select(line,full_PRS_rank, full_PRS_quantile) %>%
  dplyr::left_join(hipsci.all.hr[,c("line","overallHR_rank","overallHR_quantile")]) %>%
  dplyr::mutate(rank_change = full_PRS_rank - overallHR_rank)

p1 = full %>%
  ggplot(aes(x = full_PRS_rank, y = rank_change, label = line)) + 
  geom_point() +
  geom_text_repel() +
  theme_minimal() +
  xlab("PRSice AD PRS") + 
  ylab("PRSice AD PRS - Jeremy AD PRS")


p2 = full %>%
  ggplot(aes( y = rank_change)) + 
  geom_boxplot() +
  theme_minimal() +
  ylab("PRSice AD PRS - Jeremy AD PRS") +
  theme(axis.text = element_blank(),
        axis.title = element_blank())

pdf("../../data/prs_ad_bellenguez/change_PRS_rank_PRSice_vs_Jeremy.pdf",
    width = 7, height = 5)
p1 + p2 +
  plot_layout(widths = c(3, 1))
dev.off()

######

p1 = full %>%
  ggplot(.,aes(x=overallHR_quantile,y=full_PRS_quantile, label = line)) +
  theme_minimal() +
  geom_rect(xmin=0, xmax=0.25, ymin=0, ymax=0.25, fill="#CDF5CF", alpha=0.1) + 
  geom_rect(xmin = 0.75, xmax = 1, ymin = 0.75, ymax = 1, fill="#CDF5CF", alpha=0.1) +   
  geom_rect(xmin=0.25, xmax=0.5, ymin=0.25, ymax=0.5, fill="#CDF5CF", alpha=0.1) + 
  geom_rect(xmin = 0.5, xmax = 0.75, ymin = 0.5, ymax = 0.75, fill="#CDF5CF", alpha=0.1) +  
  geom_rect(xmin=0, xmax=0.25, ymin=0.25, ymax=0.5, fill="#FAE0A0", alpha=0.1) + 
  geom_rect(xmin = 0.25, xmax = 0.5, ymin = 0, ymax = 0.25, fill="#FAE0A0", alpha=0.1) +  
  geom_rect(xmin=0.50, xmax=0.75, ymin=0.25, ymax=0.5, fill="#FAE0A0", alpha=0.1) + 
  geom_rect(xmin = 0.25, xmax = 0.5, ymin = 0.50, ymax = 0.75, fill="#FAE0A0", alpha=0.1) +  
  geom_rect(xmin = 0.75, xmax = 1, ymin = 0.5, ymax = 0.75, fill="#FAE0A0", alpha=0.1) +   
  geom_rect(xmin = 0.5, xmax = 0.75, ymin = 0.75, ymax = 1, fill="#FAE0A0", alpha=0.1) + 
  geom_rect(xmin = 0, xmax =0.25, ymin = 0.5, ymax = 0.75, fill="#F7BBA3", alpha=0.1) +   
  geom_rect(xmin = 0.5, xmax = 0.75, ymin = 0, ymax = 0.25, fill="#F7BBA3", alpha=0.1) +  
  geom_rect(xmin = 0.25, xmax =0.5, ymin = 0.75, ymax = 1, fill="#F7BBA3", alpha=0.1) +   
  geom_rect(xmin = 0.75, xmax = 1, ymin = 0.25, ymax = 0.5, fill="#F7BBA3", alpha=0.1) + 
  geom_rect(xmin = 0, xmax =0.25, ymin = 0.75, ymax = 1, fill="#F77979", alpha=0.1) +   
  geom_rect(xmin = 0.75, xmax =1, ymin = 0, ymax = 0.25, fill="#F77979", alpha=0.1) +  
  geom_point() +
  geom_text_repel() 

pdf(paste0(outputPath,"/Prsice_vs_Jeremy_PRS_full_PRS_quantile.pdf"),width = 7,height = 7)
p1 + ggtitle("PRSice vs Jeremy full PRS quantile")
dev.off()
