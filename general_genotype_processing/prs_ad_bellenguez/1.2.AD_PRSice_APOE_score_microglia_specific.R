# Incorporating non-APOE PRSice score (from microglia-specific regions from coloc) to APOE part
# then normalising to 1000 Genomes

library(tidyverse)
library(vcfR)
library(ggrepel)
library(patchwork)
# load prsice file with format FID, IID, PRS for first set at first threshold, PRS for first set at second threshold...

for(treatment in c("untreated","LPS","IFN")){
  

hipsci_IPMAR_prsice = read_delim(paste0("../../data/prs_ad_bellenguez/PRSice/PRSice_noAPOE_microglia_",
                                        treatment,".all_score"),delim = " ") %>%
  dplyr::select(FID,IID,if_else("Pt_0.1" %in% colnames(.), "Pt_0.1",
                                if_else("Pt_0.0974501" %in% colnames(.),"Pt_0.0974501","Pt_0.0970501")))  %>% # sometimes the rounding is not quite 0.1
  dplyr::rename(polygenic_score = any_of(c("Pt_0.1","Pt_0.0970501","Pt_0.0974501")),
                line = IID) %>%
  dplyr::mutate(set = "HipSci_IPMAR",
                line = ifelse(grepl("HPS",line),
                no = str_split_i(line,pattern = "-",i = 1),
                yes = str_split_i(line,pattern = "-",i = 2)))
# we are interested in the pT 0.1 threshold that showed the best AUC in Valentina's study for AD PRS

thousand_prsice = read_delim(paste0("../../data/prs_ad_bellenguez/PRSice/1000G_PRSice_noAPOE_microglia_",
                             treatment,".all_score"),delim = " ") %>%
  dplyr::select(FID,IID,if_else("Pt_0.1" %in% colnames(.), "Pt_0.1",
                                if_else("Pt_0.0974501" %in% colnames(.),"Pt_0.0974501","Pt_0.0970501")))  %>% # sometimes the rounding is not quite 0.1
  dplyr::rename(polygenic_score = any_of(c("Pt_0.1","Pt_0.0970501","Pt_0.0974501")),
                line = IID) %>%
  dplyr::mutate(set = "1000G EUR")

# adjust PRS for PCs and take residuals

# loading PCs
pcs_hipsci_ipmar = read_tsv("../../data/prs_ad_bellenguez/microglia_samples.GRCh38.filtered.eigenvec") %>%
  dplyr::select(-`#FID`) %>%
  dplyr::rename(line = IID) %>%
  dplyr::mutate(set = "HipSci_IPMAR",
                line = ifelse(grepl("HPS",line),
                              no = str_split_i(line,pattern = "-",i = 1),
                              yes = str_split_i(line,pattern = "-",i = 2)))

pcs_1kg = read_tsv("../../data/prs_ad_bellenguez/1000G_EUR_maf05perc_biallelic_nonmissing_nochr_withrsid_checked.eigenvec") %>%
  dplyr::select(-`#FID`) %>%
  dplyr::rename(line = IID)


hipsci_IPMAR_prsice = hipsci_IPMAR_prsice %>%
  left_join(pcs_hipsci_ipmar)

thousand_prsice = thousand_prsice %>%
  left_join(pcs_1kg)
# regression

hipsci_adj = lm(formula = polygenic_score ~ PC1 + PC2 + PC3 + PC4 + PC5 ,
         data = hipsci_IPMAR_prsice)

hipsci_IPMAR_prsice = hipsci_IPMAR_prsice %>%
dplyr::mutate(residuals =  residuals.lm(hipsci_adj))

kg_adj = lm(formula = polygenic_score ~ PC1 + PC2 + PC3 + PC4 + PC5 ,
                data = thousand_prsice)

thousand_prsice = thousand_prsice %>%
  dplyr::mutate(residuals =  residuals.lm(kg_adj))

# scale to 1000 Genomes EUR ########


all_prsice = rbind(hipsci_IPMAR_prsice,thousand_prsice)


pdf(paste0("../../data/prs_ad_bellenguez/PRSice/hipsci_1000g_raw_polygenic_score_AD_Bellenguez_noAPOE_microglia_",
           treatment,".pdf"),
    width = 7, height = 5)
ggplot(all_prsice, aes(x = polygenic_score, fill = set)) +
  scale_fill_manual(values = c("HipSci_IPMAR" = "black", "1000G EUR" = "lightgrey")) + 
  geom_histogram(alpha = 0.6) + 
  theme_minimal() + 
  ggtitle("PRSice: Bellenguez AD", subtitle = "For pT = 0.1") 

# is it normal that the raw outputs are so tiny?
# apparently, yes https://github.com/choishingwan/PRSice/issues/180
dev.off()

# standardising against 1000G mean and sd
thousand_mean = mean(thousand_prsice$polygenic_score)
thousand_sd = sd(thousand_prsice$polygenic_score)

# using custom function because default scale is causing problems
scale_this <- function(x,mean,stdev){
  (x - mean) / stdev
}

all_prsice = all_prsice %>%
  group_by(set) %>%
  dplyr::mutate(prs_scaled = scale_this(polygenic_score, 
                                   mean = thousand_mean,
                                   stdev = thousand_sd)) %>%
  ungroup()


pdf(paste0("../../data/prs_ad_bellenguez/PRSice/hipsci_1000g_scaled_polygenic_score_AD_Bellenguez_noAPOE_microglia_",
           treatment,".pdf"),
    width = 7, height = 5)
ggplot(all_prsice, aes(x = prs_scaled, fill = set)) +
  geom_histogram(alpha = 0.6) + 
  scale_fill_manual(values = c("HipSci_IPMAR" = "black", "1000G EUR" = "lightgrey")) + 
  theme_minimal() + 
  ggtitle("PRSice: Bellenguez AD", subtitle = "For pT = 0.1") 

dev.off()

hipsci_ipmar_apoe_vcf = vcfR::read.vcfR("../../data/prs_ad_bellenguez/hipsci_APOE.vcf")
thousand_apoe_vcf = vcfR::read.vcfR("../../data/prs_ad_bellenguez/1000G_APOE.vcf")
variants = hipsci_ipmar_apoe_vcf@fix[,3]
hipsci_ipmar_apoe = vcfR::extract.gt(hipsci_ipmar_apoe_vcf, return.alleles = TRUE) %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(SNPs = variants) %>%
  dplyr::relocate(SNPs) %>%
  tidyr::pivot_longer(cols=-SNPs,names_to = "donor",values_to = "alleles") %>%
  # reordering heterozygotes - helpful for translating to apoe alleles later
  dplyr::mutate(alleles = case_when(alleles == "C/T" ~ "T/C",
                                   .default = alleles)) %>%

dplyr::mutate(allele1 = str_split_i(alleles,pattern = "/",i = 1),
              allele2 = str_split_i(alleles,pattern = "/",i = 2),
              line = ifelse(grepl("HPS",donor),
                            no = str_split_i(donor,pattern = "-",i = 1),
                            yes = str_split_i(donor,pattern = "-",i = 2)))
  

hipsci_ipmar_apoe = hipsci_ipmar_apoe %>%
  tidyr::pivot_wider(id_cols = line,names_from = c(SNPs),values_from = c(allele1,allele2)) %>%
  group_by(line) %>%
  dplyr::summarise(apoe_allele1 = case_when(allele1_chr19_44908684_T_C== "T" & allele1_chr19_44908822_C_T ==  "C" ~ "e3",
                                            allele1_chr19_44908684_T_C == "T" &  allele1_chr19_44908822_C_T ==  "T" ~ "e2",
                                            allele1_chr19_44908684_T_C == "C" &  allele1_chr19_44908822_C_T == "C" ~ "e4",
                                            allele1_chr19_44908684_T_C == "C" &  allele1_chr19_44908822_C_T ==  "T" ~ "E5"), # error
                   apoe_allele2 = case_when(allele2_chr19_44908684_T_C == "T" & allele2_chr19_44908822_C_T == "C" ~ "e3",
                                            allele2_chr19_44908684_T_C == "T" &  allele2_chr19_44908822_C_T == "T"~ "e2",
                                            allele2_chr19_44908684_T_C == "C" & allele2_chr19_44908822_C_T== "C" ~ "e4",
                                            allele2_chr19_44908684_T_C == "C" & allele2_chr19_44908822_C_T ==  "T" ~ "E5")) %>%
  dplyr::mutate(set = "HipSci_IPMAR")

table(hipsci_ipmar_apoe$apoe_allele1,hipsci_ipmar_apoe$apoe_allele2)
# checking allele frequency
e2 = (sum(table(hipsci_ipmar_apoe$apoe_allele1,hipsci_ipmar_apoe$apoe_allele2)[1,]) / nrow(hipsci_ipmar_apoe))/2 # ~ 8%, OK
e4 = (sum(table(hipsci_ipmar_apoe$apoe_allele1,hipsci_ipmar_apoe$apoe_allele2)[,3]) / nrow(hipsci_ipmar_apoe)) / 2 # ~ 15%, OK
e2
e4

thousand_apoe =  vcfR::extract.gt(thousand_apoe_vcf, return.alleles = TRUE) %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(SNPs = variants) %>%
  dplyr::relocate(SNPs) %>%
  tidyr::pivot_longer(cols=-SNPs,names_to = "line",values_to = "alleles") %>%
  # reordering heterozygotes - helpful for translating to apoe alleles later
  dplyr::mutate(alleles = case_when(alleles == "C/T" ~ "T/C",
                                    .default = alleles)) %>%
  
  dplyr::mutate(allele1 = str_split_i(alleles,pattern = "/",i = 1),
                allele2 = str_split_i(alleles,pattern = "/",i = 2))



thousand_apoe = thousand_apoe %>%
  tidyr::pivot_wider(id_cols = line,names_from = c(SNPs),values_from = c(allele1,allele2)) %>%
  group_by(line) %>%
  dplyr::summarise(apoe_allele1 = case_when(allele1_chr19_44908684_T_C== "T" & allele1_chr19_44908822_C_T ==  "C" ~ "e3",
                                            allele1_chr19_44908684_T_C == "T" &  allele1_chr19_44908822_C_T ==  "T" ~ "e2",
                                            allele1_chr19_44908684_T_C == "C" &  allele1_chr19_44908822_C_T == "C" ~ "e4",
                                            allele1_chr19_44908684_T_C == "C" &  allele1_chr19_44908822_C_T ==  "T" ~ "E5"), # error
                   apoe_allele2 = case_when(allele2_chr19_44908684_T_C == "T" & allele2_chr19_44908822_C_T == "C" ~ "e3",
                                            allele2_chr19_44908684_T_C == "T" &  allele2_chr19_44908822_C_T == "T"~ "e2",
                                            allele2_chr19_44908684_T_C == "C" & allele2_chr19_44908822_C_T== "C" ~ "e4",
                                            allele2_chr19_44908684_T_C == "C" & allele2_chr19_44908822_C_T ==  "T" ~ "E5")) %>%
  dplyr::mutate(set = "1000G EUR")

table(thousand_apoe$apoe_allele1,thousand_apoe$apoe_allele2)
(sum(table(thousand_apoe$apoe_allele1,thousand_apoe$apoe_allele2)[1,]) / nrow(thousand_apoe))/2 # ~ 6.5%, OK
 (sum(table(thousand_apoe$apoe_allele1,thousand_apoe$apoe_allele2)[,3]) / nrow(thousand_apoe)) / 2 # ~ 14%, OK

all_APOE = rbind(hipsci_ipmar_apoe,thousand_apoe)

# getting variants in beta score format
# weighted sum as in here: https://www.nature.com/articles/s41467-021-24082-z#Sec8
# effect sizes APOE: effect sizes (B(ε2)=−0.47 and B(ε4)=1.12) 

all_prs_with_APOE = all_prsice %>%
  dplyr::left_join(.,all_APOE) %>%
  dplyr::mutate(APOE_sum = case_when(apoe_allele1 == "e2" & apoe_allele2 == "e2" ~   -0.47 - 0.47,
                                            apoe_allele1 == "e3" & apoe_allele2 == "e2" ~   -0.47,
                                            apoe_allele1 == "e2" & apoe_allele2 == "e3" ~   -0.47,
                                            apoe_allele1 == "e4" & apoe_allele2 == "e2" ~   -0.47 + 1.12,
                                            apoe_allele1 == "e2" & apoe_allele2 == "e4" ~   -0.47 + 1.12,
                                            apoe_allele1 == "e3" & apoe_allele2 == "e4" ~   1.12,
                                            apoe_allele1 == "e4" & apoe_allele2 == "e3" ~   1.12,
                                            apoe_allele1 == "e4" & apoe_allele2 == "e4" ~   1.12 + 1.12,
                                            .default = 0))


pdf(paste0("../../data/prs_ad_bellenguez/PRSice/hipsci_1000g_raw_APOE_score_microglia_",
           treatment,".pdf"),
    width = 7, height = 5)
ggplot(all_prs_with_APOE, aes(x = APOE_sum, fill = set)) +
  scale_fill_manual(values = c("HipSci_IPMAR" = "black", "1000G EUR" = "lightgrey")) + 
  geom_histogram(alpha = 0.6) + 
  theme_minimal() + 
  ggtitle("Unscaled APOE scores") 

dev.off()

# scale APOE to 1000 Genomes
thousand= all_prs_with_APOE %>%
  dplyr::filter(set == "1000G EUR") %>%
dplyr::summarise(mean = mean(APOE_sum), sd = sd(APOE_sum))

all_prs_with_APOE = all_prs_with_APOE %>%
  dplyr::mutate(APOE_sum_scaled = scale_this(APOE_sum, 
                                   mean = thousand$mean,
                                   stdev = thousand$sd)) %>%
  dplyr::group_by(line) %>%
  dplyr::mutate(full_PRS = sum(prs_scaled,APOE_sum_scaled))



pdf(paste0("../../data/prs_ad_bellenguez/PRSice/hipsci_1000g_scaled_APOE_score_microglia_",
           treatment,".pdf"),
    width = 7, height = 5)
ggplot(all_prs_with_APOE, aes(x = APOE_sum_scaled, fill = set)) +
  scale_fill_manual(values = c("HipSci_IPMAR" = "black", "1000G EUR" = "lightgrey")) + 
  geom_histogram(alpha = 0.6) + 
  theme_minimal() + 
  ggtitle("Scaled APOE scores") 

dev.off()


medians_means = all_prs_with_APOE %>%
  group_by(set) %>%
  dplyr::summarise(median = median(full_PRS),
                   mean = mean(full_PRS)) %>%
  ungroup()

# test differences in means:
t_test_result = all_prs_with_APOE %>%
  group_by(set) %>%
  summarize(full_PRS = list(full_PRS)) %>%
  ungroup() %>%
  summarize(p_value = t.test(full_PRS[[1]], full_PRS[[2]])$p.value)

# test differences in medians:
wilcox_test_result = all_prs_with_APOE %>%
  group_by(set) %>%
  summarize(full_PRS = list(full_PRS)) %>%
  ungroup() %>%
  summarize(p_value = wilcox.test(full_PRS[[1]], full_PRS[[2]])$p.value)

# both are significant at alpha < 0.05 but not 0.01
medians_means$median[1] -medians_means$median[2] # 0.19
medians_means$mean[1] -medians_means$mean[2] # 0.24

pdf(paste0("../../data/prs_ad_bellenguez/PRSice/hipsci_full_polygenic_score_AD_Bellenguez_withAPOE_microglia_",
           treatment,".pdf"),
    width = 7, height = 5)
ggplot(all_prs_with_APOE, aes(x = full_PRS, fill = set)) +
  scale_fill_manual(values = c("HipSci_IPMAR" = "black", "1000G EUR" = "lightgrey")) + 
  geom_histogram(bins = 100,alpha = 0.8) + 
  geom_vline(xintercept = medians_means$median[1], col = "darkgrey", 
             linetype = "dashed",linewidth = 0.8)+
  geom_vline(xintercept = medians_means$median[2], col = "black",
             linetype = "dashed",linewidth = 0.8)+
  theme_minimal() + 
  ggtitle("PRSice: Bellenguez AD", subtitle = "Including APOE. For pT = 0.1")

dev.off()

# ranking by standardised polygenic score only, and including APOE

all_prs_with_APOE = all_prs_with_APOE %>%
  dplyr::group_by(set) %>%
  dplyr::arrange(prs_scaled) %>%
  dplyr::mutate(polygenic_rank = row_number()) %>%
  dplyr::arrange(full_PRS) %>%
  dplyr::mutate(full_PRS_rank = row_number()) %>%
  ungroup()


# check how much the ranks change

pdf(paste0("../../data/prs_ad_bellenguez/PRSice/hipsci_full_polygenic_score_AD_Bellenguez_withAPOE_changing_ranks_microglia_",
treatment,".pdf"),
    width = 7, height = 5)

all_prs_with_APOE %>%
  dplyr::mutate(alleles = paste0(apoe_allele1,"/",apoe_allele2)) %>%
  ggplot(aes(x = polygenic_rank, y = full_PRS_rank, col = set, shape = alleles)) +
  scale_color_manual(values = c("HipSci_IPMAR" = "black", "1000G EUR" = "darkgrey")) + 
  geom_point(alpha = 0.8) +
  theme_minimal()

dev.off()

pdf(paste0("../../data/prs_ad_bellenguez/PRSice/hipsci_full_polygenic_score_AD_Bellenguez_withAPOE_changing_ranks_our_data_only_microglia_",
           treatment,".pdf"),
    width = 7, height = 5)

all_prs_with_APOE %>%
  dplyr::filter(set == "HipSci_IPMAR") %>%
  dplyr::mutate(alleles = paste0(apoe_allele1,"/",apoe_allele2)) %>%
  ggplot(aes(x = polygenic_rank, y = full_PRS_rank, col = alleles)) +

  scale_color_brewer(palette = "Reds") + 
  geom_point(alpha = 0.8) +
  theme_minimal()
dev.off()

all_prs_with_APOE = all_prs_with_APOE %>%
  dplyr::filter(set == "HipSci_IPMAR") 

all_prs_with_APOE$polygenic_quantile = ecdf(all_prs_with_APOE$prs_scaled)(all_prs_with_APOE$prs_scaled)
all_prs_with_APOE$full_PRS_quantile = ecdf(all_prs_with_APOE$full_PRS)(all_prs_with_APOE$full_PRS)

## Binning per quartile
all_prs_with_APOE$full_PRS_quartile = rep(NA, nrow(all_prs_with_APOE))
all_prs_with_APOE[all_prs_with_APOE$full_PRS_quantile <= 0.25,"full_PRS_quartile"] = "Q1"
all_prs_with_APOE[all_prs_with_APOE$full_PRS_quantile <= 0.50 & all_prs_with_APOE$full_PRS_quantile > 0.25,"full_PRS_quartile"] = "Q2"
all_prs_with_APOE[all_prs_with_APOE$full_PRS_quantile <= 0.75 & all_prs_with_APOE$full_PRS_quantile > 0.50,"full_PRS_quartile"] = "Q3"
all_prs_with_APOE[all_prs_with_APOE$full_PRS_quantile > 0.75,"full_PRS_quartile"] = "Q4"

all_prs_with_APOE$prs_scaled_quartile = rep(NA, nrow(all_prs_with_APOE))
all_prs_with_APOE[all_prs_with_APOE$polygenic_quantile <= 0.25,"prs_scaled_quartile"] = "Q1"
all_prs_with_APOE[all_prs_with_APOE$polygenic_quantile <= 0.50 & all_prs_with_APOE$polygenic_quantile > 0.25,"prs_scaled_quartile"] = "Q2"
all_prs_with_APOE[all_prs_with_APOE$polygenic_quantile <= 0.75 & all_prs_with_APOE$polygenic_quantile > 0.50,"prs_scaled_quartile"] = "Q3"
all_prs_with_APOE[all_prs_with_APOE$polygenic_quantile > 0.75,"prs_scaled_quartile"] = "Q4"

all_prs_with_APOE %>%
write_tsv(paste0("../../data/prs_ad_bellenguez/PRSice/hipsci_polygenic_score_AD_Bellenguez_withAPOE_microglia_",
                 treatment,".tsv"))


p1 = all_prs_with_APOE %>%
  dplyr::filter(set == "HipSci_IPMAR") %>%
  ggplot(aes(x=full_PRS_rank,y=full_PRS,label = line)) +
  geom_point() + 
  theme_minimal() + 
  geom_text_repel() + 
  ggtitle("HipSci + IPMAR PRS distribution")

p2 = all_prs_with_APOE %>%
  ggplot(aes(y=full_PRS, group = set, col = set)) +
  geom_boxplot(notch = TRUE) + 
  theme_minimal() +
  scale_color_manual(values = c("HipSci_IPMAR" = "black", "1000G EUR" = "lightgrey")) + 
  theme(axis.text = element_blank(),
        axis.title = element_blank())

pdf(paste0("../../data/prs_ad_bellenguez/PRSice/hipsci_full_polygenic_score_AD_Bellenguez_withAPOE_distribution_microglia_",
           treatment,".pdf"),
    width = 7, height = 5)
p1 
dev.off()

all_prs_with_APOE %>%
  summarise(mean_PRS = mean(full_PRS), # -0.0716 incl 1000G, -0.245 for hipsci + ipmar only
            sd_PRS = sd(full_PRS), # 1.45 incl 1000G, 1.55 for hipsci + ipmar only
            over_2sd = mean_PRS + 2*sd_PRS, #  2.83 incl 1000G,  2.85 for hipsci + ipmar only
            under_2sd = mean_PRS - 2*sd_PRS) # -2.97 incl 1000G, -3.34  for hipsci + ipmar only
}


## comparison of these microglia and treatment-specific PRS with the PRS generated from the full AD GWAS set

full_AD_PRS = read_tsv("../../data/prs_ad_bellenguez/PRSice/hipsci_polygenic_score_AD_Bellenguez_withAPOE.tsv")
full_AD_PRS$group = "full_AD"

microglia_AD_PRS = list()
for(treatment in c("untreated","LPS","IFN")){
  microglia_AD_PRS[[treatment]] = read_tsv(paste0("../../data/prs_ad_bellenguez/PRSice/hipsci_polygenic_score_AD_Bellenguez_withAPOE_microglia_",
  treatment,".tsv"))
  microglia_AD_PRS[[treatment]]$group = treatment

}

microglia_AD_PRS = do.call("rbind",microglia_AD_PRS)

p1 = microglia_AD_PRS %>%
  dplyr::bind_rows(full_AD_PRS) %>%
  dplyr::group_by(group) %>%
  dplyr::arrange(full_PRS_rank) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(line = factor(line,levels = unique(line),ordered = TRUE)) %>%
  ggplot(aes(x=full_PRS_rank,y=line, color = group)) +
  geom_point() + 
  theme_minimal() + 
  ggtitle("HipSci + IPMAR PRS distribution")
p2 = microglia_AD_PRS %>%
  dplyr::bind_rows(full_AD_PRS) %>%
  dplyr::group_by(group) %>%
  dplyr::arrange(full_PRS_rank) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(line = factor(line,levels = unique(line),ordered = TRUE)) %>%
  dplyr::slice_tail(n = 30, by = group) %>%
  ggplot(aes(x=full_PRS_rank,y=line, color = group)) +
  geom_point() + 
  theme_minimal() + 
  ggtitle("HipSci + IPMAR PRS distribution: top 30")

pdf(paste0("../../data/prs_ad_bellenguez/PRSice/full_PRS_vs_microglia_treatments_ranks.pdf"),
    width = 10, height = 20)
p1 + p2
dev.off()

pdf(paste0("../../data/prs_ad_bellenguez/PRSice/full_PRS_vs_microglia_treatments_ranks_top30.pdf"),
    width = 5, height = 7)
 p2
dev.off()

pdf(paste0("../../data/prs_ad_bellenguez/PRSice/full_PRS_vs_microglia_treatments_histogram.pdf"),
    width = 4, height = 4)
microglia_AD_PRS %>%
  dplyr::bind_rows(full_AD_PRS) %>%
    ggplot(aes(x=full_PRS, fill = group)) +
  geom_histogram(alpha = 0.6) + 
  scale_fill_brewer(type = "qual",palette = "Set1") +
  theme_minimal() + 
  ggtitle("HipSci + IPMAR PRS scores")
dev.off()