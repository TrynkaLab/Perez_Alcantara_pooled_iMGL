# check coloc enrichments per treatment 
# GSEA

library(tidyverse)
library(enrichplot)
library(fgsea)
library(org.Hs.eg.db)
library(ggupset)
source("./functions.R")
treatment_cols =  c(untreated = "#8D918B", IFN = "#3A5683", LPS = "#F8766D")

outdir = "../../data/results/8.colocalisation_analysis/coloc_results/"
GWAS_datasets = c("GCST90027158", # Bellenguez AD
                  "ieu-b-7_PD_Nalls_2019", # PD AD
                  "GCST005531", # MS
                  "GCST90027164") # ALS
colocs = list()
for(gwas in GWAS_datasets){
  colocs[[gwas]] = read_tsv(paste0("../../data/results/8.colocalisation_analysis/coloc_results/",gwas,
  "/coloc_summary_single_causal_variant_500kb.txt")) 
  colocs[[gwas]] = colocs[[gwas]] %>%
    dplyr::mutate(treatment = case_when(str_detect(eQTL_comparison,"untreated") == TRUE ~ "untreated",
                                              str_detect(eQTL_comparison,"LPS") == TRUE ~ "LPS",
                                              str_detect(eQTL_comparison,"IFN") == TRUE ~ "IFN"),
                  GWAS = case_when(gwas == "GCST90027158" ~ "AD",
                                   gwas ==   "ieu-b-7_PD_Nalls_2019" ~ "PD",
                                   gwas == "GCST005531" ~ "MS",
                                   gwas == "GCST90027164" ~ "ALS"))
}

colocs = do.call("rbind",colocs)

write_csv(colocs,paste0(outdir,"/all_GWAS_colocalisations.csv"))

# unique eGenes
all = colocs %>%
  dplyr::filter(PP_H4 >= 0.70) %>%
  dplyr::select(eGene,GWAS) %>%
  distinct()

all

coloc_sign_eGenes = colocs %>%
  dplyr::group_by(treatment,GWAS) %>%
  dplyr::select(eGene, treatment,GWAS,PP_H4) %>%
  dplyr::distinct(eGene,.keep_all = TRUE) %>%
dplyr::summarize(fraction = mean(PP_H4 > 0.7),
                 eGenes =  list(eGene[PP_H4 > 0.7])) %>%
  mutate(eGenes = sapply(eGenes, function(x) paste(x, collapse = ", "))) %>%
  ungroup()

coloc_number_eGenes = colocs %>%
  dplyr::group_by(treatment,GWAS) %>%
  dplyr::select(eGene, treatment,GWAS,PP_H4) %>%
  dplyr::distinct(eGene,.keep_all = TRUE) %>%
  dplyr::summarize(number = sum(PP_H4 > 0.7),
                   eGenes =  list(eGene[PP_H4 > 0.7])) %>%
  mutate(eGenes = sapply(eGenes, function(x) paste(x, collapse = ", "))) %>%
  ungroup()

coloc_sign_regions = colocs %>%
  dplyr::select(external_data_Locus, treatment,GWAS,PP_H4) %>%
  dplyr::group_by(treatment,GWAS) %>%
  dplyr::distinct(external_data_Locus,.keep_all = TRUE) %>%
  dplyr::summarize(fraction = mean(PP_H4 > 0.7),
                   loci =  list(external_data_Locus[PP_H4 > 0.7])) %>%
  mutate(loci = sapply(loci, function(x) paste(x, collapse = ", "))) %>%
  ungroup()

write_csv(coloc_sign_regions, paste0(outdir,"/coloc_sign_regions_percent.csv"))
write_csv(coloc_sign_eGenes, paste0(outdir,"/coloc_sign_eGenes_percent.csv"))


number_of_colocs = colocs %>%
  dplyr::select(external_data_Locus, treatment,GWAS,PP_H4) %>%
  dplyr::group_by(treatment,GWAS) %>%
  dplyr::distinct(external_data_Locus, .keep_all = TRUE) %>%
  dplyr::summarize(number = sum(PP_H4 > 0.7),
                   loci =  list(external_data_Locus[PP_H4 > 0.7])) %>%
  mutate(loci = sapply(loci, function(x) paste(x, collapse = ", "))) %>%
  ungroup()
# colocalisation proportions Nikos: proportion of colocalized significant regions 
# among all identified significant regions per GWAS.

# I'm not counting GWAS loci that did not have enough SNPs in common
# I'm also overcounting matches - if the colocalisation PP_H4 value is exactly 
# the same, then it's the same region - but doesn't work that well for 0 or 1 PPH4

# upsetplot coloc genes per treatment
p = colocs %>%
  dplyr::group_by(treatment,GWAS) %>%
  dplyr::select(eGene, treatment,GWAS,PP_H4) %>%
  dplyr::filter(PP_H4 > 0.7) %>%
  dplyr::distinct(eGene,.keep_all = TRUE) %>%
  dplyr::group_by(eGene,GWAS) %>%
  summarize(treatment = list(treatment))%>%
  dplyr::ungroup() %>%
  ggplot(aes(x=treatment)) +
  geom_bar() +
  scale_x_upset(n_intersections = 20) + 
  theme_minimal()

pdf(paste0(outdir,"/sign_colocs_total_shared_eGenes_upsetplot.pdf"),
    width = 5, height = 3)
plot(p)
dev.off()

p = colocs %>%
  dplyr::select(external_data_Locus, treatment,GWAS,PP_H4) %>%
  dplyr::group_by(treatment,GWAS) %>%
  dplyr::distinct(external_data_Locus,.keep_all = TRUE) %>%
  dplyr::summarize(fraction = mean(PP_H4 > 0.7),
                   loci =  list(external_data_Locus[PP_H4 > 0.7])) %>%
  mutate(loci = sapply(loci, function(x) paste(x, collapse = ", "))) %>%
  dplyr::mutate(treatment = factor(treatment,levels = c("untreated","IFN","LPS"))) %>%
  dplyr::mutate(fraction = fraction + 0.0005) %>% # to add a tiny line in barplot for samples that are 0 
  ungroup() %>%
  ggplot(., aes(fill=treatment, y=fraction *100, x=GWAS)) + 
      geom_bar(position="dodge", stat="identity") +
      scale_fill_manual(values = treatment_cols) +
      ggtitle("Percentage of loci colocalising") +
      theme_minimal() +
      xlab("") + 
  ylab("% colocalizations")

# I'm missing loci for some non-significant loci for this plot

percent = colocs %>%
  dplyr::select(external_data_Locus, treatment,GWAS,PP_H4) %>%
  dplyr::group_by(treatment,GWAS) %>%
  dplyr::distinct(external_data_Locus,.keep_all = TRUE) %>%
  dplyr::summarize(fraction = mean(PP_H4 > 0.7),
                   loci =  list(external_data_Locus[PP_H4 > 0.7])) %>%
  mutate(loci = sapply(loci, function(x) paste(x, collapse = ", "))) %>%
  dplyr::mutate(treatment = factor(treatment,levels = c("untreated","IFN","LPS")))

# total across all treatments
percent_across_treatments = colocs %>%
  dplyr::select(external_data_Locus ,GWAS,PP_H4) %>%
  dplyr::group_by(GWAS) %>%
  dplyr::distinct(external_data_Locus,.keep_all = TRUE) %>%
  dplyr::summarize(fraction = mean(PP_H4 > 0.7),
                   loci =  list(external_data_Locus[PP_H4 > 0.7])) %>%
  mutate(loci = sapply(loci, function(x) paste(x, collapse = ", ")))

percent_across_treatments

pdf(paste0(outdir,"/sign_colocs_loci_percent_barplot.pdf"),
    width = 5, height = 3)
plot(p)

dev.off()

p = colocs %>%
  dplyr::select(external_data_Locus, treatment,GWAS,PP_H4) %>%
  dplyr::group_by(treatment,GWAS) %>%
  dplyr::distinct(external_data_Locus, .keep_all = TRUE) %>%
  dplyr::summarize(number = sum(PP_H4 > 0.7),
                   loci =  list(external_data_Locus[PP_H4 > 0.7])) %>%
  mutate(loci = sapply(loci, function(x) paste(x, collapse = ", "))) %>%
  dplyr::mutate(treatment = factor(treatment,levels = c("untreated","IFN","LPS"))) %>%
  dplyr::mutate(number = number + 0.05) %>% # to add a tiny line in barplot for samples that are 0 
  ungroup() %>%
  ggplot(., aes(fill=treatment, y=number, x=GWAS)) + 
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values = treatment_cols) +
  scale_y_continuous(breaks = seq(0,12, by = 2)) +
  ggtitle("Number of loci colocalising") +
  theme_minimal() +
  xlab("") +
  ylab("Number of colocalizations")

pdf(paste0(outdir,"/sign_colocs_loci_n_barplot.pdf"),
    width = 5, height = 3)
plot(p)
dev.off()

### summary plot eQTL eGenes and GWAS loci

gwas_num = colocs %>%
  dplyr::select(external_data_Locus,GWAS) %>%
  distinct() %>%
  group_by(GWAS) %>%
  dplyr::summarise(n_loci = n())
gwas_num # this doesn't coincide with the number of loci I have initially

gwas_num$n_loci = c(83,16,126,90)

eQTL = read_csv("../../data/results/4.Inspect_eQTL_results/tensorQTL_summary.csv") %>%
  dplyr::select(variables,starts_with("35")) %>%
  pivot_longer(cols = -variables,names_to = "group",values_to = "value") %>%
  dplyr::filter(str_detect(group,"Not_proliferating") & variables == "n_genes_05") %>%
  dplyr::mutate(treatment = str_split_i(group,pattern="_",i = 2))


p1 = eQTL %>%
  ggplot(aes(x=treatment, y = value, fill = treatment)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = treatment_cols) +
  theme_minimal() +
  guides(fill="none") +
  ylab("Number of eGenes")
p1

p2 = gwas_num %>%
  ggplot(aes(x=GWAS, y = n_loci)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ylab("Number of loci")

pdf(paste0(outdir,"/sign_loci_GWAS_summary_eGenes_35PCs.pdf"),
    width = 5, height = 3)
p1 + p2

dev.off()
