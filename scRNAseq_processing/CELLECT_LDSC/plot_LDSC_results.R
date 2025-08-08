# plot CELLECT-LDSC results
library(tidyverse)
library(pheatmap)
library(cowplot)
library(ggplotify)

outdir = "../../../data/results/CELLECT/CELLECT-LDSC/"
treatment_cols =  c(untreated = "#8D918B", IFN = "#3A5683", LPS = "#F8766D")
#https://hbctraining.github.io/publication_perfect/lessons/08_figure_specific_packages.html
# specificity IDs of interest
ids_of_interest = c(
                    "microglia-cellex",
                    #"microglia-treatment-DEG-t-bottom","microglia-treatment-DEG-t-top",
                    "microglia-treatment-eGenes-lfsr")

gwas_of_interest = c( "AD_Bellenguez" ,  "PD_Nalls", "MS_Beecham","LBD_Chia" , "ALS_van_Rheenen" , "BMI_Yengo"   )


## ggplot solution with grid ########

prioritization = read_csv("../../../data/results/CELLECT/CELLECT-LDSC/results/prioritization.csv") %>%
  dplyr::arrange(pvalue) %>%
  dplyr::filter(specificity_id %in% ids_of_interest & gwas %in% gwas_of_interest) %>%
  dplyr::group_by(gwas,specificity_id) %>%  # Tobi adjusted across annotation category (so, within [grouping by] gwas and specificity_id), because those were his independent tests
  dplyr::mutate(p_bonferroni = p.adjust(pvalue,method = "bonferroni")) %>%
  dplyr::ungroup() %>%
  dplyr::rename(treatment = annotation , GWAS = gwas, `gene list` = specificity_id) %>% ## rename categories
  dplyr::mutate(treatment = factor(treatment, levels = c("untreated","IFN","LPS")),
                minus_log10_pval = -log10(pvalue),
                GWAS = factor(str_split_i(GWAS,"_",i=1),levels = str_split_i(gwas_of_interest,"_",i=1)),
                `gene list` = case_when(`gene list` == "microglia-cellex" ~ "treatment-specific genes",
                                        `gene list` =="microglia-PRS-DEG-t-higher" ~ "AD PRS high risk DEGs",
                                        `gene list` =="microglia-PRS-DEG-t-lower" ~ "AD PRS low risk DEGs",
                                        `gene list` =="microglia-APOE-DEG-t-higher" ~ "AD APOE high risk DEGs",
                                        `gene list` =="microglia-APOE-DEG-t-lower" ~ "AD APOE low risk DEGs",
                                        `gene list` =="microglia-polygenicHR-DEG-t-higher" ~ "AD polygenic high risk DEGs",
                                        `gene list` =="microglia-polygenicHR-DEG-t-lower" ~ "AD polygenic low risk DEGs",
                                        `gene list` == "microglia-treatment-eGenes-lfsr" ~ "eQTL eGenes"))

# extracting significant treatments
sign_treat = prioritization %>%
  dplyr::filter(p_bonferroni < 0.05) %>%
  dplyr::select(treatment) %>%
  distinct() %>%
  .$treatment
conditional = read_csv("../../../data/results/CELLECT/CELLECT-LDSC/results/conditional.csv") %>%
  dplyr::arrange(pvalue) %>%
  dplyr::filter(specificity_id %in% ids_of_interest & gwas %in% gwas_of_interest & annotation %in% sign_treat) %>%
  dplyr::group_by(gwas,specificity_id) %>% 
  dplyr::mutate(p_bonferroni_conditional = p.adjust(pvalue,method = "bonferroni")) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(gwas,specificity_id,annotation) %>% # summarise p_bonferroni across conditional annotations
  dplyr::summarise(p_bonferroni_mean = mean(p_bonferroni_conditional)) %>%
  dplyr::ungroup() %>%
  dplyr::rename(treatment = annotation , GWAS = gwas, `gene list` = specificity_id) %>% ## rename categories
  dplyr::mutate(treatment = factor(treatment, levels = c("untreated","IFN","LPS")),
                GWAS = factor(str_split_i(GWAS,"_",i=1),levels = str_split_i(gwas_of_interest,"_",i=1)),
                `gene list` = case_when(`gene list` == "microglia-cellex" ~ "treatment-specific genes",
                                        `gene list` =="microglia-PRS-DEG-t-higher" ~ "AD PRS high risk DEGs",
                                        `gene list` =="microglia-PRS-DEG-t-lower" ~ "AD PRS low risk DEGs",
                                        `gene list` =="microglia-APOE-DEG-t-higher" ~ "AD APOE high risk DEGs",
                                        `gene list` =="microglia-APOE-DEG-t-lower" ~ "AD APOE low risk DEGs",
                                        `gene list` =="microglia-polygenicHR-DEG-t-higher" ~ "AD polygenic high risk DEGs",
                                        `gene list` =="microglia-polygenicHR-DEG-t-lower" ~ "AD polygenic low risk DEGs",
                                        `gene list` == "microglia-treatment-eGenes-lfsr" ~ "eQTL eGenes"))
## add significance based on conditional analysis

  prioritization = prioritization %>%
     dplyr::left_join(conditional) %>%
    dplyr::mutate(p_val_sig_plot =  c("***", "**", "*", "")[findInterval(p_bonferroni_mean, c(0.001, 0.01, 0.05)) + 1])


## plot

p = ggplot(prioritization, aes(x=treatment, y=GWAS, fill = minus_log10_pval)) + 
  geom_raster() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), # lables vertical
        strip.text.y = element_blank()) +  #remove facet bar on y 
  scale_fill_gradient(low = "white", high = "#74121D",name = bquote(-log10(P^{s-LDSC}))) +
  geom_text(aes(label = p_val_sig_plot), col = "white", size = 5) + 
  ggtitle("Heritability enrichment - stratified LDSC") +
  theme_minimal() + 
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) +
  facet_grid(cols = vars(`gene list` )) +
  theme(legend.position="bottom")



pdf(paste0(outdir,"LDSC_tileplot_eQTL_cellex.pdf"), 
    height = 4,width = 4)
plot(p)
dev.off()
png(paste0(outdir,"LDSC_tileplot_eQTL_cellex.png"), 
    height = 4.5,width = 6.5, res = 400,units = "in")
plot(p)
dev.off()


####### now for PRS AD ##########
### interested in differences between treatments within PRS categories for AD only

ids_of_interest = c(
                    "microglia-polygenicHR-DEG-t-higher",
                    "microglia-APOE-DEG-t-higher", 
                    "microglia-PRS-DEG-t-higher"
                    )

gwas_of_interest = c( "AD_Bellenguez"   )

# to check
prioritization = read_csv("../../../data/results/CELLECT/CELLECT-LDSC/results/prioritization.csv") %>%
  dplyr::arrange(pvalue) %>%
  dplyr::filter(specificity_id %in% ids_of_interest & gwas %in% gwas_of_interest) %>%
  dplyr::group_by(specificity_id,annotation) %>% # grouping by gene lists tested and treatments because tests are not totally independent (there are shared genes)
  dplyr::mutate(p_bonferroni = p.adjust(pvalue,method = "bonferroni")) %>%
  dplyr::ungroup() 


conditional = read_csv("../../../data/results/CELLECT/CELLECT-LDSC/results/conditional.csv") %>%
  dplyr::arrange(pvalue) %>%
  dplyr::filter(specificity_id %in% ids_of_interest & gwas %in% gwas_of_interest) %>%
  dplyr::group_by(specificity_id,annotation) %>% # grouping by gene lists tested and treatments because tests are not totally independent (there are shared genes)
  dplyr::mutate(p_bonferroni_conditional = p.adjust(pvalue,method = "bonferroni")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
                p_fdr_conditional = p.adjust(pvalue,method = "fdr")) %>% # too few values to correctly estimate the distribution
  dplyr::rename(treatment = annotation , `conditioned on` = conditional_annotation 
                ,GWAS = gwas, `gene list` = specificity_id) %>% ## rename categories
  dplyr::mutate(treatment = factor(treatment, levels = c("untreated","IFN","LPS")),
                `conditioned on` = factor(`conditioned on`, levels = c("untreated","IFN","LPS")),
                minus_log10_pval = -log10(pvalue),
                GWAS = factor(str_split_i(GWAS,"_",i=1),levels = str_split_i(gwas_of_interest,"_",i=1)),
                `gene list` = case_when( `gene list` =="microglia-PRS-DEG-t-higher" ~ "AD PRS high risk DEGs",
                                         `gene list` =="microglia-PRS-DEG-t-lower" ~ "AD PRS low risk DEGs",
                                         `gene list` =="microglia-APOE-DEG-t-higher" ~ "AD APOE high risk DEGs",
                                         `gene list` =="microglia-APOE-DEG-t-lower" ~ "AD APOE low risk DEGs",
                                         `gene list` =="microglia-polygenicHR-DEG-t-higher" ~ "AD polygenic high risk DEGs",
                                         `gene list` =="microglia-polygenicHR-DEG-t-lower" ~ "AD polygenic low risk DEGs"
                ))
## add significance based on conditional analysis

conditional = conditional %>%
  dplyr::mutate(p_val_sig_plot =  c("***", "**", "*", "")[findInterval(p_bonferroni_conditional, c(0.001, 0.01, 0.05)) + 1])


## plot

p = ggplot(conditional, aes(x= `conditioned on` , y=treatment, fill = minus_log10_pval)) + 
  geom_raster() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), # lables vertical
        strip.text.y = element_blank()) +  #remove facet bar on y 
  scale_fill_gradient(low = "white", high = "#74121D",name = bquote(-log10(P^{s-LDSC}))) +
  geom_text(aes(label = p_val_sig_plot), col = "white", size = 5) + 
  ggtitle("Heritability enrichment - stratified LDSC") +
  xlab("Conditioned on") + 
  ylab("Treatment") +
  theme_minimal() + 
  theme(axis.ticks = element_blank(),
        panel.grid = element_blank()) +
  facet_grid(cols = vars(`gene list` )) +
  theme(legend.position="bottom")


# LPS is significant conditioned on untreated but not IFN - transcriptional similarities
# see https://elifesciences.org/articles/55851#content fig 4 sup 1
# pdf(paste0(outdir,"LDSC_tileplot_PRS.pdf"), 
#     height = 4,width = 4)
# plot(p)
# dev.off() # fuzzy
png(paste0(outdir,"LDSC_tileplot_PRS.png"), 
    height = 4.5,width = 8.5, res = 400,units = "in")
plot(p)
dev.off()

#### without directionality


ids_of_interest = c(
  "microglia-polygenicHR-DEG-t",
  "microglia-APOE-DEG-t", 
  "microglia-PRS-DEG-t"
)

gwas_of_interest = c( "AD_Bellenguez"   )

# to check
prioritization = read_csv("../../../data/results/CELLECT/CELLECT-LDSC/results/prioritization.csv") %>%
  dplyr::arrange(pvalue) %>%
  dplyr::filter(specificity_id %in% ids_of_interest & gwas %in% gwas_of_interest) %>%
  dplyr::group_by(specificity_id,annotation) %>% # grouping by gene lists tested and treatments because tests are not totally independent (there are shared genes)
  dplyr::mutate(p_bonferroni = p.adjust(pvalue,method = "bonferroni")) %>%
  dplyr::ungroup() 


conditional = read_csv("../../../data/results/CELLECT/CELLECT-LDSC/results/conditional.csv") %>%
  dplyr::arrange(pvalue) %>%
  dplyr::filter(specificity_id %in% ids_of_interest & gwas %in% gwas_of_interest) %>%
  dplyr::group_by(specificity_id,annotation) %>% # grouping by gene lists tested and treatments because tests are not totally independent (there are shared genes)
  dplyr::mutate(p_bonferroni_conditional = p.adjust(pvalue,method = "bonferroni")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    p_fdr_conditional = p.adjust(pvalue,method = "fdr")) %>% # too few values to correctly estimate the distribution
  dplyr::rename(treatment = annotation , `conditioned on` = conditional_annotation 
                ,GWAS = gwas, `gene list` = specificity_id) %>% ## rename categories
  dplyr::mutate(treatment = factor(treatment, levels = c("untreated","IFN","LPS")),
                `conditioned on` = factor(`conditioned on`, levels = c("untreated","IFN","LPS")),
                minus_log10_pval = -log10(pvalue),
                GWAS = factor(str_split_i(GWAS,"_",i=1),levels = str_split_i(gwas_of_interest,"_",i=1)),
                `gene list` = case_when( `gene list` =="microglia-PRS-DEG-t" ~ "AD PRS DEGs",
                                         `gene list` =="microglia-APOE-DEG-t" ~ "AD APOE DEGs",
                                         `gene list` =="microglia-polygenicHR-DEG-t" ~ "AD polygenic DEGs"
                ))
## add significance based on conditional analysis

conditional = conditional %>%
  dplyr::mutate(p_val_sig_plot =  c("***", "**", "*", "")[findInterval(p_bonferroni_conditional, c(0.001, 0.01, 0.05)) + 1])


## plot

p = ggplot(conditional, aes(x= `conditioned on` , y=treatment, fill = minus_log10_pval)) + 
  geom_raster() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), # lables vertical
        strip.text.y = element_blank()) +  #remove facet bar on y 
  scale_fill_gradient(low = "white", high = "#74121D",name = bquote(-log10(P^{s-LDSC}))) +
  geom_text(aes(label = p_val_sig_plot), col = "white", size = 5) + 
  ggtitle("Heritability enrichment - stratified LDSC") +
  xlab("Conditioned on") + 
  ylab("Treatment") +
  theme_minimal() + 
  theme(axis.ticks = element_blank(),
        panel.grid = element_blank()) +
  facet_grid(cols = vars(`gene list` )) +
  theme(legend.position="bottom")


# LPS is significant conditioned on untreated but not IFN - transcriptional similarities
# see https://elifesciences.org/articles/55851#content fig 4 sup 1
# pdf(paste0(outdir,"LDSC_tileplot_PRS_no_dir.pdf"), 
#     height = 4,width = 4)
# plot(p)
# dev.off() # fuzzy
png(paste0(outdir,"LDSC_tileplot_PRS_no_dir.png"), 
    height = 4.5,width = 8.5, res = 400,units = "in")
plot(p)
dev.off()
