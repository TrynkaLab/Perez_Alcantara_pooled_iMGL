# checking rare variant burden associations
.libPaths(c('/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/',"/software/teamtrynka/ma23/R4.1/libs",
            "/software/R-4.1.0/lib/R/library"))
library(patchwork)
library(tidyverse)
library(lmerTest)
source("./functions.R")

dir = "../../data/results/2.1.rare_associations_preMAC_microglia_proliferation"

res = readRDS(paste0(dir,"/res.rds"))
res_dichotomy = readRDS(paste0(dir,"/res_dichotomy.rds"))

# summary deletereous burden
# continuous measure
p_burden = list()
coef_burden = list()
se_burden = list()
p_IFN_interaction = list()
coef_IFN_interaction = list()
se_IFN_interaction = list()
p_LPS_interaction = list()
coef_LPS_interaction = list()
se_LPS_interaction = list()
gene = list()
for(i in 1:length(res$del)){
  p_burden[[i]] = res$del[[i]]$coefficients["rare_burden","Pr(>|t|)"] # pval
  coef_burden[[i]] = res$del[[i]]$coefficients["rare_burden","Estimate"] # coefficient
  se_burden[[i]] = res$del[[i]]$coefficients["rare_burden","Std. Error"] # std error
  gene[[i]] = res$del[[i]]$gene
  p_IFN_interaction[[i]] = res$del[[i]]$coefficients["treatmentIFN:rare_burden","Pr(>|t|)"] # pval
  coef_IFN_interaction[[i]] = res$del[[i]]$coefficients["treatmentIFN:rare_burden","Estimate"] # coefficient
  se_IFN_interaction[[i]] = res$del[[i]]$coefficients["treatmentIFN:rare_burden","Std. Error"] 
  
  p_LPS_interaction[[i]] = res$del[[i]]$coefficients["treatmentLPS:rare_burden","Pr(>|t|)"] # pval
  coef_LPS_interaction[[i]] = res$del[[i]]$coefficients["treatmentLPS:rare_burden","Estimate"] # coefficient
  se_LPS_interaction[[i]] = res$del[[i]]$coefficients["treatmentLPS:rare_burden","Std. Error"] 
  
}

p_burden = unlist(p_burden)
coef_burden = unlist(coef_burden)
se_burden = unlist(se_burden)
p_IFN_interaction = unlist(p_IFN_interaction)
coef_IFN_interaction = unlist(coef_IFN_interaction)
se_IFN_interaction = unlist(se_IFN_interaction)
p_LPS_interaction = unlist(p_LPS_interaction)
coef_LPS_interaction = unlist(coef_LPS_interaction)
se_LPS_interaction = unlist(se_LPS_interaction)
gene = unlist(gene)

res_df_del = data.frame(coef_burden=coef_burden,p_burden=p_burden,se_burden = se_burden,
                        coef_IFN_interaction = coef_IFN_interaction,p_IFN_interaction = p_IFN_interaction,se_IFN_interaction = se_IFN_interaction,
                    coef_LPS_interaction = coef_LPS_interaction,p_LPS_interaction = p_LPS_interaction,se_LPS_interaction = se_LPS_interaction, 
                    gene = gene) %>%
  dplyr::relocate(gene) %>%
  distinct()


hist(res_df_del$p_burden)

# dichotomised

p_burden = list()
coef_burden = list()
se_burden = list()
p_IFN_interaction = list()
coef_IFN_interaction = list()
se_IFN_interaction = list()
p_LPS_interaction = list()
coef_LPS_interaction = list()
se_LPS_interaction = list()
gene = list()

for(i in 1:length(res_dichotomy$del)){
  p_burden[[i]] = res_dichotomy$del[[i]]$coefficients["rare_burden_dichotomyyes","Pr(>|t|)"] # pval
  coef_burden[[i]] = res_dichotomy$del[[i]]$coefficients["rare_burden_dichotomyyes","Estimate"] # coefficient
  se_burden[[i]] = res_dichotomy$del[[i]]$coefficients["rare_burden_dichotomyyes","Std. Error"] # std error
  gene[[i]] = res_dichotomy$del[[i]]$gene
  
  p_IFN_interaction[[i]] = res_dichotomy$del[[i]]$coefficients["treatmentIFN:rare_burden_dichotomyyes","Pr(>|t|)"] # pval
  coef_IFN_interaction[[i]] = res_dichotomy$del[[i]]$coefficients["treatmentIFN:rare_burden_dichotomyyes","Estimate"] # coefficient
  se_IFN_interaction[[i]] = res_dichotomy$del[[i]]$coefficients["treatmentIFN:rare_burden_dichotomyyes","Std. Error"] 
  
  p_LPS_interaction[[i]] = res_dichotomy$del[[i]]$coefficients["treatmentLPS:rare_burden_dichotomyyes","Pr(>|t|)"] # pval
  coef_LPS_interaction[[i]] = res_dichotomy$del[[i]]$coefficients["treatmentLPS:rare_burden_dichotomyyes","Estimate"] # coefficient
  se_LPS_interaction[[i]] = res_dichotomy$del[[i]]$coefficients["treatmentLPS:rare_burden_dichotomyyes","Std. Error"] 
  
}

p_burden = unlist(p_burden)
coef_burden = unlist(coef_burden)
se_burden = unlist(se_burden)
p_IFN_interaction = unlist(p_IFN_interaction)
coef_IFN_interaction = unlist(coef_IFN_interaction)
se_IFN_interaction = unlist(se_IFN_interaction)
p_LPS_interaction = unlist(p_LPS_interaction)
coef_LPS_interaction = unlist(coef_LPS_interaction)
se_LPS_interaction = unlist(se_LPS_interaction)
gene = unlist(gene)

res_df_dichotomised_del = data.frame(coef_burden=coef_burden,p_burden=p_burden,
                                     p_adjust_burden = p.adjust(p_burden,method = "fdr"),
                                     se_burden = se_burden,
                        coef_IFN_interaction = coef_IFN_interaction, p_IFN_interaction = p_IFN_interaction,
                        p_adjust_IFN_interaction = p.adjust(p_IFN_interaction,method = "fdr"),
                        se_IFN_interaction = se_IFN_interaction,
                        coef_LPS_interaction = coef_LPS_interaction,p_LPS_interaction = p_LPS_interaction,
                        p_adjust_LPS_interaction = p.adjust(p_LPS_interaction,method = "fdr"),
                        se_LPS_interaction = se_LPS_interaction, 
                        gene = gene) %>%
  dplyr::relocate(gene) %>%
  distinct()

hist(res_df_dichotomised_del$p_burden)



### checking residuals for model fit
ggcheck_the_model(res$del)

# add figures like in https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/adding-covariates-to-a-linear-model
# to check covariate need.