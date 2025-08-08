.libPaths(c('/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/',"/software/teamtrynka/ma23/R4.1/libs",
            "/software/R-4.1.0/lib/R/library"))
library(patchwork)
library(tidyverse)
source("./functions.R")

input_dir = "../../data/results/1.alluvial_plots"
output_dir = "../../data/results/2.associations_preMAC_microglia_proliferation"
dir.create(output_dir, recursive = TRUE)
data = read.table(paste0(input_dir,
                         "/pools2-11_13-15_changing_props_iPSC_preMacs_microglia_WGS_sc.txt"),
                  header = TRUE)

# read in proportion error estimates
error_estimates = read.table("../../../OTAR2065_phenotypic_QTLs/data/w/error_approximations_generic_pool.txt", header = TRUE) %>%
  dplyr::filter(coverage ==5 & genotype =="old")

# Perform the inner join based on the closest proportion values
data_with_error =  data %>%
  dplyr::rowwise() %>%
  dplyr::mutate(closest_error_proportion = find_closest_value(prop, error_estimates$w_real)) %>%
  dplyr::ungroup() %>%
  dplyr::inner_join(error_estimates, by = c("closest_error_proportion" = "w_real")) %>%
  dplyr::mutate(prop_adjusted_mean = dplyr::case_when(prop <0.1 ~ prop - mean_wdif,
                                                      prop >=0.1 ~ prop + mean_wdif),
                prop_adjusted_median = dplyr::case_when(prop <0.1 ~ prop - median_wdif,
                                                        prop >=0.1 ~ prop + median_wdif),
                prop_adjusted_max = dplyr::case_when(prop <0.1 ~ prop - max_wdif,
                                                     prop >=0.1 ~ prop + max_wdif))
# if any estimate becomes negative, make it zero:
data_with_error[data_with_error$prop_adjusted_mean<0, "prop_adjusted_mean"] = 0
data_with_error[data_with_error$prop_adjusted_median<0, "prop_adjusted_median"] = 0
data_with_error[data_with_error$prop_adjusted_max<0, "prop_adjusted_max"] = 0

# need to distinguish between those failed at point of preMACS and at point of microglia
# If Inf, then they were 0 at preMAC
# If 0, they disappear at microglia
efficiency = list()
for(p in paste0("pool",3:15)){
  if(p == "pool3"){
    efficiency[[p]] = data_with_error %>%
      dplyr::filter(stage!="iPSC" & pool=="pool3") %>%
      tidyr::pivot_wider(names_from = sample,values_from = prop_adjusted_mean,id_cols = Line) %>%
      dplyr::reframe(line=Line,
                     sc_untreated = log1p(h3_diff2_untreated) / log1p(d35_preMAC),
                     sc_IFN =  log1p(h3_diff2_IFNg) / log1p(d35_preMAC),
                     sc_LPS =  log1p(h3_diff2_LPS) / log1p(d35_preMAC),
                     phago_untreated =  log1p(`phago_s3D-1`) / log1p(d57_preMAC),
                     phago_IFN =  log1p(`phago_s3E-1`) / log1p(d57_preMAC),
                     phago_LPS =  log1p(`phago_s3F-1`) / log1p(d57_preMAC)) %>%
      # Puigdevall et al prolif. rate: not much different distribution from raw props.
      
      tidyr::pivot_longer(cols =-line,names_to = "sample",values_to = "scaled_proportion") %>%
      dplyr::mutate(pool="pool3")
    # don't have right preMAC for migration
    # don't have right preMACs for pool2
  }
  if(p == "pool4"){
    efficiency[[p]] = data_with_error %>%
      dplyr::filter(stage!="iPSC" & pool=="pool4") %>%
      tidyr::pivot_wider(names_from = sample,values_from = prop,id_cols = Line) %>%
      dplyr::reframe(line=Line,
                     sc_untreated = log1p(h4_diff2_untreated) / log1p(d35_preMAC),
                     sc_IFN =  log1p(h4_diff2_IFNg) / log1p(d35_preMAC),
                     sc_LPS =  log1p(h4_diff2_LPS) / log1p(d35_preMAC),
                     phago_untreated =  log1p(`phago_s4D-1`) / log1p(d57_preMAC),
                     phago_IFN =  log1p(`phago_s4E-1`) / log1p(d57_preMAC),
                     phago_LPS =  log1p(`phago_s4F-1`) / log1p(d57_preMAC)) %>%
      tidyr::pivot_longer(cols =-line,names_to = "sample",values_to = "scaled_proportion") %>%
      dplyr::mutate(pool="pool4")
    # don't have right preMAC for migration
    
  }
  if(p == "pool5"){
    efficiency[[p]] = data_with_error %>%
      dplyr::filter(stage!="iPSC" & pool=="pool5") %>%
      tidyr::pivot_wider(names_from = sample,values_from = prop_adjusted_mean,id_cols = Line) %>%
      dplyr::reframe(line=Line,
                     sc_untreated = log1p(P5_diff2_untreated) / log1p(d35_preMAC),
                     sc_IFN =  log1p(P5_diff2_IFN) / log1p(d35_preMAC),
                     sc_LPS =  log1p(P5_diff2_LPS) / log1p(d35_preMAC)) %>%
      tidyr::pivot_longer(cols =-line,names_to = "sample",values_to = "scaled_proportion") %>%
      dplyr::mutate(pool="pool5")
    # don't have seeded cells for migration
    # don't have right preMAC for phagocytosis
    
  }
  if(p == "pool6"){
    efficiency[[p]] = data_with_error %>%
      dplyr::filter(stage!="iPSC" & pool=="pool6") %>%
      tidyr::pivot_wider(names_from = sample,values_from = prop_adjusted_mean,id_cols = Line) %>%
      dplyr::reframe(line=Line,
                     sc_untreated = log1p(P6_diff2_untreated) / log1p(d35_preMAC),
                     sc_IFN =  log1p(P6_diff2_IFN) / log1p(d35_preMAC),
                     sc_LPS =  log1p(P6_diff2_LPS) / log1p(d35_preMAC)) %>%
      tidyr::pivot_longer(cols =-line,names_to = "sample",values_to = "scaled_proportion") %>%
      dplyr::mutate(pool="pool6")
    # don't have seeded cells for migration
    # don't have right preMAC for phagocytosis
    
  }
  if(p == "pool7"){
    efficiency[[p]] = data_with_error %>%
      dplyr::filter(stage!="iPSC" & pool=="pool7") %>%
      tidyr::pivot_wider(names_from = sample,values_from = prop_adjusted_mean,id_cols = Line) %>%
      dplyr::reframe(line=Line,
                     sc_untreated = log1p(P7_diff2_untreated) / log1p(D36_preMAC),
                     sc_IFN =  log1p(P7_diff2_IFN) / log1p(D36_preMAC),
                     sc_LPS =  log1p(P7_diff2_LPS) / log1p(D36_preMAC),
                     migr_untreated =  log1p(Mig17_Untr) / log1p(D47_preMAC),
                     migr_LPS =  log1p(Mig17_LPS) / log1p(D47_preMAC),
                     phago_untreated =  log1p(`phago_s7A-1`) / log1p(D54_preMAC),
                     phago_IFN =  log1p(`phago_s7B-1`) / log1p(D54_preMAC),
                     phago_LPS =  log1p(`phago_s7C-1`) / log1p(D54_preMAC)) %>%
      tidyr::pivot_longer(cols =-line,names_to = "sample",values_to = "scaled_proportion") %>%
      dplyr::mutate(pool="pool7")
  }
  if(p == "pool8"){
    efficiency[[p]] = data_with_error %>%
      dplyr::filter(stage!="iPSC" & pool=="pool8") %>%
      tidyr::pivot_wider(names_from = sample,values_from = prop_adjusted_mean,id_cols = Line) %>%
      dplyr::reframe(line=Line,
                     sc_untreated = log1p(P8_diff2_untreated) / log1p(D36_preMAC),
                     sc_IFN =  log1p(P8_diff2_IFN) / log1p(D36_preMAC),
                     sc_LPS =  log1p(P8_diff2_LPS) / log1p(D36_preMAC),
                     migr_untreated =  log1p(Mig16_Untr) / log1p(D40_preMAC),
                     migr_LPS =  log1p(Mig16_LPS) / log1p(D40_preMAC),
                     phago_untreated =  log1p(`phago_s8A-1`) / log1p(D54_preMAC),
                     phago_IFN =  log1p(`phago_s8A-1`) / log1p(D54_preMAC),
                     phago_LPS =  log1p(`phago_s8A-1`) / log1p(D54_preMAC)) %>%
      tidyr::pivot_longer(cols =-line,names_to = "sample",values_to = "scaled_proportion") %>%
      dplyr::mutate(pool="pool8")
  }
  if(p == "pool9"){
    efficiency[[p]] = data_with_error %>%
      dplyr::filter(stage!="iPSC" & pool=="pool9") %>%
      tidyr::pivot_wider(names_from = sample,values_from = prop_adjusted_mean,id_cols = Line) %>%
      dplyr::reframe(line=Line,
                     sc_untreated_A = log1p(P9_Diff2_Untreated_A) / log1p(D35_PreMac),
                     sc_untreated_B = log1p(P9_Diff2_Untreated_B) / log1p(D35_PreMac),
                     sc_untreated_C = log1p(P9_Diff2_Untreated_B) / log1p(D35_PreMac),
                     sc_untreated_D = log1p(P9_Diff2_Untreated_B) / log1p(D35_PreMac),
                     sc_IFN_A =  log1p(P9_Diff2_IFN_A) / log1p(D35_PreMac),
                     sc_IFN_B =  log1p(P9_Diff2_IFN_B) / log1p(D35_PreMac),
                     sc_IFN_C =  log1p(P9_Diff2_IFN_C) / log1p(D35_PreMac),
                     sc_LPS_A =  log1p(P9_Diff2_LPS_A) / log1p(D35_PreMac),
                     sc_LPS_B =  log1p(P9_Diff2_LPS_B) / log1p(D35_PreMac),
                     sc_LPS_C =  log1p(P9_Diff2_LPS_C) / log1p(D35_PreMac),
                     migr_untreated =  log1p(Untr_cells) / log1p(D39_PreMac),
                     migr_LPS =  log1p(LPS_cells) / log1p(D39_PreMac),
                     phago_untreated =  log1p(`phago_s9A-1`) / log1p(D50_PreMac),
                     phago_LPS =  log1p(`phago_s9C-1`) / log1p(D50_PreMac)) %>%
      tidyr::pivot_longer(cols =-line,names_to = "sample",values_to = "scaled_proportion") %>%
      dplyr::mutate(pool="pool9")
  }
  if(p == "pool10"){
    efficiency[[p]] = data_with_error %>%
      dplyr::filter(stage!="iPSC" & pool=="pool10") %>%
      tidyr::pivot_wider(names_from = sample,values_from = prop_adjusted_mean,id_cols = Line) %>%
      dplyr::reframe(line=Line,
                     sc_untreated_A = log1p(P10_Diff2_Untreated_A) / log1p(D36_PreMac),
                     sc_untreated_B = log1p(P10_Diff2_Untreated_B) / log1p(D36_PreMac),
                     sc_untreated_C = log1p(P10_Diff2_Untreated_B) / log1p(D36_PreMac),
                     sc_untreated_D = log1p(P10_Diff2_Untreated_B) / log1p(D36_PreMac),
                     sc_untreated_old_A = log1p(P10_Diff5_Untreated_A) / log1p(D54_PreMac),
                     sc_untreated_old_B = log1p(P10_Diff5_Untreated_B) / log1p(D54_PreMac),
                     sc_IFN_A =  log1p(P10_Diff2_IFN_A) / log1p(D36_PreMac),
                     sc_IFN_B =  log1p(P10_Diff2_IFN_B) / log1p(D36_PreMac),
                     sc_IFN_C =  log1p(P10_Diff2_IFN_C) / log1p(D36_PreMac),
                     sc_IFN_old = log1p(P10_Diff5_IFN) / log1p(D54_PreMac),
                     sc_LPS_A =  log1p(P10_Diff2_LPS_A) / log1p(D36_PreMac),
                     sc_LPS_B =  log1p(P10_Diff2_LPS_B) / log1p(D36_PreMac),
                     sc_LPS_C =  log1p(P10_Diff2_LPS_C) / log1p(D36_PreMac),
                     sc_LPS_old = log1p(P10_Diff5_LPS) / log1p(D54_PreMac),
                     migr_untreated =  log1p(Untr_cell) / log1p(D43_PreMac),
                     migr_LPS =  log1p(LPS_cells) / log1p(D43_PreMac),
                     phago_untreated =  log1p(`phago_s10A-1`) / log1p(D50_PreMac),
                     phago_IFN =  log1p(`phago_s10B-1`) / log1p(D50_PreMac),
                     phago_LPS =  log1p(`phago_s10C-1`) / log1p(D50_PreMac)) %>%
      tidyr::pivot_longer(cols =-line,names_to = "sample",values_to = "scaled_proportion") %>%
      dplyr::mutate(pool="pool10")
  }
  if(p == "pool11"){
    efficiency[[p]] = data_with_error %>%
      dplyr::filter(stage!="iPSC" & pool=="pool11") %>%
      tidyr::pivot_wider(names_from = sample,values_from = prop_adjusted_mean,id_cols = Line) %>%
      dplyr::reframe(line=Line,
                     sc_untreated_A = log1p(P11_Diff2_Untreated_A) / log1p(day35_preMAC),
                     sc_untreated_B = log1p(P11_Diff2_Untreated_B) / log1p(day35_preMAC),
                     sc_IFN =  log1p(P11_Diff2_IFN) / log1p(day35_preMAC),
                     sc_LPS =  log1p(P11_Diff2_LPS) / log1p(day35_preMAC),
                     migr_untreated =  log1p(pool11_UNTR) / log1p(day39_preMAC),
                     migr_LPS =  log1p(pool11_LPS) / log1p(day39_preMAC),
                     phago_untreated =  log1p(pool11_A1) / log1p(day49_preMAC),
                     phago_IFN =  log1p(pool11_B1) / log1p(day49_preMAC),
                     phago_LPS =  log1p(pool11_C1) / log1p(day49_preMAC)) %>%
      tidyr::pivot_longer(cols =-line,names_to = "sample",values_to = "scaled_proportion") %>%
      dplyr::mutate(pool="pool11")
  }
  if(p == "pool13"){
    efficiency[[p]] = data_with_error %>%
      dplyr::filter(stage!="iPSC" & pool=="pool13") %>%
      tidyr::pivot_wider(names_from = sample,values_from = prop_adjusted_mean,id_cols = Line) %>%
      dplyr::reframe(line=Line,
                     sc_untreated_A = log1p(P13_Diff2_Untreated_A) / log1p(day36_preMAC),
                     sc_untreated_B = log1p(P13_Diff2_Untreated_B) / log1p(day36_preMAC),
                     sc_IFN =  log1p(P13_Diff2_IFN) / log1p(day36_preMAC),
                     sc_LPS =  log1p(P13_Diff2_LPS) / log1p(day36_preMAC),
                     migr_untreated =  log1p(pool13_UNTR) / log1p(day43_preMAC),
                     migr_IFN =  log1p(pool13_IFN) / log1p(day43_preMAC),
                     migr_LPS =  log1p(pool13_LPS) / log1p(day43_preMAC),
                     phago_untreated =  log1p(pool13_A1) / log1p(day50_preMAC),
                     phago_IFN =  log1p(pool13_B1) / log1p(day50_preMAC),
                     phago_LPS =  log1p(pool13_C1) / log1p(day50_preMAC)) %>%
      tidyr::pivot_longer(cols =-line,names_to = "sample",values_to = "scaled_proportion") %>%
      dplyr::mutate(pool="pool13")
  }
  if(p == "pool14"){
    efficiency[[p]] = data_with_error %>%
      dplyr::filter(stage!="iPSC" & pool=="pool14") %>%
      tidyr::pivot_wider(names_from = sample,values_from = prop_adjusted_mean,id_cols = Line) %>%
      dplyr::reframe(line=Line,
                     sc_untreated_A = log1p(P14_Diff2_Untreated_A) / log1p(D36_PreMac),
                     sc_untreated_B = log1p(P14_Diff2_Untreated_B) / log1p(D36_PreMac),
                     sc_IFN_A =  log1p(P14_Diff2_IFN_A) / log1p(D36_PreMac),
                     sc_IFN_B =  log1p(P14_Diff2_IFN_B) / log1p(D36_PreMac),
                     sc_LPS_A =  log1p(P14_Diff2_LPS_A) / log1p(D36_PreMac),
                     sc_LPS_B =  log1p(P14_Diff2_LPS_B) / log1p(D36_PreMac),
                     migr_untreated =  log1p(Mig24_Untr) / log1p(D39_PreMac),
                     migr_LPS =  log1p(Mig24_LPS) / log1p(D39_PreMac)) %>%
      tidyr::pivot_longer(cols =-line,names_to = "sample",values_to = "scaled_proportion") %>%
      dplyr::mutate(pool="pool14")
  }
  if(p == "pool15"){
    efficiency[[p]] = data_with_error %>%
      dplyr::filter(stage!="iPSC" & pool=="pool15") %>%
      tidyr::pivot_wider(names_from = sample,values_from = prop_adjusted_mean,id_cols = Line) %>%
      dplyr::reframe(line=Line,
                     sc_untreated_A = log1p(P15_Diff2_Untreated_A) / log1p(D36_PreMac),
                     sc_untreated_B = log1p(P15_Diff2_Untreated_B) / log1p(D36_PreMac),
                     sc_IFN_A =  log1p(P15_Diff2_IFN_A) / log1p(D36_PreMac),
                     sc_IFN_B =  log1p(P15_Diff2_IFN_B) / log1p(D36_PreMac),
                     sc_LPS_A =  log1p(P15_Diff2_LPS_A) / log1p(D36_PreMac),
                     sc_LPS_B =  log1p(P15_Diff2_LPS_B) / log1p(D36_PreMac),
                     migr_untreated =  log1p(Mig23_Untr) / log1p(D39_PreMac),
                     migr_IFN =  log1p(Mig23_IFN) / log1p(D39_PreMac),
                     migr_LPS =  log1p(Mig23_LPS) / log1p(D39_PreMac)) %>%
      tidyr::pivot_longer(cols =-line,names_to = "sample",values_to = "scaled_proportion") %>%
      dplyr::mutate(pool="pool15")
  }
  
  
}
efficiency = do.call("rbind",efficiency)
summary(efficiency)
efficiency = efficiency %>%
  dplyr::filter(!scaled_proportion %in% c(Inf,NA,NaN)) %>%
  dplyr::filter(line!="sh5y5y") %>% # not interested in phagocytosis line
  dplyr::mutate(sequencing = case_when(stringr::str_detect(sample,"sc") ~ "single_cell",
                                       stringr::str_detect(sample,"phago") ~ "WGS",
                                       stringr::str_detect(sample,"migr") ~ "WGS"),
                treatment = case_when(stringr::str_detect(sample,"untreated") ~ "untreated",
                                      stringr::str_detect(sample,"IFN") ~ "IFN",
                                      stringr::str_detect(sample,"LPS") ~ "LPS")) 
summary(efficiency$scaled_proportion)
hist(efficiency$scaled_proportion, breaks = 100)
# where do the NAs come from?
# this is not normal at all
hist(log10(efficiency$scaled_proportion),breaks = 40)
# still very extreme values
# model without transforming via GLMM?

# diagnostic plots

# Run a diagnostic lot of the group variances vs group
# means (genotype x nutrient x clipping grouping).  Code
# used to produce the plot :
# https://github.com/QCBSRworkshops/workshop07/blob/main/pres-fr/data/glmm_e.r

# without logging

summary_tibble = efficiency %>%
  dplyr::group_by(line,treatment,pool) %>%
  dplyr::summarise(means =mean(scaled_proportion),
                   vars = var(scaled_proportion))
# dplyr::filter(!log_means %in% c(-Inf,Inf,NA,NaN))
vars = summary_tibble$vars
means = summary_tibble$means

# Quasi-Poisson
lm1 = lm(vars ~ means - 1)
phi.fit = coef(lm1)
# The -1 specifies a model with the intercept set to zero

# Negative binomial
lm2 = lm(vars ~ I(means^2) + offset(means) - 1)
k.fit = 1/coef(lm2)
# The offset() is used to specify that we want the group
# means added as a term with its coefficient fixed to 1

# Non-parametric loess fit
Lfit = loess(vars ~ means)

# The plot
plot(vars ~ means, xlab = "Group means", ylab = "Group variances ")
abline(a = 0, b = 1, lty = 2)
text(-4, -5, "Poisson")
curve(phi.fit * x, col = 2, add = TRUE)
# bquote() is used to substitute numeric values in
# equations with symbols
text(4, 4, bquote(paste("QP: ", sigma^2 == .(round(phi.fit,
                                                   1)) * mu)), col = 2)
curve(x * (1 + x/k.fit), col = 4, add = TRUE)
text(5, 1, paste("NB: k = ", round(k.fit, 1), sep = ""),
     col = 4)
mvec <- 0:120
lines(mvec, predict(Lfit, mvec), col = 5)
text(5, 6, "loess", col = 5)

# logging
summary_tibble = efficiency %>%
  dplyr::group_by(line,treatment,pool) %>%
  dplyr::summarise(log_means = log(mean(scaled_proportion)),
                   log_vars = log(var(scaled_proportion))) %>%
  dplyr::filter(!log_means %in% c(-Inf,Inf,NA,NaN))
log_vars = summary_tibble$log_vars
log_means = summary_tibble$log_means

# normal
lm0 = lm(log_vars ~ log_means)
fit = coef(lm0)

# Quasi-Poisson
lm1 = lm(log_vars ~ log_means - 1)
phi.fit = coef(lm1)
# The -1 specifies a model with the intercept set to zero

# Negative binomial
lm2 = lm(log_vars ~ I(log_means^2) + offset(log_means) - 1)
k.fit = 1/coef(lm2)
# The offset() is used to specify that we want the group
# means added as a term with its coefficient fixed to 1

# Non-parametric loess fit
Lfit = loess(log_vars ~ log_means)

# The plot
plot(log_vars ~ log_means, xlab = "Group means (log)", ylab = "Group variances (log)")
abline(a = 0, b = 1, lty = 2)
text(-4, -5, "Poisson")

abline(fit) 
text(-4, -11, "Normal")

curve(phi.fit * x, col = 2, add = TRUE)
# bquote() is used to substitute numeric values in
# equations with symbols
text(4, 4, bquote(paste("QP: ", sigma^2 == .(round(phi.fit,
                                                   1)) * mu)), col = 2)
curve(x * (1 + x/k.fit), col = 4, add = TRUE)
text(5, 1, paste("NB: k = ", round(k.fit, 1), sep = ""),
     col = 4)
mvec <- 0:120
lines(mvec, predict(Lfit, mvec), col = 5)
text(5, 6, "loess", col = 5)

# Quasi-Poisson, normal or loess?
#https://r.qcbs.ca/workshop07/book-en/choose-an-error-distribution.html



####### association tests
line_info=read.csv("../../../OTAR2065_phenotypic_QTLs/data/allinfo_hipsci_PD_AD_PRS_nondisease_feederfree_european.csv") %>%
  dplyr::mutate(sex = case_when(grepl(pattern = "Female",.$Sex) ~ "Female",
                                grepl(pattern = "Male",.$Sex) ~ "Male")) %>%
  dplyr::select(Line,sex) %>%
  dplyr::rename(line=Line)

efficiency = efficiency %>%
  dplyr::left_join(.,line_info)
# load minor allele dosage file for all lines

genotype = readr::read_csv("../../../OTAR2065_phenotypic_QTLs/data/genotypes/full_genotype/genotype_minor_allele_dosage_1.csv") %>%
  dplyr::relocate(rn) %>%
  dplyr::ungroup() %>%
  #dplyr::filter(rn %in% to_subset) %>%
  tidyr::pivot_longer(cols = !rn, values_to = "alt_allele_dosage", names_to = "line") %>%
  dplyr::mutate(alt_allele_dosage = dplyr::case_when(
    alt_allele_dosage == 0.5 ~ 1,
    alt_allele_dosage == 1   ~ 2,
    TRUE          ~ alt_allele_dosage  # Recode so var of allele dosage is 1, Keep 0 unchanged
  ))

length(unique(genotype$rn)) # number of SNPs after filtering

# remove genotypes under 10% MAF
snps_to_keep = genotype %>%
  dplyr::mutate(alt_allele_dosage = as.factor(alt_allele_dosage)) %>%
  count(rn,alt_allele_dosage) %>%
  dplyr::group_by(rn) %>%
  dplyr::filter(min(n) >= (0.1 * length(unique(genotype$line)))) %>%
  ungroup() %>%
  dplyr::select(rn) %>%
  distinct() %>%
  unlist()


genotype_subset = genotype %>%
  dplyr::filter(rn %in% snps_to_keep) %>%
  dplyr::mutate(rn = str_replace(rn, "chr", "")) 

length(unique(genotype_subset$rn)) # number of SNPs after filtering

if(nrow(genotype_subset_eqtl)!=0){
  
  message("Adding genotype PCs")
  
  full_genotype_pcs = read.table("../../../OTAR2065_sc_eQTL/data/genotype/plink_genotypes/all_pools.genotype.MAF05.eigenvec")
  rownames(full_genotype_pcs) = ifelse(full_genotype_pcs$V1 ==full_genotype_pcs$V2,
                                       yes = full_genotype_pcs$V1,
                                       no = paste(full_genotype_pcs$V1,full_genotype_pcs$V2,sep = "_"))
  full_genotype_pcs = full_genotype_pcs[c(2:7)] %>%
    dplyr::rename(line=V2,genotypePC1 = V3,genotypePC2 = V4,genotypePC3 = V5,genotypePC4 = V6,genotypePC5 = V7)
  
  efficiency_w_info = efficiency %>%
    dplyr::inner_join(.,genotype_subset) %>%
    dplyr::distinct() %>%
    dplyr::inner_join(.,full_genotype_pcs) %>%
    dplyr::distinct() 
  
  rm(genotype,full_genotype_pcs)
  gc()
  
  # cleaning up Inf and NAs
  efficiency_w_info = efficiency_w_info %>%
    dplyr::mutate(log_scaled_proportion = log(scaled_proportion)) %>%
    dplyr::filter(!log_scaled_proportion %in% c(-Inf,Inf,NA,NaN)) %>%
    distinct()
  # List to store the results
  res = list()
  result <- microbenchmark(
    for(i in 1:length(unique(genotype_subset$rn))){
      snp =  unique(genotype_subset$rn)[i]
      
      test = efficiency_w_info[efficiency_w_info$rn == snp, ] %>%
        dplyr::filter(!is.na(sex)) # remove missing data causing problems
      
      
      # Poisson GLMM Given the mean-variance relationship, we
      # will most likely need a model with over-dispersion.
      my_fit =  lmerTest::lmer(log_scaled_proportion ~ sex  +  genotypePC1 + genotypePC2 + genotypePC3 + genotypePC4 + genotypePC5 + sequencing + alt_allele_dosage*treatment + (1|line) + (1|pool),
                               data = test, 
                               control = lmerControl(optCtrl=list(maxfun=5000) ))
      
      sum_res = summary(my_fit)
      res[[i]] = list(
        coefficients = sum_res$coefficients,
        varcor =  sum_res$varcor,
        residuals =  sum_res$residuals,
        call =  sum_res$call,
        AICtab =  sum_res$AICtab)
      
      if(i%% 1000 ==0)  { print(paste0("i is ",i))}
    }
    
    ,times=1)
  print(result)
}

# save these as rds
saveRDS(object = res,paste0(output_dir,"/res.rds"))
# incorporate maximum proportion information in each pair, to then fit as covariate
# also associations with preMAC aging - timecourse model