# checking line proportions of migration assay
# Pools of 19 lines, migration assay with different conditions
# with media only (media_only),  with C5a (chemoattractant) (3nM_C5a), with C5a and C5a inhibitor (PMX-53 ) (3nM_C5a_PMX53), with media and PMX5 (media_PMX53)
# double-check these
# 6 well and 12 well plate
# with and without IFN
# main objective of this essay is 1) to see changes in whole population migration in different conditions 
# and 2) how the line proportion change using WGS for all the samples
# 
# Note that sample 81B and 91B went through extra cycles of PCR due to low yield from the initial library prep.

.libPaths(c( "/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/",
             "/software/teamtrynka/conda/otar2065/lib/R/library",
             "/software/teamtrynka/ma23/R4.1/libs"))

library(patchwork)
library(tidyverse)
library(plotrix)
source("./functions.R")

dir.create(file.path("../../data/results/1.check_line_proportions"),recursive = TRUE)
sample_info=read.csv("../../data/all_pools_migration_sample_info.csv")
line_info=read.csv("../../data/allinfo_hipsci_PD_AD_PRS_nondisease_feederfree_european.csv") %>%
  dplyr::mutate(sex = case_when(grepl(pattern = "Female",.$Sex) ~ "Female",
                                grepl(pattern = "Male",.$Sex) ~ "Male")) %>%
  dplyr::select(Line,sex) %>%
  dplyr::rename(line=Line)
w_estimates = read_w_estimates(sample_info,
                               "../../data/w/",
                               assay = "migration",
                               pools = paste0("pool", c(2:11,13)))
# convert slightly negative values to 0
w_estimates[w_estimates$proportion<0, "proportion"] = 0
w_estimates = w_estimates %>%
  dplyr::rename(line=donor)
write.table(w_estimates,"../../data/results/1.check_line_proportions/w_estimates.txt",
            row.names = FALSE,quote = FALSE)


# read in proportion error estimates
error_estimates = read.table("../../data/w/error_approximations_generic_pool.txt", header = TRUE) %>%
  dplyr::filter(coverage ==5 & genotype =="old")

# Perform the inner join based on the closest proportion values
w_estimates =  w_estimates %>%
  dplyr::rowwise() %>%
  dplyr::mutate(closest_error_proportion = find_closest_value(proportion, error_estimates$w_real)) %>%
  dplyr::ungroup() %>%
  dplyr::inner_join(error_estimates, by = c("closest_error_proportion" = "w_real"))

# Errors:
## If the proportion is under 0.1, it will likely be overestimated - substract the error
## If the proportion is over 0.1, it will likely be underestimated - add the error

w_estimates_adjusted = w_estimates %>%
  dplyr::mutate(prop_adjusted_mean = dplyr::case_when(proportion <0.1 ~ proportion - mean_wdif,
                                                      proportion >=0.1 ~ proportion + mean_wdif),
                prop_adjusted_median = dplyr::case_when(proportion <0.1 ~ proportion - median_wdif,
                                                        proportion >=0.1 ~ proportion + median_wdif),
                prop_adjusted_max = dplyr::case_when(proportion <0.1 ~ proportion - max_wdif,
                                                     proportion >=0.1 ~ proportion + max_wdif))
# if any estimate becomes negative, make it zero:
w_estimates_adjusted[w_estimates_adjusted$prop_adjusted_mean<0, "prop_adjusted_mean"] = 0
w_estimates_adjusted[w_estimates_adjusted$prop_adjusted_median<0, "prop_adjusted_median"] = 0
w_estimates_adjusted[w_estimates_adjusted$prop_adjusted_max<0, "prop_adjusted_max"] = 0

# plotting lines from higher to lower prop

ordered_samples = w_estimates_adjusted%>%
  dplyr::group_by(line) %>%
  dplyr::mutate(max_prop = max(proportion))

ordered_samples$line = factor(ordered_samples$line,
                              levels = rev(unique(ordered_samples[order(ordered_samples$max_prop),"line"])$line),
                              ordered = T)

length(unique(ordered_samples$line)) # 143

pdf("../../data/results/1.check_line_proportions/top_bottom_line_proportions_per_line.pdf",
    width = 13, height = 10 )
ordered_samples %>%
  ggplot( aes(x=line,y=proportion, color=top_or_bottom)) +
  geom_point(binaxis='y', stackdir='center', alpha = 0.6) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust =1)) +
  scale_color_manual(values=c("#1d3557","#EDA4BD","#9d0208")) +
  ylab("line proportions") +
  ggtitle("Proportions per well per line") +
  coord_flip()
dev.off()

# changes in line proportions
line_prop_changes = w_estimates_adjusted %>%
  dplyr::group_by(line,well) %>%
  dplyr::reframe(log_fraction_mean = log(prop_adjusted_mean[top_or_bottom == "bottom"] / prop_adjusted_mean[top_or_bottom == "top"]),
                 log_fraction_median = log(prop_adjusted_median[top_or_bottom == "bottom"] / prop_adjusted_median[top_or_bottom == "top"]),
                 log_fraction_max = log(prop_adjusted_max[top_or_bottom == "bottom"] / prop_adjusted_max[top_or_bottom == "top"]),
                 prop_unadjusted_max_value = max(proportion),
                 prop_unadjusted_min_value = min(proportion),
                 pool=pool,condition=condition,treatment=treatment,replicate = replicate) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()
pdf ("../../data/results/1.check_line_proportions/change_in_total_top_bottom_line_proportions_per_line_mean_error.pdf",
     width = 13, height = 9 )
line_prop_changes %>%
  
  ggplot( aes(x=line,y=log_fraction_mean,color = condition,
              size = prop_unadjusted_max_value)) +
  geom_point(binaxis='y', stackdir='center', alpha = 0.6) +
  theme_bw() +
  scale_color_brewer(palette  = "Set1") +
  theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust =1)) +
  ylab("log(bottom prop / top prop)") +
  ggtitle("log(bottom prop / top prop) [adjust: mean error]")

dev.off()

## add sex info
line_prop_changes = line_prop_changes %>%
  dplyr::left_join(.,line_info)


line_prop_changes %>%
  write.csv("../../data/results/1.check_line_proportions/line_prop_changes_per_well.csv",
            quote = FALSE,row.names = FALSE)
# C5a should show more migration to bottom well than media only.
# where is the migration info for all pools? How to incorporate this info to the model?

# now, averaging wells and then averaging pools, for those that are shared ####
###########

line_prop_changes %>%
  # first group by rep, and calculate rep mean and variance
  dplyr::group_by(line,condition,treatment,pool) %>%
  dplyr::filter(condition %in% c("media_only", "3nM_C5a")) %>%
  dplyr::filter(!log_fraction_mean %in% c(NA,NaN,Inf,-Inf)) %>%
  # dplyr::filter(prop_unadjusted_max_value > 0.01) %>%. # trying a less stringert filter
  dplyr::filter(prop_unadjusted_max_value > 0.005) %>% # see another filter by variance below
  # remember in log_fraction_mean, the "mean" refers to mean error used to adjust the proportions
  dplyr::summarise(mean_outcome_rep = mean(log_fraction_mean,na.rm=TRUE),
                   var_outcome_rep = var(log_fraction_mean,na.rm=TRUE),
                   mean_max_prop_rep = mean(prop_unadjusted_max_value,na.rm=TRUE),
                   var_max_prop_rep = var(prop_unadjusted_max_value,na.rm=TRUE)) %>%
  ungroup() %>%
  # remove lines with high variability between replicates
  dplyr::filter(var_outcome_rep < 1) %>% 
  # then, calculate mean and variance across pools 
  # variance is variance of the means
  dplyr::group_by(line,condition,treatment) %>% 
  # calculating mean of means, and var of means
  dplyr::mutate(mean_outcome_pool = mean(mean_outcome_rep,na.rm=TRUE),
                var_outcome_pool = var(mean_outcome_rep,na.rm=TRUE),
                mean_max_prop_pool = mean(mean_max_prop_rep,na.rm=TRUE),
                var_max_prop_pool = var(mean_max_prop_rep,na.rm=TRUE),
                n_pools = n()) %>%
  ungroup() %>%
  dplyr::left_join(.,line_info) %>%
  write.csv("../../data/results/1.check_line_proportions/line_prop_changes_averages_variances_1pct_filtered.csv",
            quote = FALSE,row.names = FALSE)


#
# checking relationship of variance from pQTL outcome (log of fraction) with min or max or average proportion in the fraction
# in the well replicates, and also in the pool replicates from shared lines
# to try to inform a reasonable threshold for prefiltering based on proportion of the line

# calculate variance between replicates
line_prop_changes_shared = line_prop_changes %>%
  dplyr::filter(condition %in% c("media_only", "3nM_C5a")) %>%
  dplyr::group_by(line,condition,treatment,pool) %>%
  dplyr::filter(!log_fraction_mean %in% c(NA,NaN,Inf,-Inf)) %>%
  # remember in log_fraction_mean, the "mean" refers to mean error used to adjust the proportions
  dplyr::summarise(mean_outcome = mean(log_fraction_mean, na.rm=TRUE),
                   var_outcome = var(log_fraction_mean, na.rm=TRUE),
                   mean_max_prop = mean(prop_unadjusted_max_value, na.rm=TRUE),
                   var_max_prop = var(prop_unadjusted_max_value, na.rm=TRUE),
                   CV_max_prop = sqrt(var(prop_unadjusted_max_value)) / mean(prop_unadjusted_max_value)) %>%
  
  ungroup()

p = line_prop_changes_shared %>%
  group_by(line) %>%
  ggplot( aes(x=log10(mean_max_prop),y=var_outcome,col=pool)) +
  geom_point(aes(size = CV_max_prop), alpha = 0.6) +
  theme_bw() +
  geom_vline(xintercept = -2, col="grey", linetype="dashed") +
  ylab("Variance of fractions") +
  xlab("log10(Mean of maximum proportion) -> larger prop. to right") +
  ggtitle("Mean(max prop) vs var(fractions) ", subtitle = "Across replicates.") +
  facet_wrap(vars(condition,treatment), ncol = 3)

pdf("../../data/results/1.check_line_proportions/mean_vs_var_all_lines_by_well_reps.pdf",
    width = 13, height = 9)
p
dev.off()
# checking problematic reps
line_prop_changes %>% dplyr::filter(line == "ouvb_2" & condition == "3nM_C5a" & treatment == "LPS" & pool == "pool8")
line_prop_changes %>% dplyr::filter(line == "aion_2" & condition == "3nM_C5a" & treatment == "IFN" & pool == "pool3")
# checking mean_max_prop supports removing lines under 0.5% (0.005) proportion
# what's with the missing values?
line_prop_changes %>% dplyr::filter(line == "aehn_22" & condition == "3nM_C5a" & treatment == "LPS" & pool == "pool6")
line_prop_changes %>% dplyr::filter(line == "aipt_33" & condition == "3nM_C5a" & treatment == "IFN" & pool == "pool6")
line_prop_changes %>% dplyr::filter(line == "miaj_4" & condition == "3nM_C5a" & treatment == "IFN" & pool == "pool8")
# these are tiny, disappearing fractions (hence NA in one rep). 
# Remove values for which there is only one well replicate because of vanishingly small proportions?

# now checking for shared lines, variation across pools
line_prop_changes_shared  = line_prop_changes %>%
  dplyr::group_by(line,condition,treatment) %>%
  dplyr::filter(!log_fraction_mean %in% c(NA,NaN,Inf,-Inf)) %>%
  # Filter for lines that are present in more than one pool
  dplyr::filter(n_distinct(pool) > 1) %>%
  dplyr::ungroup()

# first group by rep, and calculate rep mean and variance
line_prop_changes_shared =  line_prop_changes_shared %>%
  dplyr::group_by(line,condition,treatment,pool) %>%
  dplyr::filter(!log_fraction_mean %in% c(NA,NaN,Inf,-Inf)) %>%
  # remember in log_fraction_mean, the "mean" refers to mean error used to adjust the proportions
  dplyr::summarise(mean_outcome_rep = mean(log_fraction_mean, na.rm=TRUE),
                   var_outcome_rep = var(log_fraction_mean, na.rm=TRUE),
                   mean_max_prop_rep = mean(prop_unadjusted_max_value, na.rm=TRUE),
                   var_max_prop_rep = var(prop_unadjusted_max_value, na.rm=TRUE)) %>%
  ungroup() %>%
  
  # then, calculate mean and variance across pools 
  # variance is variance of the means
  dplyr::group_by(line,condition,treatment) %>% 
  # calculating mean of means, and var of means
  dplyr::summarise(mean_outcome_pool = mean(mean_outcome_rep, na.rm=TRUE),
                   var_outcome_pool = var(mean_outcome_rep, na.rm=TRUE),
                   mean_max_prop_pool = mean(mean_max_prop_rep, na.rm=TRUE),
                   var_max_prop_pool = var(mean_max_prop_rep, na.rm=TRUE),
                   n_pools = n()) %>%
  ungroup() 

p = line_prop_changes_shared %>%
  group_by(line) %>%
  ggplot( aes(x=log10(mean_max_prop_pool),y=var_outcome_pool,col=line)) +
  geom_point(aes(size=n_pools), alpha = 0.6) +
  theme_bw() +
  geom_vline(xintercept = -2, col="grey", linetype="dashed") +
  ylab("Variance of fractions") +
  xlab("log10(Mean of maximum proportion) -> larger prop. to right") +
  ggtitle("Variance of log(bottom prop/top prop) in shared lines vs maximum proportion", subtitle = "mean(mean[max prop] across rep, across pools) vs var(mean[fractions] per rep, across pools)") +
  facet_wrap(vars(condition,treatment), ncol = 3)

pdf("../../data/results/1.check_line_proportions/mean_vs_var_shared_lines_across_pools.pdf",
    width = 13, height = 9)
p
dev.off()

# check variation between pools
line_prop_changes_shared = line_prop_changes %>%
  dplyr::group_by(line,condition,treatment) %>%
  dplyr::filter(!log_fraction_mean %in% c(NA,NaN,Inf,-Inf)) %>%
  # Filter for lines that are present in more than one pool
  dplyr::filter(n_distinct(pool) > 1) %>%
  dplyr::ungroup() %>%
  # calculate mean and variance for pools and wells, per line, condition and treatment
  dplyr::group_by(line,condition,treatment,pool) %>% 
  # remember in log_fraction_mean, the "mean" refers to mean error used to adjust the proportions
  dplyr::summarise(mean_outcome = mean(log_fraction_mean, na.rm=TRUE),
                   var_outcome = var(log_fraction_mean, na.rm=TRUE),
                   mean_max_prop = mean(prop_unadjusted_max_value, na.rm=TRUE),
                   var_max_prop = var(prop_unadjusted_max_value, na.rm=TRUE),
                   CV_max_prop = sqrt(var(prop_unadjusted_max_value)) / mean(prop_unadjusted_max_value))

p = line_prop_changes_shared %>%
  group_by(line,pool) %>%
  dplyr::filter(!mean_max_prop %in% c(0,NA,Inf,-Inf)) %>%
  # remove 0 to avoid inf after log
  ggplot( aes(x=log10(mean_max_prop),y=var_outcome,col=pool)) +
  geom_point(aes(size = CV_max_prop), alpha = 0.6) +
  theme_bw() +
  ylab("Variance of fractions") +
  xlab("log10(Mean of maximum proportion) -> larger prop. to right") +
  ggtitle("Mean(max prop) vs var(fractions) ", subtitle = "for all replicates, per line and pool, in shared lines") +
  facet_wrap(vars(condition,treatment), ncol = 3)
pdf("../../data/results/1.check_line_proportions/mean_vs_var_all_wells_per_pool.pdf",
    width = 13, height = 9)
p
dev.off()

# check variation within pools, between wells
# fix the renaming of wells to replicates, it's not working now

# line_prop_changes_shared = line_prop_changes %>%
#   dplyr::group_by(line,condition,treatment) %>%
#   dplyr::filter(!log_fraction_mean %in% c(NA,Inf,-Inf)) %>%
#   # Filter for lines that are present in more than one pool
#   dplyr::filter(n_distinct(pool) > 1) %>%
#   dplyr::ungroup() %>%
#   # calculate mean and variance for pools and wells, per line, condition and treatment
#   dplyr::group_by(line,pool,treatment, condition, well) %>% 
#   dplyr::mutate(replicate_column = paste0("rep_", row_number())) %>%
#   
# 
# p = line_prop_changes_shared %>%
#   group_by(line) %>%
#   ggplot( aes(x=mean_outcome,y=var_outcome,col=line)) +
#   geom_point(aes(size = CV_max_prop), alpha = 0.6) +
#   theme_bw() +
#   ylab("Variance of fractions for pools and wells") +
#   xlab("Mean of fractions for pools and wells") +
#   ggtitle("Mean-var relationship", subtitle = "for all pools and wells, per line, in shared lines") +
#   facet_wrap(vars(condition,treatment), ncol = 3)
# p

############ from here on it's legacy code ##########
#####################################################
# 
# # Compare with errors from simulations
# # Assuming coverage of 2 (min mean depth for these samples = 2.5, shouldn't matter much in terms of max errors)
# w_estimate_sim_C2 =  read.table(
#   paste0(
#     "../../../simulation_WGS_exome_data/data/pool_19d_another/WGS/w/w_estimate_line_proportions_fixed_100_trials_WGS_C2.txt.gz"))
# 
# colnames(w_estimate_sim_C2) = c("line","w_est","trial","coverage")
# w_estimate_sim_C2$coverage = gsub("C","",w_estimate_sim_C2$coverage)
# w_real = read.table(
#   paste0(
#     "../../../simulation_WGS_exome_data/data/pool_19d_another/w/line_proportions_fixed_10_sampling_reps_100_trials.txt.gz"), header = T)
# # retain only the first sampling replicate (only one taken to simulate WGS)
# w_real = w_real[w_real$sampling_rep==1,1:ncol(w_real)-1]
# w_real = w_real %>% tidyr::pivot_longer(!trial,values_to = "w_real",names_to = "line")
# w_errors = inner_join(w_estimate_sim_C2, w_real, by = c("trial", "line"))
# w_errors$w_diff = w_errors$w_real - w_errors$w_est
# 
# w_errors = w_errors %>%
#   group_by(w_real) %>%
#   summarise(max_error = max(abs(w_diff)))
# 
# # assign error by closest w_real value from proportion
# ordered_samples$sim_error = NA
# ordered_samples$sim_error = as.double(ordered_samples$sim_error)
# for(n in 1:nrow(ordered_samples)){
#   prop = ordered_samples[n,"proportion"]
#   ordered_samples[n,"sim_error"]= w_errors[which.min(abs(w_errors$w_real-prop$proportion)),"max_error"]$max_error
# }
# 
# samples$sim_error = NA
# samples$sim_error = as.double(samples$sim_error)
# for(n in 1:nrow(samples)){
#   prop = samples[n,"proportion"]
#   samples[n,"sim_error"]= w_errors[which.min(abs(w_errors$w_real-prop$proportion)),"max_error"]$max_error
# }
# 
# 
# pdf ("../../data/results/1.check_line_proportions/change_in_total_top_bottom_line_proportions_per_line_sim_error.pdf",
#     width = 13, height = 10 )
# ordered_samples %>%
#   ggplot( aes(x=cells_format_assay,y=proportion, color=Top_or_bottom)) +
#   geom_point(aes(size = DNA_normalised_to_total), binaxis='y', stackdir='center',position = position_dodge(), alpha = 0.6) +
#   geom_errorbar(aes(ymin = proportion - sim_error, ymax = proportion + sim_error),
#                 position = position_dodge(width = 0.5)) +
#   theme_bw() +
#   theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust =1)) +
#   scale_color_manual(values=c("#1d3557","#EDA4BD","#9d0208")) +
#   ylab("line proportions") +
#   ggtitle("Proportions per well per line") +
#   coord_flip() +
#   facet_wrap(~line, ncol = 5, scales = "free_x" )
# dev.off()
# 
# pdf ("../../data/results/1.check_line_proportions/change_in_total_top_bottom_line_proportions_per_well_sim_error.pdf",
#     width = 13, height = 9 )
# samples %>%
#   ggplot( aes(x=line,y=proportion, color=Top_or_bottom)) +
#   geom_point(aes(size = DNA_normalised_to_total), binaxis='y', stackdir='center',position = position_dodge(), alpha = 0.6) +
#   geom_errorbar(aes(ymin = proportion - sim_error, ymax = proportion + sim_error),
#                 position = position_dodge(width = 0.5)) +
#   theme_bw() +
#   theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust =1)) +
#   scale_color_manual(values=c("#1d3557","#EDA4BD","#9d0208")) +
#   ylab("line proportions") +
#   ggtitle("Proportions per line per well") +
#   # coord_flip() +
#   facet_wrap(~cells_format_assay, ncol = 3, scales = "free_y" )
# 
# dev.off()
# 
# # check changes in Total proportions between IFN and no IFN
# pdf ("../../data/results/1.check_line_proportions/change_in_total_line_proportions_per_line.pdf",
#     width = 11, height = 10 )
# ordered_samples %>%
#   filter(Top_or_bottom == "Total") %>%
#   ggplot( aes(x=Cells,y=proportion, color=Assay_condition, shape = Format)) +
#   geom_point(binaxis='y', stackdir='center', alpha = 0.6, size = 5) +
#   theme_bw() +
#   theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust =1)) +
#   scale_color_manual(values=c("#ffbe0b","#003049","#8338ec","#ff006e")) +
#   ylab("line proportions") + ggtitle("Change in total line proportions per line") +
#   facet_wrap(~line, ncol = 4, scales = "free_y" )
# dev.off()
# 
# 
# 
# pdf ("../../data/results/1.check_line_proportions/change_in_total_line_proportions_per_line_sim_errors.pdf",
#     width = 11, height = 10 )
# ordered_samples %>%
#   filter(Top_or_bottom == "Total") %>%
#   ggplot( aes(x=Cells,y=proportion, color=Assay_condition, shape = Format)) +
#   geom_pointrange(aes(ymin = proportion-sim_error, ymax = proportion+sim_error),
#                   position = position_dodge(width=0.5),
#                   alpha = 0.6, size = 0.5) +
#   theme_bw() +
#   theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust =1)) +
#   scale_color_manual(values=c("#ffbe0b","#003049","#8338ec","#ff006e")) +
#   ylab("line proportions") + ggtitle("Change in total proportions per line") +
#   facet_wrap(~line, ncol = 4, scales = "free_y" )
# dev.off()
# 
# 
# # More variability in the 6 well? CV
# samples %>%
#   filter(Top_or_bottom == "Total") %>%
#   group_by(Cells, Format, line) %>%
#   summarise(sd_IFN_format = sd(proportion), mean_IFN_format = mean(proportion), n = n()) %>%
#   mutate(CV_IFN_format = sd_IFN_format / mean_IFN_format) %>%
#   ggplot( aes(x=line,y=CV_IFN_format, color=Cells, shape = Format)) +
#   geom_point(binaxis='y', stackdir='center', alpha = 0.6, size = 5) +
#   theme_bw() +
#   theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust =1)) +
#   # scale_color_manual(values=c("#1d3557","#EDA4BD","#9d0208")) +
#   ylab("Coefficient of variation")
# 
# # check coefficients of variation divided by Format and/or IFN status
# 
# samples %>%
#   filter(Top_or_bottom == "Total") %>%
#   group_by(Cells, Format, line) %>%
#   summarise(sd_IFN_format = sd(proportion), mean_IFN_format = mean(proportion), n = n()) %>%
#   mutate(CV_IFN_format = sd_IFN_format / mean_IFN_format) %>%
#   ggplot( aes(x=Format,y=CV_IFN_format, color=Cells)) +
#   geom_jitter() +
#   theme_bw() +
#   theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust =1)) +
#   # scale_color_manual(values=c("#1d3557","#EDA4BD","#9d0208")) +
#   ylab("Coefficient of variation")
# 
# # zoom in on lines that are over 1% mean proportion
# 
# samples %>%
#   filter(Top_or_bottom == "Total") %>%
#   group_by(Cells, Format, line) %>%
#   summarise(sd_IFN_format = sd(proportion), mean_IFN_format = mean(proportion), n = n()) %>%
#   mutate(CV_IFN_format = sd_IFN_format / mean_IFN_format) %>%
#   filter(mean_IFN_format > 0.01) %>%
#   ggplot( aes(x=Format,y=CV_IFN_format, color=Cells)) +
#   geom_jitter() +
#   theme_bw() +
#   theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust =1)) +
#   # scale_color_manual(values=c("#1d3557","#EDA4BD","#9d0208")) +
#   ylab("Coefficient of variation")
# # There seems to be more variability in the variation at mean proportion per line over 1%
# 
# # Measure mean and 95% CI per per "cells" (IFN/not) per line
# pdf ("../../data/results/1.check_line_proportions/mean_prop_per_line_per_IFN_per_format.pdf",
#     width = 7, height = 4 )
# samples %>%
#   filter(Top_or_bottom == "Total") %>%
#   group_by(Cells, Format, line) %>%
#   summarise(sd_IFN_format =sd(proportion) ,mean_IFN_format = mean(proportion), n = n()) %>%
#   mutate(CI95_upper = mean_IFN_format + 1.96*sd_IFN_format,
#          CI95_lower = ifelse(mean_IFN_format - 1.96*sd_IFN_format<0,0,
#                              mean_IFN_format - 1.96*sd_IFN_format) ) %>%
#   ggplot( aes(x=line,y=mean_IFN_format, color=Cells, shape=Format)) +
#   geom_point(position = position_dodge(0.4)) +
#   geom_errorbar(aes(ymin=CI95_lower, ymax=CI95_upper), width=.1, position=position_dodge(0.4)) +
#   scale_color_manual(values = c("#ff6b35","#004e89"))+
#   theme_bw() +
#   theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust =1)) +
#   # scale_color_manual(values=c("#1d3557","#EDA4BD","#9d0208")) +
#   ylab("Mean proportion per line (with 95% CI)")
# dev.off()
# 
# pdf ("../../data/results/1.check_line_proportions/mean_prop_per_line_per_IFN_per_format_se.pdf",
#     width = 7, height = 4 )
# samples %>%
#   filter(Top_or_bottom == "Total") %>%
#   group_by(Cells, Format, line) %>%
#   summarise(se_IFN_format = std.error(proportion),mean_IFN_format = mean(proportion), n = n()) %>%
#   mutate(se_upper = mean_IFN_format + se_IFN_format,
#          se_lower = ifelse(mean_IFN_format - se_IFN_format<0,0,
#                            mean_IFN_format - se_IFN_format) ) %>%
#   ggplot( aes(x=line,y=mean_IFN_format, color=Cells, shape=Format)) +
#   geom_point(position = position_dodge(0.4)) +
#   geom_errorbar(aes(ymin=se_lower, ymax=se_upper), width=.1, position=position_dodge(0.4)) +
#   scale_color_manual(values = c("#ff6b35","#004e89"))+
#   theme_bw() +
#   theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust =1)) +
#   # scale_color_manual(values=c("#1d3557","#EDA4BD","#9d0208")) +
#   ylab("Mean proportion per line (+/- SEM)")
# dev.off()
# # Differences in means seem to be driven by IFN status, not so much the 6/12 well format
# 
# 
# # compare proportions to previous differentiations
# chrmX_diff = read.table("../../../OTAR2065_NovaSeq_run42058_chromiumX_loading_test/data/w/chrmX_WGS_w_estimate.txt")
# colnames(chrmX_diff) = c("line", "proportion")
# chrmX_diff$experiment = "chrmX_diff"
# reichGM_media_comp_diff = read.table("../../../OTAR2065_NovaSeq_run41362_WGS_media_and_macPrecursors_october2021/data/w/3D_MICR2c3_160721_w_estimate.txt")
# colnames(reichGM_media_comp_diff) = c("line", "proportion")
# reichGM_media_comp_diff$experiment = "reichGM_media_comp_diff"
# merged_WGS = rbind(chrmX_diff,reichGM_media_comp_diff)
# merged_WGS[merged_WGS$proportion<0,"proportion"] = 0
# 
# new_WGS = samples %>%
#   filter(Top_or_bottom == "Total") %>%
#   select(line, proportion, cells_format_assay)
# colnames(new_WGS) = c("line", "proportion", "experiment")
# 
# 
# # 89 and 101 are exactly the same?
# samples %>% filter(Well == "89" | Well == "101")
# # not exactly
# 
# 
# # IFN stimulates macrophage migration
# # Test DNA in bottom well for IFN vs no IFN, and within same IFN status the effect of C5a + inhibitors
# pdf ("../../data/results/1.check_line_proportions/DNA_in_bottom_well_over_total_boxplot_conditions.pdf",
#     width = 7, height = 4 )
# samples %>%
#   filter(Top_or_bottom == "Bottom", line == "letw_1" ) %>%
#   ggplot( aes(x=Cells,y=DNA_normalised_to_total)) +
#   geom_boxplot(color="grey", alpha=0.2) +
#   geom_jitter(aes(color = Assay_condition),shape=16, position=position_jitter(0.2),size = 3) +
#   theme_bw() +
#   theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust =1)) +
#   scale_color_manual(values=c("#ffbe0b","#003049","#8338ec","#ff006e")) +
#   ylab("DNA in bottom well over total")
# 
# dev.off()
# # The conditions behave as expected in the IFN case (+ inhibitors or without C5a have lower bottom DNA prop)
# # But not in the no_IFN case. Why?
