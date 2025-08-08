# checking line proportions of phagocytosis assay
# Pools of 19-24 lines + sh5yhy, phagocytosis assay with different replicates (untreated, IFN, LPS)
# pre-phagocytosis assay (no cargo), and low (C5) and high (C3) mCherry (phagocytosis) bins 
# main objective of this assay:
# how do the line proportions change using WGS for all the samples?

library(patchwork)
library(tidyverse)
#library(plotrix)
source("../functions.R")

outdir="../../../data/results/migration/1.check_line_proportions"
dir.create(file.path(outdir),recursive = TRUE)
sample_info=read.csv("../../../data/all_pools_migration_sample_info.csv")
donor_info=read_csv("../../../../OTAR2065_differentiation_efficiency/data/donor_metadata_complete_with_imputed_sex.csv") %>%
  dplyr::rename(sex=Sex) %>%
  dplyr::select(donor,sex)
w_estimates = read_w_estimates(sample_info,
                               "../../../data/w/",
                               assay = "migration",
                               pools = paste0("pool", c(2:11,13:17))) # we omit pool 12 because of one donor making the majority of the pool
# adding donor
w_estimates = w_estimates %>%
  dplyr::rename(line = donor) %>%
  dplyr::mutate(donor = if_else(str_detect(line,"-"),
                                true=str_split_i(line,pattern="-",i=1), 
                                false=str_split_i(line,pattern="_",i=1))) %>%
  dplyr::relocate(line,donor)

# convert slightly negative values to 0
w_estimates[w_estimates$proportion<0, "proportion"] = 0
write.table(w_estimates,paste0(outdir,"/w_estimates.txt"),
            row.names = FALSE,quote = FALSE)
# read in proportion error estimates
error_estimates = read.table("../../../data/w/error_approximations_generic_pool.txt", header = TRUE) %>%
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
  dplyr::mutate(max_prop = max(proportion)) %>%
  dplyr::ungroup()

ordered_samples$line = factor(ordered_samples$line,
                              levels = rev(unique(ordered_samples[order(ordered_samples$max_prop),"line"])$line),
                              ordered = T)

length(unique(ordered_samples$line)) # 261


pdf("../../../data/results/migration/1.check_line_proportions/raw_top_bottom_line_proportions_per_line_incl_all_condition.pdf",
    width = 7, height = 17 )
ordered_samples %>%
  ggplot( aes(x=line,y=proportion, color=top_or_bottom)) +
  geom_point(binaxis='y', stackdir='center', alpha = 0.6) +
  theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust =1)) +
  scale_color_manual(values=c("#1d3557","#EDA4BD","#9d0208")) +
  ylab("line proportions") +
  ggtitle("Proportions per well per line") +
  coord_flip() + 
  theme_bw()+
  facet_wrap(vars(condition)) 
dev.off()

# recode media_only and samples with PMX53 as "-C5a"
# and 3nM_C5a to "+C5a" (shortening names)
w_estimates_adjusted = w_estimates_adjusted %>%
  dplyr::mutate(condition = case_when(condition == "3nM_C5a" ~ "+C5a",
                                      .default = "-C5a"))



##### The changes in proportions should all go in same direction irrespective of the pool ###
### examining pool effects in directionality and variance ####
# see 1.3.more_checks_on_batches.R

##### examining the fractions ######
line_prop_changes = w_estimates_adjusted %>%
  dplyr::group_by(line,well) %>%
  dplyr::reframe(log_fraction_mean = log(prop_adjusted_mean[top_or_bottom == "bottom"] / prop_adjusted_mean[top_or_bottom == "top"]),
                 prop_unadjusted_max_value = max(proportion),
                 prop_unadjusted_min_value = min(proportion),
                 pool=pool,condition=condition,treatment=treatment,replicate = replicate,donor=donor) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()

## normalise (center) against -C5a average value per line, pool and treatment

centering_values = line_prop_changes %>%
  dplyr::filter(condition == "-C5a") %>%
  dplyr::group_by(line,pool,treatment) %>%
  dplyr::summarise(
    centering_value = mean(log_fraction_mean)
  ) %>%
  ungroup()

line_prop_changes_scaled_well = line_prop_changes %>%
  dplyr::left_join(centering_values) %>%
  dplyr::filter(condition == "+C5a") %>%
  dplyr::mutate(log_fraction_centered_media_only = NA) %>%
  dplyr::group_by(well) %>%
  dplyr::mutate(log_fraction_centered_media_only = log_fraction_mean - centering_value) %>%
# centering to the average of the well and scaling (as in phagocytosis)
  dplyr::filter(!log_fraction_centered_media_only %in% c("-Inf","Inf","NA","NaN")) %>%
  dplyr::mutate(scaled_log_fraction = scale_this(log_fraction_centered_media_only)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!scaled_log_fraction %in% c("-Inf","Inf","NA","NaN")) %>%
  dplyr::filter(!is.na(scaled_log_fraction))

## add sex info
line_prop_changes_scaled_well = line_prop_changes_scaled_well %>%
  dplyr::left_join(.,donor_info)

## histograms

p1 = ggplot(line_prop_changes_scaled_well,aes(x=scaled_log_fraction)) + 
  geom_histogram(bins = 50) +
  theme_minimal() + 
  ggtitle("s&c (well mean) log fraction:  phagocytosis",
          subtitle = "scaled first to -C5a and then per well")

pdf (paste0(outdir,"/histogram_line_prop_changes_scaled_well.pdf"),
     width = 6, height = 6 )
p1
dev.off()


### what proportion of donors are 2sd away from the mean ###
summary_stats = line_prop_changes_scaled_well %>%
  summarize(
    mean_log = mean(scaled_log_fraction),
    sd_log = sd(scaled_log_fraction),
    mean_scaled_well = mean(scaled_log_fraction),
    sd_scaled_well = sd(scaled_log_fraction),
    
  ) %>% 
  dplyr::mutate(
    lower_threshold_log = mean_log - 3 * sd_log,
    upper_threshold_log = mean_log + 3 * sd_log,
    lower_threshold_scaled_well = mean_scaled_well - 3 * sd_scaled_well,
    upper_threshold_scaled_well = mean_scaled_well + 3 * sd_scaled_well,
  )

# Count outliers
outliers_count_log = line_prop_changes_scaled_well %>%
  mutate(
    is_outlier = log_fraction_mean < summary_stats$lower_threshold_log | log_fraction_mean > summary_stats$upper_threshold_log
  ) %>%
  filter(is_outlier == TRUE) %>%
  nrow()

outliers_count_scaled_well = line_prop_changes_scaled_well %>%
  mutate(
    is_outlier = scaled_log_fraction < summary_stats$lower_threshold_scaled_well | scaled_log_fraction > summary_stats$upper_threshold_scaled_well
  ) %>%
  filter(is_outlier == TRUE) %>%
  nrow()

# Calculate the percentage of outliers
(outliers_count_log / nrow(line_prop_changes_scaled_well)) * 100 # 0.24%
(outliers_count_scaled_well / nrow(line_prop_changes_scaled_well)) * 100. # 1.27%
# few outliers in both

### is there any particular pool that skews the values?
  
  pdf (paste0(outdir,"/scaled_centered_media_only_well_average_log_fractions_per_pool_treatment_line.pdf"),
       width = 7, height = 4 )
  ordered = line_prop_changes_scaled_well %>%
    dplyr::arrange(desc(prop_unadjusted_max_value)) 
  ordered$line = factor(ordered$line,levels = unique(ordered$line), ordered = TRUE) 
  for(lin in unique(ordered$line)){
    subset = ordered %>%
      dplyr::filter(line %in% lin)
    if(sum(is.na(subset$log_fraction_mean)) == nrow(subset)){
      next()
    }else{
      subset$pool = factor(subset$pool,levels = sort(unique(subset$pool)), ordered = TRUE)
      subset$treatment = factor(subset$treatment,levels = c("untreated","IFN","LPS"), ordered = TRUE)
      p =  subset %>%
        ggplot(aes(x = pool,y = scaled_log_fraction, col = treatment)) +
        geom_point() +
        geom_boxplot()+
        geom_hline(yintercept = 0,linetype = "dotted", col="grey") + 
        theme_minimal() + 
        ggtitle(lin)
      plot(p)
    }
  }
  dev.off()
  
  # are fractions more extreme and more variable the smaller the donor proportion?
  
  pdf (paste0(outdir,"/mean_stdev_log_fractions_per_treatment_line.pdf"),
       width = 4, height = 12 )
  ordered %>%
    dplyr::mutate(treatment = factor(treatment,levels = c("untreated","IFN","LPS"), ordered = TRUE)) %>%
    dplyr::group_by(line,treatment) %>%
    dplyr::summarise(mean_fraction = mean(log_fraction_mean),
                     st_dev_mean_fraction = sd(log_fraction_mean)) %>%
    dplyr::filter(!is.na(mean_fraction)) %>%
    dplyr::ungroup() %>%
    ggplot(aes(x=line,y=mean_fraction,size=st_dev_mean_fraction,col=treatment))+
    geom_point() +
    theme_minimal() +
    geom_hline(yintercept = 0,linetype = "dotted", col="grey") + 
    xlab("Line (--> smaller proportions)") +
    ylab("Mean fraction across pools and replicates") +
    coord_flip() 
  # lines tend to have more extreme fraction values, and these are more variable, the smaller the line proportion
  
  dev.off()
  
  p1 = ordered %>%
    dplyr::mutate(treatment = factor(treatment,levels = c("untreated","IFN","LPS"), ordered = TRUE)) %>%
                 
    dplyr::group_by(line,treatment) %>%
    dplyr::summarise(mean_fraction = mean(log_fraction_mean),
                     st_dev_mean_fraction = sd(log_fraction_mean),
                     prop_unadjusted_mean = mean(prop_unadjusted_max_value)) %>%
    dplyr::mutate( prop_bin = case_when(prop_unadjusted_mean < 0.01 ~ "under_1%",
                                        0.01 <= prop_unadjusted_mean & prop_unadjusted_mean < 0.1  ~ "1%-10%",
                                        prop_unadjusted_mean>=0.1 ~ "10%")) %>%
    dplyr::mutate(prop_bin = factor(prop_bin,levels = c("under_1%","1%-10%", "10%"))) %>%
    dplyr::filter(!is.na(mean_fraction)) %>%
    dplyr::ungroup() %>%
    ggplot(aes(x=prop_bin,y=mean_fraction,size=st_dev_mean_fraction,col=treatment))+
    geom_boxplot() +
    # geom_point(alpha = 0.1) +
    theme_minimal() +
    geom_hline(yintercept = 0,linetype = "dotted", col="grey") + 
    xlab("Mean maximum proportion in fraction") +
    ylab("Mean fraction across pools and replicates") 
  
 p2 =  ordered %>%
    dplyr::mutate(treatment = factor(treatment,levels = c("untreated","IFN","LPS"), ordered = TRUE)) %>%
    
    dplyr::group_by(line,treatment) %>%
    dplyr::summarise(mean_fraction = mean(log_fraction_mean),
                     st_dev_mean_fraction = sd(log_fraction_mean),
                     prop_unadjusted_mean = mean(prop_unadjusted_max_value)) %>%
    dplyr::mutate( prop_bin = case_when(prop_unadjusted_mean < 0.01 ~ "under_1%",
                                        0.01 <= prop_unadjusted_mean & prop_unadjusted_mean < 0.1  ~ "1%-10%",
                                        prop_unadjusted_mean>=0.1 ~ "10%")) %>%
    dplyr::mutate(prop_bin = factor(prop_bin,levels = c("under_1%","1%-10%", "10%"))) %>%
    dplyr::filter(!is.na(mean_fraction)) %>%
    dplyr::ungroup() %>%
    ggplot(aes(x=mean_fraction,y=st_dev_mean_fraction,col=prop_bin))+
     geom_point(alpha = 0.4) +
    theme_minimal() +
    geom_vline(xintercept = 0,linetype = "dotted", col="grey") + 
    xlab("Mean fraction across pools and replicates") +
    ylab("Stdev fraction across pools and replicates") +
    facet_wrap(vars(treatment))
    
 pdf (paste0(outdir,"/mean_stdev_log_fractions_per_treatment_max_prop.pdf"),
      width = 7, height = 7 )
 p1 / p2
 
 dev.off()
 
 # same for scaled
 ordered = line_prop_changes_scaled_well %>%
   dplyr::arrange(desc(prop_unadjusted_max_value)) 
 ordered$line = factor(ordered$line,levels = unique(ordered$line), ordered = TRUE) 
 
 p1 = ordered %>%
   dplyr::mutate(treatment = factor(treatment,levels = c("untreated","IFN","LPS"), ordered = TRUE)) %>%
   
   dplyr::group_by(line,treatment) %>%
   dplyr::summarise(mean_fraction = mean(scaled_log_fraction),
                    st_dev_mean_fraction = sd(scaled_log_fraction),
                    prop_unadjusted_mean = mean(prop_unadjusted_max_value)) %>%
   dplyr::mutate( prop_bin = case_when(prop_unadjusted_mean < 0.01 ~ "under_1%",
                                       0.01 <= prop_unadjusted_mean & prop_unadjusted_mean < 0.1  ~ "1%-10%",
                                       prop_unadjusted_mean>=0.1 ~ "10%")) %>%
   dplyr::mutate(prop_bin = factor(prop_bin,levels = c("under_1%","1%-10%", "10%"))) %>%
   dplyr::filter(!is.na(mean_fraction)) %>%
   dplyr::ungroup() %>%
   ggplot(aes(x=prop_bin,y=mean_fraction,size=st_dev_mean_fraction,col=treatment))+
   geom_boxplot() +
   # geom_point(alpha = 0.1) +
   theme_minimal() +
   geom_hline(yintercept = 0,linetype = "dotted", col="grey") + 
   xlab("Mean maximum proportion in fraction") +
   ylab("Mean fraction across pools and replicates") 
 
 p2 =  ordered %>%
   dplyr::mutate(treatment = factor(treatment,levels = c("untreated","IFN","LPS"), ordered = TRUE)) %>%
   
   dplyr::group_by(line,treatment) %>%
   dplyr::summarise(mean_fraction = mean(scaled_log_fraction),
                    st_dev_mean_fraction = sd(scaled_log_fraction),
                    prop_unadjusted_mean = mean(prop_unadjusted_max_value)) %>%
   dplyr::mutate( prop_bin = case_when(prop_unadjusted_mean < 0.01 ~ "under_1%",
                                       0.01 <= prop_unadjusted_mean & prop_unadjusted_mean < 0.1  ~ "1%-10%",
                                       prop_unadjusted_mean>=0.1 ~ "10%")) %>%
   dplyr::mutate(prop_bin = factor(prop_bin,levels = c("under_1%","1%-10%", "10%"))) %>%
   dplyr::filter(!is.na(mean_fraction)) %>%
   dplyr::ungroup() %>%
   ggplot(aes(x=mean_fraction,y=st_dev_mean_fraction,col=prop_bin))+
   geom_point(alpha = 0.4) +
   theme_minimal() +
   geom_vline(xintercept = 0,linetype = "dotted", col="grey") + 
   xlab("Mean fraction across pools and replicates") +
   ylab("Stdev fraction across pools and replicates") +
   facet_wrap(vars(treatment))
 
 pdf (paste0(outdir,"/mean_stdev_centered_mean_only_scaled_fractions_per_treatment_max_prop.pdf"),
      width = 7, height = 7 )
 p1 / p2
 
 dev.off()
line_prop_changes_scaled_well %>%
  write.csv(paste0(outdir,"/line_prop_changes_per_well.csv"),
            quote = FALSE,row.names = FALSE)

line_prop_changes_scaled_well %>%
  # first group by rep, and calculate rep mean and variance
  dplyr::group_by(line,treatment,pool) %>%
  dplyr::filter(!scaled_log_fraction %in% c(NA,NaN,Inf,-Inf)) %>%
  dplyr::filter(prop_unadjusted_max_value >= 0.01) %>%
  # remember in log_fraction_mean, the "mean" refers to mean error used to adjust the proportions
  dplyr::reframe(mean_outcome_rep = mean(scaled_log_fraction,na.rm = TRUE),
                   var_outcome_rep = var(scaled_log_fraction,na.rm = TRUE),
                   mean_max_prop_rep = mean(prop_unadjusted_max_value,na.rm = TRUE),
                   var_max_prop_rep = var(prop_unadjusted_max_value,na.rm = TRUE),
                 sex=sex,donor=donor) %>%
  ungroup() %>%
  
  # then, calculate mean and variance across pools 
  # variance is variance of the means
  dplyr::group_by(line,treatment) %>% 
  # calculating mean of means, and var of means
  dplyr::mutate(mean_outcome_pool = mean(mean_outcome_rep,na.rm = TRUE),
                var_outcome_pool = var(mean_outcome_rep,na.rm = TRUE),
                mean_max_prop_pool = mean(mean_max_prop_rep,na.rm = TRUE),
                var_max_prop_pool = var(mean_max_prop_rep,na.rm = TRUE),
                n_pools = n()) %>%
  ungroup() %>%
  write.csv(paste0(outdir,"/line_prop_changes_averages_variances_1pct_filtered.csv"),
            quote = FALSE,row.names = FALSE)
###########

# number of lines before filters:
length(unique(line_prop_changes_scaled_well$line)) # 216 pools 2-17
# with NA, INF filters

# with 1% filters
line_prop_changes_scaled_well %>%
  # first group by rep, and calculate rep mean and variance
  dplyr::group_by(line,replicate,treatment,pool) %>%
  dplyr::filter(prop_unadjusted_max_value >= 0.01) %>% 
  dplyr::ungroup() %>%
  summarize(unique_elements = n_distinct(line)) # 158 pools 2-17


# with 0.5% filters
line_prop_changes_scaled_well %>%
  # first group by rep, and calculate rep mean and variance
  dplyr::group_by(line,replicate,treatment,pool) %>%
  dplyr::filter(prop_unadjusted_max_value >= 0.005) %>%
  dplyr::ungroup() %>%
  summarize(unique_elements = n_distinct(line)) # 188 pools 2-17
