# checking line proportions of phagocytosis assay
# Pools of 19-24 lines + sh5yhy, phagocytosis assay with different replicates (untreated, IFN, LPS)
# pre-phagocytosis assay (no cargo), and low (C5) and high (C3) mCherry (phagocytosis) bins 
# main objective of this assay:
# how do the line proportions change using WGS for all the samples?

library(patchwork)
library(tidyverse)
#library(plotrix)
source("../functions.R")

outdir="../../../data/results/phagocytosis/1.check_line_proportions"
dir.create(file.path(outdir),recursive = TRUE)
sample_info=read.csv("../../../data/all_pools_phagocytosis_sample_info.csv")
donor_info=read_csv("../../../../OTAR2065_differentiation_efficiency/data/donor_metadata_complete_with_imputed_sex.csv") %>%
  dplyr::mutate(sex = case_when(grepl(pattern = "Female",.$Sex) ~ "Female",
                                grepl(pattern = "Male",.$Sex) ~ "Male")) %>%
  dplyr::select(donor,sex) 
w_estimates = read_w_estimates(sample_info,
                               "../../../data/w/",
                               assay = "phagocytosis",
                               pools = paste0("pool", c(3:11,13:17)))

# convert slightly negative values to 0
w_estimates[w_estimates$proportion<0, "proportion"] = 0
w_estimates = w_estimates %>%
  dplyr::rename(line=donor) %>%
  dplyr::mutate(donor = case_when(str_detect(line,"_") ~ str_split_i(line,pattern = "_",i=1),
                                  str_detect(line,"-") ~ str_split_i(line,pattern = "-",i=1),
                                  .default = line))
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

# check differences in proportions by seeding density pool 16
# also scaled fractions
# 
# seeding_dens = w_estimates_adjusted %>%
#   dplyr::filter(pool=="pool16") %>%
#   dplyr::select(line,prop_adjusted_mean,well, mCherry, seeding_density) %>%
#   tidyr::pivot_wider(names_from = c(well), values_from = prop_adjusted_mean) %>%
#   dplyr::group_by(line,top_or_bottom) %>%
#   mutate( mymean = mean(c(`406_pellet` ,`406_swab`) )) %>% 
#   arrange(mymean) %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(line=factor(line, unique(line)))
# 
# 
# ggplot(pellet_vs_swab) +
#   geom_segment( aes(x=line, xend=line, y=`406_pellet`, yend=`406_swab`), color="grey") +
#   geom_point( aes(x=line, y=`406_pellet`, color="darkblue"), size=3 ) +
#   geom_point( aes(x=line, y=`406_swab`,color="red"), size=3 ) +
#   coord_flip()+
#   theme_minimal() +
#   ggtitle("Pellet vs swab") +
#   xlab("") +
#   ylab("Adjusted donor proportion")

# plotting lines from higher to lower prop

ordered_samples = w_estimates_adjusted%>%
  dplyr::group_by(line) %>%
  dplyr::mutate(max_prop = max(proportion)) %>%
  dplyr::ungroup()

ordered_samples$line = factor(ordered_samples$line,
                              levels = rev(unique(ordered_samples[order(ordered_samples$max_prop),"line"])$line),
                              ordered = T)
R5_mCherryNeg20=rev(unique(ordered_samples[order(ordered_samples$max_prop),"line"])$line)[1:20]
pdf(paste0(outdir,"/unsorted_neg_pos_line_proportions_col_by_mCherry_R5_mCherryNeg20.pdf"),
    width = 13, height = 10 )
ordered_samples %>%
  dplyr::filter(line %in% R5_mCherryNeg20) %>%
  ggplot( aes(x=line,y=proportion, color=mCherry)) +
  geom_boxplot(binaxis='y', stackdir='center', alpha = 0.6) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust =1)) +
  scale_color_manual(values=c("#EDA4BD","#1d3557","#9d0208","red")) +
  ylab("line proportions") +
  ggtitle("Proportions per bin per line") +
  coord_flip()
dev.off()

### check with lines that join same pool and rep
ordered_samples = ordered_samples %>%
  dplyr::mutate(groups_lines=paste(pool,treatment,replicate,sep = "_"),
                groups_dots=paste(pool,mCherry))
# check sh5y5y
pdf(paste0(outdir,"/unsorted_neg_pos_line_proportions_col_by_mCherry_sh5y5y.pdf"),
    width = 10, height = 5 )

p = ordered_samples %>% 
  dplyr::filter(line == "sh5y5y") %>%
  ggplot( aes(x=groups_dots,y=proportion,color=mCherry)) +
  scale_color_manual(values=c("#EDA4BD","#1d3557","#9d0208","red")) +
  geom_boxplot( alpha = 0.6  ) +
  geom_point(alpha = 0.3) + # Adding individual dots with slight jitter
  geom_line(aes(group = groups_lines), alpha = 0.5,col = "black") + # Connecting dots from the same pool
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust =1)) +
  ylab("line proportions") +
  ggtitle("SH5Y5Y proportions per pool and bin",subtitle = "Grey lines grouped by pool,treatment,replicate.Unsorted mCherry samples have no SH5Y5Y.") +
  scale_x_discrete(breaks =unique(ordered_samples$groups_dots), labels = c(
                                                     "", "", "pool3",
                                                     "", "", "pool4",
                                                     "", "", "pool5",
                                                     "", "", "pool6",
                                                     "", "", "pool7",
                                                     "", "", "pool8",
                                                     "", "","", "pool9",
                                                     "", "","", "pool10",
                                                     "", "", "pool11",
                                                     "", "", "pool13",
                                                     "", "", "pool14", 
                                                     "", "", "pool15",
                                                     "","","pool16",
                                                     "","","pool17"
                                                     )) # Set custom labels
plot(p)

dev.off()

pdf(paste0(outdir,"/unsorted_neg_pos_line_proportions_col_by_mCherry_sh5y5y_no_weird_samples.pdf"),
    width = 10, height = 5 )

ordered_samples2= ordered_samples %>% 
  dplyr::filter(line == "sh5y5y" & mCherry!="unsorted?")
p = ordered_samples2 %>%
  ggplot( aes(x=groups_dots,y=proportion,color=mCherry)) +
  scale_color_manual(values=c("#EDA4BD","#1d3557","#9d0208")) +
  geom_boxplot( alpha = 0.6  ) +
  geom_point(alpha = 0.3) + # Adding individual dots with slight jitter
  geom_line(aes(group = groups_lines), alpha = 0.5,col = "black") + # Connecting dots from the same pool
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust =1)) +
  ylab("line proportions") +
  ggtitle("SH5Y5Y proportions per pool and bin",subtitle = "Grey lines grouped by pool,treatment,replicate.Unsorted mCherry samples have no SH5Y5Y.") +
  scale_x_discrete(breaks =unique(ordered_samples2$groups_dots), labels = c(
    "", "", "pool3",
    "", "", "pool4",
    "", "", "pool5",
    "", "", "pool6",
    "", "", "pool7",
    "", "", "pool8",
    "", "", "pool9",
    "", "", "pool10",
    "", "", "pool11",
    "", "", "pool13",
    "", "", "pool14", 
    "", "", "pool15",
    "","","pool16",
    "","","pool17"
  )) # Set custom labels
plot(p)

dev.off()

# are sh5y5y proportions within error?
sh5y5y_check = ordered_samples %>% 
  dplyr::filter(line == "sh5y5y" & mCherry == "unsorted") %>%
  dplyr::arrange(desc(proportion)) %>%
  dplyr::mutate(within_error_max = proportion - max_wdif) %>%
  dplyr::filter(within_error_max > 0)
sh5y5y_check$well
sh5y5y_check$within_error_max

# PhagoAssay7_7A-1 not within maximum estimated error for this proportion
# for all donors: multi-page pdf

pdf(paste0(outdir,"/unsorted_neg_pos_line_proportions_col_by_mCherry_one_line_per_page.pdf"),
    width = 10, height = 5 )

for(lin in levels(ordered_samples$line)){
  p1 = ordered_samples %>% 
    dplyr::filter(line == lin) %>%
    ggplot( aes(x=groups_dots,y=proportion,color=mCherry)) +
    scale_color_manual(values=c("#EDA4BD","#1d3557","#9d0208","red")) +
    geom_boxplot( alpha = 0.6  ) +
    geom_point(alpha = 0.3) + # Adding individual dots with slight jitter
    geom_line(aes(group = groups_lines), alpha = 0.5,col = "black") + # Connecting dots from the same pool
    theme_bw() +
    theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust =1)) +
    ylab("line proportions") +
    ggtitle(paste0(lin," proportions per pool and bin"),
                   subtitle = "Grey lines grouped by pool,treatment,replicate.") +
    scale_x_discrete(breaks =unique(ordered_samples$groups_dots), labels = c(
      "", "", "pool3",
      "", "", "pool4",
      "", "", "pool5",
      "", "", "pool6",
      "", "", "pool7",
      "", "", "pool8",
      "", "","", "pool9",
      "", "","", "pool10",
      "", "", "pool11",
      "", "", "pool13",
      "", "", "pool14", 
      "", "", "pool15",
      "","","pool16",
      "","","pool17"
    )) # Set custom labels
  p2 = ordered_samples %>% 
    dplyr::filter(line == lin) %>%
    ggplot( aes(x=mCherry,y=proportion,color=pool, shape=replicate)) +
    #scale_color_manual(values=c("#2D82B7","#42E2B8","#AA4586")) +
    #geom_boxplot( alpha = 0.6  ) +
    geom_point() + # Adding individual dots with slight jitter
    geom_line(aes(group = groups_lines, linetype=treatment), alpha = 0.5) + # Connecting dots from the same pool
    theme_bw() +
    theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust =1)) +
    ylab("line proportions") 
  p=p1+p2
  print(p)
  
}
dev.off()

##### The changes in proportions should all go in same direction irrespective of the pool ###
### examining pool effects in directionality and variance ####
# see 1.3.more_checks_on_batches.R

##### examining the fractions ######
line_prop_changes = w_estimates_adjusted %>%
  dplyr::group_by(line,well) %>%
  dplyr::reframe(log_fraction_mean = log(prop_adjusted_mean[mCherry == "R3_mCherryPos"] / prop_adjusted_mean[mCherry == "R5_mCherryNeg"]),
                 log_fraction_median = log(prop_adjusted_median[mCherry == "R3_mCherryPos"] / prop_adjusted_median[mCherry == "R5_mCherryNeg"]),
                 log_fraction_max = log(prop_adjusted_max[mCherry == "R3_mCherryPos"] / prop_adjusted_max[mCherry == "R5_mCherryNeg"]),
                 prop_unadjusted_max_value = max(proportion),
                 prop_unadjusted_min_value = min(proportion),
                 pool=pool,condition=condition,treatment=treatment,replicate = replicate, donor = donor) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()

# testing inverse normal transformation
# done per well, so that ranks are done per line within well and not within whole population
line_prop_changes_INT = w_estimates_adjusted %>%
  dplyr::group_by(line,well) %>%
  dplyr::reframe(fraction_mean = prop_adjusted_mean[mCherry == "R3_mCherryPos"] / prop_adjusted_mean[mCherry == "R5_mCherryNeg"],
                 fraction_median =  prop_adjusted_median[mCherry == "R3_mCherryPos"] / prop_adjusted_median[mCherry == "R5_mCherryNeg"],
                 pool=pool,condition=condition,treatment=treatment,replicate = replicate, donor = donor) %>%
  dplyr::ungroup() %>%
  dplyr::distinct() %>%
  dplyr::group_by(well) %>%
  dplyr::mutate(
                 INT_fraction_mean = qnorm((rank(fraction_mean,
                                                 na.last="keep")-0.5)/sum(!is.na(fraction_mean))),
                 INT_fraction_median = qnorm((rank(fraction_median,
                                                   na.last="keep")-0.5)/sum(!is.na(fraction_median))),
                 
                 pool=pool,condition=condition,treatment=treatment,replicate = replicate) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()

# the values to center against are the mean log(fraction) values for hegp3, per well
centering_values = line_prop_changes %>%
  dplyr::filter(line == "hegp_3") %>%
  dplyr::group_by(well) %>%
  dplyr::summarise(
    centering_value = log_fraction_mean
  ) %>%
  ungroup()

# centering
line_prop_changes_scaled_hegp3= line_prop_changes %>%
  dplyr::left_join(centering_values) %>%
  dplyr::mutate(log_fraction_centered = NA)

line_prop_changes_scaled_hegp3$log_fraction_centered = line_prop_changes_scaled_hegp3$log_fraction_mean - line_prop_changes_scaled_hegp3$centering_value
# scaling
line_prop_changes_scaled_hegp3 = line_prop_changes_scaled_hegp3 %>%
  dplyr::filter(!log_fraction_centered %in% c("-Inf","Inf","NA","NaN")) %>%
  dplyr::group_by(well) %>%
  dplyr::mutate(scaled_log_fraction = scale(log_fraction_centered,scale = TRUE, center = FALSE)) %>%
  dplyr::ungroup()

# centering to the average of the well and scaling 
line_prop_changes_scaled_well = line_prop_changes %>%
  dplyr::filter(!log_fraction_mean %in% c("-Inf","Inf","NA","NaN")) %>%
  dplyr::group_by(well) %>%
  dplyr::mutate(scaled_log_fraction = scale_this(log_fraction_mean)) %>%
  dplyr::ungroup()


# centering to the average of the pool and scaling 
line_prop_changes_scaled_pool = line_prop_changes %>%
  dplyr::filter(!log_fraction_mean %in% c("-Inf","Inf","NA","NaN")) %>%
  dplyr::group_by(pool,treatment) %>%
  dplyr::mutate(scaled_log_fraction = scale_this(log_fraction_mean)) %>%
  dplyr::ungroup()


# beware INF/-INF fractions are valid for INT and get transformed to extreme quantiles
## histograms

p1 = line_prop_changes %>%
  dplyr::filter(!log_fraction_mean %in% c("-Inf","Inf","NA","NaN")) %>%
  ggplot(.,aes(x=log_fraction_mean)) + 
  geom_histogram(bins = 50) +
  theme_minimal() + 
  ggtitle(" log fraction mean: phagocytosis ")

p2 =  line_prop_changes %>%
  dplyr::filter(!log_fraction_median %in% c("-Inf","Inf","NA","NaN")) %>%
  ggplot(.,aes(x=log_fraction_median)) + 
  geom_histogram(bins = 50) +
  theme_minimal() + 
  ggtitle(" log fraction median: phagocytosis ")

p3 = ggplot(line_prop_changes_INT,aes(x=INT_fraction_mean)) + 
  geom_histogram(bins = 50) +
  theme_minimal() + 
  ggtitle("INT fraction mean error: phagocytosis")

p4 = ggplot(line_prop_changes_INT,aes(x=INT_fraction_median)) + 
  geom_histogram(bins = 50) +
  theme_minimal() + 
  ggtitle("INT fraction median error:: phagocytosis")

p5 = ggplot(line_prop_changes_scaled_hegp3,aes(x=scaled_log_fraction)) + 
  geom_histogram(bins = 50) +
  theme_minimal() + 
  ggtitle("scaled and centered (s&c)", subtitle = "(hegp3) log fraction: phagocytosis")

p6 = ggplot(line_prop_changes_scaled_well,aes(x=scaled_log_fraction)) + 
  geom_histogram(bins = 50) +
  theme_minimal() + 
  ggtitle("s&c (well mean) log fraction:  phagocytosis")

p7 = ggplot(line_prop_changes_scaled_pool,aes(x=scaled_log_fraction)) + 
  geom_histogram(bins = 50) +
  theme_minimal() + 
  ggtitle("s&c (pool mean) log fraction: phagocytosis")

pdf (paste0(outdir,"/histograms_log_INT_fraction_mean_median_errors.pdf"),
     width = 13, height = 6 )
(p1+ p2 + p3 +p4) / (p5 + p6 + p7 + patchwork::plot_spacer())
dev.off()


### what proportion of donors are 2sd away from the mean ###
summary_stats = line_prop_changes_scaled_well %>%
  summarize(
    mean_log = mean(log_fraction_mean),
    sd_log = sd(log_fraction_mean),
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
(outliers_count_log / nrow(line_prop_changes_scaled_well)) * 100
(outliers_count_scaled_well / nrow(line_prop_changes_scaled_well)) * 100
# there are fewer outliers (defined as data points > 3 sd from mean) in the scaled data (0.88%) vs log data (1.75%)
##
pdf (paste0(outdir,"/change_in_total_R5_mCherryNeg_R3_mCherryPos_line_proportions_per_line_mean_error.pdf"),
     width = 13, height = 9 )
line_prop_changes %>%
  
  ggplot( aes(x=line,y=log_fraction_mean,color = replicate,
              size = prop_unadjusted_max_value)) +
  geom_point(binaxis='y', stackdir='center', alpha = 0.6) +
  theme_bw() +
  scale_color_brewer(palette  = "Set1") +
  theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust =1)) +
  ylab("log(R3_mCherryPos prop / R5_mCherryNeg prop)") +
  ggtitle("log(R3_mCherryPos prop / R5_mCherryNeg prop) [adjust: mean error]")

dev.off()

pdf (paste0(outdir,"/change_in_total_R5_mCherryNeg_R3_mCherryPos_line_proportions_per_line_mean_error_boxplot.pdf"),
     width = 7, height = 14 )
ordered = line_prop_changes %>%
  dplyr::arrange(desc(prop_unadjusted_max_value)) 
ordered$line = factor(ordered$line,levels = unique(ordered$line), ordered = TRUE)  
ordered %>%
  ggplot( aes(x=line,y=log_fraction_mean,col=treatment
  )) +
  geom_boxplot(binaxis='y', stackdir='center', alpha = 0.6) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust =1)) +
  ylab("log(R3_mCherryPos prop / R5_mCherryNeg prop)") +
  ggtitle("log (R3_mCherryPos prop / R5_mCherryNeg prop) [adjust: mean error]") +
  coord_flip()

dev.off()
pdf (paste0(outdir,"/change_in_total_R5_mCherryNeg_R3_mCherryPos_line_proportions_per_line_mean_error_boxplot_INT.pdf"),
     width = 7, height = 14 )
line_prop_changes_INT %>%
  
  ggplot( aes(x=line,y=INT_fraction_mean,col=treatment
              )) +
  geom_boxplot(binaxis='y', stackdir='center', alpha = 0.6) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust =1)) +
  ylab("INT(R3_mCherryPos prop / R5_mCherryNeg prop)") +
  ggtitle("INT (R3_mCherryPos prop / R5_mCherryNeg prop) [adjust: mean error]") +
  coord_flip()

dev.off()

## check INT transformation after taking the log

line_prop_changes %>%
  dplyr::group_by(well) %>%
  dplyr::mutate( INT_log_fraction_mean = qnorm((rank(log_fraction_mean,
                                    na.last="keep")-0.5)/sum(!is.na(log_fraction_mean)))) %>%
  ggplot( aes(x=line,y=INT_log_fraction_mean )) +
  geom_boxplot(binaxis='y', stackdir='center', alpha = 0.6) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust =1)) +
  ylab("INT(log(R3_mCherryPos prop / R5_mCherryNeg prop))") +
  ggtitle("INT (log(R3_mCherryPos prop / R5_mCherryNeg prop)) [adjust: mean error]") +
  coord_flip()

# same as before the log

## add sex info
line_prop_changes_scaled_well = line_prop_changes_scaled_well %>%
  dplyr::left_join(.,donor_info)

line_prop_changes_scaled_well %>%
  write.csv(paste0(outdir,"/line_prop_changes_per_well.csv"),
            quote = FALSE,row.names = FALSE)

### is there any particular pool that skews the values?

pdf (paste0(outdir,"/log_fractions_per_pool_treatment_line.pdf"),
     width = 7, height = 4 )
ordered = line_prop_changes %>%
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
      ggplot(aes(x = pool,y = log_fraction_mean, col = treatment)) +
      geom_point() +
     geom_boxplot()+
      geom_hline(yintercept = 0,linetype = "dotted", col="grey") + 
      theme_minimal() + 
      ggtitle(lin)
   plot(p)
  }
}
  dev.off()
  
  pdf (paste0(outdir,"/scaled_centered_hegp3_log_fractions_per_pool_treatment_line.pdf"),
       width = 7, height = 4 )
  ordered = line_prop_changes_scaled_hegp3 %>%
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
  
  
  pdf (paste0(outdir,"/scaled_centered_well_average_log_fractions_per_pool_treatment_line.pdf"),
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
  
  
  pdf (paste0(outdir,"/INT_fractions_per_pool_treatment_line.pdf"),
       width = 7, height = 4 )
  ordered_INT = line_prop_changes_INT %>%
    dplyr::left_join(ordered[,c("line","well","prop_unadjusted_max_value")]) %>%
    dplyr::arrange(desc(prop_unadjusted_max_value)) 
  ordered_INT$line = factor(ordered_INT$line,levels = unique(ordered_INT$line), ordered = TRUE) 
  for(lin in unique(ordered_INT$line)){
    subset = ordered_INT %>%
      dplyr::filter(line %in% lin)
    if(sum(is.na(subset$log_fraction_mean)) == nrow(subset)){
      next()
    }else{
      subset$pool = factor(subset$pool,levels = sort(unique(subset$pool)), ordered = TRUE)
      subset$treatment = factor(subset$treatment,levels = c("untreated","IFN","LPS"), ordered = TRUE)
      p =  subset %>%
        ggplot(aes(x = pool,y = INT_fraction_mean, col = treatment)) +
        geom_point() +
        geom_boxplot()+
        geom_hline(yintercept = 0,linetype = "dotted", col="grey") + 
        theme_bw() + 
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
 
 pdf (paste0(outdir,"/mean_stdev_scaled_fractions_per_treatment_max_prop.pdf"),
      width = 7, height = 7 )
 p1 / p2
 
 dev.off()

# C5a should show more phagocytosis to R3_mCherryPos well than media only.
# where is the phagocytosis info for all pools? How to incorporate this info to the model?

# now, averaging wells and then averaging pools, for those that are shared ####
###########

# number of lines before filters:
length(unique(line_prop_changes_scaled_well$line)) # 260 pools 3-17
# with NA, INF filters
line_prop_changes_scaled_well %>%
  # first group by rep, and calculate rep mean and variance
  dplyr::group_by(line,replicate,treatment,pool) %>%
  dplyr::filter(!scaled_log_fraction %in% c(NA,NaN,Inf,-Inf)) %>%
  dplyr::ungroup() %>%
  summarize(unique_elements = n_distinct(line)) # 202 pools 3-17
# with 1% filters
line_prop_changes_scaled_well %>%
  # first group by rep, and calculate rep mean and variance
  dplyr::group_by(line,replicate,treatment,pool) %>%
  dplyr::filter(!scaled_log_fraction %in% c(NA,NaN,Inf,-Inf)) %>%
  dplyr::filter(prop_unadjusted_max_value >= 0.01) %>% 
  dplyr::ungroup() %>%
  summarize(unique_elements = n_distinct(line)) # 143 pools 3-17

# with 0.5% filters
line_prop_changes_scaled_well %>%
  # first group by rep, and calculate rep mean and variance
  dplyr::group_by(line,replicate,treatment,pool) %>%
  dplyr::filter(!scaled_log_fraction %in% c(NA,NaN,Inf,-Inf)) %>%
  dplyr::filter(prop_unadjusted_max_value >= 0.005) %>%
  dplyr::ungroup() %>%
  summarize(unique_elements = n_distinct(line)) # 174 pools 3-17

line_prop_changes_scaled_well %>%
  # first group by rep, and calculate rep mean and variance
  dplyr::group_by(line,treatment,pool) %>%
  dplyr::filter(!scaled_log_fraction %in% c(NA,NaN,Inf,-Inf)) %>%
  dplyr::filter(prop_unadjusted_max_value >= 0.01) %>%
  dplyr::summarise(mean_outcome_rep = mean(scaled_log_fraction,na.rm = TRUE),
                   var_outcome_rep = var(scaled_log_fraction,na.rm = TRUE),
                   mean_max_prop_rep = mean(prop_unadjusted_max_value,na.rm = TRUE),
                   mean_min_prop_rep = mean(prop_unadjusted_min_value,na.rm = TRUE),
                   var_max_prop_rep = var(prop_unadjusted_max_value,na.rm = TRUE)) %>%
  ungroup() %>%
  
  # then, calculate mean and variance across pools 
  # variance is variance of the means
  dplyr::group_by(line,treatment) %>% 
  # calculating mean of means, and var of means
  dplyr::mutate(mean_outcome_pool = mean(mean_outcome_rep,na.rm = TRUE),
                var_outcome_pool = var(mean_outcome_rep,na.rm = TRUE),
                mean_max_prop_pool = mean(mean_max_prop_rep,na.rm = TRUE),
                mean_min_prop_pool = mean(mean_min_prop_rep,na.rm = TRUE),
                var_max_prop_pool = var(mean_max_prop_rep,na.rm = TRUE),
                n_pools = n()) %>%
  ungroup() %>%
  write.csv(paste0(outdir,"/line_prop_changes_averages_variances_1pct_filtered.csv"),
            quote = FALSE,row.names = FALSE)
# same for 0.5%
line_prop_changes_scaled_well %>%
  # first group by rep, and calculate rep mean and variance
  dplyr::group_by(line,treatment,pool) %>%
  dplyr::filter(!scaled_log_fraction %in% c(NA,NaN,Inf,-Inf)) %>%
  dplyr::filter(prop_unadjusted_max_value >= 0.005) %>%
  dplyr::summarise(mean_outcome_rep = mean(scaled_log_fraction,na.rm = TRUE),
                   var_outcome_rep = var(scaled_log_fraction,na.rm = TRUE),
                   mean_max_prop_rep = mean(prop_unadjusted_max_value,na.rm = TRUE),
                   mean_min_prop_rep = mean(prop_unadjusted_min_value,na.rm = TRUE),
                   var_max_prop_rep = var(prop_unadjusted_max_value,na.rm = TRUE)) %>%
  ungroup() %>%
  
  # then, calculate mean and variance across pools 
  # variance is variance of the means
  dplyr::group_by(line,treatment) %>% 
  # calculating mean of means, and var of means
  dplyr::mutate(mean_outcome_pool = mean(mean_outcome_rep,na.rm = TRUE),
                var_outcome_pool = var(mean_outcome_rep,na.rm = TRUE),
                mean_max_prop_pool = mean(mean_max_prop_rep,na.rm = TRUE),
                mean_min_prop_pool = mean(mean_min_prop_rep,na.rm = TRUE),
                var_max_prop_pool = var(mean_max_prop_rep,na.rm = TRUE),
                n_pools = n()) %>%
  ungroup() %>%
  write.csv(paste0(outdir,"/line_prop_changes_averages_variances_halfpct_filtered.csv"),
            quote = FALSE,row.names = FALSE)

# checking relationship of variance from pQTL outcome (log of fraction) with min or max or average proportion in the fraction
# in the well replicates, and also in the pool replicates from shared lines
# to try to inform a reasonable threshold for prefiltering based on proportion of the line

# calculate variance between replicates
line_prop_changes_shared = line_prop_changes %>%
  dplyr::group_by(line,treatment,pool) %>%
  dplyr::filter(!log_fraction_mean %in% c(NA,NaN,Inf,-Inf)) %>%
  # remember in log_fraction_mean, the "mean" refers to mean error used to adjust the proportions
  dplyr::summarise(mean_outcome = mean(log_fraction_mean),
                   var_outcome = var(log_fraction_mean),
                   mean_max_prop = mean(prop_unadjusted_max_value),
                   var_max_prop = var(prop_unadjusted_max_value),
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
  facet_wrap(vars(treatment), ncol = 3)

pdf(paste0(outdir,"/mean_vs_var_all_lines_by_well_reps.pdf"),
    width = 10, height = 5)
p
dev.off()

# now checking for shared lines, variation across pools
line_prop_changes_shared  = line_prop_changes %>%
  dplyr::group_by(line,treatment) %>%
  dplyr::filter(!log_fraction_mean %in% c(NA,NaN,Inf,-Inf)) %>%
  # Filter for lines that are present in more than one pool
  dplyr::filter(n_distinct(pool) > 1) %>%
  dplyr::ungroup()

# first group by rep, and calculate rep mean and variance
line_prop_changes_shared =  line_prop_changes_shared %>%
  dplyr::group_by(line,treatment,pool) %>%
  dplyr::filter(!log_fraction_mean %in% c(NA,NaN,Inf,-Inf)) %>%
  # remember in log_fraction_mean, the "mean" refers to mean error used to adjust the proportions
  dplyr::summarise(mean_outcome_rep = mean(log_fraction_mean),
                   var_outcome_rep = var(log_fraction_mean),
                   mean_max_prop_rep = mean(prop_unadjusted_max_value),
                   var_max_prop_rep = var(prop_unadjusted_max_value)) %>%
  ungroup() %>%
  
  # then, calculate mean and variance across pools 
  # variance is variance of the means
  dplyr::group_by(line,treatment) %>% 
  # calculating mean of means, and var of means
  dplyr::summarise(mean_outcome_pool = mean(mean_outcome_rep),
                   var_outcome_pool = var(mean_outcome_rep),
                   mean_max_prop_pool = mean(mean_max_prop_rep),
                   var_max_prop_pool = var(mean_max_prop_rep),
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
  ggtitle("Variance of log(R3_mCherryPos prop/R5_mCherryNeg prop) in shared lines vs maximum proportion", subtitle = "mean(mean[max prop] across rep, across pools) vs var(mean[fractions] per rep, across pools)") +
  facet_wrap(vars(treatment), ncol = 3)

pdf(paste0(outdir,"/mean_vs_var_shared_lines_across_pools.pdf"),
    width = 10, height = 6)
p
dev.off()

# check variation between pools
line_prop_changes_shared = line_prop_changes %>%
  dplyr::group_by(line,treatment) %>%
  dplyr::filter(!log_fraction_mean %in% c(NA,NaN,Inf,-Inf)) %>%
  # Filter for lines that are present in more than one pool
  dplyr::filter(n_distinct(pool) > 1) %>%
  dplyr::ungroup() %>%
  # calculate mean and variance for pools and wells, per line, replicate and treatment
  dplyr::group_by(line,treatment,pool) %>% 
  # remember in log_fraction_mean, the "mean" refers to mean error used to adjust the proportions
  dplyr::summarise(mean_outcome = mean(log_fraction_mean),
                   var_outcome = var(log_fraction_mean),
                   mean_max_prop = mean(prop_unadjusted_max_value),
                   var_max_prop = var(prop_unadjusted_max_value),
                   CV_max_prop = sqrt(var(prop_unadjusted_max_value)) / mean(prop_unadjusted_max_value))

p = line_prop_changes_shared %>%
  group_by(line,pool) %>%
  dplyr::filter(!mean_max_prop %in% c(0,NA,Inf,-Inf)) %>%
  # remove 0 to avoid inf after log
  ggplot( aes(x=log10(mean_max_prop),y=var_outcome,col=pool)) +
  geom_point(aes(size = CV_max_prop), alpha = 0.6) +
  geom_vline(xintercept = -2, col="grey", linetype="dashed") +
  theme_bw() +
  ylab("Variance of fractions") +
  xlab("log10(Mean of maximum proportion) -> larger prop. to right") +
  ggtitle("Mean(max prop) vs var(fractions) ", subtitle = "for all replicates, per line and pool, in shared lines") +
  facet_wrap(vars(treatment), ncol = 3)
pdf(paste0(outdir,"/mean_vs_var_all_wells_per_pool.pdf"),
    width = 10, height = 5)
p
dev.off()

# check variation within pools, between wells
# fix the renaming of wells to replicates, it's not working now

# line_prop_changes_shared = line_prop_changes %>%
#   dplyr::group_by(line,replicate,treatment) %>%
#   dplyr::filter(!log_fraction_mean %in% c(NA,Inf,-Inf)) %>%
#   # Filter for lines that are present in more than one pool
#   dplyr::filter(n_distinct(pool) > 1) %>%
#   dplyr::ungroup() %>%
#   # calculate mean and variance for pools and wells, per line, replicate and treatment
#   dplyr::group_by(line,pool,treatment, replicate, well) %>% 
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
#   facet_wrap(vars(replicate,treatment), ncol = 3)
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
# pdf ("../../data/results/4.check_donor_proportions/change_in_total_R5_mCherryNeg_R3_mCherryPos_line_proportions_per_line_sim_error.pdf ",
#     width = 13, height = 10 )
# ordered_samples %>%
#   ggplot( aes(x=cells_format_assay,y=proportion, color=mCherry)) +
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
# pdf ("../../data/results/4.check_donor_proportions/change_in_total_R5_mCherryNeg_R3_mCherryPos_line_proportions_per_well_sim_error.pdf ",
#     width = 13, height = 9 )
# samples %>%
#   ggplot( aes(x=line,y=proportion, color=mCherry)) +
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
# pdf ("../../data/results/4.check_donor_proportions/change_in_total_line_proportions_per_line.pdf ",
#     width = 11, height = 10 )
# ordered_samples %>%
#   filter(mCherry == "Total") %>%
#   ggplot( aes(x=Cells,y=proportion, color=Assay_replicate, shape = Format)) +
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
# pdf ("../../data/results/4.check_donor_proportions/change_in_total_line_proportions_per_line_sim_errors.pdf ",
#     width = 11, height = 10 )
# ordered_samples %>%
#   filter(mCherry == "Total") %>%
#   ggplot( aes(x=Cells,y=proportion, color=Assay_replicate, shape = Format)) +
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
#   filter(mCherry == "Total") %>%
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
#   filter(mCherry == "Total") %>%
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
#   filter(mCherry == "Total") %>%
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
# pdf ("../../data/results/4.check_donor_proportions/mean_prop_per_line_per_IFN_per_format.pdf ",
#     width = 7, height = 4 )
# samples %>%
#   filter(mCherry == "Total") %>%
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
# pdf ("../../data/results/4.check_donor_proportions/mean_prop_per_line_per_IFN_per_format_se.pdf ",
#     width = 7, height = 4 )
# samples %>%
#   filter(mCherry == "Total") %>%
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
#   filter(mCherry == "Total") %>%
#   select(line, proportion, cells_format_assay)
# colnames(new_WGS) = c("line", "proportion", "experiment")
# 
# 
# # 89 and 101 are exactly the same?
# samples %>% filter(well == "89" | well == "101")
# # not exactly
# 
# 
# # IFN stimulates macrophage phagocytosis
# # Test DNA in R3_mCherryPos well for IFN vs no IFN, and within same IFN status the effect of C5a + inhibitors
# pdf ("../../data/results/4.check_donor_proportions/DNA_in_R3_mCherryPos_well_over_total_boxplot_replicates.pdf ",
#     width = 7, height = 4 )
# samples %>%
#   filter(mCherry == "R3_mCherryPos", line == "letw_1" ) %>%
#   ggplot( aes(x=Cells,y=DNA_normalised_to_total)) +
#   geom_boxplot(color="grey", alpha=0.2) +
#   geom_jitter(aes(color = Assay_replicate),shape=16, position=position_jitter(0.2),size = 3) +
#   theme_bw() +
#   theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust =1)) +
#   scale_color_manual(values=c("#ffbe0b","#003049","#8338ec","#ff006e")) +
#   ylab("DNA in R3_mCherryPos well over total")
# 
# dev.off()
# # The replicates behave as expected in the IFN case (+ inhibitors or without C5a have lower R3_mCherryPos DNA prop)
# # But not in the no_IFN case. Why?
