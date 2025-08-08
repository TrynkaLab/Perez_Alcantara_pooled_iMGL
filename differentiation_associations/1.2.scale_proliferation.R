# scale data and format for associations
# proportions have been loaded first for plotting in alluvial plots
library(patchwork)
library(tidyverse)
source("./functions.R")

input_dir = "../../data/results/1.alluvial_plots"
output_dir = "../../data/results/1.2.scale_proliferation"
dir.create(output_dir, recursive = TRUE)
w_estimates = readr::read_tsv(paste0(input_dir,
                         "/pools2-11_13-17_changing_props_iPSC_preMacs_microglia_WGS_sc.txt")) %>%
  dplyr::rename(line = Line,proportion = prop) %>%
  dplyr::mutate(donor = if_else(str_detect(line,"-"),
                                true=str_split_i(line,pattern="-",i=1), 
                                false=str_split_i(line,pattern="_",i=1))) %>%
  dplyr::relocate(line,donor)
donor_info=read_csv("../../../OTAR2065_differentiation_efficiency/data/donor_metadata_complete_with_imputed_sex.csv") %>%
  dplyr::rename(sex=Sex) %>%
  dplyr::select(donor,sex)

sample_info = readr::read_csv("../../data/OTAR2065_pools_runs_samples_precursor_info.csv") %>%
  dplyr::select(sample:treatment)

w_estimates = w_estimates %>%
  dplyr::left_join(donor_info, by="donor") %>%
  dplyr::left_join(sample_info, by = c("sample","pool","stage"))
##### adjusting WGS proportions ########

# convert slightly negative values to 0
w_estimates[w_estimates$proportion<0, "proportion"] = 0
write.table(w_estimates,paste0(output_dir,"/w_estimates.txt"),
            row.names = FALSE,quote = FALSE)
# read in proportion error estimates
error_estimates = read.table("../../../OTAR2065_phenotypic_QTLs/data/w/error_approximations_generic_pool.txt", header = TRUE) %>%
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
  dplyr::mutate(prop_adjusted_mean = dplyr::case_when(proportion <0.1 & type == "WGS"  ~ proportion - mean_wdif,
                                                      proportion >=0.1  & type == "WGS"  ~ proportion + mean_wdif,
                                                      .default = proportion), # scRNA-seq don't change proportion
                prop_adjusted_median = dplyr::case_when(proportion <0.1 & type == "WGS"  ~ proportion - median_wdif,
                                                        proportion >=0.1 & type == "WGS"  ~ proportion + median_wdif,
                                                        .default = proportion),
                prop_adjusted_max = dplyr::case_when(proportion <0.1 & type == "WGS"  ~ proportion - max_wdif,
                                                     proportion >=0.1 & type == "WGS"  ~ proportion + max_wdif,
                                                     .default = proportion),
    )
# if any estimate becomes negative, make it zero:
w_estimates_adjusted[w_estimates_adjusted$prop_adjusted_mean<0, "prop_adjusted_mean"] = 0
w_estimates_adjusted[w_estimates_adjusted$prop_adjusted_median<0, "prop_adjusted_median"] = 0
w_estimates_adjusted[w_estimates_adjusted$prop_adjusted_max<0, "prop_adjusted_max"] = 0

#####  end of adjustments #####

# plotting lines from higher to lower prop

ordered_samples = w_estimates_adjusted%>%
  dplyr::group_by(line) %>%
  dplyr::mutate(max_prop = max(proportion)) %>%
  dplyr::ungroup()

ordered_samples$line = factor(ordered_samples$line,
                              levels = rev(unique(ordered_samples[order(ordered_samples$max_prop),"line"])$line),
                              ordered = T)

length(unique(ordered_samples$line)) # 261


pdf(paste0(output_dir,"/raw_line_proportions_per_line_incl_all_condition_dotplot.pdf"),
    width = 10, height = 17 )
ordered_samples %>%
  ggplot( aes(x=line,y=proportion, color=stage)) +
  geom_point(binaxis='y', stackdir='center', alpha = 0.6) +
  theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust =1)) +
  scale_color_manual(values=c("iPSC" = "lightblue",
                              "preMac" = "orange1",
                              "microglia" = "darkred")) +
  ylab("line proportions") +
  ggtitle("Proportions per line") +
  coord_flip() + 
  theme_bw() +
  facet_wrap(vars(type))
dev.off()

pdf(paste0(output_dir,"/raw_top10_line_proportions_per_line_incl_all_condition_boxplot.pdf"),
    width = 10, height = 17 )
ordered_samples %>%
  dplyr::filter(treatment %in% c(NA,"untreated") & line %in% levels(ordered_samples$line)[1:10]) %>%
  ggplot( aes(x=line,y=proportion, color=stage)) +
  geom_boxplot()+
  theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust =1)) +
  scale_color_manual(values=c("iPSC" = "lightblue",
                              "preMac" = "orange1",
                              "microglia" = "darkred")) +
  ylab("line proportions") +
  ggtitle("Proportions per line") +
  coord_flip() + 
  theme_bw() 
dev.off()



####### scaling ########
# centering to the average of the sample and scaling
### need to do this carefully per pool for preMac / iPSC, microglia / preMac, and microglia / iPSC

### preMac vs iPSC ###
line_prop_changes_premac_iPSC = w_estimates_adjusted %>%
  dplyr::filter(stage %in% c("iPSC","preMac")) %>%
  # selecting the youngest preMacs for the comparison
  dplyr::group_by(pool, stage,line,differentiation) %>%
  dplyr::slice_min(preMAC_age) %>%
 dplyr::ungroup() %>%
 dplyr::group_by(line,pool, differentiation) %>% # grouping also per differentiation because samples in pool2 have several levels of this
  dplyr::reframe(log_fraction_mean = log(prop_adjusted_mean[stage == "preMac"] / prop_adjusted_mean[stage == "iPSC"]),
                 
                 prop_unadjusted_max_value = max(proportion),
                 prop_unadjusted_min_value = min(proportion),
                 comparison = "youngest_preMac_vs_iPSC",
                 preMAC_age = na.omit(preMAC_age),
                 pool=pool,type = type, sex = sex) %>%
  dplyr::ungroup() %>%
  dplyr::distinct() %>%
  dplyr::filter(!log_fraction_mean %in% c("-Inf","Inf","NA","NaN")) %>% # removing incorrect fractions 
  # centering to the average of the pool, differentiation and preMac age groups (and scaling within same grouping)
  dplyr::group_by(pool,differentiation,preMAC_age) %>%
  dplyr::mutate(scaled_log_fraction = scale_this(log_fraction_mean)) %>%
  dplyr::ungroup()
hist(line_prop_changes_premac_iPSC$scaled_log_fraction) # normal

### microglia vs preMac ###

# for microglia, need to do this per treatment as well
# expanding preMacs categories for IFNg and LPS (treatment), 
# and for scRNAseq (type) to facilitate the comparisons
expanded_premacs_untreated = w_estimates_adjusted %>%
  dplyr::filter(stage %in% c("preMac")) %>%
  dplyr::mutate(treatment = "IFNg",type = "scRNA-seq")

expanded_premacs_IFN = w_estimates_adjusted %>%
  dplyr::filter(stage %in% c("preMac")) %>%
  dplyr::mutate(treatment = "IFNg") %>%
  dplyr::bind_rows(w_estimates_adjusted %>%
                     dplyr::filter(stage %in% c("preMac")) %>%
                     dplyr::mutate(treatment = "IFNg",type = "scRNA-seq")) # not the prettiest code
  
expanded_premacs_LPS = w_estimates_adjusted %>%
  dplyr::filter(stage %in% c("preMac")) %>%
  dplyr::mutate(treatment = "LPS") %>%
  dplyr::bind_rows(w_estimates_adjusted %>%
                     dplyr::filter(stage %in% c("preMac")) %>%
                     dplyr::mutate(treatment = "LPS",type = "scRNA-seq")) # not the prettiest code


line_prop_changes_microglia_premac = w_estimates_adjusted %>%
  dplyr::filter(stage %in% c("preMac","microglia")) %>%
  # for those microglia  pools that have several scRNA-seq samples, get the average
  dplyr::group_by(pool,preMAC_age, treatment,stage,differentiation,type,line) %>%
  dplyr::reframe(prop_adjusted_mean = mean(prop_adjusted_mean),
                 proportion = mean(proportion),
                 preMAC_age = na.omit(preMAC_age),
                 treatment = na.omit(treatment),
                 pool=pool,sex = sex) %>%
  dplyr::ungroup() %>%
  dplyr::distinct() %>%
  # add expanded preMacs
  dplyr::bind_rows(expanded_premacs_IFN) %>%
  dplyr::bind_rows(expanded_premacs_LPS) %>%
  dplyr::bind_rows(expanded_premacs_untreated) %>%
  dplyr::group_by(pool,preMAC_age,differentiation,line,treatment,type) %>%
  dplyr::reframe(log_fraction_mean = log(prop_adjusted_mean[stage == "microglia"] / prop_adjusted_mean[stage == "preMac"]),
                 
                 prop_unadjusted_max_value = max(proportion),
                 prop_unadjusted_min_value = min(proportion),
                 comparison = "microglia_vs_preMac",
                sex = unique(sex)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct() %>%
  dplyr::filter(!log_fraction_mean %in% c("-Inf","Inf","NA","NaN")) %>% # removing incorrect fractions 
  # centering to the average of the pool, differentiation, preMac age, treatment and type groups (and scaling within same grouping)
  dplyr::group_by(pool,differentiation,preMAC_age,treatment,type) %>%
  dplyr::mutate(scaled_log_fraction = scale_this(log_fraction_mean)) %>%
  dplyr::ungroup()
hist(line_prop_changes_microglia_premac$scaled_log_fraction) # normal


### microglia vs iPSC? ###
### preMac aging ###
# age continuum
line_prop_changes_aging_premac = w_estimates_adjusted %>%
  dplyr::filter(stage =="preMac") %>%
  # scale adjusted proliferation directly within sample, across all lines
  dplyr::group_by(pool,differentiation,preMAC_age) %>%
  dplyr::mutate(scaled_proportion = scale_this(prop_adjusted_mean),
                 prop_unadjusted_max_value = max(proportion),
                 prop_unadjusted_min_value = min(proportion),
                 comparison = "preMac_across_age") %>%
  dplyr::ungroup() %>%
  dplyr::distinct() %>%
  dplyr::filter(!scaled_proportion %in% c("-Inf","Inf","NA","NaN")) %>%
  # now we have a choice: if I do the following centering and scaling across premac ages, FOR ALL LINES: the extreme / very low prolif lines don't completely lose their extreme values)
  # but the effects will be driven by either the outliers or the abundance of very small changes
  # if I center and scale across premac age per line, the smaller proportion changes will be maximised and the larger will be shrunk
  # the second might be better for genome-wide analysis but the first for burden tests?
  dplyr::group_by(pool,differentiation) %>%
  dplyr::mutate(scaled_proportion = scale_this(scaled_proportion)) %>%
  dplyr::ungroup()

hist(line_prop_changes_aging_premac$scaled_proportion) # all over the place if I include "line" in the last grouping 
# for all lines looks Poisson


# old vs young
line_prop_changes_old_vs_young_premac = w_estimates_adjusted %>%
  dplyr::filter(stage =="preMac") %>%
  dplyr::group_by(pool,differentiation,line) %>%
  dplyr::slice(c(which.min(preMAC_age), which.max(preMAC_age))) %>%
  dplyr::reframe(log_fraction_mean = log(prop_adjusted_mean[preMAC_age ==max(preMAC_age)] / prop_adjusted_mean[preMAC_age ==min(preMAC_age)]),
                 
                 prop_unadjusted_max_value = max(proportion),
                 prop_unadjusted_min_value = min(proportion),
                 comparison = "preMac_old_vs_young",
                 sex = unique(sex)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct() %>%
  dplyr::filter(!log_fraction_mean %in% c("-Inf","Inf","NA","NaN")) %>% # removing incorrect fractions 
  # centering to the average of the pool, differentiation, and preMac age groups (and scaling within same grouping)
  dplyr::group_by(pool,differentiation) %>%
  dplyr::mutate(scaled_log_fraction = scale_this(log_fraction_mean)) %>%
  dplyr::ungroup()
hist(line_prop_changes_old_vs_young_premac$scaled_log_fraction) # normal


### write

write_csv(line_prop_changes_premac_iPSC,paste0(output_dir,"/line_prop_changes_premac_iPSC.csv"))
write_csv(line_prop_changes_microglia_premac,paste0(output_dir,"/line_prop_changes_microglia_premac.csv"))
write_csv(line_prop_changes_old_vs_young_premac,paste0(output_dir,"/line_prop_changes_old_vs_young_premac.csv"))
write_csv(line_prop_changes_aging_premac,paste0(output_dir,"/line_prop_changes_aging_premac.csv"))
##### examining the scaled fractions against iPSC and premac/microglia proportions ######

p1 = line_prop_changes_aging_premac %>%
  dplyr::left_join(w_estimates_adjusted) %>%
  dplyr::filter(pool %in% paste0("pool",3)) %>%
  ggplot(aes(x = line, y = scaled_proportion)) + 
    geom_jitter() +
    theme_bw()
  p2 = line_prop_changes_aging_premac %>%
    dplyr::left_join(w_estimates_adjusted) %>%
    dplyr::filter(pool %in% paste0("pool",3))  %>%
  ggplot(aes(x = line, y = prop_adjusted_mean)) + 
    geom_jitter() +
    theme_bw()

  p1 + p2

