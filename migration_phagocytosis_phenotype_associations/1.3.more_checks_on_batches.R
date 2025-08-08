# Tests on variance of fraction estimates per pool


library(patchwork)
library(tidyverse)
library(plotrix)
library(PCAtools)
library(pheatmap)
source("../functions.R")

outdir1="../../../data/results/phagocytosis/1.check_donor_proportions"
outdir2="../../../data/results/phagocytosis/1.3.more_checks_on_batches"

dir.create(file.path(outdir2),recursive = TRUE)
sample_info=read.csv("../../../data/all_pools_phagocytosis_sample_info.csv")
donor_info=read.csv("../../../data/allinfo_hipsci_PD_AD_PRS_nondisease_feederfree_european.csv") %>%
  dplyr::mutate(sex = case_when(grepl(pattern = "Female",.$Sex) ~ "Female",
                                grepl(pattern = "Male",.$Sex) ~ "Male")) %>%
  dplyr::select(Line,sex) %>%
  dplyr::rename(line=Line)
w_estimates = read_w_estimates(sample_info,
                               "../../../data/w/",
                               assay = "phagocytosis",
                               pools = paste0("pool", c(3:11,13)))

# convert slightly negative values to 0
w_estimates[w_estimates$proportion<0, "proportion"] = 0
w_estimates = w_estimates %>%
  dplyr::rename(line=donor)
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

line_prop_changes_per_well = read.csv(paste0(outdir1,"/line_prop_changes_per_well.csv")) %>%
  dplyr::mutate(pool = factor(pool,levels = paste0("pool",c(2:11,13)))) 


genotype = readr::read_csv("../../../data/genotypes/full_genotype/genotype_minor_allele_dosage.csv") %>%
  dplyr::relocate(rn) %>%
  dplyr::sample_n(100000) %>% # sample 100k positions
  tidyr::pivot_longer(cols = !rn, values_to = "alt_allele_dosage", names_to = "line") %>%
  dplyr::mutate(alt_allele_dosage = dplyr::case_when(
    alt_allele_dosage == 0.5 ~ 1,
    alt_allele_dosage == 1   ~ 2,
    TRUE          ~ alt_allele_dosage  # Recode so var of allele dosage is 1, Keep 0 unchanged
  ))

# boxplots without averaging replicates
snp = "chr20_62137982_T_C"

subset = line_prop_changes_per_well %>%
  dplyr::filter(!log_fraction_mean %in% c(NA,Inf,-Inf)) %>%
  ungroup() %>%
  left_join(.,genotype[genotype$rn ==snp,]) %>%
  dplyr::filter(treatment == "untreated")  %>%
  distinct() 
p = boxplot_pQTL_exploratory_reps(variant = snp,phenotype_with_genotype_info_df = subset,
                                  color_by="pool",
                                  phenotype = "phagocytosis")

# for donors that experience average reduction in proportions among shared donors,
# are there pool differences?

# do certain pools have lower / higher ourcomes, per genotype?

# homozygotes REF positions for shared donor hegp_3
# 
# hegp_3_fractions = line_prop_changes_per_well %>%
#   dplyr::filter(line =="hegp_3") %>% 
#   dplyr::filter(!log_fraction_mean %in% c(NA,Inf,-Inf)) %>%
#   ungroup() %>%
#   left_join(.,genotype) %>%
#   dplyr::filter(treatment == "untreated")  %>%
#   distinct() 

p1 = line_prop_changes_per_well %>%
  dplyr::filter(line =="hegp_3") %>% 
  dplyr::filter(!log_fraction_mean %in% c(NA,Inf,-Inf)) %>%
  ggplot(aes(y = log_fraction_mean, x = pool, 
             col = treatment, size = prop_unadjusted_max_value)) +
  geom_point() + 
  geom_hline(yintercept =0) +
  theme_bw() + 
  ggtitle("hegp_3 fractions in all pools and replicates.") + 
  ylab("log(mCherry+ prop/ mCherry- prop)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 = line_prop_changes_per_well %>%
  dplyr::filter(line =="hegp_3") %>% 
  dplyr::filter(!log_fraction_mean %in% c(NA,Inf,-Inf)) %>%
  ggplot(aes(y = log_fraction_mean, x = treatment, col = treatment)) +
  geom_boxplot() + 
  geom_hline(yintercept =0) +
  theme_bw() + 
  ylab("log(mCherry+ prop/ mCherry- prop)")


p3 = line_prop_changes_per_well %>%
  dplyr::filter(line =="aowh_2") %>% 
  dplyr::filter(!log_fraction_mean %in% c(NA,Inf,-Inf)) %>%
  ggplot(aes(y = log_fraction_mean, x = pool, 
             col = treatment, size = prop_unadjusted_max_value)) +
  geom_point() + 
  geom_hline(yintercept =0) +
  theme_bw() + 
  ggtitle("aowh_2 fractions in all pools and replicates.") + 
  ylab("log(mCherry+ prop/ mCherry- prop)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p4 = line_prop_changes_per_well %>%
  dplyr::filter(line =="aowh_2") %>% 
  dplyr::filter(!log_fraction_mean %in% c(NA,Inf,-Inf)) %>%
  ggplot(aes(y = log_fraction_mean, x = treatment, col = treatment)) +
  geom_boxplot() + 
  geom_hline(yintercept =0) +
  theme_bw() + 
  ylab("log(mCherry+ prop/ mCherry- prop)")

pdf(paste0(outdir2,"/fraction_changes_pool_treatment_hegp_3_aowh_2.pdf"),
    width = 12, height = 7)

(p1 + p2) / (p3 + p4)

dev.off()
# if splitting by genotype (hom REF, ALT, het) they don't change in the majority
# as expected from effects of individual donors being consistent

# for lines shared across pools, you'd expect them to have the same phagocytosis
# throughout ~ same log(mCherry+ prop/ mCherry- prop)
# this is not the case
# how does this deviation from average fraction correlate with maximum proportion?
# scale fraction and maximum proportion per donor and correlate
# should I do fraction calculations with scaled proportions per donor?
# scale with respect to the untreated prop to preserve bigger effect sizes?


# 1st scale and center proportions per pool
# checking how things would look like with hegp_3 and aowh_2
unsorted_hegp_3 =   w_estimates_adjusted %>%
  dplyr::filter(line =="hegp_3" & mCherry=="unsorted") %>%
  dplyr::mutate(value_to_center = prop_adjusted_mean) %>%
  dplyr::group_by(pool, treatment) %>%
  dplyr::select(pool,treatment,line,value_to_center) %>%
  distinct()


w_estimates_centered_scaled_hegp_3 = w_estimates_adjusted %>%
  dplyr::filter(line =="hegp_3" & mCherry!="unsorted" &  mCherry!="unsorted?") %>% 
  dplyr::group_by(pool, treatment) %>%
  dplyr::left_join(unsorted_hegp_3) %>%
  dplyr::mutate(value_to_center = case_when(
    is.na(value_to_center) ~ mean(prop_adjusted_mean),
    .default = value_to_center
  )) %>%
  dplyr::mutate(centered_prop = prop_adjusted_mean - value_to_center) %>% # trying to do it same way as scale() with minimum hassle
  dplyr::mutate(scaled_centered_prop = centered_prop / sd(centered_prop)) %>%
  ungroup()

p1 = w_estimates_centered_scaled_hegp_3 %>%
  dplyr::filter(line =="hegp_3") %>% 
  dplyr::mutate(grouping=paste0(pool,mCherry)) %>%
  ggplot(aes(y = scaled_centered_prop, x = grouping, 
             col = treatment, shape = mCherry)) +
  geom_point() + 
  geom_hline(yintercept =0) +
  theme_bw() + 
  ggtitle("hegp_3 scaled proportions in all pools and replicates.") + 
  ylab("Scaled and centered adjusted prop.") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 = w_estimates_centered_scaled_hegp_3 %>%
  dplyr::filter(line =="hegp_3") %>% 
  dplyr::mutate(grouping=paste0(treatment,mCherry)) %>%
  ggplot(aes(y = scaled_centered_prop, x =  grouping,
             col = treatment)) +
  geom_boxplot() + 
  geom_hline(yintercept =0) +
  theme_bw() + 
  ylab("Scaled and centered adjusted prop.") 

pdf(paste0(outdir2,"/scaled_centered_proportions_per_pool_hegp_3.pdf"),
    width = 15, height = 10)
p1 / p2
dev.off()

# not centering - negative values are not good for the fractions

w_estimates_scaled_hegp_3 = w_estimates_adjusted %>%
  dplyr::filter(line =="hegp_3" & mCherry!="unsorted" &  mCherry!="unsorted?") %>% 
  dplyr::group_by(pool, treatment) %>%
  dplyr::mutate(scaled_prop = prop_adjusted_mean / sqrt(sum(prop_adjusted_mean)^2/(n()-1))) #  as in scale() when no centering

p1 = w_estimates_scaled_hegp_3 %>%
  dplyr::filter(line =="hegp_3") %>% 
  dplyr::mutate(grouping=paste0(pool,mCherry)) %>%
  ggplot(aes(y = scaled_prop, x = grouping, 
             col = treatment, shape = mCherry)) +
  geom_point() + 
  geom_hline(yintercept =0) +
  theme_bw() + 
  ggtitle("hegp_3 scaled proportions in all pools and replicates.") + 
  ylab("Scaled adjusted prop.") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 = w_estimates_scaled_hegp_3 %>%
  dplyr::filter(line =="hegp_3") %>% 
  dplyr::mutate(grouping=paste0(treatment,mCherry)) %>%
  ggplot(aes(y = scaled_prop, x =  grouping,
             col = treatment)) +
  geom_boxplot() + 
  geom_hline(yintercept =0) +
  theme_bw() + 
  ylab("Scaled adjusted prop.") 

p3 = w_estimates_adjusted %>%
  dplyr::filter(line =="hegp_3" & mCherry!="unsorted" &  mCherry!="unsorted?") %>% 
  dplyr::mutate(grouping=paste0(pool,mCherry)) %>%
  ggplot(aes(y = prop_adjusted_mean, x = grouping, 
             col = treatment, shape = mCherry)) +
  geom_point() + 
  geom_hline(yintercept =0) +
  theme_bw() + 
  ggtitle("hegp_3 scaled proportions in all pools and replicates.") + 
  ylab("Adjusted prop. (no scaling)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p4 = w_estimates_adjusted %>%
  dplyr::filter(line =="hegp_3" & mCherry!="unsorted" &  mCherry!="unsorted?") %>% 
  dplyr::mutate(grouping=paste0(treatment,mCherry)) %>%
  ggplot(aes(y = prop_adjusted_mean, x =  grouping,
             col = treatment)) +
  geom_boxplot() + 
  geom_hline(yintercept =0) +
  theme_bw() + 
  ylab("Adjusted prop. (no scaling)") 

pdf(paste0(outdir2,"/scaled_proportions_per_pool_hegp_3.pdf"),
    width = 15, height = 10)
(p1 + p2) / (p3 + p4)
dev.off()

# hegp_3 fractions with 3 different methods

# usual log(mCherry+ prop/ mCherry- prop)
p1 = line_prop_changes_per_well %>%
  dplyr::filter(line =="hegp_3") %>% 
  dplyr::filter(!log_fraction_mean %in% c(NA,Inf,-Inf)) %>%
  ggplot(aes(y = log_fraction_mean, x = pool, 
             col = treatment, size = prop_unadjusted_max_value)) +
  geom_point() + 
  geom_hline(yintercept =0) +
  theme_bw() + 
  ggtitle("hegp_3 fractions in all pools and replicates.") + 
  ylab("log(mCherry+ prop/ mCherry- prop)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# log(mCherry+ scaled prop/ mCherry- scaled prop)
p2 = w_estimates_scaled_hegp_3 %>%
  dplyr::group_by(line,well) %>%
  dplyr::reframe(log_fraction_mean = log(scaled_prop[mCherry == "R3_mCherryPos"] / scaled_prop[mCherry == "R5_mCherryNeg"]),
                 prop_unadjusted_max_value = max(scaled_prop),
                 prop_unadjusted_min_value = min(scaled_prop),
                 pool=pool,condition=condition,treatment=treatment,replicate = replicate) %>%
  dplyr::ungroup() %>%
  dplyr::distinct() %>% 
  dplyr::filter(line =="hegp_3") %>% 
  dplyr::filter(!log_fraction_mean %in% c(NA,Inf,-Inf)) %>%
  dplyr::mutate(pool = factor(pool,levels = paste0("pool",c(2:11,13)))) %>%
  ggplot(aes(y = log_fraction_mean, x = pool, 
             col = treatment, size = prop_unadjusted_max_value)) +
  geom_point() + 
  geom_hline(yintercept =0) +
  theme_bw() + 
  ggtitle("hegp_3 fractions in all pools and replicates.") + 
  ylab("log(mCherry+ scaled prop/ mCherry- scaled prop)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# centering and scaling
p3 =  w_estimates_centered_scaled_hegp_3 %>%
  dplyr::group_by(line,well) %>%
  dplyr::reframe(log_fraction_mean = log(scaled_centered_prop[mCherry == "R3_mCherryPos"] / scaled_centered_prop[mCherry == "R5_mCherryNeg"]),
                 prop_unadjusted_max_value = max(scaled_centered_prop),
                 prop_unadjusted_min_value = min(scaled_centered_prop),
                 pool=pool,condition=condition,treatment=treatment,replicate = replicate) %>%
  dplyr::ungroup() %>%
  dplyr::distinct() %>% 
  dplyr::filter(line =="hegp_3") %>% 
  dplyr::filter(!log_fraction_mean %in% c(NA,Inf,-Inf)) %>%
  dplyr::mutate(pool = factor(pool,levels = paste0("pool",c(2:11,13)))) %>%
  ggplot(aes(y = log_fraction_mean, x = pool, 
             col = treatment, size = prop_unadjusted_max_value)) +
  geom_point() + 
  geom_hline(yintercept =0) +
  theme_bw() + 
  ggtitle("hegp_3 fractions in all pools and replicates.") + 
  ylab("log(mCherry+ scaled-cent prop/ mCherry- scaled-cent prop)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(paste0(outdir2,"/scaled_proportions_per_pool_then_calc_fractions_across_pools_hegp_3.pdf"),
    width = 8, height = 5)
p2/p3
dev.off()

# scaled and centered (p3) doesn't work because of the negative values produced in centering
# scaled and unscaled produce same results

# zscores (scale and center) after calculating the fractions? #######

fractions_centered_scaled_hegp_3 = line_prop_changes_per_well %>%
  dplyr::filter(line =="hegp_3") %>%
  dplyr::filter(!log_fraction_mean %in% c(NA,Inf,-Inf)) %>%
  dplyr::group_by(treatment) %>%
  dplyr::mutate(scaled_centered_fraction = scale(log_fraction_mean)) %>%
  ungroup()

p1 = line_prop_changes_per_well %>%
  dplyr::filter(line =="hegp_3") %>% 
  dplyr::filter(!log_fraction_mean %in% c(NA,Inf,-Inf)) %>%
  ggplot(aes(y = log_fraction_mean, x = pool, 
             col = treatment, size = prop_unadjusted_max_value)) +
  geom_point() + 
  geom_hline(yintercept =0) +
  theme_bw() + 
  ggtitle("hegp_3 fractions in all pools and replicates.") + 
  ylab("log(mCherry+ prop/ mCherry- prop)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# scaled log(mCherry+ prop/ mCherry- prop)
p2 = fractions_centered_scaled_hegp_3 %>%
  dplyr::mutate(pool = factor(pool,levels = paste0("pool",c(2:11,13)))) %>%
  ggplot(aes(y = scaled_centered_fraction, x = pool, 
             col = treatment, size = prop_unadjusted_max_value)) +
  geom_point() + 
  geom_hline(yintercept =0) +
  theme_bw() + 
  ggtitle("hegp_3 fractions - rescaled across pools and replicates") + 
  ylab("scaled log(mCherry+ prop/ mCherry- prop)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(paste0(outdir2,"/scaled_fractions_across_pools_hegp_3.pdf"),
    width = 8, height = 5)
p1/p2
dev.off()

# rescaling all lines and plot heatmap ####
# to see if any particular pool is an outlier
# only for lines shared across more than 3 pools - otherwise it's meaningless
# (those present only in one pool will have a value of 0)

# lines shared across 3 or more pools:
shared_lines = line_prop_changes_per_well %>%
  dplyr::filter(treatment == "untreated") 
sort(rowSums(table(shared_lines$line,shared_lines$pool)>0),decreasing = TRUE)
shared_lines = c("aowh_2",  "hegp_3",  "letw_5",  "lizq_1",
                 "laey_6",  "zaie_1",  "zaie_5")
scaled_fractions_all_lines = line_prop_changes_per_well %>%
  dplyr::filter(!log_fraction_mean %in% c(NA,Inf,-Inf)) %>%
  dplyr::group_by(treatment, line) %>%
  dplyr::mutate(scaled_centered_fraction = scale(log_fraction_mean)) %>%
  dplyr::group_by(treatment, line, pool) %>%
  dplyr::summarise(scaled_centered_fraction_mean_per_pool = abs(mean(scaled_centered_fraction))) %>%
  # abs to see better the extreme values in the heatmap
  ungroup() 

for_pheatmap=scaled_fractions_all_lines %>%
  dplyr::filter(treatment == "untreated" & line %in% shared_lines) %>%
  tidyr::pivot_wider(id_cols = c(line), 
                     names_from = pool, values_from = scaled_centered_fraction_mean_per_pool) 
lines = for_pheatmap$line
pools = colnames(for_pheatmap)[-1]
for_pheatmap = for_pheatmap %>%
  dplyr::select(-line)
mat = as.matrix(for_pheatmap)
colnames(mat) = pools
rownames(mat) = lines
# can't cluster due to amount of NAs
p = pheatmap(mat, kmeans_k = NA,cluster_rows = FALSE, cluster_cols = FALSE)

pdf(paste0(outdir2,"/abs_scaled_fractions_across_pools_shared_lines_untreated.pdf"),
    width = 8, height = 5)
print(p)
dev.off()

for_pheatmap=scaled_fractions_all_lines %>%
  dplyr::filter(treatment == "IFN" & line %in% shared_lines) %>%
  tidyr::pivot_wider(id_cols = c(line), 
                     names_from = pool, values_from = scaled_centered_fraction_mean_per_pool) 
lines = for_pheatmap$line
pools = colnames(for_pheatmap)[-1]
for_pheatmap = for_pheatmap %>%
  dplyr::select(-line)
mat = as.matrix(for_pheatmap)
colnames(mat) = pools
rownames(mat) = lines
# can't cluster due to amount of NAs
p = pheatmap(mat, kmeans_k = NA,cluster_rows = FALSE, cluster_cols = FALSE)

pdf(paste0(outdir2,"/abs_scaled_fractions_across_pools_shared_lines_IFN.pdf"),
    width = 8, height = 5)
print(p)
dev.off()

for_pheatmap=scaled_fractions_all_lines %>%
  dplyr::filter(treatment == "LPS" & line %in% shared_lines) %>%
  tidyr::pivot_wider(id_cols = c(line), 
                     names_from = pool, values_from = scaled_centered_fraction_mean_per_pool) 
lines = for_pheatmap$line
pools = colnames(for_pheatmap)[-1]
for_pheatmap = for_pheatmap %>%
  dplyr::select(-line)
mat = as.matrix(for_pheatmap)
colnames(mat) = pools
rownames(mat) = lines
# can't cluster due to amount of NAs
p = pheatmap(mat, kmeans_k = NA,cluster_rows = FALSE, cluster_cols = FALSE)

pdf(paste0(outdir2,"/abs_scaled_fractions_across_pools_shared_lines_LPS.pdf"),
    width = 8, height = 5)
print(p)
dev.off()
#############
# case_when to deal with negative values in centering? Re-scaling to positives

# check out
scales::rescale(x = c(0.05,0.1,0.5),to = c(1,10),from = range(c(0,1), na.rm = TRUE, finite = TRUE))
# https://www.theanalysisfactor.com/rescaling-variables-to-be-same/

########## thoughts
# ensure consistent directionality of proportions per line across replicates and across pools
# if not consistent across pools, remove pool with lowerproportion and test again
# if not consistent across replicates, remove altogether

# how does the bigger deviations from averge in fractions in hegp_3 can be used 
# to correct the rest?
# maybe provide size of correction per pool
# but would need to leverage with other donors that are shared across fewer pools
# and take into account directionality
# similar to weights used in RNA-seq to normalise against total library size per sample?
# https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29
# https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html#:~:text=However%2C%20DESeq2%20uses%20a%20specific,variation%20(%20Var%20%2F%20μ%20).
# similarities: variation in fraction is larger the smaller the average prop
# differences: this will be different per pool
# would need to measure this for shared donors
# need to play around with data to see this in a clearer way

### change in analysis
# mixed models? pool and line as random effects