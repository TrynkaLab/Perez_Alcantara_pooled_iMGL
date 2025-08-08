# create table of donors (actually, lines) in pools
.libPaths(c('/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/',"/software/teamtrynka/ma23/R4.1/libs",
            "/software/R-4.1.0/lib/R/library"))
library(patchwork)
library(tidyverse)
library(Seurat)
library(future)
library(janitor)
source("./helpers.R")

directory = "../../data/results/0.donor_numbers_table/"
dir.create(directory, recursive = T)

# read in line proportions from scRNA-seq

sc_ncells = load_sc_donors(doublets_unassigned_in_proportion = TRUE,
                           # if TRUE, consider doublet and unassigned in final proportion 
                           # (final proportions for individual lines will be underestimated, particularly for donors with high proportions?)
                           # fix so doublets are assigned to their individual donors
                           pools = c(paste0("pool",2:15)),
                           directory =  "../../../OTAR2065_scRNA_data_processing/data/")
sc_ncells = sc_ncells %>%
  mutate(pool = factor(pool, levels = c(paste0("pool",2:15)), ordered = TRUE),
         sample = factor(sample, levels = unique(sample), ordered = TRUE)
)

# pivot to wider format and save

sc_ncells %>% 
  pivot_wider(id_cols = Line,names_from = set,values_fill = 0,values_from = Frequency) %>%
  janitor::adorn_totals(where = c("row", "col")) %>%
  write.csv(.,paste0(directory,"donor_frequency_all_pools_with_doublets_and_unassigned.csv"),quote = FALSE, row.names = FALSE)

# per 10k cells

sc_ncells %>% 
  group_by(sample) %>%
  dplyr::filter(Line!="unassigned_plus_doublets") %>%
  mutate(freq_per_10k=(Frequency/sum(Frequency))*10000) %>%
  pivot_wider(id_cols = Line,names_from = set,values_fill = 0,values_from = freq_per_10k) %>%
  janitor::adorn_totals(where = c("row", "col")) %>%
  write.csv(.,paste0(directory,"donor_frequency_per_10k_all_pools_with_doublets_and_unassigned.csv"),quote = FALSE, row.names = FALSE)

# check only proportions of unassigned + doublets in barplot

png(paste0(directory,"unassigned_plus_doublets_vireo.png"),res = 400,
    units = "in",  width = 10, height = 6)
sc_ncells %>%
  dplyr::filter(Line=="unassigned_plus_doublets") %>%
  ggplot( aes(x=sample, y=final_prop, fill = pool)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  ggtitle("Fraction unassigned + doublets") + 
  ylab("Fraction unassigned + doublets") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
dev.off()

#singlets per pool, boxplot

png(paste0(directory,"singlet_boxplot_vireo.png"),res = 400,
    units = "in",  width = 6, height = 6)
sc_ncells %>%
  dplyr::filter(Line!="unassigned_plus_doublets") %>%
  mutate(singlet_status = case_when(!Line %in% c("unassigned","doublet") ~ "singlet")) %>%
  dplyr::filter(singlet_status == "singlet") %>%
  group_by(sample) %>%
  summarise(percent = sum(final_prop),pool=pool) %>%
  unique() %>%
  ggplot( aes(x=pool, y=percent)) +
  geom_boxplot()+
  geom_jitter(alpha=0.3)+
  ggtitle("Percentage of singlets per pool") + 
  ylab("Percentage of singlets") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),limits = c(0,1)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
dev.off()

png(paste0(directory,"singlets_unassigned_doublets_vireo.png"),res = 400,
    units = "in",  width = 10, height = 6)
sc_ncells %>%
  dplyr::filter(Line!="unassigned_plus_doublets") %>%
  mutate(singlet_status = case_when(!Line %in% c("unassigned","doublet") ~ "singlet",
                             Line == "unassigned" ~ "unassigned",
                             Line == "doublet" ~ "doublet")) %>%
  mutate(singlet_status = factor(singlet_status,levels = c("doublet","unassigned","singlet"), 
                                 ordered = TRUE)) %>%
  ggplot( aes(x=sample, y=Frequency, fill = singlet_status)) +
  geom_bar(stat="identity", color="black") +
  ggtitle("Frequency of donor identities") + 
  ylab("Number of cells") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1))
dev.off()

png(paste0(directory,"singlets_unassigned_doublets_per_treatment_vireo.png"),res = 400,
    units = "in",  width = 10, height = 6)
sc_ncells %>%
  dplyr::filter(Line!="unassigned_plus_doublets") %>%
  mutate(singlet_status = case_when(!Line %in% c("unassigned","doublet") ~ "singlet",
                                    Line == "unassigned" ~ "unassigned",
                                    Line == "doublet" ~ "doublet"),
         treatment = case_when(grepl("plst|untreated|Untreated|AG|YC|ITMG", sample) ~ "untreated",
                               grepl("IFN", sample) ~ "IFN",
                               grepl("LPS", sample)~ "LPS")) %>%
  mutate(singlet_status = factor(singlet_status,levels =c("doublet","unassigned","singlet"), ordered = TRUE)) %>%
  group_by(treatment, singlet_status) %>%
  summarise(n = sum(Frequency)) %>%
  mutate(freq = paste0(round( n / sum(n)*100,0),"%")) %>%
ggplot( aes(x=treatment, y=n, fill = singlet_status, color = singlet_status, label = freq)) +
  geom_bar(stat="identity")+
  ggtitle("Frequency of donor identities") + 
  geom_text( col = "white",position = position_stack(vjust = 0.5), size = 6) + 
  ylab("Number of cells") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
dev.off()


# now without doublets and unassigned
# review donors present in pool 11 for which we don't have genotype info right now
sc_ncells = load_sc_donors(doublets_unassigned_in_proportion = FALSE,
                           # if TRUE, consider doublet and unassigned in final proportion 
                           # (final proportions for individual lines will be underestimated, particularly for donors with high proportions?)
                           # fix so doublets are assigned to their individual donors
                           pools = c(paste0("pool",2:15)),
                           directory =  "../../../OTAR2065_scRNA_data_processing/data/")

sc_ncells = sc_ncells %>%
  mutate(pool = factor(pool, levels = c(paste0("pool",2:15)), ordered = TRUE),
         sample = factor(sample, levels = unique(sample), ordered = TRUE)
  )

# fix hegp_3 proportion in pool 8 (2.5%)
sc_ncells$starting_prop = ifelse(sc_ncells$pool == "pool8" & sc_ncells$Line == "hegp_3",
                                 0.025, sc_ncells$starting_prop)
# same in pool 11 (x0.25, ~1.075%)
sc_ncells$starting_prop = ifelse(sc_ncells$pool == "pool11" & sc_ncells$Line == "hegp_3",
                                 0.01075, sc_ncells$starting_prop)

# same in pool 15 (x0.25, ~1.013%), double check with Juliette
sc_ncells$starting_prop = ifelse(sc_ncells$pool == "pool15" & sc_ncells$Line == "hegp_3",
                                 0.01075, sc_ncells$starting_prop)

# same in pool 15 (x0.25, ~1.013%), double check with Juliette
sc_ncells$starting_prop = ifelse(sc_ncells$pool == "pool15" & sc_ncells$Line == "hegp_3",
                                 0.01075, sc_ncells$starting_prop)

### double-check and fix for other pools
# fix kolf_2 proportion in pool 9 
sc_ncells$starting_prop = ifelse(sc_ncells$pool == "pool9" & sc_ncells$Line == "kolf_2",
                                (1/72) *2, sc_ncells$starting_prop)

## I should maybe adjust the proportions of the other donors but they are so small it probably doesn't matter?

sc_ncells %>% 
  pivot_wider(id_cols = Line,names_from = set,values_fill = 0,values_from = Frequency) %>%
  janitor::adorn_totals(where = c("row", "col")) %>%
  write.csv(.,paste0(directory,"donor_frequency_all_pools.csv"),quote = FALSE, row.names = FALSE)

# per 10k
sc_ncells = sc_ncells %>% 
  group_by(sample) %>%
  mutate(freq_per_10k=(Frequency/sum(Frequency))*10000) 

# find shared donors - top 2 most frequent lines
shared_donors = sc_ncells %>% 
  group_by(Line) %>% 
  summarise(common=n()) %>% 
  filter(row_number(desc(common)) %in% 1:2) %>% 
  .$Line 

# barplots - not with facet_wrap to repeat colors

p_list = list()
for(p in unique(sc_ncells$pool)){
  filtered =  sc_ncells %>% 
    filter(pool == p) %>%
    mutate(Line = factor(Line, 
                         levels = unique(Line[order(Frequency, decreasing = TRUE)], ordered = TRUE))) 
  p_list[["Frequency_per_10k"]] = filtered %>%
    ggplot(aes(x = Line, y = freq_per_10k, fill = sample)) + 
    geom_bar(stat = "identity",position = "dodge") + 
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),breaks = c(0,10,50,100,1000,2000,5000,10000),
                       limits = c(0,10000)) + 
    theme_bw() +
    ylab("Frequency per 10k cells (pseudolog10 axis)") +
    theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
    coord_flip()
  
  p_list[["Proportion"]] = filtered %>%
    ggplot(aes(x = Line, y = final_prop, fill = sample)) + 
    geom_bar(stat = "identity",position = "dodge") + 
    theme_bw() +
    geom_hline(aes(yintercept = min(unique(starting_prop))),linetype = "dashed") + 
    ylab("Proportion of cells") + 
    theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
    coord_flip()
  plots = patchwork::wrap_plots(p_list) + 
    plot_annotation(title = paste0(p,": ",unique(filtered$n_lines)," lines")) +
    plot_layout(guides = "collect")
  
  pdf(paste0(directory,p,"_barplot_line_freq.pdf"),
      width = 10,height = 16)
  plot(plots)
  dev.off()
}

sc_ncells = sc_ncells %>%
  group_by(set) %>%
  mutate(centered_control = final_prop - mean(final_prop[Line %in% shared_donors])) %>%
  mutate(scaled_centered_control = centered_control/ sd(final_prop),
         centered_starting = final_prop - starting_prop,
         scaled_centered_starting = (final_prop - starting_prop) / sd(final_prop),
         log1p_scaled_starting = log1p(final_prop) / log1p(starting_prop),   # Puigdevall et al prolif. rate: not much different distribution from raw props.
         log1p_substract_starting_scaled = (log1p(final_prop) - log1p(starting_prop)) /  log1p(mean(final_prop[Line %in% shared_donors])) # modified puigdevall
         ) %>%
  ungroup()

sc_ncells %>%
  select(set,Line,Frequency,pool,sample,n_lines,starting_prop,final_prop,freq_per_10k,log1p_scaled_starting) %>%
  write.csv(.,paste0(directory,"donor_proliferation_for_pau.csv"),quote = FALSE, row.names = FALSE)

# boxplot all donors
png(paste0(directory,"boxplot_line_log1p_scaled_starting_prop.png"),
    width = 30,height = 6,units = "in",res = 400)
sc_ncells %>%
  ggplot(aes(x=reorder(Line,-log1p_scaled_starting,FUN = median),y=log1p_scaled_starting)) +
  geom_boxplot() +
  theme_bw() + 
  geom_hline(yintercept = 1,linetype = "dotted") +
  ggtitle("log1p(final_prop) / log1p(starting_prop)") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  +
  coord_flip()
  
dev.off()

pdf(paste0(directory,"boxplot_line_scaled_centered_starting_prop.pdf"),
    width = 6,height = 30)
sc_ncells %>%
  ggplot(aes(x=reorder(Line,-scaled_centered_starting,FUN = median),y=scaled_centered_starting)) +
  geom_boxplot() +
  theme_bw() + 
  geom_hline(yintercept = 0,linetype = "dotted") +
  ggtitle("(final_prop - starting_prop) / sd(final_prop)") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip()

dev.off()

# averaging for pools 13, also boxplots
pdf(paste0(directory,"boxplot_line_log1p_scaled_starting_prop_pool15.pdf"),
    width = 6,height = 30)
sc_ncells %>%
  dplyr::filter(pool %in% c("pool15"))%>%
  ggplot(aes(x=reorder(Line,-log1p_scaled_starting,FUN = median),y=log1p_scaled_starting)) +
  geom_boxplot() +
  theme_bw() + 
  geom_hline(yintercept = 1,linetype = "dotted") +
  ggtitle("Scaled proportions pool 15 ") + 
  xlab("Donors from pool 15")+
  ylab("log1p(final_prop) / log1p(starting_prop)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  coord_flip()

dev.off()

pdf(paste0(directory,"boxplot_summed_cell_abundance_per_treatment_pool15.pdf"),
    width = 6,height = 30)
sc_ncells %>%
  dplyr::filter(pool %in% c("pool15"))%>%
  mutate(treatment = case_when(grepl("plst|untreated|Untreated|AG|YC|ITMG", sample) ~ "untreated",
                        grepl("IFN", sample) ~ "IFN",
                        grepl("LPS", sample)~ "LPS")) %>%
  group_by(Line,pool,treatment) %>%
  summarise(Frequency = sum(Frequency)) %>%
  ungroup() %>%
  ggplot(aes(x=reorder(Line,-Frequency,FUN = median),
             y=Frequency+1)) +
  geom_boxplot() +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000)) +
  theme_bw() + 
  geom_hline(yintercept = 100,linetype = "dotted") +
  ggtitle("Number of cells (added per pool and treatment),pool 15",
          "Dotted line indicates 100 cells") + 
  xlab("Donors from pool 15")+
  ylab("total number of cells + 1 (log10 transformed)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  coord_flip()

dev.off()

# save
sc_ncells %>%
  dplyr::filter(pool %in% c("pool15"))%>%
  group_by(Line) %>%
  summarise(sum_ncells = sum(Frequency),
            average_diff_efficiency = mean(log1p_scaled_starting)) %>%
  arrange(average_diff_efficiency) %>%
  write_csv(file=paste0(directory,"averaged_diff_efficiency_sum_ncells_per_line_pool15.csv"))

################################################################################
# calculate correlation of scores for replicate donors
################################################################################
# within pools, within conditions, between replicates

for(pooll in c("pool2","pool9","pool10")){
  for(cond in c("untreated","IFN","LPS")){
  p = within_pool_within_cond_between_rep(sc_ncells,pool = pooll,condition = cond)
  
  png(paste0(directory,"within_pool_within_cond_between_rep_scaled_centered_starting_prop_",pooll,"_",cond,".png"),
      width = ceiling((length(p$grobs)+1)*2),
      height = ceiling((length(p$grobs)+1)*1.3),
      units = "in",res = 400)
  plot(p)
  dev.off()
  }
}

# between pools, within conditions

  for(cond in c("untreated","IFN","LPS")){
    message(cond)
    p =   between_pools_within_conditions(sc_ncells,condition = cond)
    
    png(paste0(directory,"between_pool_within_cond_scaled_centered_starting_prop_",cond,".png"),
        width = 6,
        height = 6,
        units = "in",res = 400)
    plot(p)
    dev.off()
    cor = between_pools_within_conditions_correlations(sc_ncells,condition = cond) 
    print(min(cor))
    print(max(cor))
    print(mean(cor))
    # results from paired correlation test in 1000 shuffles of pools assignments to pool A or B:
    # untreated: 0.62-0.68 (mean = 63)
    # IFN: 0.60-0.67 (mean = 0.60)
    # LPS: 0.65-0.73 (mean=0.66)
  }
  
###########

# zscore of these, per donor: (x - mu)/sd
# calculating mean per donor where there are repeated values
zscores = sc_ncells %>%
  group_by(Line) %>%
  summarise(log1p_scaled_starting_per_donor=mean(log1p_scaled_starting),
            scaled_centered_starting_per_donor = mean(scaled_centered_starting)) %>%
  summarise(Line= Line,
            z_log1p = (log1p_scaled_starting_per_donor - mean(log1p_scaled_starting_per_donor))/sd(log1p_scaled_starting_per_donor),
            z_centered = (scaled_centered_starting_per_donor - mean(scaled_centered_starting_per_donor))/sd(scaled_centered_starting_per_donor))

png(paste0(directory,"density_log1p_scaled_starting_prop_zscore.png"),
    width = 6,height = 6,units = "in",res = 400)
zscores %>%
  ggplot(aes(x=z_log1p)) +
  geom_density()
dev.off()

png(paste0(directory,"density_scaled_centered_starting_prop_zscore.png"),
    width = 6,height = 6,units = "in",res = 400)
zscores %>%
  ggplot(aes(x=z_centered)) +
  geom_density()
dev.off()
################################################################################
######### z scores are non-normal, so can't calculate outliers from these ########
######### get just bottom 20% and top 20% of data? ########
################################################################################

ourliers_classif = sc_ncells %>%
  group_by(Line) %>%
  summarise(scaled_centered_starting_median = median(scaled_centered_starting)) %>%
  arrange(desc(scaled_centered_starting_median)) %>%
  mutate(quantile_rank = percent_rank(scaled_centered_starting_median)) %>%
  mutate(outlier_classif = case_when(quantile_rank > 0.80 ~ "top_20",
                                     quantile_rank < 0.20 ~ "bottom_20",
                                     quantile_rank <= 0.80 & quantile_rank >= 0.20 ~ "middle_80"))

# read in file with information from Puigdevall's paper 2023
puig = read.csv("../../../resources/Puigdevall_Neuroseq_efficiency_2023/TableS1.csv", 
                header = TRUE)
puig = puig %>%
  mutate(donor_iPSC = str_replace_all(string = donor_iPSC,pattern = ".*-",replacement = "")) %>%
  dplyr::rename(Line = donor_iPSC)

ourliers_classif = merge(ourliers_classif, puig, by = "Line")

for(mut in c("all", "synonymous" ,      "missense", "deleterious",
             "missPatho" ,"ptv", "other" ,"dinuc" , "onlySnvs" ,"indels" ,
             "cnv_num_different_regions" , "cnv_length_different_regions_Mbp" ,
             "cnv_length_shared_differences_Mbp")){
  p = ggplot(ourliers_classif, aes_string(x="outlier_classif",
                                          y=mut)) +
    geom_violin()+
    geom_jitter(alpha = 0.5) + theme_bw()
  
  png(paste0(directory,"prolif_efficiency_mutation_burden_",mut,".png"),
      width = 6,
      height = 6,
      units = "in",res = 400)
  plot(p)
  dev.off()
}


plot_prop_hist(df = sc_ncells, 
               directory = directory,
               plot.name = "scaled_prop_histograms.png",
               nbins = 50,
               plot.width = 11,plot.height = 9)

plot_prop_scatter(df = sc_ncells,
                  directory = directory,
                  plot.name = "scaled_prop_scatter.png",
                  plot.width = 6,plot.height = 10)

plot_prop_density(df = sc_ncells,
                  directory = directory,
                  plot.name = "scaled_prop_density_by_set.png",
                  plot.width = 13,plot.height = 8)

plot_prop_density(df = sc_ncells,
                  directory = directory,
                  plot.name = "scaled_prop_density_by_pool.png",
                  col.by = "pool",plot.legend = TRUE,
                  plot.width = 13,plot.height = 8)

#### revise everything after this point ####



# separate histograms by pool
#### do the normalisation to shared donor for differentiation efficiency!!
## I believe it doesn't make sense for proliferation (we already normalise to total counts)
# I believe the best thing to do for selecting best and worst performers across all differentiations is 
# to center to the means of the shared donors, and scale so that they all have same stdev
# distribution of lines doesn't change much from final_prop, but lower props have moved a bit to the center
# centering around the mean of these two donors means that for example, if the mean 
# of the two is lower in IFN sample vs LPS / untreated in pool7
# we assume a less efficient differentiation in general in that sample compared to others
rm(plots,filtered, p_list)

p_list = list()
for(p in unique(sc_ncells$pool)){
  filtered =  sc_ncells %>% 
    group_by(set) %>%
    mutate(centered = final_prop - mean(final_prop[Line %in% shared_donors])) %>%
    mutate(scaled = centered/ sd(centered)) %>%
    ungroup() %>%
    filter(pool == p) %>%
    mutate(Line = factor(Line, 
                         levels = unique(Line[order(Frequency, decreasing = TRUE)], ordered = TRUE))) 
    
  p_list[["Proportion"]] = filtered %>%
    ggplot(aes(x = Line, y = final_prop, fill = sample)) + 
    geom_bar(stat = "identity",position = "dodge") + 
    theme_bw() +
    geom_hline(aes(yintercept = unique(starting_prop)),linetype = "dashed") + 
    ylab("Proportion of cells") + 
    theme(axis.text.x = element_text(angle = 45,hjust = 1)) 
  
  p_list[["Scaled"]] = filtered %>%
    ggplot(aes(x = Line, y = scaled, fill = sample)) + 
    geom_bar(stat = "identity",position = "dodge") + 
    theme_bw() +
    ylab("Scaled proportion of cells") + 
    theme(axis.text.x = element_text(angle = 45,hjust = 1)) 
  
  plots = patchwork::wrap_plots(p_list) + 
    plot_annotation(title = paste0(p,": ",unique(filtered$n_lines)," lines")) +
    plot_layout(guides = "collect")
  
  png(paste0(directory,p,"_barplot_line_scaled_prop.png"),
      width = 16,height = 6,units = "in",res = 400)
  plot(plots)
  dev.off()
}

# donor frequency AFTER FILTERING ########
########### AFTER FILTERING ########
########### AFTER FILTERING ########
########### AFTER FILTERING ########
# do per condition, extract metadata and merge
seurat_list = readRDS("../../data/results/1.QC/filtered_media_noRegression_list.rds")

  input_merged =  merge(seurat_list[[1]],
                        seurat_list[2:length(seurat_list)],
                        add.cell.ids = names(seurat_list))
  rm(seurat_list)
gc()
# get percentages per pool
sc_ncells = input_merged %>%
  group_by(orig.ident, donor_id) %>%
  summarise(Frequency = n(), sample = sample, pool=pool, treatment = treatment) %>%
  unique() %>%
  rename(Line=donor_id) %>%
  ungroup() %>%
  split(f = as.factor(.$orig.ident)) 

# fill in missing donors per pool
lines_in_pools = read.table("../../../OTAR2065_scRNA_data_processing/data/lines_in_pools.txt", 
                            header = TRUE)

for (pool in unique(input_merged$pool)) {
  message("Working on pool ", pool, "...")
  n_lines = lines_in_pools[lines_in_pools$Pool == pool,"N_lines"]
  
  lines = unlist(strsplit(lines_in_pools[lines_in_pools$Pool == pool, "Lines"], split = ";"))
  for (ident in names(sc_ncells)[grep(pool, names(sc_ncells))]) {
    sc_ncells[[ident]] =  sc_ncells[[ident]] %>%
      rows_insert(tibble(Line = lines[!lines %in% sc_ncells[[ident]]$Line])) %>%
      replace_na(list(Frequency = 0)) %>%
      mutate(orig.ident = unique(na.omit(orig.ident)), # fill in details for donors that are ==0
             sample = unique(na.omit(sample)),
             pool = unique(na.omit(pool)),
             treatment = unique(na.omit(treatment)))
    
    sc_ncells[[ident]]  = sc_ncells[[ident]] %>%
      mutate(n_lines = n_lines) %>%
      mutate(starting_prop = 1 / n_lines) %>%
      mutate(final_prop = Frequency / sum(Frequency)) %>%
      mutate(frequency_per_10k = Frequency / 10000,
             log1p_scaled_starting = log1p(final_prop) / log1p(starting_prop)   # Puigdevall et al prolif. rate: not much different distribution from raw props.
      )
    
  }
}
sc_ncells = do.call(rbind,sc_ncells)

# fix hegp_3 proportion in pool 8 (2.5%)
sc_ncells$starting_prop = ifelse(sc_ncells$pool == "pool8" & sc_ncells$Line == "hegp_3",
                                 0.025, sc_ncells$starting_prop)

# pivot to wider format and save

sc_ncells %>% 
  pivot_wider(id_cols = Line,names_from = orig.ident,values_fill = 0,values_from = Frequency) %>%
  janitor::adorn_totals(where = c("row", "col")) %>%
  write.csv(.,paste0(directory,"donor_frequency_all_pools_after_QC_filtering.csv"),quote = FALSE, row.names = FALSE)
