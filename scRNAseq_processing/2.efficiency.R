# measuring efficiency of differentiation
# pools 2, 3, 4, 5, 6
.libPaths("/software/teamtrynka/ma23/R4.1/libs") # path to libraries
library(patchwork)
library(tidyverse)
library(future)
source("./helpers.R")
directory = "../../data/results/2.efficiency/"
dir.create(directory, recursive = T)

options(future.globals.maxSize= 150097152000)


##### start ########


sc_ncells = load_sc_donors(doublets_unassigned_in_proportion = FALSE,
                           # if TRUE, consider doublet and unassigned in final proportion 
                           # (final proportions for individual lines will be underestimated, particularly for donors with high proportions?)
                           # fix so doublets are assigned to their individual donors
                           pools = c("pool2", "pool3", "pool4", "pool5", "pool6"),
                           directory =  "../../../OTAR2065_data_processing/data/")

# find shared donors - top 2 most frequent lines
shared_donors = sc_ncells %>% 
  group_by(Line) %>% 
  summarise(common=n()) %>% 
  filter(row_number(desc(common)) %in% 1:2) %>% 
  .$Line 


# center the data around the starting proportion
sc_ncells = sc_ncells %>%
  mutate(centered_prop = final_prop - starting_prop) %>%
  group_by(sample) %>% # should I scale to sd within each sample or of all samples together?
  mutate(scaled_centered_prop = (final_prop - starting_prop) / sd(final_prop), 
         log1p_scaled_prop = log1p(final_prop) / log1p(starting_prop)) %>% 
  # Puigdevall et al prolif. rate: not much different distribution from raw props.
  ungroup()

plot_prop_hist(df = sc_ncells, 
               directory = "../../data/results/2.efficiency/",
               plot.name = "scaled_prop_histograms.png",
               nbins = 50,
               plot.width = 11,plot.height = 5)

# plot scaled_centered coloring those that are 0 final prop


plot_prop_scatter(df = sc_ncells,
                  directory = "../../data/results/2.efficiency/",
                  plot.name = "scaled_prop_scatter.png",
                  plot.width = 11,plot.height = 6)


sc_ncells = sc_ncells %>%
  mutate(centered_prop = final_prop - starting_prop) %>%
  group_by(sample) %>% # should I scale to sd within each sample or of all samples together?
  mutate(scaled_centered_prop = (final_prop - starting_prop) / sd(final_prop), 
         log1p_scaled_prop = log1p(final_prop) / log1p(starting_prop)) %>% 
  # Puigdevall et al prolif. rate: not much different distribution from raw props.
  ungroup()

# I believe the scaled and centered distribution is the most appropriate
# As it's more normal and the <0 has an intuitive meaning (lower final prop. than starting prop.)

## load Rouhani 2022 data (BCOR and other mutations affecting neuron differentiation)

Rouhani = list()
files = list.files("../../../resources/Rouhani_BCOR_selection_2022/", pattern = ".csv")
Rouhani = lapply(paste0("../../../resources/Rouhani_BCOR_selection_2022/",files), function(x){
  read.csv(x, header=TRUE)
})
names(Rouhani) = gsub(pattern = ".csv",replacement = "",x = files)

# gather all BCOR-related mutations
BCOR = list()
BCOR[["S5_cancermut_DSNV"]] = Rouhani$S5_cancermut_DSNV %>%
  filter(gene == "BCOR") %>%
  select(ips,gene,type,pvalue,FATHMM) %>%
  mutate(table = "Rouhani_S5_cancermut_DSNV")

BCOR[["Rouhani_S5_cancermut_SNV"]] = Rouhani$S5_cancermut_DSNV %>%
  filter(gene == "BCOR") %>%
  select(ips,gene,type,pvalue,FATHMM) %>%
  mutate(table = "Rouhani_S5_cancermut_SNV")

BCOR[["S5_cancermut_indel"]] = Rouhani$S5_cancermut_indel %>%
  filter(Gene == "BCOR") %>%
  select(ips,Gene,pvalue) %>%
  mutate(table = "S5_cancermut_indel")

