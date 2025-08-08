# check otar overlaps from colocs
# which ones have genetic scores from credible sets over 0.2

library(tidyverse)

colocs = readr::read_csv("../../data/results/8.colocalisation_analysis/coloc_results/all_GWAS_colocalisations.csv") %>%
  dplyr::filter(PP_H4 > 0.70)

scores = list()
for(disease in c("AD","PD","ALS","MS")){
  file = list.files("../../../resources/open_targets/general_disease_genetic_association_scores/",
                    pattern = disease,full.names = TRUE)
  scores[[disease]] = readr::read_tsv(file) %>%
    dplyr::mutate(disease = disease) %>%
    dplyr::select(symbol,globalScore,gwasCredibleSets,disease)

  
}

scores = do.call("rbind",scores)

colocs = colocs %>%
  # take distinct
  dplyr::select(eGene,GWAS) %>%
  distinct() %>%
  dplyr::left_join(scores, by = c("eGene" = "symbol","GWAS" = "disease")) %>%
  dplyr::mutate(gwasCredibleSets = case_when(gwasCredibleSets == "No data" ~ "0", # if there is no data, assume there is no genetic evidence
                                             is.na(gwasCredibleSets) ~ "0", # if the gene is not in the otar scores list
                                             .default = gwasCredibleSets)) %>%
  dplyr::mutate(gwasCredibleSets = as.numeric(gwasCredibleSets))

colocs %>% 
readr::write_csv("../../data/results/8.colocalisation_analysis/coloc_results/annotated_colocs.csv")
# how many of our coloc results match high genetic evidence
table = colocs %>%
  mutate(
    above_threshold_0.2 = gwasCredibleSets > 0.2
  ) %>%
  count(GWAS, above_threshold_0.2) %>%
  group_by(GWAS) %>%
  mutate(
    percent = n / sum(n) * 100
  )

table
# GWAS  above_thresh     n percent
# <chr> <lgl>        <int>   <dbl>
#   1 AD    FALSE           15    55.6
# 2 AD    TRUE            12    44.4
# 3 ALS   TRUE             4   100  
# 4 MS    FALSE           20    48.8
# 5 MS    TRUE            21    51.2
# 6 PD    FALSE            2    28.6
# 7 PD    TRUE             5    71.4

table %>% 
  readr::write_csv("../../data/results/8.colocalisation_analysis/coloc_results/annotated_colocs_percentages.csv")
