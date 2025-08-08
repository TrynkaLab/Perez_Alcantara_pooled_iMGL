# eQTL to CELLECT
library(tidyverse)
library(patchwork)
library(annotate)
library(org.Hs.eg.db)
source("../helpers.R")
options(scipen = 999)
# between treatments
# lfsr ( local false  sign rate): condition-specific measure of significance for 
# each effect, the, which is analogous to a false discovery rate, but more stringent 
# because it requires true discoveries to be not only nonzero, but also correctly 
# signed
eGenes = read_rds("../../../../OTAR2065_sc_eQTL/data/results/8.4.Shared_eqtls_mashr/mashr_results.rds")
eGenes_lfsr = eGenes$result$lfsr
eGenes_lfsr = as.data.frame(eGenes_lfsr)
eGenes_lfsr$eGene_variant = rownames(eGenes_lfsr)

hist(-log10(eGenes_lfsr$`73_IFN_Not_proliferating`))
summary(eGenes_lfsr$`73_IFN_Not_proliferating`)
eGenes_lfsr = eGenes_lfsr %>%
  tidyr::pivot_longer(cols = -eGene_variant, 
                      names_to = "group",
                      values_to = "lfsr")
# build a sort of t-statistic from -log10(lfsr)*beta (aka osterior mean)
eGenes_postmean = eGenes$result$PosteriorMean
eGenes_postmean = as.data.frame(eGenes_postmean)
eGenes_postmean$eGene_variant = rownames(eGenes_postmean)
hist(eGenes_postmean$`73_IFN_Not_proliferating`)
summary(eGenes_postmean$`73_IFN_Not_proliferating`)

eGenes = eGenes_postmean %>%
  tidyr::pivot_longer(cols = -eGene_variant, 
                      names_to = "group",
                      values_to = "PosteriorMean") %>%
  dplyr::left_join(eGenes_lfsr,by = c("eGene_variant","group")) %>%
  # convert 0 pvalues to very small p-values
  dplyr::mutate(lfsr = case_when(lfsr == 0 ~ unique(sort(eGenes_lfsr$lfsr))[2],
                                 .default = lfsr)) %>%
  # create pseudo t -stat
  dplyr::mutate(t = -log10(lfsr)*PosteriorMean)


# getting to the right format for CELLECT
t = eGenes %>%
  as_tibble() %>%
  dplyr::mutate(group = case_when( group == '63_untreated_Not_proliferating' ~ "untreated",
                                   group =="84_LPS_Not_proliferating" ~ "LPS" ,
                                   group == "73_IFN_Not_proliferating" ~ "IFN"),
                gene = str_split_i(eGene_variant,pattern = "-",i=1)) %>%
  # some eGenes are repeated: take average min t of eQTLs
  dplyr::group_by(gene,group) %>%
  dplyr::mutate(t = mean(t)) %>%
  ungroup() %>%
  dplyr::group_by(group) %>%
  dplyr::distinct(gene,.keep_all = TRUE) %>%
  dplyr::mutate(scaled_t = center_this(t)) %>%
  dplyr::mutate(scaled_t_flipped = -center_this(t)) %>%
  dplyr::mutate(normalised_t_top_dir = (scaled_t - min(scaled_t)) / (max(scaled_t) - min(scaled_t)),
                normalised_t_bottom_dir = (scaled_t_flipped - min(scaled_t_flipped)) / (max(scaled_t_flipped) - min(scaled_t_flipped))) %>%
  # in the normalised with directionality, change values for which the t is negative (or positive for bottom dir) to 0
  # so that they are not counted by CELLECT towards the enrichment
  dplyr::mutate(normalised_t_top_dir = case_when(t < 0 ~ 0,
                                                 .default = normalised_t_top_dir),
                normalised_t_bottom_dir = case_when(t > 0 ~ 0,
                                                    .default = normalised_t_bottom_dir)
  )  %>%
  ungroup()


p1 = ggplot(t, aes(x=normalised_t_top_dir , fill = group)) + 
  geom_histogram(alpha = 0.8) +
  theme_minimal() +
  facet_wrap(~group)

p2 = ggplot(t, aes(x=normalised_t_top_dir,y=t, col = group)) + 
  geom_point() +
  theme_minimal() +
  facet_wrap(~group)

p3 = ggplot(t, aes(x=normalised_t_bottom_dir , fill = group)) + 
  geom_histogram(alpha = 0.8) +
  theme_minimal() +
  facet_wrap(~group)

p4 = ggplot(t, aes(x=normalised_t_bottom_dir,y=t, col = group)) + 
  geom_point() +
  theme_minimal() +
  facet_wrap(~group)

png("../../../data/results/CELLECT/treatments_eGenes_t_stat_for_cellect.png",
    width = 10,height = 10,units = "in",res = 400)
plot((p1 +p3) / (p2 + p4))
dev.off()

t %>%
  dplyr::select(gene, group, normalised_t_top_dir) %>%
  tidyr::pivot_wider(names_from = c(group),
                     values_from =normalised_t_top_dir) %>%
  tidyr::unnest(cols = c(untreated, LPS, IFN)) %>%
  replace_na(list(untreated = 0, LPS = 0, IFN = 0)) %>% # filling missing genes per treatment with 0 (no eQTL)
  write_csv("../../../data/results/CELLECT/treatments_eGenes_topdir_t_stat_for_cellect.csv.gz")

t %>%
  dplyr::select(gene, group, normalised_t_bottom_dir) %>%
  tidyr::pivot_wider(names_from = c(group),
                     values_from =normalised_t_bottom_dir) %>%
  tidyr::unnest(cols = c(untreated, LPS, IFN)) %>%
  replace_na(list(untreated = 0, LPS = 0, IFN = 0)) %>% # filling missing genes per treatment with 0 (no eQTL)
  write_csv("../../../data/results/CELLECT/treatments_eGenes_bottomdir_t_stat_for_cellect.csv.gz")

