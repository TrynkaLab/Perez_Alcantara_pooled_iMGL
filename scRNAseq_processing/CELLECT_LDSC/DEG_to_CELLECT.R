# converting DEG to CELLECT format
# 0-1
library(tidyverse)
library(patchwork)
library(annotate)
library(org.Hs.eg.db)
source("../helpers.R")

# between treatments
deg = list()
for(cont in c("IFNvsUntreated","LPSvsUntreated","IFNvsLPS")){
  deg[[cont]] = read_tsv(paste0("../../../data/results/5.1.Diff_expr_limma/DiffExpr_",
                                cont,
                                "all_genes_negative_lower_in_latter.txt"))
  deg[[cont]]$contrast = cont
}

deg = do.call("rbind",deg)

deg = deg %>%
  dplyr::mutate(category = case_when(contrast == "IFNvsUntreated" & logFC>0 ~ "higher_IFN_vs_Untreated",
                                     contrast == "IFNvsUntreated" & logFC<0 ~ "higher_Untreated_vs_IFN",
                                     contrast == "LPSvsUntreated" & logFC>0 ~ "higher_LPS_vs_Untreated",
                                     contrast == "LPSvsUntreated" & logFC<0 ~ "higher_Untreated_vs_LPS",
                                     contrast == "IFNvsLPS" & logFC>0 ~ "higher_IFN_vs_LPS",
                                     contrast == "IFNvsLPS" & logFC<0 ~ "higher_LPS_vs_IFN"
  ),
  fold_change_direction = case_when(logFC>0 ~ "positive",
                                    .default = "negative")) %>%
  dplyr::group_by(contrast) %>%
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


# need ensembl gene ids as "gene"

ensembl = symbol_to_ensembl(deg$symbol)
ensembl = enframe(ensembl, name = "symbol", value = "gene")
deg = deg %>%
  dplyr::left_join(.,ensembl) %>%
  distinct() %>%
  dplyr::filter(!is.na(gene))


deg %>% 
  dplyr::select(gene, contrast, normalised_t_top_dir) %>%
  tidyr::pivot_wider(names_from = c(contrast),
                     values_from =normalised_t_top_dir) %>%
  tidyr::unnest(cols = c(IFNvsUntreated, LPSvsUntreated, IFNvsLPS)) %>%
  replace_na(list(IFNvsUntreated = 0, LPSvsUntreated = 0, IFNvsLPS = 0)) %>%
  write_csv("../../../data/results/CELLECT/treatment_DEG_normalised_t_top_dir.csv.gz")


deg %>% 
  dplyr::select(gene, contrast, normalised_t_bottom_dir) %>%
  tidyr::pivot_wider(names_from = c(contrast),
                     values_from =normalised_t_bottom_dir) %>%
  tidyr::unnest(cols = c(IFNvsUntreated, LPSvsUntreated, IFNvsLPS)) %>%
  replace_na(list(IFNvsUntreated = 0, LPSvsUntreated = 0, IFNvsLPS = 0)) %>%
  write_csv("../../../data/results/CELLECT/treatment_DEG_normalised_t_bottom_dir.csv.gz")


p1 = deg %>% 
  dplyr::select(gene, contrast, normalised_t_top_dir) %>%
  
  tidyr::pivot_wider(names_from = c(contrast),
                     values_from =normalised_t_top_dir) %>%
  tidyr::unnest(cols = c(IFNvsUntreated, LPSvsUntreated, IFNvsLPS)) %>%
  replace_na(list(IFNvsUntreated = 0, LPSvsUntreated = 0, IFNvsLPS = 0)) %>%
  tidyr::pivot_longer(cols = -gene,names_to = "contrast", values_to = "normalised_t_top_dir") %>%
  ggplot(., aes(x=normalised_t_top_dir, fill=contrast)) + 
  geom_histogram(alpha = 0.8) +
  theme_minimal() +
  facet_wrap(~contrast)

p2 =  deg %>% 
  dplyr::select(gene, contrast, normalised_t_bottom_dir) %>%
  
  tidyr::pivot_wider(names_from = c(contrast),
                     values_from =normalised_t_bottom_dir) %>%
  tidyr::unnest(cols = c(IFNvsUntreated, LPSvsUntreated, IFNvsLPS)) %>%
  replace_na(list(IFNvsUntreated = 0, LPSvsUntreated = 0, IFNvsLPS = 0)) %>%
  tidyr::pivot_longer(cols = -gene,names_to = "contrast", values_to = "normalised_t_bottom_dir") %>%
  ggplot(., aes(x=normalised_t_bottom_dir, fill=contrast)) + 
  geom_histogram(alpha = 0.8) +
  theme_minimal() +
  facet_wrap(~contrast)

png("../../../data/results/CELLECT/treatments_DEG_t_stat_for_cellect.png",
    width = 10,height = 10,units = "in",res = 400)
plot(p1 / p2)
dev.off()


# PRS DEG #############

deg = list()
for(cont in c("untreated","IFN","LPS")){
  deg[[cont]] = read_csv(paste0("../../../data/results/5.1.2.Diff_expr_limma_high_low_PRS/DiffExpr_PRS_",
                                cont,".csv"))
  deg[[cont]]$contrast = cont
}

deg = do.call("rbind",deg)

deg = deg %>%
  dplyr::mutate(category = case_when( logFC>0 ~ "higher_high_PRS",
                                      logFC<0 ~ "higher_low_PRS"
  )) %>%
  dplyr::group_by(contrast) %>%
  dplyr::mutate(scaled_t = center_this(t)) %>%
  dplyr::mutate(scaled_t_flipped = -center_this(t)) %>%
  dplyr::mutate(normalised_t_higher_risk_dir = (scaled_t - min(scaled_t)) / (max(scaled_t) - min(scaled_t)),
                normalised_t_lower_risk_dir = (scaled_t_flipped - min(scaled_t_flipped)) / (max(scaled_t_flipped) - min(scaled_t_flipped))) %>%
  # in the normalised with directionality, change values for which the t is negative (or positive for bottom dir) to 0
  # so that they are not counted by CELLECT towards the enrichment
  dplyr::mutate(normalised_t_higher_risk_dir = case_when(t < 0 ~ 0,
                                                         .default = normalised_t_higher_risk_dir),
                normalised_t_lower_risk_dir = case_when(t > 0 ~ 0,
                                                        .default = normalised_t_lower_risk_dir)
  )  %>%
  ungroup()


ensembl = symbol_to_ensembl(deg$symbol)
ensembl = enframe(ensembl, name = "symbol", value = "gene")
deg = deg %>%
  dplyr::left_join(.,ensembl) %>%
  distinct() %>%
  dplyr::filter(!is.na(gene))

deg %>% 
  dplyr::select(gene, contrast, normalised_t_higher_risk_dir) %>%
  tidyr::pivot_wider(names_from = c(contrast),
                     values_from =normalised_t_higher_risk_dir) %>%
  tidyr::unnest(cols = c(untreated, IFN, LPS)) %>%
  replace_na(list(untreated = 0, IFN = 0, LPS = 0)) %>%
  write_csv("../../../data/results/CELLECT/PRS_DEG_normalised_t_higher_risk_dir.csv.gz")


deg %>% 
  dplyr::select(gene, contrast, normalised_t_lower_risk_dir) %>%
  
  tidyr::pivot_wider(names_from = c(contrast),
                     values_from =normalised_t_lower_risk_dir) %>%
  tidyr::unnest(cols = c(untreated, IFN, LPS)) %>%
  replace_na(list(untreated = 0, IFN = 0, LPS = 0)) %>%
  write_csv("../../../data/results/CELLECT/PRS_DEG_normalised_t_lower_risk_dir.csv.gz")

p1 = deg %>% 
  dplyr::select(gene, contrast, normalised_t_higher_risk_dir) %>%
  tidyr::pivot_wider(names_from = c(contrast),
                     values_from =normalised_t_higher_risk_dir) %>%
  tidyr::unnest(cols = c(untreated, IFN, LPS)) %>%
  replace_na(list(untreated = 0, IFN = 0, LPS = 0)) %>%
  tidyr::pivot_longer(cols = -gene,names_to = "contrast", values_to = "normalised_t_higher_risk_dir") %>%
  ggplot(., aes(x=normalised_t_higher_risk_dir, fill=contrast)) + 
  geom_histogram(alpha = 0.8) +
  theme_minimal() +
  facet_wrap(~contrast)

p2 =  deg %>% 
  dplyr::select(gene, contrast, normalised_t_lower_risk_dir) %>%
  
  tidyr::pivot_wider(names_from = c(contrast),
                     values_from =normalised_t_lower_risk_dir) %>%
  tidyr::unnest(cols = c(untreated, IFN, LPS)) %>%
  replace_na(list(untreated = 0, IFN = 0, LPS = 0)) %>%
  tidyr::pivot_longer(cols = -gene,names_to = "contrast", values_to = "normalised_t_lower_risk_dir") %>%
  ggplot(., aes(x=normalised_t_lower_risk_dir, fill=contrast)) + 
  geom_histogram(alpha = 0.8) +
  theme_minimal() +
  facet_wrap(~contrast)

png("../../../data/results/CELLECT/PRS_DEG_t_stat_for_cellect.png",
    width = 10,height = 10,units = "in",res = 400)
plot(p1 / p2)
dev.off()

# APOE HR (high - e4 vs medium (e3) vs low - e2) DEG #############

deg = list()
for(cont in c("untreated","IFN","LPS")){
  deg[[cont]] = read_csv(paste0("../../../data/results/5.1.2.Diff_expr_limma_high_low_PRS/DiffExpr_apoeHR_PRS_",
                                cont,".csv"))
  deg[[cont]]$contrast = cont
}

deg = do.call("rbind",deg)

deg = deg %>%
  dplyr::mutate(category = case_when( logFC>0 ~ "higher_high_APOE_PRS",
                                      logFC<0 ~ "higher_low_APOE_PRS"
  )) %>%
  dplyr::group_by(contrast) %>%
  dplyr::mutate(scaled_t = center_this(t)) %>%
  dplyr::mutate(scaled_t_flipped = -center_this(t)) %>%
  dplyr::mutate(normalised_t_higher_risk_dir = (scaled_t - min(scaled_t)) / (max(scaled_t) - min(scaled_t)),
                normalised_t_lower_risk_dir = (scaled_t_flipped - min(scaled_t_flipped)) / (max(scaled_t_flipped) - min(scaled_t_flipped))) %>%
  # in the normalised with directionality, change values for which the t is negative (or positive for bottom dir) to 0
  # so that they are not counted by CELLECT towards the enrichment
  dplyr::mutate(normalised_t_higher_risk_dir = case_when(t < 0 ~ 0,
                                                         .default = normalised_t_higher_risk_dir),
                normalised_t_lower_risk_dir = case_when(t > 0 ~ 0,
                                                        .default = normalised_t_lower_risk_dir)
  )  %>%
  ungroup()




# need ensembl gene ids as "gene"

ensembl = symbol_to_ensembl(deg$symbol)
ensembl = enframe(ensembl, name = "symbol", value = "gene")
deg = deg %>%
  dplyr::left_join(.,ensembl) %>%
  distinct() %>%
  dplyr::filter(!is.na(gene))

deg %>% 
  dplyr::select(gene, contrast, normalised_t_higher_risk_dir) %>%
  tidyr::pivot_wider(names_from = c(contrast),
                     values_from =normalised_t_higher_risk_dir) %>%
  tidyr::unnest(cols = c(untreated, IFN, LPS)) %>%
  replace_na(list(untreated = 0, IFN = 0, LPS = 0)) %>%
  write_csv("../../../data/results/CELLECT/APOE_DEG_normalised_t_higher_risk_dir.csv.gz")


deg %>% 
  dplyr::select(gene, contrast, normalised_t_lower_risk_dir) %>%
  
  tidyr::pivot_wider(names_from = c(contrast),
                     values_from =normalised_t_lower_risk_dir) %>%
  tidyr::unnest(cols = c(untreated, IFN, LPS)) %>%
  replace_na(list(untreated = 0, IFN = 0, LPS = 0)) %>%
  write_csv("../../../data/results/CELLECT/APOE_DEG_normalised_t_lower_risk_dir.csv.gz")

p1 = deg %>% 
  dplyr::select(gene, contrast, normalised_t_higher_risk_dir) %>%
  tidyr::pivot_wider(names_from = c(contrast),
                     values_from =normalised_t_higher_risk_dir) %>%
  tidyr::unnest(cols = c(untreated, IFN, LPS)) %>%
  replace_na(list(untreated = 0, IFN = 0, LPS = 0)) %>%
  tidyr::pivot_longer(cols = -gene,names_to = "contrast", values_to = "normalised_t_higher_risk_dir") %>%
  ggplot(., aes(x=normalised_t_higher_risk_dir, fill=contrast)) + 
  geom_histogram(alpha = 0.8) +
  theme_minimal() +
  facet_wrap(~contrast)

p2 = deg %>% 
  dplyr::select(gene, contrast, normalised_t_lower_risk_dir) %>%
  
  tidyr::pivot_wider(names_from = c(contrast),
                     values_from =normalised_t_lower_risk_dir) %>%
  tidyr::unnest(cols = c(untreated, IFN, LPS)) %>%
  replace_na(list(untreated = 0, IFN = 0, LPS = 0)) %>%
  tidyr::pivot_longer(cols = -gene,names_to = "contrast", values_to = "normalised_t_lower_risk_dir") %>%
  ggplot(., aes(x=normalised_t_lower_risk_dir, fill=contrast)) + 
  geom_histogram(alpha = 0.8) +
  theme_minimal() +
  facet_wrap(~contrast)

png("../../../data/results/CELLECT/APOE_DEG_t_stat_for_cellect.png",
    width = 10,height = 10,units = "in",res = 400)
plot(p1 / p2)
dev.off()
#### polygenic component only


deg = list()
for(cont in c("untreated","IFN","LPS")){
  deg[[cont]] = read_csv(paste0("../../../data/results/5.1.2.Diff_expr_limma_high_low_PRS/DiffExpr_prs_scaled_PRS_",
                                cont,".csv"))
  deg[[cont]]$contrast = cont
}

deg = do.call("rbind",deg)

deg = deg %>%
  dplyr::mutate(category = case_when( logFC>0 ~ "higher_high_polygenic_PRS",
                                      logFC<0 ~ "higher_low_polygenic_PRS"
  )) %>%
  dplyr::group_by(contrast) %>%
  dplyr::mutate(scaled_t = center_this(t)) %>%
  dplyr::mutate(scaled_t_flipped = -center_this(t)) %>%
  dplyr::mutate(normalised_t_higher_risk_dir = (scaled_t - min(scaled_t)) / (max(scaled_t) - min(scaled_t)),
                normalised_t_lower_risk_dir = (scaled_t_flipped - min(scaled_t_flipped)) / (max(scaled_t_flipped) - min(scaled_t_flipped))) %>%
  # in the normalised with directionality, change values for which the t is negative (or positive for bottom dir) to 0
  # so that they are not counted by CELLECT towards the enrichment
  dplyr::mutate(normalised_t_higher_risk_dir = case_when(t < 0 ~ 0,
                                                         .default = normalised_t_higher_risk_dir),
                normalised_t_lower_risk_dir = case_when(t > 0 ~ 0,
                                                        .default = normalised_t_lower_risk_dir)
  )  %>%
  ungroup()


# need ensembl gene ids as "gene"

ensembl = symbol_to_ensembl(deg$symbol)
ensembl = enframe(ensembl, name = "symbol", value = "gene")
deg = deg %>%
  dplyr::left_join(.,ensembl) %>%
  distinct() %>%
  dplyr::filter(!is.na(gene))



deg %>% 
  dplyr::select(gene, contrast, normalised_t_higher_risk_dir) %>%
  tidyr::pivot_wider(names_from = c(contrast),
                     values_from =normalised_t_higher_risk_dir) %>%
  tidyr::unnest(cols = c(untreated, IFN, LPS)) %>%
  replace_na(list(untreated = 0, IFN = 0, LPS = 0)) %>%
  write_csv("../../../data/results/CELLECT/polygenicHR_DEG_normalised_t_higher_risk_dir.csv.gz")


deg %>% 
  dplyr::select(gene, contrast, normalised_t_lower_risk_dir) %>%
  tidyr::pivot_wider(names_from = c(contrast),
                     values_from =normalised_t_lower_risk_dir) %>%
  tidyr::unnest(cols = c(untreated, IFN, LPS)) %>%
  replace_na(list(untreated = 0, IFN = 0, LPS = 0)) %>%
  write_csv("../../../data/results/CELLECT/polygenicHR_DEG_normalised_t_lower_risk_dir.csv.gz")



p1 = deg %>% 
  dplyr::select(gene, contrast, normalised_t_higher_risk_dir) %>%
  tidyr::pivot_wider(names_from = c(contrast),
                     values_from =normalised_t_higher_risk_dir) %>%
  tidyr::unnest(cols = c(untreated, IFN, LPS)) %>%
  replace_na(list(untreated = 0, IFN = 0, LPS = 0)) %>%
  tidyr::pivot_longer(cols = -gene,names_to = "contrast", values_to = "normalised_t_higher_risk_dir") %>%
  ggplot(., aes(x=normalised_t_higher_risk_dir, fill=contrast)) + 
  geom_histogram(alpha = 0.8) +
  theme_minimal() +
  facet_wrap(~contrast)

p2 =  deg %>% 
  dplyr::select(gene, contrast, normalised_t_lower_risk_dir) %>%
  tidyr::pivot_wider(names_from = c(contrast),
                     values_from =normalised_t_lower_risk_dir) %>%
  tidyr::unnest(cols = c(untreated, IFN, LPS)) %>%
  replace_na(list(untreated = 0, IFN = 0, LPS = 0)) %>%
  tidyr::pivot_longer(cols = -gene,names_to = "contrast", values_to = "normalised_t_lower_risk_dir") %>%
  ggplot(., aes(x=normalised_t_lower_risk_dir, fill=contrast)) + 
  geom_histogram(alpha = 0.8) +
  theme_minimal() +
  facet_wrap(~contrast)


png("../../../data/results/CELLECT/polygenicHR_DEG_t_stat_for_cellect.png",
    width = 10,height = 10,units = "in",res = 400)
plot(p1 / p2)
dev.off()


#### without directionality ########
###################################
###################################
###################################
###################################
###################################

library(tidyverse)
library(patchwork)
library(annotate)
library(org.Hs.eg.db)
source("../helpers.R")

# between treatments
deg = list()
for(cont in c("IFNvsUntreated","LPSvsUntreated","IFNvsLPS")){
  deg[[cont]] = read_tsv(paste0("../../../data/results/5.1.Diff_expr_limma/DiffExpr_",
                                cont,
                                "all_genes_negative_lower_in_latter.txt"))
  deg[[cont]]$contrast = cont
}

deg = do.call("rbind",deg)

deg = deg %>%
  dplyr::group_by(contrast) %>%
  dplyr::mutate(scaled_t = center_this(abs(t))) %>%
  dplyr::mutate(normalised_t = (scaled_t - min(scaled_t)) / (max(scaled_t) - min(scaled_t)) ) %>%
  ungroup()

# need ensembl gene ids as "gene"

ensembl = symbol_to_ensembl(deg$symbol)
ensembl = enframe(ensembl, name = "symbol", value = "gene")
deg = deg %>%
  dplyr::left_join(.,ensembl) %>%
  distinct() %>%
  dplyr::filter(!is.na(gene))


deg %>% 
  dplyr::select(gene, contrast, normalised_t) %>%
  tidyr::pivot_wider(names_from = c(contrast),
                     values_from =normalised_t) %>%
  tidyr::unnest(cols = c(IFNvsUntreated, LPSvsUntreated, IFNvsLPS)) %>%
  replace_na(list(IFNvsUntreated = 0, LPSvsUntreated = 0, IFNvsLPS = 0)) %>%
  write_csv("../../../data/results/CELLECT/treatment_DEG_normalised_t.csv.gz")


p1 = deg %>% 
  dplyr::select(gene, contrast, normalised_t) %>%
  
  tidyr::pivot_wider(names_from = c(contrast),
                     values_from =normalised_t) %>%
  tidyr::unnest(cols = c(IFNvsUntreated, LPSvsUntreated, IFNvsLPS)) %>%
  replace_na(list(IFNvsUntreated = 0, LPSvsUntreated = 0, IFNvsLPS = 0)) %>%
  tidyr::pivot_longer(cols = -gene,names_to = "contrast", values_to = "normalised_t") %>%
  ggplot(., aes(x=normalised_t, fill=contrast)) + 
  geom_histogram(alpha = 0.8) +
  theme_minimal() +
  facet_wrap(~contrast)


png("../../../data/results/CELLECT/treatments_DEG_t_stat_for_cellect_no_dir.png",
    width = 10,height = 5,units = "in",res = 400)
plot(p1 )
dev.off()


# PRS DEG #############

deg = list()
for(cont in c("untreated","IFN","LPS")){
  deg[[cont]] = read_csv(paste0("../../../data/results/5.1.2.Diff_expr_limma_high_low_PRS/DiffExpr_PRS_",
                                cont,".csv"))
  deg[[cont]]$contrast = cont
}

deg = do.call("rbind",deg)

deg = deg %>%
  dplyr::mutate(category = case_when( logFC>0 ~ "higher_high_PRS",
                                      logFC<0 ~ "higher_low_PRS"
  )) %>%
  dplyr::group_by(contrast) %>%
  dplyr::mutate(scaled_t = center_this(abs(t))) %>%
  dplyr::mutate(normalised_t = (scaled_t - min(scaled_t)) / (max(scaled_t) - min(scaled_t))) %>%
  ungroup()


ensembl = symbol_to_ensembl(deg$symbol)
ensembl = enframe(ensembl, name = "symbol", value = "gene")
deg = deg %>%
  dplyr::left_join(.,ensembl) %>%
  distinct() %>%
  dplyr::filter(!is.na(gene))

deg %>% 
  dplyr::select(gene, contrast, normalised_t) %>%
  tidyr::pivot_wider(names_from = c(contrast),
                     values_from =normalised_t) %>%
  tidyr::unnest(cols = c(untreated, IFN, LPS)) %>%
  replace_na(list(untreated = 0, IFN = 0, LPS = 0)) %>%
  write_csv("../../../data/results/CELLECT/PRS_DEG_normalised_t.csv.gz")


p1 = deg %>% 
  dplyr::select(gene, contrast, normalised_t) %>%
  tidyr::pivot_wider(names_from = c(contrast),
                     values_from =normalised_t) %>%
  tidyr::unnest(cols = c(untreated, IFN, LPS)) %>%
  replace_na(list(untreated = 0, IFN = 0, LPS = 0)) %>%
  tidyr::pivot_longer(cols = -gene,names_to = "contrast", values_to = "normalised_t") %>%
  ggplot(., aes(x=normalised_t, fill=contrast)) + 
  geom_histogram(alpha = 0.8) +
  theme_minimal() +
  facet_wrap(~contrast)


png("../../../data/results/CELLECT/PRS_DEG_t_stat_for_cellect_no_dir.png",
    width = 10,height = 5,units = "in",res = 400)
plot(p1)
dev.off()

# APOE HR (high - e4 vs medium (e3) vs low - e2) DEG #############

deg = list()
for(cont in c("untreated","IFN","LPS")){
  deg[[cont]] = read_csv(paste0("../../../data/results/5.1.2.Diff_expr_limma_high_low_PRS/DiffExpr_APOE_sum_scaled_PRS_APOE_",
                                cont,".csv"))
  deg[[cont]]$contrast = cont
}

deg = do.call("rbind",deg)

deg = deg %>%
  dplyr::mutate(category = case_when( logFC>0 ~ "higher_high_APOE_PRS",
                                      logFC<0 ~ "higher_low_APOE_PRS"
  )) %>%
  dplyr::group_by(contrast) %>%
  dplyr::mutate(scaled_t = center_this(abs(t))) %>%
  dplyr::mutate(normalised_t = (scaled_t - min(scaled_t)) / (max(scaled_t) - min(scaled_t))) %>%
  ungroup()




# need ensembl gene ids as "gene"

ensembl = symbol_to_ensembl(deg$symbol)
ensembl = enframe(ensembl, name = "symbol", value = "gene")
deg = deg %>%
  dplyr::left_join(.,ensembl) %>%
  distinct() %>%
  dplyr::filter(!is.na(gene))

deg %>% 
  dplyr::select(gene, contrast, normalised_t) %>%
  tidyr::pivot_wider(names_from = c(contrast),
                     values_from =normalised_t) %>%
  tidyr::unnest(cols = c(untreated, IFN, LPS)) %>%
  replace_na(list(untreated = 0, IFN = 0, LPS = 0)) %>%
  write_csv("../../../data/results/CELLECT/APOE_DEG_normalised_t.csv.gz")



p1 = deg %>% 
  dplyr::select(gene, contrast, normalised_t) %>%
  tidyr::pivot_wider(names_from = c(contrast),
                     values_from =normalised_t) %>%
  tidyr::unnest(cols = c(untreated, IFN, LPS)) %>%
  replace_na(list(untreated = 0, IFN = 0, LPS = 0)) %>%
  tidyr::pivot_longer(cols = -gene,names_to = "contrast", values_to = "normalised_t") %>%
  ggplot(., aes(x=normalised_t, fill=contrast)) + 
  geom_histogram(alpha = 0.8) +
  theme_minimal() +
  facet_wrap(~contrast)


png("../../../data/results/CELLECT/APOE_DEG_t_stat_for_cellect_no_dir.png",
    width = 10,height = 5,units = "in",res = 400)
plot(p1 )
dev.off()
#### polygenic component only


deg = list()
for(cont in c("untreated","IFN","LPS")){
  deg[[cont]] = read_csv(paste0("../../../data/results/5.1.2.Diff_expr_limma_high_low_PRS/DiffExpr_prs_scaled_PRS_",
                                cont,".csv"))
  deg[[cont]]$contrast = cont
}

deg = do.call("rbind",deg)

deg = deg %>%
  dplyr::mutate(category = case_when( logFC>0 ~ "higher_high_polygenic_PRS",
                                      logFC<0 ~ "higher_low_polygenic_PRS"
  )) %>%
  dplyr::group_by(contrast) %>%
  dplyr::mutate(scaled_t = center_this(abs(t))) %>%
  dplyr::mutate(normalised_t = (scaled_t - min(scaled_t)) / (max(scaled_t) - min(scaled_t))) %>%
  ungroup()


# need ensembl gene ids as "gene"

ensembl = symbol_to_ensembl(deg$symbol)
ensembl = enframe(ensembl, name = "symbol", value = "gene")
deg = deg %>%
  dplyr::left_join(.,ensembl) %>%
  distinct() %>%
  dplyr::filter(!is.na(gene))



deg %>% 
  dplyr::select(gene, contrast, normalised_t) %>%
  tidyr::pivot_wider(names_from = c(contrast),
                     values_from =normalised_t) %>%
  tidyr::unnest(cols = c(untreated, IFN, LPS)) %>%
  replace_na(list(untreated = 0, IFN = 0, LPS = 0)) %>%
  write_csv("../../../data/results/CELLECT/polygenicHR_DEG_normalised_t.csv.gz")



p1 = deg %>% 
  dplyr::select(gene, contrast, normalised_t) %>%
  tidyr::pivot_wider(names_from = c(contrast),
                     values_from =normalised_t) %>%
  tidyr::unnest(cols = c(untreated, IFN, LPS)) %>%
  replace_na(list(untreated = 0, IFN = 0, LPS = 0)) %>%
  tidyr::pivot_longer(cols = -gene,names_to = "contrast", values_to = "normalised_t") %>%
  ggplot(., aes(x=normalised_t, fill=contrast)) + 
  geom_histogram(alpha = 0.8) +
  theme_minimal() +
  facet_wrap(~contrast)



png("../../../data/results/CELLECT/polygenicHR_DEG_t_stat_for_cellect_no_dir.png",
    width = 10,height =5,units = "in",res = 400)
plot(p1)
dev.off()

