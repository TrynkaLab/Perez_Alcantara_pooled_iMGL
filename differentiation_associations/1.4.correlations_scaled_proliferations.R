# correlations between shared donors and clones
library(tidyverse)
library(ggpubr)
library(patchwork)
####

proliferation_premac_microglia = read_csv("../../data/results/1.2.scale_proliferation/line_prop_changes_microglia_premac.csv")
proliferation_premac_aging = read_csv("../../data/results/1.2.scale_proliferation/line_prop_changes_old_vs_young_premac.csv")
proliferation_ipsc_premac = read_csv("../../data/results/1.2.scale_proliferation/line_prop_changes_premac_iPSC.csv")

donor_consistency = list()
clone_consistency = list()


#### checking variation between pools for same donors 
# extract shared donors to run this calculation
proliferation_premac_microglia_shared =  proliferation_premac_microglia %>%
  distinct(treatment, line, pool) %>%   # remove duplicate rows if any
  group_by(treatment, line) %>%
  summarise(n_pools = n_distinct(pool), .groups = "drop") %>%
  filter(n_pools >= 2) %>%
  pull(line) %>% unlist() %>% unique()

# Are proliferation effects from the same line across pools more similar than those from different lines?

donor_means <- proliferation_premac_microglia %>% # not subsetting to only lines that are shared so the between-line estimates are more accurate
  group_by(treatment, line, pool,type) %>%
  summarise(mu = mean(scaled_log_fraction), .groups = "drop")
# here we are taking the means for lines because there are several experiments with WGS

within_donor_sd <- donor_means %>%
  group_by(treatment, line, type) %>%
  summarise(sd_within = sd(mu), .groups = "drop")
# sd_within = how much a given line varies across pools
# Only defined for lines in ≥2 pools

donor_overall_means <- donor_means %>%
  group_by(treatment, line, type) %>%
  summarise(mu_donor = mean(mu), .groups = "drop")
# how much each line varies  across pools

donor_consistency[["microglia_vs_premac"]] <- within_donor_sd %>%
  left_join(donor_overall_means, by = c("treatment", "line","type")) %>%
  group_by(type) %>%
  summarise(
    between_donor_sd = sd(mu_donor, na.rm = TRUE),
    mean_within_sd   = mean(sd_within, na.rm = TRUE),
    n_lines          = n_distinct(line),
    .groups = "drop"
  ) %>%
  mutate(comparison = "microglia_vs_premac")

# same for other stages

donor_means <- proliferation_ipsc_premac %>% # not subsetting to only lines that are shared so the between-line estimates are more accurate
  group_by( line, pool,type) %>%
  summarise(mu = mean(scaled_log_fraction), .groups = "drop")
# here we are taking the means for lines because there are several experiments with WGS

within_donor_sd <- donor_means %>%
  group_by(line, type) %>%
  summarise(sd_within = sd(mu), .groups = "drop")
# sd_within = how much a given line varies across pools
# Only defined for lines in ≥2 pools

donor_overall_means <- donor_means %>%
  group_by( line, type) %>%
  summarise(mu_donor = mean(mu), .groups = "drop")
# how much each line varies  across pools

donor_consistency[["premac_vs_ipsc"]] <- within_donor_sd %>%
  left_join(donor_overall_means, by = c( "line","type")) %>%
  group_by(type) %>%
  summarise(
    between_donor_sd = sd(mu_donor, na.rm = TRUE),
    mean_within_sd   = mean(sd_within, na.rm = TRUE),
    n_lines          = n_distinct(line),
    .groups = "drop"
  ) %>%
  mutate(comparison = "premac_vs_ipsc")

# and premac
donor_means <- proliferation_premac_aging %>% # not subsetting to only lines that are shared so the between-line estimates are more accurate
  mutate(type = "WGS") %>%
  group_by( line, pool,type) %>%
  summarise(mu = mean(scaled_log_fraction), .groups = "drop")
# here we are taking the means for lines because there are several experiments with WGS

within_donor_sd <- donor_means %>%
  group_by(line, type) %>%
  summarise(sd_within = sd(mu), .groups = "drop")
# sd_within = how much a given line varies across pools
# Only defined for lines in ≥2 pools

donor_overall_means <- donor_means %>%
  group_by( line, type) %>%
  summarise(mu_donor = mean(mu), .groups = "drop")
# how much each line varies  across pools

donor_consistency[["old_vs_young_premac"]] <- within_donor_sd %>%
  left_join(donor_overall_means, by = c( "line","type")) %>%
  group_by(type) %>%
  summarise(
    between_donor_sd = sd(mu_donor, na.rm = TRUE),
    mean_within_sd   = mean(sd_within, na.rm = TRUE),
    n_lines          = n_distinct(line),
    .groups = "drop"
  ) %>%
  mutate(comparison = "old_vs_young_premac")
#### now checking if lines from clones also vary less than the average

clones_of_interest <- c(
  "lizq_3","lizq1","zaie_5","zaie_1","letw_5","letw_1",
  "seru_1","seru_7","romx_1","romx_2","romx_1",
  "qonc_2","qonc_1","sebn_4","sebn_3"
)

donor_means <- proliferation_premac_microglia %>%
  group_by( line, pool, type) %>%
  summarise(mu = mean(scaled_log_fraction), .groups = "drop")

within_donor_sd_clones <- donor_means %>%
  filter(line %in% clones_of_interest) %>%   # only your selected clones
  # create clone groupings
  dplyr::mutate(clone_group = str_remove_all(line, "_.*$")) %>%
  group_by( clone_group, type) %>%
  summarise(sd_within = sd(mu), .groups = "drop")
# same as before - sd_within = how much a given line varies across pools 

clone_consistency[["microglia_vs_premac"]] <- within_donor_sd_clones %>%
  group_by( type) %>%
  summarise(
    mean_within_sd_clones   = mean(sd_within, na.rm = TRUE),
    n_lines                 = n_distinct(clone_group)*2, # 2 clones for each line
    .groups = "drop"
  ) %>%
  mutate(comparison = "microglia_vs_premac")

# premac vs iPSC
donor_means <- proliferation_ipsc_premac %>%
  group_by( line, pool, type) %>%
  summarise(mu = mean(scaled_log_fraction), .groups = "drop")

within_donor_sd_clones <- donor_means %>%
  filter(line %in% clones_of_interest) %>%   # only your selected clones
  # create clone groupings
  dplyr::mutate(clone_group = str_remove_all(line, "_.*$")) %>%
  group_by( clone_group, type) %>%
  summarise(sd_within = sd(mu), .groups = "drop")
# same as before - sd_within = how much a given line varies across pools 

clone_consistency[["premac_vs_ipsc"]] <- within_donor_sd_clones %>%
  group_by( type) %>%
  summarise(
    mean_within_sd_clones   = mean(sd_within, na.rm = TRUE),
    n_lines                 = n_distinct(clone_group)*2, # 2 clones for each line
    .groups = "drop"
  ) %>%
  mutate(comparison = "premac_vs_ipsc")

#  old vs young premac
donor_means <- proliferation_premac_aging %>%
  mutate(type = "WGS") %>%
  group_by( line, pool, type) %>%
  summarise(mu = mean(scaled_log_fraction), .groups = "drop")

within_donor_sd_clones <- donor_means %>%
  filter(line %in% clones_of_interest) %>%   # only your selected clones
  # create clone groupings
  dplyr::mutate(clone_group = str_remove_all(line, "_.*$")) %>%
  group_by( clone_group, type) %>%
  summarise(sd_within = sd(mu), .groups = "drop")
# same as before - sd_within = how much a given line varies across pools 

clone_consistency[["old_vs_young_premac"]] <- within_donor_sd_clones %>%
  group_by( type) %>%
  summarise(
    mean_within_sd_clones   = mean(sd_within, na.rm = TRUE),
    n_lines                 = n_distinct(clone_group)*2, # 2 clones for each line
    .groups = "drop"
  ) %>%
  mutate(comparison = "old_vs_young_premac")


### joining global and clone consistency
comparison <- do.call("rbind",donor_consistency) %>%
  left_join(do.call("rbind",clone_consistency), by = c( "type","comparison")) %>%
  group_by(comparison) %>%
  summarise(
    between_donor_sd_global = mean(between_donor_sd),
    mean_within_sd   = mean(mean_within_sd),
    mean_within_sd_clones   = mean(mean_within_sd_clones),
    .groups = "drop"
  ) %>%
  select(comparison,between_donor_sd_global, mean_within_sd,
         mean_within_sd_clones)

comparison
mean(comparison$mean_within_sd_clones)
mean(comparison$mean_within_sd)
mean(comparison$between_donor_sd_global)

## mean clone SD (0.71) is slightly larger than the same line SD across pools (0.61)
# as expected
# but lower than the SD of all donors across pools (0.92)

plot_df <- comparison %>%
  pivot_longer(
    cols = c(between_donor_sd_global,
             mean_within_sd,
             mean_within_sd_clones),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    dataset = ifelse(grepl("clones", metric), "Clones only", "All lines"),
    sd_type = case_when(
      metric == "between_donor_sd_global" ~ "SD across all\n donors and pools",
      metric == "mean_within_sd" ~ "mean SD within\n donors across pools",
      metric == "mean_within_sd_clones" ~ "mean SD within\n clones across pools"
    )
  ) %>%
  mutate(comparison = case_when(comparison == "microglia_vs_premac" ~ "macrophage precursor\n to iMGL",
                                comparison == "old_vs_young_premac" ~ "young to aged \n macrophage precursor",
                                comparison == "premac_vs_ipsc" ~ "iPSC to macrophage \n precursor"),
         comparison = factor(comparison, levels = c("iPSC to macrophage \n precursor",
                                                    "young to aged \n macrophage precursor",
                                                    "macrophage precursor\n to iMGL"
                                                    
                                                    ), ordered = TRUE)
  )

p_shared <- ggplot(plot_df, aes(x = comparison)) +
  
  geom_point(
    aes(y = value,
        color = sd_type
  ),size = 3) +
  
  scale_color_manual(
    values = c(
      "SD across all\n donors and pools" ="#F87575" ,
      "mean SD within\n donors across pools" = "#922D50",
      "mean SD within\n clones across pools" = "#3C1B43"
    ),
    name = ""
  ) +
  
  labs(
    x = NULL,
    y = "SD of differentiation efficiency"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "top",
    legend.title = element_text(face = "bold")
  )

p_shared

png("../../data/results/1.4.correlations_scaled_proliferation/SD_shared_donors_across_pools.png",
    width = 6,height = 3,units = "in",res = 400)
plot(p_shared)
dev.off()
### histograms
p1 = proliferation_premac_microglia %>%
ggplot(aes(x = scaled_log_fraction, fill = treatment)) +
  geom_vline(xintercept = 0,linetype = "dashed") + 
  geom_histogram(
    aes(y = after_stat(count)),
    bins = 40,
    color = "grey20",
    alpha = 0.3
  ) +
  theme_pubr() +
  xlim(-5, 5) +
  ylim(0, 320)+
  ggtitle("macrophage precursor to iMGL") + 
  xlab("Differentiation efficiency")
p1
p2 = proliferation_premac_aging %>%
  ggplot(aes(x = scaled_log_fraction)) +
    geom_vline(xintercept = 0,linetype = "dashed") + 
    
  geom_histogram(
    aes(y = after_stat(count)),
    bins = 40,
    color = "grey20",
    fill = "white",
    alpha = 0.3
  ) +
  xlim(-5, 5) +
  ylim(0, 100) +
  ggtitle("young to aged macrophage precursor") +
  theme_pubr() + 
  theme(axis.title.x = element_blank())
p2


p3 = proliferation_ipsc_premac %>%
  ggplot(aes(x = scaled_log_fraction)) +
  geom_vline(xintercept = 0,linetype = "dashed") + 
  
  geom_histogram(
    aes(y = after_stat(count)),
    bins = 40,
    color = "grey20",
    fill = "white",
    alpha = 0.3
  ) +
  xlim(-5, 5) +
  ylim(0, 100) +
  ggtitle("iPSC to macrophage precursor") +
  theme_pubr() +
  theme(axis.title.x = element_blank())


p3

png("../../data/results/1.4.correlations_scaled_proliferation/histogram_differentiation_efficiency.png",
    width = 5,height = 7,units = "in",res = 400)
(p3 / p2 / p1 ) 
dev.off()


