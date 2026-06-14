# Comparing proliferation efficiency to that of Jerber et al 2021 (dopaminergic neurons)
# Their definition of differentiation efficiency is different
# Neuronal differentiation efficiency was defined as the sum of the proportion of 
# serotonergic-like and dopaminergic neurons present on day 52. 
# The neuronal differentiation efficiency of iPSC lines was calculated as 
# the average of the efficiencies across all pools in which that cell line was included.
# Since we have a homogeneous population, even excluding the proliferative cluster 
# the differentiation efficiency for microglia would be very high

# Trying to make comparable measures of proliferation efficiency
library(tidyverse)
source("functions.R")
library(ggrepel)
library(scales)

### read in files
jerber_cells <- read.csv("../../../resources/Jerber_differentiation_efficicency/jerber_supplementary_t2.csv")

jerber_diff_eff <- read_tsv("../../../resources/Jerber_differentiation_efficicency/neuroseq_differentiation_efficiency.pool1_17_D52.tsv") %>%
  mutate(line = str_remove(donor_id, ".*-")) %>%
  select(-donor_id) %>%
  mutate(scaled_diff_efficiency = scale_this(diff_efficiency))

# ours
iPSC_premac <- read.csv("../../data/results/1.2.scale_proliferation/line_prop_changes_premac_iPSC.csv")
premac_microglia <- read.csv("../../data/results/1.2.scale_proliferation/line_prop_changes_microglia_premac.csv")

#### calculate donor proportions for day 52 from Jerber

# taking all cell types together
jerber_proportions_all_together <- jerber_cells %>%
  group_by(pool_id, time_point, treatment, cell_line) %>%
  summarise(n = sum(n_cells), .groups = "drop") %>%
  group_by(pool_id, time_point, treatment) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# only for neurons: dopaminergic DA and serotonergic Sert
jerber_proportions_neurons <- jerber_cells %>%
  mutate(time_cell = paste(time_point, celltype, sep = "-")) %>%
  filter(time_point == "D11" | time_cell %in% c("D52-DA", "D52-Sert")) %>%
  group_by(pool_id, time_point, treatment, cell_line) %>%
  summarise(n = sum(n_cells), .groups = "drop") %>%
  group_by(pool_id, time_point, treatment) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# calculating proliferation efficiency: all cell types
jerber_proliferation_all_together <- jerber_proportions_all_together %>%
  filter(treatment == "NONE") %>%
  group_by(pool_id, cell_line) %>%
  reframe(
    log_fraction_mean = log(prop[time_point == "D52"] / prop[time_point == "D11"]),
    comparison = "D52_vs_D11"
  ) %>%
  ungroup() %>%
  distinct() %>%
  filter(!log_fraction_mean %in% c("-Inf", "Inf", "NA", "NaN")) %>%
  group_by(pool_id) %>%
  mutate(scaled_log_fraction = scale_this(log_fraction_mean)) %>%
  ungroup() %>%
  mutate(line = str_remove(cell_line, ".*-")) %>%
  group_by(line) %>%
  summarise(jerber_mean_scaled_log_fraction = mean(scaled_log_fraction), .groups = "drop")

# calculating proliferation efficiency: neurons only
jerber_proliferation_neurons <- jerber_proportions_neurons %>%
  filter(treatment == "NONE") %>%
  group_by(pool_id, cell_line) %>%
  reframe(
    log_fraction_mean = log(prop[time_point == "D52"] / prop[time_point == "D11"]),
    comparison = "D52_vs_D11"
  ) %>%
  ungroup() %>%
  distinct() %>%
  filter(!log_fraction_mean %in% c("-Inf", "Inf", "NA", "NaN")) %>%
  group_by(pool_id) %>%
  mutate(scaled_log_fraction = scale_this(log_fraction_mean)) %>%
  ungroup() %>%
  mutate(line = str_remove(cell_line, ".*-")) %>%
  group_by(line) %>%
  summarise(jerber_mean_scaled_log_fraction = mean(scaled_log_fraction), .groups = "drop")

# average for microglia
mean_premac_microglia <- premac_microglia %>%
  filter(treatment == "untreated") %>%
  group_by(line) %>%
  summarise(
    microglia_mean_scaled_log_fraction = mean(scaled_log_fraction),
    .groups = "drop"
  )

# quick distributions
hist(jerber_proliferation_all_together$jerber_mean_scaled_log_fraction)
hist(jerber_proliferation_neurons$jerber_mean_scaled_log_fraction)
hist(jerber_diff_eff$diff_efficiency)
hist(jerber_diff_eff$scaled_diff_efficiency)
hist(mean_premac_microglia$microglia_mean_scaled_log_fraction)

# main plot data: microglia vs Jerber all-cell proliferation
plot_df <- mean_premac_microglia %>%
  inner_join(jerber_proliferation_all_together, by = "line") %>%
  filter(
    !is.na(microglia_mean_scaled_log_fraction),
    !is.na(jerber_mean_scaled_log_fraction)
  ) %>%
  mutate(
    x = microglia_mean_scaled_log_fraction,
    y = jerber_mean_scaled_log_fraction,
    agreement = abs(x - y),
    magnitude = abs(x) + abs(y)
  )

# overlay data: microglia vs Jerber neuron-only proliferation
plot_df_neurons_overlay <- mean_premac_microglia %>%
  inner_join(jerber_proliferation_neurons, by = "line") %>%
  filter(
    !is.na(microglia_mean_scaled_log_fraction),
    !is.na(jerber_mean_scaled_log_fraction)
  ) %>%
  mutate(
    x = microglia_mean_scaled_log_fraction,
    y = jerber_mean_scaled_log_fraction
  )

# overlay data: microglia vs Jerber scaled neuronal differentiation efficiency
plot_df_diff_eff_overlay <- mean_premac_microglia %>%
  inner_join(jerber_diff_eff, by = "line") %>%
  filter(
    !is.na(microglia_mean_scaled_log_fraction),
    !is.na(scaled_diff_efficiency)
  ) %>%
  mutate(
    x = microglia_mean_scaled_log_fraction,
    y = scaled_diff_efficiency
  )

# labels
n_to_label <- 10

label_df <- plot_df %>%
  filter(magnitude >= quantile(magnitude, 0.5, na.rm = TRUE)) %>%
  slice_min(order_by = agreement, n = n_to_label, with_ties = FALSE)

# correlations
ct_all <- cor.test(plot_df$x, plot_df$y, method = "pearson")
ct_neurons <- cor.test(plot_df_neurons_overlay$x, plot_df_neurons_overlay$y, method = "pearson")
ct_diff_eff <- cor.test(plot_df_diff_eff_overlay$x, plot_df_diff_eff_overlay$y, method = "pearson")
label_txt_all <- paste0(
  " Proliferation, all cells: Pearson r = ", round(unname(ct_all$estimate), 2),
  ", p = ", format.pval(ct_all$p.value, digits = 2, eps = 1e-3)
)

label_txt_neurons <- paste0(
  "Proliferation, neurons only: Pearson r = ", round(unname(ct_neurons$estimate), 2),
  ", p = ", format.pval(ct_neurons$p.value, digits = 2, eps = 1e-3)
)

label_txt_diff <- paste0(
  " Diff. efficiency: Pearson r = ", round(unname(ct_diff_eff$estimate), 2),
  ", p = ", format.pval(ct_diff_eff$p.value, digits = 2, eps = 1e-3)
)

# plot
p <- ggplot(plot_df, aes(x = x, y = y)) +
  geom_point(size = 2.5, alpha = 0.75, color = "grey30") +
  
  # Main regression: all Jerber cell types
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 1,
    color = "#2C7FB8",
    fill = "#A6CEE3"
  ) +
  
  # Overlay regression: neuron-only Jerber proliferation
  geom_smooth(
    data = plot_df_neurons_overlay,
    aes(x = x, y = y),
    method = "lm",
    se = FALSE,
    linewidth = 1,
    color = "#D95F02"
  ) +
  
  # Overlay regression: Jerber scaled neuronal differentiation efficiency
  geom_smooth(
    data = plot_df_diff_eff_overlay,
    aes(x = x, y = y),
    method = "lm",
    se = FALSE,
    linewidth = 1,
    color = "#1B9E77"
  ) +
  
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "grey60"
  ) +
  
  geom_label_repel(
    data = label_df,
    aes(label = line),
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.2,
    min.segment.length = 0
  ) +
  

  annotate(
    "text",
    x = -Inf, y = Inf,
    label = label_txt_all,
    hjust = -0.05, vjust = 1.5,
    size = 3.3,
    fontface = "bold",
    color = "#2C7FB8"
  ) +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = label_txt_neurons,
    hjust = -0.05, vjust = 3.0,
    size = 3.3,
    fontface = "bold",
    color = "#D95F02"
  ) +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = label_txt_diff,
    hjust = -0.05, vjust = 4.5,
    size = 3.3,
    fontface = "bold",
    color = "#1B9E77"
  ) +
  labs(
    x = "Microglia mean scaled log fraction",
    y = "Jerber mean scaled measure",
    title = "Association between microglia and Jerber proliferation / differentiation",
    subtitle = paste0(
      "Labels show ", n_to_label, " cell lines with smallest |x - y|"
    )
  ) +
  scale_y_continuous(
    limits = c(NA, max(plot_df$y, na.rm = TRUE) + 0.3)
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

plot(p)

# save
out_dir <- "../../data/results/review_compare_proliferation_to_dopaminergic_neuron_Jerber_2021"

dir.create(
  out_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

png(
  file.path(out_dir, "proliferation_efficiency_correlation.png"),
  width = 5,
  height = 5.5,
  units = "in",
  res = 1000
)
plot(p)
dev.off()

pdf(
  file.path(out_dir, "proliferation_efficiency_correlation.pdf"),
  width = 5,
  height = 5.5
)
plot(p)
dev.off()
