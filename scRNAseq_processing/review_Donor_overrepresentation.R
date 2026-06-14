# see if donors are over-represented in any cluster

library(Seurat)
library(future)
library(tidyverse)
library(lme4)
library(ggrepel)
library(ggpubr)

outdir = "../../data/results/1.6.donor_overrepresentation"
dir.create(file.path(outdir))

source("./helpers.R")
# set this option when analyzing large datasets
options(future.globals.maxSize = 110000 * 1024^2) # 110Gb
treatment_cols =  c(untreated = "#8D918B", IFN = "#3A5683", LPS = "#F8766D")

# for new seurat assay slots
options(Seurat.object.assay.version = "v5")

directory = "../../data/results/1.4.markers"

# read in annotated seurat objects
seu_subset = list()
for(treatment in c("untreated", "IFN", "LPS")) {
  
  seu_subset[[treatment]] = LoadSeuratRds(
    file = paste0(directory,"/",treatment,"_filtered_harmony/",treatment,"_filtered_harmony.Rds"))
  
  
}


##### test for every donor if their expected proportion is skewed in one cluster vs others
results_global = list()

for(t in names(seu_subset)){
  

meta <- seu_subset[[t]]@meta.data

# Total cells per donor per cluster
donor_cluster_counts <- meta %>%
  count(donor_id, merged_clusters) %>%
  rename(n_cluster = n)

# Total cells per donor (global)
donor_totals <- meta %>%
  count(donor_id) %>%
  rename(total_cells = n)

# Total cells per cluster (global)
cluster_totals <- meta %>%
  count(merged_clusters) %>%
  rename(cluster_total = n)

# All clusters
all_clusters <- unique(meta$merged_clusters)

# We want for a given donor:
#   
#   Observed: vector of counts across all clusters → n_cluster
# 
# Expected: same total cells, distributed across clusters according to global cluster sizes, weighted by total cells of donor
# 
# Only test if a donor is over-represented in a specific cluster, not under-represented.

test_donor_over_global <- function(donor, min_cells = 10) {
  
  # donor's counts across all clusters
  donor_counts <- donor_cluster_counts %>%
    filter(donor_id == donor) %>%
    complete(merged_clusters = all_clusters, fill = list(n_cluster = 0)) %>%
    arrange(merged_clusters)
  
  # skip donors with too few total cells
  total_cells <- sum(donor_counts$n_cluster)
  if(total_cells < min_cells) return(NULL)
  
  # global cluster proportions (weights)
  global_props <- cluster_totals %>%
    arrange(merged_clusters) %>%
    mutate(prop = cluster_total / sum(cluster_total)) %>%
    pull(prop)
  
  # observed for this donor
  observed <- donor_counts$n_cluster
  
  # expected counts if donor distributed globally proportional to cluster sizes
  expected <- total_cells * global_props
  
  # Chi-squared test (goodness of fit)
  test <- chisq.test(x = observed, p = global_props, rescale.p = TRUE, simulate.p.value = TRUE)
  
  # Only keep clusters where observed > expected (over-representation)
  donor_counts <- donor_counts %>%
    mutate(expected = expected,
           donor_prop = n_cluster / total_cells,
           over_rep = donor_prop > expected / total_cells)
  
  tibble(
    donor_id = donor,
    merged_clusters = donor_counts$merged_clusters,
    n_cluster = donor_counts$n_cluster,
    donor_prop = donor_counts$donor_prop,
    expected_prop = donor_counts$expected / total_cells,
    p.value = test$p.value
  )
}
# rescale.p = TRUE ensures probabilities sum to 1.
# 
# simulate.p.value = TRUE avoids issues with small counts.

all_donors <- unique(meta$donor_id)

results_global[treatment] <- map_dfr(all_donors, ~test_donor_over_global(.x))

# FDR correction
results_global[treatment] <- results_global[treatment] %>%
  mutate(p.adj = p.adjust(p.value, method = "BH"),
         treatment = treatment) %>%
  filter(donor_prop > expected_prop)  # keep only over-represented clusters
}

results_enrichment = do.call("rbind", results_global)

info = read_csv("../../data/donor_metadata_complete_with_imputed_sex.csv")

results_enrichment = results_enrichment %>%
  mutate(donor = str_replace(donor_id, "[-_].*$", "")) %>%
  left_join(info, by = "donor") %>%
  mutate("Disease status" = case_when(Disease.Status == "Control" ~ "Normal",
                                      Disease.Status ==  "EOAD" ~ "Early-onset AD",
                                      Disease.Status ==  "LOAD" ~ "Late-onset AD",
                                      .default = "Normal"))

# ------------------------------
# Step 1: select top 5 donors per cluster
# ------------------------------

top5_donors <- results_enrichment %>%
  group_by(merged_clusters,treatment) %>%
  slice_max(order_by = donor_prop, n = 5) %>%
  ungroup()

# ------------------------------

# Set a jitter object
jitter_pos <- position_jitter(width = 0.1)

p <- ggplot(results_enrichment, aes(x = merged_clusters, y = donor_prop)) +
  geom_jitter(aes(color = `Disease status`), width = 0.1, alpha = 0.6, size = 1.5) +
  # geom_text_repel with same jitter
  geom_text_repel(
    data = top5_donors,
    aes(label = donor_id, y = donor_prop, color = `Disease status`),
    size = 3,
    max.overlaps = Inf,
    segment.color = "black"
  ) +
  scale_color_manual(values = c("Early-onset AD" = "#922D50", 
                                "Late-onset AD" = "#F87575",
                                "Normal" = "#1f78b4")) +
  theme_minimal() +
  labs(
    x = "Cluster",
    y = "Fraction of cluster contributed by donor",
    title = "Donor proportion per cluster (top 5 donors labeled)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  facet_wrap(vars(treatment),scales = "free" )

p





##### other way of looking at it:   ########
# Does any donor contribute more cells to this cluster than all other donors?
results_enrichment = list() 
for(t in names(seu_subset)){
  meta <- seu_subset[[t]]@meta.data
  
# ------------------------------
# Step 1: donor × cluster counts
# ------------------------------
donor_cluster_counts <- meta %>%
  count(donor_id, merged_clusters) %>%
  rename(n_cluster = n)

# ------------------------------
# Step 2: compute donor_prop per cluster
# donor_prop = fraction of cluster contributed by this donor
# ------------------------------
donor_props <- donor_cluster_counts %>%
  group_by(merged_clusters) %>%
  mutate(donor_prop = n_cluster / sum(n_cluster)) %>%
  ungroup() %>%
  filter(n_cluster >= 10)  # filter out donors with <10 cells in that cluster

# ------------------------------
# Step 3: Wilcoxon test function per donor
# ------------------------------
test_donor_enrichment <- function(cluster_df, donor) {
  donor_prop_val <- cluster_df %>%
    filter(donor_id == donor) %>%
    pull(donor_prop)
  
  other_props <- cluster_df %>%
    filter(donor_id != donor) %>%
    pull(donor_prop)
  
  if(length(donor_prop_val) == 0 | length(other_props) < 2) return(NULL)
  
  test <- wilcox.test(donor_prop_val, other_props, alternative = "greater")
  
  tibble(
    donor_id = donor,
    merged_clusters = unique(cluster_df$merged_clusters),
    donor_prop = donor_prop_val,
    median_other_prop = median(other_props),
    p.value = test$p.value
  )
}

# ------------------------------
# Step 4: run test for all donors × clusters
# ------------------------------
all_clusters <- unique(donor_props$merged_clusters)

results_enrichment[[t]] <- map_dfr(all_clusters, function(clust) {
  cluster_df <- donor_props %>% filter(merged_clusters == clust)
  map_dfr(unique(cluster_df$donor_id), ~test_donor_enrichment(cluster_df, .x))
})

# FDR correction
results_enrichment[[t]] <- results_enrichment[[t]] %>%
  mutate(p.adj_BH = p.adjust(p.value, method = "BH"),
         p.adj_bonferroni= p.adjust(p.value, method = "bonferroni")) %>%
  arrange(p.adj_BH)

results_enrichment[[t]] = results_enrichment[[t]]  %>%
  mutate( treatment = t)


}

results_enrichment = do.call("rbind",results_enrichment)
info = read_csv("../../data/donor_metadata_complete_with_imputed_sex.csv")

results_enrichment = results_enrichment %>%
  mutate(donor = str_replace(donor_id, "[-_].*$", "")) %>%
  left_join(info, by = "donor") %>%
  mutate("Disease status" = case_when(Disease.Status == "Control" ~ "Normal",
                                      Disease.Status ==  "EOAD" ~ "Early-onset AD",
                                      Disease.Status ==  "LOAD" ~ "Late-onset AD",
                                      .default = "Normal"))

# ------------------------------
# Step 1: select top 5 donors per cluster
# ------------------------------

top5_donors <- results_enrichment %>%
  group_by(merged_clusters,treatment) %>%
  slice_max(order_by = donor_prop, n = 5) %>%
  ungroup()

# ------------------------------

# Set a jitter object
jitter_pos <- position_jitter(width = 0.1)

p <- ggplot(results_enrichment, aes(x = merged_clusters, y = donor_prop)) +
  geom_jitter(aes(color = `Disease status`), width = 0.1, alpha = 0.6, size = 1.5) +
  # geom_text_repel with same jitter
  geom_text_repel(
    data = top5_donors,
    aes(label = donor_id, y = donor_prop, color = `Disease status`),
    size = 3,
    max.overlaps = Inf,
    segment.color = "black"
  ) +
  scale_color_manual(values = c("Early-onset AD" = "#922D50", 
                                "Late-onset AD" = "#F87575",
                                "Normal" = "#1f78b4")) +
  theme_minimal() +
  labs(
    x = "Cluster",
    y = "Fraction of cluster contributed by donor",
    title = "Donor proportion per cluster (top 5 donors labeled)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  facet_wrap(vars(treatment),scales = "free" )

p

png("../../data/results/1.6.donor_overrepresentation/donor_prop_per_cluster.png",
    width = 10,height = 5,res = 400,units = "in")
plot(p)
dev.off()

write_csv(results_enrichment,
          "../../data/results/1.6.donor_overrepresentation/enrichment_results.csv")

min(results_enrichment$p.adj_BH) #0.74
min(results_enrichment$p.adj_bonferroni) #1
