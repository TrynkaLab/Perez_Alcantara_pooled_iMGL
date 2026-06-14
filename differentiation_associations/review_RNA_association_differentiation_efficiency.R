# Checking for any differential RNA to distinguish between lines at the extremes 
# of differentiation (proliferation) efficiency
# This be done at the latest stage (microglia) as we don't have RNA-seq data from premac
# we can also see if there is any association between iPSC gene expression from the Cuomo et al 2020 paper
# and proliferation from microglia to premac
# We can also check if any donor has increased signature of macrophage precursor cells
# Macrophage precursor signature: 
# https://www.frontiersin.org/journals/cellular-neuroscience/articles/10.3389/fncel.2025.1552241/full
# MRC1 (CD206) - from A2 macrophages later present in CAM
# F4/80+ from A2 macrophages, also known as Adgre1 or EMR1. Might be hard to distinguish 
# precursor cells from CAMs without direct comparison to Bian subtypes like I did in
# Washer et al.


library(arrow)
library(tidyverse)
library(limma)
library(edgeR)
library(ggpubr)
library(ggrepel)
# Loading pseudobulk for microglia
microglia_raw_counts= read.table("../../../OTAR2065_sc_eQTL/data/for_tensorQTL/pseudobulk_raw_counts.txt")

microglia_long <- microglia_raw_counts %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = c("treatment", "proliferation", "line"),
    names_pattern = "^(untreated)_(Proliferating|Not_proliferating)_(.+)$",
    values_to = "value"
  )

rm(microglia_raw_counts)
gc()

## loading metadata

donor_metadata = read_csv("../../../OTAR2065_differentiation_efficiency/data/donor_metadata_complete_with_imputed_sex.csv")
metadata_info = read_csv("../../../OTAR2065_differentiation_efficiency/data/metadata_info_hipsci_IPMAR.csv") %>%
  mutate(Disease.Status = case_when(Disease.Status == "Control" ~ "Normal",
                                    Disease.Status == "EOAD" ~ "AD",
                                    Disease.Status == "LOAD" ~ "AD",
                                    .default = "Normal")) 
merged_metadata = donor_metadata %>%
  left_join(metadata_info %>%   dplyr::select(donor,line,Pluritest.pluripotency.score, Pluritest.novelty.score), by = "donor") %>% 
  distinct() %>%
  mutate(line = case_when(donor == "Arlene" ~ "Arlene-003",
                          donor == "Bertha" ~ "Bertha-004",
                          donor == "Cindy" ~ "Cindy-005",
                          donor == "Dexter" ~ "Dexter-006",
                          donor == "Fiona" ~ "Fiona-010",
                          donor == "Gilma" ~ "Gilma-009",
                          donor == "Hector" ~ "Hector-011",
                          donor == "Imani" ~ "Imani-012",
                          donor == "Javier" ~ "Javier-013",
                          donor == "Keoni" ~ "Keoni-014",
                          donor == "Mindy" ~ "Mindy-016",
                          donor == "Nestor" ~ "Nestor-017",
                          donor == "Olaf" ~ "Olaf-018",
                          donor == "Qiana" ~ "Qiana-022",
                          donor == "eorc" ~ "eorc_2",
                          donor == "hayt" ~ "hayt_3",
                          donor == "jaqg" ~ "jaqg",
                          donor == "jotn" ~ "jotn_2",
                          donor == "nurh" ~ "nurh",
                          donor == "pitg" ~ "pitg_2",
                          donor == "qecv" ~ "qecv_2",
                          donor == "uict" ~ "uict_1",
                          .default = line),
         Disease.Status = case_when(Disease.Status == "Control" ~ "Normal",
                                    Disease.Status == "EOAD" ~ "AD",
                                    Disease.Status == "LOAD" ~ "AD",
                                    .default = "Normal"))

#### loading proliferation
ipsc_premac = read_csv("../../data/results/1.2.scale_proliferation/line_prop_changes_premac_iPSC.csv")
premac_microglia = read_csv("../../data/results/1.2.scale_proliferation/line_prop_changes_microglia_premac.csv")
## loading iPSC data from Cuomo et al and format it to pseudobulk from single cell raw counts

if(!file.exists("../../../resources/Cuomo_HipSci_2020/ipsc_pseudobulk_summed_raw_counts.csv.gz")){
ipsc = arrow::read_csv_arrow("../../../resources/Cuomo_HipSci_2020/raw_counts.csv")
ipsc_metadata = read_tsv("../../../resources/Cuomo_HipSci_2020/cell_metadata_cols.tsv") %>%
  dplyr::rename(cellname = compatible_fragment_ratio,
         line = donor_long_id) %>%
  dplyr::select(cellname,line) %>%
  dplyr::mutate(cell_line = paste(cellname,line,sep = "-"))

ipsc_metadata %>%
  group_by(line) %>%
  filter(n() > 1)

# rename lines
colnames(ipsc) <- ipsc_metadata$cell_line[ match(colnames(ipsc), ipsc_metadata$cellname) ]

## aggregating counts from matching lines

lines <- sub(".*-", "", colnames(ipsc))

ipsc_pseudobulk <- sapply(
  split(seq_len(ncol(ipsc)), lines),
  function(i) rowSums(ipsc[, i, drop = FALSE])
)

ipsc_pseudobulk = as.data.frame(ipsc_pseudobulk)
ipsc_pseudobulk$gene = unname(unlist(ipsc[,1]))

ipsc_pseudobulk <- ipsc_pseudobulk %>%
  separate_wider_delim(
    cols = gene,
    delim = "_",
    names = c("gene_id", "gene_name")
  ) %>%
  relocate(c("gene_id", "gene_name"))

write_csv(ipsc_pseudobulk,"../../../resources/Cuomo_HipSci_2020/ipsc_pseudobulk_summed_raw_counts.csv.gz")  
} else{
  ipsc_pseudobulk = read_csv("../../../resources/Cuomo_HipSci_2020/ipsc_pseudobulk_summed_raw_counts.csv.gz")
}
  
#### differential expression for proliferation iPSC -> premac ######
## average where there are different pools

ipsc_premac_average = ipsc_premac %>%
  group_by(line) %>%
  summarise(mean_proliferation = mean(scaled_log_fraction)) %>%
  ungroup() %>%
  filter(line %in% colnames(ipsc_pseudobulk)) # 46 after subsetting to shared lines

ipsc_pseudobulk_subset = ipsc_pseudobulk %>%
  dplyr::select(gene_id, gene_name, ipsc_premac_average$line)
  
# extract count matrix
counts <- ipsc_pseudobulk_subset %>%
  dplyr::select(-gene_id, -gene_name) %>%
  as.matrix()

rownames(counts) <- ipsc_pseudobulk_subset$gene_id

# subset metadata
subset_metadata = merged_metadata %>%
  dplyr::select(line, Sex, contains("Pluritest")) %>%
  dplyr::filter(line %in% colnames(ipsc_pseudobulk_subset)) %>%
  left_join(ipsc_premac_average) %>%
  as.data.frame()
rownames(subset_metadata) = (subset_metadata$line)
subset_metadata = subset_metadata[-1]

counts <- counts[, rownames(subset_metadata)]  # match metadata

all(colnames(counts) == rownames(subset_metadata)) # checking they match

subset_metadata$Sex <- factor(subset_metadata$Sex)

design <- model.matrix(
  ~ mean_proliferation + Sex + Pluritest.pluripotency.score + Pluritest.novelty.score,
  data = subset_metadata
)
 # calculating normalisation factors
y <- DGEList(counts)
y <- calcNormFactors(y)

# calculating mean variance trend and weights for fitting
v <- voom(y, design, plot = TRUE)  

### diff expr
# fit linear model
fit <- lmFit(v, design)
fit <- eBayes(fit)

# results with shrinkage
fit_treat <- treat(fit, lfc = 0.1)
check = topTable(fit_treat, number = Inf,sort.by = "P")
table(check$adj.P.Val< 0.05)

# no significant result - probably why Julie didn't do this?
# doing correlations as in Jerber
correlate_expression_with_covariate <- function(
    expr,
    metadata,
    outcome_col = "mean_proliferation",
    covariates = c(
      "Sex",
      "Pluritest.pluripotency.score",
      "Pluritest.novelty.score"
    ),
    gene_annotation = NULL,
    gene_id_col = "gene_id",
    gene_name_col = "gene_name",
    n_perm = 1000,
    seed = 1
) {
 
  
  # Expression matrix: genes x samples
  expr <- as.matrix(expr)
  
  metadata_vars <- c(outcome_col, covariates)
  
  # Keep samples with complete metadata and expression
  keep <- complete.cases(metadata[, metadata_vars, drop = FALSE]) &
    complete.cases(t(expr))
  
  metadata2 <- metadata[keep, , drop = FALSE]
  expr2 <- expr[, keep, drop = FALSE]
  
  stopifnot(ncol(expr2) == nrow(metadata2))
  
  # Reorder metadata to match expression columns
  metadata2 <- metadata2[colnames(expr2), , drop = FALSE]
  stopifnot(all(colnames(expr2) == rownames(metadata2)))
  
  # Build covariate formula
  covar_formula <- as.formula(
    paste("~", paste(covariates, collapse = " + "))
  )
  
  # Design matrix for adjustment covariates
  C <- model.matrix(covar_formula, data = metadata2)
  
  # Residual-maker matrix
  M <- diag(nrow(C)) - C %*% solve(t(C) %*% C) %*% t(C)
  
  # Residualise outcome
  x <- metadata2[[outcome_col]]
  x_resid <- as.numeric(M %*% x)
  x_resid <- as.numeric(scale(x_resid, center = TRUE, scale = FALSE))
  
  # Residualise expression across samples
  expr_resid <- expr2 %*% M
  rownames(expr_resid) <- rownames(expr2)
  colnames(expr_resid) <- colnames(expr2)
  
  expr_resid_centered <- expr_resid - rowMeans(expr_resid)
  
  # Pearson correlation and slopes
  numerator <- as.numeric(expr_resid_centered %*% x_resid)
  
  expr_ss <- rowSums(expr_resid_centered^2)
  x_ss <- sum(x_resid^2)
  
  cor_vals <- numerator / sqrt(expr_ss * x_ss)
  lm_slopes <- numerator / x_ss
  
  results <- tibble(
    gene_id = rownames(expr_resid),
    cor = cor_vals,
    slope = lm_slopes
  )
  
  # Optional gene annotation
  if (!is.null(gene_annotation)) {
    results <- results %>%
      left_join(
        gene_annotation %>%
          dplyr::select(
            gene_id = all_of(gene_id_col),
            gene_name = all_of(gene_name_col)
          ) %>%
          distinct(),
        by = "gene_id"
      )
  }
  
  # Permutation test
  set.seed(seed)
  
  n_genes <- nrow(expr_resid)
  perm_slopes <- matrix(NA_real_, nrow = n_genes, ncol = n_perm)
  
  for (i in seq_len(n_perm)) {
    x_perm <- sample(x_resid)
    perm_numerator <- as.numeric(expr_resid_centered %*% x_perm)
    perm_slopes[, i] <- perm_numerator / sum(x_perm^2)
  }
  
  obs_slopes <- results$slope[match(rownames(expr_resid), results$gene_id)]
  
  p_vals <- (rowSums(abs(perm_slopes) >= abs(obs_slopes)) + 1) / (n_perm + 1)
  
  results <- results %>%
    mutate(
      p_perm = p_vals[match(gene_id, rownames(expr_resid))],
      FDR = p.adjust(p_perm, method = "BH")
    ) %>%
    arrange(p_perm)
  
  return(results)
}

results <- correlate_expression_with_covariate(
  expr = v$E,
  metadata = subset_metadata,
  outcome_col = "mean_proliferation",
  covariates = c(
    "Sex",
    "Pluritest.pluripotency.score",
    "Pluritest.novelty.score"
  ),
  gene_annotation = ipsc_pseudobulk_subset,
  gene_id_col = "gene_id",
  gene_name_col = "gene_name",
  n_perm = 1000,
  seed = 1
)

table(results$p_perm < 0.05) # 921

# selecting positive and negatively correlated vectors

positive_genes <- results %>%
  filter(cor > 0, p_perm < 0.05) %>%
  arrange(desc(cor)) %>%
  pull(gene_name)
# https://biit.cs.ut.ee/gplink/l/a_Q3dxI-vST
negative_genes <- results %>%
  filter(cor < 0, p_perm < 0.05) %>%
  arrange(cor) %>%
  pull(gene_name)
# These negative genes are enriched in E2F3 TF binding sites 
# And other E2F- family TFs
# the genes in this list are significantly enriched for targets that are 
# normally activated by the E2F3 transcription factor. 
#  E2F3 is a key activator of genes needed for cell cycle progression (G1/S transition),
# so this finding suggests a reduction in cell proliferation or a halt in the cell cycle.
# https://pmc.ncbi.nlm.nih.gov/articles/PMC316459/
# gprofiler version e114_eg62_p19_27110d83 
# cite https://biit.cs.ut.ee/gprofiler/page/citing
# Bonferroni adjusted p-value 2x10-12
# https://biit.cs.ut.ee/gplink/l/a7vBpP9sKRs
dir.create("../../data/results/review_RNA_association_differentiation_efficiency")
write_csv(results %>% arrange(desc(cor)),"../../data/results/review_RNA_association_differentiation_efficiency/correlations_iPSC_expr_with_iPSC_to_premac_proliferation.csv")


## Histogram


plot_cor_histogram <- function(
    df,
    genes_to_label,
    cor_col = cor,
    p_col = p_perm,
    gene_col = gene_name,
    p_cutoff = 0.05,
    bins = 50
) {
  
  cor_col  <- enquo(cor_col)
  p_col    <- enquo(p_col)
  gene_col <- enquo(gene_col)
  
  # Significant cutoffs closest to zero on each side
  sig_cutoffs <- df %>%
    filter(!!p_col < p_cutoff) %>%
    summarise(
      min_positive_cor = min((!!cor_col)[!!cor_col > 0], na.rm = TRUE),
      max_negative_cor = max((!!cor_col)[!!cor_col < 0], na.rm = TRUE)
    )
  
  neg_cutoff <- sig_cutoffs$max_negative_cor
  pos_cutoff <- sig_cutoffs$min_positive_cor
  
  # Pre-compute histogram bins so we can color bars by position
  hist_df <- df %>%
    filter(!is.na(!!cor_col)) %>%
    mutate(cor_value = !!cor_col) %>%
    ggplot(aes(x = cor_value)) +
    geom_histogram(bins = bins)
  
  hist_df <- ggplot_build(hist_df)$data[[1]] %>%
    as_tibble() %>%
    mutate(
      bin_mid = (xmin + xmax) / 2,
      fill_group = if_else(
        bin_mid > neg_cutoff & bin_mid < pos_cutoff,
        "between_cutoffs",
        "outside_cutoffs"
      )
    )
  
  # Genes to label
  label_df <- df %>%
    filter(!!gene_col %in% genes_to_label) %>%
    mutate(
      label = as.character(!!gene_col),
      y_pos = 0
    )
  
  ggplot(hist_df, aes(x = x, y = count, fill = fill_group)) +
    geom_col(
      aes(width = xmax - xmin),
      color = "white"
    ) +
    scale_fill_manual(
      values = c(
        between_cutoffs = "grey90",
        outside_cutoffs = "grey70"
      ),
      guide = "none"
    ) +
    geom_vline(
      data = sig_cutoffs,
      aes(xintercept = min_positive_cor),
      linetype = "dotted",
      linewidth = 1,
      inherit.aes = FALSE
    ) +
    geom_vline(
      data = sig_cutoffs,
      aes(xintercept = max_negative_cor),
      linetype = "dotted",
      linewidth = 1,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = label_df,
      aes(x = !!cor_col, y = y_pos),
      inherit.aes = FALSE,
      size = 2
    ) +
    geom_label_repel(
      data = label_df,
      aes(x = !!cor_col, y = y_pos, label = label),
      inherit.aes = FALSE,
      nudge_y = 5,
      min.segment.length = 0,
      box.padding = 0.4,
      point.padding = 0.3,
      max.overlaps = Inf
    ) +
    labs(
      x = "Correlation",
      y = "Number of genes",
      title = "Distribution of correlations",
      subtitle = paste0(
        "Dotted lines show permuted p-value < ",
        p_cutoff
      )
    ) +
    theme_pubr()
}

p = plot_cor_histogram(results,genes_to_label = c("LRRC55",
                                              "BCOR", "MLLT4"))  
pdf("../../data/results/review_RNA_association_differentiation_efficiency/correlations_iPSC_expr_with_iPSC_to_premac_proliferation.pdf",
    height = 4,width = 5)
plot(p)
dev.off()

png("../../data/results/review_RNA_association_differentiation_efficiency/correlations_iPSC_expr_with_iPSC_to_premac_proliferation.png",
    height = 4,width = 5, res = 600,units = "in")
plot(p)
dev.off()



# lines that are in iPSC but not in premac 
ipsc_pseudobulk_alllines = read_csv("../../../resources/Cuomo_HipSci_2020/ipsc_pseudobulk_summed_raw_counts.csv.gz")
ipsc_pseudobulk_alllines = ipsc_pseudobulk_alllines %>%
  column_to_rownames("gene_id") %>%
  select(-gene_name)

plated = read.table("../../data/results/1.alluvial_plots/pools2-11_13-17_changing_props_iPSC_preMacs_microglia_WGS_sc.txt") %>%
  dplyr::select(V1) %>%
  unlist() %>%
  unique()

# which ones disappear in premac already
plated_dropped_premac = plated[!plated %in% unique(ipsc_premac$line)] # 21
plated_dropped_microglia = plated[!plated %in% unique(premac_microglia$line)] #24

table(unique(plated_dropped_premac) %in% colnames(ipsc_pseudobulk_alllines)) # 7 dropped we have iPSC expression
table(unique(plated_dropped_microglia) %in% colnames(ipsc_pseudobulk_alllines)) # same 7
# Not enough to do anything


## Now for microglia

head(microglia_long)
microglia_proliferating = microglia_long %>%
  filter(proliferation == "Proliferating" & treatment == "untreated")
microglia_not_proliferating = microglia_long %>%
  filter(proliferation == "Not_proliferating" & treatment == "untreated")


counts_prolif <- microglia_proliferating %>%
  dplyr::select(-treatment, -proliferation) %>%
  pivot_wider(id_cols = "gene", names_from = "line",values_from = "value",
              values_fill = NA) %>%
  as.data.frame() %>%
  column_to_rownames("gene") %>%
  as.matrix()

counts_not_prolif <- microglia_not_proliferating %>%
  dplyr::select(-treatment, -proliferation) %>%
  pivot_wider(id_cols = "gene", names_from = "line",values_from = "value",
              values_fill = NA) %>%
  as.data.frame() %>%
  column_to_rownames("gene") %>%
  as.matrix()

head(premac_microglia)
# 
# premac_microglia_average =  premac_microglia %>%
#   dplyr::filter(treatment == "untreated") %>%
#   group_by(line) %>%
#   summarise(mean_proliferation = mean(scaled_log_fraction)) %>%
#   ungroup() %>%
#   filter(line %in% unique(colnames(counts_not_prolif),
#                           colnames(counts_prolif))) # 164 after subsetting to shared lines
# 
# subset metadata
subset_metadata_prolif = merged_metadata %>%
  dplyr::select(line, Sex, contains("Pluritest")) %>%
  dplyr::filter(line %in% colnames(counts_prolif)) %>%
  left_join(premac_microglia_average) %>%
  as.data.frame() # 65 lines

subset_metadata_not_prolif = merged_metadata %>%
  dplyr::select(line, Sex, contains("Pluritest")) %>%
  dplyr::filter(line %in% colnames(counts_not_prolif)) %>%
  left_join(premac_microglia_average) %>%
  as.data.frame() # 171 lines

# Lines in prolif not in not prolif
table(colnames(counts_prolif) %in% colnames(counts_not_prolif)) # all present in not prolif
# Lines in "not prolif" that are not in "prolif"
table(colnames(counts_not_prolif) %in% colnames(counts_prolif)) # 112 lines

# that are in not prolif are missing in prolif
# Those lines are the ones that are mostly present in the better differentiated microglia cluster
# Not dominated by the proliferative signature

rownames(subset_metadata_prolif) = (subset_metadata_prolif$line)
subset_metadata_prolif = subset_metadata_prolif[-1]

rownames(subset_metadata_not_prolif) = (subset_metadata_not_prolif$line)
subset_metadata_not_prolif = subset_metadata_not_prolif[-1]


run_voom_correlation <- function(
    counts,
    metadata,
    outcome_col,
    covar_formula,
    n_perm = 1000,
    seed = 1,
    voom_plot = TRUE
) {

  outcome_col <- enquo(outcome_col)

  # Make sure sample order matches
  counts <- counts[, rownames(metadata), drop = FALSE]
  stopifnot(all(colnames(counts) == rownames(metadata)))

  # Keep complete samples
  vars_needed <- all.vars(covar_formula)
  outcome_name <- quo_name(outcome_col)
  metadata_vars <- c(outcome_name, vars_needed)

  keep_samples <- complete.cases(metadata[, metadata_vars]) &
    complete.cases(t(counts))

  metadata2 <- metadata[keep_samples, , drop = FALSE]
  counts2 <- counts[, rownames(metadata2), drop = FALSE]

  stopifnot(all(colnames(counts2) == rownames(metadata2)))

  # limma/voom transformation
  design_full <- model.matrix(
    as.formula(
      paste0("~ ", outcome_name, " + ", paste(vars_needed, collapse = " + "))
    ),
    data = metadata2
  )

  y <- DGEList(counts2)
  y <- calcNormFactors(y)
  v <- voom(y, design_full, plot = voom_plot)

  expr <- as.matrix(v$E)

  # Residualize expression and outcome against covariates
  C <- model.matrix(covar_formula, data = metadata2)
  M <- diag(nrow(C)) - C %*% solve(t(C) %*% C) %*% t(C)

  x <- metadata2[[outcome_name]]
  x_resid <- as.numeric(M %*% x)
  x_resid <- as.numeric(scale(x_resid, center = TRUE, scale = FALSE))

  expr_resid <- expr %*% M
  rownames(expr_resid) <- rownames(expr)
  colnames(expr_resid) <- colnames(expr)

  expr_resid_centered <- expr_resid - rowMeans(expr_resid)

  numerator <- as.numeric(expr_resid_centered %*% x_resid)
  expr_ss <- rowSums(expr_resid_centered^2)
  x_ss <- sum(x_resid^2)

  cor_vals <- numerator / sqrt(expr_ss * x_ss)
  slopes <- numerator / x_ss

  results <- tibble(
    gene_id = rownames(expr_resid),
    cor = cor_vals,
    slope = slopes
  )

  # Permutation p-values
  set.seed(seed)

  n_genes <- nrow(expr_resid)
  perm_slopes <- matrix(NA_real_, nrow = n_genes, ncol = n_perm)

  for (i in seq_len(n_perm)) {
    x_perm <- sample(x_resid)
    perm_numerator <- as.numeric(expr_resid_centered %*% x_perm)
    perm_slopes[, i] <- perm_numerator / sum(x_perm^2)
  }

  obs_slopes <- slopes

  p_perm <- (rowSums(abs(perm_slopes) >= abs(obs_slopes)) + 1) / (n_perm + 1)

  results <- results %>%
    mutate(
      p_perm = p_perm,
      FDR = p.adjust(p_perm, method = "BH")
    ) %>%
    arrange(p_perm)

  list(
    results = results,
    voom = v,
    metadata = metadata2,
    expr = expr,
    expr_resid = expr_resid,
    x_resid = x_resid
  )
}

premac_microglia_analysis_prolif = run_voom_correlation(counts = counts_prolif,
                                                        metadata = subset_metadata_prolif,
                                                        outcome_col = mean_proliferation,
                                                        covar_formula = ~ Sex + Pluritest.pluripotency.score + Pluritest.novelty.score,
                                                        n_perm = 1000,
                                                        seed = 1,
                                                        voom_plot = TRUE)

premac_microglia_analysis_not_prolif = run_voom_correlation(counts = counts_not_prolif,
                                                        metadata = subset_metadata_not_prolif,
                                                        outcome_col = mean_proliferation,
                                                        covar_formula = ~ Sex + Pluritest.pluripotency.score + Pluritest.novelty.score,
                                                        n_perm = 1000,
                                                        seed = 1,
                                                        voom_plot = TRUE)

premac_microglia_results_prolif <- premac_microglia_analysis_prolif$results

premac_microglia_results_not_prolif <- premac_microglia_analysis_not_prolif$results


table(premac_microglia_results_prolif$p_perm < 0.05)
table(premac_microglia_results_not_prolif$p_perm < 0.05)

negative_deleterious_genes_premac_microglia = c("GSPT1","CFHR3","ZAR1","OVGP1")
positive_deleterious_genes_premac_microglia = c("KIAA1239","BCAR1","SACS","PDGFRB",
                                                "IL37","INPP5D","ATP11A","CASP5",
                                                "AP4E1","PGM3","DOCK9","CIDEA",
                                                "CACNA1F","KANSL2")

premac_microglia_results_not_prolif %>%
  dplyr::filter(gene_id %in% negative_deleterious_genes_premac_microglia)

premac_microglia_results_not_prolif %>%
  dplyr::filter(gene_id %in% positive_deleterious_genes_premac_microglia)


write_csv(
  premac_microglia_results_prolif %>% arrange(desc(cor)),
  "../../data/results/review_RNA_association_differentiation_efficiency/correlations_microglia_expr_with_premac_to_microglia_proliferation_prolif_cluster.csv"
)

write_csv(
  premac_microglia_results_not_prolif %>% arrange(desc(cor)),
  "../../data/results/review_RNA_association_differentiation_efficiency/correlations_microglia_expr_with_premac_to_microglia_proliferation_not_prolif_cluster.csv"
)

p1 = plot_cor_histogram(premac_microglia_results_not_prolif %>% rename(gene_name = gene_id),genes_to_label = c(negative_deleterious_genes_premac_microglia,
                                                  positive_deleterious_genes_premac_microglia))  

p2 = plot_cor_histogram(premac_microglia_results_prolif %>% rename(gene_name = gene_id),genes_to_label = c(negative_deleterious_genes_premac_microglia,
                                                                                                              positive_deleterious_genes_premac_microglia))  


png("../../data/results/review_RNA_association_differentiation_efficiency/correlations_microglia_expr_with_premac_to_iMGL_proliferation.png",
    height = 4,width = 5, res = 600,units = "in")
plot(p1)
dev.off()

### Association with burden of deleterious variants

load("../../../OTAR2065_Daianna/output_data/04_Burden_test_proliferation/01_Burden_tests/all_Del_Burdens_per_line.Rdata")

mutBurden_del = readRDS("../../../resources/Puigdevall_Neuroseq_efficiency_2023/mutBurden_del.RDS")

ipsc_pseudobulk_bcor = ipsc_pseudobulk_subset %>%
  filter(gene_name=="BCOR") %>%
  select(-gene_id) %>%
  pivot_longer(cols = -gene_name, names_to = "line", values_to = "BCOR_expr_value") %>%
  select(-gene_name)

deleterious_bcor_vs_iPSC = mutBurden_del[c("BCOR", "MLLT4"),] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("line") %>%
  dplyr::filter(
    stringr::str_detect(
      line,
      stringr::str_c(colnames(ipsc_pseudobulk_subset), collapse = "|")
    )
  ) %>%
  mutate(line = stringr::str_remove(line, pattern = "^.*?-"),
         BCOR = if_else(BCOR > 1,1,BCOR)  ## aggregate BCOR del mutation in presence or absence
         ) %>%
  left_join(ipsc_pseudobulk_bcor) 

subset_metadata
## checking correlation divided by bcor status
# no del mutations
no_del_bcor_lines = deleterious_bcor_vs_iPSC %>% filter(BCOR == 0) %>% pull(line)
matched_lines <- intersect(no_del_bcor_lines, colnames(v$E))
matched_lines <- intersect(matched_lines, rownames(subset_metadata))

results_no_del_bcor_lines <- correlate_expression_with_covariate(
  expr = v$E[, matched_lines, drop = FALSE],
  metadata = subset_metadata[matched_lines, , drop = FALSE],
  outcome_col = "mean_proliferation",
  covariates = c(
    "Sex",
    "Pluritest.pluripotency.score",
    "Pluritest.novelty.score"
  ),
  gene_annotation = ipsc_pseudobulk_subset %>%
    dplyr::select(gene_id, gene_name) %>%
    dplyr::distinct(),
  gene_id_col = "gene_id",
  gene_name_col = "gene_name",
  n_perm = 1000,
  seed = 1
)

with_del_bcor_lines = deleterious_bcor_vs_iPSC %>% filter(BCOR == 1) %>% pull(line)
matched_lines <- intersect(with_del_bcor_lines, colnames(v$E))
matched_lines <- intersect(matched_lines, rownames(subset_metadata))

results_with_del_bcor_lines <- correlate_expression_with_covariate(
  expr = v$E[, matched_lines, drop = FALSE],
  metadata = subset_metadata[matched_lines, , drop = FALSE],
  outcome_col = "mean_proliferation",
  covariates = c(
    "Sex",
    "Pluritest.pluripotency.score",
    "Pluritest.novelty.score"
  ),
  gene_annotation = ipsc_pseudobulk_subset %>%
    dplyr::select(gene_id, gene_name) %>%
    dplyr::distinct(),
  gene_id_col = "gene_id",
  gene_name_col = "gene_name",
  n_perm = 1000,
  seed = 1
)
# no difference

wilcox_bcor_vs_iPSC <- deleterious_bcor_vs_iPSC %>%
  dplyr::filter(!is.na(BCOR), !is.na(BCOR_value)) %>%
  dplyr::mutate(BCOR = factor(BCOR, levels = c(0, 1))) %>%
  wilcox.test(BCOR_value ~ BCOR, data = .)
# not significant, despite correlation between gene expression and diff eff being significant

# subset to genes with significant burden
# only BCOR, as MLLT4 is not significant and for LRR55 we have no expression
# caveat: we don't know at what end of the differentiation the gene expression deficiency
# manifests (e.g. if LRRC55 is associated with higher premac differentiation efficiency
# but is expressed only in premac)

##### explore pathways of BCOS del present vs absent and see if there are obvious 
#### TFs, GO patthways and so on that come up
## iPSC RNA-seq

