
library(tidyverse)
library(dplyr)
library(SummarizedExperiment)
library(rlang)
library(limma)
library(edgeR)
library(ComplexHeatmap)
library(cowplot)
library(sessioninfo)

#-------------------------------------------------------------------------------
#             5.1 DGE analysis for line phagocytosis phenotypes
#-------------------------------------------------------------------------------
#  Code to run within-treatment DGE for cell line phagocytosis estimates.
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

## Set working dir
setwd("/lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_Daianna/")

## Define output data and plot dir 
outdir = paste(getwd(), "output_data", "05_Diff_expr_phagocytosis_permutations", "01_run_phagocytosis_DGE", sep = "/")
plotdir = paste(getwd(), "plots", "05_Diff_expr_phagocytosis_permutations", "01_run_phagocytosis_DGE", sep = "/")
dir.create(outdir, recursive = T)
dir.create(plotdir, recursive = T)

## Input dir
input_dir00 = paste(getwd(), "output_data", "01_Diff_expr_PRS_permutations", "00_Data_exploration_processing", sep = "/")
input_dir01 = paste(getwd(), "output_data", "05_Diff_expr_phagocytosis_permutations", "00_Data_exploration_processing", sep = "/")

## color and shape dict
load(paste0(input_dir00, "/var_colors.Rdata"), verbose = T)
load(paste0(input_dir00, "/level_colors.Rdata"), verbose = T)
load(paste0(input_dir00, "/level_shapes.Rdata"), verbose = T)
## Add new ones
var_colors["phagocytosis_mean_log_fraction"] <- "darksalmon"
var_colors["phagocytosis_mean_min_prop"] <-  "darkseagreen"
var_colors[c("genotype_PC1", "genotype_PC2", 
             "genotype_PC3", "genotype_PC4", "genotype_PC5")] <-  c("mediumpurple2", "mediumorchid1", "mediumvioletred",
                                                                    "lightslateblue", "mediumpurple4")

save(var_colors, file = paste0(paste("output_data", "01_Diff_expr_PRS_permutations", "00_Data_exploration_processing", sep = "/"), "/var_colors.Rdata"))


## load rse with phagocytosis data
load(paste0(input_dir01, "/rse_non_prolif_filt_phago_IFN.Rdata"), verbose = T)
load(paste0(input_dir01, "/rse_non_prolif_filt_phago_LPS.Rdata"), verbose = T)
load(paste0(input_dir01, "/rse_non_prolif_filt_phago_untreated.Rdata"), verbose = T)


## Run DGE within each treatment 
run_phagocytosis_DEG <- function(treatment){
  
  rse <- eval(parse_expr(paste("rse_non_prolif_filt_phago", treatment, sep = "_")))
  
  message(paste(dim(rse)[2], treatment, "samples across", length(unique(rse$line)), "unique lines of", length(unique(rse$donor)), "donors"))
  # IFN: 180 IFN samples across 139 unique lines of 138 donors
  # LPS: 182 LPS samples across 146 unique lines of 145 donors
  # untreated: 175 untreated samples across 140 unique lines of 139 donors
  
  ## Filter lowly-expressed genes (>0.1 CPM in at least 30% samples)
  expressed_genes <- rowSums(edgeR::cpm(rse) > 0.1) >= floor(0.3*ncol(rse))
  message(paste0(sum(expressed_genes), " genes out of ", nrow(rse), " with >0.1 CPM in at least 30% of samples"))
  # IFN: 12823 genes out of 18648 with >0.1 CPM in at least 30% of samples
  # LPS: 12840 genes out of 18648 with >0.1 CPM in at least 30% of samples
  # Untreated: 13170 genes out of 18648 with >0.1 CPM in at least 30% of samples
  
  rse_expr = rse[expressed_genes, ]  
  
  ## 0 zero-expressed genes
  table(rowSums(assay(rse_expr))==0)
  
  ## TMM norm factors for voom CPM
  rse_expr_norm <- calcNormFactors(rse_expr, method = "TMM")
  
  formula <- ~ phagocytosis_mean_scaled_log_fraction + Sex + phagocytosis_mean_min_prop + genotype_PC1 + genotype_PC2
  model = model.matrix(formula, data = rse_expr_norm$samples) %>%  as.data.frame()
  
  v = voom(rse_expr_norm, design = model, plot = TRUE)
  
  ## Intra-pool corr 
  cor = duplicateCorrelation(v, design = model, block = rse_expr$pool)
  
  ## Re-compute voom weights 
  v2 = voom(rse_expr_norm, design = model, plot=TRUE, block = rse_expr$pool, correlation = cor$consensus)
  
  assays(rse_expr)$logcounts <- v2$E
  
  ## Corr based on corrected weights
  cor2 = duplicateCorrelation(v2, design = model, block = rse_expr$pool)
  
  ## Fit linear model
  fit = lmFit(v2, design = model, block = rse_expr$pool, correlation = cor2$consensus)
  eBGene = eBayes(fit)
  
  top_genes = topTable(eBGene, coef = "phagocytosis_mean_scaled_log_fraction", p.value = 1, number = nrow(rse_expr), sort.by="none")
  
  ## DEGs?
  de_genes <- subset(top_genes, adj.P.Val<0.05)
  message("DGEs for phagocytosis:")
  print(paste(dim(de_genes)[1], "DEGs:", 
              dim(subset(de_genes, logFC>0))[1], "up-regulated and", 
              dim(subset(de_genes, logFC<0))[1], "down-regulated"))
  
  return(list(rse_expr, top_genes, de_genes))
}

############
## IFN DGE
############
phago_DGE_results_IFN <- run_phagocytosis_DEG("IFN")
# "4914 DEGs: 2420 up-regulated and 2494 down-regulated"

rse_non_prolif_filt_expr_phago_IFN <- phago_DGE_results_IFN[[1]]
save(rse_non_prolif_filt_expr_phago_IFN, file = paste0(outdir, "/rse_non_prolif_filt_expr_phago_IFN.Rdata"))
top_genes_phago_IFN <- phago_DGE_results_IFN[[2]]
save(top_genes_phago_IFN, file = paste0(outdir, "/top_genes_phago_IFN.Rdata"))
de_genes_phago_IFN <- phago_DGE_results_IFN[[3]]
save(de_genes_phago_IFN, file = paste0(outdir, "/de_genes_phago_IFN.Rdata"))

############
## LPS DGE
############
phago_DGE_results_LPS <- run_phagocytosis_DEG("LPS")
# "393 DEGs: 178 up-regulated and 215 down-regulated"

rse_non_prolif_filt_expr_phago_LPS <- phago_DGE_results_LPS[[1]]
save(rse_non_prolif_filt_expr_phago_LPS, file = paste0(outdir, "/rse_non_prolif_filt_expr_phago_LPS.Rdata"))
top_genes_phago_LPS <- phago_DGE_results_LPS[[2]]
save(top_genes_phago_LPS, file = paste0(outdir, "/top_genes_phago_LPS.Rdata"))
de_genes_phago_LPS <- phago_DGE_results_LPS[[3]]
save(de_genes_phago_LPS, file = paste0(outdir, "/de_genes_phago_LPS.Rdata"))

#################
## untreated DGE
##################
phago_DGE_results_untreated <- run_phagocytosis_DEG("untreated")
# "5365 DEGs: 2705 up-regulated and 2660 down-regulated"

rse_non_prolif_filt_expr_phago_untreated <- phago_DGE_results_untreated[[1]]
save(rse_non_prolif_filt_expr_phago_untreated, file = paste0(outdir, "/rse_non_prolif_filt_expr_phago_untreated.Rdata"))
top_genes_phago_untreated <- phago_DGE_results_untreated[[2]]
save(top_genes_phago_untreated, file = paste0(outdir, "/top_genes_phago_untreated.Rdata"))
de_genes_phago_untreated <- phago_DGE_results_untreated[[3]]
save(de_genes_phago_untreated, file = paste0(outdir, "/de_genes_phago_untreated.Rdata"))



## Compare phagocytosis log-fractions across treatments
data <- rbind(cbind(colData(rse_non_prolif_filt_expr_phago_IFN)[, c("sample_id", "phagocytosis_mean_scaled_log_fraction")], "treatment" = "IFN"), 
              cbind(colData(rse_non_prolif_filt_expr_phago_LPS)[, c("sample_id", "phagocytosis_mean_scaled_log_fraction")], "treatment" = "LPS"),
              cbind(colData(rse_non_prolif_filt_expr_phago_untreated)[, c("sample_id", "phagocytosis_mean_scaled_log_fraction")], "treatment" = "untreated"))

set.seed(12192024)
plot <- ggplot(data = data, mapping = aes(x = treatment,
                                          y = phagocytosis_mean_scaled_log_fraction,
                                          color = treatment)) +
  geom_boxplot(size = 0.25, width=0.32, color='black', outlier.color = NA) +
  geom_jitter(width = 0.15, alpha = 1, size = 1, show.legend = F) +
  scale_color_manual(values = level_colors[["treatment"]]) +
  labs(x = "Treatment", y = "Phagocytosis scaled log-fraction") +
  theme_bw() +
  theme(axis.title = element_text(size = (8)),
        axis.text = element_text(size = (6)),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8), 
        legend.key.height= unit(0.3, 'cm'))

ggsave(filename = paste0(plotdir, "/phagocytosis_log_fractions_per_treatment.pdf"), width = 3, height = 3)



## Function to create volcano plot for DEGs
plots_DEGs <- function(treatment){
  
  results <-  eval(parse_expr(paste0("top_genes_phago_", treatment)))
  results$symbol <- rownames(results)
  
  ## Label top 10 up and down DEGs
  top_results_up <- subset(results, adj.P.Val<0.05 & logFC>0)[order(subset(results, adj.P.Val<0.05 & logFC>0)$adj.P.Val, decreasing = F), ][1:5,]
  top_results_down <- subset(results, adj.P.Val<0.05 & logFC<0)[order(subset(results, adj.P.Val<0.05 & logFC<0)$adj.P.Val, decreasing = F), ][1:5,]
  results$DEG_symbol <- apply(results, 1, function(x){if(x["symbol"] %in% c(top_results_up$symbol, top_results_down$symbol)){x["symbol"]} else {NA}}) %>%  unlist
  
  ## Significance category
  results$signif <- apply(results, 1, function(x){if(as.numeric(x["adj.P.Val"])<0.05 & as.numeric(x["logFC"]) > 0){"Up"}
    else if(as.numeric(x["adj.P.Val"])<0.05 & as.numeric(x["logFC"]) <0){"Down"}
    else if(as.numeric(x["adj.P.Val"])>=0.05){"n.s."}})

  ## DEGs
  de_genes <- subset(results, adj.P.Val<0.05)
  de_genes_up <- subset(de_genes, logFC>0)
  de_genes_down <- subset(de_genes, logFC<0)
  
  ## Volcano plot
  set.seed(10282024)
  p <-  ggplot(data = results,
               aes(x = logFC,
                   y = -log10(adj.P.Val),
                   col = signif,
                   label = DEG_symbol)) +
    geom_point() +
    theme_minimal() +
    ggrepel::geom_text_repel(color = "black", max.overlaps = Inf, size = 2.5, min.segment.length = 0,force = 2) +
    scale_color_manual(values = c("Down"="#3A5683",
                                  "n.s."= "#A7A9AE",
                                  "Up"="#74121D")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey20") +
    ylab("-log10(adj. p-val)") + 
    labs(color = "Diff. expr.", subtitle = paste(dim(de_genes)[1], "DEGs:", 
                                                 dim(de_genes_up)[1], "Up and", 
                                                 dim(de_genes_down)[1], "Down")) +
    ggtitle(paste0("DGE for phagocytosis in ", treatment)) +
    theme(plot.title = element_text(hjust=0.5, face="bold"),
          plot.subtitle = element_text(color='gray40'))
  
  return(p)
  
}

p1 <- plots_DEGs("IFN")
p2 <- plots_DEGs("LPS")
p3 <- plots_DEGs("untreated")

plot_grid(p1, p2, p3, ncol = 3, align = "hv") 
ggsave(filename = paste0(plotdir, "/volcano_plots_DEGs_3treatments.pdf"), width = 15, height = 5)



## Function to compare t-stats between treatments  
t_vs_t <- function(treatment1, treatment2){
  
  results_treatment1 <-  eval(parse_expr(paste0("top_genes_phago_", treatment1)))
  results_treatment2 <-  eval(parse_expr(paste0("top_genes_phago_", treatment2)))
  
  ## Subset to common genes
  results_treatment1_common <-  results_treatment1[intersect(rownames(results_treatment1), rownames(results_treatment2)), ]
  results_treatment2_common <- results_treatment2[match(rownames(results_treatment1_common), rownames(results_treatment2)), ]
  identical(rownames(results_treatment1_common), rownames(results_treatment2_common))
  # [1] TRUE
  
  ## Bind 
  colnames(results_treatment1_common) <- paste0(colnames(results_treatment1_common), "_treatment1")
  colnames(results_treatment2_common) <- paste0(colnames(results_treatment2_common), "_treatment2")
  merged_results <-  cbind(results_treatment1_common, results_treatment2_common)
  merged_results$symbol <- rownames(merged_results)
  
  ## Significance category for treatment 1
  for(i in seq_along(rownames(merged_results))){
    
    if(merged_results[i, "adj.P.Val_treatment1"] < 0.05){
      if(merged_results[i, "logFC_treatment1"] > 0) { merged_results$signif_treatment1[i] <- "Up" }
      else{ merged_results$signif_treatment1[i] <- "Down" }
    }
    else{
      merged_results$signif_treatment1[i] <- "n.s." 
    }
  }
  
  ## Significance category for treatment 2
  for(i in seq_along(rownames(merged_results))){
    
    if(merged_results[i, "adj.P.Val_treatment2"] < 0.05){
      if(merged_results[i, "logFC_treatment2"] > 0) { merged_results$signif_treatment2[i] <- "Up" }
      else{ merged_results$signif_treatment2[i] <- "Down" }
    }
    else{
      merged_results$signif_treatment2[i] <- "n.s." 
    }
  }
  
  ## DGE category for both treatments
  for(i in seq_along(rownames(merged_results))){
    
    if( merged_results[i, "signif_treatment1"] == "Up" ){
      if( merged_results[i, "signif_treatment2"] == "Up" ) { merged_results$signif[i] <- "Both Up" }
      else if( merged_results[i, "signif_treatment2"] == "Down"){ merged_results$signif[i] <- paste0(treatment1, " Up; ", treatment2, " Down") }
      else { merged_results$signif[i] <-  paste0(treatment1, " Up; ", treatment2, " n.s.") }
    }
    else if( merged_results[i, "signif_treatment1"] == "Down" ){
      if( merged_results[i, "signif_treatment2"] == "Up" ) { merged_results$signif[i] <- paste0(treatment1, " Down; ", treatment2, " Up") }
      else if( merged_results[i, "signif_treatment2"] == "Down"){ merged_results$signif[i] <- "Both Down" }
      else { merged_results$signif[i] <- paste0(treatment1, " Down; ", treatment2, " n.s.") }
    }
    else{
      if( merged_results[i, "signif_treatment2"] == "Up" ) { merged_results$signif[i] <- paste0(treatment1, " n.s.; ", treatment2, " Up") }
      else if( merged_results[i, "signif_treatment2"] == "Down"){ merged_results$signif[i] <- paste0(treatment1, " n.s.; ", treatment2, " Down") }
      else { merged_results$signif[i] <- "Both n.s." }
    }
    
  }
  
  ## Lable max 10 common DEGs with concordant/distinct directionality
  Both_Up_top10_symbol <- subset(merged_results[order(merged_results$adj.P.Val_treatment1, decreasing = F), ], signif == "Both Up")$symbol[1:10]
  Both_Down_top10_symbol <- subset(merged_results[order(merged_results$adj.P.Val_treatment1, decreasing = F), ], signif == "Both Down")$symbol[1:10]
  IFN_Down_untreated_Up_symbol <- subset(merged_results[order(merged_results$adj.P.Val_treatment1, decreasing = F), ], signif == paste0(treatment1, " Down; ", treatment2, " Up"))$symbol
  IFN_Up_untreated_Down_symbol <- subset(merged_results[order(merged_results$adj.P.Val_treatment1, decreasing = F), ], signif == paste0(treatment1, " Up; ", treatment2, " Down"))$symbol
  
  ## Label only top 10 DEGs with different direction in IFN vs untreated
  if(treatment1 == "IFN" & treatment2 == "untreated"){
    IFN_Down_untreated_Up_symbol <- subset(merged_results[order(merged_results$adj.P.Val_treatment1, decreasing = F), ], signif == paste0(treatment1, " Down; ", treatment2, " Up"))$symbol[1:10]
    IFN_Up_untreated_Down_symbol <- subset(merged_results[order(merged_results$adj.P.Val_treatment1, decreasing = F), ], signif == paste0(treatment1, " Up; ", treatment2, " Down"))$symbol[1:10]
    
  }
  
  merged_results$DEG_symbol <- apply(merged_results, 1, function(x){if(x["symbol"] %in% 
                                                                       c(Both_Up_top10_symbol, Both_Down_top10_symbol, 
                                                                         IFN_Down_untreated_Up_symbol, IFN_Up_untreated_Down_symbol)){x["symbol"]} 
    else {NA}}) %>%  unlist

  ## Pearson corr between t-stats
  r <- cor(merged_results$t_treatment1, merged_results$t_treatment2, method = "pearson")
  r_anno = paste0("r = ", format(round(r, 2), nsmall = 2))
  
  ## Colors and transparency for DGE cat
  cols <- c("red3", "royalblue3", 
            "salmon1", "lightskyblue2", 
            "darksalmon", "slategray2", 
            "thistle2", "thistle2", "grey") 
  
  names(cols) <- c("Both Up",  
                   "Both Down",
                   paste0(treatment1, " Up; ", treatment2, " n.s."),
                   paste0(treatment1, " Down; ", treatment2, " n.s."),
                   paste0(treatment1, " n.s.; ", treatment2, " Up"),
                   paste0(treatment1, " n.s.; ", treatment2, " Down"),
                   paste0(treatment1, " Up; ", treatment2, " Down"), 
                   paste0(treatment1, " Down; ", treatment2, " Up"),
                   "Both n.s.")
  
  alphas <- c(1, 1, 0.8, 0.8, 0.8, 0.8, 1, 1, 0.3)  
  names(alphas) <- names(cols)
  ## Order
  merged_results$signif <-  factor(merged_results$signif, levels = names(cols)[which(names(cols) %in% unique(merged_results$signif))])
  
  cat_labels <- c(paste0("Up | Up (n = ", dim(subset(merged_results, signif == "Both Up"))[1], ")"),
                  paste0("Down | Down (n = ", dim(subset(merged_results, signif == "Both Down"))[1], ")"), 
                  paste0("Up | n.s. (n = ", dim(subset(merged_results, signif == paste0(treatment1, " Up; ", treatment2, " n.s.")))[1], ")"), 
                  paste0("Down | n.s. (n = ", dim(subset(merged_results, signif == paste0(treatment1, " Down; ", treatment2, " n.s.")))[1], ")"),
                  paste0("n.s. | Up (n = ", dim(subset(merged_results, signif == paste0(treatment1, " n.s.; ", treatment2, " Up")))[1], ")"), 
                  paste0("n.s. | Down (n = ", dim(subset(merged_results, signif == paste0(treatment1, " n.s.; ", treatment2, " Down")))[1], ")"), 
                  paste0("Up | Down (n = ", dim(subset(merged_results, signif == paste0(treatment1, " Up; ", treatment2, " Down")))[1], ")"), 
                  paste0("Down | Up (n = ", dim(subset(merged_results, signif == paste0(treatment1, " Down; ", treatment2, " Up")))[1], ")"), 
                  paste0("n.s. (n = ", dim(subset(merged_results, signif == "Both n.s."))[1], ")"))
  names(cat_labels) <- names(cols)
  
  p <- ggplot(merged_results, aes(x = t_treatment1, y = t_treatment2, 
                                  color = signif, alpha = signif, 
                                  label = DEG_symbol)) +
    geom_point(size = 1.5) +
    scale_color_manual(values = cols, labels = cat_labels, drop = FALSE) + 
    scale_alpha_manual(values = alphas, labels = cat_labels, drop=FALSE) +
    ggrepel::geom_text_repel(color = "black", max.overlaps = Inf, size = 2.5, min.segment.length = 0, fontface = "bold", force = 2) +
    labs(x = paste("t-stats in", treatment1), 
         y = paste("t-stats in", treatment2),
         title = paste("DGE effects in", treatment1, "vs.", treatment2),
         subtitle = r_anno, 
         color = paste0("DGE result (", treatment1, " | ", treatment2, ")")) +
    guides(alpha = 'none', color = guide_legend(override.aes = list(size=2))) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey20") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey20") +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 11),
          plot.subtitle = element_text(size = 9),
          plot.margin = unit(c(1,1,1,1), "cm"),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 9),
          legend.text = element_text(size=9),
          legend.title = element_text(size=10))
  
  
  return(list(p, merged_results))
  
}

p1 <- t_vs_t("IFN", "LPS")[[1]]
merged_IFN_vs_LPS <- t_vs_t("IFN", "LPS")[[2]]
ggsave(filename = paste0(plotdir, "/t_vs_t_IFN_vs_LPS.pdf"), width = 8, height = 5)
p2 <- t_vs_t("IFN", "untreated")[[1]]
merged_IFN_vs_untreated <- t_vs_t("IFN", "untreated")[[2]]
ggsave(filename = paste0(plotdir, "/t_vs_t_IFN_vs_untreated.pdf"), width = 8, height = 5)
p3 <- t_vs_t("LPS", "untreated")[[1]]
merged_LPS_vs_untreated <- t_vs_t("LPS", "untreated")[[2]]
ggsave(filename = paste0(plotdir, "/t_vs_t_LPS_vs_untreated.pdf"), width = 8, height = 5)



## Explore corr between covariates and phago log-fractions
hs <- list()
for(treatment in c("IFN", "LPS", "untreated")){
  
  rse <-  eval(parse_expr(paste0("rse_non_prolif_filt_expr_phago_", treatment)))
  formula <- ~ phagocytosis_mean_scaled_log_fraction + Sex + phagocytosis_mean_min_prop + genotype_PC1 + genotype_PC2 + (1|pool)
  C = canCorPairs(formula, colData(rse))
  
  ## Heatmap
  hs[[treatment]] <- pheatmap(
        mat = C, 
        main = treatment,
        name = "Cor",
        color = hcl.colors(50, "Reds 3", rev = TRUE),
        display_numbers = T,
        fontsize = 11,
        border_color = "black"
      )
  assign(paste0("CCA_", treatment), C)
}

par(mfrow=c(1,3))
pdf(file = paste0(plotdir, "/CCA_3treatments.pdf"), height = 5.5, width = 6)
hs[[1]]
hs[[2]]
hs[[3]]
dev.off()


## Plot relationships between phagocytosis and covariates
explore_phago_vs_covariate <- function(treatment, covariate){
  
  rse <-  eval(parse_expr(paste0("rse_non_prolif_filt_expr_phago_", treatment)))
  CCA <- eval(parse_expr(paste0("CCA_", treatment)))
  
  ## Order levels of variables
  rse$pool <- factor(rse$pool, levels = paste0("pool", c("3", "4", "5", "6", "7", "8", "9", "10",
                                                         "11", "13", "14", "15", "16", "17")))
  rse$Sex <- factor(rse$Sex, levels = c("Female", "Male"))
  
  color_var = "treatment"
  shape_var = "Sex"
  
  if(covariate == "Sex"){
    lab = "Sex"
    color_var = "Sex"
    shape_var = ""}
  else if(covariate == "pool"){
    lab = "Pool"
    color_var = "pool"
    shape_var = ""
  }
  else if(covariate == "genotype_PC1"){
    lab = "genotype PC1"}
  else if(covariate == "genotype_PC2"){
    lab = "genotype PC2"}
  else if(covariate == "phagocytosis_mean_min_prop"){
    lab = "Phagocytosis min proportion"}
    
  data <- colData(rse)
  
  ## Scatterplot for continuous variable
  if(type(data[, covariate]) ==  "double"){
    plot <- ggplot(data, aes(x = phagocytosis_mean_scaled_log_fraction,
                             y = !!rlang::sym(covariate),
                             color = !!rlang::sym(color_var))) +
      geom_point(aes(shape = !!rlang::sym(shape_var)), size=1) +
      stat_smooth (geom="line", alpha=0.6, size=0.5, span=0.25, method = lm, color='orangered3') +
      scale_color_manual(values = level_colors[[color_var]]) +
      scale_shape_manual(values = level_shapes[[shape_var]]) +
      theme_bw() +
      labs(x = "Phagocytosis scaled log-fraction", y = lab,
           title = paste("Correlation in", treatment),
           subtitle = paste("Cor:", signif(CCA["phagocytosis_mean_scaled_log_fraction", covariate], digits = 3)) ) +
      theme(plot.title = element_text(size = 8, face = "bold"),
            plot.subtitle = element_text(size = 7),
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 6),
            legend.text = element_text(size=7),
            legend.title = element_text(size=8), 
            legend.key.height= unit(0.3, 'cm'))
  }
  
  ## Boxplot for categorical variable
  else if(type(data[, covariate]) ==  "character"){
    set.seed(12192024)
    plot <- ggplot(data = data, mapping = aes(x = !!rlang::sym(covariate),
                                              y = phagocytosis_mean_scaled_log_fraction,
                                              color = !!rlang::sym(color_var))) +
      geom_boxplot(size = 0.25, width=0.32, color='black', outlier.color = NA) +
      geom_jitter(aes(shape = !!rlang::sym(shape_var)), width = 0.15, alpha = 1, size = 1) +
      scale_color_manual(values = level_colors[[color_var]]) +
      scale_shape_manual(values = level_shapes[[shape_var]]) +
      theme_bw() +
      labs(y = "Phagocytosis scaled log-fraction", x = lab, 
           title = paste("Correlation in", treatment),
           subtitle = paste("Cor:", signif(CCA["phagocytosis_mean_scaled_log_fraction", covariate], digits = 3)) )+
      theme(plot.title = element_text(size = 8, face = "bold"),
            plot.subtitle = element_text(size = 7),
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 6),
            legend.text = element_text(size = 6),
            legend.title = element_text(size = 8), 
            legend.key.height= unit(0.3, 'cm'))
  }
  
  return(plot)
}


plots <- list()
i = 1
for(treatment in c("IFN", "LPS", "untreated")){
  for(cov in c("Sex", "pool", "genotype_PC1", "genotype_PC2", "phagocytosis_mean_min_prop")){
    plots[[i]] <- explore_phago_vs_covariate(treatment, cov)
    i = i +1
  }
}

plot_grid(plotlist = plots, ncol = 5, align = "vh", rel_widths = c(0.7, 1, 0.7, 0.7, 0.7))
ggsave(filename = paste0(plotdir, "/phago_scaled_log_fractions_vs_cov_3treatments.pdf"), width = 21, height = 9)



## Plot gene expr vs phagocytosis in 3 treatments
plot_gene_expr_all_treatments <- function(gene_id){
  
  color_var <- "treatment"
  shape_var <- "Sex"
  
  rse_IFN <- rse_non_prolif_filt_expr_phago_IFN
  rse_LPS <- rse_non_prolif_filt_expr_phago_LPS
  rse_untreated <- rse_non_prolif_filt_expr_phago_untreated
  
  ## sample data per treatment
  sample_data_IFN <- colData(rse_IFN)
  sample_data_LPS <- colData(rse_LPS)
  sample_data_untreated <- colData(rse_untreated)
  
  ## gene log-cpm per treatment 
  logcounts_IFN <- assays(rse_IFN)$logcounts[gene_id, ]
  logcounts_LPS <- assays(rse_LPS)$logcounts[gene_id, ]
  logcounts_untreated <- assays(rse_untreated)$logcounts[gene_id, ]
  
  ## Add gene expr to sample data
  sample_data_IFN <-cbind(sample_data_IFN, gene_expr = logcounts_IFN) %>% as.data.frame()
  sample_data_LPS <-cbind(sample_data_LPS, gene_expr = logcounts_LPS) %>% as.data.frame()
  sample_data_untreated <-cbind(sample_data_untreated, gene_expr = logcounts_untreated) %>% as.data.frame()
  
  ## Min and max gene expr
  min <- min(c(sample_data_IFN$gene_expr, sample_data_LPS$gene_expr, sample_data_untreated$gene_expr))
  max <- max(c(sample_data_IFN$gene_expr, sample_data_LPS$gene_expr, sample_data_untreated$gene_expr))
  
  ## IFN
  plot1 <- ggplot(sample_data_IFN, aes(x = phagocytosis_mean_scaled_log_fraction,
                                       y = gene_expr,
                                       color = !!rlang::sym(color_var))) +
    geom_point(aes(shape = !!rlang::sym(shape_var)), size=1, show.legend = F) +
    stat_smooth (geom = "line", alpha=1, size=0.5, span=0.45, method = lm, color='orangered3') +
    scale_color_manual(values = level_colors[[color_var]]) +
    scale_shape_manual(values = level_shapes[[shape_var]]) +
    theme_bw() +
    ylim(min, max) +
    labs(title = paste(gene_id, "in IFN"), subtitle = paste0("adj. P = ", signif(top_genes_phago_IFN[gene_id, "adj.P.Val"], digits = 3), "      ",
                                            "logFC = ", signif(top_genes_phago_IFN[gene_id, "logFC"], digits = 3)), 
         x = "Phagocytosis scaled log-fraction", y = "log-CPM") +
    theme(plot.title = element_text(face = "bold", size = 9),
          plot.subtitle = element_text(size = 8, color = "gray40"),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 6))
  
  ## LPS
  plot2 <- ggplot(sample_data_LPS, aes(x = phagocytosis_mean_scaled_log_fraction,
                                       y = gene_expr,
                                       color = !!rlang::sym(color_var))) +
    geom_point(aes(shape = !!rlang::sym(shape_var)), size=1, show.legend = F) +
    stat_smooth (geom = "line", alpha=1, size=0.5, span=0.45, method = lm, color='orangered3') +
    scale_color_manual(values = level_colors[[color_var]]) +
    scale_shape_manual(values = level_shapes[[shape_var]]) +
    theme_bw() +
    ylim(min, max) +
    labs(title = paste(gene_id, "in LPS"), subtitle = paste0("adj. P = ", signif(top_genes_phago_LPS[gene_id, "adj.P.Val"], digits = 3), "      ",
                                            "logFC = ", signif(top_genes_phago_LPS[gene_id, "logFC"], digits = 3)), 
         x = "Phagocytosis scaled log-fraction", y = "log-CPM") +
    theme(plot.title = element_text(face = "bold", size = 9),
          plot.subtitle = element_text(size = 8, color = "gray40"),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 6))
  
  ## untreated
  plot3 <- ggplot(sample_data_untreated, aes(x = phagocytosis_mean_scaled_log_fraction,
                                             y = gene_expr,
                                             color = !!rlang::sym(color_var))) +
    geom_point(aes(shape = !!rlang::sym(shape_var)), size=1) +
    stat_smooth (geom="line", alpha=1, size=0.5, span=0.45, method = lm, color='orangered3') +
    scale_color_manual(values = level_colors[[color_var]]) +
    scale_shape_manual(values = level_shapes[[shape_var]]) +
    theme_bw() +
    ylim(min, max) +
    guides(color = "none") +
    labs(title = paste(gene_id, "in untreated"), subtitle = paste0("adj. P = ", signif(top_genes_phago_untreated[gene_id, "adj.P.Val"], digits = 3), "      ",
                                                                  "logFC = ", signif(top_genes_phago_untreated[gene_id, "logFC"], digits = 3)), 
         x = "Phagocytosis scaled log-fraction", y = "log-CPM") +
    theme(plot.title = element_text(face = "bold", size = 9),
          plot.subtitle = element_text(size = 8, color = "gray40"),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 6),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 8), 
          legend.key.height= unit(0.3, 'cm'))
  
  return(list(plot1, plot2, plot3))
  
}


## Top DEGs Up in IFN and non-DE in LPS
up_IFN_only <- subset(merged_IFN_vs_LPS, signif == "IFN Up; LPS n.s.")
up_IFN_only <- up_IFN_only[order(up_IFN_only$adj.P.Val_treatment1, decreasing = F), "symbol"][1:10]
plots <- list()
i = 1
for(gene_id in up_IFN_only){
  
  plots[i:(i+2)] <- plot_gene_expr_all_treatments(gene_id)[1:3]
  i = i + 3
}
plot_grid(plotlist = plots, nrow = 10, align = "vh")
ggsave(filename = paste0(plotdir, "/topUpDEGs_IFN_and_nonDEGs_LPS.pdf"), width = 12, height = 27)


## Top DEGs Down in IFN and non-DE in LPS
down_IFN_only <- subset(merged_IFN_vs_LPS, signif == "IFN Down; LPS n.s.")
down_IFN_only <- down_IFN_only[order(down_IFN_only$adj.P.Val_treatment1, decreasing = F), "symbol"][1:10]
plots <- list()
i = 1
for(gene_id in down_IFN_only){
  
  plots[i:(i+2)] <- plot_gene_expr_all_treatments(gene_id)[1:3]
  i = i + 3
}
plot_grid(plotlist = plots, nrow = 10, align = "vh")
ggsave(filename = paste0(plotdir, "/topDownDEGs_IFN_and_nonDEGs_LPS.pdf"), width = 12, height = 27)


## Top DEGs Up in untreated and non-DE in LPS
up_untreated_only <- subset(merged_LPS_vs_untreated, signif == "LPS n.s.; untreated Up")
up_untreated_only <- up_untreated_only[order(up_untreated_only$adj.P.Val_treatment2, decreasing = F), "symbol"][1:10]
plots <- list()
i = 1
for(gene_id in up_untreated_only){
  
  plots[i:(i+2)] <- plot_gene_expr_all_treatments(gene_id)[1:3]
  i = i + 3
}
plot_grid(plotlist = plots, nrow = 10, align = "vh")
ggsave(filename = paste0(plotdir, "/topUpDEGs_untreated_and_nonDEGs_LPS.pdf"), width = 12, height = 27)


## Top DEGs Down in untreated and non-DE in LPS
down_untreated_only <- subset(merged_LPS_vs_untreated, signif == "LPS n.s.; untreated Down")
down_untreated_only <- down_untreated_only[order(down_untreated_only$adj.P.Val_treatment2, decreasing = F), "symbol"][1:10]
plots <- list()
i = 1
for(gene_id in down_untreated_only){
  
  plots[i:(i+2)] <- plot_gene_expr_all_treatments(gene_id)[1:3]
  i = i + 3
}
plot_grid(plotlist = plots, nrow = 10, align = "vh")
ggsave(filename = paste0(plotdir, "/topDownDEGs_untreated_and_nonDEGs_LPS.pdf"), width = 12, height = 27)







## Reproducibility info
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.1 (2023-06-16)
# os       Ubuntu 22.04.4 LTS
# system   x86_64, linux-gnu
# ui       RStudio
# language (EN)
# collate  en_GB.UTF-8
# ctype    en_GB.UTF-8
# tz       Europe/Belfast
# date     2024-12-19
# rstudio  2024.04.0+735 Chocolate Cosmos (server)
# pandoc   3.1.12.3 @ /opt/view/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version     date (UTC) lib source
# abind                  1.4-5       2016-07-21 [1] CRAN (R 4.3.1)
# AnnotationDbi        * 1.64.1      2023-11-03 [1] Bioconductor
# aod                    1.3.3       2023-12-13 [1] CRAN (R 4.3.1)
# backports              1.4.1       2021-12-13 [1] CRAN (R 4.3.1)
# Biobase              * 2.62.0      2023-10-24 [1] Bioconductor
# BiocFileCache          2.10.2      2024-03-27 [1] Bioconductor 3.18 (R 4.3.1)
# BiocGenerics         * 0.48.1      2023-11-01 [1] Bioconductor
# BiocParallel         * 1.36.0      2023-10-24 [1] Bioconductor
# biomaRt              * 2.58.2      2024-01-30 [1] Bioconductor 3.18 (R 4.3.1)
# Biostrings             2.70.3      2024-03-13 [1] Bioconductor 3.18 (R 4.3.1)
# bit                    4.0.5       2022-11-15 [1] CRAN (R 4.3.1)
# bit64                  4.0.5       2020-08-30 [1] CRAN (R 4.3.1)
# bitops                 1.0-7       2021-04-24 [1] CRAN (R 4.3.1)
# blob                   1.2.4       2023-03-17 [1] CRAN (R 4.3.1)
# boot                   1.3-30      2024-02-26 [1] CRAN (R 4.3.1)
# broom                  1.0.5       2023-06-09 [1] CRAN (R 4.3.1)
# cachem                 1.0.8       2023-05-01 [1] CRAN (R 4.3.1)
# car                    3.1-2       2023-03-30 [1] CRAN (R 4.3.1)
# carData                3.0-5       2022-01-06 [1] CRAN (R 4.3.1)
# caTools                1.18.2      2021-03-28 [1] CRAN (R 4.3.1)
# circlize             * 0.4.16      2024-02-20 [1] CRAN (R 4.3.1)
# cli                    3.6.2       2023-12-11 [1] CRAN (R 4.3.1)
# clue                   0.3-65      2023-09-23 [1] CRAN (R 4.3.1)
# cluster                2.1.6       2023-12-01 [1] CRAN (R 4.3.1)
# codetools              0.2-20      2024-03-31 [1] CRAN (R 4.3.1)
# colorspace             2.1-0       2023-01-23 [1] CRAN (R 4.3.1)
# ComplexHeatmap       * 2.18.0      2023-10-24 [1] Bioconductor
# corpcor                1.6.10      2021-09-16 [1] CRAN (R 4.3.1)
# cowplot              * 1.1.3       2024-01-22 [1] CRAN (R 4.3.1)
# crayon                 1.5.2       2022-09-29 [1] CRAN (R 4.3.1)
# curl                   5.2.1       2024-03-01 [1] CRAN (R 4.3.1)
# data.table             1.15.4      2024-03-30 [1] CRAN (R 4.3.1)
# DBI                    1.2.2       2024-02-16 [1] CRAN (R 4.3.1)
# dbplyr                 2.5.0       2024-03-19 [1] CRAN (R 4.3.1)
# DelayedArray           0.28.0      2023-10-24 [1] Bioconductor
# digest                 0.6.35      2024-03-11 [1] CRAN (R 4.3.1)
# doParallel             1.0.17      2022-02-07 [1] CRAN (R 4.3.1)
# dplyr                * 1.1.4       2023-11-17 [1] CRAN (R 4.3.1)
# edgeR                * 4.0.16      2024-02-18 [1] Bioconductor 3.18 (R 4.3.1)
# EnvStats               2.8.1       2023-08-22 [1] CRAN (R 4.3.1)
# fANCOVA                0.6-1       2020-11-13 [1] CRAN (R 4.3.1)
# fansi                  1.0.6       2023-12-08 [1] CRAN (R 4.3.1)
# farver                 2.1.1       2022-07-06 [1] CRAN (R 4.3.1)
# fastmap                1.1.1       2023-02-24 [1] CRAN (R 4.3.1)
# fastmatch              1.1-4       2023-08-18 [1] CRAN (R 4.3.1)
# fgsea                * 1.28.0      2023-10-24 [1] Bioconductor
# filelock               1.0.3       2023-12-11 [1] CRAN (R 4.3.1)
# forcats              * 1.0.0       2023-01-29 [1] CRAN (R 4.3.1)
# foreach                1.5.2       2022-02-02 [1] CRAN (R 4.3.1)
# generics               0.1.3       2022-07-05 [1] CRAN (R 4.3.1)
# GenomeInfoDb         * 1.38.8      2024-03-15 [1] Bioconductor 3.18 (R 4.3.1)
# GenomeInfoDbData       1.2.11      2024-11-06 [1] Bioconductor
# GenomicRanges        * 1.54.1      2023-10-29 [1] Bioconductor
# GetoptLong             1.0.5       2020-12-15 [1] CRAN (R 4.3.1)
# ggplot2              * 3.5.1       2024-04-23 [1] CRAN (R 4.3.1)
# ggpubr                 0.6.0       2023-02-10 [1] CRAN (R 4.3.1)
# ggrepel                0.9.5       2024-01-10 [1] CRAN (R 4.3.1)
# ggsignif               0.6.4       2022-10-13 [1] CRAN (R 4.3.1)
# GlobalOptions          0.1.2       2020-06-10 [1] CRAN (R 4.3.1)
# glue                   1.7.0       2024-01-09 [1] CRAN (R 4.3.1)
# gplots                 3.1.3.1     2024-02-02 [1] CRAN (R 4.3.1)
# gtable                 0.3.4       2023-08-21 [1] CRAN (R 4.3.1)
# gtools                 3.9.5       2023-11-20 [1] CRAN (R 4.3.1)
# hms                    1.1.3       2023-03-21 [1] CRAN (R 4.3.1)
# httr                   1.4.7       2023-08-15 [1] CRAN (R 4.3.1)
# IRanges              * 2.36.0      2023-10-24 [1] Bioconductor
# iterators              1.0.14      2022-02-05 [1] CRAN (R 4.3.1)
# KEGGREST               1.42.0      2023-10-24 [1] Bioconductor
# KernSmooth             2.23-22     2023-07-10 [1] CRAN (R 4.3.1)
# labeling               0.4.3       2023-08-29 [1] CRAN (R 4.3.1)
# lattice                0.22-6      2024-03-20 [1] CRAN (R 4.3.1)
# lifecycle              1.0.4       2023-11-07 [1] CRAN (R 4.3.1)
# limma                * 3.58.1      2023-10-31 [1] Bioconductor
# lme4                   1.1-35.2    2024-03-28 [1] CRAN (R 4.3.1)
# lmerTest               3.1-3       2020-10-23 [1] CRAN (R 4.3.1)
# locfit                 1.5-9.9     2024-03-01 [1] CRAN (R 4.3.1)
# lubridate            * 1.9.3       2023-09-27 [1] CRAN (R 4.3.1)
# magrittr               2.0.3       2022-03-30 [1] CRAN (R 4.3.1)
# MASS                   7.3-60.0.1  2024-01-13 [1] CRAN (R 4.3.1)
# Matrix                 1.6-5       2024-01-11 [1] CRAN (R 4.3.1)
# MatrixGenerics       * 1.14.0      2023-10-24 [1] Bioconductor
# matrixStats          * 1.2.0       2023-12-11 [1] CRAN (R 4.3.1)
# memoise                2.0.1       2021-11-26 [1] CRAN (R 4.3.1)
# mgcv                   1.9-1       2023-12-21 [1] CRAN (R 4.3.1)
# minqa                  1.2.6       2023-09-11 [1] CRAN (R 4.3.1)
# munsell                0.5.1       2024-04-01 [1] CRAN (R 4.3.1)
# mvtnorm                1.2-4       2023-11-27 [1] CRAN (R 4.3.1)
# nlme                   3.1-164     2023-11-27 [1] CRAN (R 4.3.1)
# nloptr                 2.0.3       2022-05-26 [1] CRAN (R 4.3.1)
# numDeriv               2016.8-1.1  2019-06-06 [1] CRAN (R 4.3.1)
# org.Hs.eg.db         * 3.18.0      2024-11-06 [1] Bioconductor
# pbkrtest               0.5.2       2023-01-19 [1] CRAN (R 4.3.1)
# pillar                 1.9.0       2023-03-22 [1] CRAN (R 4.3.1)
# pkgconfig              2.0.3       2019-09-22 [1] CRAN (R 4.3.1)
# plyr                   1.8.9       2023-10-02 [1] CRAN (R 4.3.1)
# png                    0.1-8       2022-11-29 [1] CRAN (R 4.3.1)
# prettyunits            1.2.0       2023-09-24 [1] CRAN (R 4.3.1)
# progress               1.2.3       2023-12-06 [1] CRAN (R 4.3.1)
# purrr                * 1.0.2       2023-08-10 [1] CRAN (R 4.3.1)
# R6                     2.5.1       2021-08-19 [1] CRAN (R 4.3.1)
# ragg                   1.3.0       2024-03-13 [1] CRAN (R 4.3.1)
# rappdirs               0.3.3       2021-01-31 [1] CRAN (R 4.3.1)
# rbibutils              2.2.16      2023-10-25 [1] CRAN (R 4.3.1)
# RColorBrewer           1.1-3       2022-04-03 [1] CRAN (R 4.3.1)
# Rcpp                   1.0.12      2024-01-09 [1] CRAN (R 4.3.1)
# RCurl                  1.98-1.14   2024-01-09 [1] CRAN (R 4.3.1)
# Rdpack                 2.6         2023-11-08 [1] CRAN (R 4.3.1)
# readr                * 2.1.5       2024-01-10 [1] CRAN (R 4.3.1)
# remaCor                0.0.18      2024-02-08 [1] CRAN (R 4.3.1)
# reshape2               1.4.4       2020-04-09 [1] CRAN (R 4.3.1)
# RhpcBLASctl            0.23-42     2023-02-11 [1] CRAN (R 4.3.1)
# rjson                  0.2.21      2022-01-09 [1] CRAN (R 4.3.1)
# rlang                * 1.1.3       2024-01-10 [1] CRAN (R 4.3.1)
# RSQLite                2.3.6       2024-03-31 [1] CRAN (R 4.3.1)
# rstatix                0.7.2       2023-02-01 [1] CRAN (R 4.3.1)
# rstudioapi             0.16.0      2024-03-24 [1] CRAN (R 4.3.1)
# S4Arrays               1.2.1       2024-03-04 [1] Bioconductor 3.18 (R 4.3.1)
# S4Vectors            * 0.40.2      2023-11-23 [1] Bioconductor 3.18 (R 4.3.1)
# scales                 1.3.0       2023-11-28 [1] CRAN (R 4.3.1)
# sessioninfo          * 1.2.2       2021-12-06 [1] CRAN (R 4.3.1)
# shape                  1.4.6.1     2024-02-23 [1] CRAN (R 4.3.1)
# SparseArray            1.2.4       2024-02-11 [1] Bioconductor 3.18 (R 4.3.1)
# statmod                1.5.0       2023-01-06 [1] CRAN (R 4.3.1)
# stringi                1.8.3       2023-12-11 [1] CRAN (R 4.3.1)
# stringr              * 1.5.1       2023-11-14 [1] CRAN (R 4.3.1)
# SummarizedExperiment * 1.32.0      2023-10-24 [1] Bioconductor
# systemfonts            1.0.6       2024-03-07 [1] CRAN (R 4.3.1)
# textshaping            0.3.7       2023-10-09 [1] CRAN (R 4.3.1)
# tibble               * 3.2.1       2023-03-20 [1] CRAN (R 4.3.1)
# tidyr                * 1.3.1       2024-01-24 [1] CRAN (R 4.3.1)
# tidyselect             1.2.1       2024-03-11 [1] CRAN (R 4.3.1)
# tidyverse            * 2.0.0       2023-02-22 [1] CRAN (R 4.3.1)
# timechange             0.3.0       2024-01-18 [1] CRAN (R 4.3.1)
# tzdb                   0.4.0       2023-05-12 [1] CRAN (R 4.3.1)
# utf8                   1.2.4       2023-10-22 [1] CRAN (R 4.3.1)
# variancePartition    * 1.32.5      2024-02-16 [1] Bioconductor 3.18 (R 4.3.1)
# vctrs                  0.6.5       2023-12-01 [1] CRAN (R 4.3.1)
# withr                  3.0.0       2024-01-16 [1] CRAN (R 4.3.1)
# XML                    3.99-0.16.1 2024-01-22 [1] CRAN (R 4.3.1)
# xml2                   1.3.6       2023-12-04 [1] CRAN (R 4.3.1)
# XVector                0.42.0      2023-10-24 [1] Bioconductor
# zlibbioc               1.48.2      2024-03-13 [1] Bioconductor 3.18 (R 4.3.1)
# 
# [1] /opt/view/rlib/R/library
# [2] /software/hgi/installs/softpack/rstudio/.spack/linux-ubuntu22.04-x86_64_v3/gcc-11.4.0/r-4.3.1-bfwldrk76z6f52upk47zepliekn7ayqz/rlib/R/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────









