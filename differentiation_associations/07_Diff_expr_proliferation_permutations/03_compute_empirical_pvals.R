
library(tidyverse)
library(dplyr)
library(rlang)
library(ggplot2)
library(cowplot)
library(sessioninfo)

#-------------------------------------------------------------------------------
#           7.3 Compute empirical p-values per treatment and gene   
#-------------------------------------------------------------------------------
#  Code to compute empirical gene- and treatment-level p-values based on 
#  the DGE results for the 100 permuted proliferation phenotypes.
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

## Set working dir
setwd("/lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_Daianna/")

## Define output data and plot dir 
outdir = paste(getwd(), "output_data", "07_Diff_expr_proliferation_permutations", "03_compute_empirical_pvals", sep = "/")
plotdir = paste(getwd(), "plots", "07_Diff_expr_proliferation_permutations", "03_compute_empirical_pvals", sep = "/")
dir.create(outdir, recursive = T)
dir.create(plotdir, recursive = T)

## Input dir
input_dir00 = paste(getwd(), "output_data", "01_Diff_expr_PRS_permutations", "00_Data_exploration_processing", sep = "/")
input_dir01 = paste(getwd(), "output_data", "07_Diff_expr_proliferation_permutations", "01_run_proliferation_DGE", sep = "/")
input_dir02 = paste(getwd(), "output_data", "07_Diff_expr_proliferation_permutations", "02_permuted_proliferation_DGE", sep = "/")

## Load color dict
load(paste0(input_dir00, "/level_colors.Rdata"), verbose = T)

## Load data used for DGE
load(paste0(input_dir01, "/rse_non_prolif_filt_expr_proliferation_IFN.Rdata"), verbose = T)
load(paste0(input_dir01, "/rse_non_prolif_filt_expr_proliferation_LPS.Rdata"), verbose = T)
load(paste0(input_dir01, "/rse_non_prolif_filt_expr_proliferation_untreated.Rdata"), verbose = T)

## Load original DGE results
load(paste0(input_dir01, "/top_genes_proliferation_IFN.Rdata"), verbose = T)
load(paste0(input_dir01, "/top_genes_proliferation_LPS.Rdata"), verbose = T)
load(paste0(input_dir01, "/top_genes_proliferation_untreated.Rdata"), verbose = T)

## Load results for 100 permutations
permutation_results_IFN <- list()
permutation_results_LPS<- list()
permutation_results_untreated <- list()

for(i in 1:100){
  
  load(paste0(input_dir02,"/dge_results_permutation_", i, ".Rdata"), verbose = T)
  assign(paste0("dge_results_", i), dge_permutation_results)
  
  permutation_results_IFN[[i]] <- eval(parse_expr(paste("dge_results_", i, "$IFN[, c(\"logFC\", \"t\", \"P.Value\", \"adj.P.Val\")]", sep = "")))
  colnames(permutation_results_IFN[[i]]) <- paste0(colnames(permutation_results_IFN[[i]]), "_", i)
  permutation_results_LPS[[i]] <- eval(parse_expr(paste("dge_results_", i, "$LPS[, c(\"logFC\", \"t\", \"P.Value\", \"adj.P.Val\")]", sep = "")))
  colnames(permutation_results_LPS[[i]]) <- paste0(colnames(permutation_results_LPS[[i]]), "_", i)
  permutation_results_untreated[[i]] <- eval(parse_expr(paste("dge_results_", i, "$untreated[, c(\"logFC\", \"t\", \"P.Value\", \"adj.P.Val\")]", sep = "")))
  colnames(permutation_results_untreated[[i]]) <- paste0(colnames(permutation_results_untreated[[i]]), "_", i)
}

## Check num genes per treatment
lapply(permutation_results_IFN, dim) %>% unique()
# [1] 12757    4
lapply(permutation_results_LPS, dim) %>% unique()
# [1] 12823    4
lapply(permutation_results_untreated, dim) %>% unique()
# [1] 13121    4

## Bind results per treatment
permutation_results_IFN <- do.call(cbind, permutation_results_IFN)
permutation_results_LPS <- do.call(cbind, permutation_results_LPS)
permutation_results_untreated <- do.call(cbind, permutation_results_untreated)



## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#                 7.3.1 Compute empirical p-value per treatment
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

## Function to compute #DEGs for each permutation and treatment empirical p
treatment_empirical_p <- function(treatment){
  
  permutation_results <- eval(parse_expr(paste0("permutation_results_", treatment)))
  
  ## Observed DEGs
  original_DGE_results <- eval(parse_expr(paste0("top_genes_proliferation_", treatment)))
  original_DGE_results$symbol <-  rownames(original_DGE_results)
  observed_DEGs <- subset(original_DGE_results,  adj.P.Val<0.05)
  
  ## Number of DEGs per permutation
  num_DEGs <-  apply(permutation_results[, paste0("adj.P.Val_", 1:100)], 2, function(i){length(which(i<0.05))})
  
  ## N = number of permutations yielding >= observed DEGs
  message(sum(num_DEGs >= dim(observed_DEGs)[1]))
  # IFN: 0
  # LPS: 0
  # untreated: 0
  
  ## Empirical p = (N + 1)/[number of permutations + 1]
  p <- (sum(num_DEGs >= dim(observed_DEGs)[1]) + 1)/(length(num_DEGs) + 1)
  # IFN: 0.00990099
  # LPS: 0.00990099
  # untreated: 0.00990099
  
  ## Num of overlapping DEGs 
  num_overlapping_genes <- apply(permutation_results[, paste0("adj.P.Val_", 1:100)], 2, 
                                 function(i){paste0(length(intersect(names(which(i<0.05)), observed_DEGs$symbol)), "/", length(which(i<0.05)))})
  num_overlapping_genes <- num_overlapping_genes[names(which(gsub("/.*", "", num_overlapping_genes)>1))]
  
  ## IFN:
  
  # adj.P.Val_2  adj.P.Val_4  adj.P.Val_7  adj.P.Val_9 adj.P.Val_11 adj.P.Val_13 adj.P.Val_17 adj.P.Val_20 adj.P.Val_21 
  #     "13/28"       "4/11"     "97/159"        "5/6"    "141/150"      "22/39"        "2/2"      "21/26"        "6/6" 
  # adj.P.Val_22 adj.P.Val_26 adj.P.Val_30 adj.P.Val_32 adj.P.Val_36 adj.P.Val_38 adj.P.Val_41 adj.P.Val_44 adj.P.Val_49 
  #       "2/11"      "71/99"      "41/90"        "2/6"     "62/103"        "3/3"      "28/50"    "189/226"      "15/21" 
  # adj.P.Val_52 adj.P.Val_53 adj.P.Val_55 adj.P.Val_57 adj.P.Val_59 adj.P.Val_62 adj.P.Val_63 adj.P.Val_66 adj.P.Val_67 
  #      "13/19"        "3/8"        "2/2"      "34/59"        "2/6"    "111/148"      "33/66"      "58/79"      "13/28" 
  # adj.P.Val_69 adj.P.Val_75 adj.P.Val_76 adj.P.Val_77 adj.P.Val_78 adj.P.Val_79 adj.P.Val_83 adj.P.Val_84 adj.P.Val_85 
  #        "2/6"    "169/319"       "7/18"      "11/34"        "2/2"    "110/166"    "218/285"      "54/66"        "3/7" 
  # adj.P.Val_86 adj.P.Val_88 adj.P.Val_94 adj.P.Val_95 adj.P.Val_96 adj.P.Val_97 
  #    "123/205"        "2/2"     "90/111"     "74/105"      "35/49"     "50/127" 
  
  ## LPS:
  
  # adj.P.Val_1  adj.P.Val_2  adj.P.Val_3  adj.P.Val_4  adj.P.Val_6  adj.P.Val_7  adj.P.Val_8 adj.P.Val_10 adj.P.Val_11 
  #     "18/58"     "59/204"     "68/264"      "13/38"     "54/308"        "2/5"        "2/6"       "4/10"     "44/137" 
  # adj.P.Val_12 adj.P.Val_13 adj.P.Val_14 adj.P.Val_16 adj.P.Val_21 adj.P.Val_23 adj.P.Val_26 adj.P.Val_30 adj.P.Val_31 
  #        "2/7"       "7/28"      "15/33"        "2/3"      "17/36"       "8/21"       "5/33"     "49/178"      "16/45" 
  # adj.P.Val_32 adj.P.Val_33 adj.P.Val_36 adj.P.Val_44 adj.P.Val_52 adj.P.Val_53 adj.P.Val_54 adj.P.Val_55 adj.P.Val_59 
  #        "2/4"      "50/82"      "13/57"      "29/52"    "166/607"        "2/3"        "2/8"      "20/54"       "5/20" 
  # adj.P.Val_60 adj.P.Val_63 adj.P.Val_65 adj.P.Val_68 adj.P.Val_69 adj.P.Val_74 adj.P.Val_80 adj.P.Val_81 adj.P.Val_86 
  #      "12/60"        "2/7"      "20/32"        "3/6"       "9/38"       "6/27"       "5/11"     "60/244"       "5/21" 
  # adj.P.Val_87 adj.P.Val_88 adj.P.Val_90 adj.P.Val_93 adj.P.Val_95 
  #        "4/6"        "5/9"      "10/35"      "12/36"        "2/3" 
  
  ## untreated:
  
  # adj.P.Val_2  adj.P.Val_3  adj.P.Val_4  adj.P.Val_5  adj.P.Val_8  adj.P.Val_9 adj.P.Val_10 adj.P.Val_13 adj.P.Val_14 
  #       "3/3"        "5/7"      "27/30"        "4/7"      "25/31"     "90/134"    "132/177"        "2/4"     "60/119" 
  # adj.P.Val_15 adj.P.Val_18 adj.P.Val_21 adj.P.Val_23 adj.P.Val_25 adj.P.Val_26 adj.P.Val_27 adj.P.Val_29 adj.P.Val_30 
  #        "2/5"      "11/26"        "2/5"     "63/103"        "7/8"        "6/7"        "2/8"     "74/108"        "2/3" 
  # adj.P.Val_31 adj.P.Val_32 adj.P.Val_34 adj.P.Val_35 adj.P.Val_36 adj.P.Val_37 adj.P.Val_38 adj.P.Val_40 adj.P.Val_43 
  #        "4/6"       "7/13"      "39/42"    "221/431"      "26/39"       "7/14"        "2/5"        "3/3"      "18/27" 
  # adj.P.Val_44 adj.P.Val_45 adj.P.Val_46 adj.P.Val_48 adj.P.Val_51 adj.P.Val_52 adj.P.Val_53 adj.P.Val_54 adj.P.Val_56 
  #    "116/141"      "30/61"    "265/355"      "62/93"     "93/142"        "2/9"      "26/33"       "5/12"        "2/9" 
  # adj.P.Val_57 adj.P.Val_59 adj.P.Val_61 adj.P.Val_62 adj.P.Val_63 adj.P.Val_65 adj.P.Val_68 adj.P.Val_69 adj.P.Val_70 
  #      "74/86"    "287/428"        "2/4"     "62/141"    "230/257"      "24/48"       "8/29"        "3/6"        "2/5" 
  # adj.P.Val_71 adj.P.Val_76 adj.P.Val_77 adj.P.Val_78 adj.P.Val_79 adj.P.Val_82 adj.P.Val_84 adj.P.Val_85 adj.P.Val_86 
  #      "33/66"        "2/2"        "2/5"      "31/62"      "38/51"      "29/43"       "9/15"      "23/45"      "16/30" 
  # adj.P.Val_87 adj.P.Val_88 adj.P.Val_90 adj.P.Val_91 adj.P.Val_92 adj.P.Val_94 adj.P.Val_97 adj.P.Val_98 adj.P.Val_99 
  #       "7/26"      "35/45"     "70/203"       "7/16"    "139/173"       "6/22"        "3/8"        "2/7"       "5/13" 
  
  ## Density plot
  density <- density(num_DEGs, from = min(num_DEGs), to = max(num_DEGs)) 
  df <- data.frame(x = density$x, y = density$y, extreme = sapply(density$x, function(x){if(x>= dim(observed_DEGs)[1]) {"TRUE"} else {"FALSE"} }) ) 
  df$line <- dim(observed_DEGs)[1]
  
  # if(treatment == "IFN"){
  #   hjust = -0.1
  # }
  # else{
  #   hjust = 1.1
  # }
  plot <-  ggplot(data = df, aes(x = x, ymin = 0, ymax = y, fill = extreme, alpha = extreme)) +
    geom_ribbon() +
    theme_bw() +
    scale_fill_manual(values = c("TRUE" = "yellow3", "FALSE" = "seashell2")) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.75)) +
    guides(alpha = "none") +
    geom_line(inherit.aes = F, aes(x = x, y = y), linewidth = 0.2) + 
    labs(x = "Number of DEGs across permutations", y= 'Density', fill='>= observed DEGs',
         subtitle=paste0("Empirical p-val = ", signif(p, digits = 3))) +
    geom_segment(aes(x = line, y = y, xend = line, yend = max(y)), color = 'green3', linewidth = 0.4, data = df) +
    geom_text(aes(x = dim(observed_DEGs)[1], y = max(y), label = paste(dim(observed_DEGs)[1], treatment, "DEGs")), 
              colour = 'green3', hjust = 1.1, vjust = 3, fontface = 1, show.legend = FALSE) +
    coord_cartesian(ylim = c(0.0001, max(df$y))) +
    theme(plot.subtitle = element_text(size = 10, color = "gray30"), 
          axis.text = element_text(size = (10)),
          legend.text = element_text(size = 10),
          legend.title = element_text(size =11),
          axis.title.x = element_text(size = (11.5)),
          axis.title.y = element_text(size = (11.5)))
  
  return(plot)
}

## Empirical pval per treatment
# - - - - - - - - - - - - - - -  IFN p-value  - - - - - - - - - - - - - - - -
p1 <-  treatment_empirical_p("IFN")
# - - - - - - - - - - - - - - -  LPS p-value  - - - - - - - - - - - - - - - 
p2 <- treatment_empirical_p("LPS")
# - - - - - - - - - - - - -  untreated p-value  - - - - - - - - - - - - - - -
p3 <- treatment_empirical_p("untreated")

plot_grid(plotlist = list(p1, p2, p3), ncol = 3, align = "vh")
ggsave(paste0(plotdir, "/treatment_empirical_pvals.pdf"), height = 3.8, width = 17)



## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#                 6.3.2 Compute empirical p-value per gene
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

gene_empirical_p <- function(treatment){
  
  permutation_results <- eval(parse_expr(paste0("permutation_results_", treatment)))
  original_results <- eval(parse_expr(paste0("top_genes_proliferation_", treatment)))
  
  empirical_ps <- vector()
  
  for(gene_id in rownames(original_results)){
    
    ## DGE t-stats across 100 permutations for gene
    gene_permutation_tstats <- permutation_results[gene_id, paste0("t_", 1:100)]
    
    ## Confirm all gene perm t-stats are unique
    stopifnot(length(which(duplicated(gene_permutation_tstats))) == 0)
    
    ## Observed gene t-stat
    gene_observed_tstats <- original_results[gene_id, "t"]
    
    ## Universe of t-stats
    all_tstats <- cbind(gene_permutation_tstats, "t_ob" = gene_observed_tstats)
    
    ## Empirical p = [number of perm |t| >= obs |t| + 1]/[number of permutations + 1]
    empirical_p <- (sum(abs(gene_permutation_tstats) >= abs(gene_observed_tstats)) +1) / length(all_tstats)
    
    empirical_ps[gene_id] <-  empirical_p
  }
  
  return(empirical_ps)
}


## Compute empirical p-vals for all expressed genes in each treatment
empirical_ps_IFN <- gene_empirical_p("IFN")
top_genes_proliferation_IFN$empirical_p <- empirical_ps_IFN
save(top_genes_proliferation_IFN, file = paste0(outdir, "/top_genes_proliferation_IFN_with_empiricalP.Rdata"))

empirical_ps_LPS <- gene_empirical_p("LPS")
top_genes_proliferation_LPS$empirical_p <- empirical_ps_LPS
save(top_genes_proliferation_LPS, file = paste0(outdir, "/top_genes_proliferation_LPS_with_empiricalP.Rdata"))

empirical_ps_untreated <- gene_empirical_p("untreated")
top_genes_proliferation_untreated$empirical_p <- empirical_ps_untreated
save(top_genes_proliferation_untreated, file = paste0(outdir, "/top_genes_proliferation_untreated_with_empiricalP.Rdata"))



## Compare adj p and empirical p across all genes
adjP_vs_empiricalP <- function(treatment){
  
  data <- eval(parse_expr(paste0("top_genes_proliferation_", treatment)))
  
  ## Signif categories based on adj.P and empirical p
  data$signif_cat <- apply(data, 1, function(x){if(as.numeric(x["adj.P.Val"]) < 0.05 & as.numeric(x["empirical_p"]) < 0.05){"DEG; empirical p<0.05"} 
    else if(as.numeric(x["adj.P.Val"]) < 0.05 & as.numeric(x["empirical_p"]) >= 0.05) {"DEG; empirical p>0.05"}
    else if(as.numeric(x["adj.P.Val"]) >= 0.05 & as.numeric(x["empirical_p"]) < 0.05) {"non-DEG; empirical p<0.05"}
    else if(as.numeric(x["adj.P.Val"]) >= 0.05 & as.numeric(x["empirical_p"]) >= 0.05){"non-DEG; empirical p>0.05"} })
  
  labels = c("DEG; empirical p<0.05" = paste0("DEG; empirical p<0.05 (", table(data$signif_cat)["DEG; empirical p<0.05"], ")"),
             "DEG; empirical p>0.05" = paste0("DEG; empirical p>0.05 (", table(data$signif_cat)["DEG; empirical p>0.05"], ")"),
             "non-DEG; empirical p<0.05" = paste0("non-DEG; empirical p<0.05 (", table(data$signif_cat)["non-DEG; empirical p<0.05"], ")"),
             "non-DEG; empirical p>0.05" = paste0("non-DEG; empirical p>0.05 (", table(data$signif_cat)["non-DEG; empirical p>0.05"], ")"))
  
  plot <- ggplot(data, aes(x = adj.P.Val,
                           y = empirical_p,
                           color = signif_cat,
                           alpha = signif_cat)) +
    geom_point(size=1) +
    scale_color_manual(values = c("DEG; empirical p<0.05" = "magenta",
                                  "DEG; empirical p>0.05" = "palegreen3",
                                  "non-DEG; empirical p<0.05" = "lightsteelblue1",
                                  "non-DEG; empirical p>0.05" = "gray80"), label = labels) +
    scale_alpha_manual(values = c("DEG; empirical p<0.05" = 1,
                                  "DEG; empirical p>0.05" = 0.9,
                                  "non-DEG; empirical p<0.05" = 0.9,
                                  "non-DEG; empirical p>0.05" = 0.3)) +
    geom_vline(xintercept = 0.05, colour = "orangered", linewidth = 0.5, linetype = 2) +
    geom_hline(yintercept = 0.05, colour = "orangered", linewidth = 0.5, linetype = 2) +
    theme_bw() +
    guides(alpha = "none") +
    labs(title = treatment, 
         x = "adj.P.Val for true proliferation scaled log-fractions", 
         y = "empirical p-value based on permutations",
         color = "Significance") +
    theme(plot.margin=unit (c (0.4,0.1,0.4,0.1), 'cm'),
          plot.title = element_text(size = (9), face = "bold"),
          axis.title = element_text(size = (8)),
          axis.text = element_text(size = (6)),
          legend.text = element_text(size=7),
          legend.title = element_text(size=8), 
          legend.key.height= unit(0.3, 'cm'))
  
  return(plot)
  
}

p1 <-  adjP_vs_empiricalP("IFN")
p2 <-  adjP_vs_empiricalP("LPS")
p3 <-  adjP_vs_empiricalP("untreated")

plot_grid(plotlist = list(p1, p2, p3), ncol = 3, align = "vh")
ggsave(paste0(plotdir, "/treatment_adjPVal_vs_empiricalP_all_genes.pdf"), height = 3, width = 15)







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
# date     2025-01-20
# rstudio  2024.04.0+735 Chocolate Cosmos (server)
# pandoc   3.1.12.3 @ /opt/view/bin/ (via rmarkdown)
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [1] CRAN (R 4.3.1)
# Biobase              * 2.62.0    2023-10-24 [1] Bioconductor
# BiocGenerics         * 0.48.1    2023-11-01 [1] Bioconductor
# bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.3.1)
# cli                    3.6.2     2023-12-11 [1] CRAN (R 4.3.1)
# colorspace             2.1-0     2023-01-23 [1] CRAN (R 4.3.1)
# cowplot              * 1.1.3     2024-01-22 [1] CRAN (R 4.3.1)
# crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.1)
# DelayedArray           0.28.0    2023-10-24 [1] Bioconductor
# digest                 0.6.35    2024-03-11 [1] CRAN (R 4.3.1)
# dplyr                * 1.1.4     2023-11-17 [1] CRAN (R 4.3.1)
# evaluate               0.23      2023-11-01 [1] CRAN (R 4.3.1)
# fansi                  1.0.6     2023-12-08 [1] CRAN (R 4.3.1)
# farver                 2.1.1     2022-07-06 [1] CRAN (R 4.3.1)
# fastmap                1.1.1     2023-02-24 [1] CRAN (R 4.3.1)
# forcats              * 1.0.0     2023-01-29 [1] CRAN (R 4.3.1)
# generics               0.1.3     2022-07-05 [1] CRAN (R 4.3.1)
# GenomeInfoDb         * 1.38.8    2024-03-15 [1] Bioconductor 3.18 (R 4.3.1)
# GenomeInfoDbData       1.2.11    2024-11-06 [1] Bioconductor
# GenomicRanges        * 1.54.1    2023-10-29 [1] Bioconductor
# ggplot2              * 3.5.1     2024-04-23 [1] CRAN (R 4.3.1)
# glue                   1.7.0     2024-01-09 [1] CRAN (R 4.3.1)
# gtable                 0.3.4     2023-08-21 [1] CRAN (R 4.3.1)
# hms                    1.1.3     2023-03-21 [1] CRAN (R 4.3.1)
# htmltools              0.5.8     2024-03-25 [1] CRAN (R 4.3.1)
# IRanges              * 2.36.0    2023-10-24 [1] Bioconductor
# knitr                  1.45      2023-10-30 [1] CRAN (R 4.3.1)
# labeling               0.4.3     2023-08-29 [1] CRAN (R 4.3.1)
# lattice                0.22-6    2024-03-20 [1] CRAN (R 4.3.1)
# lifecycle              1.0.4     2023-11-07 [1] CRAN (R 4.3.1)
# lubridate            * 1.9.3     2023-09-27 [1] CRAN (R 4.3.1)
# magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.3.1)
# Matrix                 1.6-5     2024-01-11 [1] CRAN (R 4.3.1)
# MatrixGenerics       * 1.14.0    2023-10-24 [1] Bioconductor
# matrixStats          * 1.2.0     2023-12-11 [1] CRAN (R 4.3.1)
# munsell                0.5.1     2024-04-01 [1] CRAN (R 4.3.1)
# pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.3.1)
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.1)
# purrr                * 1.0.2     2023-08-10 [1] CRAN (R 4.3.1)
# R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.1)
# ragg                   1.3.0     2024-03-13 [1] CRAN (R 4.3.1)
# RCurl                  1.98-1.14 2024-01-09 [1] CRAN (R 4.3.1)
# readr                * 2.1.5     2024-01-10 [1] CRAN (R 4.3.1)
# rlang                * 1.1.3     2024-01-10 [1] CRAN (R 4.3.1)
# rmarkdown              2.26      2024-03-05 [1] CRAN (R 4.3.1)
# rstudioapi             0.16.0    2024-03-24 [1] CRAN (R 4.3.1)
# S4Arrays               1.2.1     2024-03-04 [1] Bioconductor 3.18 (R 4.3.1)
# S4Vectors            * 0.40.2    2023-11-23 [1] Bioconductor 3.18 (R 4.3.1)
# scales                 1.3.0     2023-11-28 [1] CRAN (R 4.3.1)
# sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.3.1)
# SparseArray            1.2.4     2024-02-11 [1] Bioconductor 3.18 (R 4.3.1)
# stringi                1.8.3     2023-12-11 [1] CRAN (R 4.3.1)
# stringr              * 1.5.1     2023-11-14 [1] CRAN (R 4.3.1)
# SummarizedExperiment * 1.32.0    2023-10-24 [1] Bioconductor
# systemfonts            1.0.6     2024-03-07 [1] CRAN (R 4.3.1)
# textshaping            0.3.7     2023-10-09 [1] CRAN (R 4.3.1)
# tibble               * 3.2.1     2023-03-20 [1] CRAN (R 4.3.1)
# tidyr                * 1.3.1     2024-01-24 [1] CRAN (R 4.3.1)
# tidyselect             1.2.1     2024-03-11 [1] CRAN (R 4.3.1)
# tidyverse            * 2.0.0     2023-02-22 [1] CRAN (R 4.3.1)
# timechange             0.3.0     2024-01-18 [1] CRAN (R 4.3.1)
# tzdb                   0.4.0     2023-05-12 [1] CRAN (R 4.3.1)
# utf8                   1.2.4     2023-10-22 [1] CRAN (R 4.3.1)
# vctrs                  0.6.5     2023-12-01 [1] CRAN (R 4.3.1)
# withr                  3.0.0     2024-01-16 [1] CRAN (R 4.3.1)
# xfun                   0.43      2024-03-25 [1] CRAN (R 4.3.1)
# XVector                0.42.0    2023-10-24 [1] Bioconductor
# yaml                   2.3.8     2023-12-11 [1] CRAN (R 4.3.1)
# zlibbioc               1.48.2    2024-03-13 [1] Bioconductor 3.18 (R 4.3.1)
# 
# [1] /opt/view/rlib/R/library
# [2] /software/hgi/installs/softpack/rstudio/.spack/linux-ubuntu22.04-x86_64_v3/gcc-11.4.0/r-4.3.1-bfwldrk76z6f52upk47zepliekn7ayqz/rlib/R/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
