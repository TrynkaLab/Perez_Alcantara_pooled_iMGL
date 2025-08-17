
library(tidyverse)
library(dplyr)
library(SummarizedExperiment)
library(VennDiagram) 
library(rlang)
library(sessioninfo)


#-------------------------------------------------------------------------------
#     7.4 Exploration of overlapping DEGs across treatments and phenotypes
#-------------------------------------------------------------------------------
#  Code to get the DEGs common between treatments and phenotypes.
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

## Set working dir
setwd("/lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_Daianna/")

## Define output data and plot dir
outdir = paste(getwd(), "output_data", "07_Diff_expr_proliferation_permutations", "04_DEGs_overlaps_explorations", sep = "/")
plotdir = paste(getwd(), "plots", "07_Diff_expr_proliferation_permutations", "04_DEGs_overlaps_explorations", sep = "/")
dir.create(outdir, recursive = T)
dir.create(plotdir, recursive = T)

## Input dirs
input_dir01 = paste(getwd(), "output_data", "01_Diff_expr_PRS_permutations", "00_Data_exploration_processing", sep = "/")
input_dir05 = paste(getwd(), "output_data", "05_Diff_expr_phagocytosis_permutations", "01_run_phagocytosis_DGE", sep = "/")
input_dir06 = paste(getwd(), "output_data", "06_Diff_expr_migration_permutations", "01_run_migration_DGE", sep = "/")
input_dir07 = paste(getwd(), "output_data", "07_Diff_expr_proliferation_permutations", "01_run_proliferation_DGE", sep = "/")

## color dict
load(paste0(input_dir01, "/level_colors.Rdata"), verbose = T)

## Load DEGs for phagocytosis
load(paste0(input_dir05, "/top_genes_phago_IFN.Rdata"), verbose = T)
load(paste0(input_dir05, "/de_genes_phago_IFN.Rdata"), verbose = T)

load(paste0(input_dir05, "/top_genes_phago_LPS.Rdata"), verbose = T)
load(paste0(input_dir05, "/de_genes_phago_LPS.Rdata"), verbose = T)

load(paste0(input_dir05, "/top_genes_phago_untreated.Rdata"), verbose = T)
load(paste0(input_dir05, "/de_genes_phago_untreated.Rdata"), verbose = T)

## Load DEGs for migration
load(paste0(input_dir06, "/top_genes_migration_IFN.Rdata"), verbose = T)
load(paste0(input_dir06, "/de_genes_migration_IFN.Rdata"), verbose = T)

load(paste0(input_dir06, "/top_genes_migration_LPS.Rdata"), verbose = T)
load(paste0(input_dir06, "/de_genes_migration_LPS.Rdata"), verbose = T)

load(paste0(input_dir06, "/top_genes_migration_untreated.Rdata"), verbose = T)
load(paste0(input_dir06, "/de_genes_migration_untreated.Rdata"), verbose = T)

## Load DEGs for proliferation
load(paste0(input_dir07, "/top_genes_proliferation_IFN.Rdata"), verbose = T)
load(paste0(input_dir07, "/de_genes_proliferation_IFN.Rdata"), verbose = T)

load(paste0(input_dir07, "/top_genes_proliferation_LPS.Rdata"), verbose = T)
load(paste0(input_dir07, "/de_genes_proliferation_LPS.Rdata"), verbose = T)

load(paste0(input_dir07, "/top_genes_proliferation_untreated.Rdata"), verbose = T)
load(paste0(input_dir07, "/de_genes_proliferation_untreated.Rdata"), verbose = T)


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#         7.4.1 Comparison of DEGs for same phenotype across treatments 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

############################## Phagocytosis DEGs ############################### 
# * IFN:
DEGs_phago_IFN <- na.omit(de_genes_phago_IFN$gene)
# * LPS:
DEGs_phago_LPS <- na.omit(de_genes_phago_LPS$gene)
# * untreated:
DEGs_phago_untreated <- na.omit(de_genes_phago_untreated$gene)

DEGs_list <- list("Phagocytosis DEGs in IFN" = DEGs_phago_IFN,
                  "Phagocytosis DEGs in LPS" = DEGs_phago_LPS,
                  "Phagocytosis DEGs in untreated" = DEGs_phago_untreated)

colors <- level_colors$treatment 
names(colors) <-  names(DEGs_list)

v_phago <- venn.diagram(DEGs_list, fill = colors, alpha = rep(0.5, length(DEGs_list)), 
                lwd = 0, margin = 0.3, cat.cex = 0.9, cex = 1, height = 20, width = 50, units = "cm", 
                cat.dist = rep(0.15, length(DEGs_list)), filename = NULL, disable.logging = T)


############################# Migration DEGs ############################### 
# * IFN:
DEGs_mig_IFN <- na.omit(de_genes_migration_IFN$gene)
# * LPS:
DEGs_mig_LPS <- na.omit(de_genes_migration_LPS$gene)
# * untreated:
DEGs_mig_untreated <- na.omit(de_genes_migration_untreated$gene)

DEGs_list <- list("Migration DEGs in IFN" = DEGs_mig_IFN,
                  "Migration DEGs in LPS" = DEGs_mig_LPS,
                  "Migration DEGs in untreated" = DEGs_mig_untreated)

colors <- level_colors$treatment 
names(colors) <-  names(DEGs_list)

v_mig <- venn.diagram(DEGs_list, fill = colors, alpha = rep(0.5, length(DEGs_list)), 
                        lwd = 0, margin = 0.3, cat.cex = 0.9, cex = 1, height = 20, width = 50, units = "cm", 
                        cat.dist = rep(0.15, length(DEGs_list)), filename = NULL, disable.logging = T)


############################# Proliferation DEGs ############################### 
# * IFN:
DEGs_prolif_IFN <- na.omit(de_genes_proliferation_IFN$gene)
# * LPS:
DEGs_prolif_LPS <- na.omit(de_genes_proliferation_LPS$gene)
# * untreated:
DEGs_prolif_untreated <- na.omit(de_genes_proliferation_untreated$gene)

DEGs_list <- list("Proliferation DEGs in IFN" = DEGs_prolif_IFN,
                  "Proliferation DEGs in LPS" = DEGs_prolif_LPS,
                  "Proliferation DEGs in untreated" = DEGs_prolif_untreated)

colors <- level_colors$treatment 
names(colors) <-  names(DEGs_list)

v_prolif <- venn.diagram(DEGs_list, fill = colors, alpha = rep(0.5, length(DEGs_list)), 
                      lwd = 0, margin = 0.3, cat.cex = 0.9, cex = 1, height = 20, width = 50, units = "cm", 
                      cat.dist = rep(0.15, length(DEGs_list)), filename = NULL, disable.logging = T)

plot_grid(plotlist = list(v_phago, v_mig, v_prolif), ncol = 3, align = "h")
ggsave(filename = paste0(plotdir, "/Venn_diagrams_compare_treatments.pdf"), width = 12, height = 4)



## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#         7.4.2 Comparison of DEGs for same treatment across phenotypes 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

################################## IFN DEGs #################################### 
# * Phagocytosis:
DEGs_phago_IFN <- na.omit(de_genes_phago_IFN$gene)
# * Migration:
DEGs_mig_IFN <- na.omit(de_genes_migration_IFN$gene)
# * Proliferation:
DEGs_prolif_IFN <- na.omit(de_genes_proliferation_IFN$gene)

DEGs_list <- list("Phagocytosis DEGs in IFN" = DEGs_phago_IFN,
                  "Migration DEGs in IFN" = DEGs_mig_IFN,
                  "Proliferation DEGs in IFN" = DEGs_prolif_IFN)

colors <- c("darkmagenta", "mediumturquoise", "royalblue")
names(colors) <-  names(DEGs_list)

v_IFN <- venn.diagram(DEGs_list, fill = colors, alpha = rep(0.5, length(DEGs_list)), 
                        lwd = 0, margin = 0.3, cat.cex = 0.9, cex = 1, height = 20, width = 50, units = "cm", 
                        cat.dist = rep(0.15, length(DEGs_list)), filename = NULL, disable.logging = T)


################################## LPS DEGs #################################### 
# * Phagocytosis:
DEGs_phago_LPS <- na.omit(de_genes_phago_LPS$gene)
# * Migration:
DEGs_mig_LPS <- na.omit(de_genes_migration_LPS$gene)
# * Proliferation:
DEGs_prolif_LPS <- na.omit(de_genes_proliferation_LPS$gene)

DEGs_list <- list("Phagocytosis DEGs in LPS" = DEGs_phago_LPS,
                  "Migration DEGs in LPS" = DEGs_mig_LPS,
                  "Proliferation DEGs in LPS" = DEGs_prolif_LPS)

names(colors) <-  names(DEGs_list)

v_LPS <- venn.diagram(DEGs_list, fill = colors, alpha = rep(0.5, length(DEGs_list)), 
                      lwd = 0, margin = 0.3, cat.cex = 0.9, cex = 1, height = 20, width = 50, units = "cm", 
                      cat.dist = rep(0.15, length(DEGs_list)), filename = NULL, disable.logging = T)


################################## untreated DEGs #################################### 
# * Phagocytosis:
DEGs_phago_untreated <- na.omit(de_genes_phago_untreated$gene)
# * Migration:
DEGs_mig_untreated <- na.omit(de_genes_migration_untreated$gene)
# * Proliferation:
DEGs_prolif_untreated <- na.omit(de_genes_proliferation_untreated$gene)

DEGs_list <- list("Phagocytosis DEGs in untreated" = DEGs_phago_untreated,
                  "Migration DEGs in untreated" = DEGs_mig_untreated,
                  "Proliferation DEGs in untreated" = DEGs_prolif_untreated)

names(colors) <-  names(DEGs_list)

v_untreated <- venn.diagram(DEGs_list, fill = colors, alpha = rep(0.5, length(DEGs_list)), 
                      lwd = 0, margin = 0.3, cat.cex = 0.9, cex = 1, height = 20, width = 50, units = "cm", 
                      cat.dist = rep(0.15, length(DEGs_list)), filename = NULL, disable.logging = T)

plot_grid(plotlist = list(v_IFN, v_LPS, v_untreated),  ncol = 3, align = "h")
ggsave(filename = paste0(plotdir, "/Venn_diagrams_compare_phenotypes.pdf"), width = 13, height = 4)



############################### IFN-only DEGs ##################################  
## DEGs for each phenotype in IFN only
# * Phagocytosis:
DEGs_phago_IFN_only <- DEGs_phago_IFN[-which(DEGs_phago_IFN %in% c(DEGs_phago_LPS, DEGs_phago_untreated))]
# * Migration:
DEGs_mig_IFN_only <- DEGs_mig_IFN[-which(DEGs_mig_IFN %in% c(DEGs_mig_LPS, DEGs_mig_untreated))]
# * Proliferation:
DEGs_prolif_IFN_only <- DEGs_prolif_IFN[-which(DEGs_prolif_IFN %in% c(DEGs_prolif_LPS, DEGs_prolif_untreated))]

DEGs_list <- list("Phagocytosis DEGs in IFN only" = DEGs_phago_IFN_only,
                  "Migration DEGs in IFN only" = DEGs_mig_IFN_only,
                  "Proliferation DEGs in IFN only" = DEGs_prolif_IFN_only)

names(colors) <-  names(DEGs_list)

v_IFN_only <- venn.diagram(DEGs_list, fill = colors, alpha = rep(0.5, length(DEGs_list)), 
                            lwd = 0, margin = 0.3, cat.cex = 0.9, cex = 1, height = 20, width = 50, units = "cm", 
                            cat.dist = rep(0.15, length(DEGs_list)), filename = NULL, disable.logging = T)


############################### LPS-only DEGs ##################################  
## DEGs for each phenotype in LPS only
# * Phagocytosis:
DEGs_phago_LPS_only <- DEGs_phago_LPS[-which(DEGs_phago_LPS %in% c(DEGs_phago_IFN, DEGs_phago_untreated))]
# * Migration:
DEGs_mig_LPS_only <- DEGs_mig_LPS[-which(DEGs_mig_LPS %in% c(DEGs_mig_IFN, DEGs_mig_untreated))]
# * Proliferation:
DEGs_prolif_LPS_only <- DEGs_prolif_LPS[-which(DEGs_prolif_LPS %in% c(DEGs_prolif_IFN, DEGs_prolif_untreated))]

DEGs_list <- list("Phagocytosis DEGs in LPS only" = DEGs_phago_LPS_only,
                  "Migration DEGs in LPS only" = DEGs_mig_LPS_only,
                  "Proliferation DEGs in LPS only" = DEGs_prolif_LPS_only)

names(colors) <-  names(DEGs_list)

v_LPS_only <- venn.diagram(DEGs_list, fill = colors, alpha = rep(0.5, length(DEGs_list)), 
                           lwd = 0, margin = 0.3, cat.cex = 0.9, cex = 1, height = 20, width = 50, units = "cm", 
                           cat.dist = rep(0.15, length(DEGs_list)), filename = NULL, disable.logging = T)


############################ (IFN | LPS)-only DEGs #############################  
## DEGs for each phenotype in IFN | LPS only
# * Phagocytosis:
DEGs_phago_IFN_or_LPS <- union(DEGs_phago_IFN, DEGs_phago_LPS) 
DEGs_phago_IFN_or_LPS <- DEGs_phago_IFN_or_LPS[-which(DEGs_phago_IFN_or_LPS %in% DEGs_phago_untreated)]

# * Migration:
DEGs_mig_IFN_or_LPS <- union(DEGs_mig_IFN, DEGs_mig_LPS) 
DEGs_mig_IFN_or_LPS <- DEGs_mig_IFN_or_LPS[-which(DEGs_mig_IFN_or_LPS %in% DEGs_mig_untreated)]

# * Proliferation:
DEGs_prolif_IFN_or_LPS <- union(DEGs_prolif_IFN, DEGs_prolif_LPS) 
DEGs_prolif_IFN_or_LPS <- DEGs_prolif_IFN_or_LPS[-which(DEGs_prolif_IFN_or_LPS %in% DEGs_prolif_untreated)]

DEGs_list <- list("Phagocytosis DEGs in (IFN | LPS) only" = DEGs_phago_IFN_or_LPS,
                  "Migration DEGs in (IFN | LPS) only" = DEGs_mig_IFN_or_LPS,
                  "Proliferation DEGs in (IFN | LPS) only" = DEGs_prolif_IFN_or_LPS)


names(colors) <-  names(DEGs_list)

v_IFN_or_LPS_only <- venn.diagram(DEGs_list, fill = colors, alpha = rep(0.5, length(DEGs_list)), 
                           lwd = 0, margin = 0.3, cat.cex = 0.9, cex = 1, height = 20, width = 50, units = "cm", 
                           cat.dist = rep(0.15, length(DEGs_list)), filename = NULL, disable.logging = T)


############################ (IFN & LPS)-only DEGs #############################  
## DEGs for each phenotype in IFN & LPS only
# * Phagocytosis:
DEGs_phago_IFN_and_LPS <- intersect(DEGs_phago_IFN, DEGs_phago_LPS) 
DEGs_phago_IFN_and_LPS <- DEGs_phago_IFN_and_LPS[-which(DEGs_phago_IFN_and_LPS %in% DEGs_phago_untreated)]

# * Migration:
DEGs_mig_IFN_and_LPS <- intersect(DEGs_mig_IFN, DEGs_mig_LPS) 
DEGs_mig_IFN_and_LPS <- DEGs_mig_IFN_and_LPS[-which(DEGs_mig_IFN_and_LPS %in% DEGs_mig_untreated)]

# * Proliferation:
DEGs_prolif_IFN_and_LPS <- intersect(DEGs_prolif_IFN, DEGs_prolif_LPS) 
DEGs_prolif_IFN_and_LPS <- DEGs_prolif_IFN_and_LPS[-which(DEGs_prolif_IFN_and_LPS %in% DEGs_prolif_untreated)]

DEGs_list <- list("Phagocytosis DEGs in (IFN & LPS) only" = DEGs_phago_IFN_and_LPS,
                  "Migration DEGs in (IFN & LPS) only" = DEGs_mig_IFN_and_LPS,
                  "Proliferation DEGs in (IFN & LPS) only" = DEGs_prolif_IFN_and_LPS)

 
names(colors) <-  names(DEGs_list)

v_IFN_and_LPS_only <- venn.diagram(DEGs_list, fill = colors, alpha = rep(0.5, length(DEGs_list)), 
                               lwd = 0, margin = 0.3, cat.cex = 0.9, cex = 1, height = 20, width = 50, units = "cm", 
                               cat.dist = rep(0.15, length(DEGs_list)), filename = NULL, disable.logging = T)


plot_grid(plotlist = list(v_IFN_only, v_LPS_only, v_IFN_or_LPS_only, v_IFN_and_LPS_only), ncol = 4, align = "h")
ggsave(filename = paste0(plotdir, "/Venn_diagrams_compare_phenotypes_DEGs_in_treatment_only.pdf"), width = 21, height = 5.4)



## Save DEGs in IFN-only, LPS-only, or in both, per phenotype
# * Phagocytosis:
top_genes_phago_treatments <- merge(top_genes_phago_IFN, top_genes_phago_LPS, all = T,
                                    by = colnames(top_genes_phago_IFN)[1:11], 
                                    suffixes = c("_IFN", "_LPS"))
rownames(top_genes_phago_treatments) <- top_genes_phago_treatments$gene

DEGs_treatment_phagocytosis <- top_genes_phago_treatments[DEGs_phago_IFN_or_LPS,]
DEGs_treatment_phagocytosis$DEG_in_IFN <- sapply(DEGs_treatment_phagocytosis$adj.P.Val_IFN, 
                                                 function(x){if(!is.na(x) & x < 0.05){TRUE} else{FALSE}})
DEGs_treatment_phagocytosis$DEG_in_LPS <- sapply(DEGs_treatment_phagocytosis$adj.P.Val_LPS, 
                                                 function(x){if(!is.na(x) &  x<0.05){TRUE} else{FALSE}})
DEGs_treatment_phagocytosis$DGE_across_treatments <- apply(DEGs_treatment_phagocytosis, 1, 
                                                            function(x){if(x["DEG_in_IFN"] == TRUE & x["DEG_in_LPS"] == FALSE){"DEG_in_IFN_only"}
                                                                        else if(x["DEG_in_IFN"] == FALSE & x["DEG_in_LPS"] == TRUE){"DEG_in_LPS_only"}
                                                                        else if(x["DEG_in_IFN"] == TRUE & x["DEG_in_LPS"] == TRUE){"DEG_in_IFN_and_LPS_only"}} ) 
save(DEGs_treatment_phagocytosis, file = paste0(outdir, "/DEGs_in_IFN_or_LPS_phagocytosis.Rdata"))

# * Migration:
top_genes_mig_treatments <- merge(top_genes_migration_IFN, top_genes_migration_LPS, all = T,
                                    by = colnames(top_genes_migration_IFN)[1:11], 
                                    suffixes = c("_IFN", "_LPS"))
rownames(top_genes_mig_treatments) <- top_genes_mig_treatments$gene

DEGs_treatment_migration <- top_genes_mig_treatments[DEGs_mig_IFN_or_LPS,]
DEGs_treatment_migration$DEG_in_IFN <- sapply(DEGs_treatment_migration$adj.P.Val_IFN, 
                                                 function(x){if(!is.na(x) & x < 0.05){TRUE} else{FALSE}})
DEGs_treatment_migration$DEG_in_LPS <- sapply(DEGs_treatment_migration$adj.P.Val_LPS, 
                                                 function(x){if(!is.na(x) &  x<0.05){TRUE} else{FALSE}})
DEGs_treatment_migration$DGE_across_treatments <- apply(DEGs_treatment_migration, 1, 
                                                           function(x){if(x["DEG_in_IFN"] == TRUE & x["DEG_in_LPS"] == FALSE){"DEG_in_IFN_only"}
                                                             else if(x["DEG_in_IFN"] == FALSE & x["DEG_in_LPS"] == TRUE){"DEG_in_LPS_only"}
                                                             else if(x["DEG_in_IFN"] == TRUE & x["DEG_in_LPS"] == TRUE){"DEG_in_IFN_and_LPS_only"}} ) 
save(DEGs_treatment_migration, file = paste0(outdir, "/DEGs_in_IFN_or_LPS_migration.Rdata"))

# * Proliferation:
top_genes_prolif_treatments <- merge(top_genes_proliferation_IFN, top_genes_proliferation_LPS, all = T,
                                  by = colnames(top_genes_proliferation_IFN)[1:11], 
                                  suffixes = c("_IFN", "_LPS"))
rownames(top_genes_prolif_treatments) <- top_genes_prolif_treatments$gene

DEGs_treatment_proliferation <- top_genes_prolif_treatments[DEGs_prolif_IFN_or_LPS,]
DEGs_treatment_proliferation$DEG_in_IFN <- sapply(DEGs_treatment_proliferation$adj.P.Val_IFN, 
                                              function(x){if(!is.na(x) & x < 0.05){TRUE} else{FALSE}})
DEGs_treatment_proliferation$DEG_in_LPS <- sapply(DEGs_treatment_proliferation$adj.P.Val_LPS, 
                                              function(x){if(!is.na(x) &  x<0.05){TRUE} else{FALSE}})
DEGs_treatment_proliferation$DGE_across_treatments <- apply(DEGs_treatment_proliferation, 1, 
                                                        function(x){if(x["DEG_in_IFN"] == TRUE & x["DEG_in_LPS"] == FALSE){"DEG_in_IFN_only"}
                                                          else if(x["DEG_in_IFN"] == FALSE & x["DEG_in_LPS"] == TRUE){"DEG_in_LPS_only"}
                                                          else if(x["DEG_in_IFN"] == TRUE & x["DEG_in_LPS"] == TRUE){"DEG_in_IFN_and_LPS_only"}} ) 
save(DEGs_treatment_proliferation, file = paste0(outdir, "/DEGs_in_IFN_or_LPS_proliferation.Rdata"))


## Aggregate across phenotype
DEGs_in_treatments_only_all_phenotypes <- rbind(cbind("gene" = DEGs_treatment_phagocytosis$gene, 
                                                      "phenotype" = "phagocytosis", 
                                                      "DGE_across_treatments" = DEGs_treatment_phagocytosis$DGE_across_treatments),
                                                cbind("gene" = DEGs_treatment_migration$gene, 
                                                      "phenotype" = "migration", 
                                                      "DGE_across_treatments" = DEGs_treatment_migration$DGE_across_treatments),
                                                cbind("gene" = DEGs_treatment_proliferation$gene, 
                                                      "phenotype" = "proliferation", 
                                                      "DGE_across_treatments" = DEGs_treatment_proliferation$DGE_across_treatments)) 
write.csv(DEGs_in_treatments_only_all_phenotypes, file = paste0(outdir, "/DEGs_in_treatments_only_all_phenotypes.csv"), 
          row.names = F, quote = F, sep = "\t")







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
# date     2025-01-24
# rstudio  2024.04.0+735 Chocolate Cosmos (server)
# pandoc   3.1.12.3 @ /opt/view/bin/ (via rmarkdown)
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version    date (UTC) lib source
# abind                  1.4-5      2016-07-21 [1] CRAN (R 4.3.1)
# aod                    1.3.3      2023-12-13 [1] CRAN (R 4.3.1)
# backports              1.4.1      2021-12-13 [1] CRAN (R 4.3.1)
# Biobase              * 2.62.0     2023-10-24 [1] Bioconductor
# BiocGenerics         * 0.48.1     2023-11-01 [1] Bioconductor
# BiocParallel         * 1.36.0     2023-10-24 [1] Bioconductor
# bitops                 1.0-7      2021-04-24 [1] CRAN (R 4.3.1)
# boot                   1.3-30     2024-02-26 [1] CRAN (R 4.3.1)
# broom                  1.0.5      2023-06-09 [1] CRAN (R 4.3.1)
# caTools                1.18.2     2021-03-28 [1] CRAN (R 4.3.1)
# circlize               0.4.16     2024-02-20 [1] CRAN (R 4.3.1)
# cli                    3.6.2      2023-12-11 [1] CRAN (R 4.3.1)
# clue                   0.3-65     2023-09-23 [1] CRAN (R 4.3.1)
# cluster                2.1.6      2023-12-01 [1] CRAN (R 4.3.1)
# codetools              0.2-20     2024-03-31 [1] CRAN (R 4.3.1)
# colorspace             2.1-0      2023-01-23 [1] CRAN (R 4.3.1)
# ComplexHeatmap       * 2.18.0     2023-10-24 [1] Bioconductor
# corpcor                1.6.10     2021-09-16 [1] CRAN (R 4.3.1)
# cowplot              * 1.1.3      2024-01-22 [1] CRAN (R 4.3.1)
# crayon                 1.5.2      2022-09-29 [1] CRAN (R 4.3.1)
# DelayedArray           0.28.0     2023-10-24 [1] Bioconductor
# digest                 0.6.35     2024-03-11 [1] CRAN (R 4.3.1)
# doParallel             1.0.17     2022-02-07 [1] CRAN (R 4.3.1)
# dplyr                * 1.1.4      2023-11-17 [1] CRAN (R 4.3.1)
# edgeR                * 4.0.16     2024-02-18 [1] Bioconductor 3.18 (R 4.3.1)
# EnvStats               2.8.1      2023-08-22 [1] CRAN (R 4.3.1)
# evaluate               0.23       2023-11-01 [1] CRAN (R 4.3.1)
# fANCOVA                0.6-1      2020-11-13 [1] CRAN (R 4.3.1)
# fansi                  1.0.6      2023-12-08 [1] CRAN (R 4.3.1)
# fastmap                1.1.1      2023-02-24 [1] CRAN (R 4.3.1)
# forcats              * 1.0.0      2023-01-29 [1] CRAN (R 4.3.1)
# foreach                1.5.2      2022-02-02 [1] CRAN (R 4.3.1)
# formatR                1.14       2023-01-17 [1] CRAN (R 4.3.1)
# futile.logger        * 1.4.3      2016-07-10 [1] CRAN (R 4.3.1)
# futile.options         1.0.1      2018-04-20 [1] CRAN (R 4.3.1)
# generics               0.1.3      2022-07-05 [1] CRAN (R 4.3.1)
# GenomeInfoDb         * 1.38.8     2024-03-15 [1] Bioconductor 3.18 (R 4.3.1)
# GenomeInfoDbData       1.2.11     2024-11-06 [1] Bioconductor
# GenomicRanges        * 1.54.1     2023-10-29 [1] Bioconductor
# GetoptLong             1.0.5      2020-12-15 [1] CRAN (R 4.3.1)
# ggplot2              * 3.5.1      2024-04-23 [1] CRAN (R 4.3.1)
# ggrepel                0.9.5      2024-01-10 [1] CRAN (R 4.3.1)
# GlobalOptions          0.1.2      2020-06-10 [1] CRAN (R 4.3.1)
# glue                   1.7.0      2024-01-09 [1] CRAN (R 4.3.1)
# gplots                 3.1.3.1    2024-02-02 [1] CRAN (R 4.3.1)
# gtable                 0.3.4      2023-08-21 [1] CRAN (R 4.3.1)
# gtools                 3.9.5      2023-11-20 [1] CRAN (R 4.3.1)
# hms                    1.1.3      2023-03-21 [1] CRAN (R 4.3.1)
# htmltools              0.5.8      2024-03-25 [1] CRAN (R 4.3.1)
# IRanges              * 2.36.0     2023-10-24 [1] Bioconductor
# iterators              1.0.14     2022-02-05 [1] CRAN (R 4.3.1)
# jtools                 2.2.2      2023-07-11 [1] CRAN (R 4.3.1)
# KernSmooth             2.23-22    2023-07-10 [1] CRAN (R 4.3.1)
# knitr                  1.45       2023-10-30 [1] CRAN (R 4.3.1)
# lambda.r               1.2.4      2019-09-18 [1] CRAN (R 4.3.1)
# lattice                0.22-6     2024-03-20 [1] CRAN (R 4.3.1)
# lifecycle              1.0.4      2023-11-07 [1] CRAN (R 4.3.1)
# limma                * 3.58.1     2023-10-31 [1] Bioconductor
# lme4                   1.1-35.2   2024-03-28 [1] CRAN (R 4.3.1)
# lmerTest               3.1-3      2020-10-23 [1] CRAN (R 4.3.1)
# locfit                 1.5-9.9    2024-03-01 [1] CRAN (R 4.3.1)
# lubridate            * 1.9.3      2023-09-27 [1] CRAN (R 4.3.1)
# magrittr               2.0.3      2022-03-30 [1] CRAN (R 4.3.1)
# MASS                   7.3-60.0.1 2024-01-13 [1] CRAN (R 4.3.1)
# Matrix                 1.6-5      2024-01-11 [1] CRAN (R 4.3.1)
# MatrixGenerics       * 1.14.0     2023-10-24 [1] Bioconductor
# matrixStats          * 1.2.0      2023-12-11 [1] CRAN (R 4.3.1)
# minqa                  1.2.6      2023-09-11 [1] CRAN (R 4.3.1)
# munsell                0.5.1      2024-04-01 [1] CRAN (R 4.3.1)
# mvtnorm                1.2-4      2023-11-27 [1] CRAN (R 4.3.1)
# nlme                   3.1-164    2023-11-27 [1] CRAN (R 4.3.1)
# nloptr                 2.0.3      2022-05-26 [1] CRAN (R 4.3.1)
# numDeriv               2016.8-1.1 2019-06-06 [1] CRAN (R 4.3.1)
# pander                 0.6.5      2022-03-18 [1] CRAN (R 4.3.1)
# pbkrtest               0.5.2      2023-01-19 [1] CRAN (R 4.3.1)
# pillar                 1.9.0      2023-03-22 [1] CRAN (R 4.3.1)
# pkgconfig              2.0.3      2019-09-22 [1] CRAN (R 4.3.1)
# plyr                   1.8.9      2023-10-02 [1] CRAN (R 4.3.1)
# png                    0.1-8      2022-11-29 [1] CRAN (R 4.3.1)
# purrr                * 1.0.2      2023-08-10 [1] CRAN (R 4.3.1)
# R6                     2.5.1      2021-08-19 [1] CRAN (R 4.3.1)
# ragg                   1.3.0      2024-03-13 [1] CRAN (R 4.3.1)
# rbibutils              2.2.16     2023-10-25 [1] CRAN (R 4.3.1)
# RColorBrewer           1.1-3      2022-04-03 [1] CRAN (R 4.3.1)
# Rcpp                   1.0.12     2024-01-09 [1] CRAN (R 4.3.1)
# RCurl                  1.98-1.14  2024-01-09 [1] CRAN (R 4.3.1)
# Rdpack                 2.6        2023-11-08 [1] CRAN (R 4.3.1)
# readr                * 2.1.5      2024-01-10 [1] CRAN (R 4.3.1)
# remaCor                0.0.18     2024-02-08 [1] CRAN (R 4.3.1)
# reshape2               1.4.4      2020-04-09 [1] CRAN (R 4.3.1)
# RhpcBLASctl            0.23-42    2023-02-11 [1] CRAN (R 4.3.1)
# rjson                  0.2.21     2022-01-09 [1] CRAN (R 4.3.1)
# rlang                * 1.1.3      2024-01-10 [1] CRAN (R 4.3.1)
# rmarkdown              2.26       2024-03-05 [1] CRAN (R 4.3.1)
# rstudioapi             0.16.0     2024-03-24 [1] CRAN (R 4.3.1)
# S4Arrays               1.2.1      2024-03-04 [1] Bioconductor 3.18 (R 4.3.1)
# S4Vectors            * 0.40.2     2023-11-23 [1] Bioconductor 3.18 (R 4.3.1)
# scales                 1.3.0      2023-11-28 [1] CRAN (R 4.3.1)
# sessioninfo          * 1.2.2      2021-12-06 [1] CRAN (R 4.3.1)
# shape                  1.4.6.1    2024-02-23 [1] CRAN (R 4.3.1)
# SparseArray            1.2.4      2024-02-11 [1] Bioconductor 3.18 (R 4.3.1)
# statmod                1.5.0      2023-01-06 [1] CRAN (R 4.3.1)
# stringi                1.8.3      2023-12-11 [1] CRAN (R 4.3.1)
# stringr              * 1.5.1      2023-11-14 [1] CRAN (R 4.3.1)
# SummarizedExperiment * 1.32.0     2023-10-24 [1] Bioconductor
# systemfonts            1.0.6      2024-03-07 [1] CRAN (R 4.3.1)
# textshaping            0.3.7      2023-10-09 [1] CRAN (R 4.3.1)
# tibble               * 3.2.1      2023-03-20 [1] CRAN (R 4.3.1)
# tidyr                * 1.3.1      2024-01-24 [1] CRAN (R 4.3.1)
# tidyselect             1.2.1      2024-03-11 [1] CRAN (R 4.3.1)
# tidyverse            * 2.0.0      2023-02-22 [1] CRAN (R 4.3.1)
# timechange             0.3.0      2024-01-18 [1] CRAN (R 4.3.1)
# tzdb                   0.4.0      2023-05-12 [1] CRAN (R 4.3.1)
# utf8                   1.2.4      2023-10-22 [1] CRAN (R 4.3.1)
# variancePartition    * 1.32.5     2024-02-16 [1] Bioconductor 3.18 (R 4.3.1)
# vctrs                  0.6.5      2023-12-01 [1] CRAN (R 4.3.1)
# VennDiagram          * 1.7.3      2022-04-12 [1] CRAN (R 4.3.1)
# withr                  3.0.0      2024-01-16 [1] CRAN (R 4.3.1)
# xfun                   0.43       2024-03-25 [1] CRAN (R 4.3.1)
# XVector                0.42.0     2023-10-24 [1] Bioconductor
# yaml                   2.3.8      2023-12-11 [1] CRAN (R 4.3.1)
# zlibbioc               1.48.2     2024-03-13 [1] Bioconductor 3.18 (R 4.3.1)
# 
# [1] /opt/view/rlib/R/library
# [2] /software/hgi/installs/softpack/rstudio/.spack/linux-ubuntu22.04-x86_64_v3/gcc-11.4.0/r-4.3.1-bfwldrk76z6f52upk47zepliekn7ayqz/rlib/R/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
