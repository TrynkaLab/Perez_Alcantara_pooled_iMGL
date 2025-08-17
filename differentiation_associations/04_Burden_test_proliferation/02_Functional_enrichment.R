
library(tidyverse)
library(rlang)
library(org.Hs.eg.db)
library(biomaRt)
library(grid)
library(circlize)
library(ComplexHeatmap)
library(fgsea)
library(cowplot)
library(sessioninfo)

#-------------------------------------------------------------------------------
#    4.2 Functional enrichment analysis for proliferation-associated genes
#-------------------------------------------------------------------------------
#  Code to perform enrichment analysis for genes whose variants associate with
#  proliferation efficiency from iPSC -> macrophage precursor -> microglia, and
#  genes implicated in macrophage survival or essential genes.
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

## Set working dir
setwd("/lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_Daianna/")

## Define output data and plot dir 
outdir = paste(getwd(), "output_data", "04_Burden_test_proliferation", "02_Functional_enrichment", sep = "/")
plotdir = paste(getwd(), "plots", "04_Burden_test_proliferation", "02_Functional_enrichment", sep = "/")
dir.create(outdir, recursive = T)
dir.create(plotdir, recursive = T)

## Subdir for sub analyses
sapply(paste0(outdir, c("/01_ORA", "/02_GSEA")), function(p){dir.create(p, recursive = T)})
sapply(paste0(plotdir, c("/01_ORA", "/02_GSEA")), function(p){dir.create(p, recursive = T)})


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
##                           4.2.0 Prepare gene sets
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

## Essential DepMap genes implicated in viability 
essential_genes_DEPMAP_2022 <- read_csv("/lustre/scratch123/hgi/projects/otar2065/resources/essential_genes_DEPMAP_2022/common_essentials.csv") %>% 
                                  as.data.frame()
head(essential_genes_DEPMAP_2022)
#               gene
# 1        AAMP (14)
# 2       AARS1 (16)
# 3 AASDHPPT (60496)
# 4       ABCB7 (22)
# 5     ABCE1 (6059)
# 6       ABCF1 (23)

essential_genes_DEPMAP_2022$gene_name <- gsub(" \\(.*", "", essential_genes_DEPMAP_2022$gene)
which(duplicated(essential_genes_DEPMAP_2022$gene_name))
# integer(0)
dim(essential_genes_DEPMAP_2022)
# [1] 1247    1


## Load mouse genes implicated in macrophage survival*
# *from Covarrubias et al. (2020), https://doi.org/10.1016/j.celrep.2020.108541

#  - Based on MAGeCK:
essential_genes_Covarrubias_2020_MAGeCK <- read_csv("/lustre/scratch123/hgi/projects/otar2065/resources/essential_genes_macrophages_Covarrubias_2020/mageck.csv") %>% 
                                                as.data.frame()
head(essential_genes_Covarrubias_2020_MAGeCK)
#             id  p-value      fdr
# 1  mmu-mir-715 1.78e-07 0.000225
# 2    NM_026638 1.78e-07 0.000225
# 3    NM_024212 1.78e-07 0.000225
# 4 NM_001289522 1.78e-07 0.000225
# 5    NM_052835 1.78e-07 0.000225
# 6    NM_016980 1.78e-07 0.000225

which(duplicated(essential_genes_Covarrubias_2020_MAGeCK$id))
# integer(0)
dim(essential_genes_Covarrubias_2020_MAGeCK)
# [1] 417   3

## All significant genes? Yes
which(essential_genes_Covarrubias_2020_MAGeCK$fdr>=0.05)
# integer(0)

## Add human orthologs
mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "mmusculus_gene_ensembl")
## Discard mmu-mir ids
essential_genes_Covarrubias_2020_MAGeCK <- essential_genes_Covarrubias_2020_MAGeCK[grep("mmu", essential_genes_Covarrubias_2020_MAGeCK$id, invert = T), ]

## Get gene name first
mouse_gene_ids <- getBM(values = essential_genes_Covarrubias_2020_MAGeCK$id,
                         mart = mart,
                         attributes = c("refseq_mrna",
                                        "ensembl_gene_id",
                                        "external_gene_name"),
                         filters = "refseq_mrna")
which(duplicated(mouse_gene_ids$refseq_mrna))
# integer(0)
which(duplicated(mouse_gene_ids$external_gene_name))
# integer(0)
mouse_gene_ids$id <- mouse_gene_ids$refseq_mrna
essential_genes_Covarrubias_2020_MAGeCK <- inner_join(essential_genes_Covarrubias_2020_MAGeCK, mouse_gene_ids, by = "id")
which(duplicated(essential_genes_Covarrubias_2020_MAGeCK$external_gene_name))
# integer(0)
dim(essential_genes_Covarrubias_2020_MAGeCK)
# [1] 409   6

## Get human ortholog
mouse_to_human_ids <- getBM(values = essential_genes_Covarrubias_2020_MAGeCK$external_gene_name,
                        mart = mart,
                        attributes = c("external_gene_name",
                                       "hsapiens_homolog_ensembl_gene",
                                       "hsapiens_homolog_associated_gene_name"),
                        filters = "external_gene_name")
## Discard mouse genes without human ortholog
mouse_to_human_ids <- subset(mouse_to_human_ids, hsapiens_homolog_associated_gene_name != "")
dim(mouse_to_human_ids)
# [1] 393   3

## Unique genes in human
which(duplicated(mouse_to_human_ids$hsapiens_homolog_associated_gene_name))
# integer(0)

## 4 mouse genes have >1 human ortholog
which(duplicated(mouse_to_human_ids$external_gene_name))
# [1]  53 254 291 344

## Keep all human genes per mouse gene
essential_genes_Covarrubias_2020_MAGeCK <- inner_join(essential_genes_Covarrubias_2020_MAGeCK, 
                                                      mouse_to_human_ids, by = "external_gene_name", multiple = "all")
dim(essential_genes_Covarrubias_2020_MAGeCK)
# [1] 393   8


#  - Based on Mann–Whitney:
essential_genes_Covarrubias_2020_Mann_Whitney <- read_csv("/lustre/scratch123/hgi/projects/otar2065/resources/essential_genes_macrophages_Covarrubias_2020/mann_whitney.csv") %>% 
                                                    as.data.frame()
head(essential_genes_Covarrubias_2020_Mann_Whitney)
#     gene MW phenotype  p-value
# 1 Gemin5    -6.619930 3.59e-11
# 2   Men1    -6.024644 1.69e-09
# 3 Fermt3    -5.768306 8.01e-09
# 4  Rplp0    -5.737613 9.60e-09
# 5   Ctc1    -5.726412 1.03e-08
# 6  Bub1b    -5.635505 1.75e-08

d <- which(duplicated(essential_genes_Covarrubias_2020_Mann_Whitney$gene))
# 8527 12938
essential_genes_Covarrubias_2020_Mann_Whitney <- essential_genes_Covarrubias_2020_Mann_Whitney[-d, ]
dim(essential_genes_Covarrubias_2020_Mann_Whitney)
# [1] 21377     3

mouse_to_human_ids <- getBM(values = essential_genes_Covarrubias_2020_Mann_Whitney$gene,
                         mart = mart,
                         attributes = c("external_gene_name",
                                      "hsapiens_homolog_ensembl_gene",
                                      "hsapiens_homolog_associated_gene_name"),
                         filters = "external_gene_name")

## Discard mouse genes without human gene
mouse_to_human_ids <- subset(mouse_to_human_ids, hsapiens_homolog_associated_gene_name != '')
dim(mouse_to_human_ids)
# [1] 16999     3

## Mouse genes with >1 human ortholog
length(which(duplicated(mouse_to_human_ids$external_gene_name)))
# [1] 1026

## Human genes with >1 mouse ortholog
length(which(duplicated(mouse_to_human_ids$hsapiens_homolog_associated_gene_name)))
# [1] 1184

mouse_to_human_ids$gene <- mouse_to_human_ids$external_gene_name
essential_genes_Covarrubias_2020_Mann_Whitney <- inner_join(essential_genes_Covarrubias_2020_Mann_Whitney, mouse_to_human_ids, 
                                                                   by = "gene", multiple = "all")
dim(essential_genes_Covarrubias_2020_Mann_Whitney)
# [1] 16999    6

length(unique(essential_genes_Covarrubias_2020_Mann_Whitney$hsapiens_homolog_associated_gene_name))
# [1] 15815

## Retain most signif p for repeated human genes (i.e. with >1 mouse gene)
repeated_human_genes <- unique(essential_genes_Covarrubias_2020_Mann_Whitney$hsapiens_homolog_associated_gene_name[which(duplicated(essential_genes_Covarrubias_2020_Mann_Whitney$hsapiens_homolog_associated_gene_name))])
for(g in repeated_human_genes){
  
  rows <- rownames(subset(essential_genes_Covarrubias_2020_Mann_Whitney, hsapiens_homolog_associated_gene_name == g))
  s <- essential_genes_Covarrubias_2020_Mann_Whitney[rows, ]
  no_most_signif <- rows[rows != rownames(s[order(s$`p-value`, decreasing = F), ][1,])]
  
  ## Keep most signif and remove remaining repeats
  essential_genes_Covarrubias_2020_Mann_Whitney <- essential_genes_Covarrubias_2020_Mann_Whitney[-which(rownames(essential_genes_Covarrubias_2020_Mann_Whitney) %in% no_most_signif), ]
}

dim(essential_genes_Covarrubias_2020_Mann_Whitney)
# [1] 15815     6

## Signif genes
essential_genes_Covarrubias_2020_Mann_Whitney_signif <- subset(essential_genes_Covarrubias_2020_Mann_Whitney, `p-value` < 0.05)
dim(essential_genes_Covarrubias_2020_Mann_Whitney_signif)
# [1] 3175    6






## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
##                   4.2.1 Overrepresentation analysis (ORA)
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

ORA <- function(universe, signif_group, gene_sets){
  
  fisher_ps <- vector()
  fisher_OR <- vector()
  intersect_size <-  vector()
  
  for(i in 1:length(gene_sets)){
    
    gene_set <- gene_sets[[i]]
    gene_set_name <- names(gene_sets)[i]
    
    ## Subset to genes in set that are in universe
    gene_set <- gene_set[gene_set %in% universe]
    
    ## Define gene groups:
    non_signif_group <- universe[!universe %in% signif_group]
    gene_set_signif_group <- signif_group[signif_group %in% gene_set]
    gene_set_non_signif_group <- non_signif_group[non_signif_group %in% gene_set]
    
    no_gene_set_signif_group <- signif_group[!signif_group %in% gene_set]
    no_gene_set_non_signif_group <- non_signif_group[!non_signif_group %in% gene_set]
    
    ## Confirm length of universe, signif group, non-signif group, and gene set
    stopifnot(length(signif_group) + length(non_signif_group) == length(universe))
    stopifnot(length(gene_set_signif_group) + length(no_gene_set_signif_group) == length(signif_group))
    stopifnot(length(gene_set_non_signif_group) + length(no_gene_set_non_signif_group) == length(non_signif_group))
    stopifnot(length(gene_set_signif_group) + length(gene_set_non_signif_group) == length(gene_set))
   
    m <- matrix(data = c(length(gene_set_signif_group), 
                         length(no_gene_set_signif_group), 
                         length(gene_set_non_signif_group), 
                         length(no_gene_set_non_signif_group)), ncol = 2)
    
    f <- fisher.test(m, alternative = "greater")
    
    fisher_ps[gene_set_name] <- f$p.value
    fisher_OR[gene_set_name] <- f$estimate
    intersect_size[gene_set_name] <- paste0(length(gene_set_signif_group), "/",  length(gene_set))
    
  }
  
  return(list("p" = fisher_ps, "OR" = fisher_OR, "intersect" = intersect_size))
  
}

  

#######################  Enrichment among SKAT-O genes  ########################

## Gene sets to test:

# - Viability genes:
DEPMAP_essential_genes <- essential_genes_DEPMAP_2022$gene_name 

# - Macrophage survival genes:
Covarrubias_MAGeCK_macro_survival <- essential_genes_Covarrubias_2020_MAGeCK$hsapiens_homolog_associated_gene_name
Covarrubias_Mann_Whitney_macro_survival <- essential_genes_Covarrubias_2020_Mann_Whitney_signif$hsapiens_homolog_associated_gene_name
Covarrubias_macro_survival <- union(Covarrubias_MAGeCK_macro_survival, Covarrubias_Mann_Whitney_macro_survival)

gene_sets <- list("DEPMAP_essential_genes" = DEPMAP_essential_genes, 
                  "Covarrubias_MAGeCK_macro_survival" = Covarrubias_MAGeCK_macro_survival, 
                  "Covarrubias_Mann_Whitney_macro_survival" = Covarrubias_Mann_Whitney_macro_survival, 
                  "Covarrubias_macro_survival" = Covarrubias_macro_survival)
save(gene_sets, file = paste0(outdir, "/Essential_and_macrophageSurvival_gene_sets.Rdata"))


# ====================== SKAT-O for deleterious variants =======================
## SKAT-O results
deleterious_SKATO_results <- read_csv(paste0("/lustre/scratch123/hgi/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/data/results/2.1.check_SKAT_WES/",
                                             "deleterious_SKATO_scaled_centered_prop_pvals.csv"))
## Subset to genes tested
deleterious_SKATO_results <- subset(deleterious_SKATO_results, ! is.na(resampling_pval))
dim(deleterious_SKATO_results)
# [1] 11410     6
## Per comparison
deleterious_SKATO_results_premac_vs_iPSC <- subset(deleterious_SKATO_results, comparison == "line_prop_changes_premac_iPSC")
deleterious_SKATO_results_old_vs_young_premac <- subset(deleterious_SKATO_results, comparison == "line_prop_changes_old_vs_young_premac")
deleterious_SKATO_results_microglia_vs_premac <- subset(deleterious_SKATO_results, comparison == "line_prop_changes_microglia_premac")

## Signif genes 
deleterious_SKATO_results_signif <- subset(deleterious_SKATO_results, resampling_pval < 0.05)
dim(deleterious_SKATO_results_signif)
# [1] 683   6
## Per comparison
deleterious_SKATO_results_premac_vs_iPSC_signif <- subset(deleterious_SKATO_results_premac_vs_iPSC, resampling_pval < 0.05)
deleterious_SKATO_results_old_vs_young_premac_signif <- subset(deleterious_SKATO_results_old_vs_young_premac, resampling_pval < 0.05)
deleterious_SKATO_results_microglia_vs_premac_signif <- subset(deleterious_SKATO_results_microglia_vs_premac, resampling_pval < 0.05)


# -- -- -- -- -- SKAT-O genes associated to any prolif. comparison -- -- -- -- -

## Universe = genes tested in SKAT-O for at least one comparison
universe <- unique(deleterious_SKATO_results$gene_name)

# Signif genes = genes associated with prolif in at least one comparison
signif_group <- unique(deleterious_SKATO_results_signif$gene_name)
universe_name <- paste0("SKAT-O (Del) any comparison (n = ", length(signif_group), " sig)")

ORA_results <- list()
ps <- vector()
ORs <- vector()
intersections <- vector()
ORA_results[[universe_name]] <- ORA(universe, signif_group, gene_sets)
ps <- as.data.frame(rbind(ps, ORA_results[[universe_name]]$p))
ORs <- as.data.frame(rbind(ORs, ORA_results[[universe_name]]$OR))
intersections <- as.data.frame(rbind(intersections, ORA_results[[universe_name]]$intersect))
rownames(ps) <- universe_name
rownames(ORs) <- universe_name
rownames(intersections) <- universe_name

# -- -- -- -- SKAT-O genes associated to prolif. from iPSC to premac -- -- -- -- 

## Universe = genes tested in SKAT-O for iPSC to premac
universe <- unique(deleterious_SKATO_results_premac_vs_iPSC$gene_name)
# Signif genes = genes associated with prolif from iPSC to premac
signif_group <- unique(deleterious_SKATO_results_premac_vs_iPSC_signif$gene_name)
universe_name <- paste0("SKAT-O (Del) Premac vs iPSC (n = ", length(signif_group), " sig)")

ORA_results[[universe_name]] <- ORA(universe, signif_group, gene_sets)
ps <- rbind(ps, ORA_results[[universe_name]]$p)
ORs <- rbind(ORs, ORA_results[[universe_name]]$OR)
intersections <- rbind(intersections, ORA_results[[universe_name]]$intersect)
rownames(ps)[dim(ps)[1]] <- rownames(ORs)[dim(ORs)[1]] <-  rownames(intersections)[dim(intersections)[1]] <-  universe_name

# -- -- -- SKAT-O genes associated to prolif. from young to old premac -- -- -- 

## Universe = genes tested in SKAT-O for young to old premac
universe <- unique(deleterious_SKATO_results_old_vs_young_premac$gene_name)
# Signif genes = genes associated with prolif from young to old premac
signif_group <- unique(deleterious_SKATO_results_old_vs_young_premac_signif$gene_name)
universe_name <- paste0("SKAT-O (Del) Old vs Young Premac (n = ", length(signif_group), " sig)")

ORA_results[[universe_name]] <- ORA(universe, signif_group, gene_sets)
ps <- rbind(ps, ORA_results[[universe_name]]$p)
ORs <- rbind(ORs, ORA_results[[universe_name]]$OR)
intersections <- rbind(intersections, ORA_results[[universe_name]]$intersect)
rownames(ps)[dim(ps)[1]] <- rownames(ORs)[dim(ORs)[1]] <-  rownames(intersections)[dim(intersections)[1]] <-  universe_name

# -- -- -- SKAT-O genes associated to prolif. from premac to microglia -- -- -- 

## Universe = genes tested in SKAT-O for premac to microglia
universe <- unique(deleterious_SKATO_results_microglia_vs_premac$gene_name)
# Signif genes = genes associated with prolif from premac to microglia
signif_group <- unique(deleterious_SKATO_results_microglia_vs_premac_signif$gene_name)
universe_name <- paste0("SKAT-O (Del) Microglia vs Premac (n = ", length(signif_group), " sig)")

ORA_results[[universe_name]] <- ORA(universe, signif_group, gene_sets)
ps <- rbind(ps, ORA_results[[universe_name]]$p)
ORs <- rbind(ORs, ORA_results[[universe_name]]$OR)
intersections <- rbind(intersections, ORA_results[[universe_name]]$intersect)
rownames(ps)[dim(ps)[1]] <- rownames(ORs)[dim(ORs)[1]] <-  rownames(intersections)[dim(intersections)[1]] <-  universe_name


# ==================== SKAT-O for non-deleterious variants =====================
## SKAT-O results
miss_non_deleterious_SKATO_results <- read_csv(paste0("/lustre/scratch123/hgi/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/data/results/2.1.check_SKAT_WES/",
                                                      "missense_non_deleterious_SKATO_scaled_centered_prop_pvals.csv"))
## Subset to genes tested
miss_non_deleterious_SKATO_results <- subset(miss_non_deleterious_SKATO_results, ! is.na(resampling_pval))
dim(miss_non_deleterious_SKATO_results)
# [1] 19756     6

## Per comparison
miss_non_deleterious_SKATO_results_premac_vs_iPSC <- subset(miss_non_deleterious_SKATO_results, comparison == "line_prop_changes_premac_iPSC")
miss_non_deleterious_SKATO_results_old_vs_young_premac <- subset(miss_non_deleterious_SKATO_results, comparison == "line_prop_changes_old_vs_young_premac")
miss_non_deleterious_SKATO_results_microglia_vs_premac <- subset(miss_non_deleterious_SKATO_results, comparison == "line_prop_changes_microglia_premac")

## Signif genes 
miss_non_deleterious_SKATO_results_signif <- subset(miss_non_deleterious_SKATO_results, resampling_pval < 0.05)
dim(miss_non_deleterious_SKATO_results_signif)
# [1] 1071    6

## Per comparison
miss_non_deleterious_SKATO_results_premac_vs_iPSC_signif <- subset(miss_non_deleterious_SKATO_results_premac_vs_iPSC, resampling_pval < 0.05)
miss_non_deleterious_SKATO_results_old_vs_young_premac_signif <- subset(miss_non_deleterious_SKATO_results_old_vs_young_premac, resampling_pval < 0.05)
miss_non_deleterious_SKATO_results_microglia_vs_premac_signif <- subset(miss_non_deleterious_SKATO_results_microglia_vs_premac, resampling_pval < 0.05)


# -- -- -- -- -- SKAT-O genes associated to any prolif. comparison -- -- -- -- -

universe <- unique(miss_non_deleterious_SKATO_results$gene_name)
signif_group <- unique(miss_non_deleterious_SKATO_results_signif$gene_name)
universe_name <- paste0("SKAT-O (non-Del) any comparison (n = ", length(signif_group), " sig)")

ORA_results[[universe_name]] <- ORA(universe, signif_group, gene_sets)
ps <- rbind(ps, ORA_results[[universe_name]]$p)
ORs <- rbind(ORs, ORA_results[[universe_name]]$OR)
intersections <- rbind(intersections, ORA_results[[universe_name]]$intersect)
rownames(ps)[dim(ps)[1]] <- rownames(ORs)[dim(ORs)[1]] <-  rownames(intersections)[dim(intersections)[1]] <-  universe_name

# -- -- -- -- SKAT-O genes associated to prolif. from iPSC to premac -- -- -- -- 

universe <- unique(miss_non_deleterious_SKATO_results_premac_vs_iPSC$gene_name)
signif_group <- unique(miss_non_deleterious_SKATO_results_premac_vs_iPSC_signif$gene_name)
universe_name <- paste0("SKAT-O (non-Del) Premac vs iPSC (n = ", length(signif_group), " sig)")

ORA_results[[universe_name]] <- ORA(universe, signif_group, gene_sets)
ps <- rbind(ps, ORA_results[[universe_name]]$p)
ORs <- rbind(ORs, ORA_results[[universe_name]]$OR)
intersections <- rbind(intersections, ORA_results[[universe_name]]$intersect)
rownames(ps)[dim(ps)[1]] <- rownames(ORs)[dim(ORs)[1]] <-  rownames(intersections)[dim(intersections)[1]] <-  universe_name

# -- -- -- SKAT-O genes associated to prolif. from young to old premac -- -- -- 

universe <- unique(miss_non_deleterious_SKATO_results_old_vs_young_premac$gene_name)
signif_group <- unique(miss_non_deleterious_SKATO_results_old_vs_young_premac_signif$gene_name)
universe_name <- paste0("SKAT-O (non-Del) Old vs Young Premac (n = ", length(signif_group), " sig)")

ORA_results[[universe_name]] <- ORA(universe, signif_group, gene_sets)
ps <- rbind(ps, ORA_results[[universe_name]]$p)
ORs <- rbind(ORs, ORA_results[[universe_name]]$OR)
intersections <- rbind(intersections, ORA_results[[universe_name]]$intersect)
rownames(ps)[dim(ps)[1]] <- rownames(ORs)[dim(ORs)[1]] <-  rownames(intersections)[dim(intersections)[1]] <-  universe_name

# -- -- -- SKAT-O genes associated to prolif. from premac to microglia -- -- -- 

universe <- unique(miss_non_deleterious_SKATO_results_microglia_vs_premac$gene_name)
signif_group <- unique(miss_non_deleterious_SKATO_results_microglia_vs_premac_signif$gene_name)
universe_name <- paste0("SKAT-O (non-Del) Microglia vs Premac (n = ", length(signif_group), " sig)")

ORA_results[[universe_name]] <- ORA(universe, signif_group, gene_sets)
ps <- rbind(ps, ORA_results[[universe_name]]$p)
ORs <- rbind(ORs, ORA_results[[universe_name]]$OR)
intersections <- rbind(intersections, ORA_results[[universe_name]]$intersect)
rownames(ps)[dim(ps)[1]] <- rownames(ORs)[dim(ORs)[1]] <-  rownames(intersections)[dim(intersections)[1]] <-  universe_name


## Heatmap
ps <- t(ps)
ORs <- t(ORs)
intersections <- t(intersections)
ORA_results_SKAT_O <- list("p" = ps, "OR" = ORs, "intersection_size" = intersections)
save(ORA_results_SKAT_O, file = paste0(outdir, "/01_ORA", "/ORA_results_enrichment_in_SKAT_O_genes.Rdata"))

h_row_title <- "Sets of viability and survival genes"
h_row_labels <- c("DEPMAP_essential_genes" = "DepMap essential genes", 
                  "Covarrubias_macro_survival" = "Covarrubias all macrophage survival genes",
                  "Covarrubias_Mann_Whitney_macro_survival" = "Covarrubias Mann Whitney macrophage survival genes",
                  "Covarrubias_MAGeCK_macro_survival" = "Covarrubias MAGeCK macrophage survival genes")

h <- ComplexHeatmap::Heatmap(-log10(as.matrix(ps)), 
                             row_title = h_row_title, 
                             row_title_gp = gpar(fontsize = 10),
                             name = "-log10(pval)", 
                             col = colorRamp2(c(0, 0.5, 1, 1.5, 2), c("white", "lavenderblush2", "lightpink", "salmon", "salmon4")), 
                             border = T, 
                             rect_gp = gpar(col = "black", lwd = 1),
                             row_names_gp = gpar(fontsize = 8.5), 
                             column_names_gp = gpar(fontsize = 6.5), 
                             row_names_max_width = unit(10, "cm"),
                             cluster_columns = F, 
                             cluster_rows = F, 
                             row_order = c("DEPMAP_essential_genes",
                                           "Covarrubias_macro_survival",
                                           "Covarrubias_Mann_Whitney_macro_survival",
                                           "Covarrubias_MAGeCK_macro_survival"),
                             row_labels = h_row_labels[rownames(ps)], 
                             cell_fun = function(j, i, x, y, width, height, fill) {
                               if(log10(as.matrix(ps))[i,j] < log10(0.05)){
                                 grid.text("*", x, y, gp = gpar(fontsize = 12.5, col = "yellow"), just = c(0.5,-1.2))
                                        grid.text(paste0(signif(ORs[i, j], digits = 3), "\n", "(", intersections[i,j], ")"), 
                                                  x, y, gp = gpar(fontsize = 8.5)) 
                                 }
                               else{
                                 grid.text(paste0(signif(ORs[i, j], digits = 3), "\n", "(", intersections[i,j], ")"), x, y, gp = gpar(fontsize = 8.5))
                               }
                             })

pdf(file = paste0(plotdir, "/01_ORA", "/ORA_results_enrichment_in_SKAT_O_genes", ".pdf"), height = 7, width = 10)
h
dev.off()



####################  Enrichment among Burden test genes  ######################

# =================== Burden tests for deleterious variants ====================
## Burden test results (all genes tested)
deleterious_burden_results <- as.data.frame(read_csv("/lustre/scratch123/hgi/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/data/results/2.2.rare_vars_vs_prolif/del_burden_scaled_centered_prop_pvals.csv"))
dim(deleterious_burden_results)
# [1] 487   5

## Per comparison
deleterious_burden_results_premac_vs_iPSC <- subset(deleterious_burden_results, comparison == "line_prop_changes_premac_iPSC")
deleterious_burden_results_old_vs_young_premac <- subset(deleterious_burden_results, comparison == "line_prop_changes_old_vs_young_premac")
deleterious_burden_results_microglia_vs_premac <- subset(deleterious_burden_results, comparison == "line_prop_changes_microglia_premac")

## Signif genes
deleterious_burden_results_signif <- subset(deleterious_burden_results, p_Bonf < 0.05)
dim(deleterious_burden_results_signif)
# [1] 22  5

## Per comparison
deleterious_burden_results_premac_vs_iPSC_signif <- subset(deleterious_burden_results_premac_vs_iPSC, p_Bonf < 0.05)
deleterious_burden_results_old_vs_young_premac_signif <- subset(deleterious_burden_results_old_vs_young_premac, p_Bonf < 0.05)
deleterious_burden_results_microglia_vs_premac_signif <- subset(deleterious_burden_results_microglia_vs_premac, p_Bonf < 0.05)


# -- -- -- -- Burden test genes associated to any prolif. comparison -- -- -- --

## Universe = genes tested in burden test for at least one comparison
universe <- unique(deleterious_burden_results$gene_name)

# Signif genes = genes associated with prolif in at least one comparison
signif_group <- unique(deleterious_burden_results_signif$gene_name)
universe_name <- paste0("Burden test (Del) any comparison (n = ", length(signif_group), " sig)")

ORA_results <- list()
ps <- vector()
ORs <- vector()
intersections <- vector()
ORA_results[[universe_name]] <- ORA(universe, signif_group, gene_sets)
ps <- as.data.frame(rbind(ps, ORA_results[[universe_name]]$p))
ORs <- as.data.frame(rbind(ORs, ORA_results[[universe_name]]$OR))
intersections <- as.data.frame(rbind(intersections, ORA_results[[universe_name]]$intersect))
rownames(ps) <- universe_name
rownames(ORs) <- universe_name
rownames(intersections) <- universe_name

# -- -- -- Burden test genes associated to prolif. from iPSC to premac -- -- -- 

universe <- unique(deleterious_burden_results_premac_vs_iPSC$gene_name)
signif_group <- unique(deleterious_burden_results_premac_vs_iPSC_signif$gene_name)
universe_name <- paste0("Burden test (Del) Premac vs iPSC (n = ", length(signif_group), " sig)")

ORA_results[[universe_name]] <- ORA(universe, signif_group, gene_sets)
ps <- rbind(ps, ORA_results[[universe_name]]$p)
ORs <- rbind(ORs, ORA_results[[universe_name]]$OR)
intersections <- rbind(intersections, ORA_results[[universe_name]]$intersect)
rownames(ps)[dim(ps)[1]] <- rownames(ORs)[dim(ORs)[1]] <-  rownames(intersections)[dim(intersections)[1]] <-  universe_name

# -- -- - Burden test genes associated to prolif. from young to old premac -- -- 

universe <- unique(deleterious_burden_results_old_vs_young_premac$gene_name)
signif_group <- unique(deleterious_burden_results_old_vs_young_premac_signif$gene_name)
universe_name <- paste0("Burden test (Del) Old vs Young Premac (n = ", length(signif_group), " sig)")

ORA_results[[universe_name]] <- ORA(universe, signif_group, gene_sets)
ps <- rbind(ps, ORA_results[[universe_name]]$p)
ORs <- rbind(ORs, ORA_results[[universe_name]]$OR)
intersections <- rbind(intersections, ORA_results[[universe_name]]$intersect)
rownames(ps)[dim(ps)[1]] <- rownames(ORs)[dim(ORs)[1]] <-  rownames(intersections)[dim(intersections)[1]] <-  universe_name

# -- -- - Burden test genes associated to prolif. from premac to microglia -- -- 

universe <- unique(deleterious_burden_results_microglia_vs_premac$gene_name)
signif_group <- unique(deleterious_burden_results_microglia_vs_premac_signif$gene_name)
universe_name <- paste0("Burden test (Del) Premac vs Microglia (n = ", length(signif_group), " sig)")

ORA_results[[universe_name]] <- ORA(universe, signif_group, gene_sets)
ps <- rbind(ps, ORA_results[[universe_name]]$p)
ORs <- rbind(ORs, ORA_results[[universe_name]]$OR)
intersections <- rbind(intersections, ORA_results[[universe_name]]$intersect)
rownames(ps)[dim(ps)[1]] <- rownames(ORs)[dim(ORs)[1]] <-  rownames(intersections)[dim(intersections)[1]] <-  universe_name


# ================= Burden tests for non-deleterious variants ==================

## Burden test results (all genes tested)
non_deleterious_burden_results <- as.data.frame(read_csv("/lustre/scratch123/hgi/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/data/results/2.2.rare_vars_vs_prolif/ptv_burden_scaled_centered_prop_pvals.csv"))
dim(non_deleterious_burden_results)
# [1] 84  5

## Per comparison
non_deleterious_burden_results_premac_vs_iPSC <- subset(non_deleterious_burden_results, comparison == "line_prop_changes_premac_iPSC")
non_deleterious_burden_results_old_vs_young_premac <- subset(non_deleterious_burden_results, comparison == "line_prop_changes_old_vs_young_premac")
non_deleterious_burden_results_microglia_vs_premac <- subset(non_deleterious_burden_results, comparison == "line_prop_changes_microglia_premac")

## Signif genes
non_deleterious_burden_results_signif <- subset(non_deleterious_burden_results, p_Bonf < 0.05)
dim(non_deleterious_burden_results_signif)
# [1] 5  5

## Per comparison
non_deleterious_burden_results_premac_vs_iPSC_signif <- subset(non_deleterious_burden_results_premac_vs_iPSC, p_Bonf < 0.05)
# non_deleterious_burden_results_old_vs_young_premac_signif <- subset(non_deleterious_burden_results_old_vs_young_premac, p_Bonf < 0.05)   zero
non_deleterious_burden_results_microglia_vs_premac_signif <- subset(non_deleterious_burden_results_microglia_vs_premac, p_Bonf < 0.05)


# -- -- -- -- Burden test genes associated to any prolif. comparison -- -- -- --

## Universe 
universe <- unique(non_deleterious_burden_results$gene_name)
# Signif genes 
signif_group <- unique(non_deleterious_burden_results_signif$gene_name)
universe_name <- paste0("Burden test (non-Del) any comparison (n = ", length(signif_group), " sig)")

ORA_results[[universe_name]] <- ORA(universe, signif_group, gene_sets)
ps <- as.data.frame(rbind(ps, ORA_results[[universe_name]]$p))
ORs <- as.data.frame(rbind(ORs, ORA_results[[universe_name]]$OR))
intersections <- as.data.frame(rbind(intersections, ORA_results[[universe_name]]$intersect))
rownames(ps)[dim(ps)[1]] <- rownames(ORs)[dim(ORs)[1]] <-  rownames(intersections)[dim(intersections)[1]] <-  universe_name

# -- -- -- Burden test genes associated to prolif. from iPSC to premac -- -- -- 

universe <- unique(non_deleterious_burden_results_premac_vs_iPSC$gene_name)
signif_group <- unique(non_deleterious_burden_results_premac_vs_iPSC_signif$gene_name)
universe_name <- paste0("Burden test (non-Del) Premac vs iPSC (n = ", length(signif_group), " sig)")

ORA_results[[universe_name]] <- ORA(universe, signif_group, gene_sets)
ps <- rbind(ps, ORA_results[[universe_name]]$p)
ORs <- rbind(ORs, ORA_results[[universe_name]]$OR)
intersections <- rbind(intersections, ORA_results[[universe_name]]$intersect)
rownames(ps)[dim(ps)[1]] <- rownames(ORs)[dim(ORs)[1]] <-  rownames(intersections)[dim(intersections)[1]] <-  universe_name

# -- -- - Burden test genes associated to prolif. from premac to microglia -- -- 

universe <- unique(non_deleterious_burden_results_microglia_vs_premac$gene_name)
signif_group <- unique(non_deleterious_burden_results_microglia_vs_premac_signif$gene_name)
universe_name <- paste0("Burden test (non-Del) Premac vs Microglia (n = ", length(signif_group), " sig)")

ORA_results[[universe_name]] <- ORA(universe, signif_group, gene_sets)
ps <- rbind(ps, ORA_results[[universe_name]]$p)
ORs <- rbind(ORs, ORA_results[[universe_name]]$OR)
intersections <- rbind(intersections, ORA_results[[universe_name]]$intersect)
rownames(ps)[dim(ps)[1]] <- rownames(ORs)[dim(ORs)[1]] <-  rownames(intersections)[dim(intersections)[1]] <-  universe_name


# =================== Burden tests for synonymous variants =====================

## Burden test results (all genes tested)
syn_burden_results <- as.data.frame(read_csv("/lustre/scratch123/hgi/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/data/results/2.2.rare_vars_vs_prolif/syn_burden_scaled_centered_prop_pvals.csv"))
dim(syn_burden_results)
# [1] 365  5

## Per comparison
syn_burden_results_premac_vs_iPSC <- subset(syn_burden_results, comparison == "line_prop_changes_premac_iPSC")
syn_burden_results_old_vs_young_premac <- subset(syn_burden_results, comparison == "line_prop_changes_old_vs_young_premac")
syn_burden_results_microglia_vs_premac <- subset(syn_burden_results, comparison == "line_prop_changes_microglia_premac")

## Signif genes
syn_burden_results_signif <- subset(syn_burden_results, p_Bonf < 0.05)
dim(syn_burden_results_signif)
# [1] 2  5

## Per comparison
# syn_burden_results_premac_vs_iPSC_signif <- subset(syn_burden_results_premac_vs_iPSC, p_Bonf < 0.05)    zero
# syn_burden_results_old_vs_young_premac_signif <- subset(syn_burden_results_old_vs_young_premac, p_Bonf < 0.05)   zero
syn_burden_results_microglia_vs_premac_signif <- subset(syn_burden_results_microglia_vs_premac, p_Bonf < 0.05)


# -- -- -- -- Burden test genes associated to any prolif. comparison -- -- -- --

universe <- unique(syn_burden_results$gene_name)
signif_group <- unique(syn_burden_results_signif$gene_name)
universe_name <- paste0("Burden test (Syn) any comparison (n = ", length(signif_group), " sig)")

ORA_results[[universe_name]] <- ORA(universe, signif_group, gene_sets)
ps <- rbind(ps, ORA_results[[universe_name]]$p)
ORs <- rbind(ORs, ORA_results[[universe_name]]$OR)
intersections <- rbind(intersections, ORA_results[[universe_name]]$intersect)
rownames(ps)[dim(ps)[1]] <- rownames(ORs)[dim(ORs)[1]] <-  rownames(intersections)[dim(intersections)[1]] <-  universe_name

# -- -- - Burden test genes associated to prolif. from premac to microglia -- -- 

universe <- unique(syn_burden_results_microglia_vs_premac$gene_name)
signif_group <- unique(syn_burden_results_microglia_vs_premac_signif$gene_name)
universe_name <- paste0("Burden test (Syn) Premac vs Microglia (n = ", length(signif_group), " sig)")

ORA_results[[universe_name]] <- ORA(universe, signif_group, gene_sets)
ps <- rbind(ps, ORA_results[[universe_name]]$p)
ORs <- rbind(ORs, ORA_results[[universe_name]]$OR)
intersections <- rbind(intersections, ORA_results[[universe_name]]$intersect)
rownames(ps)[dim(ps)[1]] <- rownames(ORs)[dim(ORs)[1]] <-  rownames(intersections)[dim(intersections)[1]] <-  universe_name


## Heatmap
ps <- t(ps)
ORs <- t(ORs)
intersections <- t(intersections)
ORA_results_Burden_tests <- list("p" = ps, "OR" = ORs, "intersection_size" = intersections)
save(ORA_results_Burden_tests, file = paste0(outdir, "/01_ORA", "/ORA_results_enrichment_in_Burden_test_genes.Rdata"))

h_row_title <- "Sets of viability and survival genes"
h_row_labels <- c("DEPMAP_essential_genes" = "DepMap essential genes", 
                  "Covarrubias_macro_survival" = "Covarrubias all macrophage survival genes",
                  "Covarrubias_Mann_Whitney_macro_survival" = "Covarrubias Mann Whitney macrophage survival genes",
                  "Covarrubias_MAGeCK_macro_survival" = "Covarrubias MAGeCK macrophage survival genes")

h <- ComplexHeatmap::Heatmap(-log10(as.matrix(ps)), 
                             row_title = h_row_title, 
                             row_title_gp = gpar(fontsize = 10),
                             name = "-log10(pval)", 
                             col = colorRamp2(c(0, 0.5, 1, 1.5, 2), c("white", "lavenderblush2", "lightpink", "salmon", "salmon4")), 
                             border = T, 
                             rect_gp = gpar(col = "black", lwd = 1),
                             row_names_gp = gpar(fontsize = 8.5), 
                             column_names_gp = gpar(fontsize = 6.5), 
                             row_names_max_width = unit(10, "cm"),
                             cluster_columns = F, 
                             cluster_rows = F, 
                             row_order = c("DEPMAP_essential_genes",
                                           "Covarrubias_macro_survival",
                                           "Covarrubias_Mann_Whitney_macro_survival",
                                           "Covarrubias_MAGeCK_macro_survival"),
                             row_labels = h_row_labels[rownames(ps)], 
                             cell_fun = function(j, i, x, y, width, height, fill) {
                               if(log10(as.matrix(ps))[i,j] < log10(0.05)){
                                 grid.text("*", x, y, gp = gpar(fontsize = 12.5, col = "yellow"), just = c(0.5,-1.2))
                                 grid.text(paste0(signif(ORs[i, j], digits = 3), "\n", "(", intersections[i,j], ")"), 
                                           x, y, gp = gpar(fontsize = 8.5)) 
                               }
                               else{
                                 grid.text(paste0(signif(ORs[i, j], digits = 3), "\n", "(", intersections[i,j], ")"), x, y, gp = gpar(fontsize = 8.5))
                               }
                             })

pdf(file = paste0(plotdir, "/01_ORA", "/ORA_results_enrichment_in_Burden_test_genes", ".pdf"), height = 7, width = 10)
h
dev.off()



###############  Enrichment among SKAT-O signif negative genes  ################

## Add gene sets of Mann-Whitney genes divided by sign
gene_sets[["Covarrubias_Mann_Whitney_macro_survival_pos"]] <- subset(essential_genes_Covarrubias_2020_Mann_Whitney_signif, `MW phenotype` > 0)$hsapiens_homolog_associated_gene_name
gene_sets[["Covarrubias_Mann_Whitney_macro_survival_neg"]] <- subset(essential_genes_Covarrubias_2020_Mann_Whitney_signif, `MW phenotype` < 0)$hsapiens_homolog_associated_gene_name

# ====================== SKAT-O for deleterious variants =======================

# -- -- -- -- -- SKAT-O genes decreasing any prolif. comparison -- -- -- -- -

## Universe = genes tested in SKAT-O for at least one comparison
universe <- unique(deleterious_SKATO_results$gene_name)

## All genes tested for burden were SKAT signif
which(!deleterious_burden_results$gene_name %in% deleterious_SKATO_results_signif$gene_name)
# integer(0)

# Signif genes = genes negatively associated with prolif in at least one comparison
neg_signif_group <- unique(subset(deleterious_burden_results, coef<0)$gene_name)
universe_name <- paste0("SKAT-O (Del) any comparison (n = ", length(neg_signif_group), " sig negative)")

ORA_results <- list()
ps <- vector()
ORs <- vector()
intersections <- vector()
ORA_results[[universe_name]] <- ORA(universe, neg_signif_group, gene_sets)
ps <- as.data.frame(rbind(ps, ORA_results[[universe_name]]$p))
ORs <- as.data.frame(rbind(ORs, ORA_results[[universe_name]]$OR))
intersections <- as.data.frame(rbind(intersections, ORA_results[[universe_name]]$intersect))
rownames(ps) <- rownames(ORs)  <- rownames(intersections)  <- universe_name

# -- -- -- -- - SKAT-O genes decreasing prolif. from iPSC to premac -- -- -- -- 

universe <- unique(deleterious_SKATO_results_premac_vs_iPSC$gene_name)
# Signif genes = genes negatively associated with prolif from iPSC to premac
neg_signif_group <-  unique(subset(deleterious_burden_results, comparison == "line_prop_changes_premac_iPSC" & coef<0)$gene_name)
which(!neg_signif_group %in% deleterious_SKATO_results_premac_vs_iPSC_signif$gene_name)
# integer(0)
universe_name <- paste0("SKAT-O (Del) Premac vs iPSC (n = ", length(neg_signif_group), " sig negative)")

ORA_results[[universe_name]] <- ORA(universe, neg_signif_group, gene_sets)
ps <- rbind(ps, ORA_results[[universe_name]]$p)
ORs <- rbind(ORs, ORA_results[[universe_name]]$OR)
intersections <- rbind(intersections, ORA_results[[universe_name]]$intersect)
rownames(ps)[dim(ps)[1]] <- rownames(ORs)[dim(ORs)[1]] <-  rownames(intersections)[dim(intersections)[1]] <-  universe_name

# -- -- -- SKAT-O genes decreasing prolif. from young to old premac -- -- -- --

universe <- unique(deleterious_SKATO_results_old_vs_young_premac$gene_name)
# Signif genes = genes negatively associated with prolif from young to old premac
neg_signif_group <-  unique(subset(deleterious_burden_results, comparison == "line_prop_changes_old_vs_young_premac" & coef<0)$gene_name)
which(!neg_signif_group %in% deleterious_SKATO_results_old_vs_young_premac_signif$gene_name)
# integer(0)
universe_name <- paste0("SKAT-O (Del) Old vs Young Premac (n = ", length(neg_signif_group), " sig negative)")

ORA_results[[universe_name]] <- ORA(universe, neg_signif_group, gene_sets)
ps <- rbind(ps, ORA_results[[universe_name]]$p)
ORs <- rbind(ORs, ORA_results[[universe_name]]$OR)
intersections <- rbind(intersections, ORA_results[[universe_name]]$intersect)
rownames(ps)[dim(ps)[1]] <- rownames(ORs)[dim(ORs)[1]] <-  rownames(intersections)[dim(intersections)[1]] <-  universe_name

# -- -- -- SKAT-O genes decreasing prolif. from premac to microglia -- -- -- --

universe <- unique(deleterious_SKATO_results_microglia_vs_premac$gene_name)
# Signif genes = genes negatively associated with prolif from premac to microglia
neg_signif_group <-  unique(subset(deleterious_burden_results, comparison == "line_prop_changes_microglia_premac" & coef<0)$gene_name)
which(!neg_signif_group %in% deleterious_SKATO_results_microglia_vs_premac_signif$gene_name)
# integer(0)
universe_name <- paste0("SKAT-O (Del) Microglia vs Premac (n = ", length(neg_signif_group), " sig negative)")

ORA_results[[universe_name]] <- ORA(universe, neg_signif_group, gene_sets)
ps <- rbind(ps, ORA_results[[universe_name]]$p)
ORs <- rbind(ORs, ORA_results[[universe_name]]$OR)
intersections <- rbind(intersections, ORA_results[[universe_name]]$intersect)
rownames(ps)[dim(ps)[1]] <- rownames(ORs)[dim(ORs)[1]] <-  rownames(intersections)[dim(intersections)[1]] <-  universe_name


## Heatmap
ps <- t(ps)
ORs <- t(ORs)
intersections <- t(intersections)
ORA_results_SKAT_O_neg <- list("p" = ps, "OR" = ORs, "intersection_size" = intersections)
save(ORA_results_SKAT_O_neg, file = paste0(outdir, "/01_ORA", "/ORA_results_enrichment_in_SKAT_O_negative_genes.Rdata"))

h_row_title <- "Sets of viability and survival genes"
h_row_labels <- c("DEPMAP_essential_genes" = "DepMap essential genes", 
                  "Covarrubias_macro_survival" = "Covarrubias all macrophage survival genes",
                  "Covarrubias_MAGeCK_macro_survival" = "Covarrubias MAGeCK macrophage survival genes",
                  "Covarrubias_Mann_Whitney_macro_survival" = "Covarrubias Mann Whitney macrophage survival genes",
                  "Covarrubias_Mann_Whitney_macro_survival_pos" = "Covarrubias Mann Whitney macrophage survival positive genes",
                  "Covarrubias_Mann_Whitney_macro_survival_neg" = "Covarrubias Mann Whitney macrophage survival negative genes")

h <- ComplexHeatmap::Heatmap(-log10(as.matrix(ps)), 
                             row_title = h_row_title, 
                             row_title_gp = gpar(fontsize = 10),
                             name = "-log10(pval)", 
                             col = colorRamp2(c(0, 0.5, 1, 1.5, 2), c("white", "lavenderblush2", "lightpink", "salmon", "salmon4")), 
                             border = T, 
                             rect_gp = gpar(col = "black", lwd = 1),
                             row_names_gp = gpar(fontsize = 8.5), 
                             column_names_gp = gpar(fontsize = 6), 
                             row_names_max_width = unit(10, "cm"),
                             cluster_columns = F, 
                             cluster_rows = F, 
                             row_order = c("DEPMAP_essential_genes",
                                           "Covarrubias_macro_survival",
                                           "Covarrubias_MAGeCK_macro_survival",
                                           "Covarrubias_Mann_Whitney_macro_survival",
                                           "Covarrubias_Mann_Whitney_macro_survival_pos",
                                           "Covarrubias_Mann_Whitney_macro_survival_neg"),
                             row_labels = h_row_labels[rownames(ps)], 
                             cell_fun = function(j, i, x, y, width, height, fill) {
                               if(log10(as.matrix(ps))[i,j] < log10(0.05)){
                                 grid.text("*", x, y, gp = gpar(fontsize = 12.5, col = "yellow"), just = c(0.5,-1.2))
                                 grid.text(paste0(signif(ORs[i, j], digits = 3), "\n", "(", intersections[i,j], ")"), 
                                           x, y, gp = gpar(fontsize = 8.5)) 
                               }
                               else{
                                 grid.text(paste0(signif(ORs[i, j], digits = 3), "\n", "(", intersections[i,j], ")"), x, y, gp = gpar(fontsize = 8.5))
                               }
                             })

pdf(file = paste0(plotdir, "/01_ORA", "/ORA_results_enrichment_in_SKAT_O_negative_genes", ".pdf"), height = 7, width = 10)
h
dev.off()



###############  Enrichment among SKAT-O signif positive genes  ################

# ====================== SKAT-O for deleterious variants =======================

# -- -- -- -- -- SKAT-O genes increasing any prolif. comparison -- -- -- -- -

## Universe = genes tested in SKAT-O for at least one comparison
universe <- unique(deleterious_SKATO_results$gene_name)

# Signif genes = genes positively associated with prolif in at least one comparison
pos_signif_group <- unique(subset(deleterious_burden_results, coef>0)$gene_name)
universe_name <- paste0("SKAT-O (Del) any comparison (n = ", length(pos_signif_group), " sig positive)")

ORA_results <- list()
ps <- vector()
ORs <- vector()
intersections <- vector()
ORA_results[[universe_name]] <- ORA(universe, pos_signif_group, gene_sets)
ps <- as.data.frame(rbind(ps, ORA_results[[universe_name]]$p))
ORs <- as.data.frame(rbind(ORs, ORA_results[[universe_name]]$OR))
intersections <- as.data.frame(rbind(intersections, ORA_results[[universe_name]]$intersect))
rownames(ps) <- rownames(ORs)  <- rownames(intersections)  <- universe_name

# -- -- -- -- - SKAT-O genes increasing prolif. from iPSC to premac -- -- -- -- 

universe <- unique(deleterious_SKATO_results_premac_vs_iPSC$gene_name)
# Signif genes = genes positively associated with prolif from iPSC to premac
pos_signif_group <-  unique(subset(deleterious_burden_results, comparison == "line_prop_changes_premac_iPSC" & coef>0)$gene_name)
which(!pos_signif_group %in% deleterious_SKATO_results_premac_vs_iPSC_signif$gene_name)
# integer(0)
universe_name <- paste0("SKAT-O (Del) Premac vs iPSC (n = ", length(pos_signif_group), " sig positive)")

ORA_results[[universe_name]] <- ORA(universe, pos_signif_group, gene_sets)
ps <- rbind(ps, ORA_results[[universe_name]]$p)
ORs <- rbind(ORs, ORA_results[[universe_name]]$OR)
intersections <- rbind(intersections, ORA_results[[universe_name]]$intersect)
rownames(ps)[dim(ps)[1]] <- rownames(ORs)[dim(ORs)[1]] <-  rownames(intersections)[dim(intersections)[1]] <-  universe_name

# -- -- -- SKAT-O genes increasing prolif. from young to old premac -- -- -- --

universe <- unique(deleterious_SKATO_results_old_vs_young_premac$gene_name)
# Signif genes = genes positively associated with prolif from young to old premac
pos_signif_group <-  unique(subset(deleterious_burden_results, comparison == "line_prop_changes_old_vs_young_premac" & coef>0)$gene_name)
which(!pos_signif_group %in% deleterious_SKATO_results_old_vs_young_premac_signif$gene_name)
# integer(0)
universe_name <- paste0("SKAT-O (Del) Old vs Young Premac (n = ", length(pos_signif_group), " sig positive)")

ORA_results[[universe_name]] <- ORA(universe, pos_signif_group, gene_sets)
ps <- rbind(ps, ORA_results[[universe_name]]$p)
ORs <- rbind(ORs, ORA_results[[universe_name]]$OR)
intersections <- rbind(intersections, ORA_results[[universe_name]]$intersect)
rownames(ps)[dim(ps)[1]] <- rownames(ORs)[dim(ORs)[1]] <-  rownames(intersections)[dim(intersections)[1]] <-  universe_name

# -- -- -- SKAT-O genes increasing prolif. from premac to microglia -- -- -- --

universe <- unique(deleterious_SKATO_results_microglia_vs_premac$gene_name)
# Signif genes = genes positively associated with prolif from premac to microglia
pos_signif_group <-  unique(subset(deleterious_burden_results, comparison == "line_prop_changes_microglia_premac" & coef>0)$gene_name)
which(!pos_signif_group %in% deleterious_SKATO_results_microglia_vs_premac_signif$gene_name)
# integer(0)
universe_name <- paste0("SKAT-O (Del) Microglia vs Premac (n = ", length(pos_signif_group), " sig positive)")

ORA_results[[universe_name]] <- ORA(universe, pos_signif_group, gene_sets)
ps <- rbind(ps, ORA_results[[universe_name]]$p)
ORs <- rbind(ORs, ORA_results[[universe_name]]$OR)
intersections <- rbind(intersections, ORA_results[[universe_name]]$intersect)
rownames(ps)[dim(ps)[1]] <- rownames(ORs)[dim(ORs)[1]] <-  rownames(intersections)[dim(intersections)[1]] <-  universe_name


## Heatmap
ps <- t(ps)
ORs <- t(ORs)
intersections <- t(intersections)
ORA_results_SKAT_O_pos <- list("p" = ps, "OR" = ORs, "intersection_size" = intersections)
save(ORA_results_SKAT_O_pos, file = paste0(outdir, "/01_ORA", "/ORA_results_enrichment_in_SKAT_O_positive_genes.Rdata"))

h_row_title <- "Sets of viability and survival genes"
h_row_labels <- c("DEPMAP_essential_genes" = "DepMap essential genes", 
                  "Covarrubias_macro_survival" = "Covarrubias all macrophage survival genes",
                  "Covarrubias_MAGeCK_macro_survival" = "Covarrubias MAGeCK macrophage survival genes",
                  "Covarrubias_Mann_Whitney_macro_survival" = "Covarrubias Mann Whitney macrophage survival genes",
                  "Covarrubias_Mann_Whitney_macro_survival_pos" = "Covarrubias Mann Whitney macrophage survival positive genes",
                  "Covarrubias_Mann_Whitney_macro_survival_neg" = "Covarrubias Mann Whitney macrophage survival negative genes")

h <- ComplexHeatmap::Heatmap(-log10(as.matrix(ps)), 
                             row_title = h_row_title, 
                             row_title_gp = gpar(fontsize = 10),
                             name = "-log10(pval)", 
                             col = colorRamp2(c(0, 0.5, 1, 1.5, 2), c("white", "lavenderblush2", "lightpink", "salmon", "salmon4")), 
                             border = T, 
                             rect_gp = gpar(col = "black", lwd = 1),
                             row_names_gp = gpar(fontsize = 8.5), 
                             column_names_gp = gpar(fontsize = 6), 
                             row_names_max_width = unit(10, "cm"),
                             cluster_columns = F, 
                             cluster_rows = F, 
                             row_order = c("DEPMAP_essential_genes",
                                           "Covarrubias_macro_survival",
                                           "Covarrubias_MAGeCK_macro_survival",
                                           "Covarrubias_Mann_Whitney_macro_survival",
                                           "Covarrubias_Mann_Whitney_macro_survival_pos",
                                           "Covarrubias_Mann_Whitney_macro_survival_neg"),
                             row_labels = h_row_labels[rownames(ps)], 
                             cell_fun = function(j, i, x, y, width, height, fill) {
                               if(log10(as.matrix(ps))[i,j] < log10(0.05)){
                                 grid.text("*", x, y, gp = gpar(fontsize = 12.5, col = "yellow"), just = c(0.5,-1.2))
                                 grid.text(paste0(signif(ORs[i, j], digits = 3), "\n", "(", intersections[i,j], ")"), 
                                           x, y, gp = gpar(fontsize = 8.5)) 
                               }
                               else{
                                 grid.text(paste0(signif(ORs[i, j], digits = 3), "\n", "(", intersections[i,j], ")"), x, y, gp = gpar(fontsize = 8.5))
                               }
                             })

pdf(file = paste0(plotdir, "/01_ORA", "/ORA_results_enrichment_in_SKAT_O_positive_genes", ".pdf"), height = 7, width = 10)
h
dev.off()



####################  Enrichment among Mann–Whitney genes  #####################

universe <- unique(essential_genes_Covarrubias_2020_Mann_Whitney$hsapiens_homolog_associated_gene_name)
signif_group <- unique(essential_genes_Covarrubias_2020_Mann_Whitney_signif$hsapiens_homolog_associated_gene_name)
universe_name <- paste0("Mann-Whitney macrophage survival genes (n = ", length(signif_group), " sig)")

## Gene sets = genes associated with prolif based on del variants 
gene_sets <- list("SKAT-O (Del) genes for any comparison" = unique(deleterious_SKATO_results_signif$gene_name),
                  "SKAT-O (Del) genes for Premac vs iPSC" = unique(deleterious_SKATO_results_premac_vs_iPSC_signif$gene_name),
                  "SKAT-O (Del) genes for Old vs Young Premac" = unique(deleterious_SKATO_results_old_vs_young_premac_signif$gene_name),
                  "SKAT-O (Del) genes for Microglia vs Premac" = unique(deleterious_SKATO_results_microglia_vs_premac_signif$gene_name), 
                  "Burden test (Del) genes for any comparison" = unique(deleterious_burden_results_signif$gene_name),
                  "Burden test (Del) genes for Premac vs iPSC" = unique(deleterious_burden_results_premac_vs_iPSC_signif$gene_name),
                  "Burden test (Del) genes for Old vs Young Premac" = unique(deleterious_burden_results_old_vs_young_premac_signif$gene_name),
                  "Burden test (Del) genes for Microglia vs Premac" = unique(deleterious_burden_results_microglia_vs_premac_signif$gene_name))

ORA_results <- list()
ps <- vector()
ORs <- vector()
intersections <- vector()
ORA_results[[universe_name]] <- ORA(universe, signif_group, gene_sets)
ps <- as.data.frame(rbind(ps, ORA_results[[universe_name]]$p))
ORs <- as.data.frame(rbind(ORs, ORA_results[[universe_name]]$OR))
intersections <- as.data.frame(rbind(intersections, ORA_results[[universe_name]]$intersect))
rownames(ps) <- rownames(ORs) <-  rownames(intersections) <- universe_name

ps <- t(ps)
ORs <- t(ORs)
intersections <- t(intersections)
ORA_results_Mann_Whitney_tests <- list("p" = ps, "OR" = ORs, "intersection_size" = intersections)
save(ORA_results_Mann_Whitney_tests, file = paste0(outdir, "/01_ORA", "/ORA_results_enrichment_in_macrophage_Mann_Whitney_test_genes.Rdata"))


h_row_title <- "Sets of proliferation-associated genes"
h <- ComplexHeatmap::Heatmap(-log10(as.matrix(ps)), 
                             row_title = h_row_title, 
                             row_title_gp = gpar(fontsize = 10),
                             name = "-log10(pval)", 
                             col = colorRamp2(c(0, 0.5, 1, 1.5, 2), c("white", "lavenderblush2", "lightpink", "salmon", "salmon4")), 
                             border = T, 
                             rect_gp = gpar(col = "black", lwd = 1),
                             row_names_gp = gpar(fontsize = 8.5), 
                             column_names_gp = gpar(fontsize = 6.5), 
                             row_names_max_width = unit(10, "cm"),
                             cluster_columns = F, 
                             cluster_rows = F, 
                             column_labels = paste0("Mann-Whitney macrophage", "\n", "survival genes (n = ", length(signif_group), " sig)"),
                             cell_fun = function(j, i, x, y, width, height, fill) {
                               if(log10(as.matrix(ps))[i,j] < log10(0.05)){
                                 grid.text("*", x, y, gp = gpar(fontsize = 12.5, col = "yellow"), just = c(0.5,-1.2))
                                 grid.text(paste0(signif(ORs[i, j], digits = 3), "\n", "(", intersections[i,j], ")"), 
                                           x, y, gp = gpar(fontsize = 8.5)) 
                               }
                               else{
                                 grid.text(paste0(signif(ORs[i, j], digits = 3), "\n", "(", intersections[i,j], ")"), x, y, gp = gpar(fontsize = 8.5))
                               }
                             })

pdf(file = paste0(plotdir, "/01_ORA", "/ORA_results_enrichment_in_macrophage_Mann_Whitney_test_genes", ".pdf"), height = 8, width = 5)
h
dev.off()




## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
##                   4.2.2 Gene Set Enrichment Analysis (GSEA)
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

###################  Enrichment among ranked SKAT-O genes  #####################

# ====================== SKAT-O for deleterious variants =======================

# -- - SKAT-O genes ranked for association with prolif. from iPSC to premac -- -

## Rank genes by decreasing -log10(p)
deleterious_SKATO_results_premac_vs_iPSC$metric <- -log10(deleterious_SKATO_results_premac_vs_iPSC$resampling_pval) 
ranked_deleterious_SKATO_results_premac_vs_iPSC <- deleterious_SKATO_results_premac_vs_iPSC[order(deleterious_SKATO_results_premac_vs_iPSC$metric, decreasing = T), ]

## Remove duplicated gene (keep most signif)
dup_gene <- ranked_deleterious_SKATO_results_premac_vs_iPSC[which(duplicated(ranked_deleterious_SKATO_results_premac_vs_iPSC$gene_name)), "gene_name"]
dup_rows <- which(ranked_deleterious_SKATO_results_premac_vs_iPSC$gene_name == as.character(dup_gene))
# [1]  351 2072
ranked_deleterious_SKATO_results_premac_vs_iPSC <- ranked_deleterious_SKATO_results_premac_vs_iPSC[-dup_rows[2], ]

ranked_genes <- ranked_deleterious_SKATO_results_premac_vs_iPSC$metric
names(ranked_genes) <- ranked_deleterious_SKATO_results_premac_vs_iPSC$gene_name


## Run fast GSEA
set.seed(12112024)
fgsea_del_SKATO_premac_vs_iPSC <- fgsea(pathways = gene_sets, 
                                         stats    = ranked_genes,
                                         minSize  = 1, # min size of set
                                         maxSize  = Inf, # max size of set
                                         eps = 0) # set to 0 for accurate pval calculation 

#                                        pathway      pval      padj    log2err        ES      NES  size  leadingEdge
#                                         <char>     <num>     <num>      <num>     <num>    <num> <int>       <list>
# 1:           Covarrubias_MAGeCK_macro_survival 0.4505495 0.4505495 0.05039643 0.4295551 1.015627    55 RABGGTA,....
# 2:     Covarrubias_Mann_Whitney_macro_survival 0.3526474 0.4231768 0.06184060 0.4004733 1.018693   545 AURKA, B....
# 3: Covarrubias_Mann_Whitney_macro_survival_neg 0.3156843 0.4231768 0.06720651 0.4019289 1.022565   509 AURKA, N....
# 4: Covarrubias_Mann_Whitney_macro_survival_pos 0.2997003 0.4231768 0.06977925 0.4712099 1.086009    36 BCOR, CY....
# 5:                  Covarrubias_macro_survival 0.3356643 0.4231768 0.06421409 0.4009615 1.020107   547 AURKA, B....
# 6:                      DEPMAP_essential_genes 0.2127872 0.4231768 0.08783126 0.4244503 1.062078   204 NOB1, SU....

plotGseaTable(pathways = gene_sets, 
              stats = ranked_genes, 
              fgseaRes = fgsea_del_SKATO_premac_vs_iPSC)
ggsave(filename = paste0(plotdir, "/02_GSEA", "/Results_fgsea_", "del_SKATO_premac_vs_iPSC", ".pdf"), height = 4, width = 10)

plots <- list()
for(i in 1:length(gene_sets)){
  plots[[i]] <- plotEnrichment(gene_sets[[i]], ranked_genes) + labs(title=names(gene_sets)[i]) 
}
plot_grid(plotlist = plots, ncol = 2)
ggsave(filename = paste0(plotdir, "/02_GSEA", "/ES_plots_", "del_SKATO_premac_vs_iPSC", ".pdf"), height = 18, width = 17)


# -- SKAT-O genes ranked for association with prolif. from young to old premac --

## Rank genes by decreasing -log10(p)
deleterious_SKATO_results_old_vs_young_premac$metric <- -log10(deleterious_SKATO_results_old_vs_young_premac$resampling_pval) 
ranked_deleterious_SKATO_results_old_vs_young_premac <- deleterious_SKATO_results_old_vs_young_premac[order(deleterious_SKATO_results_old_vs_young_premac$metric, decreasing = T), ]

## Remove duplicated gene (keep most signif)
dup_gene <- ranked_deleterious_SKATO_results_old_vs_young_premac[which(duplicated(ranked_deleterious_SKATO_results_old_vs_young_premac$gene_name)), "gene_name"]
dup_rows <- which(ranked_deleterious_SKATO_results_old_vs_young_premac$gene_name == as.character(dup_gene))
# [1]  2542 2583
ranked_deleterious_SKATO_results_old_vs_young_premac <- ranked_deleterious_SKATO_results_old_vs_young_premac[-dup_rows[2], ]

ranked_genes <- ranked_deleterious_SKATO_results_old_vs_young_premac$metric
names(ranked_genes) <- ranked_deleterious_SKATO_results_old_vs_young_premac$gene_name

fgsea_del_SKATO_old_vs_young_premac <- fgsea(pathways = gene_sets, 
                                        stats    = ranked_genes,
                                        minSize  = 1, # min size of set
                                        maxSize  = Inf, # max size of set
                                        eps = 0) # set to 0 for accurate pval calculation 

#                                        pathway      pval      padj    log2err        ES      NES  size  leadingEdge
#                                         <char>     <num>     <num>      <num>     <num>    <num> <int>       <list>
# 1:           Covarrubias_MAGeCK_macro_survival 0.1038961 0.4183816 0.13427345 0.4971476 1.178281    55 MED12, N....
# 2:     Covarrubias_Mann_Whitney_macro_survival 0.3486513 0.4183816 0.06238615 0.4043715 1.019695   537 SLC37A4,....
# 3: Covarrubias_Mann_Whitney_macro_survival_neg 0.2907093 0.4183816 0.07130530 0.4096385 1.030681   502 SLC37A4,....
# 4: Covarrubias_Mann_Whitney_macro_survival_pos 0.4255744 0.4255744 0.05302125 0.4512532 1.037819    35 SLFN13, ....
# 5:                  Covarrubias_macro_survival 0.3106893 0.4183816 0.06799226 0.4069248 1.026330   539 SLC37A4,....
# 6:                      DEPMAP_essential_genes 0.3046953 0.4183816 0.06895674 0.4193029 1.040171   201 BRCA1, M....

plotGseaTable(pathways = gene_sets, 
              stats = ranked_genes, 
              fgseaRes = fgsea_del_SKATO_old_vs_young_premac)
ggsave(filename = paste0(plotdir, "/02_GSEA", "/Results_fgsea_", "del_SKATO_old_vs_young_premac", ".pdf"), height = 4, width = 10)

plots <- list()
for(i in 1:length(gene_sets)){
  plots[[i]] <- plotEnrichment(gene_sets[[i]], ranked_genes) + labs(title=names(gene_sets)[i]) 
}
plot_grid(plotlist = plots, ncol = 2)
ggsave(filename = paste0(plotdir, "/02_GSEA", "/ES_plots_", "del_SKATO_old_vs_young_premac", ".pdf"), height = 18, width = 17)


# -- SKAT-O genes ranked for association with prolif. from premac to microglia --

## Rank genes by decreasing -log10(p)
deleterious_SKATO_results_microglia_vs_premac$metric <- -log10(deleterious_SKATO_results_microglia_vs_premac$resampling_pval) 
ranked_deleterious_SKATO_results_microglia_vs_premac <- deleterious_SKATO_results_microglia_vs_premac[order(deleterious_SKATO_results_microglia_vs_premac$metric, decreasing = T), ]

## Remove duplicated gene (keep most signif)
dup_gene <- ranked_deleterious_SKATO_results_microglia_vs_premac[which(duplicated(ranked_deleterious_SKATO_results_microglia_vs_premac$gene_name)), "gene_name"]
dup_rows <- which(ranked_deleterious_SKATO_results_microglia_vs_premac$gene_name == as.character(dup_gene))
# [1]  53 3145
ranked_deleterious_SKATO_results_microglia_vs_premac <- ranked_deleterious_SKATO_results_microglia_vs_premac[-dup_rows[2], ]

ranked_genes <- ranked_deleterious_SKATO_results_microglia_vs_premac$metric
names(ranked_genes) <- ranked_deleterious_SKATO_results_microglia_vs_premac$gene_name

fgsea_del_SKATO_microglia_vs_premac <- fgsea(pathways = gene_sets, 
                                             stats    = ranked_genes,
                                             minSize  = 1, # min size of set
                                             maxSize  = Inf, # max size of set
                                             eps = 0) # set to 0 for accurate pval calculation 

#                                         pathway      pval      padj    log2err        ES       NES  size  leadingEdge
#                                         <char>     <num>     <num>      <num>     <num>     <num> <int>       <list>
# 1:           Covarrubias_MAGeCK_macro_survival 0.1548452 0.2181818 0.10672988 0.4930035 1.1293909    55 BIRC6, N....
# 2:     Covarrubias_Mann_Whitney_macro_survival 0.1398601 0.2181818 0.11331291 0.4257319 1.0506029   546 BIRC6, I....
# 3: Covarrubias_Mann_Whitney_macro_survival_neg 0.1498501 0.2181818 0.10882013 0.4260832 1.0507229   509 BIRC6, I....
# 4: Covarrubias_Mann_Whitney_macro_survival_pos 0.5814186 0.5814186 0.03871667 0.4360675 0.9727014    37 CYFIP1, ....
# 5:                  Covarrubias_macro_survival 0.1148851 0.2181818 0.12687573 0.4276398 1.0554207   548 BIRC6, I....
# 6:                      DEPMAP_essential_genes 0.1818182 0.2181818 0.09688777 0.4395813 1.0645455   204 IWS1, PN....

plotGseaTable(pathways = gene_sets, 
              stats = ranked_genes, 
              fgseaRes = fgsea_del_SKATO_microglia_vs_premac)
ggsave(filename = paste0(plotdir, "/02_GSEA", "/Results_fgsea_", "del_SKATO_microglia_vs_premac", ".pdf"), height = 4, width = 10)

plots <- list()
for(i in 1:length(gene_sets)){
  plots[[i]] <- plotEnrichment(gene_sets[[i]], ranked_genes) + labs(title=names(gene_sets)[i]) 
}
plot_grid(plotlist = plots, ncol = 2)
ggsave(filename = paste0(plotdir, "/02_GSEA", "/ES_plots_", "del_SKATO_microglia_vs_premac", ".pdf"), height = 18, width = 17)



#################  Enrichment among ranked Burden test genes  ##################

# =================== Burden test for deleterious variants =====================

# -- -- Burden test genes ranked for association with prolif. from iPSC to premac -- -- 

## Rank genes by decreasing -log10(p_Bonf) x coef
deleterious_burden_results_premac_vs_iPSC$metric <- -log10(deleterious_burden_results_premac_vs_iPSC$p_Bonf) * deleterious_burden_results_premac_vs_iPSC$coef
ranked_deleterious_burden_results_premac_vs_iPSC <- deleterious_burden_results_premac_vs_iPSC[order(deleterious_burden_results_premac_vs_iPSC$metric, decreasing = T), ]

## No duplicated genes
which(duplicated(ranked_deleterious_burden_results_premac_vs_iPSC$gene_name))
# integer(0)
ranked_genes <- ranked_deleterious_burden_results_premac_vs_iPSC$metric
names(ranked_genes) <- ranked_deleterious_burden_results_premac_vs_iPSC$gene_name

fgsea_del_Burden_premac_vs_iPSC <- fgsea(pathways = gene_sets, 
                                        stats    = ranked_genes,
                                        minSize  = 1, # min size of set
                                        maxSize  = Inf, # max size of set
                                        eps = 0) # set to 0 for accurate pval calculation 

#                                        pathway      pval      padj    log2err         ES        NES  size  leadingEdge
#                                         <char>     <num>     <num>      <num>      <num>      <num> <int>       <list>
# 1:           Covarrubias_MAGeCK_macro_survival 0.5121951 0.6386986 0.07455008 -0.5348837 -0.8338363     3 YARS2, R....
# 2:     Covarrubias_Mann_Whitney_macro_survival 0.6386986 0.6386986 0.05922192 -0.8092370 -0.9734414    25         BCOR
# 3: Covarrubias_Mann_Whitney_macro_survival_neg 0.2035176 0.6105528 0.15419097  0.9717069  1.2570636    20 PTPRD, AURKA
# 4: Covarrubias_Mann_Whitney_macro_survival_pos 0.1131768 0.6105528 0.28780513 -0.9921260 -1.4728024     5         BCOR
# 5:                  Covarrubias_macro_survival 0.6386986 0.6386986 0.05922192 -0.8092370 -0.9734414    25         BCOR
# 6:                      DEPMAP_essential_genes 0.5965217 0.6386986 0.06321912 -0.9256198 -1.1737746    11  TAF5, YARS2
# 

plotGseaTable(pathways = gene_sets, 
              stats = ranked_genes, 
              fgseaRes = fgsea_del_Burden_premac_vs_iPSC)
ggsave(filename = paste0(plotdir, "/02_GSEA", "/Results_fgsea_", "del_Burden_premac_vs_iPSC", ".pdf"), height = 4, width = 10)

plots <- list()
for(i in 1:length(gene_sets)){
  plots[[i]] <- plotEnrichment(gene_sets[[i]], ranked_genes) + labs(title=names(gene_sets)[i]) 
}
plot_grid(plotlist = plots, ncol = 2)
ggsave(filename = paste0(plotdir, "/02_GSEA", "/ES_plots_", "del_Burden_premac_vs_iPSC", ".pdf"), height = 18, width = 17)


# -- -- Burden test genes ranked for association with prolif. from young to old premac -- -- 

deleterious_burden_results_old_vs_young_premac$metric <- -log10(deleterious_burden_results_old_vs_young_premac$p_Bonf) * deleterious_burden_results_old_vs_young_premac$coef
ranked_deleterious_burden_results_old_vs_young_premac <- deleterious_burden_results_old_vs_young_premac[order(deleterious_burden_results_old_vs_young_premac$metric, decreasing = T), ]

which(duplicated(ranked_deleterious_burden_results_old_vs_young_premac$gene_name))
# integer(0)
ranked_genes <- ranked_deleterious_burden_results_old_vs_young_premac$metric
names(ranked_genes) <- ranked_deleterious_burden_results_old_vs_young_premac$gene_name

fgsea_del_Burden_old_vs_young_premac <- fgsea(pathways = gene_sets, 
                                         stats    = ranked_genes,
                                         minSize  = 1, # min size of set
                                         maxSize  = Inf, # max size of set
                                         eps = 0) # set to 0 for accurate pval calculation 

#                                        pathway      pval      padj    log2err         ES        NES  size  leadingEdge
#                                         <char>     <num>     <num>      <num>      <num>      <num> <int>       <list>
# 1:           Covarrubias_MAGeCK_macro_survival 0.7009524 0.8411429 0.05986031  0.2857143  0.5907423     6 ALG8, GP....
# 2:     Covarrubias_Mann_Whitney_macro_survival 0.3325301 0.5000000 0.11426650 -0.9856115 -1.3427975    28         AMTN
# 3: Covarrubias_Mann_Whitney_macro_survival_neg 0.3202765 0.5000000 0.11378726 -0.9857143 -1.3730539    27         AMTN
# 4: Covarrubias_Mann_Whitney_macro_survival_pos 0.3320080 0.5000000 0.10208011 -0.8373494 -1.1323553     1       SLFN13
# 5:                  Covarrubias_macro_survival 0.3333333 0.5000000 0.11284336 -0.9855072 -1.3272038    29         AMTN
# 6:                      DEPMAP_essential_genes 0.9123506 0.9123506 0.04969014  0.1587302  0.2497668    14 ALG8, AT....

plotGseaTable(pathways = gene_sets, 
              stats = ranked_genes, 
              fgseaRes = fgsea_del_Burden_old_vs_young_premac)
ggsave(filename = paste0(plotdir, "/02_GSEA", "/Results_fgsea_", "del_Burden_old_vs_young_premac", ".pdf"), height = 4, width = 10)

plots <- list()
for(i in 1:length(gene_sets)){
  plots[[i]] <- plotEnrichment(gene_sets[[i]], ranked_genes) + labs(title=names(gene_sets)[i]) 
}
plot_grid(plotlist = plots, ncol = 2)
ggsave(filename = paste0(plotdir, "/02_GSEA", "/ES_plots_", "del_Burden_old_vs_young_premac", ".pdf"), height = 18, width = 17)


# -- -- Burden test genes ranked for association with prolif. from premac to microglia -- -- 

deleterious_burden_results_microglia_vs_premac$metric <- -log10(deleterious_burden_results_microglia_vs_premac$p_Bonf) * deleterious_burden_results_microglia_vs_premac$coef
ranked_deleterious_burden_results_microglia_vs_premac <- deleterious_burden_results_microglia_vs_premac[order(deleterious_burden_results_microglia_vs_premac$metric, decreasing = T), ]

which(duplicated(ranked_deleterious_burden_results_microglia_vs_premac$gene_name))
# integer(0)
ranked_genes <- ranked_deleterious_burden_results_microglia_vs_premac$metric
names(ranked_genes) <- ranked_deleterious_burden_results_microglia_vs_premac$gene_name

fgsea_del_Burden_microglia_vs_premac <- fgsea(pathways = gene_sets, 
                                              stats    = ranked_genes,
                                              minSize  = 1, # min size of set
                                              maxSize  = Inf, # max size of set
                                              eps = 0) # set to 0 for accurate pval calculation 

#                                        pathway      pval      padj   log2err         ES        NES  size  leadingEdge
#                                         <char>     <num>     <num>     <num>      <num>      <num> <int>       <list>
# 1:           Covarrubias_MAGeCK_macro_survival 0.2515091 0.3018109 0.1209851  0.9405405  1.3371538     3 PGM3, DNAJA3
# 2:     Covarrubias_Mann_Whitney_macro_survival 0.1131387 0.2484663 0.1429011  0.9356023  1.2208892    29 ATP11A, ....
# 3: Covarrubias_Mann_Whitney_macro_survival_neg 0.1648079 0.2484663 0.1167392  0.9277887  1.1949885    27 ATP11A, ....
# 4: Covarrubias_Mann_Whitney_macro_survival_pos 0.1656442 0.2484663 0.1541910  0.9516129  1.4100114     2 INPP5D, ....
# 5:                  Covarrubias_macro_survival 0.1131387 0.2484663 0.1429011  0.9356023  1.2208892    29 ATP11A, ....
# 6:                      DEPMAP_essential_genes 0.6870027 0.6870027 0.0772747 -0.5211507 -0.7376094    10        GSPT1

plotGseaTable(pathways = gene_sets, 
              stats = ranked_genes, 
              fgseaRes = fgsea_del_Burden_microglia_vs_premac)
ggsave(filename = paste0(plotdir, "/02_GSEA", "/Results_fgsea_", "del_Burden_microglia_vs_premac", ".pdf"), height = 4, width = 10)

plots <- list()
for(i in 1:length(gene_sets)){
  plots[[i]] <- plotEnrichment(gene_sets[[i]], ranked_genes) + labs(title=names(gene_sets)[i]) 
}
plot_grid(plotlist = plots, ncol = 2)
ggsave(filename = paste0(plotdir, "/02_GSEA", "/ES_plots_", "del_Burden_microglia_vs_premac", ".pdf"), height = 18, width = 17)



##############  Enrichment among ranked Covarrubias MAGeCK genes  ##############

## Rank genes by -log(adj. p)
essential_genes_Covarrubias_2020_MAGeCK$metric <- -log10(essential_genes_Covarrubias_2020_MAGeCK$fdr)
ranked_essential_genes_Covarrubias_2020_MAGeCK<- essential_genes_Covarrubias_2020_MAGeCK[order(essential_genes_Covarrubias_2020_MAGeCK$metric, decreasing = T), ]

## Unique human genes
which(duplicated(ranked_essential_genes_Covarrubias_2020_MAGeCK$hsapiens_homolog_associated_gene_name))
# integer(0)

ranked_genes <- ranked_essential_genes_Covarrubias_2020_MAGeCK$metric
names(ranked_genes) <- ranked_essential_genes_Covarrubias_2020_MAGeCK$hsapiens_homolog_associated_gene_name

## Gene sets
gene_sets <- list("SKAT-O (Del) genes for any comparison" = unique(deleterious_SKATO_results_signif$gene_name),
                  "SKAT-O (Del) genes for Premac vs iPSC" = unique(deleterious_SKATO_results_premac_vs_iPSC_signif$gene_name),
                  "SKAT-O (Del) genes for Old vs Young Premac" = unique(deleterious_SKATO_results_old_vs_young_premac_signif$gene_name),
                  "SKAT-O (Del) genes for Microglia vs Premac" = unique(deleterious_SKATO_results_microglia_vs_premac_signif$gene_name), 
                  "Burden test (Del) genes for any comparison" = unique(deleterious_burden_results_signif$gene_name),
                  "Burden test (Del) genes for Premac vs iPSC" = unique(deleterious_burden_results_premac_vs_iPSC_signif$gene_name),
                  "Burden test (Del) genes for Old vs Young Premac" = unique(deleterious_burden_results_old_vs_young_premac_signif$gene_name),
                  "Burden test (Del) genes for Microglia vs Premac" = unique(deleterious_burden_results_microglia_vs_premac_signif$gene_name))

fgsea_Covarrubias_MAGeCK_genes <- fgsea(pathways = gene_sets, 
                                              stats    = ranked_genes,
                                              minSize  = 1, 
                                              maxSize  = Inf, 
                                              eps = 0, 
                                              scoreType = "pos") 

#                                            pathway      pval     padj    log2err          ES         NES  size  leadingEdge
#                                             <char>     <num>    <num>      <num>       <num>       <num> <int>       <list>
# 1: Burden test (Del) genes for Microglia vs Premac 0.8781219 0.996004 0.01699711 0.107142857 0.220913108     1         PGM3
# 2:      Burden test (Del) genes for any comparison 0.8781219 0.996004 0.01699711 0.107142857 0.220913108     1         PGM3
# 3:      SKAT-O (Del) genes for Microglia vs Premac 0.7092907 0.996004 0.02921031 0.231590116 0.721454320     5 DNAJA3, ....
# 4:      SKAT-O (Del) genes for Old vs Young Premac 0.6933067 0.996004 0.03034673 0.233022944 0.762881602     6 NUBP1, H....
# 5:           SKAT-O (Del) genes for Premac vs iPSC 0.9960040 0.996004 0.00288973 0.002570694 0.007615778     4 GEMIN4, ....
# 6:           SKAT-O (Del) genes for any comparison 0.9500500 0.996004 0.01046104 0.069027654 0.282841204    13 GEMIN4, ....


plotGseaTable(pathways = gene_sets, 
              stats = ranked_genes, 
              fgseaRes = fgsea_Covarrubias_MAGeCK_genes)
ggsave(filename = paste0(plotdir, "/02_GSEA", "/Results_fgsea_", "Covarrubias_MAGeCK", ".pdf"), height = 4, width = 10)

plots <- list()
for(i in 1:length(gene_sets)){
  plots[[i]] <- plotEnrichment(gene_sets[[i]], ranked_genes) + labs(title=names(gene_sets)[i]) 
}
plot_grid(plotlist = plots, ncol = 2)
ggsave(filename = paste0(plotdir, "/02_GSEA", "/ES_plots_", "Covarrubias_MAGeCK", ".pdf"), height = 24, width = 17)




##############  Enrichment among ranked Covarrubias Mann-Whitney genes  ##############

## Rank genes by -log(p) x coef
essential_genes_Covarrubias_2020_Mann_Whitney$metric <- -log10(essential_genes_Covarrubias_2020_Mann_Whitney$`p-value`) * essential_genes_Covarrubias_2020_Mann_Whitney$`MW phenotype`
ranked_essential_genes_Covarrubias_2020_Mann_Whitney <- essential_genes_Covarrubias_2020_Mann_Whitney[order(essential_genes_Covarrubias_2020_Mann_Whitney$metric, decreasing = T), ]

## Unique human genes
which(duplicated(ranked_essential_genes_Covarrubias_2020_Mann_Whitney$hsapiens_homolog_associated_gene_name))
# integer(0)

ranked_genes <- ranked_essential_genes_Covarrubias_2020_Mann_Whitney$metric
names(ranked_genes) <- ranked_essential_genes_Covarrubias_2020_Mann_Whitney$hsapiens_homolog_associated_gene_name

fgsea_Covarrubias_Mann_Whitney_genes <- fgsea(pathways = gene_sets, 
                                        stats    = ranked_genes,
                                        minSize  = 1, 
                                        maxSize  = Inf, 
                                        eps = 0) 

#                                            pathway       pval      padj    log2err         ES        NES  size  leadingEdge
#                                             <char>      <num>     <num>      <num>      <num>      <num> <int>       <list>
# 1: Burden test (Del) genes for Microglia vs Premac 0.53039832 0.7432567 0.04513442 -0.6686144 -1.0208700    14 KANSL2, PGM3
# 2: Burden test (Del) genes for Old vs Young Premac 0.31416838 0.7432567 0.10755438  0.8494372  1.1261838     1        NEK11
# 3:      Burden test (Del) genes for Premac vs iPSC 0.05101872 0.4081497 0.32177592  0.9592740  1.3906083     2 BCOR, LRRC55
# 4:      Burden test (Del) genes for any comparison 0.69596691 0.7432567 0.03185347 -0.6038868 -0.9269830    17 KANSL2, PGM3
# 5:      SKAT-O (Del) genes for Microglia vs Premac 0.74325674 0.7432567 0.02681588 -0.5948426 -0.9644108   219 DNAJA3, ....
# 6:      SKAT-O (Del) genes for Old vs Young Premac 0.27272727 0.7432567 0.07455008 -0.6427954 -1.0410376   199 ZFR, MED....
# 7:           SKAT-O (Del) genes for Premac vs iPSC 0.64035964 0.7432567 0.03419460 -0.6080056 -0.9817346   153 FARS2, G....
# 8:           SKAT-O (Del) genes for any comparison 0.50849151 0.7432567 0.04486451 -0.6134024 -1.0020137   517 ZFR, MED....

plotGseaTable(pathways = gene_sets, 
              stats = ranked_genes, 
              fgseaRes = fgsea_Covarrubias_Mann_Whitney_genes)
ggsave(filename = paste0(plotdir, "/02_GSEA", "/Results_fgsea_", "Covarrubias_Mann_Whitney", ".pdf"), height = 4, width = 10)

plots <- list()
for(i in 1:length(gene_sets)){
  plots[[i]] <- plotEnrichment(gene_sets[[i]], ranked_genes) + labs(title=names(gene_sets)[i]) 
}
plot_grid(plotlist = plots, ncol = 2)
ggsave(filename = paste0(plotdir, "/02_GSEA", "/ES_plots_", "Covarrubias_Mann_Whitney", ".pdf"), height = 24, width = 17)







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
# date     2024-12-11
# rstudio  2024.04.0+735 Chocolate Cosmos (server)
# pandoc   3.1.12.3 @ /opt/view/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package          * version     date (UTC) lib source
# AnnotationDbi    * 1.64.1      2023-11-03 [1] Bioconductor
# Biobase          * 2.62.0      2023-10-24 [1] Bioconductor
# BiocFileCache      2.10.2      2024-03-27 [1] Bioconductor 3.18 (R 4.3.1)
# BiocGenerics     * 0.48.1      2023-11-01 [1] Bioconductor
# BiocParallel       1.36.0      2023-10-24 [1] Bioconductor
# biomaRt          * 2.58.2      2024-01-30 [1] Bioconductor 3.18 (R 4.3.1)
# Biostrings         2.70.3      2024-03-13 [1] Bioconductor 3.18 (R 4.3.1)
# bit                4.0.5       2022-11-15 [1] CRAN (R 4.3.1)
# bit64              4.0.5       2020-08-30 [1] CRAN (R 4.3.1)
# bitops             1.0-7       2021-04-24 [1] CRAN (R 4.3.1)
# blob               1.2.4       2023-03-17 [1] CRAN (R 4.3.1)
# cachem             1.0.8       2023-05-01 [1] CRAN (R 4.3.1)
# Cairo              1.6-2       2023-11-28 [1] CRAN (R 4.3.1)
# circlize         * 0.4.16      2024-02-20 [1] CRAN (R 4.3.1)
# cli                3.6.2       2023-12-11 [1] CRAN (R 4.3.1)
# clue               0.3-65      2023-09-23 [1] CRAN (R 4.3.1)
# cluster            2.1.6       2023-12-01 [1] CRAN (R 4.3.1)
# codetools          0.2-20      2024-03-31 [1] CRAN (R 4.3.1)
# colorspace         2.1-0       2023-01-23 [1] CRAN (R 4.3.1)
# ComplexHeatmap   * 2.18.0      2023-10-24 [1] Bioconductor
# cowplot          * 1.1.3       2024-01-22 [1] CRAN (R 4.3.1)
# crayon             1.5.2       2022-09-29 [1] CRAN (R 4.3.1)
# curl               5.2.1       2024-03-01 [1] CRAN (R 4.3.1)
# data.table         1.15.4      2024-03-30 [1] CRAN (R 4.3.1)
# DBI                1.2.2       2024-02-16 [1] CRAN (R 4.3.1)
# dbplyr             2.5.0       2024-03-19 [1] CRAN (R 4.3.1)
# digest             0.6.35      2024-03-11 [1] CRAN (R 4.3.1)
# doParallel         1.0.17      2022-02-07 [1] CRAN (R 4.3.1)
# dplyr            * 1.1.4       2023-11-17 [1] CRAN (R 4.3.1)
# fansi              1.0.6       2023-12-08 [1] CRAN (R 4.3.1)
# farver             2.1.1       2022-07-06 [1] CRAN (R 4.3.1)
# fastmap            1.1.1       2023-02-24 [1] CRAN (R 4.3.1)
# fastmatch          1.1-4       2023-08-18 [1] CRAN (R 4.3.1)
# fgsea            * 1.28.0      2023-10-24 [1] Bioconductor
# filelock           1.0.3       2023-12-11 [1] CRAN (R 4.3.1)
# forcats          * 1.0.0       2023-01-29 [1] CRAN (R 4.3.1)
# foreach            1.5.2       2022-02-02 [1] CRAN (R 4.3.1)
# generics           0.1.3       2022-07-05 [1] CRAN (R 4.3.1)
# GenomeInfoDb       1.38.8      2024-03-15 [1] Bioconductor 3.18 (R 4.3.1)
# GenomeInfoDbData   1.2.11      2024-11-06 [1] Bioconductor
# GetoptLong         1.0.5       2020-12-15 [1] CRAN (R 4.3.1)
# ggplot2          * 3.5.1       2024-04-23 [1] CRAN (R 4.3.1)
# ggrepel            0.9.5       2024-01-10 [1] CRAN (R 4.3.1)
# GlobalOptions      0.1.2       2020-06-10 [1] CRAN (R 4.3.1)
# glue               1.7.0       2024-01-09 [1] CRAN (R 4.3.1)
# gtable             0.3.4       2023-08-21 [1] CRAN (R 4.3.1)
# hms                1.1.3       2023-03-21 [1] CRAN (R 4.3.1)
# httr               1.4.7       2023-08-15 [1] CRAN (R 4.3.1)
# IRanges          * 2.36.0      2023-10-24 [1] Bioconductor
# iterators          1.0.14      2022-02-05 [1] CRAN (R 4.3.1)
# jtools             2.2.2       2023-07-11 [1] CRAN (R 4.3.1)
# KEGGREST           1.42.0      2023-10-24 [1] Bioconductor
# labeling           0.4.3       2023-08-29 [1] CRAN (R 4.3.1)
# lattice            0.22-6      2024-03-20 [1] CRAN (R 4.3.1)
# lifecycle          1.0.4       2023-11-07 [1] CRAN (R 4.3.1)
# lubridate        * 1.9.3       2023-09-27 [1] CRAN (R 4.3.1)
# magick             2.8.3       2024-02-18 [1] CRAN (R 4.3.1)
# magrittr           2.0.3       2022-03-30 [1] CRAN (R 4.3.1)
# Matrix             1.6-5       2024-01-11 [1] CRAN (R 4.3.1)
# matrixStats        1.2.0       2023-12-11 [1] CRAN (R 4.3.1)
# memoise            2.0.1       2021-11-26 [1] CRAN (R 4.3.1)
# munsell            0.5.1       2024-04-01 [1] CRAN (R 4.3.1)
# org.Hs.eg.db     * 3.18.0      2024-11-06 [1] Bioconductor
# pander             0.6.5       2022-03-18 [1] CRAN (R 4.3.1)
# pillar             1.9.0       2023-03-22 [1] CRAN (R 4.3.1)
# pkgconfig          2.0.3       2019-09-22 [1] CRAN (R 4.3.1)
# png                0.1-8       2022-11-29 [1] CRAN (R 4.3.1)
# prettyunits        1.2.0       2023-09-24 [1] CRAN (R 4.3.1)
# progress           1.2.3       2023-12-06 [1] CRAN (R 4.3.1)
# purrr            * 1.0.2       2023-08-10 [1] CRAN (R 4.3.1)
# R6                 2.5.1       2021-08-19 [1] CRAN (R 4.3.1)
# ragg               1.3.0       2024-03-13 [1] CRAN (R 4.3.1)
# rappdirs           0.3.3       2021-01-31 [1] CRAN (R 4.3.1)
# RColorBrewer       1.1-3       2022-04-03 [1] CRAN (R 4.3.1)
# Rcpp               1.0.12      2024-01-09 [1] CRAN (R 4.3.1)
# RCurl              1.98-1.14   2024-01-09 [1] CRAN (R 4.3.1)
# readr            * 2.1.5       2024-01-10 [1] CRAN (R 4.3.1)
# rjson              0.2.21      2022-01-09 [1] CRAN (R 4.3.1)
# rlang            * 1.1.3       2024-01-10 [1] CRAN (R 4.3.1)
# RSQLite            2.3.6       2024-03-31 [1] CRAN (R 4.3.1)
# rstudioapi         0.16.0      2024-03-24 [1] CRAN (R 4.3.1)
# S4Vectors        * 0.40.2      2023-11-23 [1] Bioconductor 3.18 (R 4.3.1)
# scales             1.3.0       2023-11-28 [1] CRAN (R 4.3.1)
# sessioninfo      * 1.2.2       2021-12-06 [1] CRAN (R 4.3.1)
# shape              1.4.6.1     2024-02-23 [1] CRAN (R 4.3.1)
# stringi            1.8.3       2023-12-11 [1] CRAN (R 4.3.1)
# stringr          * 1.5.1       2023-11-14 [1] CRAN (R 4.3.1)
# systemfonts        1.0.6       2024-03-07 [1] CRAN (R 4.3.1)
# textshaping        0.3.7       2023-10-09 [1] CRAN (R 4.3.1)
# tibble           * 3.2.1       2023-03-20 [1] CRAN (R 4.3.1)
# tidyr            * 1.3.1       2024-01-24 [1] CRAN (R 4.3.1)
# tidyselect         1.2.1       2024-03-11 [1] CRAN (R 4.3.1)
# tidyverse        * 2.0.0       2023-02-22 [1] CRAN (R 4.3.1)
# timechange         0.3.0       2024-01-18 [1] CRAN (R 4.3.1)
# tzdb               0.4.0       2023-05-12 [1] CRAN (R 4.3.1)
# utf8               1.2.4       2023-10-22 [1] CRAN (R 4.3.1)
# vctrs              0.6.5       2023-12-01 [1] CRAN (R 4.3.1)
# vroom              1.6.5       2023-12-05 [1] CRAN (R 4.3.1)
# withr              3.0.0       2024-01-16 [1] CRAN (R 4.3.1)
# XML                3.99-0.16.1 2024-01-22 [1] CRAN (R 4.3.1)
# xml2               1.3.6       2023-12-04 [1] CRAN (R 4.3.1)
# XVector            0.42.0      2023-10-24 [1] Bioconductor
# zlibbioc           1.48.2      2024-03-13 [1] Bioconductor 3.18 (R 4.3.1)
# 
# [1] /opt/view/rlib/R/library
# [2] /software/hgi/installs/softpack/rstudio/.spack/linux-ubuntu22.04-x86_64_v3/gcc-11.4.0/r-4.3.1-bfwldrk76z6f52upk47zepliekn7ayqz/rlib/R/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────





