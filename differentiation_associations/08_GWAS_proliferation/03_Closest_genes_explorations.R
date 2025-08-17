
library(tidyverse)
library(dplyr)
library(biomaRt)
library(httr)
library(cowplot)
library(sessioninfo)

#-------------------------------------------------------------------------------
#                8.3 Explorations of closest genes to lead SNPs 
#-------------------------------------------------------------------------------
#  Code to explore overlaps between the closest genes to lead SNPs for each 
#  proliferation comparison, with genes associated with proliferation 
#  according to DGE and Burden tests, and OpenTarget trait target genes.
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

## Working dir
setwd("/lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_Daianna/")

## Output data and plot dirs 
outdir = paste(getwd(), "output_data", "08_GWAS_proliferation", "03_Closest_genes_explorations", sep = "/")
input_dir = paste(getwd(), "input_data", "08_GWAS_proliferation", "03_Closest_genes_explorations", sep = "/")
plotdir = paste(getwd(), "plots", "08_GWAS_proliferation", "03_Closest_genes_explorations", sep = "/")

dir.create(outdir, recursive = T)
dir.create(input_dir, recursive = T)
dir.create(plotdir, recursive = T)

## Input dirs
input_dir00 = paste(getwd(), "output_data", "08_GWAS_proliferation", "00_Data_exploration_processing", sep = "/")
input_dir01 <- paste(getwd(), "output_data", "07_Diff_expr_proliferation_permutations", "01_run_proliferation_DGE", sep = "/")
input_dir02 <- paste(getwd(), "output_data", "08_GWAS_proliferation", "02_Check_association_results", sep = "/")


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#                           8.3.0 Data preparation
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

##############  DEGs*  ##############  
# *for proliferation from preMac -> microglia 

load(paste0(input_dir01, "/de_genes_proliferation_IFN.Rdata"), verbose = T)
load(paste0(input_dir01, "/de_genes_proliferation_LPS.Rdata"), verbose = T)
load(paste0(input_dir01, "/de_genes_proliferation_untreated.Rdata"), verbose = T)

de_genes_proliferation_IFN <- de_genes_proliferation_IFN$gene_id
de_genes_proliferation_LPS <- de_genes_proliferation_LPS$gene_id
de_genes_proliferation_untreated <- de_genes_proliferation_untreated$gene_id


##############  GWAS closest genes  ##############  

load(paste0(input_dir02, "/closest_genes_leadSNPs_premac_vs_iPSC.Rdata"), verbose = T)
closest_genes_premac_vs_iPSC <- lapply(gene_tracks_across_chrs_premac_vs_iPSC, 
                                       function(chr)lapply(chr, function(leadSNP){leadSNP[,"gene_id"]})) %>% unlist %>% unique

load(paste0(input_dir02, "/closest_genes_leadSNPs_old_vs_young_premac.Rdata"), verbose = T)
closest_genes_old_vs_young_premac <- lapply(gene_tracks_across_chrs_old_vs_young_premac, 
                                       function(chr)lapply(chr, function(leadSNP){leadSNP[,"gene_id"]})) %>% unlist %>% unique

load(paste0(input_dir02, "/closest_genes_leadSNPs_microglia_vs_premac_IFN.Rdata"), verbose = T)
closest_genes_microglia_vs_premac_IFN <- lapply(gene_tracks_across_chrs_microglia_vs_premac_IFN, 
                                            function(chr)lapply(chr, function(leadSNP){leadSNP[,"gene_id"]})) %>% unlist %>% unique

load(paste0(input_dir02, "/closest_genes_leadSNPs_microglia_vs_premac_LPS.Rdata"), verbose = T)
closest_genes_microglia_vs_premac_LPS <- lapply(gene_tracks_across_chrs_microglia_vs_premac_LPS, 
                                                function(chr)lapply(chr, function(leadSNP){leadSNP[,"gene_id"]})) %>% unlist %>% unique

load(paste0(input_dir02, "/closest_genes_leadSNPs_microglia_vs_premac_untreated.Rdata"), verbose = T)
closest_genes_microglia_vs_premac_untreated <- lapply(gene_tracks_across_chrs_microglia_vs_premac_untreated, 
                                                function(chr)lapply(chr, function(leadSNP){leadSNP[,"gene_id"]})) %>% unlist %>% unique


##############   Burden test signif. genes  ##############  

del_burden_results <- as.data.frame(read_csv("/lustre/scratch123/hgi/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/data/results/2.2.rare_vars_vs_prolif/del_burden_scaled_centered_prop_pvals.csv"))
deleterious_Burden_results_signif <- subset(del_burden_results, p_Bonf < 0.05)

## Divide by comparison
deleterious_Burden_signif_premac_vs_iPSC <- subset(deleterious_Burden_results_signif, comparison == "line_prop_changes_premac_iPSC")$gene_name
deleterious_Burden_signif_old_vs_young_premac <- subset(deleterious_Burden_results_signif, comparison == "line_prop_changes_old_vs_young_premac")$gene_name
deleterious_Burden_signif_microglia_vs_premac <- subset(deleterious_Burden_results_signif, comparison == "line_prop_changes_microglia_premac")$gene_name

## Add Ensembl IDs
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
deleterious_Burden_signif_premac_vs_iPSC_ids <- getBM(values = deleterious_Burden_signif_premac_vs_iPSC,
                                                mart = ensembl,
                                                attributes = c("ensembl_gene_id"),
                                                filters = "external_gene_name") %>% unname() %>% unlist()

deleterious_Burden_signif_old_vs_young_premac_ids <- getBM(values = deleterious_Burden_signif_old_vs_young_premac,
                                                      mart = ensembl,
                                                      attributes = c("ensembl_gene_id"),
                                                      filters = "external_gene_name") %>% unname() %>% unlist()

deleterious_Burden_signif_microglia_vs_premac_ids <- getBM(values = deleterious_Burden_signif_microglia_vs_premac,
                                                           mart = ensembl,
                                                           attributes = c("ensembl_gene_id"),
                                                           filters = "external_gene_name") %>% unname() %>% unlist()


##############   OpenTargets proliferation-related genes   ############## 

## Original instructions for OT data access from: 
##     https://platform-docs.opentargets.org/data-access/graphql-api

# Build query string to extract disease info and targets  
get_dissease_targets <- function(disease_id, n_targets){
  
  query_string = paste0("
    query targetsAssociated($efoId: String!){
      disease(efoId: $efoId) {
        id
        name
        associatedTargets(page: { index: 0, size: ", n_targets, "}) {
          count
          rows {
            target {
              id
              approvedSymbol
            }
            score
          }
        }
      }
    }
  ")
    
  # Set base URL of GraphQL API endpoint
  base_url <- "https://api.platform.opentargets.org/api/v4/graphql"
  
  # Set variables object with arguments to be passed to endpoint
  variables <- list("efoId" = disease_id)
  
  # Construct POST request body object with query string and variables
  post_body <- list(query = query_string, variables = variables)
  # POST request
  r <- POST(url = base_url, body = post_body, encode = 'json')
  
  return(r)
  
}

# ********************
#       Cancer  
# ********************

disease_id <- "MONDO_0004992"
r <- get_dissease_targets(disease_id, 20562)

cancer_genes <- do.call(rbind, lapply(content(r)$data$disease$associatedTargets$rows, function(g){unlist(g)})) %>% as.data.frame()
cancer_genes$score <- as.numeric(cancer_genes$score)
cancer_genes = cancer_genes %>% arrange(-score)
save(cancer_genes, file = paste0(outdir, "/OpenTargets_cancer_genes.Rdata"))

dim(cancer_genes)
# [1] 20562     3
which(duplicated(cancer_genes$target.id))
# integer(0)


# ****************************
#      Cell count traits 
# ****************************

# 1. Complete blood cell count: quantification of blood components

## Extract targets 
disease_id <- "EFO_0004586"
r <- get_dissease_targets(disease_id, 9980)

blood_cell_count_genes <- do.call(rbind, lapply(content(r)$data$disease$associatedTargets$rows, function(g){unlist(g)})) %>% as.data.frame()
blood_cell_count_genes$score <- as.numeric(blood_cell_count_genes$score)
blood_cell_count_genes = blood_cell_count_genes %>% arrange(-score)
save(blood_cell_count_genes, file = paste0(outdir, "/OpenTargets_complete_blood_cell_count_genes.Rdata"))

dim(blood_cell_count_genes)
# [1] 9980     3
which(duplicated(blood_cell_count_genes$target.id))
# integer(0)


# 2. Leukocyte count: number of WHITE BLOOD CELLS per unit volume in venous BLOOD
disease_id <- "EFO_0004308"
r <- get_dissease_targets(disease_id, 7818)

leukocyte_count_genes <- do.call(rbind, lapply(content(r)$data$disease$associatedTargets$rows, function(g){unlist(g)})) %>% as.data.frame()
leukocyte_count_genes$score <- as.numeric(leukocyte_count_genes$score)
leukocyte_count_genes = leukocyte_count_genes %>% arrange(-score)
save(leukocyte_count_genes, file = paste0(outdir, "/OpenTargets_leukocyte_count_genes.Rdata"))

dim(leukocyte_count_genes)
# [1] 7818     3
which(duplicated(leukocyte_count_genes$target.id))
# integer(0)


# 3. Myeloid white cell count: number of myeloid leukocytes in a specified volume of blood
disease_id <- "EFO_0007988"
r <- get_dissease_targets(disease_id, 5679)

myeloid_white_cell_count_genes <- do.call(rbind, lapply(content(r)$data$disease$associatedTargets$rows, function(g){unlist(g)})) %>% as.data.frame()
myeloid_white_cell_count_genes$score <- as.numeric(myeloid_white_cell_count_genes$score)
myeloid_white_cell_count_genes = myeloid_white_cell_count_genes %>% arrange(-score)
save(myeloid_white_cell_count_genes, file = paste0(outdir, "/OpenTargets_myeloid_white_cell_count_genes.Rdata"))

dim(myeloid_white_cell_count_genes)
# [1] 5679     3
which(duplicated(myeloid_white_cell_count_genes$target.id))
# integer(0)


# 4. Monocyte (blood-circulating myeloid leukocytes) count: number of monocytes in the blood
disease_id <- "EFO_0005091"
r <- get_dissease_targets(disease_id, 3240)

monocyte_count_genes <- do.call(rbind, lapply(content(r)$data$disease$associatedTargets$rows, function(g){unlist(g)})) %>% as.data.frame()
monocyte_count_genes$score <- as.numeric(monocyte_count_genes$score)
monocyte_count_genes = monocyte_count_genes %>% arrange(-score)
save(monocyte_count_genes, file = paste0(outdir, "/OpenTargets_monocyte_count_genes.Rdata"))

dim(monocyte_count_genes)
# [1] 3240     3
which(duplicated(monocyte_count_genes$target.id))
# integer(0)


# ************************************
#    Cell population proliferation
# ************************************
## Extract targets 
disease_id <- "GO_0008283"
r <- get_dissease_targets(disease_id, 1)

## Only one gene: TOX (ENSG00000198846)


# ****************************
#       Brain disorders
# ****************************

# 1. Parkison's disease
disease_id <- "MONDO_0005180"
r <- get_dissease_targets(disease_id, 4455)

parkinson_genes <- do.call(rbind, lapply(content(r)$data$disease$associatedTargets$rows, function(g){unlist(g)})) %>% as.data.frame()
parkinson_genes$score <- as.numeric(parkinson_genes$score)
parkinson_genes = parkinson_genes %>% arrange(-score)
save(parkinson_genes, file = paste0(outdir, "/OpenTargets_Parkinson_disease_genes.Rdata"))

dim(parkinson_genes)
# [1] 4455     3
which(duplicated(parkinson_genes$target.id))
# integer(0)


# 2. Alzheimer's disease
disease_id <- "MONDO_0004975"
r <- get_dissease_targets(disease_id, 4626)

alzheimer_genes <- do.call(rbind, lapply(content(r)$data$disease$associatedTargets$rows, function(g){unlist(g)})) %>% as.data.frame()
alzheimer_genes$score <- as.numeric(alzheimer_genes$score)
alzheimer_genes = alzheimer_genes %>% arrange(-score)
save(alzheimer_genes, file = paste0(outdir, "/OpenTargets_Alzheimer_disease_genes.Rdata"))

dim(alzheimer_genes)
# [1] 4626     3
which(duplicated(alzheimer_genes$target.id))
# integer(0)


## Overlaps between all OTAR targets
OTAR_sets <- list(blood_cell_count_genes$target.id, leukocyte_count_genes$target.id, 
                  myeloid_white_cell_count_genes$target.id, monocyte_count_genes$target.id, 
                  cancer_genes$target.id, parkinson_genes$target.id, alzheimer_genes$target.id)
names(OTAR_sets) <- paste(c("Blood cell count", "Leukocyte count", "Myeloid white cell count", "Monocyte count", 
                            "Cancer", "Parkinson's disease", "Alzheimer's disease"), 
                          paste0("(n = ", c(lapply(OTAR_sets, length) %>% unlist), ")"))

sapply(names(OTAR_sets), function(i){ 
  sapply(names(OTAR_sets), function(j){
    print( length(intersect(OTAR_sets[[i]], OTAR_sets[[j]])) ) } )
  })

#                                    Blood cell count (n = 9980)   Leukocyte count (n = 7818)  Myeloid white cell count (n = 5679)  Monocyte count (n = 3240)
# Blood cell count (n = 9980)                               9980                         7818                                 5679                       3240
# Leukocyte count (n = 7818)                                7818                         7818                                 5679                       3240
# Myeloid white cell count (n = 5679)                       5679                         5679                                 5679                       2355
# Monocyte count (n = 3240)                                 3240                         3240                                 2355                       3240
# Cancer (n = 20562)                                        9339                         7359                                 5367                       3065
# Parkinson's disease (n = 4455)                            2561                         2020                                 1515                        922
# Alzheimer's disease (n = 4626)                            2524                         1987                                 1461                        881
#                                     Cancer (n = 20562)   Parkinson's disease (n = 4455)   Alzheimer's disease (n = 4626)
# Blood cell count (n = 9980)                       9339                             2561                             2524
# Leukocyte count (n = 7818)                        7359                             2020                             1987
# Myeloid white cell count (n = 5679)               5367                             1515                             1461
# Monocyte count (n = 3240)                         3065                              922                              881
# Cancer (n = 20562)                               20562                             4372                             4413
# Parkinson's disease (n = 4455)                    4372                             4455                             2119
# Alzheimer's disease (n = 4626)                    4413                             2119                             4626


## Assess enrichment of trait 1 genes among trait 2
overlaps_between_OTAR_traits <- function(OT_trait1, OT_trait2){
  
  OTAR_genes_trait1 <- get(paste0(OT_trait1, "_genes"))
  OTAR_genes_trait2 <- get(paste0(OT_trait2, "_genes"))
  
  OT_trait_labs <- c("blood_cell_count" = "complete blood cell count",
                     "leukocyte_count" = "leukocyte count",
                     "myeloid_white_cell_count" = "myeloid white cell count",
                     "monocyte_count" = "monocyte count",
                     "cancer" = "cancer",
                     "alzheimer" = "Alzheimer's disease",
                     "parkinson" = "Parkinson's disease")
  
  ## Subset trait 1 genes to same X num of closest genes x phenotype
  closest_genes_num <- sapply(paste0("closest_genes_", c("premac_vs_iPSC", "old_vs_young_premac", "microglia_vs_premac_IFN", 
                                                             "microglia_vs_premac_LPS", "microglia_vs_premac_untreated")), function(l){length(get(l))})
  overlaps <- vector()
  for(i in 1:length(closest_genes_num)){
    
    ## Top X targets for trait 1
    OTAR_genes_trait1_top <- OTAR_genes_trait1 %>% arrange(-score) %>% .[1:closest_genes_num[i], "target.id"]
      
      ## Enrichments among trait 2 targets above threshold score 
      thresholds <- seq(from = min(OTAR_genes_trait2$score), to = max(OTAR_genes_trait2$score), length.out = 20)
      
      for(threshold in thresholds){
        OTAR_genes_trait2_thresholded <- OTAR_genes_trait2[OTAR_genes_trait2$score >= threshold, "target.id"]
        
        ## Fisher test taking ~40k genes in the genome as universe
        ## (all existing genes with/without OTAR evidence)
        n_trait1 <- length(OTAR_genes_trait1_top)
        n_trait2 <- length(OTAR_genes_trait2_thresholded)
        
        n_trait2_and_1 =  length(intersect(OTAR_genes_trait1_top, OTAR_genes_trait2_thresholded))
        n_trait2_and_not1 = n_trait2 - n_trait2_and_1
        n_trait1_and_not2 = n_trait1 - n_trait2_and_1
        rest = 40000 - (n_trait2_and_1 + n_trait2_and_not1 + n_trait1_and_not2)
        
        m <- matrix(c(n_trait2_and_1, n_trait2_and_not1, n_trait1_and_not2, rest), ncol = 2)
        f <- fisher.test(m, alternative = "greater")
        p <- f$p.value
        
        overlaps <- rbind(overlaps, 
                          c("threshold" = threshold, 
                            "num_genes_trait1" = n_trait1, 
                            "num_genes_trait2" = n_trait2, 
                            "overlap_size" = n_trait2_and_1, 
                            "p" = p))
      
        }
      
  }
  
  overlaps <- as.data.frame(overlaps)
  overlaps[,1:5] <- apply(overlaps[,1:5], 2, as.numeric)
  
  plot <- ggplot(overlaps, aes(x = threshold, y = overlap_size)) +
    facet_grid(rows = vars(num_genes_trait1)) +
    geom_bar(stat = "identity", fill = "mistyrose1", color = "black", linewidth = 0.3, width = 0.03) +
    theme_classic() +
    labs(title = paste("Enrichment of",  OT_trait_labs[OT_trait1], "genes among", OT_trait_labs[OT_trait2], "genes"),
         x = paste("Association threshold for",  OT_trait_labs[OT_trait2], "genes"), 
         y = "Size of overlap") +
    scale_y_continuous(expand = c(0,0), limits = c(0, max(overlaps$overlap_size) + 15)) +
    geom_text(data = subset(overlaps, p < 0.05),
              aes(x = threshold, y = overlap_size + (0.45 * (max(overlaps$overlap_size) + 1)/10), label = "*"), 
              size = 3, fontface = "bold", color = "red") +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5), 
          axis.title = element_text(size = 9),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 8))

  return(plot)
}

p1 <- overlaps_between_OTAR_traits("leukocyte_count", "blood_cell_count")
p2 <- overlaps_between_OTAR_traits("myeloid_white_cell_count", "blood_cell_count")
p3 <- overlaps_between_OTAR_traits("monocyte_count", "blood_cell_count")
p4 <- overlaps_between_OTAR_traits("myeloid_white_cell_count", "leukocyte_count")
p5 <- overlaps_between_OTAR_traits("monocyte_count", "leukocyte_count")
p6 <- overlaps_between_OTAR_traits("monocyte_count", "myeloid_white_cell_count")

plot_grid(plotlist = list(p1, p2, p3, p4, p5, p6), nrow = 2)
ggsave(paste0(plotdir, "/Overlaps_between_OTAR_cell_count_trait_targets.pdf"), height = 19, width = 20)



## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#     8.3.1 Overlaps between closest genes and proliferation-related genes
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

## Function to get overlaps with OT targets at multiple association cutoffs
OT_overlaps_at_cutoff <- function(prolif_pheno, OT_trait){
  
  closest_genes <- get(paste0("closest_genes_", prolif_pheno)) 
  OT_genes_data <- get(paste0(OT_trait, "_genes"))
  
  prolif_pheno_labs <- switch (prolif_pheno,
                               "premac_vs_iPSC" = c("iPSC", "mac. precursor"),
                               "old_vs_young_premac" = c("young", "old mac. precursor"),
                               "microglia_vs_premac_IFN" = c("mac. precursor", "microglia in IFN"),
                               "microglia_vs_premac_LPS" = c("mac. precursor", "microglia in LPS"),
                               "microglia_vs_premac_untreated" = c("mac. precursor", "microglia in untreated")
                               )
  OT_trait_lab <- switch(OT_trait,
                         "blood_cell_count" = "complete blood cell count",
                         "leukocyte_count" = "leukocyte count",
                         "myeloid_white_cell_count" = "myeloid white cell count",
                         "monocyte_count" = "monocyte count",
                         "cancer" = "cancer",
                         "alzheimer" = "Alzheimer's disease",
                         "parkinson" = "Parkinson's disease"
                         )
  
  thresholds <- seq(from = min(OT_genes_data$score), 
                    to = max(OT_genes_data$score), 
                    length.out = 20)
  
  overlaps <- vector()
  for(threshold in thresholds){
     
    OT_genes <- OT_genes_data[OT_genes_data$score >= threshold, "target.id"]
    
    overlapping_genes <- intersect(closest_genes, OT_genes)
    overlapping_genes_names <- OT_genes_data[match(overlapping_genes, OT_genes_data$target.id), ] %>% 
         arrange(-score) %>% .[, "target.approvedSymbol"]
    
    ## Fisher test taking ~40k (protein coding + lncRNAs + miRNAs) genes in the genome as universe
    ## (these are the genes for which we exmained proximity to lead SNPs and with/without OTAR evidence)
    n_closest <- length(closest_genes)
    n_OT_genes <- length(OT_genes)
    
    n_closest_OT =  length(overlapping_genes)
    n_closest_noOT = n_closest - n_closest_OT
    n_non_closest_OT = n_OT_genes - n_closest_OT
    n_non_closest_noOT = 40000 - c(n_closest_OT + n_closest_noOT + n_non_closest_OT)
    
    m <- matrix(c(n_closest_OT, n_closest_noOT, n_non_closest_OT, n_non_closest_noOT), ncol = 2)
    f <- fisher.test(m, alternative = "greater")
    p <- f$p.value
    
    overlaps <- rbind(overlaps, 
                      c("threshold" = threshold, 
                        "num_OT_genes" = n_OT_genes, 
                        "overlap_size" = n_closest_OT, 
                        "p" = p, 
                        "overlapping_genes" = paste(overlapping_genes_names[1:min(n_closest_OT, 28)], collapse = "\n")))
    
  }
  
  overlaps <- as.data.frame(overlaps)
  overlaps[,1:4] <- apply(overlaps[,1:4], 2, as.numeric)
  
  plot <- ggplot(overlaps, aes(x = threshold, y = overlap_size)) +
    geom_bar(stat = "identity", fill = "mistyrose1", color = "black", linewidth = 0.3, width = 0.03) +
    theme_classic() +
    labs(title = paste("Overlap between closest genes to lead SNPs for prolif. from", "\n",
                       prolif_pheno_labs[1], "to", prolif_pheno_labs[2], 
                       paste0("(n = ", n_closest, ")"),
                       "and target genes for", OT_trait_lab),
         x = paste("Gene-Trait association threshold", "\n", "(n of associated genes above threshold)"), 
         y = "Size of overlap with GWAS closest genes") +
    scale_x_continuous(expand = c(0.012, 0), breaks = signif(overlaps$threshold, digits = 2), 
                       labels = paste(signif(overlaps$threshold, digits = 2), "\n", 
                                      paste0("(n = ", overlaps$num_OT_genes, ")"))) +
    scale_y_continuous(expand = c(0,0), limits = c(0, max(overlaps$overlap_size)*(1 + (1/10)) )) +
    geom_text(aes(x = threshold, y = overlap_size + (0.25 * (max(overlaps$overlap_size) + 1)/10), 
                  label = overlap_size), size = 2.5) +
    geom_text(data = subset(overlaps, p < 0.05),
              aes(x = threshold, y = overlap_size + (0.45 * (max(overlaps$overlap_size) + 1)/10), label = "*"), 
              size = 3, fontface = "bold", color = "red") +
    geom_text(aes(x = 0.65, y = (1 - 0.4)*(max(overlaps$overlap_size)), 
                  label = paste("All overlapping genes:", "\n", 
                              overlaps[1,"overlapping_genes"])), size = 2.5, fontface = "italic") +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5), 
          axis.title = element_text(size = 9),
          axis.text.x = element_text(size = 7, angle = 55, hjust = 0.8),
          axis.text.y = element_text(size = 8))
  
  return(plot)
}


#~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ 
##   Potential target genes of (iPSC -> young preMac) proliferation lead SNPs
#~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ 

## Overlaps with Burden signif genes
length(intersect(closest_genes_premac_vs_iPSC, deleterious_Burden_signif_premac_vs_iPSC_ids)) # 0
length(intersect(closest_genes_premac_vs_iPSC, deleterious_Burden_signif_old_vs_young_premac_ids)) # 0
length(intersect(closest_genes_premac_vs_iPSC, deleterious_Burden_signif_microglia_vs_premac_ids)) # 0

## Overlaps with DEGs
length(intersect(closest_genes_premac_vs_iPSC, de_genes_proliferation_IFN)) / length(closest_genes_premac_vs_iPSC) * 100
# 23.94% (17/71 - 6579 DEGs)
length(intersect(closest_genes_premac_vs_iPSC, de_genes_proliferation_LPS)) / length(closest_genes_premac_vs_iPSC) * 100
# 12.67% (9/71 - 3035 DEGs)
length(intersect(closest_genes_premac_vs_iPSC, de_genes_proliferation_untreated)) / length(closest_genes_premac_vs_iPSC) * 100
# 21.12% (15/71 - 7454 DEGs)

## Overlaps with OT genes at variable association thresholds
prolif_pheno = "premac_vs_iPSC"
p1 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "blood_cell_count")
p2 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "leukocyte_count")
p3 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "myeloid_white_cell_count")
p4 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "monocyte_count")
p5 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "cancer")
p6 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "parkinson")
p7 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "alzheimer")

plot_grid(p1, p2, p3, p4, p5, p6, p7, nrow = 2, align = "vh")
ggsave(filename = paste0(plotdir, "/Overlaps_closest_genes_for_", prolif_pheno, "_and_OT_genes.pdf"), width = 27, height = 12)


#~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ 
##   Potential target genes of (yound -> old preMac) proliferation lead SNPs
#~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ 

intersect(closest_genes_old_vs_young_premac, deleterious_Burden_signif_premac_vs_iPSC_ids) # 0
intersect(closest_genes_old_vs_young_premac, deleterious_Burden_signif_old_vs_young_premac_ids) # 0
intersect(closest_genes_old_vs_young_premac, deleterious_Burden_signif_microglia_vs_premac_ids) # 0

length(intersect(closest_genes_old_vs_young_premac, de_genes_proliferation_IFN)) / length(closest_genes_old_vs_young_premac) * 100
# 10.23% (9/88 - 6579 DEGs)
length(intersect(closest_genes_old_vs_young_premac, de_genes_proliferation_LPS)) / length(closest_genes_old_vs_young_premac) * 100
# 6.81% (6/88 - 3035 DEGs)
length(intersect(closest_genes_old_vs_young_premac, de_genes_proliferation_untreated)) / length(closest_genes_old_vs_young_premac) * 100
# 12.5% (11/88 - 7454 DEGs)

## Overlaps with OT genes 
prolif_pheno = "old_vs_young_premac"
p1 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "blood_cell_count")
p2 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "leukocyte_count")
p3 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "myeloid_white_cell_count")
p4 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "monocyte_count")
p5 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "cancer")
p6 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "parkinson")
p7 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "alzheimer")

plot_grid(p1, p2, p3, p4, p5, p6, p7, nrow = 2, align = "vh")
ggsave(filename = paste0(plotdir, "/Overlaps_closest_genes_for_", prolif_pheno, "_and_OT_genes.pdf"), width = 27, height = 12)


#~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ 
##   Potential target genes of (preMac -> microglia in IFN) proliferation lead SNPs
#~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ 

intersect(closest_genes_microglia_vs_premac_IFN, deleterious_Burden_signif_premac_vs_iPSC_ids) # 0
intersect(closest_genes_microglia_vs_premac_IFN, deleterious_Burden_signif_old_vs_young_premac_ids) # 0
intersect(closest_genes_microglia_vs_premac_IFN, deleterious_Burden_signif_microglia_vs_premac_ids) 
# "ENSG00000182223" - ZAR1 (1/106 - 18 Burden signif. genes)

length(intersect(closest_genes_microglia_vs_premac_IFN, de_genes_proliferation_IFN)) / length(closest_genes_microglia_vs_premac_IFN) * 100
# 13.207% (14/106 - 6579 DEGs)
length(intersect(closest_genes_microglia_vs_premac_IFN, de_genes_proliferation_LPS)) / length(closest_genes_microglia_vs_premac_IFN) * 100
# 4.71% (5/106 - 3035 DEGs)
length(intersect(closest_genes_microglia_vs_premac_IFN, de_genes_proliferation_untreated)) / length(closest_genes_microglia_vs_premac_IFN) * 100
# 15.09% (16/106 - 7454 DEGs)

## Overlaps with OT genes 
prolif_pheno = "microglia_vs_premac_IFN"
p1 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "blood_cell_count")
p2 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "leukocyte_count")
p3 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "myeloid_white_cell_count")
p4 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "monocyte_count")
p5 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "cancer")
p6 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "parkinson")
p7 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "alzheimer")

plot_grid(p1, p2, p3, p4, p5, p6, p7, nrow = 2, align = "vh")
ggsave(filename = paste0(plotdir, "/Overlaps_closest_genes_for_", prolif_pheno, "_and_OT_genes.pdf"), width = 27, height = 12)


#~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ 
##   Potential target genes of (preMac -> microglia in LPS) proliferation lead SNPs
#~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ 

intersect(closest_genes_microglia_vs_premac_LPS, deleterious_Burden_signif_premac_vs_iPSC_ids) # 0
intersect(closest_genes_microglia_vs_premac_LPS, deleterious_Burden_signif_old_vs_young_premac_ids) # 0
intersect(closest_genes_microglia_vs_premac_LPS, deleterious_Burden_signif_microglia_vs_premac_ids) 
# "ENSG00000182223" - ZAR1 (1/123 - 18 Burden signif. genes)

length(intersect(closest_genes_microglia_vs_premac_LPS, de_genes_proliferation_IFN)) / length(closest_genes_microglia_vs_premac_LPS) * 100
# 21.13% (26/123 - 6579 DEGs)
length(intersect(closest_genes_microglia_vs_premac_LPS, de_genes_proliferation_LPS)) / length(closest_genes_microglia_vs_premac_LPS) * 100
# 8.13% (10/123 - 3035 DEGs)
length(intersect(closest_genes_microglia_vs_premac_LPS, de_genes_proliferation_untreated)) / length(closest_genes_microglia_vs_premac_LPS) * 100
# 19.51% (24/123 - 7454 DEGs)

## Overlaps with OT genes 
prolif_pheno = "microglia_vs_premac_LPS"
p1 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "blood_cell_count")
p2 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "leukocyte_count")
p3 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "myeloid_white_cell_count")
p4 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "monocyte_count")
p5 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "cancer")
p6 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "parkinson")
p7 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "alzheimer")

plot_grid(p1, p2, p3, p4, p5, p6, p7, nrow = 2, align = "vh")
ggsave(filename = paste0(plotdir, "/Overlaps_closest_genes_for_", prolif_pheno, "_and_OT_genes.pdf"), width = 27, height = 12)


#~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ 
##   Potential target genes of (preMac -> microglia in untreated) proliferation lead SNPs
#~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ 

intersect(closest_genes_microglia_vs_premac_untreated, deleterious_Burden_signif_premac_vs_iPSC_ids) # 0
intersect(closest_genes_microglia_vs_premac_untreated, deleterious_Burden_signif_old_vs_young_premac_ids) # 0
intersect(closest_genes_microglia_vs_premac_untreated, deleterious_Burden_signif_microglia_vs_premac_ids)
# "ENSG00000168918" - INPP5D (1/51 - 18 Burden signif. genes)

length(intersect(closest_genes_microglia_vs_premac_untreated, de_genes_proliferation_IFN)) / length(closest_genes_microglia_vs_premac_untreated) * 100
# 7.84% (4/51 - 6579 DEGs)
length(intersect(closest_genes_microglia_vs_premac_untreated, de_genes_proliferation_LPS)) / length(closest_genes_microglia_vs_premac_untreated) * 100
# 3.92% (2/51 - 3035 DEGs)
length(intersect(closest_genes_microglia_vs_premac_untreated, de_genes_proliferation_untreated)) / length(closest_genes_microglia_vs_premac_untreated) * 100
# 11.76% (6/51 - 7454 DEGs)

## Overlaps with OT genes 
prolif_pheno = "microglia_vs_premac_untreated"
p1 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "blood_cell_count")
p2 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "leukocyte_count")
p3 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "myeloid_white_cell_count")
p4 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "monocyte_count")
p5 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "cancer")
p6 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "parkinson")
p7 <- OT_overlaps_at_cutoff(prolif_pheno, OT_trait = "alzheimer")

plot_grid(p1, p2, p3, p4, p5, p6, p7, nrow = 2, align = "vh")
ggsave(filename = paste0(plotdir, "/Overlaps_closest_genes_for_", prolif_pheno, "_and_OT_genes.pdf"), width = 27, height = 12)



## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#            8.3.2 Further explorations for overlapping genes 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

## Plot allele dosage vs prolif
allele_dosage_vs_prolif <- function(variant_rsID, chunk, gene_id, phenotype){
  
  ## Association results for variant 
  load(paste0(input_dir02, "/summary_stats_", phenotype, ".Rdata"), verbose = T)
  variant_data <- subset(summary_stats, rsID == variant_rsID)
  
  variant_id <- variant_data$variant_id
  effect_size <- signif(variant_data$Estimate, digits = 3)
  t <- signif(variant_data$`t value`, digits = 3)
  pval <- signif(variant_data$`Pr(>|t|)`, digits = 3)
  
  ## Allele dosage and prolif data
  line_prop_changes <-  load(paste0(input_dir00, "/", phenotype, "_data/", 
                          "line_prop_changes_", phenotype, "_genotype_", chunk, ".Rdata"), verbose = T)
  line_prop_changes <- eval(parse_expr(line_prop_changes))
  
  allele_dosage_prolif <- line_prop_changes[,  c("line", "mean_scaled_log_fraction", "npool", "sex", variant_id)] %>% rename(genotype = variant_id)
  allele_dosage_prolif$genotype_char <- as.character(allele_dosage_prolif$genotype)
  
  p <- ggplot(allele_dosage_prolif, aes(x = genotype_char, y = mean_scaled_log_fraction)) +
    geom_point(color = "black") +
    geom_violin(aes(x = genotype_char, y = mean_scaled_log_fraction), 
                color = "black", fill = NA, size = 0.25,  width=0.32) +
    stat_smooth(aes(group = 1), geom = "line", alpha = 1, size = 1, method = lm, color='firebrick2') +
    theme_bw() +
    labs(title = paste0(variant_rsID, " (", gsub("_", "-", variant_id), ") close to ", gene_id), x = "Minor allele dosage", 
         y = paste("Proliferation mean scaled log-fraction", "\n", 
                   paste0("(", strsplit(phenotype, "_")[[1]][1], " in ", 
                          strsplit(phenotype, "_")[[1]][4], " vs ",  strsplit(phenotype, "_")[[1]][3], ")") )) +
    geom_text(aes(x = 1.5, y = 1, label = paste("t =", t, "\n", "p =", pval)), size = 3) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 9),
          axis.text = element_text(size = 8),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 9))
  
  return(p)
}


## Plot gene del burden vs prolif

## Num of del variants x gene x donor 
mutBurden_del <- read_rds("/lustre/scratch123/hgi/projects/otar2065/resources/Puigdevall_Neuroseq_efficiency_2023/mutBurdenTabs/mutBurden_del.RDS") %>% as.data.frame()
mutBurden_del = mutBurden_del %>% t() %>% as.data.frame()
mutBurden_del$line <-  gsub(".*-", "", rownames(mutBurden_del))

del_burden_vs_prolif <- function(gene_id, phenotype){
  
  ## Prolif data
  line_prop_changes <-  load(paste0(input_dir00, "/", phenotype, "_data/", 
                                    "line_prop_changes_", phenotype, "_genotype_", 187, ".Rdata"), verbose = T)
  line_prop_changes <- eval(parse_expr(line_prop_changes))[, 1:4]
  
  ## Add del burdens for gene
  line_prop_changes <- left_join(line_prop_changes, mutBurden_del[, c("line", gene_id)], by = "line")
  
  ## Gene stats
  del_burden_results <- as.data.frame(read_csv("/lustre/scratch123/hgi/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/data/results/2.2.rare_vars_vs_prolif/del_burden_scaled_centered_prop_pvals.csv"))
  del_burden_results <- subset(del_burden_results, gene_name == gene_id & comparison == unique(del_burden_results$comparison) %>% .[grep(strsplit(phenotype, "_")[[1]][1], .)])
  coef = signif(del_burden_results$coef, digits = 3)
  padj = signif(del_burden_results$p_Bonf, digits = 3)
  
  ## Remove lines with NA burdens
  line_prop_changes <- line_prop_changes[-which(is.na(line_prop_changes[, gene_id])), ]
  
  p <- ggplot(line_prop_changes, aes(x = as.character(get(gene_id)), y = mean_scaled_log_fraction)) +
    geom_point(color = "black") +
    geom_violin(aes(x = as.character(get(gene_id)), y = mean_scaled_log_fraction), color = "black", 
                fill = NA, size = 0.25,  width=0.32) +
    stat_smooth(aes(group = 1), geom = "line", alpha = 1, size = 1, method = lm, color='firebrick2') +
    theme_bw() +
    labs(title = paste0(gene_id), x = "Burden of deleterious variants", 
         y = paste("Proliferation mean scaled log-fraction", "\n", 
                   paste0("(", strsplit(phenotype, "_")[[1]][1], " vs ",  
                          strsplit(phenotype, "_")[[1]][3], ")") )) +
    geom_text(aes(x = 3, y = -1, label = paste("Effect size =", coef, "\n", "adj. p =", padj)), size = 3) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 9),
          axis.text = element_text(size = 8),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 9))
  
  return(p)
}


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
##                           ZAR1 - ENSG00000182223
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

## Lead SNP in region overlapping ZAR1 for preMac -> microglia in IFN
gene_tracks_microglia_vs_premac_IFN <-  lapply(gene_tracks_across_chrs_microglia_vs_premac_IFN, function(chr) 
                                            { if(length(chr)>0){do.call(rbind, chr) %>% rownames_to_column("leadSNP")} }) %>% do.call(rbind, .)

subset(gene_tracks_microglia_vs_premac_IFN, gene_id == "ENSG00000182223")
#       leadSNP seqnames    start      end  width  strand          gene_id  gene_name    gene_biotype  seq_coord_system
#  rs17656010.6        4 48490252 48494389   4138       +  ENSG00000182223       ZAR1  protein_coding        chromosome
#                                         description    gene_id_version   canonical_transcript   symbol entrezid tx_seq_start tx_seq_end
#  zygote arrest 1 [Source:HGNC Symbol;Acc:HGNC:20436]  ENSG00000182223.7       ENST00000327939     ZAR1   326340     48490252   48494389

p1 <- allele_dosage_vs_prolif("rs17656010", chunk = 1627, gene_id = "ZAR1", phenotype = "microglia_vs_premac_IFN")
p2 <- del_burden_vs_prolif("ZAR1", phenotype = "microglia_vs_premac_IFN")
plot_grid(p1, p2, ncol = 2, align = "h")
ggsave(filename = paste0(plotdir, "/ZAR1_rs17656010_effects_on_prolif_microglia_vs_premac_IFN.pdf"), width = 10, height = 5)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## Lead SNP in region overlapping ZAR1 for preMac -> microglia in LPS
gene_tracks_microglia_vs_premac_LPS <-  lapply(gene_tracks_across_chrs_microglia_vs_premac_LPS, function(chr) 
{ if(length(chr)>0){do.call(rbind, chr) %>% rownames_to_column("leadSNP")} }) %>% do.call(rbind, .)

subset(gene_tracks_microglia_vs_premac_LPS, gene_id == "ENSG00000182223")
#      leadSNP  seqnames     start      end  width strand          gene_id  gene_name   gene_biotype seq_coord_system
#  rs6447641.7         4  48490252 48494389   4138      +  ENSG00000182223       ZAR1 protein_coding       chromosome
#                                          description    gene_id_version  canonical_transcript symbol entrezid tx_seq_start tx_seq_end
#  zygote arrest 1 [Source:HGNC Symbol;Acc:HGNC:20436]  ENSG00000182223.7       ENST00000327939   ZAR1   326340     48490252   48494389

p3 <- allele_dosage_vs_prolif("rs6447641", chunk = 1626, gene_id = "ZAR1", phenotype = "microglia_vs_premac_LPS")
p4 <- del_burden_vs_prolif("ZAR1", phenotype = "microglia_vs_premac_LPS")
plot_grid(p3, p4, ncol = 2, align = "h")
ggsave(filename = paste0(plotdir, "/ZAR1_rs6447641_effects_on_prolif_microglia_vs_premac_LPS.pdf"), width = 10, height = 5)


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
##                          INPP5D - ENSG00000168918
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

## Lead SNP in region overlapping INPP5D for preMac -> microglia in untreated
gene_tracks_microglia_vs_premac_untreated <-  lapply(gene_tracks_across_chrs_microglia_vs_premac_untreated, function(chr) 
{ if(length(chr)>0){do.call(rbind, chr) %>% rownames_to_column("leadSNP")} }) %>% do.call(rbind, .)

subset(gene_tracks_microglia_vs_premac_untreated, gene_id == "ENSG00000168918")
#      leadSNP  seqnames      start       end   width  strand          gene_id  gene_name   gene_biotype seq_coord_system
#  rs4663490.1         2  233059967 233207903  147937       +  ENSG00000168918     INPP5D protein_coding       chromosome
#                                                               description     gene_id_version   canonical_transcript  symbol entrezid
#  inositol polyphosphate-5-phosphatase D [Source:HGNC Symbol;Acc:HGNC:6079]  ENSG00000168918.14       ENST00000445964  INPP5D     3635
#  tx_seq_start tx_seq_end
#    233059967  233207903

p5 <- allele_dosage_vs_prolif("rs4663490", chunk = 1009, gene_id = "INPP5D", phenotype = "microglia_vs_premac_untreated")
p6 <- del_burden_vs_prolif("INPP5D", phenotype = "microglia_vs_premac_untreated")
plot_grid(p5, p6, ncol = 2, align = "h")
ggsave(filename = paste0(plotdir, "/INPP5D_rs4663490_effects_on_prolif_microglia_vs_premac_untreated.pdf"), width = 10, height = 5)







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
# date     2025-02-25
# rstudio  2024.04.0+735 Chocolate Cosmos (server)
# pandoc   3.1.12.3 @ /opt/view/bin/ (via rmarkdown)
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package          * version     date (UTC) lib source
# AnnotationDbi      1.64.1      2023-11-03 [1] Bioconductor
# Biobase            2.62.0      2023-10-24 [1] Bioconductor
# BiocFileCache      2.10.2      2024-03-27 [1] Bioconductor 3.18 (R 4.3.1)
# BiocGenerics       0.48.1      2023-11-01 [1] Bioconductor
# biomaRt          * 2.58.2      2024-01-30 [1] Bioconductor 3.18 (R 4.3.1)
# Biostrings         2.70.3      2024-03-13 [1] Bioconductor 3.18 (R 4.3.1)
# bit                4.0.5       2022-11-15 [1] CRAN (R 4.3.1)
# bit64              4.0.5       2020-08-30 [1] CRAN (R 4.3.1)
# bitops             1.0-7       2021-04-24 [1] CRAN (R 4.3.1)
# blob               1.2.4       2023-03-17 [1] CRAN (R 4.3.1)
# cachem             1.0.8       2023-05-01 [1] CRAN (R 4.3.1)
# cli                3.6.2       2023-12-11 [1] CRAN (R 4.3.1)
# colorspace         2.1-0       2023-01-23 [1] CRAN (R 4.3.1)
# cowplot          * 1.1.3       2024-01-22 [1] CRAN (R 4.3.1)
# crayon             1.5.2       2022-09-29 [1] CRAN (R 4.3.1)
# curl               5.2.1       2024-03-01 [1] CRAN (R 4.3.1)
# DBI                1.2.2       2024-02-16 [1] CRAN (R 4.3.1)
# dbplyr             2.5.0       2024-03-19 [1] CRAN (R 4.3.1)
# digest             0.6.35      2024-03-11 [1] CRAN (R 4.3.1)
# dplyr            * 1.1.4       2023-11-17 [1] CRAN (R 4.3.1)
# evaluate           0.23        2023-11-01 [1] CRAN (R 4.3.1)
# fansi              1.0.6       2023-12-08 [1] CRAN (R 4.3.1)
# farver             2.1.1       2022-07-06 [1] CRAN (R 4.3.1)
# fastmap            1.1.1       2023-02-24 [1] CRAN (R 4.3.1)
# filelock           1.0.3       2023-12-11 [1] CRAN (R 4.3.1)
# forcats          * 1.0.0       2023-01-29 [1] CRAN (R 4.3.1)
# generics           0.1.3       2022-07-05 [1] CRAN (R 4.3.1)
# GenomeInfoDb       1.38.8      2024-03-15 [1] Bioconductor 3.18 (R 4.3.1)
# GenomeInfoDbData   1.2.11      2025-02-07 [1] Bioconductor
# ggplot2          * 3.5.1       2024-04-23 [1] CRAN (R 4.3.1)
# ggrepel            0.9.5       2024-01-10 [1] CRAN (R 4.3.1)
# glue               1.7.0       2024-01-09 [1] CRAN (R 4.3.1)
# gtable             0.3.4       2023-08-21 [1] CRAN (R 4.3.1)
# hms                1.1.3       2023-03-21 [1] CRAN (R 4.3.1)
# htmltools          0.5.8       2024-03-25 [1] CRAN (R 4.3.1)
# httr             * 1.4.7       2023-08-15 [1] CRAN (R 4.3.1)
# IRanges            2.36.0      2023-10-24 [1] Bioconductor
# jsonlite           1.8.8       2023-12-04 [1] CRAN (R 4.3.1)
# jtools             2.2.2       2023-07-11 [1] CRAN (R 4.3.1)
# KEGGREST           1.42.0      2023-10-24 [1] Bioconductor
# knitr              1.45        2023-10-30 [1] CRAN (R 4.3.1)
# labeling           0.4.3       2023-08-29 [1] CRAN (R 4.3.1)
# lifecycle          1.0.4       2023-11-07 [1] CRAN (R 4.3.1)
# lubridate        * 1.9.3       2023-09-27 [1] CRAN (R 4.3.1)
# magrittr           2.0.3       2022-03-30 [1] CRAN (R 4.3.1)
# memoise            2.0.1       2021-11-26 [1] CRAN (R 4.3.1)
# munsell            0.5.1       2024-04-01 [1] CRAN (R 4.3.1)
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
# Rcpp               1.0.12      2024-01-09 [1] CRAN (R 4.3.1)
# RCurl              1.98-1.14   2024-01-09 [1] CRAN (R 4.3.1)
# readr            * 2.1.5       2024-01-10 [1] CRAN (R 4.3.1)
# rlang              1.1.3       2024-01-10 [1] CRAN (R 4.3.1)
# rmarkdown          2.26        2024-03-05 [1] CRAN (R 4.3.1)
# RSQLite            2.3.6       2024-03-31 [1] CRAN (R 4.3.1)
# rstudioapi         0.16.0      2024-03-24 [1] CRAN (R 4.3.1)
# S4Vectors          0.40.2      2023-11-23 [1] Bioconductor 3.18 (R 4.3.1)
# scales             1.3.0       2023-11-28 [1] CRAN (R 4.3.1)
# sessioninfo      * 1.2.2       2021-12-06 [1] CRAN (R 4.3.1)
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
# xfun               0.43        2024-03-25 [1] CRAN (R 4.3.1)
# XML                3.99-0.16.1 2024-01-22 [1] CRAN (R 4.3.1)
# xml2               1.3.6       2023-12-04 [1] CRAN (R 4.3.1)
# XVector            0.42.0      2023-10-24 [1] Bioconductor
# yaml               2.3.8       2023-12-11 [1] CRAN (R 4.3.1)
# zlibbioc           1.48.2      2024-03-13 [1] Bioconductor 3.18 (R 4.3.1)
# 
# [1] /opt/view/rlib/R/library
# [2] /software/hgi/installs/softpack/rstudio/.spack/linux-ubuntu22.04-x86_64_v3/gcc-11.4.0/r-4.3.1-bfwldrk76z6f52upk47zepliekn7ayqz/rlib/R/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────


