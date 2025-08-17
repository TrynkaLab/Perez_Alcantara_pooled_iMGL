
library(tidyverse)
library(dplyr)
library(biomaRt)
library(httr)
library(sessioninfo)


######################################################################################################
##  11. Extract Open Targets GWAS and additional genetic studies to support coloc genes associations
######################################################################################################

#-------------------------------------------------------------------------------
#     11.1 Extract GWAS publications used in Open Targets association scores
#-------------------------------------------------------------------------------
#  Code to extract the GWAS publications from which gene-disease genetics-based 
#  association scores were calculated in Open Targets. 
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

## Define dirs
setwd("/lustre/scratch125/humgen/teams_v2/trynka/otar2065/OTAR2065_Daianna/")
input_dir = paste("input_data", "11_genetics_evidence_for_coloc_genes", "01_extract_OTAR_GWAS_publications", sep = "/")
outdir = paste("output_data", "11_genetics_evidence_for_coloc_genes", "01_extract_OTAR_GWAS_publications", sep = "/")
dir.create(outdir, recursive = T)
dir.create(input_dir, recursive = T)

## Input 79 coloc genes (77 unique)
coloc_genes <- readxl::read_xls(paste0(input_dir, "/Suppl_table_10_OTAR_annotated_signif_colocs.xls"))

## Add Ensembl ID to gene symbols
ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
ensembl_ids <- getBM(values = unique(coloc_genes$eGene),
                     mart = ensembl,
                     attributes = c("external_gene_name",
                                    "ensembl_gene_id"),
                     filters = "external_gene_name")

## Keep 1st occurrence only
ensembl_ids <- ensembl_ids[-which(duplicated(ensembl_ids$external_gene_name)), ]
colnames(ensembl_ids) <- c("eGene", "ensembl_gene_id")

## Add Ensembl ID of NSF used in OTAR
ensembl_ids[which(ensembl_ids$eGene == "NSF"), "ensembl_gene_id"] <- "ENSG00000073969"

coloc_genes <- left_join(coloc_genes, ensembl_ids, by = "eGene")

## Add disease EFO ID
coloc_genes <- coloc_genes %>% 
  mutate(efo_id = case_when(
    GWAS == "AD" ~ "MONDO_0004975",
    GWAS == "PD" ~ "MONDO_0021095",
    GWAS == "MS" ~ "MONDO_0005301",
    GWAS == "ALS" ~ "MONDO_0004976")
  )


# Build query string to extract GWAS publications for x gene x disease 
get_gene_disease_GWAS <- function(gene_id, disease_id, n_targets = 1000){
  
  query_string = paste0("
    query GwasQuery($ensemblId: String!, $efoId: String!) {
        disease(efoId: $efoId) {
          gwasCredibleSets: evidences(
            ensemblIds: [$ensemblId]
            enableIndirect: true
            datasourceIds: [\"gwas_credible_sets\"]
            size:", n_targets, 
                        ") {
            count
            rows {
              credibleSet {
                study {
                  id
                  publicationFirstAuthor
                  publicationDate
                  pubmedId
                }
              }
            }
          }
        }
    }"
  )
  
  
  # Set base URL of GraphQL API endpoint
  base_url <- "https://api.platform.opentargets.org/api/v4/graphql"
  
  # Set variables object with arguments to be passed to endpoint
  variables <- list("efoId" = disease_id, "ensemblId" = gene_id)
  
  # Construct POST request body object with query string and variables
  post_body <- list(query = query_string, variables = variables)
  
  # POST request
  r <- POST(url = base_url, body = post_body, encode = 'json')
  
  return(r)
}


for(i in 1:dim(coloc_genes)[1]){
  
  gene_id <- coloc_genes[i, "ensembl_gene_id"] %>% as.character()
  disease_id <- coloc_genes[i, "efo_id"] %>% as.character()
  
  r <- get_gene_disease_GWAS(gene_id, disease_id)
  contents <- content(r)
  
  if(contents$data$disease$gwasCredibleSets$count > 0){
    
    ## FINNGEN studies have other format (add provided data)
    
    ## Extract 1st author surname 
    studies <- do.call(rbind,  lapply(contents$data$disease$gwasCredibleSets$rows, unlist)) %>% as.data.frame()
    surname <- gsub(" [A-Z]{1,2}$", "", studies$credibleSet.study.publicationFirstAuthor)
    if(length(surname)>0){studies$surname <- surname} else {studies$surname <- ""}
    
    ## Extract publication year
    year <- str_split_i(studies$credibleSet.study.publicationDate, "-", 1)
    if(length(year)>0){studies$year <- year} else {studies$year <- ""}
    
    ## Add PubMed ID
    pubmedId <- studies$credibleSet.study.pubmedId
    if(length(pubmedId)>0){studies$pubmedId <- pubmedId} else {studies$pubmedId <- ""}
    
    studies$info_to_add <- paste0(studies$surname, " ", studies$year, " (", studies$pubmedId, ")")
    studies[grep("FINNGEN", studies$credibleSet.study.id), "info_to_add"] <- studies[grep("FINNGEN", studies$credibleSet.study.id), "credibleSet.study.id"]
    
    coloc_genes[i, "GWAS_in_gwasCredibleSets_score"] <- paste(sort(unique(studies$info_to_add)), collapse = ', ')
    
  } else{
    coloc_genes[i, "GWAS_in_gwasCredibleSets_score"] <- "NA"
  }
  
}

## Export table                                                                                                                                                                                                                                                                                                        6 PILRA
write_excel_csv(coloc_genes, file = paste0(outdir, "/coloc_genes_with_OTRA_GWAS.csv"), col_names = T)



## Scores exported on May 27, 2025:
AD_scores <- read_tsv(file = paste0(input_dir, "/OT-MONDO_0004975-associated-targets-27_05_2025-v25_03-2.tsv")) %>% as.data.frame()
PD_scores <- read_tsv(file = paste0(input_dir, "/OT-MONDO_0021095-associated-targets-27_05_2025-v25_03.tsv")) %>% as.data.frame()
MS_scores <- read_tsv(file = paste0(input_dir, "/OT-MONDO_0005301-associated-targets-27_05_2025-v25_03.tsv")) %>% as.data.frame()
ALS_scores <- read_tsv(file = paste0(input_dir, "/OT-MONDO_0004976-associated-targets-27_05_2025-v25_03.tsv")) %>% as.data.frame()

## Subset to coloc genes (not all have OTAR scores)
AD_coloc_genes_previous_scores <- subset(coloc_genes, GWAS == "AD")
AD_scores_subset <- AD_scores[AD_scores$symbol %in% AD_coloc_genes_previous_scores$eGene, 1:3] %>% rename(eGene = symbol)
AD_coloc_genes_previous_scores <- AD_coloc_genes_previous_scores %>% 
  left_join(AD_scores_subset, by = "eGene")

PD_coloc_genes_previous_scores <- subset(coloc_genes, GWAS == "PD")
PD_scores_subset <- PD_scores[PD_scores$symbol %in% PD_coloc_genes_previous_scores$eGene, 1:3] %>% rename(eGene = symbol)
PD_coloc_genes_previous_scores <- PD_coloc_genes_previous_scores %>% 
  left_join(PD_scores_subset, by = "eGene")

MS_coloc_genes_previous_scores <- subset(coloc_genes, GWAS == "MS")
MS_scores_subset <- MS_scores[MS_scores$symbol %in% MS_coloc_genes_previous_scores$eGene, 1:3] %>% rename(eGene = symbol)
MS_coloc_genes_previous_scores <- MS_coloc_genes_previous_scores %>% 
  left_join(MS_scores_subset, by = "eGene")

ALS_coloc_genes_previous_scores <- subset(coloc_genes, GWAS == "ALS")
ALS_scores_subset <- ALS_scores[ALS_scores$symbol %in% ALS_coloc_genes_previous_scores$eGene, 1:3] %>% rename(eGene = symbol)
ALS_coloc_genes_previous_scores <- ALS_coloc_genes_previous_scores %>% 
  left_join(ALS_scores_subset, by = "eGene")


## Compare previous to new scores
AD_global_score_diffs <- as.numeric(AD_coloc_genes_previous_scores$globalScore) - as.numeric(AD_coloc_genes_previous_scores$OTAR_globalScore)
AD_gwas_score_diffs <- as.numeric(AD_coloc_genes_previous_scores$gwasCredibleSets) - as.numeric(AD_coloc_genes_previous_scores$OTAR_gwasCredibleSets)

PD_global_score_diffs <- as.numeric(PD_coloc_genes_previous_scores$globalScore) - as.numeric(PD_coloc_genes_previous_scores$OTAR_globalScore)
PD_gwas_score_diffs <- as.numeric(PD_coloc_genes_previous_scores$gwasCredibleSets) - as.numeric(PD_coloc_genes_previous_scores$OTAR_gwasCredibleSets)

MS_global_score_diffs <- as.numeric(MS_coloc_genes_previous_scores$globalScore) - as.numeric(MS_coloc_genes_previous_scores$OTAR_globalScore)
MS_gwas_score_diffs <- as.numeric(MS_coloc_genes_previous_scores$gwasCredibleSets) - as.numeric(MS_coloc_genes_previous_scores$OTAR_gwasCredibleSets)

ALS_global_score_diffs <- as.numeric(ALS_coloc_genes_previous_scores$globalScore) - as.numeric(ALS_coloc_genes_previous_scores$OTAR_globalScore)
ALS_gwas_score_diffs <- as.numeric(ALS_coloc_genes_previous_scores$gwasCredibleSets) - as.numeric(ALS_coloc_genes_previous_scores$OTAR_gwasCredibleSets)


## Max diff
na.omit(c(AD_global_score_diffs, PD_global_score_diffs, MS_global_score_diffs, ALS_global_score_diffs)) %>% abs %>% max
# [1] 0.1839296
na.omit(c(AD_gwas_score_diffs, PD_gwas_score_diffs, MS_gwas_score_diffs, ALS_gwas_score_diffs)) %>% abs %>% max
# [1] 0.1079528


## Use new score or previous if not reported !

AD_coloc_genes_previous_scores = 
  AD_coloc_genes_previous_scores %>%   
  mutate(globalScore_toReport = 
           case_when(
             globalScore == "No data" | is.na(globalScore) | globalScore == "NA" ~ as.character(OTAR_globalScore), 
             .default = as.character(globalScore))
  ) %>% 
  mutate(gwasCredibleSets_toReport = 
           case_when(
             gwasCredibleSets == "No data" | is.na(gwasCredibleSets) | gwasCredibleSets == "NA" ~ as.character(OTAR_gwasCredibleSets), 
             .default = as.character(gwasCredibleSets))
  )


PD_coloc_genes_previous_scores = 
  PD_coloc_genes_previous_scores %>%   
  mutate(globalScore_toReport = 
           case_when(
             globalScore == "No data" | is.na(globalScore) | globalScore == "NA" ~ as.character(OTAR_globalScore), 
             .default = as.character(globalScore))
  ) %>% 
  mutate(gwasCredibleSets_toReport = 
           case_when(
             gwasCredibleSets == "No data" | is.na(gwasCredibleSets) | gwasCredibleSets == "NA" ~ as.character(OTAR_gwasCredibleSets), 
             .default = as.character(gwasCredibleSets))
  )


MS_coloc_genes_previous_scores = 
  MS_coloc_genes_previous_scores %>%   
  mutate(globalScore_toReport = 
           case_when(
             globalScore == "No data" | is.na(globalScore) | globalScore == "NA" ~ as.character(OTAR_globalScore), 
             .default = as.character(globalScore))
  ) %>% 
  mutate(gwasCredibleSets_toReport = 
           case_when(
             gwasCredibleSets == "No data" | is.na(gwasCredibleSets) | gwasCredibleSets == "NA" ~ as.character(OTAR_gwasCredibleSets), 
             .default = as.character(gwasCredibleSets))
  )


ALS_coloc_genes_previous_scores = 
  ALS_coloc_genes_previous_scores %>%   
  mutate(globalScore_toReport = 
           case_when(
             globalScore == "No data" | is.na(globalScore) | globalScore == "NA" ~ as.character(OTAR_globalScore), 
             .default = as.character(globalScore))
  ) %>% 
  mutate(gwasCredibleSets_toReport = 
           case_when(
             gwasCredibleSets == "No data" | is.na(gwasCredibleSets) | gwasCredibleSets == "NA" ~ as.character(OTAR_gwasCredibleSets), 
             .default = as.character(gwasCredibleSets))
  )


coloc_genes_toReport <- rbind(AD_coloc_genes_previous_scores, PD_coloc_genes_previous_scores, 
                              MS_coloc_genes_previous_scores, ALS_coloc_genes_previous_scores) %>% as.data.frame()

coloc_genes <- as.data.frame(coloc_genes)
rownames(coloc_genes) <- paste(coloc_genes$eGene, coloc_genes$GWAS, sep = "-")
rownames(coloc_genes_toReport) <- paste(coloc_genes_toReport$eGene, coloc_genes_toReport$GWAS, sep = "-")
coloc_genes_toReport <- coloc_genes_toReport[rownames(coloc_genes), ] 

## Add Disease column as abb (EFO ID)   
coloc_genes_toReport$Disease <- paste0(coloc_genes_toReport$GWAS, " (", coloc_genes_toReport$efo_id, ")")

write_excel_csv(coloc_genes_toReport, file = paste0(outdir, "/coloc_genes_with_OTRA_GWAS_toReport.csv"), col_names = T)



## % of coloc genes x disease above >0.2 in newly-defined gwasCredibleSets score
by(coloc_genes_toReport, coloc_genes_toReport[, "GWAS"], function(x){table(x[, "gwasCredibleSets_toReport"] >0.2)})

# coloc_genes_toReport[, "GWAS"]: AD
# FALSE  TRUE 
#    15    12 
# ------------------------------------------------------------------------------------------------ 
#   coloc_genes_toReport[, "GWAS"]: ALS
# TRUE 
#    4 
# ------------------------------------------------------------------------------------------------ 
#   coloc_genes_toReport[, "GWAS"]: MS
# FALSE  TRUE 
#    20    21 
# ------------------------------------------------------------------------------------------------ 
#   coloc_genes_toReport[, "GWAS"]: PD
# FALSE  TRUE 
#     2     5 

table(coloc_genes_toReport$gwasCredibleSets_toReport>0.2)
# FALSE  TRUE 
#    37    42 







## Reproducibility info
options(width = 120)
session_info()
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.4.3 (2025-02-28)
# os       macOS Sequoia 15.4.1
# system   x86_64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       Europe/London
# date     2025-05-28
# rstudio  2024.12.1+563 Kousa Dogwood (desktop)
# pandoc   NA
# quarto   1.5.57 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# ! package          * version date (UTC) lib source
# AnnotationDbi      1.68.0  2024-10-29 [1] Bioconductor 3.20 (R 4.4.1)
# Biobase            2.66.0  2024-10-29 [1] Bioconductor 3.20 (R 4.4.1)
# BiocFileCache      2.14.0  2024-10-29 [1] Bioconductor 3.20 (R 4.4.1)
# BiocGenerics       0.52.0  2024-10-29 [1] Bioconductor 3.20 (R 4.4.1)
# BiocManager      * 1.30.25 2024-08-28 [1] CRAN (R 4.4.1)
# biomaRt          * 2.62.1  2025-01-30 [1] Bioconductor 3.20 (R 4.4.2)
# Biostrings         2.74.1  2024-12-16 [1] Bioconductor 3.20 (R 4.4.2)
# bit                4.6.0   2025-03-06 [1] CRAN (R 4.4.1)
# bit64              4.6.0-1 2025-01-16 [1] CRAN (R 4.4.1)
# blob               1.2.4   2023-03-17 [1] CRAN (R 4.4.0)
# cachem             1.1.0   2024-05-16 [1] CRAN (R 4.4.0)
# callr              3.7.6   2024-03-25 [1] CRAN (R 4.4.0)
# cellranger         1.1.0   2016-07-27 [1] CRAN (R 4.4.0)
# cli                3.6.4   2025-02-13 [1] CRAN (R 4.4.1)
# codetools          0.2-20  2024-03-31 [2] CRAN (R 4.4.3)
# colorspace         2.1-1   2024-07-26 [1] CRAN (R 4.4.0)
# crayon             1.5.3   2024-06-20 [1] CRAN (R 4.4.0)
# curl               6.2.2   2025-03-24 [1] CRAN (R 4.4.1)
# DBI                1.2.3   2024-06-02 [1] CRAN (R 4.4.0)
# dbplyr             2.5.0   2024-03-19 [1] CRAN (R 4.4.0)
# desc               1.4.3   2023-12-10 [1] CRAN (R 4.4.0)
# devtools           2.4.5   2022-10-11 [1] CRAN (R 4.4.0)
# digest             0.6.37  2024-08-19 [1] CRAN (R 4.4.1)
# dplyr            * 1.1.4   2023-11-17 [1] CRAN (R 4.4.0)
# ellipsis           0.3.2   2021-04-29 [1] CRAN (R 4.4.0)
# fastmap            1.2.0   2024-05-15 [1] CRAN (R 4.4.0)
# filelock           1.0.3   2023-12-11 [1] CRAN (R 4.4.0)
# forcats          * 1.0.0   2023-01-29 [1] CRAN (R 4.4.0)
# fs                 1.6.5   2024-10-30 [1] CRAN (R 4.4.1)
# generics           0.1.3   2022-07-05 [1] CRAN (R 4.4.0)
# GenomeInfoDb       1.42.3  2025-01-27 [1] Bioconductor 3.20 (R 4.4.2)
# GenomeInfoDbData   1.2.13  2025-05-27 [1] Bioconductor
# ggplot2          * 3.5.1   2024-04-23 [1] CRAN (R 4.4.0)
# V glue               1.8.0   2025-05-27 [1] Github (jimhester/fstrings@a3f80d6) (on disk 1.8.0.9000)
# gtable             0.3.6   2024-10-25 [1] CRAN (R 4.4.1)
# hms                1.1.3   2023-03-21 [1] CRAN (R 4.4.0)
# htmltools          0.5.8.1 2024-04-04 [1] CRAN (R 4.4.0)
# htmlwidgets        1.6.4   2023-12-06 [1] CRAN (R 4.4.0)
# httpuv             1.6.15  2024-03-26 [1] CRAN (R 4.4.0)
# httr             * 1.4.7   2023-08-15 [1] CRAN (R 4.4.0)
# httr2              1.1.2   2025-03-26 [1] CRAN (R 4.4.1)
# IRanges            2.40.1  2024-12-05 [1] Bioconductor 3.20 (R 4.4.2)
# jsonlite           2.0.0   2025-03-27 [1] CRAN (R 4.4.1)
# KEGGREST           1.46.0  2024-10-29 [1] Bioconductor 3.20 (R 4.4.1)
# later              1.4.2   2025-04-08 [1] CRAN (R 4.4.1)
# lifecycle          1.0.4   2023-11-07 [1] CRAN (R 4.4.0)
# lubridate        * 1.9.4   2024-12-08 [1] CRAN (R 4.4.1)
# magrittr           2.0.3   2022-03-30 [1] CRAN (R 4.4.0)
# memoise            2.0.1   2021-11-26 [1] CRAN (R 4.4.0)
# mime               0.13    2025-03-17 [1] CRAN (R 4.4.1)
# miniUI             0.1.1.1 2018-05-18 [1] CRAN (R 4.4.0)
# munsell            0.5.1   2024-04-01 [1] CRAN (R 4.4.0)
# pillar             1.10.2  2025-04-05 [1] CRAN (R 4.4.1)
# pkgbuild           1.4.7   2025-03-24 [1] CRAN (R 4.4.1)
# pkgconfig          2.0.3   2019-09-22 [1] CRAN (R 4.4.0)
# pkgload            1.4.0   2024-06-28 [1] CRAN (R 4.4.0)
# png                0.1-8   2022-11-29 [1] CRAN (R 4.4.0)
# prettyunits        1.2.0   2023-09-24 [1] CRAN (R 4.4.0)
# processx           3.8.6   2025-02-21 [1] CRAN (R 4.4.1)
# profvis            0.4.0   2024-09-20 [1] CRAN (R 4.4.1)
# progress           1.2.3   2023-12-06 [1] CRAN (R 4.4.0)
# promises           1.3.2   2024-11-28 [1] CRAN (R 4.4.1)
# ps                 1.9.0   2025-02-18 [1] CRAN (R 4.4.1)
# purrr            * 1.0.4   2025-02-05 [1] CRAN (R 4.4.1)
# R6                 2.6.1   2025-02-15 [1] CRAN (R 4.4.1)
# rappdirs           0.3.3   2021-01-31 [1] CRAN (R 4.4.0)
# Rcpp               1.0.14  2025-01-12 [1] CRAN (R 4.4.1)
# readr            * 2.1.5   2024-01-10 [1] CRAN (R 4.4.0)
# readxl             1.4.5   2025-03-07 [1] CRAN (R 4.4.1)
# remotes            2.5.0   2024-03-17 [1] CRAN (R 4.4.0)
# rlang            * 1.1.5   2025-01-17 [1] CRAN (R 4.4.1)
# RSQLite            2.3.11  2025-05-04 [1] CRAN (R 4.4.1)
# rstudioapi         0.17.1  2024-10-22 [1] CRAN (R 4.4.1)
# S4Vectors          0.44.0  2024-10-29 [1] Bioconductor 3.20 (R 4.4.1)
# scales             1.3.0   2023-11-28 [1] CRAN (R 4.4.0)
# sessioninfo      * 1.2.3   2025-02-05 [1] CRAN (R 4.4.1)
# shiny              1.10.0  2024-12-14 [1] CRAN (R 4.4.1)
# stringi            1.8.7   2025-03-27 [1] CRAN (R 4.4.1)
# stringr          * 1.5.1   2023-11-14 [1] CRAN (R 4.4.0)
# tibble           * 3.2.1   2023-03-20 [1] CRAN (R 4.4.0)
# tidyr            * 1.3.1   2024-01-24 [1] CRAN (R 4.4.0)
# tidyselect         1.2.1   2024-03-11 [1] CRAN (R 4.4.0)
# tidyverse        * 2.0.0   2023-02-22 [1] CRAN (R 4.4.0)
# timechange         0.3.0   2024-01-18 [1] CRAN (R 4.4.0)
# tzdb               0.5.0   2025-03-15 [1] CRAN (R 4.4.1)
# UCSC.utils         1.2.0   2024-10-29 [1] Bioconductor 3.20 (R 4.4.1)
# urlchecker         1.0.1   2021-11-30 [1] CRAN (R 4.4.0)
# usethis            3.1.0   2024-11-26 [1] CRAN (R 4.4.1)
# utf8               1.2.4   2023-10-22 [1] CRAN (R 4.4.0)
# vctrs              0.6.5   2023-12-01 [1] CRAN (R 4.4.0)
# vroom              1.6.5   2023-12-05 [1] CRAN (R 4.4.0)
# withr              3.0.2   2024-10-28 [1] CRAN (R 4.4.1)
# xml2               1.3.8   2025-03-14 [1] CRAN (R 4.4.1)
# xtable             1.8-4   2019-04-21 [1] CRAN (R 4.4.1)
# XVector            0.46.0  2024-10-29 [1] Bioconductor 3.20 (R 4.4.1)
# zlibbioc           1.52.0  2024-10-29 [1] Bioconductor 3.20 (R 4.4.1)
# 
# [1] /Users/dg40/Library/R/x86_64/4.4/library
# [2] /Library/Frameworks/R.framework/Versions/4.4-x86_64/Resources/library
# 
# * ── Packages attached to the search path.
# V ── Loaded and on-disk version mismatch.
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

