library(tidyverse)
library(ensembldb)
library(locuszoomr)
library(patchwork)
library(AnnotationHub)
source("functions.R")

# Global options
options(future.globals.maxSize = 90000 * 1024^2) # 90Gb
options(ucscChromosomeNames = FALSE) # Allow arbitrary chromosome identifiers
options(scipen = 999) # Prevent small numbers from showing as 0
readRenviron(path = "/nfs/users/nfs_m/ma23/.Renviron") # Use LDlinkR token

# Load data
dirs <- list.dirs("../../data/results/8.colocalisation_analysis/coloc_results")[-1]
genes <- read_csv("/lustre/scratch123/hgi/teams/trynka/resources/biomart/Homo_sapiens.GRCh38.111.genes.csv", col_names = TRUE) %>%
  select(gene_id, gene_name)
dbsnp <- read_delim("../../data/genotype/dbsnp/dbsnp.156.tensorqtl.subset.tsv", col_names = FALSE, delim = " ") %>%
  rename(chrom = X1, pos = X2, rsid = X3, ref = X4, alt = X5) %>%
  mutate(variant_id = paste(chrom, pos, ref, alt, sep = "_")) %>%
  distinct(variant_id, .keep_all = TRUE) %>%
  select(rsid, variant_id)
ah <- AnnotationHub()
ensDb_v106 <- ah[["AH100643"]] # GRCh38 genome build

# Helper functions
plot_and_save <- function(condition, locus_name, eGene, genes, dbsnp, outDir, gwas_dataset = NULL, nudge_y = 1) {
  pdf_file <- paste0(outDir, "/", condition, "_", str_replace(locus_name, "/", "-"), "_", eGene, "_selected_coloc_plots.pdf")
  pdf(file = pdf_file, width = 8, height = 8)
  
  plist <- plot_selected_coloc_locuszoom(condition, locus_name, eGene, genes, dbsnp, nudge_y=1)
  plot(plist$p)
  
  if (!is.null(gwas_dataset)) {
    eqtl_index_snp <- plist$eQTL %>% dplyr::filter(pval_nominal == min(pval_nominal)) %>% pull(rsid)
    plist_gwas <- plot_selected_GWAS_locuszoom(gwas_dataset, gene = eGene, locus_name, eqtl_index_snp, nudge_y)
    plot(plist_gwas$p)
    rm(plist_gwas)
    gc()
  }
  
  dev.off()
}

process_locus <- function(outDir, locus_name, eGenes, coloc_res, genes, dbsnp, gwas_datasets = NULL) {
  coloc_res <- read_rds(paste0(outDir, "/coloc_results_single_causal_variant_", 500000 / 1000, "kb.rds"))
  
  if (length(coloc_res) == 0) {
    stop("Error: coloc_res is empty!")
  }
  
  message("Processing locus: ", locus_name)  # Debugging message
  
  # First loop: Iterate over each condition in coloc_res
  for (condition in names(coloc_res)) {
    message("Processing condition: ", condition)  # Debugging message
    
    # Second loop: Iterate over each eGene
    for (eGene in eGenes) {
      # Check if gwas_datasets exists for the current condition
      gwas_dataset <- if (!is.null(gwas_datasets) && condition %in% names(gwas_datasets)) gwas_datasets[[condition]] else NULL
      
      message("Processing eGene: ", eGene, " | GWAS dataset: ", ifelse(is.null(gwas_dataset), "None", gwas_dataset))  # Debugging message
      if(is.null(coloc_res[[condition]][[locus_name]][[eGene]])){
        message("eGene ",eGene,
                " is not present in the coloc results for locus ", locus_name,
                " and condition ", condition,". Skipping...")
        next()
      }
      # Call the plotting function
      if(is.null(gwas_dataset)){
        plot_and_save(condition, locus_name, eGene, genes, dbsnp, outDir)
        
      }else{
        plot_and_save(condition, locus_name, eGene, genes, dbsnp, outDir, gwas_dataset)
        
      }
    }
  }
}


### AD GWAS loci
# Process TREM2 locus
outDir <- "../../data/results/8.colocalisation_analysis/coloc_results/GCST90027158" # Bellenguez AD
locus_name <- "6_41181270_TREM2"
eGene <- "TREM2"
gwas_datasets <- list("63_untreated_Not_proliferating" = "GCST90027158") # AD Bellenguez

process_locus(outDir, locus_name, eGene, coloc_res, genes, dbsnp, gwas_datasets)

# Process PILRA/B locus
locus_name <- "7_100334426_ZCWPW1/NYAP1"
eGenes <- c("PILRA", "PILRB")
gwas_datasets <- list("73_IFN_Not_proliferating" = "GCST90027158") # AD Bellenguez

process_locus(outDir, locus_name, eGenes, coloc_res, genes, dbsnp, gwas_datasets)

# Process EPHA1 locus
locus_name <- "7_143413669_EPHA1"
eGenes <- c("GSTK1", "FAM131B", "ZYX")
gwas_datasets <- list("63_untreated_Not_proliferating" = "GCST90027158") # AD Bellenguez

process_locus(outDir, locus_name, eGenes, coloc_res, genes, dbsnp, gwas_datasets)

# Process LRRK2 PD GWAS
outDir <- "../../data/results/8.colocalisation_analysis/coloc_results/ieu-b-7_PD_Nalls_2019/" # Nalls PD
locus_names <- c("12_40220632_LRRK2", "12_40340400_LRRK2")
eGenes <- c("LRRK2")

purrr::walk(locus_names, ~ {
  process_locus(outDir, .x, eGenes, coloc_res, genes, dbsnp)
})

