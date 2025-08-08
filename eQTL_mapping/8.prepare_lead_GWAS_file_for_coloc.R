# prepare lead GWAS file for 8.1.prepare_harmonised_GWAS_for_coloc.R

library(tidyverse)
library(tidyr)
library(stringr)
library(purrr)

options(scipen = 999) # prevents from showing very small numbers as 0

args = commandArgs(trailingOnly=TRUE)

# Checking arguments
if (length(args)<1) {
  stop("1 argument must be supplied", call. = FALSE)
}

GWAS_folder = as.character(args[1]) # folder that contains the harmonised summary stats from selected GWAS trait/study

## test
GWAS_folder="/lustre/scratch123/hgi/teams/trynka/resources/summary_statistics/public/GCST006900/"

###
outDir = "../../data/results/8.colocalisation_analysis/subsets_for_coloc/"
dir.create(paste0(outDir,basename(GWAS_folder)))


GWAS_lead_variants_path = list.files(as.character(GWAS_folder),pattern = "*GRCh38_loci*",full.names = TRUE)
if(!is_empty(GWAS_lead_variants_path)){
  

gwas_lead =  readr::read_tsv(GWAS_lead_variants_path)
# column names must be:
#   chr   snp_pos variant_id      locus_name

columns_to_rename = c("chr", "snp_pos", "variant_id", "locus_name")
if ("Gene" %in% colnames(gwas_lead)) {
  if (sum(colnames(gwas_lead) %in% columns_to_rename) == 4) {
    gwas_lead = gwas_lead %>%
      dplyr::mutate(locus_name = dplyr::case_when(locus_name == "New" ~ Gene,
                                                  .default = locus_name))
    write_tsv(gwas_lead, paste0(outDir, basename(GWAS_folder), "/GRCh38_loci.txt"))
    
  } 
  if("Minor/Major_allele" %in% colnames(gwas_lead)) {
    # change colnames to fit
    gwas_lead = gwas_lead %>%
      dplyr::rename_with( ~ tolower(.), matches(paste0("(?i)", columns_to_rename)))  # fix any uppercase issues
    
    # if loci don't have proper names, rename
    
    gwas_lead = gwas_lead %>%
      dplyr::mutate(locus_name = dplyr::case_when(locus_name == "New" ~ Gene,
                                                  .default = locus_name)) %>%
      # change variant ids
      dplyr::mutate(variant_id = paste0(
        chr,
        "_",
        snp_pos,
        "_",
        str_replace(`Minor/Major_allele`, "([^/]+)/([^/]+)", "\\2_\\1")
      ))
    
    write_tsv(gwas_lead, paste0(outDir, basename(GWAS_folder), "/GRCh38_loci.txt"))
    
  }
} 
if("DATE ADDED TO CATALOG"  %in% colnames(gwas_lead)){ # GWAS catalog format
  #   chr   snp_pos variant_id      locus_name
  
  gwas_lead = gwas_lead %>%
    dplyr::rename(chr = CHR_ID,
                  snp_pos = CHR_POS,
                  locus_name = MAPPED_GENE) %>%
    dplyr::select(chr,snp_pos,locus_name,SNPS)
  
  # load summary stats data and rename to chr_pos_ref_alt
  summary_file = list.files(paste0(GWAS_folder,"/harmonised/"),pattern = "*h.tsv.gz")
  # omit cellect files
  summary_file = summary_file[str_detect(summary_file,pattern = "cellect",negate = TRUE)]
  summary = read_tsv(paste0(GWAS_folder,"/harmonised/",summary_file))
  
  gwas_lead = gwas_lead %>%
    dplyr::left_join(summary[,c("hm_variant_id", "hm_rsid")],by = join_by(SNPS == hm_rsid)) %>%
    dplyr::filter(!is.na(hm_variant_id)) # some variants might disappear - GWAS catalog error?
  
  gwas_lead %>%
    dplyr::rename( variant_id= hm_variant_id) %>%
    dplyr::select(chr ,  snp_pos, variant_id  ,    locus_name) %>%
    write_tsv(., paste0(outDir, basename(GWAS_folder), "/GRCh38_loci.txt"))
  
} else {
  gwas_lead = gwas_lead %>%
    dplyr::rename_with( ~ tolower(.), matches(paste0("(?i)", columns_to_rename)))  # fix any uppercase issues
  
  
  write_tsv(gwas_lead, paste0(outDir, basename(GWAS_folder), "/GRCh38_loci.txt"))
}
} else{
  # create lead associations GWAS file 
  summary_path = list.files(paste0(as.character(GWAS_folder),"/harmonised"),pattern = "*.h.tsv.gz",full.names = TRUE)

  summary = read_tsv(summary_path)
  # filter to genome-wide significant associations
  summary = summary %>%
    dplyr::filter(p_value < 5*1e-8) # bonferroni correction for 1million indp. tests: 0.05/1000000
  
  # get LD file and retain variants NOT in LD with any other
  # defining LD as R > 0.7
  
  ld_file =list()
  for(chr in 1:22){
    ld_file[[chr]] = read_tsv(paste0("/lustre/scratch123/hgi/teams/trynka/resources/1000g/1000G/pairwise_ld/grch38_R_over_0.7/chr",
                               chr,"_GRCh38.EUR.ld.gz"))
  }

  ld_file = do.call("rbind",ld_file)
  
  # get the minimum value within independent variants
  summary = summary %>%
    dplyr::mutate(interval = floor(hm_pos / 250) * 250) %>%
    dplyr::group_by(interval) %>%
    dplyr::filter(p_value == min(p_value)) %>%
    dplyr::ungroup()
  
}
