# for coloc - running on my own GWAS
# against external GWAS
library(tidyverse)
library(tidyr)
library(stringr)
library(purrr)

options(scipen = 999) # prevents from showing very small numbers as 0


### fix this to run on server
 args = commandArgs(trailingOnly=TRUE)
# 
 # Checking arguments
 if (length(args)<3) {
   stop("3 arguments must be supplied", call. = FALSE)
 }
# 
 external_GWAS_folder = as.character(args[1]) # folder that contains the harmonised summary stats from selected GWAS trait/study
 my_GWAS_file= as.character(args[2]) # file path of our GWAS results of interest 
 distance_threshold=as.numeric(args[3]) # Define the distance threshold around lead variant in bp (e.g., 500,000 bp)
##########

outDir = "../../../data/results/migration/3.coloc/subsets_for_coloc/"
dir.create(outDir, recursive = TRUE)

  # for testing
# external_GWAS_folder = as.character("/lustre/scratch123/hgi/teams/trynka/resources/summary_statistics/public/GCST90027158/") 
# my_GWAS_file = as.character("../../../data/results/phagocytosis/2.check_association_results/lm_1pct_filtered_deflated/untreated/untreated_GWAS_GRCh38.tsv.gz")
# distance_threshold = 500000

# GWAS data
  treat = stringr::str_split(my_GWAS_file,"/")[[1]][9]
  external_GWAS_path = list.files(paste0(external_GWAS_folder,"/harmonised"),pattern = "*.h.tsv.gz",full.names = TRUE) #there should be just one file with that pattern
  # omit cellect version
  external_GWAS_path = external_GWAS_path[str_detect(string =external_GWAS_path, "cellect",negate = TRUE)]  
# GWAS data
gwas = readr::read_tsv(external_GWAS_path) %>%
  dplyr::rename(standard_error = if_else(condition = "se" %in% colnames(.),true = "se", false = "standard_error")) %>%
  dplyr::select(starts_with("hm_var"), hm_beta,p_value,hm_effect_allele, hm_chrom,hm_pos,
                standard_error,standard_error,effect_allele_frequency) %>%
  dplyr::rename(variant_id = starts_with("hm_var"),
                beta = hm_beta, # otherwise beta = log(OR) https://github.com/chr1swallace/coloc/issues/149
                pval_nominal = p_value,
                effect_allele = hm_effect_allele,
                maf=effect_allele_frequency) %>%
  dplyr::mutate(chr = as.numeric(hm_chrom),
                pos=as.numeric(hm_pos), # ensuring types are correct
                varbeta = standard_error^2) %>%
  # no need to adjust for N, see https://github.com/chr1swallace/coloc/issues/14

  dplyr::select(chr, pos,variant_id,effect_allele, beta,varbeta,pval_nominal)


message("Printing colnames of GWAS file:")
print( colnames(gwas))

# loading lead variants file
# prepared in 8.prepare_lead_GWAS_variant....from eQTL folder
### should be GRCh38 and have chr, snp_pos, variant_id and locus_name 
GWAS_lead_variants_path = list.files(as.character(paste0("../../../../OTAR2065_sc_eQTL/data/results/8.colocalisation_analysis/subsets_for_coloc/",
                                                         basename(external_GWAS_folder))),pattern = "*GRCh38_loci*",full.names = TRUE)

gwas_lead =  readr::read_tsv(GWAS_lead_variants_path)

# check all lead variants are in the summary file
if(sum(!gwas_lead$variant_id %in% gwas$variant_id)==0){
  message("All variants from lead GWAS file are present in the GWAS summary file")
} else{
  message("WARNING: Some variants from lead GWAS file are not present in the GWAS summary file")
  message(gwas_lead$variant_id[!gwas_lead$variant_id %in% gwas$variant_id])
  message("Swapping the REF and ALT alleles in those variants")
  gwas_lead = gwas_lead %>%
    dplyr::mutate(variant_id = case_when(!variant_id %in% gwas$variant_id ~ paste0(chr,"_",snp_pos,"_",str_replace(`Minor/Major_allele`, "([^/]+)/([^/]+)", "\\1_\\2")),
                                         .default = variant_id))
  if(sum(!gwas_lead$variant_id %in% gwas$variant_id)==0){
    message("All variants from lead GWAS file are now present in the GWAS summary file")
  }else{
    message("Warning: Swapping the REF and ALT alleles did not work in all cases")
    
  }
}

## our GWAS data

message("Reading in GWAS results from:")
message(my_GWAS_file)
my_GWAS = readr::read_tsv(my_GWAS_file)
message("Head of my GWAS file:")
print(head(my_GWAS))
dim(my_GWAS) # ~ 3 million
### checking flipped alleles
my_GWAS = my_GWAS %>%
  dplyr::mutate(swapped_allele = paste(chr,snp_pos,ALT,REF,sep = "_")) %>%
  dplyr::mutate(variant_id = case_when(swapped_allele %in% gwas$variant_id ~ swapped_allele,
                                       .default = variant_id),
                coef =  case_when(swapped_allele %in% gwas$variant_id ~ -coef,
                                  .default = coef),
                swapped = case_when(swapped_allele %in% gwas$variant_id ~ TRUE,
                                    .default = FALSE))

  my_gwas_for_coloc = my_GWAS %>%
    dplyr::rename(beta = coef,pval_nominal= p_value) %>% # the coef is the beta
    dplyr::mutate(varbeta = se^2,
                  chr = as.numeric(chr),
                  pos=as.numeric(snp_pos)) %>%
    dplyr::select(chr, pos, variant_id , beta,varbeta,pval_nominal) %>%
    dplyr::arrange(chr,pos) 
  
  message("Printing colnames of my GWAS file:")
  print( colnames(my_gwas_for_coloc))
  
  ### gathering my GWAS and GWAS variants from external data that are close enough #####
  
  # subset GWAS summary at same location
  # taken from lead external GWAS loci
  
  subset_gwas = list()
  for(n_locus in 1:nrow(gwas_lead)){
    locus = paste(gwas_lead$chr[n_locus],gwas_lead$snp_pos[n_locus], gwas_lead$locus_name[n_locus], sep = "_")
    
    subset_gwas[[locus]] = gwas %>%
      dplyr::filter(chr == gwas_lead[n_locus,]$chr & pos <  (gwas_lead[n_locus,]$snp_pos + distance_threshold) & pos >  (gwas_lead[n_locus,]$snp_pos - distance_threshold))
    subset_gwas[[locus]] = subset_gwas[[locus]]  %>%
      dplyr::filter(!is.na(beta )) 
    
    
  }
  # Subset all my GWAS positions within +/- 500kb external GWAS leads
  subset_my_gwas = list()

  # gwas_lead must have the following columns:
  # chr   snp_pos variant_id      locus_name
  for (n_locus in 1:nrow(gwas_lead)) {
    # subset to shared variants
    locus = paste(gwas_lead$chr[n_locus],gwas_lead$snp_pos[n_locus], gwas_lead$locus_name[n_locus], sep = "_")
    
    subset_my_gwas[[locus]] = my_gwas_for_coloc %>%
      dplyr::filter(variant_id %in% subset_gwas[[locus]]$variant_id)
    subset_gwas[[locus]] = subset_gwas[[locus]] %>%
      dplyr::filter(variant_id %in% subset_my_gwas[[locus]]$variant_id)
    
      
          check_rows <- function(df1, df2) {
            if (nrow(df1) != nrow(df2)) {
              stop("Error: The eQTL and GWAS data frames don't have an identical number of rows.")
            }
          }
          tryCatch({
            check_rows(subset_gwas[[locus]], subset_my_gwas[[locus]])
          }, error = function(e) {
            cat("Error occurred:", e$message, "\n")
          })
          
          message(
            "There are ",
            sum(
              subset_gwas[[locus]]$variant_id %in% unique(subset_my_gwas[[locus]]$variant_id)
            ),
            " shared variants between external GWAS and my GWAS at locus ",
            locus
          )
        }
      
   
  # write objects
  dir.create(paste0(outDir,basename(external_GWAS_folder)))
  
  readr::write_rds(x = subset_gwas,file = paste0(outDir,basename(external_GWAS_folder),"/",treat,"_external_GWAS_subsets.rds"))
  readr::write_rds(x = subset_my_gwas,file = paste0(outDir,basename(external_GWAS_folder),"/",treat,"_my_GWAS_subsets.rds"))


