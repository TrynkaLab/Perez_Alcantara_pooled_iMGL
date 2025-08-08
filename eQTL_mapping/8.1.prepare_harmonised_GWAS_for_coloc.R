# colocalization analysis
# check relevant datasets here
# https://docs.google.com/spreadsheets/d/1jcB7UG1wx9SBtBKgjyE9CedRofMyaiUWBRJv0HKBVg4/edit#gid=0
# https://cran.r-project.org/web/packages/coloc/vignettes/a03_enumeration.html
library(tidyverse)
library(tidyr)
library(stringr)
library(purrr)

options(scipen = 999) # prevents from showing very small numbers as 0



args = commandArgs(trailingOnly=TRUE)

# Checking arguments
if (length(args)<4) {
  stop("4 arguments must be supplied", call. = FALSE)
}

GWAS_folder = as.character(args[1]) # folder that contains the harmonised summary stats from selected GWAS trait/study
eQTL_path = as.character(args[2]) # file path of the eQTL results of interest (created in 4.Inspect_eQTL_results.R)
eQTL_nominal_path=as.character(args[3]) # the nominal eQTL results for same conditions (see Snakefile)
distance_threshold=as.numeric(args[4]) # Define the distance threshold around lead variant in bp (e.g., 500,000 bp)



outDir = "../../data/results/8.colocalisation_analysis/subsets_for_coloc/"
dir.create(outDir)
# 
# # test data
# GWAS_folder = as.character("/lustre/scratch123/hgi/teams/trynka/resources/summary_statistics/public/GCST90027158")
# eQTL_path = as.character("../../data/results/4.Inspect_eQTL_results/tensorQTL_variant_gene_35PCs.csv")
# eQTL_nominal_path=as.character("../../data/results/tensorqtl/35/")
# distance_threshold = 500000

# GWAS data
GWAS_path = list.files(paste0(GWAS_folder,"/harmonised"),pattern = "*.h.tsv.gz",full.names = TRUE) #there should be just one file with that pattern
# omit cellect version
GWAS_path = GWAS_path[str_detect(string =GWAS_path, "cellect",negate = TRUE)]
  
# GWAS data
gwas = readr::read_tsv(GWAS_path) %>%
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
# prepared in 8.prepare_lead_GWAS_variant....
### should be GRCh38 and have chr, snp_pos, variant_id and locus_name 
GWAS_lead_variants_path = list.files(as.character(paste0(outDir,basename(GWAS_folder))),pattern = "*GRCh38_loci*",full.names = TRUE)

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

## eQTL data
message("Reading in eQTL results from:")
message(eQTL_path)
eqtl_full = readr::read_csv(eQTL_path)
message("Head of the eQTL result file:")
print(head(eqtl_full))

# run per condition
for(comp in unique(eqtl_full$group)){
  
  eqtl = eqtl_full %>%
    dplyr::rename(beta = slope) %>% # the slope is the beta
    tidyr::separate_wider_delim(cols =  position, names = c("chr", "pos"), delim = "_") %>%
    dplyr::mutate(varbeta = slope_se^2,
                  chr = as.numeric(chr),
                  pos=as.numeric(pos)) %>%
    dplyr::filter(group == comp) %>%
    dplyr::select(chr, pos,gene_id, gene_name, variant_id , beta,varbeta,pval_nominal, qval) %>%
    dplyr::filter(complete.cases(.)) %>% # will remove rows with NA, which are those from the X chr
    #group_by(variant_id) %>% # select lead associations per variant
    #dplyr::filter(pval_nominal == min(pval_nominal)) %>%
    #ungroup() %>%
    dplyr::arrange(chr,pos) 
  
  message("Printing colnames of eqtl file:")
  print( colnames(eqtl))
  
  # nominal results
  # read in from each condition (treatment x cluster)
  eqtl_nominal=read_tsv(paste0(eQTL_nominal_path,
                               "sum_sizefactorsNorm_log2_scaled_centered_",
                               str_split_fixed(comp,"_",2)[2],
                               "_common_500kb_window_tensorQTL_nominal.txt")) %>%
    dplyr::rename(beta = slope, # the slope is the beta
                  gene_id = phenotype_id
    ) %>%
    tidyr::separate_wider_delim(cols =  variant_id, names = c("chr", "pos","minor_allele","major_allele"), delim = "_",cols_remove = FALSE) %>%
    dplyr::mutate(varbeta = slope_se^2,
                  chr = as.numeric(chr),
                  pos=as.numeric(pos)) %>%
    dplyr::select(chr, pos,minor_allele,major_allele, gene_id, variant_id , beta,varbeta,pval_nominal) %>%
    dplyr::filter(complete.cases(.)) %>% # will remove rows with NA, which are those from the X chr
    #group_by(variant_id) %>% # select lead associations per variant
    #plyr::filter(pval_nominal == min(pval_nominal)) %>%
    #ungroup() %>%
    dplyr::arrange(chr,pos)
  
  ### gathering eqtl and lead GWAS variants that are close enough #####
  
  # from previous paper: https://www.nature.com/articles/s41588-022-01066-3#Sec9
  # For each trait–cell type pair, we applied colocalization to any locus where a 
  # lead variant for a significant eQTL (q value < 0.1) was located within 100kb (100,000 bp) 
  # and in high LD (r2 > 0.5) with a significant GWAS variant (i.e., any GWAS variant
  # with nominal P value < 1 ×10−5, which enabled us to capture suggestive association 
  # signals). In addition, we required at least 50 variants to be available for 
  # testing at each candidate locus. 
  
  # from jeremy's paper:
  # colocalisation tests between GWAS and eQTL signals where the lead variants were within 500 kb of each other, 
  # and passed to coloc all variants within 200 kb of each lead variant.
  
  # subset GWAS summary at same location, and for shared variants with eqtl
  
  subset_gwas = list()
  for(n_locus in 1:nrow(gwas_lead)){
    locus = paste(gwas_lead$chr[n_locus],gwas_lead$snp_pos[n_locus], gwas_lead$locus_name[n_locus], sep = "_")
    
    subset_gwas[[locus]] = gwas %>%
      dplyr::filter(chr == gwas_lead[n_locus,]$chr & pos <  (gwas_lead[n_locus,]$snp_pos + distance_threshold) & pos >  (gwas_lead[n_locus,]$snp_pos - distance_threshold))
    subset_gwas[[locus]] = subset_gwas[[locus]]  %>%
      dplyr::filter(!is.na(beta )) 
    
    
  }
  # For every eQTL lead variant close enough to GWAS locus, subset full eqtl summary based on positions that 
  # are within the distance threshold to GWAS lead
  subset_eqtl = list()
  subset_eqtl_split=list()
  subset_gwas_split = list()
  
  # gwas_lead must have the following columns:
  # chr   snp_pos variant_id      locus_name
  for (n_locus in 1:nrow(gwas_lead)) {
    locus = paste(gwas_lead$chr[n_locus],
                  gwas_lead$snp_pos[n_locus],
                  gwas_lead$locus_name[n_locus],
                  sep = "_")
    subset_eqtl[[locus]] = eqtl %>%
      dplyr::group_by(gene_name) %>% # group by eQTL locus and extract the lead variants that reach significance
      dplyr::filter(qval == min(qval)) %>% # doesn't do any filtering because cis_map already gives the lead eQTL variants per gene
      dplyr::filter(qval < 0.05) %>% # filter for significance
      dplyr::filter(
        chr == gwas_lead[n_locus,]$chr &
          pos <  (gwas_lead[n_locus,]$snp_pos + distance_threshold) &
          pos >  (gwas_lead[n_locus,]$snp_pos - distance_threshold)
      ) %>%
      dplyr::ungroup()
    
    if (nrow(subset_eqtl[[locus]]) == 0) {
      message(
        "There are no significant eQTL variants 500kb around the lead GWAS variant in locus ",
        locus
      )
      subset_eqtl[[locus]] = "There are no significant eQTL variants 500kb around the lead GWAS variant"
      
    } else{
      # now subset using variants from the nominal tests, per gene that remains significant from the previous step
      # also at specified distance from gwas variant
      
      sign_egenes = subset_eqtl[[locus]]$gene_id
      names(sign_egenes) = subset_eqtl[[locus]]$gene_name
      ####
      
      subset_eqtl[[locus]] = subset_eqtl[[locus]] %>%
        dplyr::right_join(eqtl_nominal) %>% ###
        dplyr::filter(gene_id %in% sign_egenes) %>% # It's unclear why I'm losing some significant eQTL variants that are missing in nominal - maybe nominal results and significant are different versions of analyses
        dplyr::mutate(gene_name = names(sign_egenes[match(gene_id, sign_egenes)])) %>%
        dplyr::left_join(eqtl[c("variant_id", "qval", "gene_name")],
                         by = c("variant_id", "gene_name"),
                         keep = TRUE) %>%
        dplyr::rename(
          qval_lead_variant_eqtl = qval.y,
          variant_id = variant_id.x,
          gene_name = gene_name.x
        ) %>%
        dplyr::select(!c("variant_id.y", "gene_name.y", "qval.x"))
      
      if (nrow(subset_eqtl[[locus]]) == 0) {
        message(
          "Warning:  significant eGenes (",
          sign_egenes,
          ") are missing in nominal for locus ",
          locus
        )
        subset_eqtl[[locus]] = paste0(
          "Warning:  significant eGenes (",
          sign_egenes,
          ") are missing in nominal for locus ",
          locus
        )
        
      } else{
        # split by eGene to test correct betas
        subset_eqtl_split[[locus]] = subset_eqtl[[locus]] %>%
          dplyr::group_split(gene_name)
        
        # naming sub-lists
        gene_gene_names = c()
        for (n in 1:length(subset_eqtl_split[[locus]])) {
          gene_gene_names = append(gene_gene_names, unique(subset_eqtl_split[[locus]][[n]]$gene_name))
          subset_gwas_split[[locus]][[n]] = subset_gwas[[locus]]
        }
        names(subset_eqtl_split[[locus]]) = gene_gene_names
        names(subset_gwas_split[[locus]]) = gene_gene_names
        
        ### pre-check if there are swapped REF/ALTs and if there are, swap variant_id and betas
        # from eQTL
        for (n in names(subset_eqtl_split[[locus]])) {
          
        temp=subset_eqtl_split[[locus]][[n]] %>%
          dplyr::mutate(swapped_variant_id = paste(chr,pos,
                                                   str_split_i(variant_id,pattern="_",i=4),str_split_i(variant_id,pattern="_",i=3),
                                                   sep = "_")) %>%
          dplyr::mutate(variant_id = case_when(swapped_variant_id %in% subset_gwas_split[[locus]][[n]]$variant_id ~ swapped_variant_id,
                                               .default = variant_id),
                        swapped = case_when(swapped_variant_id %in% subset_gwas_split[[locus]][[n]]$variant_id ~ TRUE,
                                               .default = FALSE),
                        beta = case_when(swapped_variant_id %in% subset_gwas_split[[locus]][[n]]$variant_id ~ -beta,
                                               .default = beta))
        message("Table of swapped alleles: ")
        table(temp$swapped)
        
        subset_eqtl_split[[locus]][[n]] = temp %>%
          dplyr::select(-swapped_variant_id,-swapped)
        }
        
        # check which dataframe is bigger, GWAS or eQTL, then subset to shared variants
        for (n in names(subset_eqtl_split[[locus]])) {
          if (nrow(subset_eqtl_split[[locus]][[n]]) > nrow(subset_gwas_split[[locus]][[n]])) {
            subset_gwas_split[[locus]][[n]] = subset_gwas_split[[locus]][[n]] %>%
              dplyr::filter(variant_id %in% subset_eqtl_split[[locus]][[n]]$variant_id)
            subset_eqtl_split[[locus]][[n]] = subset_eqtl_split[[locus]][[n]] %>%
              dplyr::filter(variant_id %in% subset_gwas_split[[locus]][[n]]$variant_id)
          } else{
            subset_eqtl_split[[locus]][[n]] = subset_eqtl_split[[locus]][[n]] %>%
              dplyr::filter(variant_id %in% subset_gwas_split[[locus]][[n]]$variant_id)
            subset_gwas_split[[locus]][[n]] = subset_gwas_split[[locus]][[n]] %>%
              dplyr::filter(variant_id %in% subset_eqtl_split[[locus]][[n]]$variant_id)
          }
          check_rows <- function(df1, df2) {
            if (nrow(df1) != nrow(df2)) {
              stop("Error: The eQTL and GWAS data frames don't have an identical number of rows.")
            }
          }
          tryCatch({
            check_rows(subset_gwas_split[[locus]][[n]], subset_eqtl_split[[locus]][[n]])
          }, error = function(e) {
            cat("Error occurred:", e$message, "\n")
          })
          
          message(
            "There are ",
            sum(
              subset_gwas_split[[locus]][[n]]$variant_id %in% unique(subset_eqtl_split[[locus]][[n]]$variant_id)
            ),
            " shared variants between GWAS and eQTL at GWAS locus ",
            locus,
            " for eGene ",
            n
          )
        }
      }
    }
  }
  # subset GWAS list to contain only loci for which we have eQTL results
  subset_gwas_split = subset_gwas_split[names(subset_gwas_split) %in% names(subset_eqtl_split)]
  
  # write objects
  dir.create(paste0(outDir,basename(GWAS_folder)))
  
  readr::write_rds(x = subset_gwas_split,file = paste0(outDir,basename(GWAS_folder),"/",comp,"_GWAS_subsets.rds"))
  readr::write_rds(x = subset_eqtl_split,file = paste0(outDir,basename(GWAS_folder),"/",comp,"_eQTL_subsets.rds"))
}
