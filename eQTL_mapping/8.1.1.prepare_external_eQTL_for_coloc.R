# prepare external eQTL data for coloc
# output file must be an rds containing a list with the eQTL tested categories 
# and within those, another list of loci from eQTL data that overlaps our loci
# and within those, a list of shared eGenes, and each of them is a tibble of variants
# with header:
# chr       pos variant_id       effect_allele    beta   varbeta pval_nominal
# <dbl>     <dbl> <chr>            <chr>           <dbl>     <dbl>        <dbl>

library(tidyverse)
library(tidyr)
library(stringr)
library(purrr)
library(qvalue)
source("./functions.R")
options(scipen = 999) # prevents from showing very small numbers as 0
options(future.globals.maxSize = 150000 * 1024^2) # 150Gb


args = commandArgs(trailingOnly=TRUE)

# Checking arguments
if (length(args)<3) {
  stop("3 arguments must be supplied", call. = FALSE)
}
external_eQTL_folder = as.character(args[1]) # folder that contains the harmonised summary stats from selected external eQTL trait/study
eQTL_path = as.character(args[2]) # file path of the eQTL results of interest (created in 4.Inspect_eQTL_results.R)
eQTL_nominal_path=as.character(args[3]) # the nominal eQTL results for same conditions (see Snakefile)


outDir = "../../data/results/8.colocalisation_analysis/subsets_for_coloc/"
dir.create(outDir)

# test data
# external_eQTL_folder = as.character("/lustre/scratch123/hgi/projects/otar2065/resources/macromap_Panousis_2023/eQTL/macromap_Panousis_2023_IFNG_6//") 
# eQTL_path = as.character("../../data/results/4.Inspect_eQTL_results/tensorQTL_variant_gene_60PCs.csv")
# eQTL_nominal_path=as.character("../../data/results/tensorqtl/60/")

external_eQTL_lead_path = list.files(external_eQTL_folder,pattern = "*permuted.tsv.gz",full.names = TRUE) 

# external eQTL data
external_eQTL_path = list.files(external_eQTL_folder,pattern = "*.all.tsv.gz",full.names = TRUE) 
# exclude index, if present
external_eQTL_path = external_eQTL_path[!grepl("tbi",external_eQTL_path)]

### read in
# external eQTL data - nominal
external_eqtl =  readr::read_tsv(external_eQTL_path) %>%
  
  dplyr::rename(standard_error = if_else(condition = "se" %in% colnames(.),true = "se", false = "standard_error")) %>%
  dplyr::rename(ref = REF, alt = ALT) %>%
  dplyr::mutate(variant = paste(chromosome,position,ref,alt,sep = "_")) %>%
  dplyr::select(ensembl_id,variant, beta,pval_nominal,alt, chromosome,position,
                standard_error) %>%
  dplyr::rename(variant_id = variant,
                effect_allele = alt) %>%
  dplyr::mutate(chr = chromosome,
                pos=position, # ensuring types are correct
                varbeta = standard_error^2) %>%
  # no need to adjust for N, see https://github.com/chr1swallace/coloc/issues/14
  
  dplyr::select(ensembl_id,chr, pos,variant_id,effect_allele, beta,varbeta,pval_nominal)

# loading lead variants file - if present
### should be GRCh38 and have chr, snp_pos, variant_id and locus_name 

if(!is_empty(external_eQTL_lead_path)){
  

external_eQTL_lead =  readr::read_tsv(external_eQTL_lead_path) %>%
  dplyr::rename(variant_id = variant,chr = chromosome,snp_pos = position) 
# the locus_name will be the gene name, since it's a cis analysis
genes = ensembl_to_gene_name(external_eQTL_lead$ensembl_id,IDTo = "SYMBOL",IDFrom = "GENEID") 

genes = genes$map
genes = genes %>%
  dplyr::rename(ensembl_id = From, locus_name = To)

external_eQTL_lead = external_eQTL_lead %>%
  dplyr::left_join(.,genes) %>%
  dplyr::mutate(locus_name = case_when(is.na(locus_name) ~ ensembl_id,
                                       .default = locus_name))  %>%
  # omit genes without significant eQTLs
  dplyr::filter(!is.na(variant_id)) %>%
  dplyr::filter(!is.na(p_beta)) %>%
  # calculating q values
  dplyr::mutate(qval =  qvalue(p_beta)$qvalues) %>%
  dplyr::filter(qval <0.05)

# subset full eQTL results to significant eGenes
external_eqtl = external_eqtl %>%
  dplyr::filter(ensembl_id %in% external_eQTL_lead$ensembl_id)
gc()
} else {
  external_eQTL_lead = readr::read_tsv(external_eQTL_path) %>%
    dplyr::rename(standard_error = if_else(condition = "se" %in% colnames(.),true = "se", false = "standard_error")) %>%
    dplyr::rename(ref = REF, alt = ALT) %>%
    dplyr::mutate(variant = paste(chromosome,position,ref,alt,sep = "_")) %>%
    dplyr::select(ensembl_id,variant, beta,pvalue,alt, chromosome,position,
                  standard_error,maf,significant_by_2step_FDR) %>%
    dplyr::rename(variant_id = variant,
                  pval_nominal = pvalue,
                  effect_allele = alt) %>%
    dplyr::mutate(chr = chromosome,
                  pos=position, # ensuring types are correct
                  varbeta = standard_error^2) %>%
    # no need to adjust for N, see https://github.com/chr1swallace/coloc/issues/14
    
    dplyr::select(ensembl_id,chr, pos,variant_id,effect_allele, beta,varbeta,pval_nominal,significant_by_2step_FDR) %>%
    dplyr::filter(significant_by_2step_FDR == "Yes")
  
  external_eQTL_lead = external_eQTL_lead %>%
    dplyr::group_by(ensembl_id) %>%
    dplyr::arrange(pval_nominal) %>%
    dplyr::slice_head(n=1) # take the lead variant per eGene as in tensorQTL
  
  # the locus_name will be the gene name, since it's a cis analysis
  genes = ensembl_to_gene_name(external_eQTL_lead$ensembl_id,IDTo = "SYMBOL",IDFrom = "GENEID") 
  
  genes = genes$map
  genes = genes %>%
    dplyr::rename(ensembl_id = From, locus_name = To)
  
  external_eQTL_lead = external_eQTL_lead %>%
    dplyr::left_join(.,genes) %>%
    dplyr::mutate(locus_name = case_when(is.na(locus_name) ~ ensembl_id,
                                         .default = locus_name))  %>%
    dplyr::rename(snp_pos = pos)
  
  # subset to significant eGenes
  external_eqtl = external_eqtl %>%
    dplyr::filter(ensembl_id %in% external_eQTL_lead$ensembl_id) %>%
    dplyr::distinct() 

}


message("Printing colnames of external eQTL file:")
print( colnames(external_eqtl))

message("There are ",length(unique(external_eQTL_lead$ensembl_id)), " unique significant eGenes") # 899
message("There are ",length(unique(external_eQTL_lead$variant_id)), " unique significant eQTLs") # 893

# check all lead variants are in the summary file
if(sum(!external_eQTL_lead$variant_id %in% external_eqtl$variant_id)==0){
  message("All variants from lead external eQTL file are present in the external eQTL summary file")
} else{
  message("WARNING: Some variants from lead external eQTL file are not present in the external eQTL summary file")
  message(external_eQTL_lead$variant_id[!external_eQTL_lead$variant_id %in% external_eqtl$variant_id])
  message("Swapping the REF and ALT alleles in those variants")
  if("Minor/Major_allele" %in% colnames(external_eQTL_lead)){
    external_eQTL_lead = external_eQTL_lead %>%
      dplyr::mutate(variant_id = case_when(!variant_id %in% external_eqtl$variant_id ~ paste0(chr,"_",snp_pos,"_",str_replace(`Minor/Major_allele`, "([^/]+)/([^/]+)", "\\1_\\2")),
                                           .default = variant_id))
  } else{
    external_eQTL_lead = external_eQTL_lead %>%
      dplyr::mutate(variant_id = case_when(!variant_id %in% external_eqtl$variant_id ~ paste0(chr,"_",snp_pos,"_",ALT,"_",REF),
                                           .default = variant_id))
  }
  
  if(sum(!external_eQTL_lead$variant_id %in% external_eqtl$variant_id)==0){
    message("All variants from lead external eQTL file are now present in the external eQTL summary file")
  }else{
    message("Warning: Swapping the REF and ALT alleles did not work in all cases")
    message("Subsetting to shared variants")
    external_eQTL_lead = external_eQTL_lead[external_eQTL_lead$variant_id %in% external_eqtl$variant_id,]
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
                               "/sum_sizefactorsNorm_log2_scaled_centered_",
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
  
  ### gathering eqtl and lead external eQTL variants from the same loci (same eGenes) ########
  
  # subset external eQTL summary at same location, and for shared variants with eqtl
  
  subset_external_eqtl =  list()
  for(n_locus in 1:nrow(external_eQTL_lead)){
    locus = paste(external_eQTL_lead$chr[n_locus],external_eQTL_lead$snp_pos[n_locus], external_eQTL_lead$locus_name[n_locus], sep = "_")
    
    subset_external_eqtl[[locus]] = external_eqtl %>%
      dplyr::filter(ensembl_id==external_eQTL_lead[n_locus,"ensembl_id"]$ensembl_id)
    subset_external_eqtl[[locus]] = subset_external_eqtl[[locus]]  %>%
      dplyr::filter(!is.na(beta )) 
    
    
  }
  # For every eQTL lead variant close enough to external eQTL locus, subset full eqtl summary based on positions that 
  # are within the distance threshold to external eQTL lead
  subset_eqtl = list()
  subset_eqtl_split=list()
  subset_external_eQTL_split = list()
  
  # external_eQTL_lead must have the following columns:
  # chr   snp_pos variant_id      locus_name
  for (n_locus in 1:nrow(external_eQTL_lead)) {
    locus = paste(external_eQTL_lead$chr[n_locus],
                  external_eQTL_lead$snp_pos[n_locus],
                  external_eQTL_lead$locus_name[n_locus],
                  sep = "_")
      subset_eqtl[[locus]] = eqtl %>%
        dplyr::filter(gene_id == external_eQTL_lead[n_locus,"ensembl_id"]$ensembl_id) %>% # prefilter for significanrt genes
        dplyr::group_by(gene_name) %>% # group by eQTL locus and extract the lead variants that reach significance
        dplyr::filter(qval == min(qval)) %>% # doesn't do any filtering because cis_map already gives the lead eQTL variants per gene
        dplyr::filter(qval < 0.05) %>% # filter for significance
        dplyr::ungroup()
    
    
    if (nrow(subset_eqtl[[locus]]) == 0) {
      message(
        "There are no significant eQTL variants linked to a significant external eGene in locus ",
        locus
      )
      subset_eqtl[[locus]] = "There are no significant eQTL variants 500kb around the lead external eQTL variant"
      
    } else{
      # now subset using variants from the nominal tests, per gene that remains significant from the previous step
      # also at specified distance from external eQTL variant
      
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
        
        # formatting external QTL variants
        subset_external_eqtl[[locus]] = subset_external_eqtl[[locus]] %>%
          dplyr::mutate(variant_id = str_replace(variant_id,pattern="chr",replacement=""))
        
        # naming sub-lists
        gene_names = c()
        for (n in 1:length(subset_eqtl_split[[locus]])) {
          gene_names = append(gene_names, unique(subset_eqtl_split[[locus]][[n]]$gene_name))
          subset_external_eQTL_split[[locus]][[n]] = subset_external_eqtl[[locus]]
        }
        names(subset_eqtl_split[[locus]]) = gene_names
        names(subset_external_eQTL_split[[locus]]) = gene_names
        

        # check which dataframe is bigger, external eQTL or eQTL, then 
        ###### subset to shared variants ########
        for (n in names(subset_eqtl_split[[locus]])) {
          #### check if any swapped REF ALT are present
          test = subset_eqtl_split[[locus]][[n]] %>%
            dplyr::mutate(swapped = paste(chr,pos,major_allele,minor_allele,sep="_"))
          if (sum(test$swapped %in% subset_external_eQTL_split[[locus]][[n]]$variant_id) > 0) {
            message(
              "Warning:  there are  ",
              sum(test$swapped %in% subset_external_eQTL_split[[locus]][[n]]$variant_id),
              " swapped alleles in locus ",
              locus
            )
            
          }
          
          if (nrow(subset_eqtl_split[[locus]][[n]]) > nrow(subset_external_eQTL_split[[locus]][[n]])) {
            subset_external_eQTL_split[[locus]][[n]] = subset_external_eQTL_split[[locus]][[n]] %>%
              dplyr::filter(variant_id %in% subset_eqtl_split[[locus]][[n]]$variant_id)  %>%
              distinct()
            subset_eqtl_split[[locus]][[n]] = subset_eqtl_split[[locus]][[n]] %>%
              dplyr::filter(variant_id %in% subset_external_eQTL_split[[locus]][[n]]$variant_id)  %>%
              distinct()
          } else{
            subset_eqtl_split[[locus]][[n]] = subset_eqtl_split[[locus]][[n]] %>%
              dplyr::filter(variant_id %in% subset_external_eQTL_split[[locus]][[n]]$variant_id)  %>%
              distinct()
            subset_external_eQTL_split[[locus]][[n]] = subset_external_eQTL_split[[locus]][[n]] %>%
              dplyr::filter(variant_id %in% subset_eqtl_split[[locus]][[n]]$variant_id) %>%
              distinct()
          }
          check_rows <- function(df1, df2) {
            if (nrow(df1) != nrow(df2)) {
              stop("Error: The eQTL and external eQTL data frames don't have an identical number of rows.")
            }
          }
          tryCatch({
            check_rows(subset_external_eQTL_split[[locus]][[n]], subset_eqtl_split[[locus]][[n]])
          }, error = function(e) {
            cat("Error occurred:", e$message, "\n")
          })
          
          message(
            "There are ",
            sum(
              subset_external_eQTL_split[[locus]][[n]]$variant_id %in% unique(subset_eqtl_split[[locus]][[n]]$variant_id)
            ),
            " shared variants between external eQTL and eQTL at external eQTL locus ",
            locus,
            " for eGene ",
            n
          )
        }
      }
    }
  }
  # subset external eQTL list to contain only loci for which we have eQTL results
  subset_external_eQTL_split = subset_external_eQTL_split[names(subset_external_eQTL_split) %in% names(subset_eqtl_split)]
  
  # write objects
  dir.create(paste0(outDir,"primary_eQTLs/",basename(external_eQTL_folder)),recursive = TRUE)
  
  readr::write_rds(x = subset_external_eQTL_split,file = paste0(outDir,"primary_eQTLs/",basename(external_eQTL_folder),"/",comp,"_external_eQTL_subsets.rds"))
  readr::write_rds(x = subset_eqtl_split,file = paste0(outDir,"primary_eQTLs/",basename(external_eQTL_folder),"/",comp,"_eQTL_subsets.rds"))
rm(subset_external_eQTL_split,subset_eqtl_split)
gc()
  }
