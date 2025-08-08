# coloc pQTL vs eQTL

# colocalization analysis
# https://cran.r-project.org/web/packages/coloc/vignettes/a03_enumeration.html
.libPaths(c('/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/',"/software/teamtrynka/ma23/R4.1/libs",
            "/software/R-4.1.0/lib/R/library"))
library(tidyverse)
library(tidyr)
library(coloc)
library(stringr)
library(purrr)

options(scipen = 999) # prevents from showing very small numbers as 0

outDir = "../../data/results/8.colocalisation_analysis"
dir.create(outDir)

args = commandArgs(trailingOnly=TRUE)

# Checking arguments
if (length(args)<4) {
  stop("5 arguments must be supplied", call. = FALSE)
}

GWAS_path = as.character(args[1]) # file path of the harmonised summary stats from selected GWAS trait/study
GWAS_lead_variants_path = as.character(args[2]) # file path of the lead variants from selected GWAS trait/study
eQTL_path = as.character(args[3]) # file path of the eQTL results of interest (created in 4.Inspect_eQTL_results.R)
eQTL_nominal_path=as.character(args[4]) # the nominal eQTL results for same conditions (see Snakefile)
distance_threshold=as.numeric(args[5]) # Define the distance threshold around lead variant in bp (e.g., 500,000 bp)
# test data
GWAS_path = as.character("../../../resources/harmonised_GWAS/GCST90012877_Schwartzentruber_2021/33589840-GCST90012877-EFO_0000249.h.tsv.gz") 
GWAS_lead_variants_path = as.character("../../../resources/Jeremy_medrXiv_AD_loci_GRCh38.txt")  # re-do analysis with the real file
eQTL_path = as.character("../../data/results/4.Inspect_eQTL_results/tensorQTL_variant_gene_60PCs.csv")
eQTL_nominal_path=as.character("../../data/results/tensorqtl/60/")
distance_threshold = 500000

# GWAS data
gwas = read_tsv(GWAS_path) %>%
  dplyr::select(hm_variant_id, hm_beta,p_value,hm_effect_allele, hm_chrom,hm_pos,standard_error) %>%
  dplyr::rename(variant_id = hm_variant_id,
                beta = hm_beta,
                pval_nominal = p_value,
                effect_allele = hm_effect_allele) %>%
  dplyr::mutate(chr = as.numeric(hm_chrom),
                pos=as.numeric(hm_pos), # ensuring types are correct
                varbeta = standard_error^2) %>%
  dplyr::select(chr, pos,variant_id,effect_allele, beta,varbeta,pval_nominal)

message("Printing colnames of GWAS file:")
print( colnames(gwas))

### AD loci have been converted from GRCh37 to GRCh38 in 6.1.LiftOver_Jeremy_AD_GWAS.R 
gwas_lead = read_tsv(GWAS_lead_variants_path)

# check all lead variants are in the summary file
if(sum(!gwas_lead$variant_id %in% gwas$variant_id)==0){
  message("All variants from lead GWAS file are present in the GWAS summary file")
} else{
  message("WARNING: Some variants from lead GWAS file are not present in the GWAS summary file")
  
}

## eQTL data
eqtl_full = read_csv(eQTL_path)

# run per condition
coloc_res=list()
for(comp in unique(eqtl_full$group)){
  
  eqtl = eqtl_full %>%
    dplyr::rename(beta = slope, # the slope is the beta
                  gene_id = phenotype_id
    ) %>%
    tidyr::separate_wider_delim(cols =  position, names = c("chr", "pos"), delim = "_") %>%
    dplyr::mutate(varbeta = slope_se^2,
                  chr = as.numeric(chr),
                  pos=as.numeric(pos)) %>%
    dplyr::filter(group == comp) %>%
    dplyr::select(chr, pos,gene_id, symbol, variant_id ,maf, beta,varbeta,pval_nominal, qval) %>%
    dplyr::filter(complete.cases(.)) %>% # will remove rows with NA, which are those from the X chr
    #group_by(variant_id) %>% # select lead associations per variant
    #dplyr::filter(pval_nominal == min(pval_nominal)) %>%
    #ungroup() %>%
    dplyr::arrange(chr,pos) 
  
  message("Printing colnames of eqtl file:")
  print( colnames(eqtl))
  
  # nominal results
  
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
    dplyr::select(chr, pos,minor_allele,major_allele, gene_id, variant_id ,maf, beta,varbeta,pval_nominal) %>%
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
  
  
  
  
  # For every eQTL lead variant close enough to GWAS locus, subset full eqtl summary based on positions that 
  # are within the distance threshold to GWAS lead
  subset_eqtl = list()
  for(locus in gwas_lead$Locus_name){
    
    subset_eqtl[[locus]] = eqtl %>%
      group_by(symbol) %>% # group by eQTL locus and extract the lead variants that reach significance
      dplyr::filter(qval == min(qval)) %>% # doesn't do any filtering because cis_map already gives the lead eQTL variants per gene
      dplyr::filter(qval < 0.05) %>% # filter for significance
      dplyr::filter(chr == gwas_lead[gwas_lead$Locus_name == locus,]$Chr & abs(pos -  gwas_lead[gwas_lead$Locus_name == locus,]$SNP_pos) <= distance_threshold) %>%
      ungroup()
    # now subset using variants from the nominal tests, per gene that remains significant from the previous step
    # also at specified distance from gwas variant
    
    sign_egenes = subset_eqtl[[locus]]$gene_id
    ####
    
    subset_eqtl[[locus]] = subset_eqtl[[locus]] %>%
      right_join(eqtl_nominal) %>%
      dplyr::filter(gene_id %in% sign_egenes)
    
    subset_eqtl[[locus]] = subset_eqtl[[locus]] %>%
      dplyr::filter(subset_eqtl[[locus]]$chr == gwas_lead[gwas_lead$Locus_name == locus,]$Chr & abs(pos -  gwas_lead[gwas_lead$Locus_name == locus,]$SNP_pos) <= distance_threshold) 
    
    # add qval info to the lead eqtl variant
    # also symbol
    
  }
  
  # subset GWAS summary at same location, and for shared variants with eqtl
  
  subset_gwas = list()
  for(locus in gwas_lead$Locus_name){
    
    subset_gwas[[locus]] = gwas %>%
      dplyr::filter(gwas$chr == gwas_lead[gwas_lead$Locus_name == locus,]$Chr & abs(pos -  gwas_lead[gwas_lead$Locus_name == locus,]$SNP_pos) <= distance_threshold) 
    
    message("There are ",sum(subset_gwas[[locus]]$variant_id %in% unique(subset_eqtl[[locus]]$variant_id)),
            " shared variants between GWAS and eQTL at GWAS locus ", locus)
    
    #########################################
    ##### TO-DO: check for swapped REF/ALT alleles ######
    ##########################################
    ##########!!!!!!!!!!!!!#################
    
    
    subset_gwas[[locus]]= subset_gwas[[locus]] %>%
      dplyr::filter(subset_gwas[[locus]]$variant_id %in% subset_eqtl[[locus]]$variant_id) 
    
    # subset eqtl
    message("Subsetting eQTL variants")
    subset_eqtl[[locus]]= subset_eqtl[[locus]] %>%
      dplyr::filter(subset_eqtl[[locus]]$variant_id %in% subset_gwas[[locus]]$variant_id) 
    
  }
  
  
  ### format correctly for coloc and test
  
  coloc_res[[comp]] = list()
  for(locus in gwas_lead$Locus_name){
    
    message("Building coloc input around lead GWAS variant at locus ", locus)
    
    
    # eQTL
    # sdY= 1 because the expression was standardized to have a variance of 1 
    # (was scaled and centered in 2.Filter_seurat_for_eQTL.R)
    
    #test per significant eQTL gene
    coloc_res[[comp]] [[locus]] = list()
    for(gene in unique(subset_eqtl[[locus]]$gene_id)){
      subset_eqt_at_gene =  subset_eqtl[[locus]] %>%
        dplyr::filter(gene_id == gene)
      
      
      subset_eqtl_coloc = list(beta = subset_eqt_at_gene$beta,
                               varbeta = subset_eqt_at_gene$varbeta,
                               snp = subset_eqt_at_gene$variant_id,
                               position=subset_eqt_at_gene$pos,
                               type="quant",
                               sdY= 1 
      )
      if(is.null(check_dataset(subset_eqtl_coloc)) & length(subset_eqt_at_gene$variant_id)>50){
        message("eQTL list for coloc is valid")
      } else{
        message("WARNING: eQTL list for coloc is NOT valid")
        
      }
      
      
      # Convert each column to a named list element
      subset_gwas_variants = subset_gwas[[locus]] %>%
        dplyr::filter(variant_id %in% subset_eqt_at_gene$variant_id)
      
      if(nrow(subset_gwas_variants) == nrow(subset_eqt_at_gene)){
        message("eQTL and GWAS subsets have the same number of variants")
      } else{
        message("WARNING: eQTL and GWAS subsets DO NOT have the same number of variants")
        
      }
      
      subset_gwas_coloc = list(beta = subset_gwas_variants$beta,
                               varbeta = subset_gwas_variants$varbeta,
                               snp = subset_gwas_variants$variant_id,
                               position=subset_gwas_variants$pos,
                               type="cc"
      )
      
      if(is.null(check_dataset(subset_gwas_coloc)) & length(subset_gwas[[locus]]$variant_id)>50){
        message("GWAS list for coloc is valid")
      } else{
        message("WARNING: GWAS list for coloc is NOT valid")
        
      }
      
      message("Running coloc")
      
      if(is.null(check_dataset(subset_gwas_coloc)) & length(subset_gwas[[locus]]$variant_id)>50 & is.null(check_dataset(subset_eqtl_coloc)) & length(subset_eqtl[[locus]]$variant_id)>50 ){
        coloc_res[[comp]] [[locus]][[gene]]= coloc.abf(dataset1=subset_gwas_coloc,
                                                       dataset2=subset_eqtl_coloc,
                                                       p12=1e-6)
      } else{
        coloc_res[[comp]] [[locus]][[gene]]= NULL
      }
      # gives posterior probability (Bayesian: not p-vals!!) for each of the following hypotheses:
      # H0: neither trait has a genetic association in the region
      # H1: only trait 1 has a genetic association in the region
      # H2: only trait 2 has a genetic association in the region
      # H3: both traits are associated, but with different causal variants
      # H4: both traits are associated and share a single causal variant
    }
  }
}

# A sensitivity analysis can be used, post-hoc, to determine the range of 
# prior probabilities for which a conclusion is still supported

# 
# locus="BIN1" 
# sensitivity(coloc_res[[comp]] [["BIN1"]]$ENSG00000136717,rule="H4 > 0.5") # pass 99%
# sensitivity(coloc_res[[comp]] [["BIN1"]]$ENSG00000236682,rule="H4 > 0.5") # fail
# 
# coloc_res[[comp]] [["BIN1"]]$ENSG00000136717$summary[[6]]
# 
# # check why many are not centered around gwas peak
# 
# sensitivity(coloc_res[[comp]] [["CCDC6"]]$ENSG00000108091,rule="H4 > 0.5") 
# ENSG00000108091 == CCDC6

# summarise results
coloc_res_summary = list()
for(comp in names(coloc_res)){
  for(nam in names(coloc_res[[comp]] )){
    if(length(coloc_res[[comp]][[nam]])==0){
      next
      
    }else{
      for(gene in names(coloc_res[[comp]][[nam]])){
        coloc_res_summary[[comp]][[nam]][[gene]] = data.frame(GWAS_Locus = nam,
                                                              eGene = gene, 
                                                              eQTL_comparison=comp,
                                                              PP_H4 = coloc_res[[comp]][[nam]][[gene]]$summary[[6]] # extract PP H4
        )
        
      }
      coloc_res_summary[[comp]][[nam]] = do.call("rbind",coloc_res_summary[[comp]][[nam]])
      
    }
    
  }
  
  if(!is.null(coloc_res_summary[[comp]])){
    coloc_res_summary[[comp]] = do.call("rbind",coloc_res_summary[[comp]])
    
  }
  
  
}
#Convert the list to a data frame
coloc_res_summary = do.call("rbind",coloc_res_summary)

sensitivity(coloc_res$`60_Not_proliferating_untreated`$BIN1$ENSG00000136717,"H4 > 0.75") # 0.99
sensitivity(coloc_res$`60_Not_proliferating_LPS`$PILRA$ENSG00000121716,"H4 > 0.75") # 0.99
sensitivity(coloc_res$`60_Not_proliferating_untreated`$PILRA$ENSG00000085514,"H4 > 0.75") # 0.92
sensitivity(coloc_res$`60_Not_proliferating_LPS`$ADAMTS4$ENSG00000158869,"H4 > 0.75") # 0.83
sensitivity(coloc_res$`60_Not_proliferating_untreated`$CCDC6$ENSG00000108091,"H4 > 0.75") #0.75

coloc_res_summary %>%
  dplyr::arrange(desc(PP_H4)) %>%
  write_tsv(.,paste0(outDir,"/coloc_summary_single_causal_variant_",distance_threshold/1000,"kb.txt"))

write_rds(coloc_res,paste0(outDir,"/coloc_results_single_causal_variant_",distance_threshold/1000,"kb.rds"))

### add lead variants info,symbol, variants that colocalize
# make prettier plots



### more than one causal variant ###
# coloc.abf makes the simplifying assumption that each trait has at most one causal
# variant in the region under consideration. This means it does not need to know 
# about correlation (LD) between SNPs, because no model under consideration has 
# two or more SNPs in it. coloc.susie extends the functionality by allowing more 
# than one causal variant, but this means it needs to know the LD between the SNPs. 
# You may be able to get this from the study data, but if not, reference data 
# such as 1000 Genomes can be used.

# LD matrix already. This will be a square numeric matrix of dimension equal to 
# the number of SNPs, with dimnames corresponding to the SNP ids
str(D1$LD)

# plot colocs
