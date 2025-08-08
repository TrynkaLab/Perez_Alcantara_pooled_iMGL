# run colocalisation analysis with coloc
# https://cran.r-project.org/web/packages/coloc/vignettes/a03_enumeration.html

library(tidyverse)
library(tidyr)
library(coloc)
library(stringr)
library(purrr)

options(scipen = 999) # prevents from showing very small numbers as 0

# 
# args = commandArgs(trailingOnly=TRUE)
# 
# # Checking arguments
# if (length(args)<4) {
#   stop("4 arguments must be supplied", call. = FALSE)
# }
# 


#### fix to run on pipeline
# external_data_folder = as.character(args[1]) # folder that contains the harmonised summary stats from selected external_data trait/study
# GWAS_path = as.character(args[2]) # file path of the GWAS results of interest (created in 4.Inspect_GWAS_results.R)
# GWAS_nominal_path=as.character(args[3]) # the nominal GWAS results for same conditions (see Snakefile)
# distance_threshold=as.numeric(args[4]) # Define the distance threshold around lead variant in bp (e.g., 500,000 bp)


######
outDir = "../../../data/results/phagocytosis/2.check_association_results/lm_1pct_filtered_deflated/coloc_results/"
dir.create(outDir)


# test data
external_data_folder = as.character("../../../data/results/phagocytosis/2.check_association_results/lm_1pct_filtered_deflated/subsets_for_coloc/GCST90027158/")
distance_threshold = 500000

external_data_id = basename(external_data_folder)
if(str_detect(external_data_folder,"primary_GWASs")){
  external_data_files = list.files(external_data_folder,pattern = "*external_GWAS_subsets.rds")
  
} else {
  external_data_files = list.files(external_data_folder,pattern = "*external_GWAS_subsets.rds")
  
}

external_data = list()
GWAS = list()

  for(comp in gsub(x = external_data_files,pattern = "_external_GWAS_subsets.rds",replacement = "")){
    external_data[[comp]] = readr::read_rds(file = paste0(external_data_folder,"/",comp,"_external_GWAS_subsets.rds"))
    GWAS[[comp]] = readr::read_rds(file = paste0(external_data_folder,"/",comp,"_my_GWAS_subsets.rds"))
    
  }


### format correctly for coloc and test
coloc_res= list()

for(comp in names(external_data)){
  coloc_res[[comp]] = list()
  
  for(locus in names(external_data[[comp]])){
    coloc_res[[comp]][[locus]] = list()
      
      message("Building coloc input around lead external_data variant at locus ", locus)

      
      # GWAS
      # sdY= 1 because the expression was standardized to have a variance of 1 
      # (was scaled and centered in 2.Filter_seurat_for_GWAS.R)
      
      #### remove duplicated SNPs
      GWAS[[comp]][[locus]] = GWAS[[comp]][[locus]] %>%
        distinct(variant_id,.keep_all = TRUE)
      
      #test per significant GWAS gene
      subset_GWAS_coloc = list(beta = GWAS[[comp]][[locus]]$beta,
                               varbeta = GWAS[[comp]][[locus]]$varbeta,
                               snp = GWAS[[comp]][[locus]]$variant_id,
                               position=GWAS[[comp]][[locus]]$pos,
                               type="quant",
                               sdY= 1 
      )
      if(is.null(coloc::check_dataset(subset_GWAS_coloc)) & length(GWAS[[comp]][[locus]]$variant_id)>50){
        message("GWAS list for coloc is valid")
      } else{
        message("WARNING: GWAS list for coloc is NOT valid")
        
      }
      
      
      
      if(nrow( external_data[[comp]][[locus]]) == nrow(GWAS[[comp]][[locus]])){
        message("GWAS and external_data subsets have the same number of variants")
      } else{
        message("WARNING: GWAS and external_data subsets DO NOT have the same number of variants")
        message("WARNING: Checking and removing duplicated and extra variants")
        external_data[[comp]][[locus]] = external_data[[comp]][[locus]] %>%
          dplyr::filter(!duplicated(variant_id))
        GWAS[[comp]][[locus]] = GWAS[[comp]][[locus]] %>%
          dplyr::filter(!duplicated(variant_id))
        
        
      }
      

        # double-check all the GWAS I use are like this (CC)
        # my phagocytosis GWAS wouldn't be, for example
        subset_external_data_coloc = list(beta = external_data[[comp]][[locus]]$beta,
                                          varbeta = external_data[[comp]][[locus]]$varbeta,
                                          snp = external_data[[comp]][[locus]]$variant_id,
                                          position=external_data[[comp]][[locus]]$pos,
                                          type="cc"
        )
      
      
      if(is.null(check_dataset(subset_external_data_coloc)) & length(external_data[[comp]][[locus]]$variant_id)>50){
        message("external_data list for coloc is valid")
      } else{
        message("WARNING: external_data list for coloc is NOT valid")
        
      }
      
      message("Running coloc")
      
      if(is.null(coloc::check_dataset(subset_external_data_coloc)) & length(external_data[[comp]][[locus]]$variant_id)>50 & is.null(coloc::check_dataset(subset_GWAS_coloc)) & length(GWAS[[comp]][[locus]]$variant_id)>50 ){
        coloc_res[[comp]][[locus]]= coloc.abf(dataset1=subset_external_data_coloc,
                                                      dataset2=subset_GWAS_coloc,
                                                      p12=1e-5)
      } else{
        coloc_res[[comp]][[locus]]= NULL
      }
      # gives posterior probability (Bayesian: not p-vals!!) for each of the following hypotheses:
      # H0: neither trait has a genetic association in the region
      # H1: only trait 1 has a genetic association in the region
      # H2: only trait 2 has a genetic association in the region
      # H3: both traits are associated, but with different causal variants
      # H4: both traits are associated and share a single causal variant
    
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
# # check why many are not centered around external_data peak
# 
# sensitivity(coloc_res[[comp]] [["CCDC6"]]$ENSG00000108091,rule="H4 > 0.5") 
# ENSG00000108091 == CCDC6

# plot_dataset(GWAS[["35_LPS_Not_proliferating"]][["6_41161469_TREM2"]][["TREM2"]])
# sensitivity(coloc_res[[comp]][[locus]],rule = "H4>0.7" )

# summarise results
coloc_res_summary = list()
for(comp in names(coloc_res)){
  coloc_res_summary[[comp]] = list()
  for(nam in names(coloc_res[[comp]] )){
    if(length(coloc_res[[comp]][[nam]])==0){
      next
      
    }else{
        coloc_res_summary[[comp]][[nam]] = data.frame(external_data_Locus = nam,
                                                              GWAS_comparison=comp,
                                                              PP_H4 = coloc_res[[comp]][[nam]]$summary[[6]] # extract PP H4
        )
        
      }
      
    }
  coloc_res_summary[[comp]] = do.call("rbind",coloc_res_summary[[comp]])
  
  }
  

#Convert the list to a data frame
coloc_res_summary = do.call("rbind",coloc_res_summary)

# sensitivity(coloc_res$`60_Not_proliferating_untreated`$BIN1$ENSG00000136717,"H4 > 0.75") # 0.99
# sensitivity(coloc_res$`60_Not_proliferating_LPS`$PILRA$ENSG00000121716,"H4 > 0.75") # 0.99
# sensitivity(coloc_res$`60_Not_proliferating_untreated`$PILRA$ENSG00000085514,"H4 > 0.75") # 0.92
# sensitivity(coloc_res$`60_Not_proliferating_LPS`$ADAMTS4$ENSG00000158869,"H4 > 0.75") # 0.83
# sensitivity(coloc_res$`60_Not_proliferating_untreated`$CCDC6$ENSG00000108091,"H4 > 0.75") #0.75

dir.create(paste0(outDir,"/",external_data_id))

coloc_res_summary %>%
  dplyr::arrange(desc(PP_H4)) %>%
  readr::write_tsv(.,paste0(outDir,"/",external_data_id,"/coloc_summary_single_causal_variant_",distance_threshold/1000,"kb.txt"))

readr::write_rds(coloc_res,paste0(outDir,"/",external_data_id,"/coloc_results_single_causal_variant_",distance_threshold/1000,"kb.rds"))


### add lead variants info,symbol, variants that colocalize
# make prettier plots
#sensitivity(coloc_res$`60_Not_proliferating_LPS`$`7_100334426_ZCWPW1/NYAP1`$PILRB,"H4 > 0.75") # 0.99


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
#str(D1$LD)

# plot colocs
