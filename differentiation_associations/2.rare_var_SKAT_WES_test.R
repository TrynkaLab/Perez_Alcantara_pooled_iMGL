# SKAT analysis
# testing regressing pool, type and treatment where found
# to fit pool replicates

library(vcfR)
library(tidyverse)
library(SKAT)
## input arguments
args = commandArgs(trailingOnly=TRUE)

if (length(args)<7) {
  stop("You need to detail 6 input and 1 output file paths.n",
       call. = FALSE)
} else if (length(args) == 7) {
  output_directory = args[1]
  WES_missense_del_gt_path = args[2]
  WES_missense_del_cleaninfo_path = args[3]
  line_prop_changes_path=args[4]
  genotype_pc_path=args[5]
  line_info_path=args[6]
  n_gene=as.numeric(args[7])
}

### test 
output_directory = "../../data/results/2.rare_var_SKAT_WES/with_correction/"
WES_missense_del_gt_path ="../../../hipsci_genotype_processing/data/WES/VEP/missense_clean_GRCh38.vcf.gz"
WES_missense_del_cleaninfo_path ="../../../hipsci_genotype_processing/data/WES/VEP/missense_clean_GRCh38.txt"
line_prop_changes_path = "../../data/results/1.2.scale_proliferation/line_prop_changes_microglia_premac.csv"
genotype_pc_path = "../../../OTAR2065_sc_eQTL/data/genotype/plink_genotypes/all_pools.no_outliers.genotype.MAF05.eigenvec"
line_info_path="../../../OTAR2065_differentiation_efficiency/data/donor_metadata_complete_with_imputed_sex.csv"
n_gene=178
######

comparison = str_remove(str_split_i(line_prop_changes_path,pattern = "/",i=6),".csv")
dir.create(paste0(output_directory,comparison), recursive = T)

# read in data ################
# general metadata
line_info = read_csv(line_info_path) %>%
  dplyr::select(donor,Sex)
# PCs
geno_PCs = read_tsv(genotype_pc_path) %>%
  dplyr::select(-`#FID`) %>%
  dplyr::rename(
    line = IID,
    genotypePC1 = PC1,
    genotypePC2 = PC2,
    genotypePC3 = PC3,
    genotypePC4 = PC4,
    genotypePC5 = PC5
  ) %>%
  dplyr::mutate(donor = if_else(str_detect(line,"-"),
                                str_split_i(line,pattern="-",1),
                                str_split_i(line,pattern="_",1))) %>%
  dplyr:::select(donor,line:genotypePC5)

# add to line info

line_info = line_info %>%
  dplyr::left_join(geno_PCs)

# line proportion changes 
line_prop_changes = read_csv(line_prop_changes_path) %>%
  dplyr::mutate(donor = if_else(str_detect(line,"-"),
                                str_split_i(line,pattern="-",1),
                                str_split_i(line,pattern="_",1))) %>%
  dplyr::filter(!scaled_log_fraction %in% c(NA,Inf,-Inf)) %>%
  dplyr::select(any_of(c("scaled_log_fraction","line","donor","differentiation","treatment","pool","prop_unadjusted_min_value","type"))) %>%
  distinct()

message("There are ", length(unique(line_prop_changes$line)), " lines in the analysis after Na and Inf filters") # 238 for premac vs iPSC
# 241 for microglia vs premac

## agregate all phenotype and covariate information

full_pheno_info = line_prop_changes %>%
  dplyr::left_join(line_info) %>%
  # omit donors that are outliers (no genotype PC info) and sh5y5y
  dplyr::filter(!is.na(genotypePC1))
colnames(full_pheno_info)

# clean VEP info ###
WES_missense_del_info = read_tsv(WES_missense_del_cleaninfo_path) %>%
  dplyr::mutate(CADD_PHRED = case_when(CADD_PHRED == "." ~ NA,
                                       .default = as.double(CADD_PHRED))) %>%
  # defining as deleterious - union of LoF and missense pathogenic mutations
  # LoF -  frameshift, stop-gain, transcript_ablation, splice acceptor or splice donor variants
  # see https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
  # missense pathogenic - annotated as missense or start loss with a CADD Phred score cutoff > 15 (as Pau did)
  dplyr::mutate(category = case_when((str_detect(Consequence,pattern="missense") | 
                                        str_detect(Consequence,pattern="start_lost") | 
                                        str_detect(Consequence,pattern="protein_altering")) & CADD_PHRED >=15 ~ "deleterious",
                                     (str_detect(Consequence,pattern="frameshift")|
                                        str_detect(Consequence,pattern="stop_gained")|
                                        str_detect(Consequence,pattern="transcript_ablation")|
                                        str_detect(Consequence,pattern="splice_acceptor")|
                                        str_detect(Consequence,pattern="splice_donor")
                                     ) ~ "deleterious",
                                     .default = "missense_non_deleterious"),
                ID = paste(CHROM,POS,REF,ALT,sep = "_")) %>%
  dplyr::filter(!duplicated(ID))

length(unique(WES_missense_del_info$Gene)) # 18,588 genes, 18,581 gene IDs (SYMBOL)
table(WES_missense_del_info$category) # 69,223 deleterious, 62,022 missense non-deleterious variants

# clean VCF ###
WES_missense_del_genotype = vcfR::read.vcfR(WES_missense_del_gt_path)

# INFO section
info = WES_missense_del_genotype@fix[,1:5] %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(ID = paste(CHROM,POS,REF,ALT,sep = "_"))
# GT section
WES_missense_del_genotype = vcfR::extract_gt_tidy(WES_missense_del_genotype)

WES_missense_del_genotype = WES_missense_del_genotype %>%
  dplyr::mutate(genotype = case_when(gt_GT == "1/1" ~ 2,
                                     gt_GT %in% c("0/1","1/0") ~ 1,
                                     gt_GT == "0/0" ~ 0,
                                     .default = NA),
                line = if_else(str_detect(Indiv,"HPSI"), 
                               str_split_i(Indiv,"-",i=2),
                               NA)) %>%
  dplyr::mutate(genotype = as.numeric(genotype),
                donor = str_split_i(line,"_",i=1)) 



table(info$ID %in% WES_missense_del_info$ID) # lost a handful after filters
table(WES_missense_del_info$ID %in% info$ID) 

# join genotype info
WES_missense_del_genotype = cbind(WES_missense_del_genotype,info)

### subset genes based on n_genes
length(unique(WES_missense_del_info$Gene)) # 18,588
gene_id = unique(WES_missense_del_info$Gene)[n_gene]

###### SKAT analysis
# needs:
# Z: genotype matrix with samples as rows and SNPs as columns
# y.c: vector of phenotype of interest of length = number of samples
# X: covariate matrix, with samples as rows and first column as all 0s (intercept)
# with subsequent columns as covariates to adjust for (eg. sex, genotype PCs)

for(mut_type in c("missense_non_deleterious","deleterious")){
  message("Working on ", mut_type, " variants.")
  message("\n")
  
  # subset right variants
  vars_within_mut_type = WES_missense_del_info %>%
    dplyr::filter(category == mut_type & Gene == gene_id) %>%
    dplyr::select(ID,Gene,SYMBOL)
  
  
  if("treatment" %in% colnames(full_pheno_info)){
    # for(treat in c("untreated","LPS","IFN")){
    #   message("Working on ", treat, " treatment")
    #   message("\n")
    # correct fr treatment
      

    message("subset right metadata") 
    
    meta_subset = full_pheno_info %>%
      # filter objects to shared lines
      dplyr::filter((line %in% WES_missense_del_genotype$line)) %>%
      dplyr::mutate(sex_numeric = case_when(Sex == "Male" ~ 1, # female as reference (0)
                                            .default = 0)) 
    
    # pre-adjust for pool, treatment and differentiation
    
    # lmm_result = lmerTest::lmer(scaled_log_fraction ~ pool + treatment +  (pool | treatment) + (pool | differentiation),
    #                            data = meta_subset)
    # residuals_lmm = resid(lmm_result)
    lm_result = lm(scaled_log_fraction ~ pool + treatment + differentiation,
                   data = meta_subset)
    meta_subset$residuals = lm_result$residuals
    
    # The residuals after fitting a mixed model and a regular linear model are extremely similar
    
    # so choosing a linear model because it's much faster
    # plot(residuals_lmm,lm_result$residuals)
    message("creating vector of phenotype") 
    
    y.c = meta_subset$residuals
    
      ## creating covariate matrix ###
      X = as.matrix(data.frame(#pool=meta_subset$pool, # I can't include pool, type or treatment because every element of matrix needs to be same type
        # would need to pre-adjust somehow (pre-adjust then fit residuals)
        sex=meta_subset$sex_numeric,
        prop_unadjusted_min_value= meta_subset$prop_unadjusted_min_value,
        genotypePC1 = meta_subset$genotypePC1,  genotypePC2 = meta_subset$genotypePC2))
      
      skat_res = list()
      
      #### important check ##
      # drop genes that don't amount to at 3 lines across all variants in at least 2 categories 
      # (often heterozygote or ALT homozygote don't have them bc deleterious vars are rare)
      WES_missense_del_genotype_subset = WES_missense_del_genotype %>%
        # filter objects to shared lines
        dplyr::filter((ID %in% vars_within_mut_type$ID) & (line %in% meta_subset$line))
      
      check_min_lines = sum(colSums(table(WES_missense_del_genotype_subset$ID,WES_missense_del_genotype_subset$genotype))>3)>2
      # if the tibble has < 2 rows (0 or just 1 variant for that gene and mut status), save NA object
      # if it fails the minimum line number across variants, save NA object
      if(dim(vars_within_mut_type)[1]<2 | check_min_lines == FALSE  ){
        # not saving empty gene results so I don't generate thousands of files
        # if(mut_type == "missense_non_deleterious"){
        #   message("The ", mut_type, " tibble has < 2 rows (0 or just 1 variant for that gene and mut status. Not saving file ... ")
        #   
        #   next()
        # }else{
        #   message("The ", mut_type, " tibble has < 2 rows (0 or just 1 variant for that gene and mut status. Not saving file ... ")
        #   
        #   break()
        #   
        # }
        skat_res[[gene_id]]$p.value = NA
        skat_res[[gene_id]]$resampling_pval = NA
        skat_res[[gene_id]]$gene_name = unique(WES_missense_del_info$SYMBOL)[n_gene]
        skat_res[[gene_id]]$gene_id = unique(WES_missense_del_info$Gene)[n_gene]
        
      }else{
        message("creating genotype matrix") 
        
        Z = WES_missense_del_genotype_subset %>%
          dplyr::select(ID,line,genotype) %>%
          tidyr::pivot_wider(id_cols = line,names_from = ID,values_from = genotype) %>%
          dplyr::arrange(match(line, meta_subset$line)) 
        # need to expand to all line repeats
        # doing it by joining with meta_data matrix to ensure right order
        
        Z = meta_subset %>%
          dplyr::select(line,donor) %>%
          dplyr::left_join(Z) %>%
          dplyr::select(-donor)
        
        if(identical(Z$line,meta_subset$line)){
          message("Pool repeat genotype information and metadata are in correct order")
        }else{
          message("ERROR: Pool repeat genotype information and metadata are in incorrect order!!")
          break()
          
        }
        
        # subset to needed columns
        Z = as.matrix(Z[,2:ncol(Z)])
        
        chr = str_split_i(vars_within_mut_type$ID[1],pattern = "_",i=1)
        
        message("running SKAT") 
        message("\n") 
        
        ####### running SKAT
        # see https://cran.r-project.org/web/packages/SKAT/vignettes/SKAT.pdf
        if(chr !="chrX"){
          ########### running null model with permutations to control for FWER
          obj=SKAT_Null_Model(y.c ~ X, out_type="C",n.Resampling=1000, type.Resampling="bootstrap") 
          skat_res[[gene_id]] = SKAT(Z, obj, method="SKATO") # including pool reps, ENSG00000104133 SPG11 has pval of 8.052946e-24
          # 8.760352e-06 averaging across pools and reps
          skat_res[[gene_id]]$gene_name = unique(vars_within_mut_type$SYMBOL)
          skat_res[[gene_id]]$gene_id = unique(vars_within_mut_type$Gene)
          skat_res[[gene_id]]$resampling_pval = Get_Resampling_Pvalue(skat_res[[gene_id]])$p.value # ENSG00000104133 resampling pval = 9e-4
          skat_res[[gene_id]]$genotype = WES_missense_del_genotype_subset %>%
            dplyr::select(ID,line,genotype) %>%
            tidyr::pivot_wider(id_cols = line,names_from = ID,values_from = genotype) %>%
            dplyr::arrange(match(line, meta_subset$line)) 
          skat_res[[gene_id]]$metadata = meta_subset
          
          # resampling p-value is the same even if p-value is smaller
          
        } else {
          sex = X[,"sex"] 
          sex = try(case_match(sex, 0 ~ 2, 1 ~ 1, .default = sex)) #  They should be either 1 (=male) or 2(=female)!
          prop_unadjusted_min_value = X[,"prop_unadjusted_min_value"]
          genotypePC1 = X[,"genotypePC1"]
          genotypePC2 = X[,"genotypePC2"]
          obj=SKAT_Null_Model_ChrX(y.c ~ sex + prop_unadjusted_min_value + genotypePC1 + genotypePC2, 
                                   out_type="C",
                                   n.Resampling=1000, type.Resampling="bootstrap",
                                   SexVar = "sex")
          skat_res[[gene_id]] = SKAT_ChrX(Z, obj, method="SKATO") 
          skat_res[[gene_id]]$gene_name = unique(vars_within_mut_type$SYMBOL)
          skat_res[[gene_id]]$gene_id = unique(vars_within_mut_type$Gene)
          skat_res[[gene_id]]$resampling_pval = Get_Resampling_Pvalue(skat_res[[gene_id]])$p.value
          skat_res[[gene_id]]$genotype = WES_missense_del_genotype_subset %>%
            dplyr::select(ID,line,genotype) %>%
            tidyr::pivot_wider(id_cols = line,names_from = ID,values_from = genotype) %>%
            dplyr::arrange(match(line, meta_subset$line)) 
          skat_res[[gene_id]]$metadata = meta_subset
          
        }
      }
    
  }else{
      
      message("subset right metadata") 
      
      meta_subset = full_pheno_info %>%
        # filter objects to shared lines
        dplyr::filter((line %in% WES_missense_del_genotype$line)) %>%
        dplyr::mutate(sex_numeric = case_when(Sex == "Male" ~ 1, # female as reference (0)
                                              .default = 0)) 
      
      # pre-adjust for pool and differentiation
      
      # lmm_result = lmerTest::lmer(scaled_log_fraction ~ pool + (pool | differentiation),
      #                            data = meta_subset)
      # residuals_lmm = resid(lmm_result)
      # The residuals after fitting a mixed model and a regular linear model are identical
      # because differentiation and pool are heavily confounded (for all pools except pool 2)
      # so choosing a linear model because it's much faster
      lm_result = lm(scaled_log_fraction ~ pool + differentiation,
                                 data = meta_subset)
      meta_subset$residuals = lm_result$residuals
      
      message("creating vector of phenotype") 
      
      y.c = meta_subset$residuals
      ## creating covariate matrix ###
      X = as.matrix(data.frame(#pool=meta_subset$pool, # I can't include pool because every element of matrix needs to be same type
        # would need to pre-adjust, maybe by taking the residuals after adjusting for pool as random effect?
        sex=meta_subset$sex_numeric,
        prop_unadjusted_min_value = meta_subset$prop_unadjusted_min_value,
        genotypePC1 = meta_subset$genotypePC1,  genotypePC2 = meta_subset$genotypePC2))
      
      skat_res = list()
      
      #### important check ##
      # drop genes that don't amount to at 3 lines across all variants in at least 2 categories 
      # (often heterozygote or ALT homozygote don't have them bc deleterious vars are rare)
      WES_missense_del_genotype_subset = WES_missense_del_genotype %>%
        # filter objects to shared lines
        dplyr::filter((ID %in% vars_within_mut_type$ID) & (line %in% meta_subset$line))
      
      check_min_lines = sum(colSums(table(WES_missense_del_genotype_subset$ID,WES_missense_del_genotype_subset$genotype))>3)>2
      # if the tibble has < 2 rows (0 or just 1 variant for that gene and mut status), save NA object
      # if it fails the minimum line number across variants, save NA object
      if(dim(vars_within_mut_type)[1]<2 | check_min_lines == FALSE  ){
        # not saving empty gene results so I don't generate thousands of files
        # if(mut_type == "missense_non_deleterious"){
        #   message("The ", mut_type, " tibble has < 2 rows (0 or just 1 variant for that gene and mut status. Not saving file ... ")
        #   
        #   next()
        # }else{
        #   message("The ", mut_type, " tibble has < 2 rows (0 or just 1 variant for that gene and mut status. Not saving file ... ")
        #   
        #   break()
        #   
        # }
        skat_res[[gene_id]]$p.value = NA
        skat_res[[gene_id]]$resampling_pval = NA
        skat_res[[gene_id]]$gene_name = unique(WES_missense_del_info$SYMBOL)[n_gene]
        skat_res[[gene_id]]$gene_id = unique(WES_missense_del_info$Gene)[n_gene]
        
      }else{
        message("creating genotype matrix") 
        
        Z = WES_missense_del_genotype_subset %>%
          dplyr::select(ID,line,genotype) %>%
          tidyr::pivot_wider(id_cols = line,names_from = ID,values_from = genotype) %>%
          dplyr::arrange(match(line, meta_subset$line)) 
        # need to expand to all line repeats
        # doing it by joining with meta_data matrix to ensure right order
        
        Z = meta_subset %>%
          dplyr::select(line,donor) %>%
          dplyr::left_join(Z) %>%
          dplyr::select(-donor)
        
        if(identical(Z$line,meta_subset$line)){
          message("Pool repeat genotype information and metadata are in correct order")
        }else{
          message("ERROR: Pool repeat genotype information and metadata are in incorrect order!!")
          break()
          
        }
               
        
        # subset to needed columns
        Z = as.matrix(Z[,2:ncol(Z)])
        
        chr = str_split_i(vars_within_mut_type$ID[1],pattern = "_",i=1)
        
        message("running SKAT") 
        message("\n") 
        
        ####### running SKAT
        # see https://cran.r-project.org/web/packages/SKAT/vignettes/SKAT.pdf
        if(chr !="chrX"){
          ########### running null model with permutations to control for FWER
          obj=SKAT_Null_Model(y.c ~ X, out_type="C",n.Resampling=1000, type.Resampling="bootstrap") 
          skat_res[[gene_id]] = SKAT(Z, obj, method="SKATO") # including pool reps, ENSG00000104133 SPG11 has pval of 8.052946e-24
          # 8.760352e-06 averaging across pools and reps
          skat_res[[gene_id]]$gene_name = unique(vars_within_mut_type$SYMBOL)
          skat_res[[gene_id]]$gene_id = unique(vars_within_mut_type$Gene)
          skat_res[[gene_id]]$resampling_pval = Get_Resampling_Pvalue(skat_res[[gene_id]])$p.value # ENSG00000104133 resampling pval = 9e-4
          skat_res[[gene_id]]$genotype = WES_missense_del_genotype_subset %>%
            dplyr::select(ID,line,genotype) %>%
            tidyr::pivot_wider(id_cols = line,names_from = ID,values_from = genotype) %>%
            dplyr::arrange(match(line, meta_subset$line)) 
          skat_res[[gene_id]]$metadata = meta_subset
          
        } else {
          sex = X[,"sex"] 
          sex = try(case_match(sex, 0 ~ 2, 1 ~ 1, .default = sex)) #  They should be either 1 (=male) or 2(=female)!
          prop_unadjusted_min_value = X[,"prop_unadjusted_min_value"]
          genotypePC1 = X[,"genotypePC1"]
          genotypePC2 = X[,"genotypePC2"]
          obj=SKAT_Null_Model_ChrX(y.c ~ sex + prop_unadjusted_min_value + genotypePC1 + genotypePC2, 
                                   out_type="C",
                                   n.Resampling=1000, type.Resampling="bootstrap",
                                   SexVar = "sex")
          skat_res[[gene_id]] = SKAT_ChrX(Z, obj, method="SKATO") 
          skat_res[[gene_id]]$gene_name = unique(vars_within_mut_type$SYMBOL)
          skat_res[[gene_id]]$gene_id = unique(vars_within_mut_type$Gene)
          skat_res[[gene_id]]$resampling_pval = Get_Resampling_Pvalue(skat_res[[gene_id]])$p.value
          skat_res[[gene_id]]$genotype = WES_missense_del_genotype_subset %>%
            dplyr::select(ID,line,genotype) %>%
            tidyr::pivot_wider(id_cols = line,names_from = ID,values_from = genotype) %>%
            dplyr::arrange(match(line, meta_subset$line)) 
          skat_res[[gene_id]]$metadata = meta_subset
          
        }
     
      }


  }
  message("Saving file...")
  write_rds(skat_res,paste0(output_directory,comparison,"/",n_gene,"_",mut_type,"_SKATO_result.RDS"))
  rm(skat_res)
  gc()
}

