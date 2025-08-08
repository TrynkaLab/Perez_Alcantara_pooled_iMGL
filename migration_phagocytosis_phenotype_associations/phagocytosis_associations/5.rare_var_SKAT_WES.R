# SKAT analysis

library(vcfR)
library(tidyverse)
library(SKAT)
## input arguments
args = commandArgs(trailingOnly=TRUE)

if (length(args)<8) {
  stop("You need to detail 7 input and 1 output file paths.n",
       call. = FALSE)
} else if (length(args) == 8) {
  output_directory = args[1]
  WES_missense_del_gt_path = args[2]
  WES_missense_del_cleaninfo_path = args[3]
  line_prop_changes_path=args[4]
  genotype_pc_path=args[5]
  line_info_path=args[6]
  n_gene=as.numeric(args[7])
  scaled_proliferation_path=args[8]
}

### test 
# output_directory = "../../../data/results/phagocytosis/5.2.rare_var_SKAT_WES/"
# WES_missense_del_gt_path ="../../../../hipsci_genotype_processing/data/WES/VEP/missense_clean_GRCh38.vcf.gz"
# WES_missense_del_cleaninfo_path ="../../../../hipsci_genotype_processing/data/WES/VEP/missense_clean_GRCh38.txt"
# line_prop_changes_path = "../../../data/results/phagocytosis/1.check_line_proportions/line_prop_changes_per_well.csv"
# genotype_pc_path = "../../../../OTAR2065_sc_eQTL/data/genotype/plink_genotypes/all_pools.all_donors.genotype.MAF05.eigenvec"
# line_info_path="../../../../OTAR2065_differentiation_efficiency/data/donor_metadata_complete_with_imputed_sex.csv"
# n_gene=9993
# scaled_proliferation_path = "../../../../OTAR2065_WGS_iPSC_premac_micro/data/results/1.2.scale_proliferation/line_prop_changes_microglia_premac.csv"
######

dir.create(output_directory, recursive = T)

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

# line proportion changes (phagocytosis)
line_prop_changes = readr::read_csv(line_prop_changes_path) %>%
  dplyr::mutate(donor = if_else(str_detect(line,"-"),
                                str_split_i(line,pattern="-",1),
                                str_split_i(line,pattern="_",1))) %>%
  dplyr::filter(!scaled_log_fraction %in% c(NA,Inf,-Inf)) %>%
  dplyr::select(scaled_log_fraction,replicate,line,donor,condition,treatment,pool,prop_unadjusted_min_value) %>%
  distinct()

message("There are ", length(unique(line_prop_changes$line)), " lines in the analysis after Na and Inf filters") # 202

# scaled proliferation info from premac- microglia (average across all reps) ########
scaled_proliferation = readr::read_csv(scaled_proliferation_path) %>%
  dplyr::mutate(treatment = case_when(treatment == "IFNg" ~ "IFN",
                                      .default = treatment)) %>%
  dplyr::filter(!scaled_log_fraction %in% c(NA,Inf,-Inf)) %>%
  dplyr::group_by(treatment,line) %>%
  dplyr::summarise(scaled_proliferation = mean(scaled_log_fraction,na.rm = TRUE)) %>%
  dplyr::ungroup()
  

## agregate all phenotype and covariate information

full_pheno_info = line_prop_changes %>%
  dplyr::left_join(line_info) %>%
  dplyr::left_join(scaled_proliferation) %>%
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
nrow(unique(WES_missense_del_info[WES_missense_del_info$category == "deleterious","Gene"])) # 16,348  genes with del mutations
nrow(unique(WES_missense_del_info[WES_missense_del_info$category == "missense_non_deleterious","Gene"])) # 15,523  genes with del mutations
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
  

  
  for(treat in c("untreated","LPS","IFN")){
    skat_res =list()
    
    if(dim(vars_within_mut_type)[1]==0){
      skat_res[[gene_id]]$p.value = NA
      skat_res[[gene_id]]$resampling_pval = NA
      skat_res[[gene_id]]$gene_name = unique(WES_missense_del_info$SYMBOL)[n_gene]
      skat_res[[gene_id]]$gene_id = unique(WES_missense_del_info$Gene)[n_gene]
      
    } else{
      
    
    
    message("Working on ", treat, " treatment")
    message("\n")
    
    message("subset right metadata") 
    
    meta_subset = full_pheno_info %>%
      # filter objects to shared lines
      dplyr::filter((treatment == treat) & (line %in% WES_missense_del_genotype$line)) %>%
      dplyr::mutate(sex_numeric = case_when(Sex == "Male" ~ 1, # female as reference (0)
                                            .default = 0)) %>%
      dplyr::group_by(line) %>%
      # summarising per line, pool, across all reps
      dplyr::reframe(scaled_log_fraction_pool = mean(scaled_log_fraction),
                     prop_unadjusted_min_value_mean = mean(prop_unadjusted_min_value),
                     sex_numeric = sex_numeric, genotypePC1= genotypePC1, genotypePC2 =genotypePC2,
                     scaled_proliferation = scaled_proliferation) %>%
      dplyr::select(line,scaled_log_fraction_pool,prop_unadjusted_min_value_mean,scaled_proliferation,sex_numeric, genotypePC1, genotypePC2) %>%
      distinct() %>%
      tidyr::drop_na()
    
   
    #### important check ##
    # drop genes that don't amount to at 3 lines across all variants in at least 2 categories 
    # (often heterozygote or ALT homozygote don't have them bc deleterious vars are rare)
    WES_missense_del_genotype_subset = WES_missense_del_genotype %>%
      # filter objects to shared lines
      dplyr::filter((ID %in% vars_within_mut_type$ID) & (line %in% meta_subset$line))
    
    #### check both matrices have same lines ####
    
    if(identical(sort(unique(WES_missense_del_genotype_subset$line)),sort(unique(meta_subset$line)))){
      message("Genotype and metadata have same lines")
    } else{
      message("Genotype and metadata don't have same lines!")
      skat_res[[gene_id]]$p.value = NA
      skat_res[[gene_id]]$resampling_pval = NA
      skat_res[[gene_id]]$gene_name = unique(WES_missense_del_info$SYMBOL)[n_gene]
      skat_res[[gene_id]]$gene_id = unique(WES_missense_del_info$Gene)[n_gene]
      next()
    }
    
    
    check_min_lines = sum(colSums(table(WES_missense_del_genotype_subset$ID,WES_missense_del_genotype_subset$genotype))>3)>2
    # if the tibble has < 2 rows (0 or just 1 variant for that gene and mut status), save NA object
    # if it fails the minimum line number across variants, save NA object
    if(dim(vars_within_mut_type)[1]<2 | check_min_lines == FALSE  ){
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
      
      
      #### check again both matrices have same lines ####
      
      if(identical(sort(Z$line),sort(meta_subset$line))){
        message("Genotype and metadata have same lines")
      } else{
        message("Genotype and metadata don't have same lines!")
        skat_res[[gene_id]]$p.value = NA
        skat_res[[gene_id]]$resampling_pval = NA
        skat_res[[gene_id]]$gene_name = unique(WES_missense_del_info$SYMBOL)[n_gene]
        skat_res[[gene_id]]$gene_id = unique(WES_missense_del_info$Gene)[n_gene]
        next()
      }
      
      # subset to needed columns
      Z = as.matrix(Z[,2:ncol(Z)])
      
      message("creating vector of phenotype") 
      
      y.c = meta_subset$scaled_log_fraction_pool
      ## creating covariate matrix ###
      X = as.matrix(data.frame(#pool=meta_subset$pool, # I can't include pool because every element of matrix needs to be same type
        # would need to pre-adjust somehow
        sex=meta_subset$sex_numeric,
        scaled_proliferation = meta_subset$scaled_proliferation,
        prop_unadjusted_min_value_mean = meta_subset$prop_unadjusted_min_value_mean,
        genotypePC1 = meta_subset$genotypePC1,  genotypePC2 = meta_subset$genotypePC2))
      
      skat_res = list()
      

      
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
        prop_unadjusted_min_value_mean = X[,"prop_unadjusted_min_value_mean"]
        scaled_proliferation = X[,"scaled_proliferation"]
        genotypePC1 = X[,"genotypePC1"]
        genotypePC2 = X[,"genotypePC2"]
        obj=SKAT_Null_Model_ChrX(y.c ~ sex + prop_unadjusted_min_value_mean + scaled_proliferation + genotypePC1 + genotypePC2, 
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
    write_rds(skat_res,paste0(output_directory,n_gene,"_",treat,"_",mut_type,"_SKATO_result.RDS"))
    rm(skat_res)
    gc()
  }
  
}
