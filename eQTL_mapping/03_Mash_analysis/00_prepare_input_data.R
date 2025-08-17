
library(tidyverse)
library(reshape2)
library(sessioninfo)


################################################################################
##            3. Multivariate adaptive shrinkage analysis (mash)
##       Analysis of shared eQTL effects in microglia and macrophages*
################################################################################
#  Code to run mash with eQTL effects in stimulated iPSC-derived macrophages to 
#  examine shared effects with microglia.
# * Data from Panousis et al., 2023 (https://doi.org/10.1101/2023.05.29.542425).
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

#-------------------------------------------------------------------------------
#                      3.0 Prepare input data for mash
#-------------------------------------------------------------------------------

## Set working dir
setwd("/lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_Daianna/")

## Define input, output, and plot dirs
inputdir = paste(getwd(), "input_data", "03_Mash_analysis", "00_prepare_input_data", sep = "/")
outdir = paste(getwd(), "output_data", "03_Mash_analysis", "00_prepare_input_data", sep = "/")
plotdir = paste(getwd(), "plots", "03_Mash_analysis", "00_prepare_input_data", sep = "/")
dir.create(inputdir, recursive = T)
dir.create(outdir, recursive = T)
dir.create(plotdir, recursive = T)


# ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ 
##  Step 1: define strong eQTL set across all microglia and macrophage conditions
# ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ 

# ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ 
#  Step 1.1: find significant lead eQTLs per cell type
# ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ 

###########################  Microglia lead eQTLs  #############################

## Load microglia lead eQTL results across all treatments
eqtl_lead_microglia = as.data.frame(readr::read_csv("../OTAR2065_sc_eQTL/data/results/4.Inspect_eQTL_results/tensorQTL_variant_gene.csv")) 
dim(eqtl_lead_microglia)
# [1] 53353    30

## Lead gene-variant pairs 
head(eqtl_lead_microglia, 3)
#           gene_id gene_name num_var beta_shape1 beta_shape2   true_df pval_true_df      variant_id start_distance end_distance ma_samples ma_count        af pval_nominal     slope
# 1 ENSG00000170291      ELP5    1014   1.0196383    69.66986  97.42944 8.478686e-52  17_7242662_T_C          -8755        -8755        111      134 0.3602151 6.539263e-63  1.234492
# 2 ENSG00000109861      CTSC    1843   1.0984722    62.24490  82.98165 9.018975e-42 11_88352106_A_G          87036        87036         97      115 0.3108108 6.848099e-54  1.130912
# 3 ENSG00000240344     PPIL3     810   0.9814995    49.33042 100.22078 6.580654e-45 2_200877622_C_T           6714         6714         72       82 0.2204301 5.042867e-53 -1.353377
#
#     slope_se   pval_perm    pval_beta         qval                        name                          group nPCs           cluster treatment                       toselect    position
# 1 0.03647488 0.000999001 6.324232e-51 1.984476e-47 untreated_Not_proliferating 63_untreated_Not_proliferating   63 Not_proliferating untreated 63-untreated-Not_proliferating  17_7242662
# 2 0.03797634 0.000999001 7.330265e-44 2.615268e-40       IFN_Not_proliferating       73_IFN_Not_proliferating   73 Not_proliferating       IFN       73-IFN-Not_proliferating 11_88352106
# 3 0.04962319 0.000999001 1.998529e-42 1.567793e-39 untreated_Not_proliferating 63_untreated_Not_proliferating   63 Not_proliferating untreated 63-untreated-Not_proliferating 2_200877622
#
#   swapped_tensor_alleles swapped_tensor_in_VCF swapped_slope     gene_variant_pair
# 1         17_7242662_C_T                 FALSE            NA   ELP5-17_7242662_T_C
# 2        11_88352106_G_A                 FALSE            NA  CTSC-11_88352106_A_G
# 3        2_200877622_T_C                 FALSE            NA PPIL3-2_200877622_C_T


## Define ensemblID-varID (some gene names are repeated!!!)
eqtl_lead_microglia$ensembl_id_var <- paste(eqtl_lead_microglia$gene_id, eqtl_lead_microglia$variant_id, sep = "-")

## Subset to non-prolif cluster
eqtl_lead_microglia_np <- subset(eqtl_lead_microglia, cluster == 'Not_proliferating')
dim(eqtl_lead_microglia_np)
# [1] 32067    31

## All gene-variant pairs are unique within each treatment
by(eqtl_lead_microglia_np$ensembl_id_var, eqtl_lead_microglia_np$treatment, function(x){which(duplicated(x))})

# eqtl_microglia_np$treatment: IFN
# integer(0)
# --------------------------------------------------------------------------------------------------------------------------------------------- 
# eqtl_microglia_np$treatment: LPS
# integer(0)
# --------------------------------------------------------------------------------------------------------------------------------------------- 
# eqtl_microglia_np$treatment: untreated
# integer(0)

## Num variants per gene per treatment
## * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#  Warning:                     !!!                              *
#  We don't want genes with more than 1 variant within the same  *
#  condition. These variants may be correlated (i.e. in LD).     *
#  There were not such cases.                                    *
## * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
## The top (lead) variant per gene per treatment
table(eqtl_lead_microglia_np$gene_id, eqtl_lead_microglia_np$treatment) %>% table()
#    0     1 
# 1629 32067 

## Subset to signif eQTLs 
eqtl_lead_microglia_signif <- eqtl_lead_microglia_np %>% filter(qval<0.05)
save(eqtl_lead_microglia_signif, file = paste0(outdir, "/lead_signif_eQTLs_microglia_all_treatments.Rdata"))
dim(eqtl_lead_microglia_signif)
# [1] 14053    31

## Unique gene-variant pairs within each treatment
by(eqtl_lead_microglia_signif$ensembl_id_var, eqtl_lead_microglia_signif$treatment, function(x){which(duplicated(x))})

# eqtl_microglia_signif$treatment: IFN
# integer(0)
# --------------------------------------------------------------------------------------------------------------------------------------------- 
#   eqtl_microglia_signif$treatment: LPS
# integer(0)
# --------------------------------------------------------------------------------------------------------------------------------------------- 
#   eqtl_microglia_signif$treatment: untreated
# integer(0)


## Confirm 1 variant per gene per treatment
table(eqtl_lead_microglia_signif$gene_id, eqtl_lead_microglia_signif$treatment) %>% table()
#    0     1 
# 7166 14053 

## signif lead eQTLs per treatment
table(eqtl_lead_microglia_signif$treatment)
#  IFN       LPS untreated 
# 4454      4464      5135 

## Number of eQTLs unique/shared across treatments
rowSums(table(eqtl_lead_microglia_signif$ensembl_id_var, eqtl_lead_microglia_signif$treatment)) %>% table
# Num treatments:    1     2     3 
# Num of eQTLs:  11217  1013   270 

## Lead signif eQTLs in at least one condition in microglia
signif_lead_gene_var_pairs_microglia <- unique(eqtl_lead_microglia_signif$ensembl_id_var)
length(signif_lead_gene_var_pairs_microglia)
# [1] 12500

## Chr and pos of variants in strong eQTLs
variants_strong_eQTLs_microglia <- t(sapply(signif_lead_gene_var_pairs_microglia, function(x){strsplit(strsplit(x, "-")[[1]][2], "_")[[1]][1:2]})) %>% as.data.frame()
rownames(variants_strong_eQTLs_microglia) <- colnames(variants_strong_eQTLs_microglia) <- NULL
write.table(variants_strong_eQTLs_microglia, file = paste0(outdir, "/variants_strong_eQTLs_microglia.txt"), row.names = F, quote = F, sep = "\t")

## rsIDs and ALL alleles for strong eQTLs variants (output from find_variant_rsIDs.sh) 
strong_eQTLs_rsIDs_microglia <- read_table(paste0(inputdir, "/dbsnp.156.variants_strong_eQTLs_microglia.subset.tsv"))
dim(strong_eQTLs_rsIDs_microglia)
# 68955     5     (taking all Alt alleles)
colnames(strong_eQTLs_rsIDs_microglia) <- c("chr", "pos", "rsID", "Ref", "Alt") 
strong_eQTLs_rsIDs_microglia$varID_alleles <- paste(strong_eQTLs_rsIDs_microglia$chr, 
                                                    strong_eQTLs_rsIDs_microglia$pos,
                                                    strong_eQTLs_rsIDs_microglia$Ref,
                                                    strong_eQTLs_rsIDs_microglia$Alt, sep = "_")


## Annotate rsIDs for strong eQTL variants in microglia
genes <- sapply(signif_lead_gene_var_pairs_microglia, function(x){strsplit(x, "-")[[1]][1]}) %>% as.vector()
varID_alleles <- sapply(signif_lead_gene_var_pairs_microglia, function(x){strsplit(x, "-")[[1]][2]}) %>% as.vector()
chr <- sapply(varID_alleles, function(x){strsplit(x, "_")[[1]][1]}) %>% as.vector()
pos <- sapply(varID_alleles, function(x){strsplit(x, "_")[[1]][2]}) %>% as.vector()

strong_eQTLs_alleles_rsIDs_microglia <- as.data.frame(cbind("eQTL_id_alleles" = signif_lead_gene_var_pairs_microglia, 
                                                            "ensmebl_id" = genes, 
                                                            "chr" = chr,
                                                            "pos" = pos,
                                                            "varID_alleles" = varID_alleles))
strong_eQTLs_alleles_rsIDs_microglia$rsID <- strong_eQTLs_rsIDs_microglia[match(strong_eQTLs_alleles_rsIDs_microglia$varID_alleles, strong_eQTLs_rsIDs_microglia$varID_alleles), "rsID"] %>% unlist
dim(strong_eQTLs_alleles_rsIDs_microglia)
# [1] 12500     6

## Discard eQTLs with no rsID available (not findable in macrophage data)
strong_eQTLs_alleles_rsIDs_microglia <- strong_eQTLs_alleles_rsIDs_microglia[-which(is.na(strong_eQTLs_alleles_rsIDs_microglia$rsID)), ]
strong_eQTLs_alleles_rsIDs_microglia$eQTL_id_rsID <- paste(strong_eQTLs_alleles_rsIDs_microglia$ensmebl_id, 
                                                           strong_eQTLs_alleles_rsIDs_microglia$chr,
                                                           strong_eQTLs_alleles_rsIDs_microglia$pos,
                                                           strong_eQTLs_alleles_rsIDs_microglia$rsID, sep = "_")
head(strong_eQTLs_alleles_rsIDs_microglia, 3)
#                   eQTL_id_alleles      ensmebl_id chr       pos   varID_alleles       rsID                           eQTL_id_rsID
# 1  ENSG00000170291-17_7242662_T_C ENSG00000170291  17   7242662  17_7242662_T_C   rs222843    ENSG00000170291_17_7242662_rs222843
# 2 ENSG00000109861-11_88352106_A_G ENSG00000109861  11  88352106 11_88352106_A_G rs11019479 ENSG00000109861_11_88352106_rs11019479
# 3 ENSG00000240344-2_200877622_C_T ENSG00000240344   2 200877622 2_200877622_C_T  rs2136600  ENSG00000240344_2_200877622_rs2136600

## Only one Alt allele per rsID? No
table(strong_eQTLs_alleles_rsIDs_microglia$eQTL_id_rsID) %>%  table
#    1    2 
# 9160    1 

## Delete eQTLs sharing rsID = ambiguous in macromap
r <- which(strong_eQTLs_alleles_rsIDs_microglia$eQTL_id_rsID == which(table(strong_eQTLs_alleles_rsIDs_microglia$eQTL_id_rsID) > 1) %>% names)
strong_eQTLs_alleles_rsIDs_microglia[r, ]
#                              eQTL_id_alleles      ensmebl_id chr       pos               varID_alleles       rsID
#  ENSG00000138378-2_191165320_A_AACTCTCGGAGCG ENSG00000138378   2 191165320 2_191165320_A_AACTCTCGGAGCG rs11269783
#  ENSG00000138378-2_191165320_A_AACTCACGGAGCG ENSG00000138378   2 191165320 2_191165320_A_AACTCACGGAGCG rs11269783
#                            eQTL_id_rsID
#  ENSG00000138378_2_191165320_rs11269783
#  ENSG00000138378_2_191165320_rs11269783

strong_eQTLs_alleles_rsIDs_microglia <- strong_eQTLs_alleles_rsIDs_microglia[-r, ]

## Final set: unique
dim(strong_eQTLs_alleles_rsIDs_microglia)
# [1] 9160    7
save(strong_eQTLs_alleles_rsIDs_microglia, file = paste0(outdir, "/strong_eQTLs_alleles_rsIDs_microglia.Rdata"))
which(duplicated(strong_eQTLs_alleles_rsIDs_microglia$eQTL_id_alleles))
# integer(0)


###########################  Macrophage lead eQTLs  ############################
## Data copied from: /lustre/scratch123/hgi/projects/otar2065/resources/macromap_Panousis_2023/eQTL/permuted_summary/

## Treatment x hr lead eQTLs
eqtl_lead_macro <- list()

for (treatment in c("IFNG", "sLPS", "Ctrl")){
  
  for(hr in c("6", "24")){
    lead_eQTLs <- as.data.frame(read.table(file = paste0(inputdir, "/macromap_Panousis_2023_data/", treatment, "_", hr, ".permuted.txt"), header = T))
    ## Chr num
    lead_eQTLs$chr <- gsub("chr", "", lead_eQTLs$Chromosome_phe)
    ## Ensemble ID
    lead_eQTLs$ensembl_id <- gsub("\\.[0-9]*", "", lead_eQTLs$Phenotype_ID)
    ## Gene-varID
    lead_eQTLs$gene_var_id <- paste(lead_eQTLs$ensembl_id, lead_eQTLs$chr, lead_eQTLs$Pos_variant, lead_eQTLs$Pos_variant_end, 
                                    lead_eQTLs$Best_variant_in_cis, sep = "-")
    
    eqtl_lead_macro[[paste(treatment, hr, sep = "_")]] <- lead_eQTLs
  }
}

## Same genes across all conditions 
lapply(eqtl_lead_macro, dim) %>%  unique
# [1] 14034    22

## Unique gene-var pairs within each condition? Yes
lapply(eqtl_lead_macro, function(x){which(duplicated(x[,"gene_var_id"]))}) %>%  unique
# integer(0)

## Top variant per gene within each condition? Yes!
lapply(eqtl_lead_macro, function(x){table(table(x[,"Phenotype_ID"]))}) %>%  unique
#     1 
# 14034 

## Signif eQTLs in each treatment x hr
eqtl_lead_macro_signif <- lapply(eqtl_lead_macro, function(x){subset(x, Corrected_pvalue < 0.05)[, c("Phenotype_ID", "ensembl_id", "chr", "Pos_variant", 
                                                                                                     "Pos_variant_end", "Best_variant_in_cis", "gene_var_id")]})

## Signif lead eQTLs per condition
lapply(eqtl_lead_macro_signif, function(l){dim(l)[1]}) %>% unlist
# IFNG_6 IFNG_24  sLPS_6 sLPS_24  Ctrl_6 Ctrl_24 
#   3953    4166    3704    3920    3824    4148 

## Num signif eQTLs in at least one condition
all_eqtl_lead_macro_signif <- do.call("rbind", eqtl_lead_macro_signif) %>% unique() 
dim(all_eqtl_lead_macro_signif)
# [1] 19902     7

## Chr and position of variants in all lead signif eQTLs
variants_strong_eQTLs_macrophages <- do.call("rbind", lapply(eqtl_lead_macro_signif, function(x){unique(x[, c("chr", "Pos_variant")])}) ) %>% unique() 
variants_strong_eQTLs_macrophages[,2] <- as.character(variants_strong_eQTLs_macrophages[,2])
rownames(variants_strong_eQTLs_macrophages) <- colnames(variants_strong_eQTLs_macrophages) <- NULL
dim(variants_strong_eQTLs_macrophages)
# [1] 19379     2
write.table(variants_strong_eQTLs_macrophages, file = paste0(outdir, "/variants_strong_eQTLs_macrophages.txt"), row.names = F, quote = F, sep = "\t")

## Alleles of strong eQTL variants (output from find_variant_alleles.sh)
strong_eQTLs_alleles_macrophage <- read_table(paste0(inputdir, "/dbsnp.156.variants_lead_signif_macromap_Panousis_2023.subset.tsv"))
colnames(strong_eQTLs_alleles_macrophage) <- c("chr", "pos", "rsID", "Ref", "Alt") 
strong_eQTLs_alleles_macrophage$chr_pos_rsID <- paste(strong_eQTLs_alleles_macrophage$chr, 
                                                      strong_eQTLs_alleles_macrophage$pos,
                                                      strong_eQTLs_alleles_macrophage$rsID, sep = "-")
dim(unique(strong_eQTLs_alleles_macrophage[,1:2]))
# [1] 19370     2

## % multiallelic variants
length(which(table(strong_eQTLs_alleles_macrophage$chr_pos_rsID) >1)) / 19370 * 100
# [1] 66.05059

## Expand multiallelic variants in eQTLs
all_eqtl_lead_macro_signif$chr_pos_rsID <- paste(all_eqtl_lead_macro_signif$chr, 
                                                 all_eqtl_lead_macro_signif$Pos_variant, 
                                                 all_eqtl_lead_macro_signif$Best_variant_in_cis, sep = "-")

expanded_eqtl_lead_macro_signif <- vector()
for(i in 1:dim(all_eqtl_lead_macro_signif)[1]){
  
  eQTL <- all_eqtl_lead_macro_signif[i, ]
  
  ## Expand to all Alt alleles for each chr-pos-rsID per eQTL  
  alleles <- strong_eQTLs_alleles_macrophage[which(strong_eQTLs_alleles_macrophage$chr_pos_rsID == eQTL$chr_pos_rsID), c("Ref", "Alt")]
  
  if(dim(alleles)[1] == 0){alleles <- data.frame("Ref" = NA, "Alt"= NA)}
  
  expanded_eqtl_lead_macro_signif <- rbind(expanded_eqtl_lead_macro_signif, cbind(eQTL, alleles))
}

dim(expanded_eqtl_lead_macro_signif)
# [1] 36609    10

## Remove eQTLs without allele data
expanded_eqtl_lead_macro_signif <- expanded_eqtl_lead_macro_signif[-which(is.na(expanded_eqtl_lead_macro_signif$Ref) | is.na(expanded_eqtl_lead_macro_signif$Alt)), ]
dim(expanded_eqtl_lead_macro_signif)
# [1] 35338    12

expanded_eqtl_lead_macro_signif$var_id <- paste(expanded_eqtl_lead_macro_signif$chr,
                                                expanded_eqtl_lead_macro_signif$Pos_variant,
                                                expanded_eqtl_lead_macro_signif$Ref,
                                                expanded_eqtl_lead_macro_signif$Alt, sep = "_")

## Unique eQTLs!
expanded_eqtl_lead_macro_signif$ensembl_id_var <- paste(expanded_eqtl_lead_macro_signif$ensembl_id, expanded_eqtl_lead_macro_signif$var_id, sep = "-")
which(duplicated(expanded_eqtl_lead_macro_signif$ensembl_id_var))
# integer(0)
save(expanded_eqtl_lead_macro_signif, file = paste0(outdir, "/expanded_eqtl_lead_macro_signif.Rdata"))


## eQTLs from multiallelic variants
table(table(expanded_eqtl_lead_macro_signif$gene_var_id))
#    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22   23   25   26   29   30 
# 8025 6566 3374  239  163  102   53   31   18   12    8    7    4    5    3    1    3    2    5    2    1    2    1    1    1    1    1 


## Keep alleles found in microglia nominal data only 
all_microglia_eQTLs = vector()
for(treatment in c("IFN","LPS", "untreated")){
  
  eQTLs = readr::read_delim(paste0("/lustre/scratch123/hgi/projects/otar2065/OTAR2065_sc_eQTL/data/results/tensorqtl/best_results/",
                                              "sum_sizefactorsNorm_log2_scaled_centered_",
                                              treatment, "_Not_proliferating_common_500kb_window_tensorQTL_nominal.txt")) %>% 
    dplyr::mutate(ensembl_id_var = paste(phenotype_id, variant_id, sep = "-")) %>% 
    dplyr::select(ensembl_id_var) %>% as.data.frame() %>% unlist
  
  all_microglia_eQTLs <- append(all_microglia_eQTLs, as.vector(eQTLs))
}

all_microglia_eQTLs <- unique(all_microglia_eQTLs)
length(all_microglia_eQTLs)
# [1] 12430020

strong_eQTLs_alleles_rsIDs_macrophages <- subset(expanded_eqtl_lead_macro_signif, ensembl_id_var %in% all_microglia_eQTLs)
dim(strong_eQTLs_alleles_rsIDs_macrophages)
# [1] 9973   12

## eQTLs with >1 allele present in microglia
table(strong_eQTLs_alleles_rsIDs_macrophages$gene_var_id) %>%  table
#    1    2 
# 9937   18

## Remove eQTLs sharing rsID
ambiguous_rsIDs <- which(table(strong_eQTLs_alleles_rsIDs_macrophages$gene_var_id) >1) %>% names()
strong_eQTLs_alleles_rsIDs_macrophages <- strong_eQTLs_alleles_rsIDs_macrophages[-which(strong_eQTLs_alleles_rsIDs_macrophages$gene_var_id %in% ambiguous_rsIDs), ]
dim(strong_eQTLs_alleles_rsIDs_macrophages)
# [1] 9937   12
strong_eQTLs_alleles_rsIDs_macrophages$eQTL_id_alleles <- strong_eQTLs_alleles_rsIDs_macrophages$ensembl_id_var
strong_eQTLs_alleles_rsIDs_macrophages$eQTL_id_rsID <- paste(strong_eQTLs_alleles_rsIDs_macrophages$ensembl_id, 
                                                             gsub("-", "_", strong_eQTLs_alleles_rsIDs_macrophages$chr_pos_rsID), sep = "_")
save(strong_eQTLs_alleles_rsIDs_macrophages, file = paste0(outdir, "/strong_eQTLs_alleles_rsIDs_macrophages.Rdata"))


## Strong lead eQTLs in microglia and macrophages
strong_eQTLs <- rbind(strong_eQTLs_alleles_rsIDs_microglia[, c("eQTL_id_alleles", "eQTL_id_rsID")], 
                      strong_eQTLs_alleles_rsIDs_macrophages[, c("eQTL_id_alleles", "eQTL_id_rsID")]) %>% as.data.frame()
strong_eQTLs <- unique(strong_eQTLs)
dim(strong_eQTLs)
# [1] 18491     2

## No diff rsID for same alleles
which(duplicated(strong_eQTLs$eQTL_id_alleles))
# integer(0)
## No diff alleles for same rsID
which(duplicated(strong_eQTLs$eQTL_id_rsID))
# integer(0)

save(strong_eQTLs, file = paste0(outdir, "/strong_eQTLs_alleles_rsIDs.Rdata"))


# ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ 
#  Step 1.2: extract effects and se for strong eQTLs across all conditions and cell types
# ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ ~ __ 

###########################  Microglia nominal eQTLs  #############################
microglia_nominal_all = list()

for(treatment in c("IFN","LPS", "untreated")){
  
  microglia_nominal_all[[treatment]] = readr::read_delim(paste0("/lustre/scratch123/hgi/projects/otar2065/OTAR2065_sc_eQTL/data/results/tensorqtl/best_results/",
                                                             "sum_sizefactorsNorm_log2_scaled_centered_",
                                                             treatment, "_Not_proliferating_common_500kb_window_tensorQTL_nominal.txt")) %>% 
    dplyr::mutate(treatment = treatment) %>%
    dplyr::mutate(eQTL_id_alleles = paste(phenotype_id, variant_id, sep = "-")) %>%
    dplyr::mutate(ensembl_id = phenotype_id) %>% 
    dplyr::mutate(effect = slope) %>% 
    dplyr::mutate(se = slope_se) %>% 
    dplyr::select(ensembl_id, eQTL_id_alleles, effect, se) 
}

## All eQTLs with effects and se in 3 treatments
intersect <- intersect(intersect(microglia_nominal_all$IFN$eQTL_id_alleles, microglia_nominal_all$LPS$eQTL_id_alleles), 
                 microglia_nominal_all$untreated$eQTL_id_alleles)
length(intersect)
# [1] 11251961

microglia_nominal_joint <- inner_join(microglia_nominal_all$IFN, microglia_nominal_all$LPS, 
                                     by = c("ensembl_id", "eQTL_id_alleles"), suffix = c("_IFN", "_LPS")) %>% 
                          inner_join(., microglia_nominal_all$untreated, by = c("ensembl_id", "eQTL_id_alleles")) %>% as.data.frame()

colnames(microglia_nominal_joint)[7:8] <- paste0(colnames(microglia_nominal_joint)[7:8], "_untreated")
## Unique eQTLs
which(duplicated(microglia_nominal_joint$eQTL_id_alleles))
# integer(0)
dim(microglia_nominal_joint)
# [1] 11251961        8

## Remove eQTLs from multiallelic variants
microglia_nominal_joint$gene_chr_pos_Ref <- gsub("_[ATCG]*$", "", microglia_nominal_joint$eQTL_id_alleles)
table(microglia_nominal_joint$gene_chr_pos_Ref) %>% table
#        1        2        3        4        5        6 
# 11021081    97674    10416      997       52        6 
mutiallelic_variant_eQTLs <- which(table(microglia_nominal_joint$gene_chr_pos_Ref) > 1) %>% names()
length(mutiallelic_variant_eQTLs)
# [1] 109145
microglia_nominal_joint <- microglia_nominal_joint[-which(microglia_nominal_joint$gene_chr_pos_Ref %in% mutiallelic_variant_eQTLs), ]
dim(microglia_nominal_joint)
# [1] 11021081        9

## Add rsIDs 
microglia_nominal_joint$varID_alleles <- gsub(".*-", "", microglia_nominal_joint$eQTL_id_alleles)
microglia_nominal_joint$chr <- gsub("_.*", "", gsub(".*-", "", microglia_nominal_joint$eQTL_id_alleles))
microglia_nominal_joint$pos <- gsub("_.*", "", gsub(".*-[0-9]*_", "", microglia_nominal_joint$eQTL_id_alleles))
variants_all_microglia_eQTLs <- unique(microglia_nominal_joint[, c("chr", "pos")])
colnames(variants_all_microglia_eQTLs) <- NULL
dim(variants_all_microglia_eQTLs)
# [1] 3478933       2
write.table(variants_all_microglia_eQTLs, file = paste0(outdir, "/variants_all_microglia_eQTLs.txt"), row.names = F, quote = F, sep = "\t")

## Output from find_variant_rsIDs.sh
rsIDs_all_microglia <- read_table(paste0(inputdir, "/dbsnp.156.variants_all_microglia_eQTLs.subset.tsv"))
dim(rsIDs_all_microglia)
# [1] 10036481        5
colnames(rsIDs_all_microglia) <- c("chr", "pos", "rsID", "Ref", "Alt") 
rsIDs_all_microglia$varID_alleles <- paste(rsIDs_all_microglia$chr, 
                                           rsIDs_all_microglia$pos,
                                           rsIDs_all_microglia$Ref,
                                           rsIDs_all_microglia$Alt, sep = "_")
rsIDs_all_microglia$varID_rsID <- paste(rsIDs_all_microglia$chr, 
                                        rsIDs_all_microglia$pos,
                                        rsIDs_all_microglia$rsID, sep = "_")

microglia_nominal_joint$varID_rsID <- unlist(rsIDs_all_microglia[match(microglia_nominal_joint$varID_alleles, rsIDs_all_microglia$varID_alleles), "varID_rsID"]) 

## Discard eQTLs without rsID
microglia_nominal_joint <- microglia_nominal_joint[-which(is.na(microglia_nominal_joint$varID_rsID)), ]
dim(microglia_nominal_joint)
# [1] 9688468      13

## Multiple alleles of same variant left? No! 
microglia_nominal_joint$eQTL_id_rsID <- paste(microglia_nominal_joint$ensembl_id, microglia_nominal_joint$varID_rsID, sep = "-")
table(microglia_nominal_joint$eQTL_id_rsID) %>% table
#       1 
# 9688468 

## Save
save(microglia_nominal_joint, file = paste0(outdir, "/microglia_nominal_joint.Rdata"))


###########################  Macrophage nominal eQTLs  #############################
## Data copied from: /lustre/scratch123/hgi/projects/otar2065/resources/macromap_Panousis_2023/eQTL/nominal/
## Search all strong eQTLs in each condition nominal data in macrophages

macrophage_nominal_all = list()

for(treatment in c("IFNG","sLPS", "Ctrl")){
  for(hr in c("6", "24")){
    macrophage_nominal_all[[paste(treatment, hr, sep = "_")]] = readr::read_delim(paste0(inputdir, "/macromap_Panousis_2023_data/",
                                                                      treatment, "_", hr, "_1MB_all.stderr.txt")) %>% 
      dplyr::mutate(chr = gsub("chr", "", Chromosome_phe)) %>%
      dplyr::mutate(pos = Pos_variant)%>%
      dplyr::mutate(ensembl_id = gsub("\\.[0-9]*", "", Phenotype_ID)) %>% 
      dplyr::mutate(rsID = gsub("_.*", "", Variant_in_cis)) %>%
      dplyr::mutate(varID_rsID = paste(chr, pos, rsID, sep = "_")) %>%
      dplyr::mutate(eQTL_id_rsID = paste(ensembl_id, varID_rsID, sep = "-")) %>%
      dplyr::mutate(effect = Beta_regression) %>% 
      dplyr::mutate(se = std.err) %>% 
      dplyr::select(ensembl_id, eQTL_id_rsID, chr, pos, rsID, varID_rsID, effect, se)
  }
  
}

save(macrophage_nominal_all, file = paste0(outdir, "/macrophage_nominal_all.Rdata"))

## All eQTLs with effects and se in 6 conditions
intersect1 <- intersect(macrophage_nominal_all$IFNG_6$eQTL_id_rsID, macrophage_nominal_all$IFNG_24$eQTL_id_rsID)
intersect2 <- intersect(macrophage_nominal_all$sLPS_6$eQTL_id_rsID, macrophage_nominal_all$sLPS_24$eQTL_id_rsID)
intersect3 <- intersect(macrophage_nominal_all$Ctrl_6$eQTL_id_rsID, macrophage_nominal_all$Ctrl_24$eQTL_id_rsID)
intersect4 <- intersect(intersect(intersect1, intersect2), intersect3)

length(intersect4)
# [1] 63348409

macrophage_nominal_joint <- inner_join(macrophage_nominal_all$IFNG_6, macrophage_nominal_all$IFNG_24, 
                                      by = c("eQTL_id_rsID", "ensembl_id", "chr", "pos", "rsID", "varID_rsID"), 
                                      suffix = c("_IFNG_6", "_IFNG_24")) %>% 
                            inner_join(., macrophage_nominal_all$sLPS_6, by = c("eQTL_id_rsID", "ensembl_id", "chr", "pos", "rsID", "varID_rsID")) %>% 
                            inner_join(., macrophage_nominal_all$sLPS_24, by = c("eQTL_id_rsID", "ensembl_id", "chr", "pos", "rsID", "varID_rsID")) %>% 
                            inner_join(., macrophage_nominal_all$Ctrl_6, by = c("eQTL_id_rsID", "ensembl_id", "chr", "pos", "rsID", "varID_rsID")) %>% 
                            inner_join(., macrophage_nominal_all$Ctrl_24, by = c("eQTL_id_rsID", "ensembl_id", "chr", "pos", "rsID", "varID_rsID")) %>% as.data.frame()
  

colnames(macrophage_nominal_joint)[11:12] <- c("effect_sLPS_6", "se_sLPS_6")
colnames(macrophage_nominal_joint)[13:14] <- c("effect_sLPS_24", "se_sLPS_24")
colnames(macrophage_nominal_joint)[15:16] <- c("effect_Ctrl_6", "se_Ctrl_6")
colnames(macrophage_nominal_joint)[17:18] <- c("effect_Ctrl_24", "se_Ctrl_24")

## Unique eQTLs
which(duplicated(macrophage_nominal_joint$eQTL_id_rsID))
# integer(0)
dim(macrophage_nominal_joint)
# [1] 63348409        18
save(macrophage_nominal_joint, file = paste0(outdir, "/macrophage_nominal_joint.Rdata"))



## Search microglia eQTLs in macrophage data
macrophage_nominal_joint$pos <- as.character(macrophage_nominal_joint$pos)
complete_microglia_macrophage_eQTLs <- inner_join(macrophage_nominal_joint, microglia_nominal_joint, by = c("eQTL_id_rsID", "ensembl_id", 
                                                                                                            "chr", "pos", "varID_rsID"))
dim(complete_microglia_macrophage_eQTLs)
# [1] 8719455      27

## One microglia allele --> one macrophage rsID
which(duplicated(complete_microglia_macrophage_eQTLs$eQTL_id_rsID))
# integer(0)
which(duplicated(complete_microglia_macrophage_eQTLs$gene_chr_pos_Ref))
# integer(0)
save(complete_microglia_macrophage_eQTLs, file = paste0(outdir, "/complete_microglia_macrophage_eQTLs.Rdata"))



## Extract effect and se of strong eQTLs
B_hat_strong_set <- subset(complete_microglia_macrophage_eQTLs, eQTL_id_alleles %in% strong_eQTLs$eQTL_id_alleles)[, c("eQTL_id_alleles", 
                                                                                                                       "eQTL_id_rsID", 
                                                                                                                       "ensembl_id", 
                                                                                                                       paste("effect", c("IFN", "LPS", "untreated", 
                                                                                                                                       "IFNG_6", "IFNG_24", 
                                                                                                                                       "sLPS_6", "sLPS_24", 
                                                                                                                                       "Ctrl_6", "Ctrl_24"), sep = "_"))]
dim(B_hat_strong_set)
# [1] 16427     12
save(B_hat_strong_set, file = paste0(outdir, "/B_hat_strong_set.Rdata"))


S_strong_set <- subset(complete_microglia_macrophage_eQTLs, eQTL_id_alleles %in% strong_eQTLs$eQTL_id_alleles)[, c("eQTL_id_alleles", 
                                                                                                                       "eQTL_id_rsID", 
                                                                                                                       "ensembl_id", 
                                                                                                                       paste("se", c("IFN", "LPS", "untreated", 
                                                                                                                                         "IFNG_6", "IFNG_24", 
                                                                                                                                         "sLPS_6", "sLPS_24", 
                                                                                                                                         "Ctrl_6", "Ctrl_24"), sep = "_"))]
dim(S_strong_set)
# [1] 16427     12
save(S_strong_set, file = paste0(outdir, "/S_strong_set.Rdata"))



# ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ 
##  Step 2: extract random sample from nominal eQTL data
# ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ 

## Genes of strong eQTLs
strong_genes <- unique(gsub("-.*", "", strong_eQTLs$eQTL_id_alleles))
## In B_hat and S
strong_genes <- strong_genes %>% subset(. %in% B_hat_strong_set$ensembl_id)
length(strong_genes)
# [1] 6433

## Select 5 random eQTLs from strong eQTL genes
random_eQTLs_per_gene <- sapply(strong_genes, function(gene){ sample(which(complete_microglia_macrophage_eQTLs$ensembl_id == gene), size = 5, replace = F) }) %>% 
                                  unlist() %>% as.vector()
BS_random_set <- complete_microglia_macrophage_eQTLs[random_eQTLs_per_gene, ] 
dim(BS_random_set)
# [1] 32165    27

B_random_set <- BS_random_set[, c("eQTL_id_alleles", "eQTL_id_rsID", "ensembl_id", 
                                  paste("effect", c("IFN", "LPS", "untreated", "IFNG_6", "IFNG_24", 
                                                  "sLPS_6", "sLPS_24", "Ctrl_6", "Ctrl_24"), sep = "_"))] 
S_random_set <- BS_random_set[, c("eQTL_id_alleles", "eQTL_id_rsID", "ensembl_id", 
                                  paste("se", c("IFN", "LPS", "untreated", "IFNG_6", "IFNG_24", 
                                                    "sLPS_6", "sLPS_24", "Ctrl_6", "Ctrl_24"), sep = "_"))] 
dim(B_random_set)
# [1] 32165    12

save(B_random_set, file = paste0(outdir, "/B_random_set.Rdata"))
save(S_random_set, file = paste0(outdir, "/S_random_set.Rdata"))




# ------------------------------------------------------------------------------
## Explore random and strong eQTL effects and std errors

## B_hat with 32K random eQTLs
B_hat_random <- B_random_set[, -c(1:3)] %>% as.matrix()
## S_hat with 32K random eQTLs
S_hat_random <- S_random_set[, -c(1:3)] %>% as.matrix()

pdf(file = paste0(plotdir, "/Bhat_random_set.pdf"), height = 10, width = 15)
par(mfrow=c(3, 4))

labels = c(paste("Microglia", c("IFN", "LPS", "untreated")), paste("Macrophage", c("IFN 6hr", "IFN 24hr", 
                                                                                   "sLPS 6hr", "sLPS 24hr", 
                                                                                   "Ctrl 6hr", "Ctrl 24hr")))
names(labels) <- colnames(B_hat_random)
## Effects of random eQTLs
hist(B_hat_random,  main = "All microglia and macrophage conditions", xlab = "random eQTL effects")
sapply(colnames(B_hat_random), function(col){hist(B_hat_random[, col], main = labels[col], xlab = "random eQTL effects")})
dev.off()

## Std errors of random eQTLs
pdf(file = paste0(plotdir, "/Shat_random_set.pdf"), height = 10, width = 15)
par(mfrow=c(3, 4))

labels = c(paste("Microglia", c("IFN", "LPS", "untreated")), paste("Macrophage", c("IFN 6hr", "IFN 24hr", 
                                                                                   "sLPS 6hr", "sLPS 24hr", 
                                                                                   "Ctrl 6hr", "Ctrl 24hr")))
names(labels) <- colnames(S_hat_random)

hist(S_hat_random,  main = "All microglia and macrophage conditions", xlab = "random eQTL std errors")
sapply(colnames(S_hat_random), function(col){hist(S_hat_random[, col], main = labels[col], xlab = "random eQTL std errors")})
dev.off()



## B_hat with 16K strong eQTLs
B_hat_strong <- as.matrix(B_hat_strong_set[, -c(1:3)])
## S_hat with 16K strong eQTLs
S_hat_strong <- as.matrix(S_strong_set[, -c(1:3)])

## Effects of strong eQTLs
pdf(file = paste0(plotdir, "/Bhat_strong_set.pdf"), height = 10, width = 15)
par(mfrow=c(3, 4))

labels = c(paste("Microglia", c("IFN", "LPS", "untreated")), paste("Macrophage", c("IFN 6hr", "IFN 24hr", 
                                                                                   "sLPS 6hr", "sLPS 24hr", 
                                                                                   "Ctrl 6hr", "Ctrl 24hr")))
names(labels) <- colnames(B_hat_strong)

hist(B_hat_strong,  main = "All microglia and macrophage conditions", xlab = "strong eQTL effects")
sapply(colnames(B_hat_strong), function(col){hist(B_hat_strong[, col], main = labels[col], xlab = "strong eQTL effects")})
dev.off()

## Std errors of strong eQTLs
pdf(file = paste0(plotdir, "/Shat_strong_set.pdf"), height = 10, width = 15)
par(mfrow=c(3, 4))

labels = c(paste("Microglia", c("IFN", "LPS", "untreated")), paste("Macrophage", c("IFN 6hr", "IFN 24hr", 
                                                                                   "sLPS 6hr", "sLPS 24hr", 
                                                                                   "Ctrl 6hr", "Ctrl 24hr")))
names(labels) <- colnames(S_hat_strong)

hist(S_hat_strong,  main = "All microglia and macrophage conditions", xlab = "strong eQTL std errors")
sapply(colnames(S_hat_strong), function(col){hist(S_hat_strong[, col], main = labels[col], xlab = "strong eQTL std errors")})
dev.off()







## Reproducibility info
options(width = 120)
session_info()
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.1 (2023-06-16)
# os       Ubuntu 22.04.4 LTS
# system   x86_64, linux-gnu
# ui       RStudio
# language (EN)
# collate  en_GB.UTF-8
# ctype    en_GB.UTF-8
# tz       Europe/Belfast
# date     2024-11-25
# rstudio  2024.04.0+735 Chocolate Cosmos (server)
# pandoc   3.1.12.3 @ /opt/view/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [1] CRAN (R 4.3.1)
# ashr                 * 2.2-63    2023-08-21 [1] CRAN (R 4.3.1)
# assertthat             0.2.1     2019-03-21 [1] CRAN (R 4.3.1)
# Biobase              * 2.62.0    2023-10-24 [1] Bioconductor
# BiocGenerics         * 0.48.1    2023-11-01 [1] Bioconductor
# bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.3.1)
# cli                    3.6.2     2023-12-11 [1] CRAN (R 4.3.1)
# colorspace             2.1-0     2023-01-23 [1] CRAN (R 4.3.1)
# cowplot              * 1.1.3     2024-01-22 [1] CRAN (R 4.3.1)
# crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.1)
# DelayedArray           0.28.0    2023-10-24 [1] Bioconductor
# dplyr                * 1.1.4     2023-11-17 [1] CRAN (R 4.3.1)
# edgeR                  4.0.16    2024-02-18 [1] Bioconductor 3.18 (R 4.3.1)
# fansi                  1.0.6     2023-12-08 [1] CRAN (R 4.3.1)
# forcats              * 1.0.0     2023-01-29 [1] CRAN (R 4.3.1)
# generics               0.1.3     2022-07-05 [1] CRAN (R 4.3.1)
# GenomeInfoDb         * 1.38.8    2024-03-15 [1] Bioconductor 3.18 (R 4.3.1)
# GenomeInfoDbData       1.2.11    2024-11-10 [1] Bioconductor
# GenomicRanges        * 1.54.1    2023-10-29 [1] Bioconductor
# ggplot2              * 3.5.1     2024-04-23 [1] CRAN (R 4.3.1)
# ggrepel                0.9.5     2024-01-10 [1] CRAN (R 4.3.1)
# glue                   1.7.0     2024-01-09 [1] CRAN (R 4.3.1)
# gtable                 0.3.4     2023-08-21 [1] CRAN (R 4.3.1)
# hms                    1.1.3     2023-03-21 [1] CRAN (R 4.3.1)
# invgamma               1.1       2017-05-07 [1] CRAN (R 4.3.1)
# IRanges              * 2.36.0    2023-10-24 [1] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [1] CRAN (R 4.3.1)
# lattice                0.22-6    2024-03-20 [1] CRAN (R 4.3.1)
# lifecycle              1.0.4     2023-11-07 [1] CRAN (R 4.3.1)
# limma                * 3.58.1    2023-10-31 [1] Bioconductor
# locfit                 1.5-9.9   2024-03-01 [1] CRAN (R 4.3.1)
# lubridate            * 1.9.3     2023-09-27 [1] CRAN (R 4.3.1)
# magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.3.1)
# mashr                * 0.2.79    2023-10-18 [1] CRAN (R 4.3.1)
# Matrix                 1.6-5     2024-01-11 [1] CRAN (R 4.3.1)
# MatrixGenerics       * 1.14.0    2023-10-24 [1] Bioconductor
# matrixStats          * 1.2.0     2023-12-11 [1] CRAN (R 4.3.1)
# mixsqp                 0.3-54    2023-12-20 [1] CRAN (R 4.3.1)
# munsell                0.5.1     2024-04-01 [1] CRAN (R 4.3.1)
# mvtnorm                1.2-4     2023-11-27 [1] CRAN (R 4.3.1)
# pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.3.1)
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.1)
# plyr                   1.8.9     2023-10-02 [1] CRAN (R 4.3.1)
# purrr                * 1.0.2     2023-08-10 [1] CRAN (R 4.3.1)
# R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.1)
# Rcpp                   1.0.12    2024-01-09 [1] CRAN (R 4.3.1)
# RCurl                  1.98-1.14 2024-01-09 [1] CRAN (R 4.3.1)
# readr                * 2.1.5     2024-01-10 [1] CRAN (R 4.3.1)
# reshape2             * 1.4.4     2020-04-09 [1] CRAN (R 4.3.1)
# rlang                * 1.1.3     2024-01-10 [1] CRAN (R 4.3.1)
# rmeta                  3.0       2018-03-20 [1] CRAN (R 4.3.1)
# rstudioapi             0.16.0    2024-03-24 [1] CRAN (R 4.3.1)
# S4Arrays               1.2.1     2024-03-04 [1] Bioconductor 3.18 (R 4.3.1)
# S4Vectors            * 0.40.2    2023-11-23 [1] Bioconductor 3.18 (R 4.3.1)
# scales                 1.3.0     2023-11-28 [1] CRAN (R 4.3.1)
# sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.3.1)
# SparseArray            1.2.4     2024-02-11 [1] Bioconductor 3.18 (R 4.3.1)
# SQUAREM                2021.1    2021-01-13 [1] CRAN (R 4.3.1)
# statmod                1.5.0     2023-01-06 [1] CRAN (R 4.3.1)
# stringi                1.8.3     2023-12-11 [1] CRAN (R 4.3.1)
# stringr              * 1.5.1     2023-11-14 [1] CRAN (R 4.3.1)
# SummarizedExperiment * 1.32.0    2023-10-24 [1] Bioconductor
# tibble               * 3.2.1     2023-03-20 [1] CRAN (R 4.3.1)
# tidyr                * 1.3.1     2024-01-24 [1] CRAN (R 4.3.1)
# tidyselect             1.2.1     2024-03-11 [1] CRAN (R 4.3.1)
# tidyverse            * 2.0.0     2023-02-22 [1] CRAN (R 4.3.1)
# timechange             0.3.0     2024-01-18 [1] CRAN (R 4.3.1)
# truncnorm              1.0-9     2023-03-20 [1] CRAN (R 4.3.1)
# tzdb                   0.4.0     2023-05-12 [1] CRAN (R 4.3.1)
# utf8                   1.2.4     2023-10-22 [1] CRAN (R 4.3.1)
# vctrs                  0.6.5     2023-12-01 [1] CRAN (R 4.3.1)
# withr                  3.0.0     2024-01-16 [1] CRAN (R 4.3.1)
# XVector                0.42.0    2023-10-24 [1] Bioconductor
# zlibbioc               1.48.2    2024-03-13 [1] Bioconductor 3.18 (R 4.3.1)
# 
# [1] /opt/view/rlib/R/library
# [2] /software/hgi/installs/softpack/rstudio/.spack/linux-ubuntu22.04-x86_64_v3/gcc-11.4.0/r-4.3.1-bfwldrk76z6f52upk47zepliekn7ayqz/rlib/R/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────


                  