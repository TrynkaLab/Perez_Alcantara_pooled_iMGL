
library(tidyverse)
library(rlang)
library(lmerTest)
library(cowplot)
library(ggrepel)
library(sessioninfo)


################################################################################
##            4. Burden tests for line proliferation efficiency
################################################################################

#-------------------------------------------------------------------------------
#                  4.1 Perform burden tests for proliferation
#-------------------------------------------------------------------------------
#  Code to perform burden tests to asssess the effect of deleterious, protein-
#  truncating, and synonymous variant burdens on proliferation efficiency.
#  - Note: analysis based on non-prolif cells.
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

## Set working dir
setwd("/lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_Daianna/")

## Define input, output, and plot dirs
inputdir = paste(getwd(), "input_data", "04_Burden_test_proliferation", "01_Burden_tests", sep = "/")
outdir = paste(getwd(), "output_data", "04_Burden_test_proliferation", "01_Burden_tests", sep = "/")
plotdir = paste(getwd(), "plots", "04_Burden_test_proliferation", "01_Burden_tests", sep = "/")
dir.create(inputdir, recursive = T)
dir.create(outdir, recursive = T)
dir.create(plotdir, recursive = T)

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#                  4.1.0 Explore and process proliferation data
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

## Log-fractions per donor per pool
## - From iPSC -> young preMac:
line_prop_changes_premac_vs_iPSC <- as.data.frame(read_csv("/lustre/scratch123/hgi/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/data/results/1.2.scale_proliferation/line_prop_changes_premac_iPSC.csv"))
## - From young preMac -> old preMac:
line_prop_changes_old_vs_young_premac <- as.data.frame(read_csv("/lustre/scratch123/hgi/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/data/results/1.2.scale_proliferation/line_prop_changes_old_vs_young_premac.csv"))
## - From preMac (diff ages) -> microglia (in IFN/LPS/untreated):
line_prop_changes_microglia_vs_premac <- as.data.frame(read_csv("/lustre/scratch123/hgi/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/data/results/1.2.scale_proliferation/line_prop_changes_microglia_premac.csv"))


##################  iPSC -> young preMac proportion changes  ###################

## Line x pool (x batch for some lines)
dim(line_prop_changes_premac_vs_iPSC)
# [1] 402  11
dim(unique(line_prop_changes_premac_vs_iPSC[, c("line", "pool", "differentiation")]))
# [1] 402   3

## All for same comparison
table(line_prop_changes_premac_vs_iPSC$comparison)
# youngest_preMac_vs_iPSC 
# 402 

## 238 unique lines and 15 pools
length(unique(line_prop_changes_premac_vs_iPSC$line))
# [1] 238
length(unique(line_prop_changes_premac_vs_iPSC$pool))
# [1] 15

## Lines per pool
table(line_prop_changes_premac_vs_iPSC$pool)
# pool10 pool11 pool13 pool14 pool15 pool16 pool17  pool2  pool3  pool4  pool5  pool6  pool7  pool8  pool9 
#     55     23     19     23     16     24     32     35     20     24     19     23     26     22     41 

## Scaled log prop fractions
head(line_prop_changes_premac_vs_iPSC, 3)
#         line   pool differentiation log_fraction_mean prop_unadjusted_max_value prop_unadjusted_min_value              comparison preMAC_age type    sex
# 1 Arlene-003 pool15       IPMAR15-I         0.1341518                0.11758165                0.10243138 youngest_preMac_vs_iPSC         36  WGS Female
# 2 Bertha-004 pool15       IPMAR15-I         0.3976784                0.01699674                0.01197973 youngest_preMac_vs_iPSC         36  WGS Female
# 3  Cindy-005 pool15       IPMAR15-I        -0.5454670                0.09772526                0.05582467 youngest_preMac_vs_iPSC         36  WGS Female
#   scaled_log_fraction
# 1           0.2073989
# 2           0.6956558
# 3          -1.0517854


###############  young preMac -> old preMac proportion changes  ################

## Line x pool (x batch)
dim(line_prop_changes_old_vs_young_premac)
# [1] 355   9
dim(unique(line_prop_changes_old_vs_young_premac[, c("line", "pool", "differentiation")]))
# [1] 355   3

## All for same comparison
table(line_prop_changes_old_vs_young_premac$comparison)
# preMac_old_vs_young 
# 355 

## 220 unique lines and 15 pools
length(unique(line_prop_changes_old_vs_young_premac$line))
# [1] 220
length(unique(line_prop_changes_old_vs_young_premac$pool))
# [1] 15

## Lines per pool
table(line_prop_changes_old_vs_young_premac$pool)
# pool10 pool11 pool13 pool14 pool15 pool16 pool17  pool2  pool3  pool4  pool5  pool6  pool7  pool8  pool9 
#     27     22     19     23     16     24     31     33     19     23     19     23     23     22     31 

## Scaled log prop fractions
head(line_prop_changes_old_vs_young_premac, 3)
#     pool differentiation   line log_fraction_mean prop_unadjusted_max_value prop_unadjusted_min_value          comparison    sex scaled_log_fraction
# 1 pool10      hipsci10-I aizi_3        -1.4581901               0.007771699               0.003366767 preMac_old_vs_young Female          -0.3970394
# 2 pool10      hipsci10-I babk_2         0.6393642               0.070396477               0.037606358 preMac_old_vs_young Female           1.5077873
# 3 pool10      hipsci10-I civh_1         0.7855631               0.254489145               0.116272296 preMac_old_vs_young Female           1.6405531


########### (young - old) preMac -> microglia (in IFN/LPS/untreated) ###########

dim(line_prop_changes_microglia_vs_premac)
# [1] 2344   12

## Same comparison but microglia in the 3 treatments
table(line_prop_changes_microglia_vs_premac$comparison)
# microglia_vs_preMac 
# 2344 
table(line_prop_changes_microglia_vs_premac$treatment)
# IFNg       LPS untreated 
#  823       957       564 

## Diff preMac ages
table(line_prop_changes_microglia_vs_premac$preMAC_age)
# 28  35  36  39  40  42  43  46  47  49  50  54  57  60 
# 15 300 438 340  65  60 201 164  48 160 274 155 105  19 

## Line x pool (x batch) x treatment x preMac age
dim(unique(line_prop_changes_microglia_vs_premac[, c("line", "pool", "differentiation", "treatment", "preMAC_age")]))
# [1] 2344    5

## e.g. prolif of debk_9 in pool10, batch hipsci10-I from preMac 36 -> microglia IFN, etc ... 
subset(line_prop_changes_microglia_vs_premac, line  == "debk_9" & pool == "pool10" & differentiation == "hipsci10-I") 
#       pool preMAC_age differentiation   line treatment      type log_fraction_mean prop_unadjusted_max_value prop_unadjusted_min_value          comparison
# 19  pool10         36      hipsci10-I debk_9      IFNg scRNA-seq        1.23130655               0.008798442               0.004384436 microglia_vs_preMac
# 20  pool10         36      hipsci10-I debk_9       LPS scRNA-seq        0.89415801               0.006280353               0.004384436 microglia_vs_preMac
# 125 pool10         43      hipsci10-I debk_9       LPS       WGS       -0.05135493               0.004173264               0.004055267 microglia_vs_preMac
# 126 pool10         43      hipsci10-I debk_9 untreated       WGS       -0.26619979               0.004173264               0.003622356 microglia_vs_preMac
# 178 pool10         50      hipsci10-I debk_9      IFNg       WGS        0.85555717               0.004582224               0.003130386 microglia_vs_preMac
# 179 pool10         50      hipsci10-I debk_9       LPS       WGS        0.57440237               0.003904272               0.003130386 microglia_vs_preMac
# 180 pool10         50      hipsci10-I debk_9 untreated       WGS        0.61305368               0.003986564               0.003130386 microglia_vs_preMac
#        sex scaled_log_fraction
# 19  Female           0.1455149
# 20  Female           0.2157142
# 125 Female           0.8921817
# 126 Female           0.6631039
# 178 Female           1.3773139
# 179 Female           1.1332734
# 180 Female           1.3574247


## Add genotype PC1-10
genotype_PCs = read.table("../OTAR2065_sc_eQTL/data/genotype/plink_genotypes/all_pools.no_outliers.genotype.MAF05.eigenvec")
genotype_PCs$V1 <- NULL
colnames(genotype_PCs) <- c("line", paste0("PC", 1:(dim(genotype_PCs)[2]-1)))
genotype_first_PCs <- genotype_PCs[, 1:11]

line_prop_changes_premac_vs_iPSC <- inner_join(line_prop_changes_premac_vs_iPSC, genotype_first_PCs, by = "line")
save(line_prop_changes_premac_vs_iPSC, file = paste0(outdir, "/line_prop_changes_premac_vs_iPSC.Rdata"))
line_prop_changes_old_vs_young_premac <- inner_join(line_prop_changes_old_vs_young_premac, genotype_first_PCs, by = "line")
save(line_prop_changes_old_vs_young_premac, file = paste0(outdir, "/line_prop_changes_old_vs_young_premac.Rdata"))
line_prop_changes_microglia_vs_premac <- inner_join(line_prop_changes_microglia_vs_premac, genotype_first_PCs, by = "line")
save(line_prop_changes_microglia_vs_premac, file = paste0(outdir, "/line_prop_changes_microglia_vs_premac.Rdata"))



## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#                    4.1.1 Compute genetic burdens per line  
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
##                  4.1.1.1 Global genetic burdens per line
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

## Num of del/ptv/syn variants per gene per donor 
mutBurden_del <- read_rds("/lustre/scratch123/hgi/projects/otar2065/resources/Puigdevall_Neuroseq_efficiency_2023/mutBurdenTabs/mutBurden_del.RDS") %>% as.data.frame()
mutBurden_ptv <- read_rds("/lustre/scratch123/hgi/projects/otar2065/resources/Puigdevall_Neuroseq_efficiency_2023/mutBurdenTabs/mutBurden_ptv.RDS") %>% as.data.frame()
mutBurden_syn <- read_rds("/lustre/scratch123/hgi/projects/otar2065/resources/Puigdevall_Neuroseq_efficiency_2023/mutBurdenTabs/mutBurden_synonymous.RDS") %>% as.data.frame()

#######################  Burden of deleterious variants  #######################
dim(mutBurden_del)
# [1] 19653 genes,   832 lines

mutBurden_del[1:5, 1:5]
#            HPSI0114i-bezi_1 HPSI0114i-bezi_3 HPSI0114i-eipl_1 HPSI0114i-eipl_3 HPSI0114i-fikt_3
# OR4F5                     1                1                1                1                0
# SAMD11                    1                2                1                1                1
# AL645608.1                0                0                0                0                0
# NOC2L                     1                1                1                1                1
# KLHL17                    0                0                0                0                0

## Unique genes and lines
which(duplicated(colnames(mutBurden_del)))
# integer(0)
which(duplicated(rownames(mutBurden_del)))
# integer(0)

## Global burden (across all genes)
mutBurden_del <- rbind(mutBurden_del, "Global_burden" = apply(mutBurden_del, 2, function(c){sum(as.numeric(c))}))

## Cell line
mutBurden_del <- rbind(mutBurden_del, "line" = gsub(".*-", "", colnames(mutBurden_del)))
which(duplicated(mutBurden_del["line", ]))
# integer(0)

####################  Burden of protein-truncating variants  ###################
dim(mutBurden_ptv)
# [1] 19653 genes,   832 lines

mutBurden_ptv[1:5, 1:5]
#            HPSI0114i-bezi_1 HPSI0114i-bezi_3 HPSI0114i-eipl_1 HPSI0114i-eipl_3 HPSI0114i-fikt_3
# OR4F5                     0                0                0                0                0
# SAMD11                    0                1                0                0                0
# AL645608.1                0                0                0                0                0
# NOC2L                     0                0                0                0                0
# KLHL17                    0                0                0                0                0

## Global burden 
mutBurden_ptv <- rbind(mutBurden_ptv, "Global_burden" = apply(mutBurden_ptv, 2, function(c){sum(as.numeric(c))}))

## Cell line
mutBurden_ptv <- rbind(mutBurden_ptv, "line" = gsub(".*-", "", colnames(mutBurden_ptv)))
which(duplicated(mutBurden_ptv["line", ]))
# integer(0)

########################  Burden of synonymous variants  #######################
dim(mutBurden_syn)
# [1] 19653 genes,   832 lines

mutBurden_syn[1:5, 1:5]
#           HPSI0114i-bezi_1 HPSI0114i-bezi_3 HPSI0114i-eipl_1 HPSI0114i-eipl_3 HPSI0114i-fikt_3
# OR4F5                     0                0                0                0                0
# SAMD11                    0                0                0                0                0
# AL645608.1                0                0                0                0                0
# NOC2L                     3                3                2                2                3
# KLHL17                    2                2                2                2                1

## Global burden 
mutBurden_syn <- rbind(mutBurden_syn, "Global_burden" = apply(mutBurden_syn, 2, function(c){sum(as.numeric(c))}))

## Cell line
mutBurden_syn <- rbind(mutBurden_syn, "line" = gsub(".*-", "", colnames(mutBurden_syn)))
which(duplicated(mutBurden_syn["line", ]))
# integer(0)


## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
##            4.1.1.2 Genetic burdens per line in SKAT-O genes*
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# * significant genes in SKAT-O tests for deleterious variants

## Results
deleterious_SKATO_results <- read_csv(paste0("/lustre/scratch123/hgi/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/data/results/2.1.check_SKAT_WES/",
                                              "deleterious_SKATO_scaled_centered_prop_pvals.csv"))

## Signif genes 
deleterious_SKATO_results_signif <- subset(deleterious_SKATO_results, resampling_pval < 0.05)

## Divide by comparison
deleterious_SKATO_signif_premac_vs_iPSC <- subset(deleterious_SKATO_results_signif, comparison == "line_prop_changes_premac_iPSC")
deleterious_SKATO_signif_old_vs_young_premac <- subset(deleterious_SKATO_results_signif, comparison == "line_prop_changes_old_vs_young_premac")
deleterious_SKATO_signif_microglia_vs_premac <- subset(deleterious_SKATO_results_signif, comparison == "line_prop_changes_microglia_premac")

## Genes in mutBurden data
deleterious_SKATO_signif_premac_vs_iPSC <- subset(deleterious_SKATO_signif_premac_vs_iPSC, gene_name %in% rownames(mutBurden_del))
# dim: 180   6
deleterious_SKATO_signif_old_vs_young_premac <- subset(deleterious_SKATO_signif_old_vs_young_premac, gene_name %in% rownames(mutBurden_del))
# dim: 233   6
deleterious_SKATO_signif_microglia_vs_premac <- subset(deleterious_SKATO_signif_microglia_vs_premac, gene_name %in% rownames(mutBurden_del))
# dim: 263   6

#######################  Burden of deleterious variants  #######################

# - In signif genes for iPSC -> preMac
mutBurden_del <- rbind(mutBurden_del, "SKAT0_del_premac_vs_iPSC_burden" = apply(mutBurden_del[deleterious_SKATO_signif_premac_vs_iPSC$gene_name, ], 2, function(c){sum(as.numeric(c))}))
# - In signif genes for young -> old preMac
mutBurden_del <- rbind(mutBurden_del, "SKAT0_del_old_vs_young_premac_burden" = apply(mutBurden_del[deleterious_SKATO_signif_old_vs_young_premac$gene_name, ], 2, function(c){sum(as.numeric(c))}))
# - In signif genes for preMac -> microglia
mutBurden_del <- rbind(mutBurden_del, "SKAT0_del_microglia_vs_premac_burden" = apply(mutBurden_del[deleterious_SKATO_signif_microglia_vs_premac$gene_name, ], 2, function(c){sum(as.numeric(c))}))

####################  Burden of protein-truncating variants  ###################

# - In signif genes for iPSC -> preMac
mutBurden_ptv <- rbind(mutBurden_ptv, "SKAT0_del_premac_vs_iPSC_burden" = apply(mutBurden_ptv[deleterious_SKATO_signif_premac_vs_iPSC$gene_name, ], 2, function(c){sum(as.numeric(c))}))
# - In signif genes for young -> old preMac
mutBurden_ptv <- rbind(mutBurden_ptv, "SKAT0_del_old_vs_young_premac_burden" = apply(mutBurden_ptv[deleterious_SKATO_signif_old_vs_young_premac$gene_name, ], 2, function(c){sum(as.numeric(c))}))
# - In signif genes for preMac -> microglia
mutBurden_ptv <- rbind(mutBurden_ptv, "SKAT0_del_microglia_vs_premac_burden" = apply(mutBurden_ptv[deleterious_SKATO_signif_microglia_vs_premac$gene_name, ], 2, function(c){sum(as.numeric(c))}))

########################  Burden of synonymous variants  #######################

# - In signif genes for iPSC -> preMac
mutBurden_syn <- rbind(mutBurden_syn, "SKAT0_del_premac_vs_iPSC_burden" = apply(mutBurden_syn[deleterious_SKATO_signif_premac_vs_iPSC$gene_name, ], 2, function(c){sum(as.numeric(c))}))
# - In signif genes for young -> old preMac
mutBurden_syn <- rbind(mutBurden_syn, "SKAT0_del_old_vs_young_premac_burden" = apply(mutBurden_syn[deleterious_SKATO_signif_old_vs_young_premac$gene_name, ], 2, function(c){sum(as.numeric(c))}))
# - In signif genes for preMac -> microglia
mutBurden_syn <- rbind(mutBurden_syn, "SKAT0_del_microglia_vs_premac_burden" = apply(mutBurden_syn[deleterious_SKATO_signif_microglia_vs_premac$gene_name, ], 2, function(c){sum(as.numeric(c))}))


## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
##       4.1.1.3 Genetic burdens per line in signif genes* of burden tests
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# * significant genes of burden tests for deleterious variants

## Signif genes from burden test 
del_burden_scaled_centered_prop_pvals <- as.data.frame(read_csv("/lustre/scratch123/hgi/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/data/results/2.2.rare_vars_vs_prolif/del_burden_scaled_centered_prop_pvals.csv"))
deleterious_Burden_results_signif <- subset(del_burden_scaled_centered_prop_pvals, p_Bonf < 0.05)

## Divide by comparison
deleterious_Burden_signif_premac_vs_iPSC <- subset(deleterious_Burden_results_signif, comparison == "line_prop_changes_premac_iPSC")
deleterious_Burden_signif_old_vs_young_premac <- subset(deleterious_Burden_results_signif, comparison == "line_prop_changes_old_vs_young_premac")
deleterious_Burden_signif_microglia_vs_premac <- subset(deleterious_Burden_results_signif, comparison == "line_prop_changes_microglia_premac")

## Genes in mutBurden data
deleterious_Burden_signif_premac_vs_iPSC <- subset(deleterious_Burden_signif_premac_vs_iPSC, gene_name %in% rownames(mutBurden_del))
# dim: 3  5
deleterious_Burden_signif_old_vs_young_premac <- subset(deleterious_Burden_signif_old_vs_young_premac, gene_name %in% rownames(mutBurden_del))
# dim: 1  5
deleterious_Burden_signif_microglia_vs_premac <- subset(deleterious_Burden_signif_microglia_vs_premac, gene_name %in% rownames(mutBurden_del))
# dim: 18  5

## Further divide by gene effect direction
deleterious_Burden_signif_premac_vs_iPSC_pos <- subset(deleterious_Burden_signif_premac_vs_iPSC, coef>0)
deleterious_Burden_signif_premac_vs_iPSC_neg <- subset(deleterious_Burden_signif_premac_vs_iPSC, coef<0)

deleterious_Burden_signif_old_vs_young_premac_pos <- subset(deleterious_Burden_signif_old_vs_young_premac, coef>0)
# deleterious_Burden_signif_old_vs_young_premac_neg <- subset(deleterious_Burden_signif_old_vs_young_premac, coef<0) 0 genes

deleterious_Burden_signif_microglia_vs_premac_pos <- subset(deleterious_Burden_signif_microglia_vs_premac, coef>0)
deleterious_Burden_signif_microglia_vs_premac_neg <- subset(deleterious_Burden_signif_microglia_vs_premac, coef<0)


#######################  Burden of deleterious variants  #######################

# - In signif genes for iPSC -> preMac
mutBurden_del <- rbind(mutBurden_del, "Burden_del_premac_vs_iPSC_burden" = apply(mutBurden_del[deleterious_Burden_signif_premac_vs_iPSC$gene_name, ], 2, function(c){sum(as.numeric(c))}))
# - In signif genes for young -> old preMac
mutBurden_del <- rbind(mutBurden_del, "Burden_del_old_vs_young_premac_burden" = apply(mutBurden_del[deleterious_Burden_signif_old_vs_young_premac$gene_name, ], 2, function(c){sum(as.numeric(c))}))
# - In signif genes for preMac -> microglia
mutBurden_del <- rbind(mutBurden_del, "Burden_del_microglia_vs_premac_burden" = apply(mutBurden_del[deleterious_Burden_signif_microglia_vs_premac$gene_name, ], 2, function(c){sum(as.numeric(c))}))

#-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
#  -  In signif genes for iPSC -> preMac with + beta
mutBurden_del <- rbind(mutBurden_del, "Burden_del_premac_vs_iPSC_burden_pos" = apply(mutBurden_del[deleterious_Burden_signif_premac_vs_iPSC_pos$gene_name, ], 2, function(c){sum(as.numeric(c))}))
#  -  In signif genes for iPSC -> preMac with - beta
mutBurden_del <- rbind(mutBurden_del, "Burden_del_premac_vs_iPSC_burden_neg" = apply(mutBurden_del[deleterious_Burden_signif_premac_vs_iPSC_neg$gene_name, ], 2, function(c){sum(as.numeric(c))}))

#-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
#  -  In signif genes for young -> old preMac with + beta
mutBurden_del <- rbind(mutBurden_del, "Burden_del_old_vs_young_premac_burden_pos" = apply(mutBurden_del[deleterious_Burden_signif_old_vs_young_premac_pos$gene_name, ], 2, function(c){sum(as.numeric(c))}))

#-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
#  -  In signif genes for preMac -> microglia with + beta
mutBurden_del <- rbind(mutBurden_del, "Burden_del_microglia_vs_premac_burden_pos" = apply(mutBurden_del[deleterious_Burden_signif_microglia_vs_premac_pos$gene_name, ], 2, function(c){sum(as.numeric(c))}))
#  -  In signif genes for preMac -> microglia with - beta
mutBurden_del <- rbind(mutBurden_del, "Burden_del_microglia_vs_premac_burden_neg" = apply(mutBurden_del[deleterious_Burden_signif_microglia_vs_premac_neg$gene_name, ], 2, function(c){sum(as.numeric(c))}))


####################  Burden of protein-truncating variants  ###################

# - In signif genes for iPSC -> preMac
mutBurden_ptv <- rbind(mutBurden_ptv, "Burden_del_premac_vs_iPSC_burden" = apply(mutBurden_ptv[deleterious_Burden_signif_premac_vs_iPSC$gene_name, ], 2, function(c){sum(as.numeric(c))}))
# - In signif genes for young -> old preMac
mutBurden_ptv <- rbind(mutBurden_ptv, "Burden_del_old_vs_young_premac_burden" = apply(mutBurden_ptv[deleterious_Burden_signif_old_vs_young_premac$gene_name, ], 2, function(c){sum(as.numeric(c))}))
# - In signif genes for preMac -> microglia
mutBurden_ptv <- rbind(mutBurden_ptv, "Burden_del_microglia_vs_premac_burden" = apply(mutBurden_ptv[deleterious_Burden_signif_microglia_vs_premac$gene_name, ], 2, function(c){sum(as.numeric(c))}))

#-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
#  -  In signif genes for iPSC -> preMac with + beta
mutBurden_ptv <- rbind(mutBurden_ptv, "Burden_del_premac_vs_iPSC_burden_pos" = apply(mutBurden_ptv[deleterious_Burden_signif_premac_vs_iPSC_pos$gene_name, ], 2, function(c){sum(as.numeric(c))}))
#  -  In signif genes for iPSC -> preMac with - beta
mutBurden_ptv <- rbind(mutBurden_ptv, "Burden_del_premac_vs_iPSC_burden_neg" = apply(mutBurden_ptv[deleterious_Burden_signif_premac_vs_iPSC_neg$gene_name, ], 2, function(c){sum(as.numeric(c))}))

#-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
#  -  In signif genes for young -> old preMac with + beta
mutBurden_ptv <- rbind(mutBurden_ptv, "Burden_del_old_vs_young_premac_burden_pos" = apply(mutBurden_ptv[deleterious_Burden_signif_old_vs_young_premac_pos$gene_name, ], 2, function(c){sum(as.numeric(c))}))

#-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
#  -  In signif genes for preMac -> microglia with + beta
mutBurden_ptv <- rbind(mutBurden_ptv, "Burden_del_microglia_vs_premac_burden_pos" = apply(mutBurden_ptv[deleterious_Burden_signif_microglia_vs_premac_pos$gene_name, ], 2, function(c){sum(as.numeric(c))}))
#  -  In signif genes for preMac -> microglia with - beta
mutBurden_ptv <- rbind(mutBurden_ptv, "Burden_del_microglia_vs_premac_burden_neg" = apply(mutBurden_ptv[deleterious_Burden_signif_microglia_vs_premac_neg$gene_name, ], 2, function(c){sum(as.numeric(c))}))


########################  Burden of synonymous variants  #######################

# - In signif genes for iPSC -> preMac
mutBurden_syn <- rbind(mutBurden_syn, "Burden_del_premac_vs_iPSC_burden" = apply(mutBurden_syn[deleterious_Burden_signif_premac_vs_iPSC$gene_name, ], 2, function(c){sum(as.numeric(c))}))
# - In signif genes for young -> old preMac
mutBurden_syn <- rbind(mutBurden_syn, "Burden_del_old_vs_young_premac_burden" = apply(mutBurden_syn[deleterious_Burden_signif_old_vs_young_premac$gene_name, ], 2, function(c){sum(as.numeric(c))}))
# - In signif genes for preMac -> microglia
mutBurden_syn <- rbind(mutBurden_syn, "Burden_del_microglia_vs_premac_burden" = apply(mutBurden_syn[deleterious_Burden_signif_microglia_vs_premac$gene_name, ], 2, function(c){sum(as.numeric(c))}))

#-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
#  -  In signif genes for iPSC -> preMac with + beta
mutBurden_syn <- rbind(mutBurden_syn, "Burden_del_premac_vs_iPSC_burden_pos" = apply(mutBurden_syn[deleterious_Burden_signif_premac_vs_iPSC_pos$gene_name, ], 2, function(c){sum(as.numeric(c))}))
#  -  In signif genes for iPSC -> preMac with - beta
mutBurden_syn <- rbind(mutBurden_syn, "Burden_del_premac_vs_iPSC_burden_neg" = apply(mutBurden_syn[deleterious_Burden_signif_premac_vs_iPSC_neg$gene_name, ], 2, function(c){sum(as.numeric(c))}))

#-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
#  -  In signif genes for young -> old preMac with + beta
mutBurden_syn <- rbind(mutBurden_syn, "Burden_del_old_vs_young_premac_burden_pos" = apply(mutBurden_syn[deleterious_Burden_signif_old_vs_young_premac_pos$gene_name, ], 2, function(c){sum(as.numeric(c))}))

#-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
#  -  In signif genes for preMac -> microglia with + beta
mutBurden_syn <- rbind(mutBurden_syn, "Burden_del_microglia_vs_premac_burden_pos" = apply(mutBurden_syn[deleterious_Burden_signif_microglia_vs_premac_pos$gene_name, ], 2, function(c){sum(as.numeric(c))}))
#  -  In signif genes for preMac -> microglia with - beta
mutBurden_syn <- rbind(mutBurden_syn, "Burden_del_microglia_vs_premac_burden_neg" = apply(mutBurden_syn[deleterious_Burden_signif_microglia_vs_premac_neg$gene_name, ], 2, function(c){sum(as.numeric(c))}))



## All del/ptv/syn burdens per line
all_Del_Burdens_per_line <- as.data.frame(t(mutBurden_del[19654:19666, ]))
all_Ptv_Burdens_per_line <- as.data.frame(t(mutBurden_ptv[19654:19666, ]))
all_Syn_Burdens_per_line <- as.data.frame(t(mutBurden_syn[19654:19666, ]))

save(all_Del_Burdens_per_line, file = paste0(outdir, "/all_Del_Burdens_per_line.Rdata"))
save(all_Ptv_Burdens_per_line, file = paste0(outdir, "/all_Ptv_Burdens_per_line.Rdata"))
save(all_Syn_Burdens_per_line, file = paste0(outdir, "/all_Syn_Burdens_per_line.Rdata"))



## Explore burden tests results for del/ptv/syn variants
del_burden_results <- as.data.frame(read_csv("/lustre/scratch123/hgi/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/data/results/2.2.rare_vars_vs_prolif/del_burden_scaled_centered_prop_pvals.csv"))
ptv_burden_results<- as.data.frame(read_csv("/lustre/scratch123/hgi/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/data/results/2.2.rare_vars_vs_prolif/ptv_burden_scaled_centered_prop_pvals.csv"))
syn_burden_results <- as.data.frame(read_csv("/lustre/scratch123/hgi/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/data/results/2.2.rare_vars_vs_prolif/syn_burden_scaled_centered_prop_pvals.csv"))

## Add col for variant type and merge
del_burden_results$variant_type <- "Deleterious"
ptv_burden_results$variant_type <- "Missense non-deleterious"
syn_burden_results$variant_type <- "Synonymous"

burden_results <- rbind(del_burden_results, ptv_burden_results, syn_burden_results)

## Row labels
comparison_labs <- c("iPSC to macrophage prec.", "young to old macrophage prec.", "macrophage prec. to microglia")
names(comparison_labs) <- c("line_prop_changes_premac_iPSC", "line_prop_changes_old_vs_young_premac", "line_prop_changes_microglia_premac")
burden_results$comparison <- factor(burden_results$comparison, levels = names(comparison_labs))
  
## Signif genes
burden_results$signif_dir <- apply(burden_results, 1, function(x){if(as.numeric(x["p_Bonf"])<0.05){ if(as.numeric(x["coef"]) < 0){"Decreases proliferation"} else{"Increases proliferation"} }
                                                                      else{"n.s."}})
burden_results$gene_label <- apply(burden_results, 1, function(x){ if(as.numeric(x["p_Bonf"])<0.05){x["gene_name"]} else{NA} })

set.seed(06122024)
ggplot(burden_results, aes(x = coef, y = -log10(p_Bonf), color = signif_dir, alpha = signif_dir, label = gene_label)) +
  geom_point(size = 1, show.legend = TRUE) +
  scale_color_manual(values = c("Decreases proliferation" = "blue", "Increases proliferation" = "red", "n.s." = "gray")) + 
  scale_alpha_manual(values = c("Decreases proliferation" = 1, "Increases proliferation" = 1, "n.s." = 0.5)) + 
  geom_text_repel(size = 2, show.legend = FALSE, min.segment.length = 0, color = "black", box.padding = 0.09) +
  facet_grid(rows = vars(comparison), cols = vars(variant_type), 
             labeller = labeller(comparison = comparison_labs)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey20", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey20", linewidth = 0.4) +
  labs(x = "Effect size", y = "-log10(adj. p)", color = "Effect in proliferation") + 
  guides(alpha = "none") +
  theme_void() +
  theme(legend.key.height = unit(0.25,"cm"),
        axis.title = element_text(size = (9)),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 8),
        strip.background.x = element_rect(fill = "gray90", linetype = "solid",
                                          color = "black", linewidth = 0.6),
        strip.background.y = element_rect(fill = "gray90", linetype = "solid",
                                          color = "black", linewidth = 0.6))

ggsave(filename = paste0(plotdir, "/volcano_plots_gene_burdens.pdf"), height = 6.5, width = 9)

## Pretier plot for all burdens
ggplot(burden_results, aes(x = coef, y = -log10(p_Bonf), color = signif_dir, alpha = signif_dir, label = gene_label)) +
  geom_point(size = 1, show.legend = TRUE) +
  scale_color_manual(values = c("Decreases proliferation" = "blue", "Increases proliferation" = "red", "n.s." = "gray")) + 
  scale_alpha_manual(values = c("Decreases proliferation" = 1, "Increases proliferation" = 1, "n.s." = 0.5)) + 
  geom_text_repel(size = 2.5, show.legend = FALSE, min.segment.length = 0, color = "black", box.padding = 0.2) +
  facet_grid(rows = vars(comparison), cols = vars(variant_type), 
             labeller = labeller(comparison = comparison_labs), axes = "all", axis.labels = "margins") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey20", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey20", linewidth = 0.4) +
  labs(x = "Fold change", y = expression(-log[10](P)), color = "Effect in proliferation") + 
  guides(color = "none", alpha = "none") +
  theme_classic() +
  theme(axis.title = element_text(size = 11),
        strip.text.x = element_text(size = 10.5),
        strip.text.y = element_text(size = 7),
        strip.background.x = element_blank(),
        strip.background.y = element_blank())

ggsave(filename = paste0(plotdir, "/volcano_plots_gene_burdens_prettier.pdf"), height = 5.5, width = 6.4)

## Prettier plot for del burden only
burden_results_del <- subset(burden_results, variant_type == "Deleterious")
ggplot(burden_results_del, aes(x = coef, y = -log10(p_Bonf), color = signif_dir, alpha = signif_dir, label = gene_label)) +
  geom_point(size = 1, show.legend = TRUE) +
  scale_color_manual(values = c("Decreases proliferation" = "blue", "Increases proliferation" = "red", "n.s." = "gray")) + 
  scale_alpha_manual(values = c("Decreases proliferation" = 1, "Increases proliferation" = 1, "n.s." = 0.5)) + 
  geom_text_repel(size = 2.5, show.legend = FALSE, min.segment.length = 0, color = "black", box.padding = 0.2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey20", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey20", linewidth = 0.4) +
  facet_grid(rows = vars(comparison), 
             labeller = labeller(comparison = comparison_labs), axes = "all_x") +
  labs(title = "Deleterious", x = "Fold change", y = expression(-log[10](P))) + 
  guides(color = "none", alpha = "none") +
  theme_classic() +
  # ggpubr::theme_pubr() +
  theme(plot.title = element_text(hjust = 0.5, size = 11.5) , 
        axis.title = element_text(size = (11)),
        strip.text.y = element_text(size = 7),
        strip.background.x = element_blank(),
        strip.background.y = element_blank())

ggsave(filename = paste0(plotdir, "/volcano_plots_gene_Del_burdens_prettier.pdf"), height = 5.5, width = 2.5)


## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
##              4.1.1.4 Add line burdens to proliferation data
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

##################  iPSC -> young preMac proportion changes  ###################

#~~~~~~~~~~~~~~~~~~~~~~ Burdens of deleterious variants ~~~~~~~~~~~~~~~~~~~~~~~~

## A. Global (across all genes)
line_prop_changes_premac_vs_iPSC_burdens <- inner_join(line_prop_changes_premac_vs_iPSC,  all_Del_Burdens_per_line[, c("line", "Global_burden")], by = "line")
colnames(line_prop_changes_premac_vs_iPSC_burdens)[dim(line_prop_changes_premac_vs_iPSC_burdens)[2]] <- "Global_Burden_Del"
dim(line_prop_changes_premac_vs_iPSC_burdens)
# [1] 339  22

## B. In SKAT-O genes associated with iPSC -> preMac based on their del variants
line_prop_changes_premac_vs_iPSC_burdens <- inner_join(line_prop_changes_premac_vs_iPSC_burdens,  all_Del_Burdens_per_line[, c("line", "SKAT0_del_premac_vs_iPSC_burden")], by = "line")
colnames(line_prop_changes_premac_vs_iPSC_burdens)[dim(line_prop_changes_premac_vs_iPSC_burdens)[2]] <- "SKAT0_del_premac_vs_iPSC_Burden_Del"

## C. In genes associated with iPSC -> preMac based on their del burdens
line_prop_changes_premac_vs_iPSC_burdens <- inner_join(line_prop_changes_premac_vs_iPSC_burdens,  all_Del_Burdens_per_line[, c("line", "Burden_del_premac_vs_iPSC_burden")], by = "line")
colnames(line_prop_changes_premac_vs_iPSC_burdens)[dim(line_prop_changes_premac_vs_iPSC_burdens)[2]] <- "Burden_del_premac_vs_iPSC_Burden_Del"
##     i. Genes with + beta
line_prop_changes_premac_vs_iPSC_burdens <- inner_join(line_prop_changes_premac_vs_iPSC_burdens,  all_Del_Burdens_per_line[, c("line", "Burden_del_premac_vs_iPSC_burden_pos")], by = "line")
colnames(line_prop_changes_premac_vs_iPSC_burdens)[dim(line_prop_changes_premac_vs_iPSC_burdens)[2]] <- "Burden_del_premac_vs_iPSC_pos_Burden_Del"
##     ii. Genes with - beta
line_prop_changes_premac_vs_iPSC_burdens <- inner_join(line_prop_changes_premac_vs_iPSC_burdens,  all_Del_Burdens_per_line[, c("line", "Burden_del_premac_vs_iPSC_burden_neg")], by = "line")
colnames(line_prop_changes_premac_vs_iPSC_burdens)[dim(line_prop_changes_premac_vs_iPSC_burdens)[2]] <- "Burden_del_premac_vs_iPSC_neg_Burden_Del"

#~~~~~~~~~~~~~~~~~~~ Burdens of protein-truncating variants ~~~~~~~~~~~~~~~~~~~~

## A. Global (across all genes)
line_prop_changes_premac_vs_iPSC_burdens <- inner_join(line_prop_changes_premac_vs_iPSC_burdens,  all_Ptv_Burdens_per_line[, c("line", "Global_burden")], by = "line")
colnames(line_prop_changes_premac_vs_iPSC_burdens)[dim(line_prop_changes_premac_vs_iPSC_burdens)[2]] <- "Global_Burden_Ptv"
dim(line_prop_changes_premac_vs_iPSC_burdens)
# [1] 339  27

## B. In SKAT-O genes associated with iPSC -> preMac based on their del variants
line_prop_changes_premac_vs_iPSC_burdens <- inner_join(line_prop_changes_premac_vs_iPSC_burdens,  all_Ptv_Burdens_per_line[, c("line", "SKAT0_del_premac_vs_iPSC_burden")], by = "line")
colnames(line_prop_changes_premac_vs_iPSC_burdens)[dim(line_prop_changes_premac_vs_iPSC_burdens)[2]] <- "SKAT0_del_premac_vs_iPSC_Burden_Ptv"

## C. In genes associated with iPSC -> preMac based on their del burdens
line_prop_changes_premac_vs_iPSC_burdens <- inner_join(line_prop_changes_premac_vs_iPSC_burdens,  all_Ptv_Burdens_per_line[, c("line", "Burden_del_premac_vs_iPSC_burden")], by = "line")
colnames(line_prop_changes_premac_vs_iPSC_burdens)[dim(line_prop_changes_premac_vs_iPSC_burdens)[2]] <- "Burden_del_premac_vs_iPSC_Burden_Ptv"
##     i. Genes with + beta
line_prop_changes_premac_vs_iPSC_burdens <- inner_join(line_prop_changes_premac_vs_iPSC_burdens,  all_Ptv_Burdens_per_line[, c("line", "Burden_del_premac_vs_iPSC_burden_pos")], by = "line")
colnames(line_prop_changes_premac_vs_iPSC_burdens)[dim(line_prop_changes_premac_vs_iPSC_burdens)[2]] <- "Burden_del_premac_vs_iPSC_pos_Burden_Ptv"
##     ii. Genes with - beta
line_prop_changes_premac_vs_iPSC_burdens <- inner_join(line_prop_changes_premac_vs_iPSC_burdens,  all_Ptv_Burdens_per_line[, c("line", "Burden_del_premac_vs_iPSC_burden_neg")], by = "line")
colnames(line_prop_changes_premac_vs_iPSC_burdens)[dim(line_prop_changes_premac_vs_iPSC_burdens)[2]] <- "Burden_del_premac_vs_iPSC_neg_Burden_Ptv"

#~~~~~~~~~~~~~~~~~~~~~~~ Burdens of synonymous variants ~~~~~~~~~~~~~~~~~~~~~~~~

## A. Global (across all genes)
line_prop_changes_premac_vs_iPSC_burdens <- inner_join(line_prop_changes_premac_vs_iPSC_burdens,  all_Syn_Burdens_per_line[, c("line", "Global_burden")], by = "line")
colnames(line_prop_changes_premac_vs_iPSC_burdens)[dim(line_prop_changes_premac_vs_iPSC_burdens)[2]] <- "Global_Burden_Syn"
dim(line_prop_changes_premac_vs_iPSC_burdens)
# [1] 339  32

## B. In SKAT-O genes associated with iPSC -> preMac based on their del variants
line_prop_changes_premac_vs_iPSC_burdens <- inner_join(line_prop_changes_premac_vs_iPSC_burdens,  all_Syn_Burdens_per_line[, c("line", "SKAT0_del_premac_vs_iPSC_burden")], by = "line")
colnames(line_prop_changes_premac_vs_iPSC_burdens)[dim(line_prop_changes_premac_vs_iPSC_burdens)[2]] <- "SKAT0_del_premac_vs_iPSC_Burden_Syn"

## C. In genes associated with iPSC -> preMac based on their del burdens
line_prop_changes_premac_vs_iPSC_burdens <- inner_join(line_prop_changes_premac_vs_iPSC_burdens,  all_Syn_Burdens_per_line[, c("line", "Burden_del_premac_vs_iPSC_burden")], by = "line")
colnames(line_prop_changes_premac_vs_iPSC_burdens)[dim(line_prop_changes_premac_vs_iPSC_burdens)[2]] <- "Burden_del_premac_vs_iPSC_Burden_Syn"
##     i. Genes with + beta
line_prop_changes_premac_vs_iPSC_burdens <- inner_join(line_prop_changes_premac_vs_iPSC_burdens,  all_Syn_Burdens_per_line[, c("line", "Burden_del_premac_vs_iPSC_burden_pos")], by = "line")
colnames(line_prop_changes_premac_vs_iPSC_burdens)[dim(line_prop_changes_premac_vs_iPSC_burdens)[2]] <- "Burden_del_premac_vs_iPSC_pos_Burden_Syn"
##     ii. Genes with - beta
line_prop_changes_premac_vs_iPSC_burdens <- inner_join(line_prop_changes_premac_vs_iPSC_burdens,  all_Syn_Burdens_per_line[, c("line", "Burden_del_premac_vs_iPSC_burden_neg")], by = "line")
colnames(line_prop_changes_premac_vs_iPSC_burdens)[dim(line_prop_changes_premac_vs_iPSC_burdens)[2]] <- "Burden_del_premac_vs_iPSC_neg_Burden_Syn"


##################  young -> old preMac proportion changes  ###################

#~~~~~~~~~~~~~~~~~~~~~~ Burdens of deleterious variants ~~~~~~~~~~~~~~~~~~~~~~~~

## A. Global (across all genes)
line_prop_changes_old_vs_young_premac_burdens <- inner_join(line_prop_changes_old_vs_young_premac,  all_Del_Burdens_per_line[, c("line", "Global_burden")], by = "line")
colnames(line_prop_changes_old_vs_young_premac_burdens)[dim(line_prop_changes_old_vs_young_premac_burdens)[2]] <- "Global_Burden_Del"
dim(line_prop_changes_old_vs_young_premac_burdens)
# [1] 294  20

## B. In SKAT-O genes associated with young -> old preMac based on their del variants
line_prop_changes_old_vs_young_premac_burdens <- inner_join(line_prop_changes_old_vs_young_premac_burdens,  all_Del_Burdens_per_line[, c("line", "SKAT0_del_old_vs_young_premac_burden")], by = "line")
colnames(line_prop_changes_old_vs_young_premac_burdens)[dim(line_prop_changes_old_vs_young_premac_burdens)[2]] <- "SKAT0_del_old_vs_young_premac_Burden_Del"

## C. In genes associated with young -> old preMac based on their del burdens
line_prop_changes_old_vs_young_premac_burdens <- inner_join(line_prop_changes_old_vs_young_premac_burdens,  all_Del_Burdens_per_line[, c("line", "Burden_del_old_vs_young_premac_burden")], by = "line")
colnames(line_prop_changes_old_vs_young_premac_burdens)[dim(line_prop_changes_old_vs_young_premac_burdens)[2]] <- "Burden_del_old_vs_young_premac_Burden_Del"
##     i. Genes with + beta
line_prop_changes_old_vs_young_premac_burdens <- inner_join(line_prop_changes_old_vs_young_premac_burdens,  all_Del_Burdens_per_line[, c("line", "Burden_del_old_vs_young_premac_burden_pos")], by = "line")
colnames(line_prop_changes_old_vs_young_premac_burdens)[dim(line_prop_changes_old_vs_young_premac_burdens)[2]] <- "Burden_del_old_vs_young_premac_pos_Burden_Del"

#~~~~~~~~~~~~~~~~~~~ Burdens of protein-truncating variants ~~~~~~~~~~~~~~~~~~~~

## A. Global (across all genes)
line_prop_changes_old_vs_young_premac_burdens <- inner_join(line_prop_changes_old_vs_young_premac_burdens,  all_Ptv_Burdens_per_line[, c("line", "Global_burden")], by = "line")
colnames(line_prop_changes_old_vs_young_premac_burdens)[dim(line_prop_changes_old_vs_young_premac_burdens)[2]] <- "Global_Burden_Ptv"
dim(line_prop_changes_old_vs_young_premac_burdens)
# [1] 294  24

## B. In SKAT-O genes associated with young -> old preMac based on their del variants
line_prop_changes_old_vs_young_premac_burdens <- inner_join(line_prop_changes_old_vs_young_premac_burdens,  all_Ptv_Burdens_per_line[, c("line", "SKAT0_del_old_vs_young_premac_burden")], by = "line")
colnames(line_prop_changes_old_vs_young_premac_burdens)[dim(line_prop_changes_old_vs_young_premac_burdens)[2]] <- "SKAT0_del_old_vs_young_premac_Burden_Ptv"

## C. In genes associated with young -> old preMac based on their del burdens
line_prop_changes_old_vs_young_premac_burdens <- inner_join(line_prop_changes_old_vs_young_premac_burdens,  all_Ptv_Burdens_per_line[, c("line", "Burden_del_old_vs_young_premac_burden")], by = "line")
colnames(line_prop_changes_old_vs_young_premac_burdens)[dim(line_prop_changes_old_vs_young_premac_burdens)[2]] <- "Burden_del_old_vs_young_premac_Burden_Ptv"
##     i. Genes with + beta
line_prop_changes_old_vs_young_premac_burdens <- inner_join(line_prop_changes_old_vs_young_premac_burdens,  all_Ptv_Burdens_per_line[, c("line", "Burden_del_old_vs_young_premac_burden_pos")], by = "line")
colnames(line_prop_changes_old_vs_young_premac_burdens)[dim(line_prop_changes_old_vs_young_premac_burdens)[2]] <- "Burden_del_old_vs_young_premac_pos_Burden_Ptv"

#~~~~~~~~~~~~~~~~~~~~~~~ Burdens of synonymous variants ~~~~~~~~~~~~~~~~~~~~~~~~

## A. Global (across all genes)
line_prop_changes_old_vs_young_premac_burdens <- inner_join(line_prop_changes_old_vs_young_premac_burdens,  all_Syn_Burdens_per_line[, c("line", "Global_burden")], by = "line")
colnames(line_prop_changes_old_vs_young_premac_burdens)[dim(line_prop_changes_old_vs_young_premac_burdens)[2]] <- "Global_Burden_Syn"
dim(line_prop_changes_old_vs_young_premac_burdens)
# [1] 294  28

## B. In SKAT-O genes associated with young -> old preMac based on their del variants
line_prop_changes_old_vs_young_premac_burdens <- inner_join(line_prop_changes_old_vs_young_premac_burdens,  all_Syn_Burdens_per_line[, c("line", "SKAT0_del_old_vs_young_premac_burden")], by = "line")
colnames(line_prop_changes_old_vs_young_premac_burdens)[dim(line_prop_changes_old_vs_young_premac_burdens)[2]] <- "SKAT0_del_old_vs_young_premac_Burden_Syn"

## C. In genes associated with young -> old preMac based on their del burdens
line_prop_changes_old_vs_young_premac_burdens <- inner_join(line_prop_changes_old_vs_young_premac_burdens,  all_Syn_Burdens_per_line[, c("line", "Burden_del_old_vs_young_premac_burden")], by = "line")
colnames(line_prop_changes_old_vs_young_premac_burdens)[dim(line_prop_changes_old_vs_young_premac_burdens)[2]] <- "Burden_del_old_vs_young_premac_Burden_Syn"
##     i. Genes with + beta
line_prop_changes_old_vs_young_premac_burdens <- inner_join(line_prop_changes_old_vs_young_premac_burdens,  all_Syn_Burdens_per_line[, c("line", "Burden_del_old_vs_young_premac_burden_pos")], by = "line")
colnames(line_prop_changes_old_vs_young_premac_burdens)[dim(line_prop_changes_old_vs_young_premac_burdens)[2]] <- "Burden_del_old_vs_young_premac_pos_Burden_Syn"


##################  preMac -> microglia proportion changes  ####################

#~~~~~~~~~~~~~~~~~~~~~~ Burdens of deleterious variants ~~~~~~~~~~~~~~~~~~~~~~~~

## A. Global (across all genes)
line_prop_changes_microglia_vs_premac_burdens <- inner_join(line_prop_changes_microglia_vs_premac,  all_Del_Burdens_per_line[, c("line", "Global_burden")], by = "line")
colnames(line_prop_changes_microglia_vs_premac_burdens)[dim(line_prop_changes_microglia_vs_premac_burdens)[2]] <- "Global_Burden_Del"
dim(line_prop_changes_microglia_vs_premac_burdens)
# [1] 1940   23

## B. In SKAT-O genes associated with young -> old preMac based on their del variants
line_prop_changes_microglia_vs_premac_burdens <- inner_join(line_prop_changes_microglia_vs_premac_burdens,  all_Del_Burdens_per_line[, c("line", "SKAT0_del_microglia_vs_premac_burden")], by = "line")
colnames(line_prop_changes_microglia_vs_premac_burdens)[dim(line_prop_changes_microglia_vs_premac_burdens)[2]] <- "SKAT0_del_microglia_vs_premac_Burden_Del"

## C. In genes associated with young -> old preMac based on their del burdens
line_prop_changes_microglia_vs_premac_burdens <- inner_join(line_prop_changes_microglia_vs_premac_burdens,  all_Del_Burdens_per_line[, c("line", "Burden_del_microglia_vs_premac_burden")], by = "line")
colnames(line_prop_changes_microglia_vs_premac_burdens)[dim(line_prop_changes_microglia_vs_premac_burdens)[2]] <- "Burden_del_microglia_vs_premac_Burden_Del"
##     i. Genes with + beta
line_prop_changes_microglia_vs_premac_burdens <- inner_join(line_prop_changes_microglia_vs_premac_burdens,  all_Del_Burdens_per_line[, c("line", "Burden_del_microglia_vs_premac_burden_pos")], by = "line")
colnames(line_prop_changes_microglia_vs_premac_burdens)[dim(line_prop_changes_microglia_vs_premac_burdens)[2]] <- "Burden_del_microglia_vs_premac_pos_Burden_Del"
##     ii. Genes with - beta
line_prop_changes_microglia_vs_premac_burdens <- inner_join(line_prop_changes_microglia_vs_premac_burdens,  all_Del_Burdens_per_line[, c("line", "Burden_del_microglia_vs_premac_burden_neg")], by = "line")
colnames(line_prop_changes_microglia_vs_premac_burdens)[dim(line_prop_changes_microglia_vs_premac_burdens)[2]] <- "Burden_del_microglia_vs_premac_neg_Burden_Del"

#~~~~~~~~~~~~~~~~~~~ Burdens of protein-truncating variants ~~~~~~~~~~~~~~~~~~~~

## A. Global (across all genes)
line_prop_changes_microglia_vs_premac_burdens <- inner_join(line_prop_changes_microglia_vs_premac_burdens,  all_Ptv_Burdens_per_line[, c("line", "Global_burden")], by = "line")
colnames(line_prop_changes_microglia_vs_premac_burdens)[dim(line_prop_changes_microglia_vs_premac_burdens)[2]] <- "Global_Burden_Ptv"
dim(line_prop_changes_microglia_vs_premac_burdens)
# [1] 1940   28

## B. In SKAT-O genes associated with young -> old preMac based on their del variants
line_prop_changes_microglia_vs_premac_burdens <- inner_join(line_prop_changes_microglia_vs_premac_burdens,  all_Ptv_Burdens_per_line[, c("line", "SKAT0_del_microglia_vs_premac_burden")], by = "line")
colnames(line_prop_changes_microglia_vs_premac_burdens)[dim(line_prop_changes_microglia_vs_premac_burdens)[2]] <- "SKAT0_del_microglia_vs_premac_Burden_Ptv"

## C. In genes associated with young -> old preMac based on their del burdens
line_prop_changes_microglia_vs_premac_burdens <- inner_join(line_prop_changes_microglia_vs_premac_burdens,  all_Ptv_Burdens_per_line[, c("line", "Burden_del_microglia_vs_premac_burden")], by = "line")
colnames(line_prop_changes_microglia_vs_premac_burdens)[dim(line_prop_changes_microglia_vs_premac_burdens)[2]] <- "Burden_del_microglia_vs_premac_Burden_Ptv"
##     i. Genes with + beta
line_prop_changes_microglia_vs_premac_burdens <- inner_join(line_prop_changes_microglia_vs_premac_burdens,  all_Ptv_Burdens_per_line[, c("line", "Burden_del_microglia_vs_premac_burden_pos")], by = "line")
colnames(line_prop_changes_microglia_vs_premac_burdens)[dim(line_prop_changes_microglia_vs_premac_burdens)[2]] <- "Burden_del_microglia_vs_premac_pos_Burden_Ptv"
##     ii. Genes with - beta
line_prop_changes_microglia_vs_premac_burdens <- inner_join(line_prop_changes_microglia_vs_premac_burdens,  all_Ptv_Burdens_per_line[, c("line", "Burden_del_microglia_vs_premac_burden_neg")], by = "line")
colnames(line_prop_changes_microglia_vs_premac_burdens)[dim(line_prop_changes_microglia_vs_premac_burdens)[2]] <- "Burden_del_microglia_vs_premac_neg_Burden_Ptv"

#~~~~~~~~~~~~~~~~~~~~~~~ Burdens of synonymous variants ~~~~~~~~~~~~~~~~~~~~~~~~

## A. Global (across all genes)
line_prop_changes_microglia_vs_premac_burdens <- inner_join(line_prop_changes_microglia_vs_premac_burdens,  all_Syn_Burdens_per_line[, c("line", "Global_burden")], by = "line")
colnames(line_prop_changes_microglia_vs_premac_burdens)[dim(line_prop_changes_microglia_vs_premac_burdens)[2]] <- "Global_Burden_Syn"
dim(line_prop_changes_microglia_vs_premac_burdens)
# [1] 1940   33

## B. In SKAT-O genes associated with young -> old preMac based on their del variants
line_prop_changes_microglia_vs_premac_burdens <- inner_join(line_prop_changes_microglia_vs_premac_burdens,  all_Syn_Burdens_per_line[, c("line", "SKAT0_del_microglia_vs_premac_burden")], by = "line")
colnames(line_prop_changes_microglia_vs_premac_burdens)[dim(line_prop_changes_microglia_vs_premac_burdens)[2]] <- "SKAT0_del_microglia_vs_premac_Burden_Syn"

## C. In genes associated with young -> old preMac based on their del burdens
line_prop_changes_microglia_vs_premac_burdens <- inner_join(line_prop_changes_microglia_vs_premac_burdens,  all_Syn_Burdens_per_line[, c("line", "Burden_del_microglia_vs_premac_burden")], by = "line")
colnames(line_prop_changes_microglia_vs_premac_burdens)[dim(line_prop_changes_microglia_vs_premac_burdens)[2]] <- "Burden_del_microglia_vs_premac_Burden_Syn"
##     i. Genes with + beta
line_prop_changes_microglia_vs_premac_burdens <- inner_join(line_prop_changes_microglia_vs_premac_burdens,  all_Syn_Burdens_per_line[, c("line", "Burden_del_microglia_vs_premac_burden_pos")], by = "line")
colnames(line_prop_changes_microglia_vs_premac_burdens)[dim(line_prop_changes_microglia_vs_premac_burdens)[2]] <- "Burden_del_microglia_vs_premac_pos_Burden_Syn"
##     ii. Genes with - beta
line_prop_changes_microglia_vs_premac_burdens <- inner_join(line_prop_changes_microglia_vs_premac_burdens,  all_Syn_Burdens_per_line[, c("line", "Burden_del_microglia_vs_premac_burden_neg")], by = "line")
colnames(line_prop_changes_microglia_vs_premac_burdens)[dim(line_prop_changes_microglia_vs_premac_burdens)[2]] <- "Burden_del_microglia_vs_premac_neg_Burden_Syn"



## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#                        4.1.2 Fit linear mixed models
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

burden_test_lmm <- function(stage_comparison, variant_type, genes_option, genes_beta_sign){
  
  if(stage_comparison == "premac_vs_iPSC"){
    y_axis_lab <- "iPSC -> young preMac"
    data <- line_prop_changes_premac_vs_iPSC_burdens
    fixed_effects <- c("sex", "prop_unadjusted_min_value", "PC1", "PC2")
    random_effects <- c("(1 | pool)")
  }
  
  else if(stage_comparison == "old_vs_young_premac"){
    y_axis_lab <- "young -> old preMac"
    data <- line_prop_changes_old_vs_young_premac_burdens
    fixed_effects <- c("sex", "prop_unadjusted_min_value", "PC1", "PC2")
    random_effects <- c("(1 | pool)")
  }
  
  else if(stage_comparison == "microglia_vs_premac"){
    y_axis_lab <- "preMac -> microglia"
    data <- line_prop_changes_microglia_vs_premac_burdens
    fixed_effects <- c("sex", "prop_unadjusted_min_value", "PC1", "PC2")
    random_effects <- c("(1 | pool)", "(1 | treatment)")
  }
  
  
  ## Burden variable 
  if(genes_option == "Global"){ coef <- paste0(genes_option, "_Burden_", variant_type) }
  else{
    if(!is.null(genes_beta_sign)){coef <- paste0(genes_option, "_del_", stage_comparison, "_", genes_beta_sign, "_Burden_", variant_type) }
    else{coef <- paste0(genes_option, "_del_", stage_comparison, "_Burden_", variant_type) }
    }
   
  data[, coef] <- as.numeric(data[, coef])
  
  ## Model
  formula <- as.formula(paste("scaled_log_fraction", " ~ ", paste(c(coef, fixed_effects, random_effects), collapse = " + ")))
  ## Fit LMM
  fit <- lmerTest::lmer(formula, data)
  sum_fit = summary(fit)
  
  
  ## Plot burden vs log-fraction
  variant_type_names <- c("Del" = "deleterious", 
                          "Ptv" = "protein-truncating", 
                          "Syn" = "synonymous")
  if(genes_option == "Global"){ x_axis_lab <- paste0("Global burden of ", variant_type_names[variant_type], " variants") }
  else if(genes_option == "SKAT0"){ x_axis_lab <- paste0("Burden of ", variant_type_names[variant_type], " variants in SKAT-O signif. genes") }
  else{
    if(is.null(genes_beta_sign)){x_axis_lab <- paste0("Burden of ", variant_type_names[variant_type], " variants in burden test signif. genes") }
    else{ x_axis_lab <- paste0("Burden of ", variant_type_names[variant_type], " variants in burden test signif. ", genes_beta_sign, " genes") }
  }
  
  plot <- jtools::effect_plot(fit, 
                      pred = !!rlang::sym(coef), 
                      interval = FALSE, partial.residuals = TRUE,
                      jitter = 0, plot.points = TRUE)  + 
    ggtitle(paste("Effect of mutational burden on", y_axis_lab), 
            subtitle = paste("Beta:", signif(sum_fit$coefficients[coef, "Estimate"], digits = 3), '    ', 
                                   'p:', signif(sum_fit$coefficients[coef, "Pr(>|t|)"], digits = 3) )) + 
    xlab(x_axis_lab) + ylab("Scaled log-fraction of line abundance") +
    theme(plot.title = element_text(size = 10, face = "bold"), 
          plot.subtitle = element_text(size = 9.5, color = "gray40"),
          axis.title = element_text(size = (10)))
  
  
  print(formula)
  print(sum_fit$coefficients[coef, ])
  return(list(fit, sum_fit, plot))
}


## Fit LMM to proportion changes from iPSC -> preMac
# - Deleterious variants
lmm_premac_vs_iPSC_burden_Del_Global <- burden_test_lmm("premac_vs_iPSC", "Del", "Global", NULL)
lmm_premac_vs_iPSC_burden_Del_SKAT0 <- burden_test_lmm("premac_vs_iPSC", "Del", "SKAT0", NULL)
lmm_premac_vs_iPSC_burden_Del_Burden <- burden_test_lmm("premac_vs_iPSC", "Del", "Burden", NULL)
lmm_premac_vs_iPSC_burden_Del_Burden_pos <- burden_test_lmm("premac_vs_iPSC", "Del", "Burden", "pos")
lmm_premac_vs_iPSC_burden_Del_Burden_neg <- burden_test_lmm("premac_vs_iPSC", "Del", "Burden", "neg")

# - Ptv variants
lmm_premac_vs_iPSC_burden_Ptv_Global <- burden_test_lmm("premac_vs_iPSC", "Ptv", "Global", NULL)
lmm_premac_vs_iPSC_burden_Ptv_SKAT0 <- burden_test_lmm("premac_vs_iPSC", "Ptv", "SKAT0", NULL)
lmm_premac_vs_iPSC_burden_Ptv_Burden <- burden_test_lmm("premac_vs_iPSC", "Ptv", "Burden", NULL)
lmm_premac_vs_iPSC_burden_Ptv_Burden_pos <- burden_test_lmm("premac_vs_iPSC", "Ptv", "Burden", "pos")
lmm_premac_vs_iPSC_burden_Ptv_Burden_neg <- burden_test_lmm("premac_vs_iPSC", "Ptv", "Burden", "neg")

# - Syn variants
lmm_premac_vs_iPSC_burden_Syn_Global <- burden_test_lmm("premac_vs_iPSC", "Syn", "Global", NULL)
lmm_premac_vs_iPSC_burden_Syn_SKAT0 <- burden_test_lmm("premac_vs_iPSC", "Syn", "SKAT0", NULL)
lmm_premac_vs_iPSC_burden_Syn_Burden <- burden_test_lmm("premac_vs_iPSC", "Syn", "Burden", NULL)
lmm_premac_vs_iPSC_burden_Syn_Burden_pos <- burden_test_lmm("premac_vs_iPSC", "Syn", "Burden", "pos")
lmm_premac_vs_iPSC_burden_Syn_Burden_neg <- burden_test_lmm("premac_vs_iPSC", "Syn", "Burden", "neg")

plot_grid(plotlist  = list(lmm_premac_vs_iPSC_burden_Del_Global[[3]], lmm_premac_vs_iPSC_burden_Del_SKAT0[[3]], 
                           lmm_premac_vs_iPSC_burden_Del_Burden[[3]], lmm_premac_vs_iPSC_burden_Del_Burden_pos[[3]],
                           lmm_premac_vs_iPSC_burden_Del_Burden_neg[[3]], 
                           lmm_premac_vs_iPSC_burden_Ptv_Global[[3]], lmm_premac_vs_iPSC_burden_Ptv_SKAT0[[3]], 
                           lmm_premac_vs_iPSC_burden_Ptv_Burden[[3]], lmm_premac_vs_iPSC_burden_Ptv_Burden_pos[[3]],
                           lmm_premac_vs_iPSC_burden_Ptv_Burden_neg[[3]],
                           lmm_premac_vs_iPSC_burden_Syn_Global[[3]], lmm_premac_vs_iPSC_burden_Syn_SKAT0[[3]], 
                           lmm_premac_vs_iPSC_burden_Syn_Burden[[3]], lmm_premac_vs_iPSC_burden_Syn_Burden_pos[[3]],
                           lmm_premac_vs_iPSC_burden_Syn_Burden_neg[[3]]), 
          nrow = 3, align = "vh")
ggsave(filename = paste0(plotdir, "/Burden_effects_on_premac_vs_iPSC.pdf"), height = 14, width = 27)



## Fit LMM to proportion changes from young -> old preMac
# - Deleterious variants
lmm_old_vs_young_premac_burden_Del_Global <- burden_test_lmm("old_vs_young_premac", "Del", "Global", NULL)
lmm_old_vs_young_premac_burden_Del_SKAT0 <- burden_test_lmm("old_vs_young_premac", "Del", "SKAT0", NULL)
lmm_old_vs_young_premac_burden_Del_Burden <- burden_test_lmm("old_vs_young_premac", "Del", "Burden", NULL)
lmm_old_vs_young_premac_burden_Del_Burden_pos <- burden_test_lmm("old_vs_young_premac", "Del", "Burden", "pos")

# - Ptv variants
lmm_old_vs_young_premac_burden_Ptv_Global <- burden_test_lmm("old_vs_young_premac", "Ptv", "Global", NULL)
lmm_old_vs_young_premac_burden_Ptv_SKAT0 <- burden_test_lmm("old_vs_young_premac", "Ptv", "SKAT0", NULL)
lmm_old_vs_young_premac_burden_Ptv_Burden <- burden_test_lmm("old_vs_young_premac", "Ptv", "Burden", NULL)
lmm_old_vs_young_premac_burden_Ptv_Burden_pos <- burden_test_lmm("old_vs_young_premac", "Ptv", "Burden", "pos")

# - Syn variants
lmm_old_vs_young_premac_burden_Syn_Global <- burden_test_lmm("old_vs_young_premac", "Syn", "Global", NULL)
lmm_old_vs_young_premac_burden_Syn_SKAT0 <- burden_test_lmm("old_vs_young_premac", "Syn", "SKAT0", NULL)
lmm_old_vs_young_premac_burden_Syn_Burden <- burden_test_lmm("old_vs_young_premac", "Syn", "Burden", NULL)
lmm_old_vs_young_premac_burden_Syn_Burden_pos <- burden_test_lmm("old_vs_young_premac", "Syn", "Burden", "pos")

plot_grid(plotlist  = list(lmm_old_vs_young_premac_burden_Del_Global[[3]], lmm_old_vs_young_premac_burden_Del_SKAT0[[3]], 
                           lmm_old_vs_young_premac_burden_Del_Burden[[3]], lmm_old_vs_young_premac_burden_Del_Burden_pos[[3]],
                           lmm_old_vs_young_premac_burden_Ptv_Global[[3]], lmm_old_vs_young_premac_burden_Ptv_SKAT0[[3]], 
                           lmm_old_vs_young_premac_burden_Ptv_Burden[[3]], lmm_old_vs_young_premac_burden_Ptv_Burden_pos[[3]],
                           lmm_old_vs_young_premac_burden_Syn_Global[[3]], lmm_old_vs_young_premac_burden_Syn_SKAT0[[3]], 
                           lmm_old_vs_young_premac_burden_Syn_Burden[[3]], lmm_old_vs_young_premac_burden_Syn_Burden_pos[[3]]), 
          nrow = 3, align = "vh")
ggsave(filename = paste0(plotdir, "/Burden_effects_on_old_vs_young_premac.pdf"), height = 14, width = 20)


## Fit LMM to proportion changes from preMac -> microglia
# - Deleterious variants
lmm_microglia_vs_premac_burden_Del_Global <- burden_test_lmm("microglia_vs_premac", "Del", "Global", NULL)
lmm_microglia_vs_premac_burden_Del_SKAT0 <- burden_test_lmm("microglia_vs_premac", "Del", "SKAT0", NULL)
lmm_microglia_vs_premac_burden_Del_Burden <- burden_test_lmm("microglia_vs_premac", "Del", "Burden", NULL)
lmm_microglia_vs_premac_burden_Del_Burden_pos <- burden_test_lmm("microglia_vs_premac", "Del", "Burden", "pos")
lmm_microglia_vs_premac_burden_Del_Burden_neg <- burden_test_lmm("microglia_vs_premac", "Del", "Burden", "neg")

# - Ptv variants
lmm_microglia_vs_premac_burden_Ptv_Global <- burden_test_lmm("microglia_vs_premac", "Ptv", "Global", NULL)
lmm_microglia_vs_premac_burden_Ptv_SKAT0  <- burden_test_lmm("microglia_vs_premac", "Ptv", "SKAT0", NULL)
lmm_microglia_vs_premac_burden_Ptv_Burden <- burden_test_lmm("microglia_vs_premac", "Ptv", "Burden", NULL)
lmm_microglia_vs_premac_burden_Ptv_Burden_pos <- burden_test_lmm("microglia_vs_premac", "Ptv", "Burden", "pos")
lmm_microglia_vs_premac_burden_Ptv_Burden_neg <- burden_test_lmm("microglia_vs_premac", "Ptv", "Burden", "neg")

# - Syn variants
lmm_microglia_vs_premac_burden_Syn_Global <- burden_test_lmm("microglia_vs_premac", "Syn", "Global", NULL)
lmm_microglia_vs_premac_burden_Syn_SKAT0 <- burden_test_lmm("microglia_vs_premac", "Syn", "SKAT0", NULL)
lmm_microglia_vs_premac_burden_Syn_Burden <- burden_test_lmm("microglia_vs_premac", "Syn", "Burden", NULL)
lmm_microglia_vs_premac_burden_Syn_Burden_pos <- burden_test_lmm("microglia_vs_premac", "Syn", "Burden", "pos")
lmm_microglia_vs_premac_burden_Syn_Burden_neg <- burden_test_lmm("microglia_vs_premac", "Syn", "Burden", "neg")

plot_grid(plotlist  = list(lmm_microglia_vs_premac_burden_Del_Global[[3]], lmm_microglia_vs_premac_burden_Del_SKAT0[[3]], 
                           lmm_microglia_vs_premac_burden_Del_Burden[[3]], lmm_microglia_vs_premac_burden_Del_Burden_pos[[3]],
                           lmm_microglia_vs_premac_burden_Del_Burden_neg[[3]],
                           lmm_microglia_vs_premac_burden_Ptv_Global[[3]], lmm_microglia_vs_premac_burden_Ptv_SKAT0[[3]], 
                           lmm_microglia_vs_premac_burden_Ptv_Burden[[3]], lmm_microglia_vs_premac_burden_Ptv_Burden_pos[[3]],
                           lmm_microglia_vs_premac_burden_Ptv_Burden_neg[[3]],
                           lmm_microglia_vs_premac_burden_Syn_Global[[3]], lmm_microglia_vs_premac_burden_Syn_SKAT0[[3]], 
                           lmm_microglia_vs_premac_burden_Syn_Burden[[3]], lmm_microglia_vs_premac_burden_Syn_Burden_pos[[3]],
                           lmm_microglia_vs_premac_burden_Syn_Burden_neg[[3]]), 
          nrow = 3, align = "vh")
ggsave(filename = paste0(plotdir, "/Burden_effects_on_microglia_vs_premac.pdf"), height = 14, width = 27)



## Improved plots for del burdens
nicer_plots_for_del_burdens <- function(stage_comparison){
  
  if(stage_comparison == "premac_vs_iPSC"){
    comparison_label <- "iPSC to young macrophage precursor"
    data <- line_prop_changes_premac_vs_iPSC_burdens
  }
  
  else if(stage_comparison == "old_vs_young_premac"){
    comparison_label <- "young to old macrophage precursor"
    data <- line_prop_changes_old_vs_young_premac_burdens
  }
  
  else if(stage_comparison == "microglia_vs_premac"){
    comparison_label <- "macrophage precursor to microglia"
    data <- line_prop_changes_microglia_vs_premac_burdens
  }
  
  ## Results for del burdens
  lmm_results_burden_Del_Global <- eval(parse_expr(paste0("lmm_", stage_comparison, "_burden_Del_Global")))
  lmm_results_burden_Del_SKAT0 <- eval(parse_expr(paste0("lmm_", stage_comparison, "_burden_Del_SKAT0")))
  lmm_results_burden_Del_Burden_pos <- eval(parse_expr(paste0("lmm_", stage_comparison, "_burden_Del_Burden_pos")))
  if(stage_comparison != "old_vs_young_premac"){
    lmm_results_burden_Del_Burden_neg <- eval(parse_expr(paste0("lmm_", stage_comparison, "_burden_Del_Burden_neg")))
  }
  
  ## Scatterplots for global and SKAT-O burdens
  coef <- "Global_Burden_Del"
  data_new <- jtools::make_predictions(lmm_results_burden_Del_Global[[1]], 
                                       pred = coef, 
                                       partial.residuals = TRUE)$data %>%  as.data.frame()
  
  p_global <- ggplot(data_new, aes(x = !!rlang::sym(coef), y = scaled_log_fraction, colour = !!rlang::sym(coef))) + 
                  geom_jitter(width = 0.1, alpha = 0.7, size = 0.9) +
                  scale_color_gradient2(low = "gray87", 
                                        midpoint = min(data_new[,coef]) + ((max(data_new[,coef]) - min(data_new[,coef])) /2), 
                                        mid = "gray62", high = "gray1") +
                  guides(colour = "none") +
                  stat_smooth(geom="line", alpha = 1, size = 1, method = lm, aes(group=1), color='black') +
                  labs(title = "", 
                       subtitle = paste("Beta:", signif(lmm_results_burden_Del_Global[[2]]$coefficients[coef, "Estimate"], digits = 3), '    ', 
                                        'p:', signif(lmm_results_burden_Del_Global[[2]]$coefficients[coef, "Pr(>|t|)"], digits = 3) )) + 
                  xlab("Global burden") + ylab(paste0("Proliferation from ", comparison_label)) +
                  theme_classic() +
                  theme(plot.title = element_text(size = 10, face = "bold"), 
                        plot.subtitle = element_text(size = 9.5, color = "gray40"),
                        axis.title = element_text(size = (10)))
    
    
  coef <- paste0("SKAT0_del_", stage_comparison, "_Burden_Del")
  data_new <- jtools::make_predictions(lmm_results_burden_Del_SKAT0[[1]], 
                                       pred = coef, 
                                       partial.residuals = TRUE)$data %>%  as.data.frame()
  
  p_SKATO <- ggplot(data_new, aes(x = !!rlang::sym(coef), y = scaled_log_fraction, colour = !!rlang::sym(coef))) + 
                  geom_jitter(width = 0.1, alpha = 0.7, size = 0.9) +
                  scale_color_gradient2(low = "gray87", 
                                        midpoint = min(data_new[,coef]) + ((max(data_new[,coef]) - min(data_new[,coef])) /2), 
                                        mid = "gray62", high = "gray1") +
                  guides(colour = "none") +
                  stat_smooth(geom="line", alpha = 1, size = 1, method = lm, aes(group=1), color='black') +
                  labs(title = "", 
                       subtitle = paste("Beta:", signif(lmm_results_burden_Del_SKAT0[[2]]$coefficients[coef, "Estimate"], digits = 3), '    ', 
                                        'p:', signif(lmm_results_burden_Del_SKAT0[[2]]$coefficients[coef, "Pr(>|t|)"], digits = 3) )) + 
                  xlab("Burden in SKAT-O genes") + ylab(paste0("Proliferation from ", comparison_label)) +
                  theme_classic() +
                  theme(plot.title = element_text(size = 10, face = "bold"), 
                        plot.subtitle = element_text(size = 9.5, color = "gray40"),
                        axis.title = element_text(size = (10)))
                
  
  ## Violin plots for del burdens in burden test genes
  set.seed(12062024)
  ## Color gradient for burdens
  colfunc <- colorRampPalette(c("gray87", "gray62", "gray1"))
  burden_colors <- colfunc(20)[2:20]
  names(burden_colors) <- 0:18
  
  ## - Positive genes
  coef <- paste0("Burden_del_", stage_comparison, "_pos_Burden_Del")
  ## Use scaled log-fractions with covariates regressed out
  data_new <- jtools::make_predictions(lmm_results_burden_Del_Burden_pos[[1]], 
                                        pred = coef, 
                                        partial.residuals = TRUE)$data %>%  as.data.frame()
  
  data_new[, coef] <- factor(data[, coef], levels = sort(unique(data_new[, coef])))
  
  v_burden_pos <- ggplot(data_new, aes(x = !!rlang::sym(coef), y = scaled_log_fraction, color = !!rlang::sym(coef))) + 
    geom_jitter(width = 0.1, alpha = 0.9, size = 0.9) +
    geom_violin(fill = NA, color = "black", linewidth = 0.4) +
    scale_color_manual(values = burden_colors) +
    stat_smooth(geom="line", alpha = 1, size = 1, method = lm, aes(group=1), color='black') +
    theme_classic() +
    guides(color = "none") +
    labs(title = "", 
        subtitle = paste("Beta:", signif(lmm_results_burden_Del_Burden_pos[[2]]$coefficients[coef, "Estimate"], digits = 3), '    ', 
                      'p:', signif(lmm_results_burden_Del_Burden_pos[[2]]$coefficients[coef, "Pr(>|t|)"], digits = 3) ), 
        x = "Burden in genes whose del burden increase proliferation", y = paste0("Proliferation from ", comparison_label)) +
    theme(plot.title = element_text(size = 10, face = "bold"), 
          plot.subtitle = element_text(size = 9.5, color = "gray40"),
          axis.title = element_text(size = (10)))
  
  
  ## - Negative genes
  if (stage_comparison != "old_vs_young_premac"){
    
    coef <- paste0("Burden_del_", stage_comparison, "_neg_Burden_Del")
    data_new <- jtools::make_predictions(lmm_results_burden_Del_Burden_neg[[1]], 
                                         pred = coef, 
                                         partial.residuals = TRUE)$data %>%  as.data.frame()
    
    data_new[, coef] <- as.factor(data[, coef])
    
    v_burden_neg <- ggplot(data_new, aes(x = !!rlang::sym(coef), y = scaled_log_fraction, color = !!rlang::sym(coef))) + 
      geom_jitter(width = 0.1, alpha = 0.7, size = 0.9) +
      geom_violin(fill = NA, color = "black", linewidth = 0.4) +
      scale_color_manual(values = burden_colors) +
      stat_smooth(geom="line", alpha = 1, size = 1, method = lm, aes(group=1), color='black') +
      theme_classic() +
      guides(color = "none") +
      labs(title = "", 
           subtitle = paste("Beta:", signif(lmm_results_burden_Del_Burden_neg[[2]]$coefficients[coef, "Estimate"], digits = 3), '    ', 
                            'p:', signif(lmm_results_burden_Del_Burden_neg[[2]]$coefficients[coef, "Pr(>|t|)"], digits = 3) ), 
           x = "Burden in genes whose del burden decrease proliferation", y = paste0("Proliferation from ", comparison_label)) +
      theme(plot.title = element_text(size = 10, face = "bold"), 
            plot.subtitle = element_text(size = 9.5, color = "gray40"),
            axis.title = element_text(size = (10)))
  }
  
  else { v_burden_neg <- NULL }
  
  return(list(p_global, p_SKATO, v_burden_pos, v_burden_neg))
  
}


l1 <- nicer_plots_for_del_burdens("premac_vs_iPSC")
l2 <- nicer_plots_for_del_burdens("old_vs_young_premac")
l3 <- nicer_plots_for_del_burdens("microglia_vs_premac")


plot_grid(plotlist = c(l1, l2, l3), ncol = 4, align = "h")
ggsave(filename = paste0(plotdir, "/Del_Burden_effects_3_comparisons.pdf"), height = 13, width = 16)







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
# date     2024-12-04
# rstudio  2024.04.0+735 Chocolate Cosmos (server)
# pandoc   3.1.12.3 @ /opt/view/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package     * version    date (UTC) lib source
# bit           4.0.5      2022-11-15 [1] CRAN (R 4.3.1)
# bit64         4.0.5      2020-08-30 [1] CRAN (R 4.3.1)
# boot          1.3-30     2024-02-26 [1] CRAN (R 4.3.1)
# cli           3.6.2      2023-12-11 [1] CRAN (R 4.3.1)
# colorspace    2.1-0      2023-01-23 [1] CRAN (R 4.3.1)
# cowplot     * 1.1.3      2024-01-22 [1] CRAN (R 4.3.1)
# crayon        1.5.2      2022-09-29 [1] CRAN (R 4.3.1)
# digest        0.6.35     2024-03-11 [1] CRAN (R 4.3.1)
# dplyr       * 1.1.4      2023-11-17 [1] CRAN (R 4.3.1)
# fansi         1.0.6      2023-12-08 [1] CRAN (R 4.3.1)
# farver        2.1.1      2022-07-06 [1] CRAN (R 4.3.1)
# forcats     * 1.0.0      2023-01-29 [1] CRAN (R 4.3.1)
# generics      0.1.3      2022-07-05 [1] CRAN (R 4.3.1)
# ggplot2     * 3.5.1      2024-04-23 [1] CRAN (R 4.3.1)
# glue          1.7.0      2024-01-09 [1] CRAN (R 4.3.1)
# gtable        0.3.4      2023-08-21 [1] CRAN (R 4.3.1)
# hms           1.1.3      2023-03-21 [1] CRAN (R 4.3.1)
# jtools        2.2.2      2023-07-11 [1] CRAN (R 4.3.1)
# labeling      0.4.3      2023-08-29 [1] CRAN (R 4.3.1)
# lattice       0.22-6     2024-03-20 [1] CRAN (R 4.3.1)
# lifecycle     1.0.4      2023-11-07 [1] CRAN (R 4.3.1)
# lme4        * 1.1-35.2   2024-03-28 [1] CRAN (R 4.3.1)
# lmerTest    * 3.1-3      2020-10-23 [1] CRAN (R 4.3.1)
# lubridate   * 1.9.3      2023-09-27 [1] CRAN (R 4.3.1)
# magrittr      2.0.3      2022-03-30 [1] CRAN (R 4.3.1)
# MASS          7.3-60.0.1 2024-01-13 [1] CRAN (R 4.3.1)
# Matrix      * 1.6-5      2024-01-11 [1] CRAN (R 4.3.1)
# minqa         1.2.6      2023-09-11 [1] CRAN (R 4.3.1)
# munsell       0.5.1      2024-04-01 [1] CRAN (R 4.3.1)
# nlme          3.1-164    2023-11-27 [1] CRAN (R 4.3.1)
# nloptr        2.0.3      2022-05-26 [1] CRAN (R 4.3.1)
# numDeriv      2016.8-1.1 2019-06-06 [1] CRAN (R 4.3.1)
# pander        0.6.5      2022-03-18 [1] CRAN (R 4.3.1)
# pillar        1.9.0      2023-03-22 [1] CRAN (R 4.3.1)
# pkgconfig     2.0.3      2019-09-22 [1] CRAN (R 4.3.1)
# purrr       * 1.0.2      2023-08-10 [1] CRAN (R 4.3.1)
# R6            2.5.1      2021-08-19 [1] CRAN (R 4.3.1)
# ragg          1.3.0      2024-03-13 [1] CRAN (R 4.3.1)
# Rcpp          1.0.12     2024-01-09 [1] CRAN (R 4.3.1)
# readr       * 2.1.5      2024-01-10 [1] CRAN (R 4.3.1)
# rlang       * 1.1.3      2024-01-10 [1] CRAN (R 4.3.1)
# rstudioapi    0.16.0     2024-03-24 [1] CRAN (R 4.3.1)
# scales        1.3.0      2023-11-28 [1] CRAN (R 4.3.1)
# sessioninfo * 1.2.2      2021-12-06 [1] CRAN (R 4.3.1)
# stringi       1.8.3      2023-12-11 [1] CRAN (R 4.3.1)
# stringr     * 1.5.1      2023-11-14 [1] CRAN (R 4.3.1)
# systemfonts   1.0.6      2024-03-07 [1] CRAN (R 4.3.1)
# textshaping   0.3.7      2023-10-09 [1] CRAN (R 4.3.1)
# tibble      * 3.2.1      2023-03-20 [1] CRAN (R 4.3.1)
# tidyr       * 1.3.1      2024-01-24 [1] CRAN (R 4.3.1)
# tidyselect    1.2.1      2024-03-11 [1] CRAN (R 4.3.1)
# tidyverse   * 2.0.0      2023-02-22 [1] CRAN (R 4.3.1)
# timechange    0.3.0      2024-01-18 [1] CRAN (R 4.3.1)
# tzdb          0.4.0      2023-05-12 [1] CRAN (R 4.3.1)
# utf8          1.2.4      2023-10-22 [1] CRAN (R 4.3.1)
# vctrs         0.6.5      2023-12-01 [1] CRAN (R 4.3.1)
# vroom         1.6.5      2023-12-05 [1] CRAN (R 4.3.1)
# withr         3.0.0      2024-01-16 [1] CRAN (R 4.3.1)
# 
# [1] /opt/view/rlib/R/library
# [2] /software/hgi/installs/softpack/rstudio/.spack/linux-ubuntu22.04-x86_64_v3/gcc-11.4.0/r-4.3.1-bfwldrk76z6f52upk47zepliekn7ayqz/rlib/R/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# 
