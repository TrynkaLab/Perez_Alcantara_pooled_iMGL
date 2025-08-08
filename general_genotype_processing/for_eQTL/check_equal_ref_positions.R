# Check positions that do/don't have the same reference nucleotide
.libPaths(c("/software/teamtrynka/conda/otar2065/lib/R/library","/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/"))
library(vcfR)
library(data.table)
.libPaths(c("/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/"))
library(stringr)
options(stringsAsFactors = FALSE)

# read in
larger = vcfR::read.vcfR(
  "../../data/all_pools_no_curn_iukl/all_pools_no_curn_iukl.genotype.MAF01.hg38.vcf.gz",
  nrows = 200000
)
# extract imputation info to check distribution of scores
larger_impute_info=larger@fix[,"INFO"]
larger_impute_info = stringr::str_split_i(larger_impute_info,pattern=";",i=3)
larger_impute_info = as.numeric(stringr::str_split_i(larger_impute_info,pattern="=",i=2))
hist(larger_impute_info)
# reduce to position and variant info
larger = as.data.table(larger@fix[, 1:5])
gc()
# same with others
curn =  vcfR::read.vcfR("../../data/all_pools_no_curn_iukl/curn_3.genotype_wgs.vcf.gz",
                        nrows = 200000)
curn = as.data.table(curn@fix[, 1:5])
gc()
iukl =  vcfR::read.vcfR("../../data/all_pools_no_curn_iukl/iukl_1.genotype_wgs.vcf.gz",
                        nrows = 200000)
iukl = as.data.table(iukl@fix[, 1:5])
gc()

# subset to same positions
larger$newPos = paste(larger$CHROM,larger$POS,sep = "_")
curn$newPos = paste(curn$CHROM,curn$POS,sep = "_")
iukl$newPos = paste(iukl$CHROM,iukl$POS,sep = "_")
names(larger_impute_info) = larger$newPos
merged_df = merge(larger, curn, by = "newPos", all = FALSE,suffixes = c("larger","curn"))
merged_df = merge(merged_df, iukl, by = "newPos", all = FALSE)

# I'm losing 3/4th of positions here -
# maybe I should use the full genomic for the smaller sets and just call the REF?

# how many positions have same / different REF?
merged_df$newIdLarger = paste(merged_df$newPos,merged_df$REFlarger,sep = "_")
merged_df$newIdCurn = paste(merged_df$newPos,merged_df$REFcurn,sep = "_")
merged_df$newIdIukl = paste(merged_df$newPos,merged_df$REF,sep = "_")
head(merged_df)
gc()

# compare refs
table(merged_df$newIdLarger %in% merged_df$newIdCurn )[1]/(table(merged_df$newIdLarger %in% merged_df$newIdCurn )[1] + table(merged_df$newIdLarger %in% merged_df$newIdCurn )[2])
table(merged_df$newIdLarger %in% merged_df$newIdIukl )[1]/(table(merged_df$newIdLarger %in% merged_df$newIdIukl )[1] + table(merged_df$newIdLarger %in% merged_df$newIdIukl )[2])
table(merged_df$newIdCurn %in% merged_df$newIdIukl )[1]/(table(merged_df$newIdCurn %in% merged_df$newIdIukl )[1] + table(merged_df$newIdCurn %in% merged_df$newIdIukl )[2])
# 0.3 - 0.4 % differences in REF

different_ref = !(merged_df$newIdLarger %in% merged_df$newIdCurn)
names(different_ref) = merged_df$newPos
# different ref - imputation values histogram
hist(larger_impute_info[names(different_ref[different_ref==TRUE])])
# same ref
hist(larger_impute_info[names(different_ref[different_ref==FALSE])])
# positions with different ref have inflated lower imputation values

# compare alts
merged_df$newAltLarger = paste(merged_df$newPos,merged_df$ALTlarger,sep = "_")
merged_df$newAltCurn = paste(merged_df$newPos,merged_df$ALTcurn,sep = "_")
merged_df$newAltIukl = paste(merged_df$newPos,merged_df$ALT,sep = "_")
table(merged_df$newAltLarger %in% merged_df$newAltCurn )[1]/(table(merged_df$newAltLarger %in% merged_df$newAltCurn )[1] + table(merged_df$newAltLarger %in% merged_df$newAltCurn )[2])
table(merged_df$newAltLarger %in% merged_df$newAltIukl )[1]/(table(merged_df$newAltLarger %in% merged_df$newAltIukl )[1] + table(merged_df$newAltLarger %in% merged_df$newAltIukl )[2])
table(merged_df$newAltCurn %in% merged_df$newAltIukl )[1]/(table(merged_df$newAltCurn %in% merged_df$newAltIukl )[1] + table(merged_df$newAltCurn %in% merged_df$newAltIukl )[2])
# 0.7 - 0.8 % differences in ALT
merged_df[323,]
merged_df[842,]
# I wonder how many donors have non-ref alleles, but the ALTs differ (multiallelic positions)?
# "True multiallelic sites are not observed very frequently unless you look at very large cohorts,
# so they are often taken as a sign of a noisy region where artifacts are likely."
# from GATK: https://gatk.broadinstitute.org/hc/en-us/articles/360035890771-Biallelic-vs-Multiallelic-sites

# check imputation scores for positions with different ALT
different_alt = !(merged_df$newAltLarger %in% merged_df$newAltCurn)
names(different_alt) = merged_df$newPos
# different ref - imputation values histogram
hist(larger_impute_info[names(different_alt[different_alt==TRUE])])
# same ref
hist(larger_impute_info[names(different_alt[different_alt==FALSE])])
# here it's harder to see differences.
# overall loss
table(larger_impute_info[merged_df$newPos]>=0.8)
7699 / (43322 + 7699) # 15% loss
43322 / 200000 # 79% overall loss
# check gVCF
curn =  vcfR::read.vcfR("../../data/ungenotyped_donors/nextflow_results/variant_calling/deepvariant/S1/S1.deepvariant.g.vcf.gz",
                        nrows = 200000)
curn = as.data.table(curn@fix[, 1:5])
curn$CHROM=paste0("chr",curn$CHROM)
larger = vcfR::read.vcfR(
  "../../data/all_pools_no_curn_iukl/all_pools_no_curn_iukl.genotype.MAF01.hg38.vcf.gz",
  nrows = 200000
)
larger = as.data.table(larger@fix[, 1:5])


larger$newPos = paste(larger$CHROM,larger$POS,sep = "_")
curn$newPos = paste(curn$CHROM,curn$POS,sep = "_")

merged_df = merge(larger, curn, by = "newPos", all = FALSE,suffixes = c("larger","curn"))
merged_df$newIdLarger = paste(merged_df$newPos,merged_df$REFlarger,sep = "_")
merged_df$newIdCurn = paste(merged_df$newPos,merged_df$REFcurn,sep = "_")

table(merged_df$newIdLarger %in% merged_df$newIdCurn )[1]/(table(merged_df$newIdLarger %in% merged_df$newIdCurn )[1] + table(merged_df$newIdLarger %in% merged_df$newIdCurn )[2])

head(curn@fix[curn@fix[,2]=="1000731"])
curn@gt[which(curn@fix[,2]=="1000731"),]
curn@gt[which(curn@fix[,2]=="1000112"),]

test = curn@gt[which(curn@fix[,2] %in% merged_df$POSlarger),]
# They are all coded as ALT but didn't make it to the final deepvariant VCF?
# have low quality
