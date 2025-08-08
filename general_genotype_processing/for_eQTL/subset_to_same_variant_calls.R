# subset files to be added to positions of equal REF and ALT
# after filtering to high imputation quality
.libPaths(c("/software/teamtrynka/conda/otar2065/lib/R/library","/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/"))
library(vcfR)
library(data.table)
.libPaths(c("/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/"))
library(stringr)
options(future.globals.maxSize= 30097152000) # 30Gb maxsize, 30000*1024^2


args = commandArgs(trailingOnly=TRUE)
larger_path = args[1]
curn_path = args[2]
iukl_path = args[3]
output_path = args[4]

# read in
larger = vcfR::read.vcfR(larger_path)

# extract imputation info to check distribution of scores
larger_impute_info=larger@fix[,"INFO"]
larger_impute_info = stringr::str_split_i(larger_impute_info,pattern=";",i=3)
larger_impute_info = as.numeric(stringr::str_split_i(larger_impute_info,pattern="=",i=2))
if (min(larger_impute_info) < 0.8) {
  stop("Minimum value of imputation 'INFO' is less than 0.8. Script stopped.")
}

# reduce to position and variant info
larger = as.data.table(larger@fix[, 1:5])
gc()
curn = vcfR::read.vcfR(curn_path)

curn = as.data.table(curn@fix[, 1:5])
gc()
iukl =  vcfR::read.vcfR(iukl_path)
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
# 0.2 - 0.3 % differences in REF

# compare alts
merged_df$newAltLarger = paste(merged_df$newPos,merged_df$ALTlarger,sep = "_")
merged_df$newAltCurn = paste(merged_df$newPos,merged_df$ALTcurn,sep = "_")
merged_df$newAltIukl = paste(merged_df$newPos,merged_df$ALT,sep = "_")
table(merged_df$newAltLarger %in% merged_df$newAltCurn )[1]/(table(merged_df$newAltLarger %in% merged_df$newAltCurn )[1] + table(merged_df$newAltLarger %in% merged_df$newAltCurn )[2])
table(merged_df$newAltLarger %in% merged_df$newAltIukl )[1]/(table(merged_df$newAltLarger %in% merged_df$newAltIukl )[1] + table(merged_df$newAltLarger %in% merged_df$newAltIukl )[2])
table(merged_df$newAltCurn %in% merged_df$newAltIukl )[1]/(table(merged_df$newAltCurn %in% merged_df$newAltIukl )[1] + table(merged_df$newAltCurn %in% merged_df$newAltIukl )[2])
# 0.7 - 0.8 % differences in ALT

# "True multiallelic sites are not observed very frequently unless you look at very large cohorts, 
# so they are often taken as a sign of a noisy region where artifacts are likely."
# from GATK: https://gatk.broadinstitute.org/hc/en-us/articles/360035890771-Biallelic-vs-Multiallelic-sites

# subsetting to identical REF values
# Filter rows where all three columns have identical values
merged_df = merged_df[newIdLarger == newIdCurn & newIdCurn == newIdIukl]
# same for ALTs
merged_df = merged_df[newAltLarger == newAltCurn & newAltCurn == newAltIukl]

# check again
table(merged_df$newIdLarger %in% merged_df$newIdCurn )[1]/(table(merged_df$newIdLarger %in% merged_df$newIdCurn )[1] + table(merged_df$newIdLarger %in% merged_df$newIdCurn )[2])
table(merged_df$newIdLarger %in% merged_df$newIdIukl )[1]/(table(merged_df$newIdLarger %in% merged_df$newIdIukl )[1] + table(merged_df$newIdLarger %in% merged_df$newIdIukl )[2])
table(merged_df$newIdCurn %in% merged_df$newIdIukl )[1]/(table(merged_df$newIdCurn %in% merged_df$newIdIukl )[1] + table(merged_df$newIdCurn %in% merged_df$newIdIukl )[2])

table(merged_df$newAltLarger %in% merged_df$newAltCurn )[1]/(table(merged_df$newAltLarger %in% merged_df$newAltCurn )[1] + table(merged_df$newAltLarger %in% merged_df$newAltCurn )[2])
table(merged_df$newAltLarger %in% merged_df$newAltIukl )[1]/(table(merged_df$newAltLarger %in% merged_df$newAltIukl )[1] + table(merged_df$newAltLarger %in% merged_df$newAltIukl )[2])
table(merged_df$newAltCurn %in% merged_df$newAltIukl )[1]/(table(merged_df$newAltCurn %in% merged_df$newAltIukl )[1] + table(merged_df$newAltCurn %in% merged_df$newAltIukl )[2])

#save coordinates for filtering VCF with bcftools
line="##fileformat=VCFv4.3"
write(line,file=output_path,append=FALSE)

merged_df = as.data.frame(merged_df[,c("CHROM","POS","REF","ALT")])
merged_df = merged_df %>%
  dplyr::mutate(ID = ".",
                QUAL=".",
                FILTER=".",
                INFO=".") %>%
dplyr::rename('#CHROM'=CHROM) 
merged_df = merged_df %>%
  dplyr::relocate('#CHROM',POS,ID,REF,ALT,QUAL,FILTER,INFO)
write.table(merged_df[,c('#CHROM',"POS","ID","REF","ALT","QUAL","FILTER","INFO")],
            file =  output_path,
            col.names = TRUE,row.names = FALSE,sep = "\t", quote = FALSE,
            append = TRUE)

# beware: iukl has different calls at same positions:
# chr10   1403654 .       G       A       46      PASS    .       GT:GQ:DP:AD:VAF:PL      1/1:46:27:0,26:0.962963:46,57,0
# chr10   1403654 .       T       C       45.3    PASS    .       GT:GQ:DP:AD:VAF:PL      0/1:45:27:12,15:0.555556:45,0,64

# chr10   1416766 .       G       A       56.4    PASS    .       GT:GQ:DP:AD:VAF:PL      1/1:47:17:0,17:1:56,47,0
# chr10   1416766 .       A       C       58.9    PASS    .       GT:GQ:DP:AD:VAF:PL      1/1:49:19:0,19:1:58,49,0
