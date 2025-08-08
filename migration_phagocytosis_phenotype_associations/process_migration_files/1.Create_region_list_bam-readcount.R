# Extract from vcf the chr, start column, and add an end column
# Filter VCF to only the variants that are different across all donors
# Input: VCF of non-identical SNPs for all donors
# Output: for input SNPs, just chr, pos and pos (start and end are the same, it's 1bp)
.libPaths(c("/software/teamtrynka/conda/otar2065/lib/R/library","/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/"))

library(vcfR)
library(data.table)
library(tidyverse)
options(stringsAsFactors = FALSE)
gc()


# For server
args = commandArgs(trailingOnly=TRUE)


if (length(args)<2) {
  stop("You need to detail input and output file path.n", call.=FALSE)
} else if (length(args)==2) {
  vcf_path = args[1]
  output_path = args[2]

}

# Function to convert gt vcf to minor allele dosages

message("...Reading VCF...")


vcf = vcfR::read.vcfR(vcf_path)

message("...Getting right columns...")

subset = cbind(vcf@fix[,"CHROM"],vcf@fix[,"POS"])
subset=as.data.frame(subset)
print(head(subset))
subset$V2 = as.numeric(subset$V2)
subset$V3 = subset$V2
subset = as.data.table(subset)
fwrite(subset,file = output_path,
       sep = "\t", col.names = FALSE)
