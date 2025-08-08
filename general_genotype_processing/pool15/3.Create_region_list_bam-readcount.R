# Get the region list in right format for bam-readcount
# Extract from vcf the chr, start column, and add an end column
# Input: VCF of non-identical SNPs for all donors
# Output: for input SNPs, just chr, pos and pos (start and end are the same, it's 1bp)
.libPaths("/lustre/scratch123/hgi/projects/otar2065/bin/R/R-4.0.0/lib")

library(vcfR)
library(data.table)
library(tidyverse)
options(stringsAsFactors = FALSE)


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
