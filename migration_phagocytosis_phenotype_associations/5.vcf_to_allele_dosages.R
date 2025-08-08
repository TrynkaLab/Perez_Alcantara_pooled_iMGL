# convert VCF to genotype allele dosages
message("Entering Rscript")
.libPaths(c("/software/teamtrynka/conda/otar2065/lib/R/library",
            "/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/"))
library(vcfR)
library(data.table)
library(tidyverse)
source("./functions.R")
options(stringsAsFactors = FALSE)

#For server
args = commandArgs(trailingOnly=TRUE)


if (length(args)<2) {
  stop("You need to detail 1 input and 1 output file paths.n",
       call. = FALSE)
} else if (length(args) == 2) {
  input_genotype_path = args[1]
  genotype_minor_allele_dos_path = args[2]
}


vcf = vcfR::read.vcfR(input_genotype_path)
vcf_reduced = as.data.table(vcf@fix[,1:5])
original_vcf_rows = nrow(vcf_reduced)
# Eliminate SNPs where ALT alleles are not A,C,T or G
vcf_reduced = vcf_reduced[ALT %in% c("A","C","T","G") ]

vcf_reduced$C_POS_REF_ALT = paste(vcf_reduced$CHROM,vcf_reduced$POS,vcf_reduced$REF,vcf_reduced$ALT,
                              sep = "_")
vcf_reduced = unique(vcf_reduced, by="C_POS_REF_ALT")

message("...Calculating minor allele dosages...")

ma_dosages = process_gt_dosages_ref_alt(vcf)

if(nrow(ma_dosages)<nrow(vcf_reduced)){
  vcf_reduced = vcf_reduced[vcf_reduced$C_POS_REF_ALT %in% ma_dosages$rn,]
}else{
  ma_dosages = ma_dosages[ma_dosages$rn %in% vcf_reduced$C_POS_REF_ALT,]
  
}

fwrite(ma_dosages,genotype_minor_allele_dos_path)

# 2 min 100k variants
# 15 min 1M variants
