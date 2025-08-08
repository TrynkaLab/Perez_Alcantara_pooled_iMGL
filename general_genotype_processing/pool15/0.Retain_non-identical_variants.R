# Filter VCF to only the variants that are different across all donors
#Sys.setenv(LD_LIBRARY_PATH="/software/teamtrynka/conda/otar2065/lib:/usr/lib/oracle/12.1/client64/lib:/software/lsf-farm5/10.1/linux3.10-glibc2.17-x86_64/lib")
#dyn.load("/software/teamtrynka/conda/otar2065/lib/libiconv.so")
#print(is.loaded("libiconv"))
.libPaths(c("/software/teamtrynka/conda/otar2065/lib/R/library","/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/"))

library(vcfR)
library(data.table)
library(tidyverse)
options(stringsAsFactors = FALSE)
source("functions.R")
gc()


# For server
args = commandArgs(trailingOnly=TRUE)


if (length(args)<2) {
  stop("You need to detail input and output file path.n", call.=FALSE)
} else if (length(args)==2) {
  vcf_path = args[1]
  output_path = args[2]
  
}


message("...Reading VCF...")


vcf = vcfR::read.vcfR(vcf_path)

message("...Calculating minor allele dosages...")

gt2 = process_gt_dosages(vcf)

# Get rows of non-identical alleles

cols <-c(setdiff(colnames(gt2),"rn"))

non_identical_rows = apply(gt2[,..cols], 1,function(x) length(unique(x))>1)

original_length = nrow(gt2)
gt2 = gt2[non_identical_rows,]
non_identical_genotype_length = nrow(gt2)

## subset vcf by identifying remaining snps
ids = paste0(vcf@fix[,"CHROM"],"_",vcf@fix[,"POS"],"_",vcf@fix[,"REF"],"_",vcf@fix[,"ALT"])
tosubset = ids %in% gt2$rn
vcf@gt = vcf@gt[tosubset,]
vcf@fix = vcf@fix[tosubset,]
vcf@meta = append(vcf@meta,"##R_custom_script=0.retain_non-identical_variants.R")

# subset original vcf to those rows and save
write.table(data.frame(original_length = original_length, non_identical_genotype_length = non_identical_genotype_length),
              file = paste0(vcf_path,"_Nvariants_retained.txt"),quote = F, row.names = F)
write.vcf(vcf,output_path)
