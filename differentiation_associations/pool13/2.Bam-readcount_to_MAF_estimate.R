# Convert bam-readcount output to minor allele frequency of SNP (b vector)
# Input: bam-readcount output and genotype vcf
# Output: minor allele frequency estimate (b vector)
.libPaths(c("/software/teamtrynka/conda/otar2065/lib/R/library","/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/"))

library(vcfR)
library(data.table)
library(tidyverse)
options(stringsAsFactors = FALSE)


# For server
args = commandArgs(trailingOnly=TRUE)


if (length(args)<4) {
  stop("You need to detail at least 2 input and 2 output file paths.n", call.=FALSE)
} else if (length(args)==4) {
  bam_readcount_path = args[1]
  vcf_path = args[2]
  b_estimate_path = args[3]
  genotype_minor_allele_dos_path = args[4]
}else if(length(args)==5){
  bam_readcount_path = args[1]
  vcf_path = args[2]
  b_estimate_path = args[3]
  genotype_minor_allele_dos_path = args[4]
  variants_retained_path = args[5]
}

# Reading in only 10 columns - consider selecting all if I look into non-ACGT single bp minor alleles
bam_readcount= data.table::fread( cmd = paste("zcat", bam_readcount_path,"| awk -F '\t' '{print $1 , $2 , $3 , $4, $5, $6, $7, $8, $9, $10}'"))

colnames(bam_readcount) = c("CHROM","POS","REF","TOTAL_READS","DEL","A","C","G","T","N")
bam_readcount[["REF"]] = toupper(bam_readcount[["REF"]]) 

vcf = vcfR::read.vcfR(vcf_path)
vcf_reduced = as.data.table(vcf@fix[,1:5])
original_vcf_rows = nrow(vcf_reduced)
# Eliminate SNPs where ALT alleles are not A,C,T or G
vcf_reduced = vcf_reduced[ALT %in% c("A","C","T","G") ]

vcf_reduced$C_POS_REF = paste(vcf_reduced$CHROM,vcf_reduced$POS,vcf_reduced$REF,sep = "_")

bam_readcount$C_POS_REF = paste(bam_readcount$CHROM,bam_readcount$POS,bam_readcount$REF,sep = "_")

message("...Removing duplicates and subsetting data tables to common variants...")

# Remove duplicates - there are a few in the bam-readcount file (why?)
bam_readcount = unique(bam_readcount, by="C_POS_REF")
vcf_reduced = unique(vcf_reduced, by="C_POS_REF")

# Subset datasets to common SNPs for comparison
bam_readcount = subset(bam_readcount, C_POS_REF %in% vcf_reduced$C_POS_REF )
vcf_reduced = subset(vcf_reduced, C_POS_REF %in% bam_readcount$C_POS_REF )

# Change into minor allele dosages and subset vcf to those that remain
process_gt_dosages = function(vcf){
  # Excluding indels (subsetting only to SNPs) - we may want to keep this info for calculating the proportions?
  # vcf2 <- extract.indels(vcf)
  
  ## Extracting GT
  gt <- as.data.table(extract.gt(vcf))
  # first make new names from chromosome, position, ref and alt (so there are no duplicates)
  names = paste0(vcf@fix[,"CHROM"],"_",vcf@fix[,"POS"],"_",vcf@fix[,"REF"])
  # There are rows for same SNP but different minor alleles: those will be removed with the previous line (does not include "ALT")
  # there are still some duplicates so remove those
  gt=gt[!duplicated(names),]
  message("Removing ",sum(duplicated(names)), " duplicated rows from vcf")
  
  ## removing those from the names
  names=names[!duplicated(names)]
  gt[, ("rn") := names] 
  
  # setkey(gt,rn) # this sorts the data.table, be careful
  object.size(gt)
  
  ### Some meanings: ################
  # / : genotype unphased (e.g. 0/0)
  #  | : genotype phased (e.g. 0|0)
  # Phased data are ordered along one chromosome and so from these data you know the haplotype
  # i.e. the variants are known to be on the same homologous chromosome because some reads were found to carry both variants
  # Unphased data are simply the genotypes without regard to which one of the pair of chromosomes holds that allele
  ################
  
  
  ## converting GT to minor allele score
  
  gt2 =  gt %>%
    dplyr::mutate_all(funs(str_replace_all(., "0\\|0", "0")))
  
  gt2 =  gt2 %>%
    dplyr::mutate_all(funs(str_replace_all(., "0\\/0", "0")))
  
  gt2 = gt2 %>%
    dplyr::mutate_all(funs(str_replace_all(., "0\\|1", "0.5")))
  
  gt2 = gt2 %>%
    dplyr::mutate_all(funs(str_replace_all(., "1\\|0", "0.5")))
  
  gt2 =  gt2 %>%
    dplyr::mutate_all(funs(str_replace_all(., "0\\/1", "0.5")))
  gt2 =  gt2 %>%
    dplyr::mutate_all(funs(str_replace_all(., "1\\/0", "0.5")))
  
  gt2 = gt2 %>%
    dplyr::mutate_all(funs(str_replace_all(., "1\\|1", "1")))
  gt2 = gt2 %>%
    dplyr::mutate_all(funs(str_replace_all(., "1\\/1", "1")))
  
  # change donor columns to numeric
  cols = setdiff(colnames(gt2),"rn")
  gt2[ ,(cols) := lapply(.SD, as.numeric),.SDcols = cols] 
  
  # hist(gt2)
  ## Some donors have NAs at some SNP positions
  anyNA(gt2)
  gt2[which(rowSums(is.na(gt2)) != 0)[1], ]
  gt[which(rowSums(is.na(gt2)) != 0)[1], ]
  
  # Removing SNPs with NAs in at least one donor - do I need to do this? Can I ignore the NAs somehow?
  gt2 = gt2[which(rowSums(is.na(gt2)) == 0), ]
  return(gt2)
}

message("...Calculating minor allele dosages...")

ma_dosages = process_gt_dosages(vcf)

if(nrow(ma_dosages)<nrow(vcf_reduced)){
  vcf_reduced = vcf_reduced[vcf_reduced$C_POS_REF %in% ma_dosages$rn,]
}else{
  ma_dosages = ma_dosages[ma_dosages$rn %in% vcf_reduced$C_POS_REF,]
  
}


# Are bam_readcount and vcf_reduced of same length and in the same order?
 if(!identical(vcf_reduced$C_POS_REF,bam_readcount$C_POS_REF)) stop("VCF and bam-readcount file are not identical after filtering.\n")
 if(!identical(ma_dosages$rn,bam_readcount$C_POS_REF)) stop("Dosage file and bam-readcount file are not identical after filtering.\n")


estimate_b_from_bam_readcount = function(bam_readcount,vcf){
  
  count_dt = data.table(
    total_reads = rep(0, nrow(bam_readcount)),
    A = rep(0, nrow(bam_readcount)),
    C = rep(0, nrow(bam_readcount)),
    `T` = rep(0, nrow(bam_readcount)),
    G = rep(0, nrow(bam_readcount)),
    b_estimate = as.numeric(rep(NA, nrow(bam_readcount)))
  )
  for(column in c("A","C","T","G")){
    sub_match = gsub("^([^:]+:[^:]+).*", "\\1", unlist(bam_readcount[[column]]))
    count_dt[[column]] = as.numeric(gsub(".*:","",sub_match))
    
  }
  count_dt[["total_reads"]] = bam_readcount[["TOTAL_READS"]]
  count_dt[["ALT"]] = vcf[["ALT"]]
  
  for(base in c("A","C", "G", "T")){
    count_dt[count_dt$ALT %in% base, "b_estimate"] = as.double(unlist(count_dt[count_dt$ALT %in% base,..base]) / unlist(count_dt[count_dt$ALT %in% base, "total_reads"]))
    
  }
  
  
  
  return(count_dt)

}

count_dt = estimate_b_from_bam_readcount(bam_readcount = bam_readcount, vcf = vcf_reduced)
message("The proportion of rows with b_estimate from the input vcf was ", 
        sum(table(count_dt$b_estimate)) / original_vcf_rows, " \n")

if(ncol(count_dt)!=7) stop("Error: The output file does not have the expected number of columns.\n")
if(sum(is.na(count_dt$b_estimate)) == length(count_dt$b_estimate)) stop("Error: All the minor allele frequency estimates are NA.\n")

message("Saving minor allele frequency estimates...\n")
fwrite(count_dt, file = b_estimate_path)

## Check again that minor allele dosage file and b_estimate file have same snps in same order

ma_dosages = ma_dosages[ma_dosages$rn %in% vcf_reduced$C_POS_REF,]
if(!identical(ma_dosages$rn,bam_readcount$C_POS_REF)) stop("Error: Dosage file and bam-readcount file are not identical at the end of script.\n")



if(exists("variants_retained_path")){
write.table(data.frame(original_length = original_vcf_rows, final_length = nrow(ma_dosages)),
            file = variants_retained_path,quote = F, row.names = F, col.names = T)
}
fwrite(ma_dosages,genotype_minor_allele_dos_path)

