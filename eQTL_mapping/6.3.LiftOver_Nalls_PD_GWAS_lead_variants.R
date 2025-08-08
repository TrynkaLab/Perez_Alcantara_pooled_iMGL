# liftover Nalls PD GWAS lead variants
.libPaths(c('/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/',"/software/teamtrynka/ma23/R4.1/libs",
            "/software/R-4.1.0/lib/R/library"))
library(GenomicRanges)
library(rtracklayer)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

if (length(args)<3) {
  stop("At least 3 arguments must be supplied", call. = FALSE)
}

GWAS_variants_path = as.character(args[1]) # file path of the lead variants from selected GWAS trait/study (GRCh37)
liftover_file_path = as.character(args[2]) # file path of the lead variants from selected GWAS trait/study (GRCh37)
output_path = as.character(args[3]) # file path of the lifted GWAS coordinates (GRCh38)

# to test
GWAS_variants_path = as.character("/lustre/scratch123/hgi/teams/trynka/resources/summary_statistics/public/ieu-b-7_PD_Nalls_2019/Nalls_2019_gwas_risk_variants_GRCh37.csv")
liftover_file_path = as.character("/lustre/scratch123/hgi/projects/otar2065/resources/for_liftOver/chain/hg19ToHg38.over.chain") # file path of the lead variants from selected GWAS trait/study (GRCh37)
output_path = as.character("/lustre/scratch123/hgi/teams/trynka/resources/summary_statistics/public/ieu-b-7_PD_Nalls_2019/Nalls_2019_gwas_risk_variants_GRCh38_loci.tsv")

gwas = read.csv(GWAS_variants_path) %>%
  dplyr::rename(base_pair = BP,
                chromosome = chr)


# Load the chain file for GRCh37 to GRCh38 conversion
chainFile = rtracklayer::import.chain(liftover_file_path)

# Define the input coordinates in GRCh37
# fix to have some kind of identifier for SNPs

inputCoordinates = gwas %>%
  dplyr::mutate(seqnames = paste0("chr",chromosome),
                seqinfo = paste(chromosome,base_pair,sep = "_"),
                start = base_pair, 
                end =  base_pair,
                variant_id = paste(chromosome,base_pair,REF,ALT,sep = "_")) %>%
dplyr::select(seqnames,start,end,seqinfo,REF,ALT,variant_id,Gene) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# Convert coordinates from GRCh37 to GRCh38
outputCoordinates = rtracklayer::liftOver(as(inputCoordinates, "GRanges"), chainFile)  %>%
  as_tibble() %>%
  dplyr::mutate(chr = as.numeric(gsub("chr","",seqnames)),
                snp_pos=as.numeric(start)) %>%
  dplyr::mutate(variant_id = paste(chr,snp_pos,inputCoordinates$REF,inputCoordinates$ALT,sep = "_")) %>%
  dplyr::rename(locus_name = Gene) %>%
  dplyr::select(chr ,  snp_pos, variant_id,      locus_name)


outputCoordinates %>%
  readr::write_tsv(.,output_path)
R.utils::gzip(output_path)
