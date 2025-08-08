# to liftOver Nalls PD GWAS
.libPaths(c('/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/',"/software/teamtrynka/ma23/R4.1/libs",
            "/software/R-4.1.0/lib/R/library"))
library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
outDir = "../../../resources/"

args = commandArgs(trailingOnly=TRUE)

if (length(args)<3) {
  stop("At least 3 arguments must be supplied", call. = FALSE)
}

GWAS_variants_path = as.character(args[1]) # file path of the lead variants from selected GWAS trait/study (GRCh37)
liftover_file_path = as.character(args[2]) # file path of the lead variants from selected GWAS trait/study (GRCh37)
output_path = as.character(args[3]) # file path of the lifted GWAS coordinates (GRCh38)

# to test
# GWAS_variants_path = as.character("../../../resources/harmonised_GWAS/GCST009325_Nalls_et_al_2019_PD/GCST009325.tsv")  # file path of the lead variants from selected GWAS trait/study (GRCh37)
# liftover_file_path = as.character("/lustre/scratch123/hgi/projects/otar2065/resources/for_liftOver/chain/hg19ToHg38.over.chain") # file path of the chain file for liftover
# output_path = as.character("../../../resources/harmonised_GWAS/GCST009325_Nalls_et_al_2019_PD/GCST009325_GRCh38_manually_formated.tsv")

gwas = readr::read_tsv(GWAS_variants_path) 


# Load the chain file for GRCh37 to GRCh38 conversion
chainFile = rtracklayer::import.chain(liftover_file_path)

# Define the input coordinates in GRCh37
# fix to have some kind of identifier for SNPs

inputCoordinates = gwas %>%
  dplyr::mutate(seqnames = paste0("chr",chromosome),
                seqinfo = paste(chromosome,base_pair_location,sep = "_"),
                start = base_pair_location, 
                end =  base_pair_location,
                variant_id = paste(chromosome,base_pair,other_allele,effect_allele,sep = "_")) %>%
  dplyr::rename(pval_nominal = p_value) %>%
  dplyr::select(seqnames,start,end,seqinfo,other_allele,effect_allele,beta,
                standard_error,pval_nominal,effect_allele_frequency,N_cases,N_controls,variant_id) %>%
GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# Convert coordinates from GRCh37 to GRCh38
outputCoordinates = rtracklayer::liftOver(as(inputCoordinates, "GRanges"), chainFile)  %>%
  as_tibble() %>%
  dplyr::mutate(chr = as.numeric(gsub("chr","",seqnames)),
                varbeta = as.numeric(standard_error)^2*(N_cases + N_controls),
                pos=as.numeric(start))



outputCoordinates %>%
  readr::write_tsv(.,output_path)
gzip(output_path)
