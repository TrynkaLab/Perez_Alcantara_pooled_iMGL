# to liftOver Jeremy AD GWAS lead SNP coordinates
.libPaths(c('/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/',"/software/teamtrynka/ma23/R4.1/libs",
            "/software/R-4.1.0/lib/R/library"))
library(GenomicRanges)
library(rtracklayer)

outDir = "../../../resources/"

args = commandArgs(trailingOnly=TRUE)

if (length(args)<3) {
  stop("At least 3 arguments must be supplied", call. = FALSE)
}

GWAS_lead_variants_path = as.character(args[1]) # file path of the lead variants from selected GWAS trait/study (GRCh37)
liftover_file_path = as.character(args[2]) # file path of the lead variants from selected GWAS trait/study (GRCh37)
output_path = as.character(args[3]) # file path of the lifted GWAS coordinates (GRCh38)

# to test
GWAS_lead_variants_path = as.character("../../../resources/Jeremy_medrXiv_AD_loci.txt")
liftover_file_path = as.character("/lustre/scratch123/hgi/projects/otar2065/resources/for_liftOver/chain/hg19ToHg38.over.chain") # file path of the lead variants from selected GWAS trait/study (GRCh37)
output_path = as.character("../../../resources/Jeremy_medrXiv_AD_loci_GRCh38.txt")

gwas_lead = read_tsv(GWAS_lead_variants_path) 


# Load the chain file for GRCh37 to GRCh38 conversion
chainFile = rtracklayer::import.chain(liftover_file_path)

# Define the input coordinates in GRCh37
inputCoordinates = GRanges(
  seqnames = paste0("chr",gwas_lead$Chr),
  ranges = IRanges(start = gwas_lead$SNP_pos, end =  gwas_lead$SNP_pos)
)

# Convert coordinates from GRCh37 to GRCh38
outputCoordinates = rtracklayer::liftOver(as(inputCoordinates, "GRanges"), chainFile) %>%
  as_tibble() %>%
  mutate(chr = gsub("chr","",seqnames))

### I checked they are in the same order, but be careful

gwas_lead %>%
  mutate(SNP_pos = outputCoordinates$start) %>%
  mutate(variant_id = paste(Chr,SNP_pos,
                            str_split_fixed(SNP_Lead, "_", n = 3)[, 2],
                            str_split_fixed(SNP_Lead, "_", n = 3)[, 3],
                            sep="_")) %>%
  select(Chr,SNP_pos,variant_id,Locus_name) %>%
  write_tsv(.,output_path)

