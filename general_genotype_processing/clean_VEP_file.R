
# clean VEP file from imputed genotype
library(tidyverse)
library(vcfR)

vep = list()
for(chr in c(1:22,"X")){
  vep[[chr]] = readr::read_delim(paste0("../../data/consequence_calling/missense_dirty.",chr,".txt"),delim = "\t")

}

vep = do.call("rbind",vep)

table(vep$Consequence)
table(vep$PolyPhen)

vep = vep %>%
  dplyr::mutate(CADD_PHRED = case_when(CADD_PHRED == "." ~ NA,
                                       .default = as.double(CADD_PHRED))) %>%
  # defining as deleterious - union of LoF and missense pathogenic mutations
  # LoF -  frameshift, stop-gain, transcript_ablation, splice acceptor or splice donor variants
  # see https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
  # missense pathogenic - annotated as missense or start loss with a CADD Phred score cutoff > 15 (as Pau did)
  dplyr::mutate(category = case_when((str_detect(Consequence,pattern="missense_variant")| str_detect(Consequence,pattern="start_lost")) & CADD_PHRED >=15 ~ "deleterious",
                                     (str_detect(Consequence,pattern="frameshift")|
                                        str_detect(Consequence,pattern="stop_gained")|
                                        str_detect(Consequence,pattern="transcript_ablation")|
                                        str_detect(Consequence,pattern="splice_acceptor_variant")|
                                        str_detect(Consequence,pattern="splice_donor_variant")) ~ "deleterious",
                                     .default = "missense_non_deleterious"))

table(vep$category) # 11804 deleterious
length(unique(vep$Gene)) # 15071 genes

vep %>%
  dplyr::filter(SYMBOL == "SPG11") %>%
  dplyr::filter(category == "deleterious")
# just one variant detected - Pau saw up to 5 (somatic) deleterious mutations on this gene
# because they were somatic - detected on WES but are not SNPs
# so any VEP called from imputed genotypes will only test the effects of known 
# variants that can be imputed from 1KGenomes etc.
# need to do VEP on WES files from /lustre/scratch123/hgi/projects/hipsci/releases/data/exomeseq
# and /lustre/scratch123/hgi/projects/hipsci/releases/data/deep_exomeseq

vep %>%
  dplyr::filter(category == "deleterious") %>%
  dplyr::group_by(Gene) %>%
    dplyr::reframe(number_of_mutations = n(),SYMBOL = SYMBOL) %>%
  dplyr::distinct() %>%
  dplyr::arrange(by=desc(number_of_mutations)) 
# still there are genes with tons of mutations detected
# for very large genes it seems

# load genotype information

for(chr in c(1:22,"X")){
genotype = vcfR::read.vcfR(paste0("../../data/consequence_calling/microglia_samples.GRCh38.filtered.with_X.",chr,".vep.vcf.gz"))
tidy_gt = vcfR::extract_gt_tidy(genotype,alleles = FALSE)
fix = as.data.frame(genotype@fix[,1:5]) %>%
  tibble::rownames_to_column("Key") %>%
  dplyr::mutate(Key = as.numeric(Key)) %>%
  dplyr::left_join(tidy_gt[,c( "Key", "Indiv", "gt_GT")])
rm(tidy_gt)
gc()
fix = fix %>%
  # changing to number of ALT alleles
  # 0, 1 ,2
  # beware for X chromosome males will only have 0 or 1
  # and fix line names
  dplyr::mutate(gt_GT = case_when(gt_GT %in% c("0|1" , "1|0") ~ "1",
                                  gt_GT == "0|0"  ~ "0",
                                  gt_GT == "1|1"  ~ "2",
                                  .default = gt_GT),
                Indiv = case_when(str_detect(Indiv,"HPS") ~ str_split_i(Indiv,"-",2),
                                  .default = Indiv)) %>%
  dplyr::mutate(gt_GT = as.numeric(gt_GT),
                POS = as.double(POS)) %>%
# change to wide format
  tidyr::pivot_wider(id_cols = c( "Key", "CHROM", "POS" ,"ID", "REF", "ALT"),
                     names_from = "Indiv",values_from = "gt_GT") %>%
  dplyr::select(-Key)
  
vep = vep %>%
  dplyr::left_join(fix)
}

# save vep with genotype
write_csv(vep,"../../data/consequence_calling/deleterious_missense_VEP_from_imputed_genotypes.csv")
