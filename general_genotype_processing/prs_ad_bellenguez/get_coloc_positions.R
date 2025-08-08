# get positions from coloc
library(tidyverse)
library(GenomicRanges)

options(stringsAsFactors=F)

args = commandArgs(trailingOnly = TRUE)
coloc_results_path = args[1]
gwas_positions_path = args[2]
gwas_path=args[3]
  
# for testing
#coloc_results_path="/lustre/scratch123/hgi/projects/otar2065/OTAR2065_sc_eQTL/data/results/8.colocalisation_analysis/coloc_results/all_GWAS_colocalisations.csv"
#gwas_positions_path="/lustre/scratch123/hgi/projects/otar2065/OTAR2065_sc_eQTL/data/results/8.colocalisation_analysis/coloc_results/coloc_positions_for_PRSice/AD_GWAS_coloc_pos_LPS.csv"
# gwas_path="/lustre/scratch123/hgi/teams/trynka/resources/summary_statistics/public/GCST90027158/harmonised/35379992-GCST90027158-MONDO_0004975.h.tsv.gz"
dir.create(dirname(gwas_positions_path),recursive = TRUE)

coloc_results = readr::read_csv(coloc_results_path) %>%
  dplyr::filter(GWAS == "AD" & PP_H4 >= 0.70) %>%
  dplyr::select(treatment,chr,snp_pos,external_data_Locus, eGene,locus_name) %>%
  dplyr::mutate(pos_low = snp_pos - 250000, pos_high = snp_pos + 250000) %>%
  dplyr::group_split(treatment)

coloc_results = purrr::map(coloc_results, distinct)
coloc_results = purrr::map(coloc_results, ~ .x %>% dplyr::distinct(locus_name, .keep_all = TRUE) %>%
  dplyr::arrange(chr,pos_low))

gwas = readr::read_tsv(gwas_path)
for(i in seq_along(coloc_results)){

# Convert the data to GenomicRanges
gr = GRanges(
  seqnames = Rle(as.character(coloc_results[[i]]$chr)),  # Use 'chr' as sequence names
  ranges = IRanges(start = coloc_results[[i]]$pos_low, 
                   end = coloc_results[[i]]$pos_high,
                   treatment = coloc_results[[i]]$treatment,
                   external_data_Locus = coloc_results[[i]]$external_data_Locus ,
                   eGene= coloc_results[[i]]$eGene,
                   locus_name = coloc_results[[i]]$locus_name)
)

# Reduce the ranges to merge overlapping/adjacent ones (so the regions may be larger than 500Kb)
gr_merged = reduce(gr)


# Convert the merged GRanges object back into a tibble
bedtools_results = tibble(
  CHROM = as.integer(as.character(seqnames(gr_merged))),
  BEG = start(gr_merged),
  END = end(gr_merged)
)

# subset GWAS file to coloc positions
reduced_gwas = list()
for(n in 1:nrow(bedtools_results)){
  reduced_gwas[[n]] = gwas %>%
    dplyr::filter(hm_chrom == bedtools_results$CHROM[n] & hm_pos>=bedtools_results$BEG[n] & hm_pos<=bedtools_results$END[n])
}

reduced_gwas = do.call("rbind",reduced_gwas) %>%
  dplyr::distinct()

  name = coloc_results[[i]] %>%
    pull(treatment) %>%
    unique() %>%
    paste(collapse = "_")
  
  readr::write_delim(coloc_results[[i]], delim = "\t",col_names = FALSE,
                     paste0(dirname(gwas_positions_path),"/AD_GWAS_coloc_pos_", name, ".tsv"))
  
  readr::write_tsv(reduced_gwas,
  paste0("../../data/prs_ad_bellenguez/AD_coloc_subset_",name,".tsv.gz"))


}
