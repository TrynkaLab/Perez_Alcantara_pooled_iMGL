# annotate 1000 genomes file wirh IDs like the hipsci file
library(tidyverse)
library(stringr)
outdir = "../../data/genotype/plink_genotypes/"
thousand = read_table(paste0(outdir,
                             "1000G_maf05perc_biallelic_nonmissing_nochr_withrsid_checked.bim"), 
                      col_names = c("chr", "variant_id", "gd", "pos", "a1", "a2"),
                      col_types = cols(
                        chr = col_character(), 
                        variant_id = col_character(),
                        gd = col_double(),
                        pos = col_double(),
                        a1 = col_character(),
                        a2 = col_character()
                      )
)
hipsci =read_table(paste0(outdir,
                            "all_pools.all_donors.genotype.MAF05.bim"), 
                   col_names = c("chr", "variant_id", "gd", "pos", "a1", "a2"),
                   col_types = cols(
                     chr = col_character(), # there's an X
                     variant_id = col_character(),
                     gd = col_double(),
                     pos = col_double(),
                     a1 = col_character(),
                     a2 = col_character()
                   ))

hipsci = hipsci %>%
  mutate(new_variant_id = paste(chr,pos,a2,a1, sep = ":"))
thousand = thousand %>%
  mutate(new_variant_id = paste(chr,pos,a2,a1, sep = ":"))

table(hipsci$new_variant_id %in% thousand$new_variant_id)       

thousand_annotated = thousand %>%
  left_join(hipsci,by =join_by(new_variant_id) ) %>%
  dplyr::rename(variant_id = variant_id.x) %>%
  mutate(variant_id = case_when(is.na(variant_id.y) == FALSE ~ variant_id.y,
                                is.na(variant_id.y) == TRUE ~ variant_id))

thousand_annotated %>%
  dplyr::filter(!is.na(variant_id.y))

####### remove duplicated positions (multiallelic)


thousand_annotated %>%
  select(c(chr.x, variant_id,  gd.x,  pos.x, a1.x,  a2.x)) %>%
  write_delim(.,
               file  = paste0(outdir,"1000G_maf05perc_biallelic_nonmissing_nochr_withrsid_checked_annotated.bim"),
              col_names = FALSE,delim = "\t")
missingth = thousand_annotated %>%
  filter(is.na(variant_id.y)) %>%
  select(c(variant_id))
missing = hipsci %>%
  filter(!variant_id %in% thousand_annotated$variant_id) %>%
  select(c(variant_id))

missing = unique(rbind(missing,missingth))

missing$variant_id %>%
  write_lines(.,
              file  = paste0(outdir,"1000G_maf05perc_biallelic_nonmissing_nochr_withrsid_checked_annotated.missingsnp"))           
