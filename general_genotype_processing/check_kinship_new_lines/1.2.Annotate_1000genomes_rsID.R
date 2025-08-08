# annotate 1000 genomes file wirh rsIDs like the hipsci file
.libPaths(c('/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/',"/software/teamtrynka/ma23/R4.1/libs",
            "/software/R-4.1.0/lib/R/library"))
library(tidyverse)
library(stringr)
outdir = "../../data/check_kinship_new_lines/genotype/plink_genotypes/"
thousand = readr::read_table(paste0(outdir,
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
hipsci =readr::read_table(paste0(outdir,
                            "tocheck.genotype.MAF05.bim"), 
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

table(hipsci$new_variant_id %in% thousand$variant_id)       

thousand_annotated = thousand %>%
  left_join(hipsci,by =join_by(variant_id == new_variant_id) ) %>%
  mutate(variant_id = case_when(is.na(variant_id.y) == FALSE ~ variant_id.y, # for same positions, change variant ID to that of the hispci
                                is.na(variant_id.y) == TRUE ~ variant_id))

thousand_annotated %>%
  dplyr::filter(!is.na(variant_id.y))

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

table(missing$variant_id %in% thousand_annotated$variant_id)
table(missing$variant_id %in% hipsci$variant_id)

missing$variant_id %>%
  write_lines(.,
              file  = paste0(outdir,"1000G_maf05perc_biallelic_nonmissing_nochr_withrsid_checked_annotated.missingsnp"))           
