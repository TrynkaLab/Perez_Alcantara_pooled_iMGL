# reformat LD info file
.libPaths(c('/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/',"/software/teamtrynka/ma23/R4.1/libs",
            "/software/R-4.1.0/lib/R/library"))
library(tidyverse)
library(stringr)

outDir = "../../data/results/4.Inspect_eQTL_results"
dir.create(outDir)

ld = read_table("../../data/genotype/plink_genotypes/all_pools.genotype.MAF05.filtered.tags.list")

tensorqtl_all = read_csv(paste0(outDir,"/tensorQTL_variant_gene_60PCs.csv"))
tensorqtl_all = tensorqtl_all %>% 
  dplyr::filter(qval < 0.05) 

# every row is a variant ID at a position that showed significant effects in any 
# tensorQTL results
# the TAGS column shows all SNPs that are tagged (in LD with r2 above 0.8)

all_tags_column = unique(unlist(lapply(ld$TAGS,function(x) str_split(string = x,
                                                                     pattern = "\\|"))))
# check if significant variant IDs are also in tagged column 
# meaning that the script worked as intended AND
# there are significant SNPs that are not independent from each other

table(ld$SNP %in% all_tags_column)
# 1699 SNPs from 15509 (~40%) are in r2>0.8 with each other

# adding variant info to ld file ########
# as in 5.subset_variant_id_for_plink_LD.R
# Load the bim file
bim_file = readr::read_table("../../data/genotype/plink_genotypes/all_pools.genotype.MAF05.filtered.bim", 
                             col_names = c("chr", "variant_id", "gd", "pos", "a1", "a2"),
                             col_types = cols(
                               chr = col_character(), # there's an X
                               variant_id = col_character(),
                               gd = col_double(),
                               pos = col_double(),
                               a1 = col_character(),
                               a2 = col_character()
                             ))


#########
#### filtering tensor variants #######

reformated_bim = bim_file %>%
  mutate(tensor_variant_id = paste(chr,pos,a2,a1, sep = "_"), # tensor variants with respect to bim file are chr_pos_A2_A1 (weird)
         wrong_id = paste(chr,pos,a1,a2, sep = "_")) # to check

reformated_filtered_bim = reformated_bim %>%
  dplyr::filter(tensor_variant_id %in% unique(tensorqtl_all$variant_id)) %>%
  dplyr::rename(rsID=variant_id,variant_id = tensor_variant_id ) %>%
  dplyr::select(rsID,variant_id)
rm(reformated_bim)
gc()
# add variant info

tensorqtl_all = tensorqtl_all %>%
  left_join(.,reformated_filtered_bim,by = join_by(variant_id))
anyNA(tensorqtl_all$rsID) # some variants without rsID?
table(is.na(tensorqtl_all$rsID)) # quite a few - investigate

table(tensorqtl_all$rsID %in% ld$SNP) # those without rsID are not present in LD file
table(tensorqtl_all$rsID %in% all_tags_column) # 3644 in LD with others

# add that column to tensor results table
tensorqtl_all = tensorqtl_all%>%
  left_join(.,ld[c("SNP","TAGS")], by=join_by(rsID == SNP))

# Leave only the values from TAGS that are found in rsID
tensorqtl_all = tensorqtl_all %>%
  mutate(TAGS = str_split(TAGS, "\\|")) %>%
  mutate(TAGS = map(TAGS, ~ intersect(., rsID))) %>%
  mutate(TAGS = map_chr(TAGS, ~ str_c(.x, collapse = ",")))

tensorqtl_all[c("rsID","TAGS")]

table(tensorqtl_all$rsID %in% unique(unlist(lapply(tensorqtl_all$TAGS,
                                                   function(x) str_split(string = x, pattern = "\\,")))))

# we're only interested in those variants in LD that are associated with the expression of the same gene
# add to each TAGS column that is not empty the rsID (for easier filtering - best to create new column)
tensorqtl_all = tensorqtl_all %>%
  dplyr::mutate(tags_for_filtering = if_else(str_count(TAGS)==0,
                                             true = NA,
                                             false = paste0(TAGS,",",rsID) )) 
# grep in rsIDs or tags, then assign id from largest effect
tensorqtl_all %>%
  separate_rows(tags_for_filtering, sep = ",") %>% #reorder in alphabetical order
  arrange(tags_for_filtering) %>%
  group_by(row_number()) %>%
  summarise(tags_for_filtering = toString(tags_for_filtering)) %>%
  ungroup()

#count duplicates
na.omit(tensorqtl_all$tags_for_filtering)

#save the result table with the "raw" eQTL results and the variants in LD, and then 
# present results figures e.g. an expression-genotype boxplot of that shared eGene-variant for all conditions 
# picking one variant name (the "top" one) and show the different slopes