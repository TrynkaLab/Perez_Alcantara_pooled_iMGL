# standardise kinship matrix for LMM and build correct dimensions
library(tidyverse)

kinship = read_tsv("../../data/kinship/all_pools.genotype.MAF01.hg38.king", col_names = FALSE)
id = read_tsv("../../data/kinship/all_pools.genotype.MAF01.hg38.king.id")
colnames(kinship) = id$`IID`
hist(kinship$aehn_22)
summary(kinship$aehn_22)
summary(unlist(kinship))
# no need to scale?

write_tsv(kinship,"../../data/kinship/kinship_for_LMM.tsv")
