# filter VCF/BED file to remove variants that are not covered by at least 
# 5 donors in at least 2 genotype categories
.libPaths(c('/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/',"/software/teamtrynka/ma23/R4.1/libs",
            "/software/R-4.1.0/lib/R/library"))
library(stringr)
library(tidyverse)
library(data.table)
library(genio)
library(bedr)
source("./functions.R")
directory = "../../data/genotype/"


# filter bed file directly? try

plink_data = genio::read_plink(paste0(directory,"/plink_genotypes/all_pools.genotype.MAF05"))
plink_data$X[1:10, 1:10] # 0 is homozygote ref, 1 homozygote alt, 2 het
plink_data$bim[1:10,] # some are swapped with respect to the VCF, as I knew already
# e.g. rs140337953, because PLINK places the most infrequent allele in the population as the minor allele
bed_tibble = plink_data$X %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "snp") %>%
  pivot_longer(cols=-snp,names_to = "line",values_to = "genotype")

# test = bed_tibble[1:10000,]

rm(plink_data)
gc()
bed_tibble = bed_tibble %>%
  group_by(snp) %>%
  mutate(genotype = as.factor(genotype))%>%
  summarize(max_count = sum(genotype == names(which.max(table(genotype)))),
            min_count = sum(genotype == names(which.min(table(genotype)))))

# test %>%
#   group_by(snp) %>%
#   mutate(genotype = as.factor(genotype)) %>%
#   summarize(max_count = sum(genotype == names(which.max(table(genotype)))),
#  count = sum(genotype == names(which.min(table(genotype)))))
# 
# test2 = bed_tibble %>% dplyr::filter(snp =="rs114111569")
# 
# test2 %>% group_by(rowname) %>%
#   mutate(genotype = as.factor(genotype)) %>%
#   summarize(max_count = sum(genotype == names(which.max(table(genotype)))),
#             count = sum(genotype == names(which.min(table(genotype)))))

snps_to_retain = bed_tibble %>%
  dplyr::filter(min_count>5) %>%
  .$snp

rm(bed_tibble)
gc()

# subset
plink_data = genio::read_plink(paste0(directory,"/plink_genotypes/all_pools.genotype.MAF05"))
plink_data$X = plink_data$X[rownames(plink_data$X) %in% snps_to_retain,]
plink_data$bim = plink_data$bim[plink_data$bim$id %in% snps_to_retain,]

genio::write_plink(paste0(directory,"/plink_genotypes/all_pools.genotype.MAF05.filtered"), 
                   plink_data$X , plink_data$bim, plink_data$fam)
