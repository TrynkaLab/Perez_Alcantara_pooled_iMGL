
library(vcfR)
library(tidyverse)

## checking af before imputation
# http://samtools.github.io/bcftools/howtos/plugin.af-dist.html
af_preqc = read.table("../../data/all_pools_no_curn_iukl/concat.af.dist.txt") %>%
  dplyr::mutate(af = factor(paste(V2,V3,sep="-"),ordered = TRUE),
                category = "pre_qc") %>%
  dplyr::rename(freq = V4) %>%
  dplyr::select(freq,af,category)
  
af_posqc = read.table("../../data/all_pools_no_curn_iukl/toimpute.af.dist.txt") %>%
  dplyr::mutate(af = factor(paste(V2,V3,sep="-"),ordered = TRUE),
                category = "post_qc") %>%
  dplyr::rename(freq = V4) %>%
  dplyr::select(freq,af,category) %>%
  dplyr::bind_rows(af_preqc)

  
pdf("../../data/all_pools_no_curn_iukl/toimpute.af.dist.pdf", width = 5, height = 5)

af_posqc %>%
  ggplot(aes(af,y=freq,col = category)) +
  geom_point() +
  theme_minimal() +
  xlab("Genotype probability (given 1k Genomes allele freq.)")

dev.off()

# genotype probabilities are what they should be, even larger peak at high probability
# The program checks non-reference genotypes (i.e. not necessarily all markers are looked at), 
# calculates how likely it is to observe the 0/1 or 1/1 genotype given the allele frequency (the former is 2*AF*(1-AF), the latter is AF*AF), and then outputs the counts across the probability bins. 
# one wants to see few unlikely genotypes.

af = read.table("../../data/all_pools_no_curn_iukl/toimpute.nor2filter.af.dist.txt")
pdf("../../data/all_pools_no_curn_iukl/toimpute.nor2filter.af.dist.pdf", width = 5, height = 5)
af %>%
  dplyr::mutate(af = factor(paste(V2,V3,sep="-"),ordered = TRUE)) %>%
  dplyr::rename(freq = V4) %>%
  ggplot(aes(af,y=freq)) +
  geom_point() +
  theme_minimal()

dev.off()
# before filtering very high probabilities are not enriched
# obs because I filtered to perfect genotypes (INFO=1) so then 0.9-1 probs would be higher

# 1K_10K reference - checking chromosome 22
input_genotype_path = "../../data/all_pools_no_curn_iukl/1K_10K/19.vcf.gz"
vcf = vcfR::read.vcfR(input_genotype_path,nrows = 100000)
input_qual = vcf@fix[,8]
position = vcf@fix[,2]
# extracting INFO
input_qual_mat = stringr::str_split_fixed(string = input_qual,
                                  ";", 5)
# fixing issue with unequal number of elements in info pushing INFO to 5th position
input_qual = data.frame(V1=input_qual_mat[,4],
                        V2=input_qual_mat[,5]) %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(V1=case_when(V2!="" ~ V2,
                             .default = V1))
input_qual = stringr::str_split_i(string = input_qual$V1,
                                  "=", 2)

toplot = data.frame(tenk = as.double(input_qual),
                    position = position)


p1 = ggplot(toplot,aes(x=tenk)) + 
  geom_histogram( ) +
  theme_minimal() + 
  ggtitle("1KGenomes + UK10K imputation quality")
p1
 # looks OK

### comparing WGS positions for shared donors 
wgs = vcfR::read.vcfR("../../data/ungenotyped_donors/deepvariant_nextflow/AD_DV_output_hg19.vcf.gz")
wgs = cbind(wgs@fix,wgs@gt) %>%
  dplyr::as_tibble() %>%
  dplyr::filter(CHROM == 19) %>%
  tidyr::pivot_longer(cols =!c(CHROM,POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT),
                      names_to = "donor",values_to = "genotype") %>%
  dplyr::mutate(genotype = str_split(genotype, ":", simplify = TRUE)[,1],
                CHROM = as.numeric(CHROM),
                POS=as.numeric(POS))  %>%
  dplyr::filter(donor %in% c("aowh_2","kolf2_1s")) %>%
  dplyr::mutate(donor = case_when(donor=="kolf2_1s" ~ "kolf_2",
                                  .default = donor)) %>%
  dplyr::filter(genotype %in% c("0/0","0/1","1/0","1/1")) %>%
  dplyr::mutate(allele1 = as.numeric(str_split_i(genotype, "\\/", i=1)),
                allele2=as.numeric(str_split_i(genotype, "\\/", i=2))) %>%
  dplyr::mutate(category = "wgs")
  
vcf = cbind(vcf@fix,vcf@gt) %>%
  dplyr::as_tibble() %>%
  tidyr::pivot_longer(cols =!c(CHROM,POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT),
                      names_to = "donor",values_to = "genotype") %>%
  dplyr::mutate(genotype = str_split(genotype, ":", simplify = TRUE)[,1],
                CHROM = as.numeric(CHROM),
                POS=as.numeric(POS)) %>%
  dplyr::filter(donor %in% c("aowh_2","kolf2")) %>%
  dplyr::mutate(allele1 = as.numeric(str_split_i(genotype, "\\|", i=1)),
                allele2=as.numeric(str_split_i(genotype, "\\|", i=2))) %>%
  dplyr::mutate(ID = paste(CHROM,POS,REF,ALT,sep = "_")) %>%
  dplyr::mutate(category = "imputed_old") 

merged = rbind(wgs, vcf)
merged = merged %>%
  dplyr::filter(POS %in% c(vcf$POS)) %>%
  dplyr::arrange(POS,category)



# HRC reference - checking chromosome 22
input_genotype_path = "../../data/all_pools_no_curn_iukl/HRC/22.vcf.gz"
vcf = vcfR::read.vcfR(input_genotype_path)
input_qual = vcf@fix[,8]
position = vcf@fix[,2]
# extracting INFO
input_qual_mat = stringr::str_split_fixed(string = input_qual,
                                          ";", 5)
# fixing issue with unequal number of elements in info pushing INFO to 5th position
input_qual = data.frame(V1=input_qual_mat[,4],
                        V2=input_qual_mat[,5]) %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(V1=case_when(V2!="" ~ V2,
                             .default = V1))
input_qual = stringr::str_split_i(string = input_qual$V1,
                                  "=", 2)

toplot2 = data.frame(HRC = as.double(input_qual),
                    position = as.numeric(position))
summary(toplot2$HRC) # median imputation quality from HRC is lower, but in theory more accurate
summary(toplot$tenk)

p3 = ggplot(toplot2,aes(x=HRC)) + 
  geom_histogram( ) +
  theme_minimal() + 
  ggtitle("HRC imputation quality")

p4 = ggplot(toplot2,aes(x=position, y = HRC)) + 
  geom_point(alpha=0.6) + 
  theme_minimal() + 
  ggtitle("HRC quality along chr 22")

# some positions have poor quality
maf = vcfR::maf(vcf, element = 2)
problem_pos = toplot2[toplot2$HRC<0.2,"position"]
rows = vcf@fix[vcf@fix[,2] %in% problem_pos,]

chrom = vcfR::create.chromR(name="HRC", vcf=vcf, verbose=TRUE)

chrom = vcfR::masker(chrom, min_QUAL=0, min_DP=350, max_DP=650, min_MQ=59.5, max_MQ=60.5)
chrom = vcfR::proc.chromR(chrom, verbose = TRUE)

############# IPMAR ##########################
# 1K_10K reference - checking chromosome 22
input_genotype_path = "../../data/all_pools_no_curn_iukl/1K_10K/22.vcf.gz"
vcf = vcfR::read.vcfR(input_genotype_path)
input_qual = vcf@fix[,8]
position = vcf@fix[,2]
# extracting INFO
input_qual_mat = stringr::str_split_fixed(string = input_qual,
                                          ";", 5)
# fixing issue with unequal number of elements in info pushing INFO to 5th position
input_qual = data.frame(V1=input_qual_mat[,4],
                        V2=input_qual_mat[,5]) %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(V1=case_when(V2!="" ~ V2,
                             .default = V1))
input_qual = stringr::str_split_i(string = input_qual$V1,
                                  "=", 2)

toplot = data.frame(tenk = as.double(input_qual),
                    position = position)


p1 = ggplot(toplot,aes(x=tenk)) + 
  geom_histogram( ) +
  theme_minimal() + 
  ggtitle("1KGenomes + UK10K imputation quality")

### comparing shared donors from WGS genotype call (no imputation) to old file imputed
new_genotype = readr::read_tsv("../../data/all_pools_no_curn_iukl/merged_phased_imputed_UK10K_1KG_HRC_hg38_Bellenguez_checked.txt") %>%
  dplyr::select(c("CHROM","POS","REF","ALT","kolf_2","aowh_2"))%>%
  tidyr::pivot_longer(cols =!c(CHROM,POS, REF, ALT),
  names_to = "donor",values_to = "genotype") %>%
  dplyr::mutate(genotype = str_split(genotype, ":", simplify = TRUE)[,1]) %>%
  dplyr::mutate(allele1 = as.numeric(str_split_i(genotype, "\\|", i=1)),
                allele2=as.numeric(str_split_i(genotype, "\\|", i=2))) %>%
  dplyr::mutate(gtCount = allele1 + allele2) %>%
  dplyr::select(CHROM,POS,REF,ALT,donor,gtCount) %>%
  dplyr::mutate(ID = paste(CHROM,POS,REF,ALT,sep="_")) %>%
  dplyr::mutate(type = "new") %>%
  dplyr::relocate(ID)

imputation_info = new_genotype = readr::read_tsv("../../data/all_pools_no_curn_iukl/merged_phased_imputed_UK10K_1KG_HRC_hg38_Bellenguez_checked.txt") %>%
  dplyr::select(c("CHROM","POS","REF","ALT","INFO"))
  

wgs_no_imp = readr::read_tsv("/lustre/scratch123/hgi/projects/otar2065/hipsci_genotype_processing/data/ungenotyped_donors/deepvariant_nextflow/AD_Bellenguez_GRCh38.txt") %>%
  dplyr::mutate(CHROM = paste0("chr",CHROM)) %>%
  dplyr::rename(kolf_2 = kolf2_1s) %>%
  dplyr::select(c("CHROM","POS","REF","ALT","kolf_2","aowh_2")) %>%
  dplyr::mutate(ID = paste(CHROM,POS,REF,ALT,sep="_")) %>%
  dplyr::relocate(ID) %>%
  dplyr::filter(ID %in% new_genotype$ID) %>%
  tidyr::pivot_longer(cols =!c(ID,CHROM,POS, REF, ALT),
                      names_to = "donor",values_to = "genotype") %>%
  dplyr::mutate(genotype = str_split(genotype, ":", simplify = TRUE)[,1]) %>%
  dplyr::mutate(allele1 = as.numeric(str_split_i(genotype, "\\/", i=1)),
                allele2=as.numeric(str_split_i(genotype, "\\/", i=2))) %>%
  dplyr::mutate(gtCount = allele1 + allele2) %>%
  dplyr::select(ID,CHROM,POS,REF,ALT,donor,gtCount) %>%
  dplyr::mutate(type = "wgs")

merged = rbind(wgs_no_imp,new_genotype) %>%
  dplyr::arrange(CHROM,POS,donor,type)

are_gt_identical = merged %>%
  dplyr::group_by(CHROM,POS,REF,ALT,donor) %>%
  dplyr::filter(n()>1) %>%
  dplyr::summarise(are_gt_identical = all(gtCount == first(gtCount))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(CHROM,POS,REF,ALT) %>%
  dplyr::count(are_gt_identical) %>%
  dplyr::mutate(prop_identical = n/sum(n)) %>%
  dplyr::filter(are_gt_identical == TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(.,imputation_info) %>% # adding new imputation info
  dplyr::distinct()

table(are_gt_identical$prop_identical)
13/((40*2)+13) # 14% of all positions and donors tested don't match
