# compare vcf files from old hipsci genotype chip (GRCh37 lifter over) and new WGS
# for kolf_2 and aowh_2

.libPaths(c('/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/',"/software/teamtrynka/ma23/R4.1/libs",
            "/software/R-4.1.0/lib/R/library"))

library(vcfR)
library(tidyverse)
library(patchwork)
outdir="../../data/kolf_aowh/comparison/"
dir.create(outdir)

# compare positions present in old, not in new
# and viceversa
# and of those present in both, how many calls are different


old = vcfR::read.vcfR("../../data/kolf_aowh/aowh_2.genotype.hg38.vcf.gz")
new = vcfR::read.vcfR("../../data/ungenotyped_donors/deepvariant_nextflow/aowh_2.deepvariant.vcf.gz")

varnumber = c(nrow(old),nrow(new))
names(varnumber) = c("Nvar_old","Nvar_new")

tidy_old= as_tibble(old@fix)


tidy_old %>%
  dplyr::tally(is.na(QUAL)) # all NA in QUAL for old genotype
tidy_old %>%
  dplyr::tally(FILTER =="PASS") # all deemed to PASS filters

# check out GT info
old_gt = vcfR::extract_gt_tidy(old)

old_gt%>%
  group_by(gt_GT) %>%
  summarize(sum_gt = n())

p1 =   old_gt %>%
  ggplot(aes(x = as.factor(gt_GT))) + 
  geom_bar() +
  scale_y_log10(limits=c(1,2300000),breaks = c(10,1000,10000,1e6,2.3e6))+
  ggtitle("xx")+
  theme_bw()

p2 = old_gt %>%
  ggplot(aes(x = as.numeric(gt_GQ))) + 
  geom_density() +
  geom_vline(xintercept = 0.02, col = "red")+
  geom_vline(xintercept = 20, col = "red")+
  annotate("text", x = Inf, y = Inf, label = paste0("Nvar:",nrow(different)), hjust = 1.2, vjust = 1.5) +
  ggtitle("xx")+
  theme_bw()
######### intersecting with new WGS ######
rm(old)
gc()

tidy_new = as_tibble(new@fix) %>%
  dplyr::mutate(CHROM=paste0("chr",CHROM))

# present in both
intersect = inner_join(tidy_old, tidy_new, by = c("CHROM", "POS"))
intersect_varnumber = nrow(intersect)

# How many calls are different
table(intersect$REF.x == intersect$REF.y & intersect$ALT.x == intersect$ALT.y)
# check wgs quality of those that are different
different = intersect[!(intersect$REF.x == intersect$REF.y & intersect$ALT.x == intersect$ALT.y),]

bcf_intersect = vcfR::read.vcfR("../../data/kolf_aowh/intersect/aowh_2/0003.vcf.gz")
tidy_bcf = as_tibble(bcf_intersect@fix)

nrow(intersect) - nrow(different) # should be same number of variants as intersected with bcftools
# 3,730,919
nrow(tidy_bcf)
# 3,730,716 # but it's not. bcftools doing additional filters in intersection?
# might have excluded some low quality calls

# checking quality
# QUAL is phred score of the probability of there being ALT allele (so 0 just means REF). Probabilities above 0.05 and under 20 might be problematic
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531872-Phred-scaled-quality-scores
# a 0 QUAL means 100% chance REF; 0.02 means 99.5% chance REF; a 20 QUAL means 1% chance REF (so 99% chance ALT)
# 
10^-(0.02/10)
10^-(20/10)
# GQ tells you how confident we are that the genotype we assigned to a particular sample is correct - important
# PL: "Normalized" Phred-scaled likelihoods of the possible genotypes. 
# The PL values are "normalized" so that the PL of the most likely genotype (assigned in the GT field) is 0 in the Phred scale.
# so a 1/1 witll be something,something,0
# a 0/1 will be something,0,something
# check here: https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format

p1 = tidy_bcf %>%
  ggplot(aes(x = as.numeric(QUAL))) + 
  geom_density() +
  geom_vline(xintercept = 0.02, col = "red")+
  geom_vline(xintercept = 20, col = "red")+
  annotate("text", x = Inf, y = Inf, label = paste0("Nvar:",nrow(tidy_bcf)), hjust = 1.2, vjust = 1.5) +
  ggtitle("bcftools intersect old genotype vs WGS calls")+
  theme_bw()

p1

p2 = tidy_new %>%
  ggplot(aes(x = as.numeric(QUAL))) + 
  geom_density() +
  geom_vline(xintercept = 0.02, col = "red")+
  geom_vline(xintercept = 20, col = "red")+
  annotate("text", x = Inf, y = Inf, label = paste0("Nvar:",nrow(tidy_new)),  hjust = 1.2, vjust = 1.5) +
  ggtitle("full file WGS calls")+
  theme_bw()

p2

p3 = intersect %>%
  ggplot(aes(x = as.numeric(QUAL.y))) + 
  geom_density() +
  geom_vline(xintercept = 0.02, col = "red")+
  geom_vline(xintercept = 20, col = "red")+
  annotate("text", x = Inf, y = Inf, label = paste0("Nvar:",nrow(intersect)), hjust = 1.2, vjust = 1.5) +
  ggtitle("tidy intersect old genotype vs WGS calls")+
  theme_bw()

p3

p4 = different %>%
  ggplot(aes(x = as.numeric(QUAL.y))) + 
  geom_density() +
  geom_vline(xintercept = 0.02, col = "red")+
  geom_vline(xintercept = 20, col = "red")+
  annotate("text", x = Inf, y = Inf, label = paste0("Nvar:",nrow(different)), hjust = 1.2, vjust = 1.5) +
  ggtitle("subset of the intersect with different REF or ALT")+
  theme_bw()

p4

png(paste0(outdir,"QUAL_dist_intersect_WGS_old_geno.png"),res = 400,units = "in",
    width = 12,height = 12)
plot((p1+p2) / (p3+p4))
dev.off()
# for intersecting positions that don't have different REF/ALT ("different")
# how many genotype calls (0/0,0/1,1/1) are different between old and new?


# are the intersect most 0/1 or 1/1, as  the distribution of QUAL suggests?
bcf_isect_gt = vcfR::extract_gt_tidy(bcf_intersect)
# GQ - individual genotype quality.  Phred-scaled confidence that the genotype assignment (GT) is correct
#  the GQ is the difference between the PL of the second most likely genotype, and the PL of the most likely genotype.
bcf_isect_gt %>%
  group_by(gt_GT) %>%
  summarize(sum_gt = n())

p1 =   bcf_isect_gt %>%
  ggplot(aes(x = as.factor(gt_GT))) + 
  geom_bar() +
  scale_y_log10(limits=c(1,2300000),breaks = c(10,1000,10000,1e6,2.3e6))+
  ggtitle("xx")+
  theme_bw()


# why are so many of the intersect ALT variant calls?
# checking same info in full WGS variant calls from deepvariant
new_gt =  vcfR::extract_gt_tidy(new)
new_gt %>%
  group_by(gt_GT) %>%
  summarize(sum_gt = n())
# many multiallelic variants
# (0/2...N, and 1/2...N)
# also many NAs - probably low quality

p2 =   new_gt %>%
  ggplot(aes(x = as.factor(gt_GT))) + 
  geom_bar() +
  scale_y_log10(limits=c(1,2300000),breaks = c(10,1000,10000,1e6,2.3e6))+
  ggtitle("xx")+
  theme_bw()

# why are there NAs?

p3 = bcf_isect_gt %>%
  ggplot(aes(x = as.numeric(gt_GQ))) + 
  geom_density() +
  geom_vline(xintercept = 0.02, col = "red")+
  geom_vline(xintercept = 20, col = "red")+
  annotate("text", x = Inf, y = Inf, label = paste0("Nvar:",nrow(different)), hjust = 1.2, vjust = 1.5) +
  ggtitle("xx")+
  theme_bw()
# most are high quality (>20, >99% chance correct call compared to other allele)

png(paste0(outdir,"QUAL_dist_intersect_WGS_old_geno.png"),res = 400,units = "in",
    width = 12,height = 12)
plot((p1+p2) / (p3+p4))
dev.off()

# getting variants that are different from wgs calls


p1 =   bcf_isect_gt %>%
  ggplot(aes(x = as.factor(gt_GT))) + 
  geom_bar() +
  scale_y_log10(limits=c(1,2300000),breaks = c(10,1000,10000,1e6,2.3e6))+
  ggtitle("subset of the intersect with different REF or ALT - genotype calls")+
  theme_bw()

bcf_isect_gt %>%
  group_by(gt_GT) %>%
  summarize(sum_gt = n())
# why are there NAs?
# obvs the subset of the intersect with different REF or ALT will have way more 0/1 and 1/1
# maybe the genotypes at that position were imputed in the old, and there are multiple possible SNPs at that place
# probably something to do with multiallelic variants because there are a non-trivial number, as shown above

p2 = bcf_isect_gt %>%
  ggplot(aes(x = as.numeric(gt_GQ))) + 
  geom_density() +
  geom_vline(xintercept = 0.02, col = "red")+
  geom_vline(xintercept = 20, col = "red")+
  annotate("text", x = Inf, y = Inf, label = paste0("Nvar:",nrow(different)), hjust = 1.2, vjust = 1.5) +
  ggtitle("subset of the intersect with different REF or ALT - GQ")+
  theme_bw()
# most are high quality (>20, >99% chance correct call compared to other allele)


png(paste0(outdir,"QUAL_dist_intersect_WGS_old_geno.png"),res = 400,units = "in",
    width = 12,height = 12)
plot((p1+p2) / (p3+p4))
dev.off()

hist(as.numeric(different$QUAL.y))
table(different$FILTER.x)
table(different$FILTER.y)
table(different$FILTER.y == "PASS" & different$QUAL.y<20)

dplyr::filter(QUAL > 0.02 & QUAL < 20)

# WGS: how many of the positions that are different match better the calls from 
# bam-readcount than the old genotype?
# and how many are ./.

# check what happens in the calculation of MAF (bam-readcount, b steps) in the deconvolution if the 
# don't match the ALT allele or the REF, and what happens when it matches one but not the other