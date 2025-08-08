#!/bin/bash
# check the lines present in the full vcf (with merged info)
/software/teamtrynka/conda/trynka-base/bin/bcftools view ../../data/full_genotype/chr22/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.20180102.genotypes.chr22.vcf.gz head -n 66 | tail -n 1 > ../../data/full_genotype/vcf_header.txt
# transpose
cat ../../data/full_genotype/vcf_header.txt | datamash -W transpose > ../../data/full_genotype/vcf_header_2.txt
sed -n '10,$p' ../../data/full_genotype/vcf_header_2.txt > ../../data/full_genotype/vcf_header.txt
rm ../../data/full_genotype/vcf_header_2.txt
