#!/bin/bash
#BSUB -G teamtrynka
#BSUB -q "oversubscribed"
#BSUB -n 1
#BSUB -M 40000
#BSUB -R "select[mem>40000] rusage[mem=40000]"
#BSUB -o "output%J.log"
#BSUB -e "error%J.log"

## Commands based on Marta's code

module load HGI/softpack/groups/otar2065/otar2065_9/42

## Find rsIDs of all variants in proliferation GWAS 
bcftools annotate --rename-chrs /lustre/scratch123/hgi/projects/otar2065/OTAR2065_sc_eQTL/data/genotype/dbsnp/rename_chrs.txt  \
/lustre/scratch123/hgi/projects/otar2065/OTAR2065_sc_eQTL/data/genotype/dbsnp/dbsnp_156_GRCh38.vcf.gz | bcftools view -T ../../output_data/08_GWAS_proliferation/02_Check_association_results/variants_for_GWAS.txt -Oz -o ../../input_data/08_GWAS_proliferation/02_Check_association_results/dbsnp.156.variants_for_GWAS.subset.vcf.gz

bcftools index ../../input_data/08_GWAS_proliferation/02_Check_association_results/dbsnp.156.variants_for_GWAS.subset.vcf.gz

## Split multiallelic variants 
bcftools norm -m-any ../../input_data/08_GWAS_proliferation/02_Check_association_results/dbsnp.156.variants_for_GWAS.subset.vcf.gz | bcftools query -f '%CHROM %POS %ID %REF %ALT\n' > ../../input_data/08_GWAS_proliferation/02_Check_association_results/dbsnp.156.variants_for_GWAS.subset.tsv
