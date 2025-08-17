
#!/bin/bash
#BSUB -G teamtrynka
#BSUB -q "oversubscribed"
#BSUB -n 1
#BSUB -M 40000
#BSUB -R "select[mem>40000] rusage[mem=40000]"
#BSUB -o "output%J.log"
#BSUB -e "error%J.log"

## Commands based on Marta's code

module load HGI/softpack/groups/otar2065/otar2065_9/34

bcftools annotate --rename-chrs /lustre/scratch123/hgi/projects/otar2065/OTAR2065_sc_eQTL/data/genotype/dbsnp/rename_chrs.txt  \
/lustre/scratch123/hgi/projects/otar2065/OTAR2065_sc_eQTL/data/genotype/dbsnp/dbsnp_156_GRCh38.vcf.gz | bcftools view -T ../../output_data/03_Mash_analysis/00_prepare_input_data/variants_lead_signif_macro.txt -Oz -o ../../input_data/03_Mash_analysis/00_prepare_input_data/dbsnp.156.variants_lead_signif_macromap_Panousis_2023.subset.vcf.gz
bcftools index ../../input_data/03_Mash_analysis/00_prepare_input_data/dbsnp.156.variants_lead_signif_macromap_Panousis_2023.subset.vcf.gz
# split multiallelic variants so each alt gets their own row, and subset to first 5 columns
bcftools norm -m-any ../../input_data/03_Mash_analysis/00_prepare_input_data/dbsnp.156.variants_lead_signif_macromap_Panousis_2023.subset.vcf.gz | bcftools query -f '%CHROM %POS %ID %REF %ALT\n' > ../../input_data/03_Mash_analysis/00_prepare_input_data/dbsnp.156.variants_lead_signif_macromap_Panousis_2023.subset.tsv