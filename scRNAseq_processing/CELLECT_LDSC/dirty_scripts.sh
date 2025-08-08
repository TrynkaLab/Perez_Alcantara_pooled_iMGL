conda activate munge_ldsc

#/software/hgi/envs/conda/teamtrynka/ma23/CELLECT/ldsc/mtag_munge.py \
# AD Bellenguez
/software/team152/oa3/CELLECT/ldsc/mtag_munge.py \
--sumstats /lustre/scratch123/hgi/teams/trynka/resources/summary_statistics/public/GCST90027158/harmonised/35379992-GCST90027158-MONDO_0004975-Build38.f.tsv.gz \
--a1 other_allele \
--a2 effect_allele \
--snp variant_id \
--N-cas-col n_cas \
--N-con-col n_con \
--merge-alleles /software/team152/oa3/CELLECT/data/ldsc/w_hm3.snplist \
--keep-pval \
--p p_value \
--out ../../../data/results/CELLECT/GCST90027158
# AD Jeremy
/software/team152/oa3/CELLECT/ldsc/mtag_munge.py \
--sumstats /lustre/scratch123/hgi/teams/trynka/resources/summary_statistics/public/GCST90012877/harmonised/33589840-GCST90012877-EFO_0000249-Build37.f.tsv.gz \
--a1 other_allele \
--a2 effect_allele \
--snp variant_id \
--N-cas 53042 \
--N-con 355900 \
--merge-alleles /software/team152/oa3/CELLECT/data/ldsc/w_hm3.snplist \
--keep-pval \
--p p_value \
--out ../../../data/results/CELLECT/GCST90012877
# PD Nalls
/software/team152/oa3/CELLECT/ldsc/mtag_munge.py \
--sumstats /lustre/scratch123/hgi/projects/otar2065/resources/harmonised_GWAS/ieu-b-7_PD_Nalls_2019/harmonised/ieu-b-7_manually_fixed_for_harmonisation_b38.for_cellect.h.tsv.gz \
--a1 other_allele \
--a2 effect_allele \
--snp hm_rsid \
--N-cas 33674 \
--N-con 449056 \
--merge-alleles /software/team152/oa3/CELLECT/data/ldsc/w_hm3.snplist \
--keep-pval \
--p p_value \
--out ../../../data/results/CELLECT/ieu-b-7_PD_Nalls_2019
# ALS van Rheenen
/software/team152/oa3/CELLECT/ldsc/mtag_munge.py \
--sumstats /lustre/scratch123/hgi/teams/trynka/resources/summary_statistics/public/GCST90027163/harmonised/34873335-GCST90027163-MONDO_0004976.for_cellect.h.tsv.gz \
--a1 other_allele \
--a2 effect_allele \
--snp rsid \
--N-cas 27205 \
--N-con 110881 \
--merge-alleles /software/team152/oa3/CELLECT/data/ldsc/w_hm3.snplist \
--keep-pval \
--p p_value \
--out ../../../data/results/CELLECT/GCST90027163

# Lewy body dementia Chia
/software/team152/oa3/CELLECT/ldsc/mtag_munge.py \
--sumstats /lustre/scratch123/hgi/teams/trynka/resources/summary_statistics/public/GCST90001390/harmonised/33589841-GCST90001390-EFO_0006792.for_cellect.h.tsv.gz \
--a1 other_allele \
--a2 effect_allele \
--snp rsid \
--N-cas 2591 \
--N-con 4027 \
--merge-alleles /software/team152/oa3/CELLECT/data/ldsc/w_hm3.snplist \
--keep-pval \
--p p_value \
--out ../../../data/results/CELLECT/GCST90001390

# Height Yengo
/software/team152/oa3/CELLECT/ldsc/mtag_munge.py \
--sumstats /lustre/scratch123/hgi/teams/trynka/resources/summary_statistics/public/GCST006901/height_2018_30124842_b37.tsv.gz \
--a1 other_allele \
--a2 effect_allele \
--snp variant_id \
--N 693529 \
--merge-alleles /software/team152/oa3/CELLECT/data/ldsc/w_hm3.snplist \
--keep-pval \
--p p_value \
--out ../../../data/results/CELLECT/GCST006901
### cellex specificity scores from pseudobulks per donor, pool and treatment
# ../../../data/results/CELLECT/cellex_treatment.csv.gz
# create the gene annotation file (GRCh38)
# NOT NEEDED, WORKS WITH grch37 ANNOTATIONSX
#awk -F, '{print $9 "\t" $1 "\t" $4 "\t" $5 "\t" $7 "\t" $10}' ../../../data/results/CELLECT/gene_annotations/Homo_sapiens.GRCh38.111.genes.csv > ../../../data/results/CELLECT/gene_annotations/Homo_sapiens.GRCh38.111.genes.tsv
# remove any header columns

# then run the snakefile with the config.yaml properly configured
