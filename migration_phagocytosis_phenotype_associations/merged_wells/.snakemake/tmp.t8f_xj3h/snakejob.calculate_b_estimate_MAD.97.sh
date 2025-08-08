#!/bin/sh
# properties = {"type": "single", "rule": "calculate_b_estimate_MAD", "local": false, "input": ["../../../data/bam-readcount/P3_146_bam-readcount.txt.gz", "../../../data/genotype/pool_3.merged.genotype.nonIdentical.hg38.vcf.gz", "../../../data/genotype/pool_4.merged.genotype.nonIdentical.hg38.vcf.gz"], "output": ["../../../data/b/P3_146_b_estimate.csv", "../../../data/genotypes/P3_146_genotype_minor_allele_dosage.csv", "../../../data/genotypes/P3_146_Nvariants_retained_reduced_genotype.txt"], "wildcards": {"sample": "P3_146"}, "params": {"group": "-G teamtrynka", "queue": "-q long", "threads": "-n 32", "memory": "-M50000 -R'span[hosts=1] select[mem>50000] rusage[mem=50000]'", "jobname": "-o ../../../logs/log_calculate_b_estimate.P3_146.%J.%I", "error": "-e ../../../errors/error_calculate_b_estimate.P3_146.%J.%I"}, "log": [], "threads": 1, "resources": {}, "jobid": 97, "cluster": {}}
cd /lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_run45219_migration_WGS_june2022/code/OTAR2065_run45219_migration_WGS_june2022/merged_wells && \
/software/teamtrynka/conda/trynka-base/bin/python \
-m snakemake ../../../data/b/P3_146_b_estimate.csv --snakefile /lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_run45219_migration_WGS_june2022/code/OTAR2065_run45219_migration_WGS_june2022/merged_wells/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_run45219_migration_WGS_june2022/code/OTAR2065_run45219_migration_WGS_june2022/merged_wells/.snakemake/tmp.t8f_xj3h ../../../data/bam-readcount/P3_146_bam-readcount.txt.gz ../../../data/genotype/pool_3.merged.genotype.nonIdentical.hg38.vcf.gz ../../../data/genotype/pool_4.merged.genotype.nonIdentical.hg38.vcf.gz --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules calculate_b_estimate_MAD --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch "/lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_run45219_migration_WGS_june2022/code/OTAR2065_run45219_migration_WGS_june2022/merged_wells/.snakemake/tmp.t8f_xj3h/97.jobfinished" || (touch "/lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_run45219_migration_WGS_june2022/code/OTAR2065_run45219_migration_WGS_june2022/merged_wells/.snakemake/tmp.t8f_xj3h/97.jobfailed"; exit 1)

