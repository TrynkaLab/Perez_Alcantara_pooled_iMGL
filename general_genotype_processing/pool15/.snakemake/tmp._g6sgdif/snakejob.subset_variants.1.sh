#!/bin/sh
# properties = {"type": "single", "rule": "subset_variants", "local": false, "input": ["../../data/pool_15/pool_15.merged.genotype.hg38.vcf.gz"], "output": ["../../data/pool_15/pool_15.merged.genotype.nonIdentical.hg38.vcf.gz"], "wildcards": {"pool": "pool_15"}, "params": {"group": "-G teamtrynka", "queue": "-q normal", "threads": "-n 1", "memory": "-M160000 -R'span[hosts=1] select[mem>160000] rusage[mem=160000]'", "jobname": "-o ../../logs/log_pool_15_subset_variants_WGS.%J.%I", "error": "-e ../../errors/error_pool_15_subset_variants_WGS.%J.%I"}, "log": [], "threads": 1, "resources": {}, "jobid": 1, "cluster": {}}
cd /lustre/scratch123/hgi/mdt1/projects/otar2065/hipsci_genotype_processing/code/pool15 && \
/software/teamtrynka/conda/trynka-base/bin/python \
-m snakemake ../../data/pool_15/pool_15.merged.genotype.nonIdentical.hg38.vcf.gz --snakefile /lustre/scratch123/hgi/mdt1/projects/otar2065/hipsci_genotype_processing/code/pool15/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /lustre/scratch123/hgi/mdt1/projects/otar2065/hipsci_genotype_processing/code/pool15/.snakemake/tmp._g6sgdif ../../data/pool_15/pool_15.merged.genotype.hg38.vcf.gz --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules subset_variants --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch "/lustre/scratch123/hgi/mdt1/projects/otar2065/hipsci_genotype_processing/code/pool15/.snakemake/tmp._g6sgdif/1.jobfinished" || (touch "/lustre/scratch123/hgi/mdt1/projects/otar2065/hipsci_genotype_processing/code/pool15/.snakemake/tmp._g6sgdif/1.jobfailed"; exit 1)

