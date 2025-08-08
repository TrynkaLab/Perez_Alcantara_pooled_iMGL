#!/bin/sh
# properties = {"type": "single", "rule": "get_region_list", "local": false, "input": ["../../../data/genotypes/pool13/pool13.merged.genotype.nonIdentical.hg38.vcf.gz"], "output": ["../../../data/genotypes/pool13/pool13_genotype.txt"], "wildcards": {"pool": "pool13"}, "params": {"group": "-G teamtrynka", "queue": "-q normal", "threads": "-n 2", "memory": "-M40000 -R'span[hosts=1] select[mem>40000] rusage[mem=40000]'", "jobname": "-o ../../../logs/log_get_region_list.%J.%I", "error": "-e ../../../errors/error_get_region_list.%J.%I"}, "log": [], "threads": 1, "resources": {}, "jobid": 13, "cluster": {}}
cd /lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/code/OTAR2065_WGS_iPSC_premac_micro/pool13 && \
/software/teamtrynka/conda/trynka-base/bin/python \
-m snakemake ../../../data/genotypes/pool13/pool13_genotype.txt --snakefile /lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/code/OTAR2065_WGS_iPSC_premac_micro/pool13/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/code/OTAR2065_WGS_iPSC_premac_micro/pool13/.snakemake/tmp.1pl133ri ../../../data/genotypes/pool13/pool13.merged.genotype.nonIdentical.hg38.vcf.gz --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules get_region_list --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch "/lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/code/OTAR2065_WGS_iPSC_premac_micro/pool13/.snakemake/tmp.1pl133ri/13.jobfinished" || (touch "/lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/code/OTAR2065_WGS_iPSC_premac_micro/pool13/.snakemake/tmp.1pl133ri/13.jobfailed"; exit 1)

