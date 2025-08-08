#!/bin/sh
# properties = {"type": "single", "rule": "subset_and_index_bam", "local": false, "input": ["../../../data/bams/pool15/merged_P15_D36_PreMac_040823.bam"], "output": ["../../../data/bams/pool15/P15_D36_PreMac_040823.bam.bai"], "wildcards": {"pool": "pool15", "sample": "P15_D36_PreMac_040823"}, "params": {"group": "-G teamtrynka", "queue": "-q normal", "threads": "-n 32", "memory": "-M20000 -R'span[hosts=1] select[mem>20000] rusage[mem=20000]'", "jobname": "-o ../../../logs/log_subset_and_index_bam.P15_D36_PreMac_040823.%J.%I", "error": "-e ../../../errors/error_subset_and_index_bam.P15_D36_PreMac_040823.%J.%I"}, "log": [], "threads": 1, "resources": {}, "jobid": 20, "cluster": {}}
cd /lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/code/OTAR2065_WGS_iPSC_premac_micro/pool15 && \
/software/teamtrynka/conda/trynka-base/bin/python \
-m snakemake ../../../data/bams/pool15/P15_D36_PreMac_040823.bam.bai --snakefile /lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/code/OTAR2065_WGS_iPSC_premac_micro/pool15/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/code/OTAR2065_WGS_iPSC_premac_micro/pool15/.snakemake/tmp.xzle60oo ../../../data/bams/pool15/merged_P15_D36_PreMac_040823.bam --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules subset_and_index_bam --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch "/lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/code/OTAR2065_WGS_iPSC_premac_micro/pool15/.snakemake/tmp.xzle60oo/20.jobfinished" || (touch "/lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/code/OTAR2065_WGS_iPSC_premac_micro/pool15/.snakemake/tmp.xzle60oo/20.jobfailed"; exit 1)

