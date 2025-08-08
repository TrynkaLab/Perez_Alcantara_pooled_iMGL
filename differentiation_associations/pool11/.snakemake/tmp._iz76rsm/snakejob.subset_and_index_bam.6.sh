#!/bin/sh
# properties = {"type": "single", "rule": "subset_and_index_bam", "local": false, "input": ["../../../data/bams/pool11/merged_pool11_day49_preMAC.bam"], "output": ["../../../data/bams/pool11/pool11_day49_preMAC.bam.bai"], "wildcards": {"pool": "pool11", "sample": "pool11_day49_preMAC"}, "params": {"group": "-G teamtrynka", "queue": "-q normal", "threads": "-n 32", "memory": "-M20000 -R'span[hosts=1] select[mem>20000] rusage[mem=20000]'", "jobname": "-o ../../../logs/log_subset_and_index_bam.pool11_day49_preMAC.%J.%I", "error": "-e ../../../errors/error_subset_and_index_bam.pool11_day49_preMAC.%J.%I"}, "log": [], "threads": 1, "resources": {}, "jobid": 6, "cluster": {}}
cd /lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/code/OTAR2065_WGS_iPSC_premac_micro/pool11 && \
/software/teamtrynka/conda/trynka-base/bin/python \
-m snakemake ../../../data/bams/pool11/pool11_day49_preMAC.bam.bai --snakefile /lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/code/OTAR2065_WGS_iPSC_premac_micro/pool11/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/code/OTAR2065_WGS_iPSC_premac_micro/pool11/.snakemake/tmp._iz76rsm ../../../data/bams/pool11/merged_pool11_day49_preMAC.bam --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules subset_and_index_bam --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch "/lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/code/OTAR2065_WGS_iPSC_premac_micro/pool11/.snakemake/tmp._iz76rsm/6.jobfinished" || (touch "/lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/code/OTAR2065_WGS_iPSC_premac_micro/pool11/.snakemake/tmp._iz76rsm/6.jobfailed"; exit 1)

