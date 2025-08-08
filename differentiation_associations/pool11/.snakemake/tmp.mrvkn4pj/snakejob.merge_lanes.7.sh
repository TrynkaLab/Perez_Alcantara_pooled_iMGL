#!/bin/sh
# properties = {"type": "single", "rule": "merge_lanes", "local": false, "input": ["../../../data/bams/pool11/pool11_day39_preMAC_L001.bam", "../../../data/bams/pool11/pool11_day39_preMAC_L002.bam"], "output": ["../../../data/bams/pool11/merged_pool11_day39_preMAC.bam"], "wildcards": {"pool": "pool11", "sample": "pool11_day39_preMAC"}, "params": {"group": "-G teamtrynka", "queue": "-q normal", "threads": "-n 16", "memory": "-M100000 -R'span[hosts=1] select[mem>100000] rusage[mem=100000]'", "jobname": "-o ../../../logs/log_merge_lanes.pool11_day39_preMAC.%J.%I", "error": "-e ../../../errors/error_merge_lanes.pool11_day39_preMAC.%J.%I", "samtools_threads": 16}, "log": [], "threads": 1, "resources": {}, "jobid": 7, "cluster": {}}
cd /lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/code/OTAR2065_WGS_iPSC_premac_micro/pool11 && \
/software/teamtrynka/conda/trynka-base/bin/python \
-m snakemake ../../../data/bams/pool11/merged_pool11_day39_preMAC.bam --snakefile /lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/code/OTAR2065_WGS_iPSC_premac_micro/pool11/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/code/OTAR2065_WGS_iPSC_premac_micro/pool11/.snakemake/tmp.mrvkn4pj ../../../data/bams/pool11/pool11_day39_preMAC_L001.bam ../../../data/bams/pool11/pool11_day39_preMAC_L002.bam --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules merge_lanes --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch "/lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/code/OTAR2065_WGS_iPSC_premac_micro/pool11/.snakemake/tmp.mrvkn4pj/7.jobfinished" || (touch "/lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/code/OTAR2065_WGS_iPSC_premac_micro/pool11/.snakemake/tmp.mrvkn4pj/7.jobfailed"; exit 1)

