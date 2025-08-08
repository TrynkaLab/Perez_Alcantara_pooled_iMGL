#!/bin/sh
# properties = {"type": "single", "rule": "cram_to_bam", "local": false, "input": ["../../../data/crams/pool11/pool11_day35_preMAC_L001.cram", "../../../data/crams/pool11/pool11_day35_preMAC_L002.cram"], "output": ["../../../data/bams/pool11/pool11_day35_preMAC_L001.bam", "../../../data/bams/pool11/pool11_day35_preMAC_L002.bam"], "wildcards": {"pool": "pool11", "sample": "pool11_day35_preMAC"}, "params": {"group": "-G teamtrynka", "queue": "-q normal", "threads": "-n 32", "memory": "-M100000 -R'span[hosts=1] select[mem>100000] rusage[mem=100000]'", "jobname": "-o ../../../logs/log_cram_to_bam.pool11_day35_preMAC.%J.%I", "error": "-e ../../../errors/error_cram_to_bam.pool11_day35_preMAC.%J.%I", "reference": "/lustre/scratch125/core/sciops_repository/references/Homo_sapiens/GRCh38_15_plus_hs38d1/all/fasta/Homo_sapiens.GRCh38_15_plus_hs38d1.fa"}, "log": [], "threads": 1, "resources": {}, "jobid": 12, "cluster": {}}
cd /lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/code/OTAR2065_WGS_iPSC_premac_micro/pool11 && \
/software/teamtrynka/conda/trynka-base/bin/python \
-m snakemake ../../../data/bams/pool11/pool11_day35_preMAC_L001.bam --snakefile /lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/code/OTAR2065_WGS_iPSC_premac_micro/pool11/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/code/OTAR2065_WGS_iPSC_premac_micro/pool11/.snakemake/tmp.r3yulrw1 ../../../data/crams/pool11/pool11_day35_preMAC_L001.cram ../../../data/crams/pool11/pool11_day35_preMAC_L002.cram --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules cram_to_bam --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch "/lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/code/OTAR2065_WGS_iPSC_premac_micro/pool11/.snakemake/tmp.r3yulrw1/12.jobfinished" || (touch "/lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/code/OTAR2065_WGS_iPSC_premac_micro/pool11/.snakemake/tmp.r3yulrw1/12.jobfailed"; exit 1)

