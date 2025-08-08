#!/bin/bash
# activate correct  module with
# module load HGI/softpack/groups/otar2065/otar2065_9/40
snakemake  --rerun-incomplete --keep-going --rerun-triggers mtime --jobs 5000 --cluster "bsub {params.group} {params.queue} {params.threads} {params.memory} {params.jobname} {params.error}"
# --use-envmodules doesn't seem to work with bsub
