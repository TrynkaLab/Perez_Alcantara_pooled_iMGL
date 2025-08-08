#!/bin/bash
snakemake --rerun-incomplete --printshellcmds --jobs 100 --cluster "bsub {params.group} {params.queue} {params.threads} {params.memory} {params.jobname} {params.error}"
