#!/bin/bash
snakemake --rerun-incomplete --rerun-triggers mtime --jobs 2000 --keep-going --cluster "bsub {params.group} {params.queue} {params.threads} {params.memory} {params.jobname} {params.error}"
