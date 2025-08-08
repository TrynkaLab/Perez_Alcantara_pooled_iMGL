#!/bin/bash
snakemake --rerun-incomplete --keep-going --rerun-triggers mtime --jobs 1000 --cluster "bsub {params.group} {params.queue} {params.threads} {params.memory} {params.jobname} {params.error}"
