#!/bin/bash
/software/teamtrynka/conda/trynka-base/bin/snakemake --rerun-incomplete --jobs 100 --cluster "bsub {params.group} {params.queue} {params.threads} {params.memory} {params.jobname} {params.error}"
