#!/bin/bash
/software/teamtrynka/conda/otar2065/bin/snakemake --jobs 100 --cluster "bsub {params.group} {params.queue} {params.threads} {params.memory} {params.jobname} {params.error}"
