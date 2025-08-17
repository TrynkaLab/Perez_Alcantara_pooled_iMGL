#!/bin/bash
# activate conda with: conda activate /software/hgi/envs/conda/teamtrynka/ma23/snakemake
# activate correct module with: module load HGI/softpack/groups/otar2065/otar2065_no_snakemake-1

################################################################################
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!    Important    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
################################################################################
# Snakemake wasn't run as it didn't work properly when run in an interactive
# session.
################################################################################

## Start Is with: bsub -n2 -q oversubscribed -G teamtrynka -M30000 -R'span[hosts=1] select[mem>30000] rusage[mem=30000]' -Is bash
## from there run:
snakemake  --rerun-incomplete --keep-going --rerun-triggers mtime --jobs 10000 