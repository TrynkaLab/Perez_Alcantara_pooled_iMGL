#!/bin/bash
# activate conda with 
#  conda activate /software/hgi/envs/conda/teamtrynka/ma23/snakemake
# activate correct module with
# module load HGI/softpack/groups/otar2065/otar2065_9/42
snakemake  --rerun-incomplete --keep-going --rerun-triggers mtime --jobs 1000  
# --use-envmodules doesn't seem to work with bsub
