#!/bin/bash
# first run needs to download and install packages from the python 3 and 2 environments within CELLECT
# conda activate snakemake_cellect (with python v 3.7)
# edit /software/hgi/envs/conda/teamtrynka/ma23/CELLECT/envs/cellectpy3.yml so it has python v.3.7 instead of 3.6
# remove the .snakemake folder and install it again if there are dependency problems, but do it with a job and not the head node
# otherwise there are incompatibilities (weird)
# also make sure there are no other python versions in the path (e.g. from modules) with 
# which python
# run this with bsub -o cellect.o -e cellect.e -M 60000 -R'select[mem>60000] rusage[mem=60000] span[hosts=1]' -n10 -- 'bash run_cellect.sh'
snakemake --use-conda -j -s /software/hgi/envs/conda/teamtrynka/ma23/CELLECT/cellect-magma.snakefile --configfile config.yml
