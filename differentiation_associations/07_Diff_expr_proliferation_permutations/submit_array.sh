#!/bin/bash
# Submit an array of 100 jobs of DGE, one per permuted proliferation phenotypes 

module load HGI/softpack/groups/otar2065/otar2065_9/40 
bsub  -J "array_DGE_jobs[1-100]"  -q week-chkpt -R "select[mem>40000] rusage[mem=40000]" -M 40000 -G teamtrynka -o "./DGE_jobs_out_err/out.%J.%I" -e  "./DGE_jobs_out_err/err.%J.%I" sh ./run_DGE_Rscript_per_job.sh

