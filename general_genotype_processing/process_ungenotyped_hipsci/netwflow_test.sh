#!/bin/bash
# needs nextflow v23

set -e

export HTTP_PROXY='http://wwwcache.sanger.ac.uk:3128'
export HTTPS_PROXY='http://wwwcache.sanger.ac.uk:3128'
export NXF_ANSI_LOG=false
export NXF_OPTS="-Xms14G -Xmx14G -Dnxf.pool.maxThreads=2000"

# General Nextflow variables #

# Nextflow version, change depending on installation you use
#export NXF_VER=23.04.2.5870

# singularity cache
export NXF_SINGULARITY_CACHEDIR=/lustre/scratch123/hgi/projects/otar2065/hipsci_genotype_processing/data/ungenotyped_donors/resources/nextflow/cache/singularity


OUTDIR=../../data/ungenotyped_donors/test
SAMPLESHEET=../../data/ungenotyped_donors/samplesheet.csv
mkdir -p ${OUTDIR}
/software/teamtrynka/nextflow_v23/nextflow run /software/teamtrynka/nextflow_pipelines/nf-core/sarek/workflow/main.nf \
--input ${SAMPLESHEET} \
--outdir ${OUTDIR} \
--genome null \
--igenomes_ignore \
--fasta /lustre/scratch123/hgi/teams/trynka/resources/GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz \
--save_reference \
--dbsnp /lustre/scratch123/hgi/teams/trynka/resources/dbsnp/151/b38/All_20180418.vcf.gz \
-profile singularity \
-c ./nextflow.config
-resume

# to use custom fasta and gtf I must set --genome null and --igenomes_ignore
# --save-reference will save the indices, interval files etc in the results
# directory for you to move and store in a more central location for re-use with future pipeline runs
