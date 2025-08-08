#!/usr/bin/env bash

set -o errexit
set -o pipefail
set -o nounset

USAGE="\n\nUSAGE:\n\t $0 <input file(vcf format)>\n\n\t      The absolute path should be resolved and the output\n\t      will use the original file name but with vep.vcf in\n\t      place of the vcf in the input filename\n\t      filename must contain \"vcf\"\n\n"

if [[ $HOSTNAME =~ farm5-head[0-9] ]]
then
        echo -e "\n\n\tYour trying to run this on a head node, that is not possible,\n\ttry running via bsub\n\teg. bsub -q normal -R \"select[mem>8000] rusage[mem=8000]\" -M8000 -o vep.out -e vep.err $0\n\n"
        exit 1
fi


if [ -z $@ ]
then
	echo -e $USAGE
        exit 1
else
	INPUT_VCF=$(realpath $1)
fi
VCF_DIR=$(dirname $INPUT_VCF)
VCF_FILE=$(basename $INPUT_VCF)
OUTPUT_VCF=$(echo $VCF_FILE | sed 's/vcf/vep.vcf/')

if  [[ ! $VCF_FILE =~ .*vcf.* ]]
then
	echo -e $USAGE
	exit 1
elif [ ! -e $INPUT_VCF ]
then
	echo -e "\n\n\t   File $1 does not exist\n\tCheck $1 is not a SymLink\n\n\n"
        exit 1
fi

#echo $INPUT_VCF
#echo $VCF_DIR
#echo $VCF_FILE
#echo $OUTPUT_VCF


#module load cellgen/singularity
singularity exec \
--bind $VCF_DIR:/opt/vcf \
--bind /lustre/scratch125/humgen/resources/ensembl/vep/GRCh37/vep_data:/opt/vep/.vep \
--bind /lustre/scratch125/humgen/resources/ensembl/vep/GRCh37/vep_data/Plugins:/opt/vep/.vep/Plugins \
--bind /lustre/scratch125/humgen/resources/gnomAD/release-2.1.1/exomes \
--bind /lustre/scratch125/humgen/resources/cadd_scores/20201027-GRCh37_v1.6 \
--bind /lustre/scratch125/humgen/resources/SpliceAI_data_files \
--bind /lustre/scratch123/hgi/teams/trynka/resources/CADD/v1.7_GRCh37 \
/lustre/scratch125/humgen/resources/ensembl/vep/singularity_containers/vep_111.0.sif \
vep \
--cache \
--dir_cache /opt/vep/.vep/ \
--assembly GRCh37 \
--fasta /opt/vep/.vep/homo_sapiens/102_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
--offline \
--format vcf \
--dir_plugins /opt/vep/.vep/Plugins \
-i /opt/vcf/$VCF_FILE \
--plugin SpliceRegion,Extended \
--plugin GeneSplicer,/opt/vep/.vep/Plugins/GeneSplicer/bin/linux/genesplicer,/opt/vep/.vep/Plugins/GeneSplicer/human \
--plugin UTRannotator,/opt/vep/.vep/Plugins/uORF_5UTR_GRCh37_PUBLIC.txt \
--plugin CADD,/lustre/scratch123/hgi/teams/trynka/resources/CADD/v1.7_GRCh37/whole_genome_SNVs.tsv.gz,/lustre/scratch123/hgi/teams/trynka/resources/CADD/v1.7_GRCh37/gnomad.genomes-exomes.r4.0.indel.tsv.gz \
--fork 4 \
--everything \
--plugin LoF,loftee_path:/opt/vep/.vep/Plugins,human_ancestor_fa:/opt/vep/.vep/Plugins/grch37_human_ancestor.fa.gz,conservation_file:/opt/vep/.vep/Plugins/phylocsf_gerp.sql  \
--plugin REVEL,/opt/vep/.vep/Plugins/grch37_tabbed_revel.tsv.gz \
--plugin SpliceAI,snv=/lustre/scratch125/humgen/resources/SpliceAI_data_files/spliceai_scores.raw.snv.hg19.vcf.gz,indel=/lustre/scratch125/humgen/resources/SpliceAI_data_files/spliceai_scores.raw.indel.hg19.vcf.gz \
--vcf \
-o /opt/vcf/$OUTPUT_VCF \
--compress_output bgzip \
--allele_number \
--verbose 
