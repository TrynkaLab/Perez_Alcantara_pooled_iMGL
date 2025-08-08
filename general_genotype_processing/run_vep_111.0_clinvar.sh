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


module load cellgen/singularity
singularity exec \
--bind $VCF_DIR:/opt/vcf \
--bind /lustre/scratch125/humgen/resources/ensembl/vep/GRCh38/vep_data:/opt/vep/.vep \
--bind /lustre/scratch125/humgen/resources/ensembl/vep/GRCh38/Plugins:/opt/vep/.vep/Plugins \
--bind /lustre/scratch125/humgen/resources/gnomAD/release-2.1.1/exomes \
--bind /lustre/scratch125/humgen/resources/cadd_scores/20201027-GRCh38_v1.6 \
--bind /lustre/scratch125/humgen/resources/SpliceAI_data_files \
/lustre/scratch125/humgen/resources/ensembl/vep/singularity_containers/vep_111.0.sif \
vep \
--cache \
--dir_cache /opt/vep/.vep/ \
--assembly GRCh38 \
--fasta /opt/vep/.vep/homo_sapiens/111_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
--offline \
--format vcf \
--dir_plugins /opt/vep/.vep/Plugins \
-i /opt/vcf/$VCF_FILE \
--plugin SpliceRegion,Extended \
--plugin GeneSplicer,/opt/vep/.vep/Plugins/GeneSplicer/bin/linux/genesplicer,/opt/vep/.vep/Plugins/GeneSplicer/human \
--plugin UTRannotator,/opt/vep/.vep/Plugins/uORF_starts_ends_GRCh38_PUBLIC.txt \
--plugin CADD,/lustre/scratch125/humgen/resources/cadd_scores/20201027-GRCh38_v1.6/whole_genome_SNVs.tsv.gz,/lustre/scratch125/humgen/resources/cadd_scores/20201027-GRCh38_v1.6/gnomad.genomes.r3.0.indel.tsv.gz \
--plugin dbNSFP,/opt/vep/.vep/Plugins/dbNSFP4.5a_grch38.gz,Ensembl_transcriptid,Uniprot_acc,VEP_canonical,LRT_pred,SIFT_pred,MutationTaster_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,AlphaMissense_score,AlphaMissense_rankscore,AlphaMissense_pred,PrimateAI_score,PrimateAI_rankscore,PrimateAI_pred \
--plugin DisGeNET,file=/opt/vep/.vep/Plugins/all_variant_disease_pmid_associations_final.tsv.gz \
--custom /opt/vep/.vep/Plugins/clinvar_20240611.vcf.gz,ClinVar,vcf,exact,0,CLNDN,CLNDNINCL,CLNDISDB,CLNDISDBINCL,AF_ESP,AF_EXAC,AF_TGP,ALLELEID,CLNDN,CLNDNINCL,CLNDISDB,CLNDISDBINCL,CLNHGVS,CLNREVSTAT,CLNSIG,CLNSIGCONF,CLNSIGINCL,CLNVC,CLNVCSO,CLNVI,DBVARID,GENEINFO,MC,ORIGIN,RS,SSR \
--fork 4 \
--everything \
--plugin Phenotypes,file=/opt/vep/.vep/Plugins/phenotypes.gff.gz,include_types=Variation \
--plugin Conservation,/opt/vep/.vep/Plugins/90_mammals.gerp_conservation_scores.homo_sapiens.GRCh38.bw \
--plugin LoF,loftee_path:/opt/vep/.vep/Plugins,human_ancestor_fa:/opt/vep/.vep/Plugins/GRCh38_human_ancestor.fa.gz,conservation_file:/opt/vep/.vep/Plugins/loftee.sql,gerp_bigwig:/opt/vep/.vep/Plugins/90_mammals.gerp_conservation_scores.homo_sapiens.GRCh38.bw \
--plugin REVEL,/opt/vep/.vep/Plugins/grch38_tabbed_revel.tsv.gz \
--plugin SpliceAI,snv=/lustre/scratch125/humgen/resources/SpliceAI_data_files/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/lustre/scratch125/humgen/resources/SpliceAI_data_files/spliceai_scores.raw.indel.hg38.vcf.gz \
--plugin EVE,file=/opt/vep/.vep/Plugins/EVE_data/eve_plugin/eve_merged.vcf.gz \
--plugin Enformer,file=/opt/vep/.vep/Plugins/enformer_grch38.vcf.gz \
--plugin MaveDB,file=/opt/vep/.vep/Plugins/MaveDB_variants.tsv.gz,cols=all \
--plugin RiboseqORFs,file=/opt/vep/.vep/Plugins/Ribo-seq_ORFs.bed.gz \
--plugin AncestralAllele,/opt/vep/.vep/Plugins/homo_sapiens_ancestor_GRCh38.fa.gz \
--plugin OpenTargets,file=/opt/vep/.vep/Plugins/OTGenetics.tsv.gz,cols=all \
--plugin LoFtool,/opt/vep/.vep/Plugins/LoFtool_scores.txt \
--plugin VARITY,file=/opt/vep/.vep/Plugins/varity_all_predictions_38.tsv.gz \
--plugin BayesDel,file=/opt/vep/.vep/Plugins/BayesDel_170824_addAF_all_scores_GRCh38_sorted.txt.gz \
--vcf \
-o /opt/vcf/$OUTPUT_VCF \
--compress_output bgzip \
--allele_number \
--verbose

