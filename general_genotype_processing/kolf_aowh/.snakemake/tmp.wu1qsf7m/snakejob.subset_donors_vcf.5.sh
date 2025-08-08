#!/bin/sh
# properties = {"type": "single", "rule": "subset_donors_vcf", "local": false, "input": ["../../data/kolf_aowh/kolf_aowh_donorIds.txt", "../../data/full_genotype/chr3/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.20180102.genotypes.chr3.vcf.gz"], "output": ["../../data/kolf_aowh/chr3/kolf_aowh.genotype.vcf.gz", "../../data/kolf_aowh/chr3/kolf_aowh.genotype.vcf.gz.csi"], "wildcards": {"pool": "kolf_aowh", "chromosomes": "3"}, "params": {"group": "-G teamtrynka", "queue": "-q normal", "threads": "-n 32", "memory": "-M100000 -R'span[hosts=1] select[mem>100000] rusage[mem=100000]'", "jobname": "-o ../../logs/log_kolf_aowh_3_subset_donors_vcf.%J.%I", "error": "-e ../../errors/error_kolf_aowh_3_subset_donors_vcf.%J.%I"}, "log": [], "threads": 1, "resources": {}, "jobid": 5, "cluster": {}}
cd /lustre/scratch123/hgi/mdt1/projects/otar2065/hipsci_genotype_processing/code/kolf_aowh && \
/software/teamtrynka/conda/trynka-base/bin/python \
-m snakemake ../../data/kolf_aowh/chr3/kolf_aowh.genotype.vcf.gz --snakefile /lustre/scratch123/hgi/mdt1/projects/otar2065/hipsci_genotype_processing/code/kolf_aowh/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /lustre/scratch123/hgi/mdt1/projects/otar2065/hipsci_genotype_processing/code/kolf_aowh/.snakemake/tmp.wu1qsf7m ../../data/kolf_aowh/kolf_aowh_donorIds.txt ../../data/full_genotype/chr3/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.20180102.genotypes.chr3.vcf.gz --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules subset_donors_vcf --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch "/lustre/scratch123/hgi/mdt1/projects/otar2065/hipsci_genotype_processing/code/kolf_aowh/.snakemake/tmp.wu1qsf7m/5.jobfinished" || (touch "/lustre/scratch123/hgi/mdt1/projects/otar2065/hipsci_genotype_processing/code/kolf_aowh/.snakemake/tmp.wu1qsf7m/5.jobfailed"; exit 1)

