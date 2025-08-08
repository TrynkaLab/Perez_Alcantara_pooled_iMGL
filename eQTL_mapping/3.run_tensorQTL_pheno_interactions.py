# Running tensorQTL with phenotypic interactions
import os
import os.path
import getopt
import sys
import time
import re # provides regular expression matching operations similar to those found in Perl
import subprocess #keeps command output
import tensorqtl
from tensorqtl import genotypeio, cis, trans
import numpy as np
import pandas as pd

# Load covariate and phenotype
expression_bed = snakemake.input.phenotype
covariates_file = snakemake.input.covariate
interactions_file = snakemake.input.interaction
plink_genotypes_loc = snakemake.params.plink_files
tensor_out = snakemake.params.out_file

print("Phenotype file: ")
print(expression_bed)
print("Covariate file: ")
print(covariates_file)
print("Output file: ")
print(tensor_out)

phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T
interaction_df = pd.read_csv(interactions_file, sep='\t', index_col=0)
print("Phenotype matrix shape: ")
print(phenotype_df.columns)
print("Covariates matrix shape: ")
print(covariates_df.index)
print("Interactions matrix shape: ")
print(interaction_df.index)
# Load genotypes
pr = genotypeio.PlinkReader(plink_genotypes_loc)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos', 'a0', 'a1']]

variant_df_str = variant_df
variant_df_str['chrom'] = variant_df['chrom'].astype(str)
variant_df_str['pos'] = variant_df['pos'].astype(str)
variant_df['hg38'] = variant_df_str[['chrom', 'pos', 'a1', 'a0']].apply(lambda x: '_'.join(x), axis=1)
del variant_df_str

variant_df.index = variant_df.hg38
genotype_df.index = variant_df.index

#remove donors not present in the expression data
genotype_df_subset = genotype_df[list(phenotype_df.columns)]

# Changing position from string to float
variant_df=variant_df.astype({'pos': 'float64'})


# Extract everything up to the last "/"
path = os.path.dirname(tensor_out)

# Extract everything after the last "/", excluding the file extension
prefix = os.path.splitext(os.path.basename(tensor_out))[0]

# Print the extracted values
print("Directory:", path)
print("Prefix:", prefix)

cis.map_nominal(genotype_df=genotype_df_subset,
variant_df=variant_df,
phenotype_df=phenotype_df, phenotype_pos_df=phenotype_pos_df,
prefix=prefix,
covariates_df=covariates_df,
output_dir=path,
interaction_df=interaction_df,
maf_threshold_interaction=0.05,
                run_eigenmt=True, # True throws deprecation error
               write_top=True, 
               write_stats=True) # True causes error in chr2 - potentially because of wrong sorting, should be fixed in new tensorqtl version
                

