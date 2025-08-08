import numpy as np 
import pandas as pd 
import cellex

data = pd.read_csv("../../../data/results/5.prepare_differential_expression_pseudobulks/combined_pseudobulk_raw_counts_treatmentxprolifxidxpool.txt",  sep='\t')
metadata = pd.read_csv("../../../data/results/5.prepare_differential_expression_pseudobulks/metadata_treatmentxprolifxidxpool.txt",  sep='\t')

# Set "gene" column as index (row names)
data.set_index('gene', inplace=True)
# rename metadata
metadata.rename(columns={'cols_names': 'cell_id', 'treatment': 'cell_type'}, inplace=True)
metadata.set_index('cell_id', inplace=True)
# Drop all other columns except 'cell_type'
metadata = metadata[[ 'cell_type']]

eso = cellex.ESObject(data=data, annotation=metadata, verbose=True)

eso.compute(verbose=True)

eso.results["esmu"]

# CELLECT requires that genes are in the Human Ensembl Gene ID format
cellex.utils.mapping.human_symbol_to_human_ens(eso.results["esmu"], drop_unmapped=True, verbose=True)

eso.results["esmu"].to_csv("../../../data/results/CELLECT/cellex_treatment.csv.gz")

