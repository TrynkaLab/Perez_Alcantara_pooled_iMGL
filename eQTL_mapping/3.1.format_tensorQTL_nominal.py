#format nominal pass
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

tensor_in = snakemake.input.nominal
tensor_out = snakemake.output.formated

# Extract everything up to the last "/"
path = os.path.dirname(tensor_out)

# Extract everything after the last "/", excluding the file extension
prefix = os.path.splitext(os.path.basename(tensor_out))[0]
# Print the extracted values
print("Directory:", path)
print("Prefix:", prefix)

# load results
print("Loading results  ")

chromosomes = range(1, 23)  # Iterate through chromosome numbers 1 to 22

dataframes = []  # List to store DataFrames for each chromosome

for chromosome in chromosomes:
    file_path = f"{path}/{prefix}.cis_qtl_pairs.{chromosome}.parquet"  # Construct the file path for each chromosome

    # Read the parquet file and append the DataFrame to the list
    dataframe = pd.read_parquet(file_path)
    dataframes.append(dataframe)

# Concatenate all DataFrames into a single DataFrame
pairs_df = pd.concat(dataframes)
pairs_df.to_csv(tensor_out, sep='\t')
