#!/bin/bash
## Run DGE per migration permutation column index

## Call env variable to iterate
INDEX=${LSB_JOBINDEX}

## Run script per index 
Rscript ./02_permuted_migration_DGE.R ${INDEX}
