#!/bin/bash
## Run DGE per phagocytosis permutation column index

## Call env variable to iterate
INDEX=${LSB_JOBINDEX}

## Run script per index 
Rscript ./02_permuted_phagocytosis_DGE.R ${INDEX}
