#!/bin/bash
## Run DGE per proliferation permutation column index

## Call env variable to iterate
INDEX=${LSB_JOBINDEX}

## Run script per index 
Rscript ./02_permuted_proliferation_DGE.R ${INDEX}
