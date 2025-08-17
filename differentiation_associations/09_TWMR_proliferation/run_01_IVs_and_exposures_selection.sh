#!/bin/bash
## Run 01_IVs_and_exposures_selection.R script to create TWMR models for 10 focal genes x treatment x job

## Call env variable to iterate
INDEX=${LSB_JOBINDEX}

## Run script per index 
Rscript ./01_IVs_and_exposures_selection.R ${INDEX}
