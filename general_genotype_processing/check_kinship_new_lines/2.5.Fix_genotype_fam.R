# recoding genotype fam file
.libPaths(c('/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/',"/software/teamtrynka/ma23/R4.1/libs",
            "/software/R-4.1.0/lib/R/library"))
library(tidyverse)
library(stringr)
args = commandArgs(trailingOnly=TRUE)
input.path = args[1]
output.path = args[2]

input = read_table(input.path,col_names = FALSE) %>%
  dplyr::mutate(X1=stringr::str_split_fixed(X2,pattern = "_",n=2)[,1]) 

write.table(input,
            file = output.path,
            row.names = FALSE,col.names = FALSE,sep = "\t", quote = FALSE)

