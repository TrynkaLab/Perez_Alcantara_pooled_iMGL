## adding treatment interactions for tensorQTL
# merging metadata and expression information, pairwise for untreated vs LPS or IFN
# and creating dataframe with interaction (treatment) information


library(readr)
library(tidyverse)
library(stringr)

output_dir = "../../data/for_tensorQTL/"
dir.create(output_dir)

args = commandArgs(trailingOnly=TRUE)
message(length(args)," arguments provided")
if (length(args)<2) {
  stop("You need to detail 2 arguments: treatment and PC number.n",
       call. = FALSE)
} else if (length(args) == 2) {
  treatment_combination = args[1]
  pcnumber = as.numeric(args[2])
}


# pcnumber=15
# treatment_combination = "IFN_vs_untreated"

scaled = list()
metadata = list()
for(condition in c("Not_proliferating","Proliferating")){
  metadata[[condition]] = read.table( 
                paste0(output_dir,"/",pcnumber,"/tensorQTL_metadata_sum_sizefactorsNorm_log2_",treatment,"_",condition,".txt"), 
                header = TRUE)
  
  scaled[[condition]] = read.table(paste0(output_dir,"/expr_sum_sizefactorsNorm_log2_scaled_centered_",
                                          treatment, "_", condition,".bed.gz"))
  colnames(scaled[[condition]]) = c("chr","start","end","gene_id",colnames(metadata[[condition]]))
  
  # same number of donors
  ncol(metadata[[condition]]) == ncol(scaled[[condition]][-1:-4])
  
}


