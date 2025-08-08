# add needed metadata for tensorQTL

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
  treatment = args[1]
  pcnumber = as.numeric(args[2])
}

# pcnumber=15
# treatment = "IFN"

#### load genotype PCs and subset to donors present in each condition
message("Adding genotype PCs")

full_genotype_pcs = read.table("../../../OTAR2065_sc_eQTL/data/genotype/plink_genotypes/all_pools.no_outliers.genotype.MAF05.eigenvec")
rownames(full_genotype_pcs) = full_genotype_pcs$V2 # should be line names, or donor names if line was not available

full_genotype_pcs = full_genotype_pcs[c(-1,-2)]
colnames(full_genotype_pcs) = paste0("genotypePC",1:ncol(full_genotype_pcs))

if(treatment !="IFN"){
  

# metadata
metadata = list()
scaled = list()
for(condition in c("Not_proliferating","Proliferating")){
  metadata[[condition]] = read.table(paste0(output_dir,"/metadata_noPCs_",
                                            treatment, "_", condition,".txt"),
                                     header = TRUE,sep = "\t") %>%
    dplyr::mutate(one_over_ncells = 1/count)
  
  print(metadata[[condition]]$donor_id)



# reading expression

  scaled[[condition]] = read.table(paste0(output_dir,"/expr_sum_sizefactorsNorm_log2_scaled_centered_",
                                          treatment, "_", condition,".bed.gz"))
  colnames(scaled[[condition]]) = c("chr","start","end","gene_id",metadata[[condition]]$donor_id)
  



# same number of donors:
nrow(metadata[[condition]]) == ncol(scaled[[condition]][-1:-4])


####### calculating PCS #####
message("Calculating expression PCs")

  # Subset SingleCellExperiment to contain only genes expressed
  # Calculate PCAs on those genes and samples
  numeric_df = scaled[[condition]] %>%
    dplyr::select(!c(chr,start,end,gene_id)) %>%
    as.data.frame() %>%
    mutate_all(as.numeric)
  
  pcs=prcomp(numeric_df,center = FALSE) #already centered
  metadata[[condition]] = t(as.data.frame(metadata[[condition]]))
  
  if(pcnumber > ncol(pcs$rotation)){
    warning("pcnumber is larger than the number of columns in pcs$rotation.")
    
    break
  }else{
    metadata[[condition]] =rbind(t(pcs$rotation[,1:pcnumber]),metadata[[condition]])
  }


head(metadata[[condition]])


  reduced_genotype_pc = full_genotype_pcs[match(colnames(metadata[[condition]]),rownames(full_genotype_pcs)),]
  metadata[[condition]] = rbind(metadata[[condition]],t(reduced_genotype_pc))
  

message("Saving metadata")

  ## save metadata (all info)
  write.table(metadata[[condition]],paste0(output_dir,"/",pcnumber,"/full_metadata_scaled_",
                                           treatment, "_", condition,".txt"), sep = "\t", quote = F, col.names = T, row.names = T)
  # reduced for tensorQTL
  rows_to_keep = rownames(metadata[[condition]]) %in%  c(paste0("PC",1:pcnumber),
                                                        paste0("genotypePC",1:5),
                                                        "one_over_ncells")
    # low variability among donors, setting to first 5 to avoid overcorrection
  write.table(  metadata[[condition]][rows_to_keep,],
                paste0(output_dir,"/",pcnumber,"/tensorQTL_metadata_sum_sizefactorsNorm_log2_",treatment,"_",condition,".txt"), 
                sep = "\t", quote = F, col.names = T, row.names = T)
  
}
}else{
  
  # metadata
  metadata = list()
  scaled = list()
  for(condition in c("Not_proliferating")){
    metadata[[condition]] = read.table(paste0(output_dir,"/metadata_noPCs_",
                                              treatment, "_", condition,".txt"),
                                       header = TRUE,sep = "\t") %>%
      dplyr::mutate(one_over_ncells = 1/count)
    
    print(metadata[[condition]]$donor_id)
    
    
    
    # reading expression
    
    scaled[[condition]] = read.table(paste0(output_dir,"/expr_sum_sizefactorsNorm_log2_scaled_centered_",
                                            treatment, "_", condition,".bed.gz"))
    colnames(scaled[[condition]]) = c("chr","start","end","gene_id",metadata[[condition]]$donor_id)
    
    
    
    
    # same number of donors:
    nrow(metadata[[condition]]) == ncol(scaled[[condition]][-1:-4])
    
    
    ####### calculating PCS #####
    message("Calculating expression PCs")
    
    # Subset SingleCellExperiment to contain only genes expressed
    # Calculate PCAs on those genes and samples
    numeric_df = scaled[[condition]] %>%
      dplyr::select(!c(chr,start,end,gene_id)) %>%
      as.data.frame() %>%
      mutate_all(as.numeric)
    
    pcs=prcomp(numeric_df,center = FALSE) #already centered
    metadata[[condition]] = t(as.data.frame(metadata[[condition]]))
    
    if(pcnumber > ncol(pcs$rotation)){
      warning("pcnumber is larger than the number of columns in pcs$rotation.")
      
      break
    }else{
      metadata[[condition]] =rbind(t(pcs$rotation[,1:pcnumber]),metadata[[condition]])
    }
    
    
    head(metadata[[condition]])
    
    
    reduced_genotype_pc = full_genotype_pcs[match(colnames(metadata[[condition]]),rownames(full_genotype_pcs)),]
    metadata[[condition]] = rbind(metadata[[condition]],t(reduced_genotype_pc))
    
    
    message("Saving metadata")
    
    ## save metadata (all info)
    write.table(metadata[[condition]],paste0(output_dir,"/",pcnumber,"/full_metadata_scaled_",
                                             treatment, "_", condition,".txt"), sep = "\t", quote = F, col.names = T, row.names = T)
    # reduced for tensorQTL
    rows_to_keep = rownames(metadata[[condition]]) %in%  c(paste0("PC",1:pcnumber),
                                                           paste0("genotypePC",1:5),
                                                           "one_over_ncells")
    # low variability among donors, setting to first 5 to avoid overcorrection
    write.table(  metadata[[condition]][rows_to_keep,],
                  paste0(output_dir,"/",pcnumber,"/tensorQTL_metadata_sum_sizefactorsNorm_log2_",treatment,"_",condition,".txt"), 
                  sep = "\t", quote = F, col.names = T, row.names = T)
    
  }
}

## Of note: they don't seem to explicitly account for pool in the covariate matrix, but this variability should be
# accounted for by including the PCs (as explained on the paper)


