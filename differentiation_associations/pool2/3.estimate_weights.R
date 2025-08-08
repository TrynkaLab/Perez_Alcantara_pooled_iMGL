# Estimate weights
# Input: b_estimate file and genotype vcf (reduced, from step 2)
# Output: estimated donor proportions (w_estimate)

library(pracma)
library(tidyverse)
library(data.table)
options(stringsAsFactors = FALSE)




# For server
args = commandArgs(trailingOnly=TRUE)


if (length(args)<3) {
  stop("You need to detail at least 2 input and 1 output file paths.n",
       call. = FALSE)
} else if (length(args) == 3) {
  b_estimate_path = args[1]
  genotype_minor_allele_dos_path = args[2]
  w_estimate_path = args[3]
}

message("...Reading minor allele dosage file...\n")

genotype = fread(genotype_minor_allele_dos_path)

message("...Reading b_estimate file ...\n")
b_estimate=fread(b_estimate_path)


estimate_weights = function(b, A=gt){
  ## Solve for w in Aw = b
  ## b in our case is a vector b_est (b in best case scenario, calculated directly from the genotype file and known weights), 
  # for each trial and coverage
  ## A is a matrix of rows = SNPs (m) and columns = donors (k)
  # a[mk] is 0, 0.5 or 1 (minor allele dosage)
  # w is a vector of proportions with length k, where the sum equals 1 
  # tic("Function starts ... ")
  
  # check dimensions
  length(b) == nrow(A)
  # ideally check with snp names if they are identical()
  
  # remove rows from A where all genotypes are identical - those are where variance per row equals 0
  cols = setdiff(colnames(A),"rn")
  
  is_identical_genotype <- function(x, ...) {
    # rowVar() == 0, faster version of apply(A[,..cols], MARGIN = 1, FUN = function(x) var(x) == 0)
    (rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)) == 0
  }
  to_remove =  is_identical_genotype(A[,..cols])
  # message(paste0("We are removing ", sum(to_remove), " identical rows out of ", nrow(A), 
  #                ": the ", floor((sum(to_remove)/nrow(A))*100), "%" ))
  A = A[!to_remove,]
  
  # remove same SNPs from b
  b = b[!to_remove]
  
  # Removing SNPs with NAs - not sampled (not "sequenced") in the tested coverage for this trial
  to_remove = is.na(b)
  # message(paste0("We are removing ", sum(to_remove), " unsampled SNPs out of the remaining ", length(b)))
  b = b[!to_remove]
  A = A[!to_remove,]
  # solving for w
  #  Aeqw=beq  where Aeq  is a 1×k matrix of ones and be is a length m -vector, also of ones, contrains the sum of weights to 1
  #  lb=0 and ub=1  constraints our solution to wi∈[0,1]
  
  w=lsqlincon(C=as.matrix(A[,..cols]),
              d=b,
              Aeq=matrix(1,ncol=length(cols),nrow=1),
              beq=1,
              lb=0,ub=1)
  
 
  # toc()
  return(w)
}

w = estimate_weights(b=b_estimate$b_estimate, A=genotype)

names(w) = setdiff(colnames(genotype),"rn")

write.table(w,file = w_estimate_path, col.names = F)

