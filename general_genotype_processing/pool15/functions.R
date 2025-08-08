# all functions

# Function to convert gt vcf to minor allele dosages
process_gt_dosages = function(vcf){
  # Excluding indels (subsetting only to SNPs) - we may want to keep this info for calculating the proportions?
  # vcf2 <- extract.indels(vcf)
  
  ## Extracting GT
  gt <- as.data.table(extract.gt(vcf))
  # first make new names from chromosome, position, ref and alt (so there are no duplicates)
  names = paste0(vcf@fix[,"CHROM"],"_",vcf@fix[,"POS"],"_",vcf@fix[,"REF"],"_",vcf@fix[,"ALT"])
  
  # there are still some duplicates so remove those
  gt=gt[!duplicated(names),]
  message("Removing ",sum(duplicated(names)), " duplicated rows from vcf")
  
  ## removing those from the names
  names=names[!duplicated(names)]
  gt[, ("rn") := names] 
  
  # setkey(gt,rn) # this sorts the data.table, be careful
  object.size(gt)
  
  ### Some meanings: ################
  # / : genotype unphased (e.g. 0/0)
  #  | : genotype phased (e.g. 0|0)
  # Phased data are ordered along one chromosome and so from these data you know the haplotype
  # i.e. the variants are known to be on the same homologous chromosome because some reads were found to carry both variants
  # Unphased data are simply the genotypes without regard to which one of the pair of chromosomes holds that allele
  ################
  
  
  ## converting GT to minor allele score
  
  gt2 =  gt %>%
    mutate_all(funs(str_replace_all(., "0\\|0", "0")))
  
  gt2 =  gt2 %>%
    mutate_all(funs(str_replace_all(., "0\\/0", "0")))
  
  gt2 = gt2 %>%
    mutate_all(funs(str_replace_all(., "0\\|1", "0.5")))
  
  gt2 = gt2 %>%
    mutate_all(funs(str_replace_all(., "1\\|0", "0.5")))
  
  gt2 =  gt2 %>%
    mutate_all(funs(str_replace_all(., "0\\/1", "0.5")))
  gt2 =  gt2 %>%
    mutate_all(funs(str_replace_all(., "1\\/0", "0.5")))
  
  gt2 = gt2 %>%
    mutate_all(funs(str_replace_all(., "1\\|1", "1")))
  gt2 = gt2 %>%
    mutate_all(funs(str_replace_all(., "1\\/1", "1")))
  
  # change donor columns to numeric
  cols = setdiff(colnames(gt2),"rn")
  gt2[ ,(cols) := lapply(.SD, as.numeric),.SDcols = cols] 
  
  # hist(gt2)
  ## Some donors have NAs at some SNP positions
  anyNA(gt2)
  gt2[which(rowSums(is.na(gt2)) != 0)[1], ]
  gt[which(rowSums(is.na(gt2)) != 0)[1], ]
  
  # Removing SNPs with NAs in at least one donor - do I need to do this? Can I ignore the NAs somehow?
  gt2 = gt2[which(rowSums(is.na(gt2)) == 0), ]
  return(gt2)
}


########### Populating list of reads at every SNP and every value of coverage, for given proportion of donors in pool  #######

sample_reads = function(genotype_df ,
                        depth_df,
                        proportion_df) {
  # The minor allele frequency for SNP i (b[i]) in the population can be calculated as:
  # b[i] = Σ[j=1 to j=ndonors](a[ij]*w[j])
  # where a[ij] is the genotype at SNPi for donor j (coded for minor allele: 0, 0.5 or 1, where the last is homoz. for minor)
  # w[j] is the proportion of donor j in the population
  # For every SNP, looks at the calculated number of reads sampled at that position (depth_df).
  # If the read at i == 0 , skips (unsequenced position)
  # If number of reads at i > 0, for every read:
  # This function calculates and returns a list of
  # b
  # reads: minor (1) allele with probability b[i], otherwise major allele (0) for every SNP, for every set of reads sampled per SNP
  # b_estimate: proportion of minor allele in the sampled reads per SNP - the estimate of b we will use in the real-case scenario
  
  
  ## calculate b for every trial of donor proportions
  message("...Calculating b...")
  
  sum_df =  vector(mode = "list", length = (ncol(genotype_df)-1))
  
  for (j in setdiff(colnames(genotype_df),"rn")) {
    sum_df[[j]] = genotype_df[[j]] * proportion_df[grep(j, colnames(genotype_df))]
    
  }
  length(sum_df[[j]])
  sum_df = as.data.table(do.call("cbind", sum_df))

  # Σ[j=1 to j=ndonors]
  b = sum_df[, rowSums(.SD)]
  
  # sample number of reads detailed in depth_df with probability b
  ## need to do this for every SNP, trial and coverage
  # size is the number of trials, so binom(1, 10, 0.2) has the same expected value as sum(rbinom(10, 1, 0.2))
  
  message("...Sampling reads...")
  
  # b_estimate_subset = vector(mode = "list", length = ncol(depth_df))
  # 
  # message("......Sampling from binomial distribution...")
  # 
  # for(C in setdiff(colnames(depth_df),"rn")) {
  #   ## fix rbinom so that n=0 is converted from integer(0) to NA
  #   rbinom_NA = function(x , y) {
  #     result = rbinom(n = x,
  #                     ## one integer, per C per snp e.g. depth_df[snp,C]
  #                     size = 1,
  #                     prob = y) # one numeric,  per SNP e.g. b[snp]
  #     
  #     if (length(result) == 0)
  #       return(NA)
  #     return(result)
  #   }
  #   # run for all trials and coverages
  #   reads_per_SNP = vector(mode = "list", length = nrow(depth_df))
  #   reads_per_SNP = mcmapply(
  #     depth_df[[C]],
  #     b,
  #     FUN = rbinom_NA, 
  #     mc.cores = detectCores()-1) # this will not work on windows for mc.cores>1
  #   
  #   names(reads_per_SNP) = depth_df$rn
  #   
  #   ## calculate b_estimate: proportion of minor allele in the sampled population per SNP
  #   # It's the arithmetic mean
  # 
  #   b_estimate_subset[[C]] = vapply(reads_per_SNP, mean, FUN.VALUE = 1)
  # }
  # b_estimate = do.call("cbind", b_estimate_subset)
  # b_estimate = as.data.table(b_estimate,keep.rownames = T)
  # 
  #optimization
  coverages=setdiff(colnames(depth_df),"rn")
  get_rbinom = function(C){
    return(rbinom(length(b), depth_df[[C]], b)/depth_df[[C]])
  }
  b_estimate = mclapply(X=coverages,FUN=get_rbinom, mc.cores = min(length(coverages), detectCores() - 1))
  b_estimate = do.call("cbind",b_estimate)
  colnames(b_estimate) = coverages
  # original simulation for optimization
  # fSim <- function(lambda, n){
  #   b <- runif(n)
  #   d <- rpois(n, lambda)
  #   return(rbinom(length(b), d, b)/d)
  # }
  # 
  # lambdas <- c(0.1, 5,10000)
  # system.time(final_result <- mcmapply(fSim, lambdas, 8e6, mc.cores = min(length(lambdas), detectCores() - 1)))
  # 
  gc()
  
  # time for 1 million snps:
  
  message("...Collating results...")
  
  sampled_reads = list(b = b,
                       b_estimate = b_estimate)

  return(sampled_reads)
}

##### generating proportions of donors for every trial

simulate_n_proportions_balance = function(genotype ) {
  ## Calculate a matrix of proportions - these will be the true proportions we are trying to infer later
  # Assuming we sequence a total of [seqcells], we randomly sample from the uniform distribution
  # where the minimum is (seqcells / ndonors) / 2
  # and the maximum is (seqcells/ndonors) * 1.5
  # This ensures that the proportions are reasonably balanced between donors, irrespective of the number of donors in the pool
  seqcells = 10000
  mat_prop = matrix(floor(runif(
    (ncol(genotype)-1),
    min = (seqcells / (ncol(genotype)-1)) / 2,
    max = (seqcells / (ncol(genotype)-1)) * 1.5
  )),
  ncol = (ncol(genotype)-1))
  
  colnames(mat_prop) = setdiff(colnames(genotype),"rn")
  message("There is a mean of ", mean(rowSums(mat_prop)), " cells sequenced")
  mat_prop = proportions(mat_prop, margin = 1)
  message("The average proportion of cells is ", mean(mat_prop))
  
  return(mat_prop)
}

simulate_n_proportions_imbalance = function(genotype) {
  ## Calculate a matrix of proportions - these will be the true proportions we are trying to infer later
  # Starting from a total of [seqcells] (this number doesn't matter that much), we randomly sample from the uniform distribution
  ## With imbalance: so 3 donors take up from 10 to 90% of cells approx
  ## then random sample from uniform distribution from remaining cells as before
  # where the minimum is (seqcells / ndonors) / 2
  # and the maximum is (seqcells/ndonors) * 1.5
  # This ensures that the proportions are reasonably balanced between donors, irrespective of the number of donors in the pool
  seqcells = 10000
  prop_list = list()
  cols = ncol(genotype)-1

    # Fix a proportion between 0.1 and 0.9
    imbalanced = floor(runif(3,
                             min = seqcells * 0.1,
                             max = seqcells * 0.9))
    ## rest at random from distribution of remaining cells
    other = floor(runif(
      (cols - 3),
      min = ((seqcells - imbalanced) / (cols - 2)) / 2,
      max = ((seqcells - imbalanced) / (cols - 2)) * 1.5
    ))
    ## shuffle
    mat_prop= sample(c(imbalanced , other))

  names(mat_prop) = setdiff(colnames(genotype),"rn")
  message("There is a mean of ", mean(mat_prop), " cells sequenced")
  mat_prop = proportions(mat_prop)
  message("The average proportion of cells is ", mean(mat_prop))
  
  return(mat_prop)
  
}



simulate_proportions_fixed_steps_Nsamples_3large= function(ndonors=30,
                                                           start_prop = 0.0005, 
                                                           small_step = 0.002, 
                                                           large_step=0.1){
  prop_list = list()
   
  for (n in 1:ndonors) {
    
    if (n == 1){
      prop_list[[n]] = start_prop
    } 
    if(n >1 & n<(ndonors-2)){
      prop_list[[n]]  = prop_list[[n-1]]+ small_step
    }
    if(n >= (ndonors-2) & n < ndonors){
      prop_list[[n]]  = prop_list[[n-1]]+ large_step
    }
    if(n == ndonors){
      prop_list[[n]]  = 1 - sum(unlist(prop_list))
    }
  }
  if (sum(unlist(prop_list))>1){ 
    stop("The cumulative proportion is > 1. Consider smaller steps or different number of donors.")
  }
  
  return(unlist(prop_list))
  
    
  }
  
# estimate weights
estimate_weights = function(b, A=gt){
  ## Solve for w in Aw = b
  ## b in our case is a vector b_est (b in best case scenario, calculated directly from the genotype file and known weights), 
  # for each trial and coverage
  ## A is a matrix of rows = SNPs (m) and columns = donors (k)
  # a[mk] is 0, 0.5 or 1 (minor allele dosage)
  # w is a vector of proportions with length k, where the sum equals 1 

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
  
  
  return(w)
}

# estimate w per trial and coverage
calculate_w_trial_coverage = function(b_estimate, .gt = gt) {
  w_trials = vector(mode = "list", length = trials)
  for (t in 1:trials) {
    
    message(paste0("....... Working on trial ",t, " .........."))
    cols = setdiff(colnames(b_estimate),c("trial","rn"))
    b = b_estimate[trial == t,..cols]
    
    w_trials[[t]]= apply(b, 2, FUN = function(x) estimate_weights(b=x))
    
    # reformat
    w_trials[[t]] =  data.frame(w_trials[[t]])
    rownames(w_trials[[t]]) = setdiff(colnames(gt),c("rn"))
    w_trials[[t]] = as.data.frame(t(w_trials[[t]] )) # transpose for longer format after unlist
    w_trials[[t]]$coverage =  cols
    w_trials[[t]]$trial = t
    gc() # Force garbage collection
    # takes 20 seconds per trial for all coverages
  }
  w_trials = do.call("rbind", w_trials)
  return(w_trials)
}