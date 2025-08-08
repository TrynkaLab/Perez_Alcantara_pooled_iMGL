# all functions
# Function to convert gt vcf to minor allele dosages
process_gt_dosages = function(vcf){
  # Excluding indels (subsetting only to SNPs) - we may want to keep this info for calculating the proportions?
  # vcf2 <- extract.indels(vcf)
  
  ## Extracting GT
  gt <- as.data.table(vcfR::extract.gt(vcf))
  # first make new names from chromosome, position, ref  (so there are no duplicates)
  names = paste0(vcf@fix[,"CHROM"],"_",vcf@fix[,"POS"],"_",vcf@fix[,"REF"])
  
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
    dplyr::mutate_all(funs(str_replace_all(., "0\\|0", "0")))
  
  gt2 =  gt2 %>%
    dplyr::mutate_all(funs(str_replace_all(., "0\\/0", "0")))
  
  gt2 = gt2 %>%
    dplyr::mutate_all(funs(str_replace_all(., "0\\|1", "0.5")))
  
  gt2 = gt2 %>%
    dplyr::mutate_all(funs(str_replace_all(., "1\\|0", "0.5")))
  
  gt2 =  gt2 %>%
    dplyr::mutate_all(funs(str_replace_all(., "0\\/1", "0.5")))
  gt2 =  gt2 %>%
    dplyr::mutate_all(funs(str_replace_all(., "1\\/0", "0.5")))
  
  gt2 = gt2 %>%
    dplyr::mutate_all(funs(str_replace_all(., "1\\|1", "1")))
  gt2 = gt2 %>%
    dplyr::mutate_all(funs(str_replace_all(., "1\\/1", "1")))
  
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

# same as above, but returns SNPs named as ref alt (there might be two SNPs at same position)
process_gt_dosages_ref_alt = function(vcf){
  # Excluding indels (subsetting only to SNPs) - we may want to keep this info for calculating the proportions?
  # vcf2 <- extract.indels(vcf)
  
  # solve possible issue with duplicated SNP IDs
myID = getID(vcf)
 if((length(unique(myID, incomparables = NA)) == length(myID)) ==FALSE){
  message("WARNING: Fixing issue with duplicated IDs")
   message("Example of duplicated ID")
   dups = which(duplicated(vcf@fix[,3]))
   if(length(dups) >=1){
   print(vcf@fix[vcf@fix[,3]==vcf@fix[dups[1],3],])
   }
  vcf3 = vcf[!duplicated(myID, incomparables = NA), ]
  myID <- getID(vcf3)
  length(unique(myID, incomparables = NA)) == length(myID)
  vcf = vcf3

 }
## Extracting GT

gt <- as.data.table(vcfR::extract.gt(vcf))


  

  # first make new names from chromosome, position, ref, alt (there might be duplicates)
  names = paste(vcf@fix[,"CHROM"],vcf@fix[,"POS"],vcf@fix[,"REF"], vcf@fix[,"ALT"],sep = "_")
  
  # are there more duplicates?
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
    dplyr::mutate_all(funs(str_replace_all(., "0\\|0", "0")))
  
  gt2 =  gt2 %>%
    dplyr::mutate_all(funs(str_replace_all(., "0\\/0", "0")))
  
  gt2 = gt2 %>%
    dplyr::mutate_all(funs(str_replace_all(., "0\\|1", "0.5")))
  
  gt2 = gt2 %>%
    dplyr::mutate_all(funs(str_replace_all(., "1\\|0", "0.5")))
  
  gt2 =  gt2 %>%
    dplyr::mutate_all(funs(str_replace_all(., "0\\/1", "0.5")))
  gt2 =  gt2 %>%
    dplyr::mutate_all(funs(str_replace_all(., "1\\/0", "0.5")))
  
  gt2 = gt2 %>%
    dplyr::mutate_all(funs(str_replace_all(., "1\\|1", "1")))
  gt2 = gt2 %>%
    dplyr::mutate_all(funs(str_replace_all(., "1\\/1", "1")))
  
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

#read w estimates per pool
# fix for migration
read_w_estimates = function(sample_info,path,assay="migration",pools = paste0("pool",2:6)){
  w_estimates = list()
  for(pool in pools){
    message("Reading in ",pool)
    for(s in sample_info[sample_info$pool==pool,"sample_name"]){
      message("...sample ",s)
      
      w_estimates[[paste(pool,s,sep = "_")]] = read.table(paste0(path,pool,"/",assay,"/",s,"_w_estimate.txt"))
      colnames(w_estimates[[paste(pool,s,sep = "_")]]) = c("donor","proportion")
      w_estimates[[paste(pool,s,sep = "_")]]$pool = pool
      w_estimates[[paste(pool,s,sep = "_")]]$sample_name = s
      
    }
  }
  w_estimates = do.call("rbind",w_estimates)
  w_estimates = dplyr::inner_join(w_estimates,sample_info, by = c("sample_name","pool"))
  return(w_estimates)
}

# Find the closest value in Y to a given value in X
find_closest_value = function(x, y) {
  y[which.min(abs(y - x))]
}

# phenotypic QTL boxplots
boxplot_pQTL = function(variant,phenotype_with_genotype_info_df,qval,coef,color_by="pool", phenotype="migration"){
  
  # color_by = "pool" or "well" or "none"
  message("Plotting variant ",variant)
  
  toplot = phenotype_with_genotype_info_df %>%
    dplyr::filter(rn == variant) %>%
    dplyr::mutate(REF = str_split_fixed(rn,"_",n = 4)[,3],ALT = str_split_fixed(rn,"_", n=4)[,4]) %>%
    dplyr::mutate(genotype = case_when(alt_allele_dosage == 0.0 ~ paste0(REF,"|", REF),
                                       alt_allele_dosage == 0.5 ~ paste0(REF,"|", ALT),
                                       alt_allele_dosage == 1.0 ~ paste0(ALT,"|", ALT))) %>%
    dplyr::mutate(genotype = factor(genotype,levels = c(unique(paste0(REF,"|", REF)),
                                                        unique(paste0(REF,"|", ALT)),
                                                        unique(paste0(ALT,"|", ALT))), 
                                    ordered = TRUE))
  
  
  if(phenotype=="migration"){
    ylab_message = "log(bottom proportion / top proportion)"
  }else{
    ylab_message = "log(mCherry+ proportion/mCherry- proportion)"
  }
  # Create the plot with box plot and points
  
  if(color_by == "pool"){
    p = toplot %>% 
      ggplot(aes(x = genotype, y = mean_outcome_pool)) +
      geom_boxplot() + # Create the box plot
      geom_jitter(aes( col=pool), alpha = 0.6,width = 0.2) + # Create points with size based on 'size_var'
      theme_classic() + 
      ylab(ylab_message)
    
    p = p + patchwork::plot_annotation(title = paste0("Allelic effects at SNP ", variant),
                                       subtitle = paste0("In ", unique(toplot$treatment), 
                                                         " + ",unique(toplot$condition),
                                                         ". Coefficient = ",round(coef,2),
                                                         ". Adjusted p-val = ",signif(qval,3))
    )
    
  } 
  if(color_by == "well"){
    p = toplot %>% 
      ggplot(aes(x = genotype, y = mean_outcome_pool)) +
      geom_boxplot() + # Create the box plot
      geom_jitter(aes( col=well, size = prop_unadjusted_max_value), alpha = 0.6,width = 0.2) + # Create points with size based on 'size_var'
      theme_classic() + 
      ylab(ylab_message)
    
    p = p + patchwork::plot_annotation(title = paste0("Allelic effects at SNP ", variant),
                                       subtitle = paste0("In ", unique(toplot$treatment), 
                                                         " + ",unique(toplot$condition),
                                                         ". Coefficient = ",round(coef,2),
                                                         ". Adjusted p-val = ",signif(qval,3))
    )
    
  } 
  if(color_by == "none"){
    
    p = toplot %>% 
      ggplot(aes(x = genotype, y = mean_outcome_pool)) +
      geom_boxplot() + # Create the box plot
      geom_jitter(alpha = 0.6,width = 0.2) + # Create points with size based on 'size_var'
      theme_classic() + 
      ylab(paste0("log(bottom proportion / top proportion) "))
    
    p = p + patchwork::plot_annotation(title = paste0("Allelic effects at SNP ", variant),
                                       subtitle = paste0("In ", unique(toplot$treatment), 
                                                         " + ",unique(toplot$condition),
                                                         ". Coefficient = ",round(coef,2),
                                                         ". Adjusted p-val = ",signif(qval,3))
    )
  }
  
  
  return(p)
  
}

boxplot_pQTL_noaverage = function(variant,phenotype_with_genotype_info_df,qval,coef,color_by="pool",assay=NA){
  # color_by = "pool" or "well" or "none"
  message("Plotting variant ",variant)
  
  toplot = phenotype_with_genotype_info_df %>%
    dplyr::filter(rn == variant) %>%
    dplyr::mutate(REF = str_split_fixed(rn,"_",n = 4)[,3],ALT = str_split_fixed(rn,"_", n=4)[,4]) %>%
    dplyr::mutate(genotype = case_when(alt_allele_dosage == 0.0 ~ paste0(REF,"|", REF),
                                       alt_allele_dosage == 0.5 ~ paste0(REF,"|", ALT),
                                       alt_allele_dosage == 1.0 ~ paste0(ALT,"|", ALT))) %>%
    dplyr::mutate(genotype = factor(genotype,levels = c(unique(paste0(REF,"|", REF)),
                                                        unique(paste0(REF,"|", ALT)),
                                                        unique(paste0(ALT,"|", ALT))), 
                                    ordered = TRUE))
  
  # Create the plot with box plot and points
  
  if(color_by == "pool"){
    p = toplot %>% 
      ggplot(aes(x = genotype, y = log_fraction_mean)) +
      geom_boxplot() + # Create the box plot
      geom_jitter(aes( col=pool, size = prop_unadjusted_max_value), alpha = 0.6,width = 0.2) + # Create points with size based on 'size_var'
      theme_classic() + 
      ylab(paste0("log(bottom proportion / top proportion) "))
    
    p = p + patchwork::plot_annotation(title = paste0("Allelic effects at SNP ", variant),
                                       subtitle = paste0(assay,":  ", unique(toplot$treatment), 
                                                        
                                                         ". Coefficient = ",round(coef,2),
                                                         ". Adjusted p-val = ",signif(qval,3))
    )
    
  } 
  if(color_by == "well"){
    p = toplot %>% 
      ggplot(aes(x = genotype, y = log_fraction_mean)) +
      geom_boxplot() + # Create the box plot
      geom_jitter(aes( col=well, size = prop_unadjusted_max_value), alpha = 0.6,width = 0.2) + # Create points with size based on 'size_var'
      theme_classic() + 
      ylab(paste0("log(bottom proportion / top proportion) "))
    
    p = p + patchwork::plot_annotation(title = paste0("Allelic effects at SNP ", variant),
                                       subtitle = paste0("In ", unique(toplot$treatment), 
                                                         " + ",unique(toplot$condition),
                                                         ". Coefficient = ",round(coef,2),
                                                         ". Adjusted p-val = ",signif(qval,3))
    )
    
  } 
  if(color_by == "none"){
    
    p = toplot %>% 
      ggplot(aes(x = genotype, y = log_fraction_mean)) +
      geom_boxplot() + # Create the box plot
      geom_jitter(alpha = 0.6,width = 0.2) + # Create points with size based on 'size_var'
      theme_classic() + 
      ylab(paste0("log(bottom proportion / top proportion) "))
    
    p = p + patchwork::plot_annotation(title = paste0("Allelic effects at SNP ", variant),
                                       subtitle = paste0("In ", unique(toplot$treatment), 
                                                         " + ",unique(toplot$condition),
                                                         ". Coefficient = ",round(coef,2),
                                                         ". Adjusted p-val = ",signif(qval,3))
    )
  }
  
  
  return(p)
  
}

boxplot_pQTL_exploratory_reps = function(variant,phenotype_with_genotype_info_df,color_by="pool", phenotype="migration"){
  
  # color_by = "pool" or "well" or "none"
  message("Plotting variant ",variant)
  
  toplot = phenotype_with_genotype_info_df %>%
    dplyr::filter(rn == variant) %>%
    dplyr::mutate(REF = str_split_fixed(rn,"_",n = 4)[,3],ALT = str_split_fixed(rn,"_", n=4)[,4]) %>%
    dplyr::mutate(genotype = case_when(alt_allele_dosage == 0.0 ~ paste0(REF,"|", REF),
                                       alt_allele_dosage == 0.5 ~ paste0(REF,"|", ALT),
                                       alt_allele_dosage == 1.0 ~ paste0(ALT,"|", ALT))) %>%
    dplyr::mutate(genotype = factor(genotype,levels = c(unique(paste0(REF,"|", REF)),
                                                        unique(paste0(REF,"|", ALT)),
                                                        unique(paste0(ALT,"|", ALT))), 
                                    ordered = TRUE))
  
  
  if(phenotype=="migration"){
    ylab_message = "log(bottom proportion / top proportion)"
  }else{
    ylab_message = "log(mCherry+ proportion/mCherry- proportion)"
  }
  # Create the plot with box plot and points
  
  if(color_by == "pool"){
    p = toplot %>% 
      ggplot(aes(x = genotype, y = log_fraction_mean)) +
      geom_boxplot() + # Create the box plot
      geom_jitter(aes( col=pool), alpha = 0.6,width = 0.2) + # Create points with size based on 'size_var'
      theme_classic() + 
      ylab(ylab_message)
    
    p = p + patchwork::plot_annotation(title = paste0("Allelic effects at SNP ", variant),
                                       subtitle = paste0("In ", unique(toplot$treatment), 
                                                         " ." )
    )
    
  } 
  if(color_by == "well"){
    p = toplot %>% 
      ggplot(aes(x = genotype, y = log_fraction_mean)) +
      geom_boxplot() + # Create the box plot
      geom_jitter(aes( col=well, size = prop_unadjusted_max_value), alpha = 0.6,width = 0.2) + # Create points with size based on 'size_var'
      theme_classic() + 
      ylab(ylab_message)
    
    p = p + patchwork::plot_annotation(title = paste0("Allelic effects at SNP ", variant),
                                       subtitle = paste0("In ", unique(toplot$treatment), 
                                                         " + ",unique(toplot$condition),
                                                         ". " )
    )
    
  } 
  if(color_by == "none"){
    
    p = toplot %>% 
      ggplot(aes(x = genotype, y = log_fraction_mean)) +
      geom_boxplot() + # Create the box plot
      geom_jitter(alpha = 0.6,width = 0.2) + # Create points with size based on 'size_var'
      theme_classic() + 
      ylab(paste0("log(bottom proportion / top proportion) "))
    
    p = p + patchwork::plot_annotation(title = paste0("Allelic effects at SNP ", variant),
                                       subtitle = paste0("In ", unique(toplot$treatment), 
                                                         " + ",unique(toplot$condition),
                                                         "." )
    )
  }
  
  
  return(p)
  
}

nicer_heatmap = function(cor_res){
  row.order = hclust(dist(cor_res$r))$order
  col.order = hclust(dist(t(cor_res$r)))$order
  
  cor_res = cor_res$r[row.order,col.order]
  # color palettes
  col_fun = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  
  ht = Heatmap(cor_res, name = "mat", col = col_fun,
               row_names_gp = gpar(fontsize = 8),
               column_names_gp = gpar(fontsize = 8),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(sprintf("%.1f", cor_res[i, j]), x, y, gp = gpar(fontsize = 10))
               }, cluster_rows = F, cluster_columns = F)
  
  p = grid.grabExpr(draw(ht, padding = unit(c(8, 5, 5, 5), "mm"))) # to save object to grid
  return(p)
  
}

swap_REF_ALT_alelle_names = function(variant_name){
  
  # Split the input string by the underscore character
  string_parts = strsplit(variant_name, "_")[[1]]
  
  # Extract the last element and swap the order of the last two letters
  last_element = string_parts[length(string_parts)]
  new_last_element =string_parts[length(string_parts)-1]
  
  # Combine the parts back together with the modified last element
  new_variant_name = paste0(paste0(string_parts[1:2], collapse = "_"), "_", last_element,"_", new_last_element)
  
  # Print the new string
  return(new_variant_name)
  
}

myManhattan <- function(df, graph.title = "", highlight = NULL, highlight.col = "green",
                        col = c("lightblue", "navy"), even.facet = FALSE, chrom.lab = NULL,
                        suggestiveline = 1e-05, suggestivecolor = "blue",
                        genomewideline = 5e-08, genomewidecolor = "red",
                        font.size = 12, axis.size = 0.5, significance = NULL, report = FALSE,
                        inf.corr = 0.95, y.step = 2, point.size = 1){
  # modified from https://github.com/alfonsosaera/myManhattan
  # altering p-values that are 0 to be x% smaller than the smallest non-zero value
  myMin <- min(df$P[df$P != 0]) * inf.corr
  df$P[df$P == 0] <- myMin
  
  ### here I fixed BP because it was plotting together basepairs from different chromosomes!!!
  
  df = df %>%
    dplyr::arrange(CHR,BP) %>%
    dplyr::mutate(row = 1:nrow(df))
  
  require(ggplot2)
  require(stats)
  # checking colors 
  y.title <- expression(-log[10](italic(p)))
  if (length(col) > length(unique(df$CHR))){
    chrom.col <- col[1:length(unique(df$CHR))]
  } else if (!(length(col) > length(unique(df$CHR)))){
    chrom.col <- rep(col, length(unique(df$CHR))/length(col))
    if (length(chrom.col) < length(unique(df$CHR))){
      dif <- length(unique(df$CHR)) - length(chrom.col)
      chrom.col <- c(chrom.col, col[1:dif])
    }
  }
  # checking min p-val for y axis plotting
  y.max <- floor(max(-log10(df$P))) + 1
  if (y.max %% 2 != 0){
    y.max <- y.max + 1
  }
  # checking chromosome labels if present
  if (!is.null(chrom.lab)){
    if (length(unique(df$CHR)) != length(chrom.lab)){
      warning("Number of chrom.lab different of number of chromosomes in dataset, argument ignored.")
    } else {
      df$CHR <- factor(df$CHR, levels = unique(df$CHR), labels=chrom.lab)
    }
  }
  

  g <- ggplot(df) +
    geom_point(aes(row, -log10(P), colour = as.factor(CHR)), size = point.size)
  if (!is.null(significance)){
    if (is.numeric(significance)){
      genomewideline <- significance
      suggestiveline <- genomewideline / 0.005
    } else if (significance == "Bonferroni"){
      BFlevel <- 0.05 / length(df$SNP)
      cat("Bonferroni correction significance level:", BFlevel, "\n")
      genomewideline <- BFlevel
      suggestiveline <- BFlevel / 0.005
    } else if (significance == "FDR"){
      df$fdr <- p.adjust(df$P, "fdr")
      genomewideline <- 0.05
      suggestiveline <- FALSE
      y.title <- expression(-log[10](italic(q)))
      g <- ggplot(df) +
        geom_point(aes(row, -log10(fdr), colour = as.factor(CHR)), size = point.size)
      if (!is.null(highlight)) {
        if (is.numeric(highlight)){
          highlight <- as.character(df$SNP[df$P < highlight])
        }
        if (any(!(highlight %in% df$SNP))){
          warning("Cannot highlight SNPs not present in the dataset. Argument is ignored.")
        } else {
          g <- g + geom_point(data = df[which(df$SNP %in% highlight), ],
                              aes(row, -log10(fdr), group=SNP, colour=SNP),
                              color = highlight.col, size = point.size)
          highlight <- NULL
          y.max <- floor(max(-log10(df$fdr))) + 1
          if (y.max %% 2 != 0){
            y.max <- y.max + 1
          }
        }
      }
    }
  }
  if (even.facet){
    g <- g + facet_grid(.~CHR, scale = "free_x", switch = "x")
  } else {
    g <- g + facet_grid(.~CHR, scale = "free_x", space = "free_x", switch = "x")
  }
  g <- g + scale_colour_manual(values = chrom.col) +
    scale_y_continuous(expand = c(0, 0), limit = c(0, y.max),
                       breaks = seq(from = 0, to = y.max, by = y.step)) +
    scale_x_continuous() +
    ggpubr::theme_pubr()+
    theme(strip.background = element_blank(), legend.position = "none",
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.spacing.x=unit(0.1, "lines"),
          axis.line.y = element_line(linewidth = axis.size, color = "black"),
          axis.ticks.y = element_line(linewidth = axis.size, color = "black"),
          axis.ticks.length = unit(axis.size * 10, "points"),
          plot.title = element_text(hjust = (0.5), size = font.size + 8),
          axis.title.y = element_text(size = font.size + 5),
          axis.title.x = element_text(size = font.size + 5),
          axis.text = element_text(size = font.size),
          strip.text.x = element_text(size = font.size))+
    labs(title = graph.title, x = "Chromosome", y = y.title)
  if (!is.null(highlight)) {
    if (is.numeric(highlight)){
      highlight <- as.character(df$SNP[df$P < highlight])
    }
    if (any(!(highlight %in% df$SNP))){
      warning("Cannot highlight SNPs not present in the dataset. Argument is ignored.")
    } else {
      g <- g + geom_point(data = df[which(df$SNP %in% highlight), ],
                          aes(row, -log10(P), group=SNP, colour=SNP),
                          color = highlight.col, size = point.size)
    }
  }
  if (suggestiveline){
    g <- g + geom_hline(yintercept = -log10(suggestiveline), color = suggestivecolor)
  }
  if (genomewideline){
    g <- g + geom_hline(yintercept = -log10(genomewideline), color = genomewidecolor)
  }
  if (report){
    if (significance == "FDR"){
      rep <- df[df$fdr < 0.05, ]
    } else if (significance == "Bonferroni"){
      rep <- df[df$P < BFlevel, ]
    } else if (is.numeric(significance)){
      rep <- df[df$P < significance, ]
    } else {
      cat("using default significance level, 5e-8")
      rep <- df[df$P < 5e-8, ]
    }
    print(rep)
  }
  return(g)
}

read_b_estimates = function(sample_info,path,assay="migration",pools = paste0("pool",2:6)){
  b_estimates = list()
  for(pool in pools){
    message("Reading in ",pool)
    for(s in sample_info[sample_info$pool==pool,"sample_name"]){
      message("...sample ",s)
      
      b_estimates[[paste(pool,s,sep = "_")]] = readr::read_csv(paste0(path,pool,"/",assay,"/",s,"_b_estimate.csv"))
      b_estimates[[paste(pool,s,sep = "_")]]$chr_pos_ref = readr::read_csv(paste0(path,"../genotypes/",pool,"/",assay,"/",s,
                                                                                  "_genotype_minor_allele_dosage.csv"
                                                                                  ),
                                                                           col_types = cols_only(rn = 'c'))$rn
      b_estimates[[paste(pool,s,sep = "_")]] = b_estimates[[paste(pool,s,sep = "_")]] %>%
        dplyr::mutate(variant_id = paste(chr_pos_ref, ALT, sep = "_")) %>%
        dplyr::rename(ALT_allele_estimate = b_estimate)
      
      b_estimates[[paste(pool,s,sep = "_")]]$pool = pool
      b_estimates[[paste(pool,s,sep = "_")]]$sample_name = s
      
    }
  }
  b_estimates = do.call("rbind",b_estimates)
  b_estimates = dplyr::inner_join(b_estimates,sample_info, by = c("sample_name","pool"))
  return(b_estimates)
}

# modified pca from PCAtools to avoid useNames = NA error
pca_modified = function (mat, metadata = NULL, center = TRUE, scale = FALSE, 
          rank = NULL, removeVar = NULL, transposed = FALSE, BSPARAM = BiocSingular::ExactParam()) 
{
  if (is.data.frame(mat)) {
    mat <- as.matrix(mat)
  }
  if (!transposed) {
    mat <- t(mat)
  }
  if (!is.null(metadata)) {
    if (!identical(rownames(mat), rownames(metadata))) {
      stop("'colnames(mat)' is not identical to 'rownames(metadata)'")
    }
  }
  .center <- if (center) 
    NULL
  else 0
  vars <- matrixStats::colVars(mat, center = .center)
  if (!is.null(removeVar)) {
    message("-- removing the lower ", removeVar * 100, "% of variables based on variance")
    varorder <- order(vars, decreasing = TRUE)
    keep <- head(varorder, max(1, ncol(mat) * (1 - removeVar)))
    mat <- mat[, keep, drop = FALSE]
    vars <- vars[keep]
  }
  if (is.null(rank)) {
    if (is(BSPARAM, "ExactParam")) {
      rank <- min(dim(mat))
    }
    else {
      stop("'rank' must be specified for approximate PCA methods")
    }
  }
  pcaobj <- BiocSingular::runPCA(mat, center = center, scale = scale, rank = rank, 
                   BSPARAM = BSPARAM)
  if (scale) {
    total.var <- length(vars)
  }
  else {
    total.var <- sum(vars)
  }
  proportionvar <- (pcaobj$sdev^2)/total.var * 100
  pcaobj <- list(rotated = data.frame(pcaobj$x), loadings = data.frame(pcaobj$rotation), 
                 variance = proportionvar, sdev = pcaobj$sdev, metadata = metadata, 
                 xvars = colnames(mat), yvars = rownames(mat), components = colnames(pcaobj$x))
  rownames(pcaobj$rotated) <- pcaobj$yvars
  rownames(pcaobj$loadings) <- pcaobj$xvars
  names(pcaobj$variance) <- pcaobj$components
  class(pcaobj) <- "pca"
  return(pcaobj)
}


fix_burden_matrices <- function(df) {
  df %>%
    tibble::as_tibble(rownames = "gene") %>%
    dplyr::rename_with(~str_extract(., "(?<=-)\\w+"), .cols = contains("H"))
}

extract_burden_results <- function(df) {
  res = data.frame("p_val"= df$coefficients[2,"Pr(>|t|)"], ### careful here to take right covariate
                   "estimate" = df$coefficients[2,"Estimate"],
                   "gene" = df$gene) %>%
    pivot_longer(cols=-gene,names_to = "coefficients",values_to = "value")
  
}

# subset column names to shared lines, including "gene"
subset_burden_shared_lines = function(x,y) { 
  x %>%
    dplyr::select(c("gene",sort(lubridate::intersect(colnames(.),y)))) 
}

filter_rare_variant_genes = function(tb,fraction) { 
  message("The fraction corresponds to ",floor((ncol(tb)-1)*fraction), " lines")
  tb %>%
    tibble::column_to_rownames(var="gene") %>%
    dplyr::filter(rowSums(. > 0) >= floor(ncol(.)*fraction)) %>%
    dplyr::filter(rowSums(. == 0) >= floor(ncol(.)*fraction)) %>%
    tibble::rownames_to_column(var="gene")
  
}


rare_variant_to_long = function(tbl_list) { 
  # returns a tibble
  for(i in names(tbl_list)){
    tbl_list[[i]] = tbl_list[[i]] %>%
      tidyr::pivot_longer(.,cols = !gene,names_to = "line",values_to = "rare_burden") %>%
      dplyr::mutate(rare_mutation_type = i)
  }
  return(do.call("rbind",tbl_list))
  
}

scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

center_this <- function(x){
  (x - mean(x, na.rm=TRUE)) 
}

## gene symbol to Emsemblr ID ############
symbol_to_ensembl <- function(x){
  ensembl = mapIds(org.Hs.eg.db, x, "ENSEMBL","SYMBOL")
  return(ensembl)
}

# Ensembl ID to symbol
ensembl_to_symbol <- function(x){
  symbol = mapIds(org.Hs.eg.db, x, "SYMBOL","ENSEMBL")
  return(symbol)
}

symbol_to_entrez <- function(x){
  symbol = mapIds(org.Hs.eg.db, x, "ENTREZID","SYMBOL")
  return(symbol)
}


extract_SKATO_results <- function(dlist) {
  dlist = dlist[[1]] ### careful here
  
  df = data.frame("p_val"= dlist$p.value, 
                  "gene_name" = dlist$gene_name,
                  "gene_id"=dlist$gene_id)
  
  if(!is.null(dlist$resampling_pval)){
  
  # calculating proportion of permuted pval < actual pval
  # a.k.a permuted p-val
  df$resampling_pval = dlist$resampling_pval
  }else{
    df$resampling_pval = NA
  }
  return(df)
}

## qqplot
gg_qqplot = function(pvals) {
  df = data.frame(
    observed = -log10(sort(pvals)),
    expected = -log10(ppoints( length(pvals)))
  )
  ggplot(df) +
    geom_point(aes(expected, observed), shape = 1, size = 3, col = "grey40") +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    theme_bw() + 
    xlab( expression(paste("Expected -log"[10], plain(P)))) +
    ylab(expression(paste("Observed -log"[10], plain(P))))
}

plot_selected_coloc_locuszoom = function(coloc_res,treatment,locus_name,eGene,genes,dbsnp,nudge_y = 0){
  
  
  eQTL = coloc_res[[treatment]][[locus_name]]$results %>%
    dplyr::select(snp,V.df1,z.df1, r.df1 ,lABF.df1, SNP.PP.H4) %>%
    dplyr::rename(variant_id = snp) %>%
    tidyr::separate_wider_delim(cols = variant_id,delim = "_",names = c("chrom","pos","other_allele","effect_allele"),cols_remove = FALSE) %>%
    dplyr::mutate(Position_GRCh38 = paste0("chr",chrom,":",pos),
                  chrom = as.numeric(chrom),
                  pos = as.numeric(pos),
                  gene_name = eGene,
                  pval_nominal = 10 ^(( pnorm(-abs(z.df1),log.p=TRUE) + log(2) ) / log(10))) %>% # transforming this bc coloc uses p-vals derived from z score if betas + se from eQTL / GWAS are available
    dplyr::left_join(genes) %>%
    dplyr::left_join(dbsnp, by = "variant_id") 
  
  top_coloc_var = eQTL %>%
    dplyr::filter(gene_name == eGene) %>% # filter to eGene, beware this removes variants without rsId which I hope are not the lead (make a warning for that)
    dplyr::slice_min(pval_nominal) %>% # if there are ties by pval, sort by rsid
    dplyr::slice_min(rsid,na_rm = FALSE) %>% # if there are NAs in rsid, but not in all rows, take minimum rsid (not natural order)
    dplyr::distinct(rsid,.keep_all = TRUE)
  
  eQTL_subset = eQTL %>%
    dplyr::filter(gene_name == eGene & !is.na(rsid)) %>% # filter to eGene, beware this removes variants without rsId which I hope are not the lead (make a warning for that)
    data.frame()
  
  GWAS = coloc_res[[treatment]][[locus_name]][["results"]] %>%
    dplyr::select(snp,V.df2,z.df2, r.df2 ,lABF.df2, SNP.PP.H4) %>%
    dplyr::rename(variant_id = snp) %>%
    tidyr::separate_wider_delim(cols = variant_id,delim = "_",names = c("chrom","pos","other_allele","effect_allele"),cols_remove = FALSE) %>%
    dplyr::mutate(Position_GRCh38 = paste0("chr",chrom,":",pos),
                  chrom = as.numeric(chrom),
                  pos = as.numeric(pos),
                  gene_name = eGene,
                  pval_nominal = 10 ^(( pnorm(-abs(z.df2),log.p=TRUE) + log(2) ) / log(10))) %>%
    dplyr::left_join(genes) %>%
    dplyr::left_join(dbsnp, by = "variant_id") 
  
  GWAS_subset = GWAS %>%
    dplyr::filter(gene_name == eGene & !is.na(rsid)) %>% # filter to eGene, beware this removes variants without rsId which I hope are not the lead (make a warning for that)
    data.frame()
  
 if (is.na(top_coloc_var$rsid)){
   message("Warning: top coloc var has no rsID - Using top GWAS var ...")
   # locus info GWAS
   loc = locuszoomr::locus(GWAS_subset, gene = eGene, flank = 2.5e5,
                           ens_db = ensDb_v106, p = "pval_nominal", labs = "rsid", pos = "pos")
   
   # locus info eQTL
   
   loc_eqtl = locuszoomr::locus(eQTL_subset, gene = eGene, flank = 2.5e5,
                                ens_db = ensDb_v106, 
                                p = "pval_nominal", labs = "rsid", pos = "pos",
                                index_snp = loc$index_snp)

   
 } else{
   
   # locus info eQTL
   
   loc_eqtl = locuszoomr::locus(eQTL_subset, gene = eGene, flank = 2.5e5,
                                ens_db = ensDb_v106, 
                                p = "pval_nominal", labs = "rsid", pos = "pos",
                                index_snp = top_coloc_var$rsid)
   # locus info GWAS
   loc = locuszoomr::locus(GWAS_subset, gene = eGene, flank = 2.5e5,
                           ens_db = ensDb_v106, p = "pval_nominal", labs = "rsid", pos = "pos",
                           index_snp = top_coloc_var$rsid) # specifying the top coloc var
   
  
}
  # ld data eQTL
  
  #loc_eqtl = locuszoomr::link_LD(loc_eqtl, token = Sys.getenv("LDLINK_TOKEN"), pop="CEU",r2d = "r2",genome_build = "grch38_high_coverage")
  
  # getting LD through parent function and not locuszoomr version because it throws esoteric errors
  lddata = LDlinkR::LDmatrix(loc_eqtl$data$rsid,token = Sys.getenv("LDLINK_TOKEN"),
                             pop="CEU",
                             r2d = "r2",
                             genome_build = "grch38_high_coverage")
  labs <- loc_eqtl$labs
  index_snp <- loc_eqtl$index_snp
  
  ld <- lddata[, index_snp]
  loc_eqtl$data$ld <- ld[match(loc_eqtl$data[, labs], lddata$RS_number)]
  
  if(all(is.na( loc_eqtl$data$ld))){
    message("Warning: r2 of top coloc variant with other variants is all NA. Using GWAS index variant instead... ")
    
    # locus info GWAS
    loc = locuszoomr::locus(GWAS_subset, gene = eGene, flank = 2.5e5,
                            ens_db = ensDb_v106, p = "pval_nominal", labs = "rsid", pos = "pos") 
    # index snp will be min pval
    
    # locus info eQTL
    loc_eqtl = locuszoomr::locus(eQTL_subset, gene = eGene, flank = 2.5e5,
                                 ens_db = ensDb_v106, 
                                 p = "pval_nominal", labs = "rsid", pos = "pos",
                                 index_snp = loc$index_snp)
    
    
    lddata = LDlinkR::LDmatrix(loc_eqtl$data$rsid,token = Sys.getenv("LDLINK_TOKEN"),
                               pop="CEU",
                               r2d = "d",
                               genome_build = "grch38_high_coverage")
    labs <- loc_eqtl$labs
    index_snp <- loc_eqtl$index_snp

    ld <- lddata[, index_snp]
    loc_eqtl$data$ld <- ld[match(loc_eqtl$data[, labs], lddata$RS_number)]

# GWAS ld
    lddata = LDlinkR::LDmatrix(loc$data$rsid,token = Sys.getenv("LDLINK_TOKEN"),
                               pop="CEU",
                               r2d = "r2",
                               genome_build = "grch38_high_coverage")
    
    labs <- loc$labs
    index_snp <- loc$index_snp
    
    ld <- lddata[, index_snp]
    loc$data$ld <- ld[match(loc$data[, labs], lddata$RS_number)]
    
  }else{
    # loc = locuszoomr::link_LD(loc, token = Sys.getenv("LDLINK_TOKEN"), pop="CEU",
    #                           r2d = "r2",genome_build = "grch38_high_coverage")
    
    # getting LD through parent function and not locuszoomr version because it throws esoteric errors
    lddata = LDlinkR::LDmatrix(loc$data$rsid,token = Sys.getenv("LDLINK_TOKEN"),
                               pop="CEU",
                               r2d = "r2",
                               genome_build = "grch38_high_coverage")
    
    labs <- loc$labs
    index_snp <- loc$index_snp
    
    ld <- lddata[, index_snp]
    loc$data$ld <- ld[match(loc$data[, labs], lddata$RS_number)]
    
  }
  
  p_eqtl = locuszoomr::gg_scatter(loc_eqtl,labels = "index", 
                                  nudge_y = nudge_y) + ggtitle("eQTL")
  
  
  p_gwas = locuszoomr::gg_scatter(loc,
                                  pcutoff = 5e-100, # setting to NULL doesn't work so I give high limit
                                  legend_pos = NULL,
                                  labels = "index", nudge_y=nudge_y) + ggtitle("GWAS")
  g = locuszoomr::gg_genetracks(loc,gene_col = "grey40",exon_col = "grey40",exon_border=NA,
                                filter_gene_biotype = c("protein_coding") )
  p = (p_gwas / p_eqtl / g) +
    patchwork::plot_annotation(title = paste0(locus_name, " locus , ", eGene, " eGene - chr",
                                                                        top_coloc_var$chrom, " ",top_coloc_var$pos, " ",
                                                                        top_coloc_var$other_allele," \u2192 ", # unicode arrow
                                                                        top_coloc_var$effect_allele ),
                                                         subtitle =  paste0("PP for colocalization = ", 
                                                                            floor(coloc_res[[treatment]][[locus_name]]$summary[[6]] * 100), "%")) + 
    patchwork::plot_layout(ncol = 1,nrow = 3,heights = c(2,2,0.5))
  # if top coloc var info at top of plot doesn't match index snp, it's
  # because we had to take r2 from the lead GWAS snp due to LD issues
  
  return(list("p" = p, "eQTL"= eQTL,"GWAS"=GWAS))    
  
}

