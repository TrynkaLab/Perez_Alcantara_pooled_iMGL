# enrichment in interesting genes/loci
.libPaths(c('/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/',"/software/teamtrynka/ma23/R4.1/libs",
            "/software/R-4.1.0/lib/R/library"))
library(tidyverse)

outDir = "../../data/results/7.eQTL_results_gene_enrichment"
dir.create(outDir)


jeremy_candidate_genes = read.table("../../../resources/Jeremy_medrXiv_AD_candidate_genes.txt", header = TRUE)

set3 = read.csv("../../../resources/AD_PD_gene_sets_Andrew/Set3_manually_curated.csv")
  

sets=list()
sets[["AD_candidates"]] = c(jeremy_candidate_genes$symbol)

sets[["PD_candidates"]] = set3 %>%
  dplyr::filter(grepl("PD", disease)) %>%
  .$gene_name

sets[["ALS_candidates"]] = set3 %>%
  dplyr::filter(grepl("ALS", disease)) %>%
  .$gene_name

# tensorQTL file - all results

tensorqtl_all = readr::read_csv("../../data/results/4.Inspect_eQTL_results/tensorQTL_variant_gene_60PCs.csv")

# Hypergeometric test

# hypergeometric test with permutations

# Total balls in the urn: background dataset (all RNA-seq genes)
# x --> number of white balls drawn without replacement from an urn which contains both black and white balls (GWAS genes in my tested dataset)
# m --> White balls in the urn: characteristic I'm testing for enrichment (GWAS genes)
# n --> Black balls in the urn: total-white balls (non-GWAS genes)
# k --> Number of balls drawn from the urn: depends on the size of the dataset tested 

# I want the p-value of getting (number of GWAS genes in my DE dataset) or more white balls in a sample 
# of size (number of genes in my DE dataset) from an urn with (number of GWAS genes) white balls and 
# (total genes RNAseq-GWAS genes in that set) black balls.

# 1-phyper(q, m, n, k)

# q in this case is x-1 (diff exp dataset -1), the probability of x or more. It is a quantile defined as the smallest value x such that F(x) ≥ p, 
# where F is the distribution function.
# 


# comparison of probabilities

phyperRandom <- function(myGeneList, myGeneSet, genome){
  myRandomGS <- sample( genome,size=length(myGeneSet) )  # take random sample of genes from whole genome, same size as set tested
  myX <- length(which(myGeneList %in% myRandomGS))  # number of GWAS in random set
  myM <- length(which(myGeneList %in% genome)) #  GWAS genes in my pool
  myN <- length(genome) - length(which(myGeneList %in% genome))  # genes in genome not GWAS 
  myK <- length(myGeneSet) # size of set tested
  return(list(prob=1-phyper(q=myX-1,
                            m=myM, n=myN, k=myK), n.overlap=myX))
}

df_list <- list()
overlap_list <- list()  # to write down number of GWAS genes in diff expr results


for(set in names(sets)){
  
  prob_results=list()
  sums_probs=list()
  b_test <- list()
  overlap <- list()
  
  full_tensor_res = tensorqtl_all %>%
    group_split(group) %>%
    set_names(unique(tensorqtl_all$group)) %>%
    map(~ pull(., gene_name))
  sign_tensor_res = tensorqtl_all %>%
    dplyr::filter(pval_beta<0.05) %>%
    group_split(group) %>%
    set_names(unique(tensorqtl_all$group)) %>%
    map(~ pull(., gene_name))
  
  for(comparison in names(sign_tensor_res)){
    
    # p-value for enrichment in GWAS genes using my differentially expressed genes for each stage
    prob_results[[comparison]]=1-phyper(q=length(which(sets[[set]] %in% sign_tensor_res[[comparison]]))-1,
                                        m=length(which(sets[[set]] %in% full_tensor_res[[comparison]])),
                                        n=length(full_tensor_res[[comparison]]) - length(which(sets[[set]] %in% full_tensor_res[[comparison]])),
                                        k=length(sign_tensor_res[[comparison]]))
    gene_list=sets[[set]][which(sets[[set]] %in% sign_tensor_res[[comparison]])]
    
    overlap[[paste0("total_in_background:",comparison)]] = length(which(sets[[set]] %in% full_tensor_res[[comparison]]))
    
    overlap[[comparison]] = length(which(sets[[set]] %in% sign_tensor_res[[comparison]]))
    names(overlap[[comparison]]) = "GWAS_genes_in_enrichment_universe"
    # same, for random set of genes of same size. random gene set distribution of p-values
    random=list()
    random_overlaps=numeric()
    probab_random=numeric()
    
    for(i in 1:10000){
      random[[i]] <- as.data.frame(phyperRandom( sets[[set]], sign_tensor_res[[comparison]], full_tensor_res[[comparison]]))  # get probs from random gene sets, for each diff expr set tested
      probab_random[i] <- random[[i]]$prob
      random_overlaps[i] <- random[[i]]$n.overlap
      
    }
    
    
    # empirical p-value 
    
    sums_probs[[comparison]]$pval=((sum(probab_random<=prob_results[[comparison]])+1)/(10000+1)) # fraction of probabilities from random draws that are more extreme (smaller) than the one from my tested set
    # +1 in case sum(probab_random<=prob_results[[comparison]]) ==0
    # that would be the permuted? p-value
    
    # calculate the 95% confidence interval from the binomial distribution
    b_test[[comparison]]=binom.test(sum(probab_random<=prob_results[[comparison]]),10000)  # If those fractions are my "successes", then the binomial test can give me a confidence interval
    sums_probs[[comparison]]$conf.int.low=b_test[[comparison]]$conf.int[1]
    sums_probs[[comparison]]$conf.int.up=b_test[[comparison]]$conf.int[2]
    
    
    print(set)
    
  }
  overlap_list[[set]] =  do.call(rbind, overlap)
  
  df_list[[set]] = do.call(rbind, sums_probs)
  
}

overlap_list=  as.data.frame(do.call(rbind, Map(cbind, set = names(overlap_list), overlap_list)))
overlap_list$gene_set = rownames(overlap_list)
df_list= as.data.frame(do.call(rbind, Map(cbind, set = names(df_list), df_list))) 
df_list$comparison = rownames(df_list)

write.table(overlap_list,paste0(outDir,"/hypergeometric_enrichment_test_GWAS_n_genes.txt"), 
            col.names = T, row.names = F, quote = F)
df_list = apply(df_list,2,as.character)
write.table(df_list,paste0(outDir,"/hypergeometric_enrichment_test_GWAS_results.txt"), 
            col.names = T, row.names = F, quote = F)

# write genes that overlap

gene_list = list()
for(n in names(sign_tensor_res)){
  gene_list[[n]]$comparison = n
  gene_list[[n]]$genes_AD = paste0(sign_tensor_res[[n]][sign_tensor_res[[n]] %in% sets$AD_candidates],collapse = ",")
  gene_list[[n]]$genes_PD = paste0(sign_tensor_res[[n]][sign_tensor_res[[n]] %in% sets$PD_candidates],collapse = ",")
  gene_list[[n]]$genes_ALS = paste0(sign_tensor_res[[n]][sign_tensor_res[[n]] %in% sets$ALS_candidates],collapse = ",")
  
}
gene_list = as.data.frame(do.call(rbind, gene_list))

# fix this
gene_list = as.data.frame(df_list) %>%
  left_join(as.data.frame(gene_list))
gene_list = apply(gene_list,2,as.character)

write.table(gene_list,"../../data/results/2.1.Diff_expr_gene_check/hypergeometric_enrichment_test_GWAS_genes_in_sign_diff_expr_comparisons.txt", 
            col.names = T, row.names = F, quote = F)
