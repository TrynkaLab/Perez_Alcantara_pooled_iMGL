# interaction QTLs with LMM
library(lme4)
library(lmerTest)
library(interactions)
library(tidyverse)
library(variancePartition)
library(ggrepel)
source("./functions.R")
options(future.globals.maxSize = (250000 * 1024 ^ 2)) #(~250 Gb)

# try first with significant gene-variant pairs from chromosome 1
## LMM - lm() or something more resistant to outliers
# interaction plots with partial residuals
output_dir = "../../data/results/3.run_LMM_pheno_interactions"
dir.create(output_dir, recursive = TRUE)

phenotype = "migration"
treat = "untreated"
prolif = "Not_proliferating"
donor_blacklist="letw_5,lizq_3,zaie_1,romx_2,seru_7,qonc_2,iukl_1,curn_3,boqx_2,garx_2,sojd_3,yoch_6"
pcnumber = 15

donor_blacklist = stringr::str_split(donor_blacklist,pattern = ",")[[1]]

####

for(treat in c("untreated", "IFN", "LPS")){
# for(treat in c( "untreated", "IFN")){
# load files

# interaction pheno data
mean_pheno = readr::read_csv(paste0("../../../OTAR2065_phenotypic_QTLs/data/results/",
                                    phenotype,
                                    "/1.check_line_proportions/line_prop_changes_per_well.csv")) %>%
  dplyr::filter(!log_fraction_mean %in% c(NA,Inf, -Inf)) %>%
  dplyr::filter(treatment  ==  treat) %>%
  dplyr::group_by(line, pool) %>%
  dplyr::summarise(scaled_phenotype = mean(scaled_log_fraction)) %>%
  # dplyr::rename(scaled_phenotype = scaled_log_fraction) %>%
  dplyr::mutate(line_pool = paste(line,pool,sep = "_"))

# covariates
metadata = readr::read_tsv("../../../OTAR2065_differentiation_efficiency/data/results/5.prepare_differential_expression_pseudobulks/metadata_treatmentxprolifxidxpool.txt") %>%
  dplyr::mutate(line_pool = paste(donor_id,pool,sep = "_")) %>%
  dplyr::filter(treatment == treat & proliferation_status == prolif & line_pool %in% mean_pheno$line_pool) %>%
  dplyr::rename(line=donor_id)
# scaled_phenotype
pseudobulk = readr::read_tsv("../../../OTAR2065_differentiation_efficiency/data/results/5.prepare_differential_expression_pseudobulks/combined_pseudobulk_raw_counts_treatmentxprolifxidxpool.txt") %>%
  dplyr::select(c("gene",metadata$cols_names)) %>%
  column_to_rownames("gene")

## kinship 
kinship = read_tsv("../../data/kinship/kinship_for_LMM.tsv")
# from ana cuomo: Gower normalisation                    kinship_mat *= (kinship_mat.shape[0] - 1) / (kinship_mat.trace() - kinship_mat.mean(0).sum())
# then to vector with economic_qs
# from limix https://github.com/limix/numpy-sugar/blob/main/numpy_sugar/linalg/qs.py
# kinship_vectors = as.data.frame(eigen(kinship, symmetric = TRUE)$vectors)
kinship_pc = as.data.frame(prcomp(kinship)$rotation)
rownames(kinship_pc) = colnames(kinship)
colnames(kinship_pc) = paste0("kinshipPC",1:ncol(kinship_pc))
# expand to all replicates
kinship_pc = kinship_pc %>%
  tibble::rownames_to_column(var = "line") %>%
  dplyr::select(c("line",paste0("kinshipPC",1:5)))

p1 = ggplot(kinship_pc,aes(x = kinshipPC1, y = kinshipPC2, label = line)) +
  geom_point() + 
  geom_text_repel() + 
  theme_bw()

metadata = metadata %>%
  dplyr::left_join(kinship_pc)
## genotype PCs

full_genotype_pcs = read.table("../../data/genotype/plink_genotypes/all_pools.genotype.MAF05.eigenvec")
rownames(full_genotype_pcs) = full_genotype_pcs$V2 # should be line names, or donor names if line was not available


full_genotype_pcs = full_genotype_pcs[c(-1,-2)]
colnames(full_genotype_pcs) = paste0("genotypePC",1:ncol(full_genotype_pcs))
# 
# p2 = full_genotype_pcs %>%
#   tibble::rownames_to_column("line") %>%
#   ggplot(.,aes(x = genotypePC1, y = genotypePC2, label = line)) +
#   geom_point() + 
#   geom_text_repel() + 
#   theme_bw()
# 
# pdf(file = paste0(output_dir,"/kinship_PC_genotype_PC.pdf"),
#     width = 6, height = 5)
# p1 + p2
# dev.off()

genotype_pcs = full_genotype_pcs %>%
  tibble::rownames_to_column(var = "line") %>%
  dplyr::select(c("line",paste0("genotypePC",1:5)))

metadata = metadata %>%
  dplyr::left_join(genotype_pcs)
# join scaled_phenotype metadata and phenotype
metadata = metadata %>%
  dplyr::left_join(mean_pheno[,c("line_pool","scaled_phenotype")],by = "line_pool")

## filters as in regular tensorQTL

message("Filtering by number of cells: remove donors with fewer than 100 cells")
retain = metadata$count >=100
metadata = metadata[retain,]
pseudobulk = pseudobulk[,metadata$cols_names]

# transform counts
metadata = metadata %>%
  dplyr::mutate(one_over_ncells = 1/count)


message("Excluding donors from the blacklist")

retain = !(metadata$line  %in% donor_blacklist)
metadata = metadata[retain,]
pseudobulk = pseudobulk[,metadata$cols_names]

# rename columns of pseudobulk by line_pool_rep
colnames(pseudobulk) = metadata$line_pool


#log2 the summed counts applying size factors
message("Normalising to lib size, calculating log2(counts + 1)")

size_factors = colSums(pseudobulk)/mean(colSums(pseudobulk))
log_sum_counts = as.data.frame(log2(t(t(pseudobulk)/size_factors) + 1))
cpm =  apply(pseudobulk,2, function(x) (x/sum(x))*1000000)

# retain genes whose mean counts across donors are > 1CPM
# and retain genes that are present in at least 30% of the donors with at least 1CPM (rounding down)
filtered = cpm %>%
  as.data.frame() %>%
  dplyr::filter(rowMeans(.) >= 1 & (floor((rowSums(. > 0)/ncol(.))*100))>=30) 
#anyDuplicated(colnames(cpm))
#cpm[,duplicated(colnames(cpm))]
log_sum_counts = log_sum_counts[rownames(filtered),]
message("There are ",nrow(log_sum_counts)," genes left after filtering")


message("Mean-SD plots")

# mean-sd plots

  p1 =  vsn::meanSdPlot(as.matrix(log_sum_counts), plot = FALSE)$gg +theme_bw() +
    ggtitle("Libsize-scaled + log2 sum of raw counts (ranks)")
  
  pdf(file = paste0(output_dir,"/mean_sd_sum_pseudobulk_noreps.pdf"),
      width = 6, height = 5)
  plot(p1)
  dev.off()

  message("Annotating ensembl ids and positions")
  
  ensembl = read_csv("../../../resources/ENSEMBL_human_v111_2023/Homo_sapiens.GRCh38.111.genes.csv") %>%
    dplyr::select(seqname,start,end,gene_id,gene_name)
  
  log_sum_counts = log_sum_counts %>%
      tibble::rownames_to_column(var = "gene_name") %>%
      dplyr::left_join(.,ensembl,by="gene_name") %>%
      dplyr::filter(!is.na(gene_id)) %>% # removing genes that are not found, around 900 from 12.5k
      dplyr::filter(seqname %in% c(1:22)) %>%  # Remove chromosomes not in 1:22
      dplyr::relocate(seqname,start,end,gene_id,gene_name) %>%
      dplyr::rename(chr = seqname)
  message("There are ",nrow(log_sum_counts)," genes left after removing unmatched ensembl ids")
  
  
    fixed_cols = log_sum_counts %>%
      dplyr::select("chr","start","end","gene_id","gene_name")

    numeric_df = log_sum_counts %>%
      dplyr::select(!c(chr,start,end,gene_id, gene_name)) %>%
      t() %>% # transposing because scales per column (and we want to scale per gene)
      as.data.frame() %>%
      mutate_all(as.numeric)
    scaled = base::scale(numeric_df) %>%
      t() %>%
      as.data.frame() %>%
      mutate_all(as.numeric) 
    
    
    scaled = cbind(fixed_cols,scaled)
    colnames(scaled) = c("chr","start","end","gene_id","gene_name",colnames(pseudobulk))
    
    
    ####### calculating PCS #####
    message("Calculating expression PCs")
    
    numeric_df = scaled %>%
      dplyr::select(!c(chr,start,end,gene_id,gene_name)) %>%
      as.data.frame() %>%
      mutate_all(as.numeric)
    
    pcs=prcomp(numeric_df,center = FALSE) #already centered

    if(pcnumber > ncol(pcs$rotation)){
      warning("pcnumber is larger than the number of columns in pcs$rotation.")
      
      break
    }else{
      pcs = pcs$rotation[,1:pcnumber] %>%
        as.data.frame() %>%
        rownames_to_column(var = "line_pool")
      metadata = metadata %>%
        dplyr::left_join(pcs, by = "line_pool")
    }
    
    

    
message("Read in genotype file with significant eQTL results")

significant_tensor_eQTL = read_csv("../../data/results/4.Inspect_eQTL_results/tensorQTL_variant_gene_60PCs.csv") %>%
  dplyr::filter(treatment == treat & cluster == prolif) %>%
  dplyr::mutate(qval = qvalue::qvalue(pval_beta)$qvalues) %>% # qvalue should already be there?
  dplyr::filter(qval < 0.01) 

message("Also retain interesting genes from Sam's phagocytosis screen, and AD and PD candidate genes")

sam_genes = read_tsv("../../../CRISPR/OTAR2065_phagocytosis_CRISPR/data/2024_07_sam_screen_results.tsv") %>%
  pull(id)

AD_candidate_genes = read_csv("../../../resources/AD_PD_gene_sets_Andrew/Set1_Jeremy_candidates.csv") %>%
  dplyr::pull(gene_sym)

PD_candidate_genes = read_csv("../../../resources/AD_PD_gene_sets_Andrew/Set3_manually_curated.csv") %>%
  dplyr::filter(str_detect(disease,"PD")) %>%
  dplyr::pull(gene_name)

all_genes = unique(c(sam_genes,AD_candidate_genes, PD_candidate_genes)) # 167 additional genes

more_tensor_eQTL = read_csv("../../data/results/4.Inspect_eQTL_results/tensorQTL_variant_gene_60PCs.csv") %>%
  dplyr::filter(treatment == treat & cluster == prolif) %>%
  dplyr::mutate(qval = qvalue::qvalue(pval_beta)$qvalues) %>%
  dplyr::filter(gene_name %in% all_genes) # 133 found

significant_tensor_eQTL = significant_tensor_eQTL %>%
  dplyr::rows_append(more_tensor_eQTL) %>%
  distinct() # 2059 rows

nrow(significant_tensor_eQTL)
message("Read in chunked genotype with allele dosages")
genotype = list()
for(n in 1:48){
  print(n)
  genotype[[n]] = readr::read_csv(paste0("../../../OTAR2065_phenotypic_QTLs/data/genotypes/full_genotype/genotype_minor_allele_dosage_",
                                         n,
                                         ".csv")) %>%
    dplyr::relocate(rn) %>%
    dplyr::rename(variant_id = rn) %>%
    dplyr::mutate(variant_id = str_replace(variant_id,pattern="chr",replacement = "")) %>%
    dplyr::filter(variant_id %in% significant_tensor_eQTL$variant_id)   # filter to shared
  
}
genotype = do.call("rbind",genotype)
saveRDS(genotype,paste0(output_dir,"/pheno_genotype_",treat,"_",prolif,"_",phenotype,".rds")) # to save a bit of time
genotype = readRDS(paste0(output_dir,"/pheno_genotype_",treat,"_",prolif,"_",phenotype,".rds"))
nrow(genotype)
nrow(significant_tensor_eQTL) # why are some missing?

## variance partition #######
# var ="2_119355944_A_G"
# test = genotype %>%
#   dplyr::filter(variant_id == var) %>%
#   tidyr::pivot_longer(!variant_id,names_to = "line", values_to = "genotype") %>%
#   dplyr::mutate(genotype = case_when(genotype == 0.5 ~ 1,
#                                      genotype == 1 ~ 2,
#                                      .default = genotype))
# # add metadata
# test = test %>%
#   dplyr::left_join(metadata[,c("line","treatment","pool","line_pool",
#                                paste0("genotypePC",1:5),paste0("kinshipPC",1:5),
#                                paste0("PC",1:pcnumber),
#                                "scaled_phenotype","one_over_ncells")]) %>%
#   dplyr::filter(!is.na(treatment))
# 
# # add scaled_phenotype###
# # subset to significant gene from tested variant
# sign_gene = significant_tensor_eQTL %>%
#   dplyr::filter(variant_id == var)
# gen = unique(sign_gene$gene_name)
# 
# subset_scaled = scaled %>%
#   dplyr::filter(gene_name == gen) %>%
#   tidyr::pivot_longer(!c("chr","start","end","gene_id","gene_name"),
#                       names_to = "line_pool",values_to = "scaled_expression")
# 
# test = test %>%
#   dplyr::left_join(subset_scaled)
# 
# form = scaled_expression ~  genotypePC1 + genotypePC2 + genotypePC3 + genotypePC4 + genotypePC5 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 +
#   genotype*scaled_phenotype + one_over_ncells + line +  (1|pool) + (1|kinshipPC1) + (1|kinshipPC2) + (1|kinshipPC3) + (1|kinshipPC4) + (1|kinshipPC5)
# 
# # canonical correlation analysis
# canon =  canCorPairs(form, test)
# # Plot correlation matrix
# # between all pairs of variables
# pdf(paste0(output_dir,"/canon_correlation_LMM_fit_DBI_2_119355944_A_G",treat,"_",phenotype,".pdf"), width = 6, height = 6)
# plotCorrMatrix(canon,sort = FALSE)
# dev.off()

# genotype PCs are correlated
# one over ncells and pool are correlated with scaled_phenotype PCs
# are genotype PCs uncorrelated when there are no repeated donors??
form = ~ genotypePC1 + genotypePC2 + genotypePC3 + genotypePC4 + genotypePC5  + (1|kinshipPC1) + (1|kinshipPC2) + (1|kinshipPC3) + (1|kinshipPC4) + (1|kinshipPC5)

rownames(kinship_pc) = kinship_pc$line
kinship_pc = kinship_pc[-1]

canon =  canCorPairs(form, merge(full_genotype_pcs,kinship_pc))
# Plot correlation matrix
# between all pairs of variables
pdf(paste0(output_dir,"/canon_correlation_orig_genotype_PCs.pdf"), width = 6, height = 6)
plotCorrMatrix(canon,sort = FALSE)
dev.off()

form = ~ kinshipPC1 + kinshipPC2 + kinshipPC3 + kinshipPC4 + kinshipPC5

canon =  canCorPairs(form, kinship_pc)
# Plot correlation matrix
# between all pairs of variables
pdf(paste0(output_dir,"/canon_correlation_orig_kinship_pc.pdf"), width = 6, height = 6)
plotCorrMatrix(canon,sort = FALSE)
dev.off()
# they are slightly correlated (not PCs)
####### LMM #########

my_fit = list()
data = list()
for(var in genotype$variant_id){
  
  # subset and fix genotype
  test = genotype %>%
    dplyr::filter(variant_id == var) %>%
    tidyr::pivot_longer(!variant_id,names_to = "line", values_to = "genotype") %>%
    dplyr::mutate(genotype = case_when(genotype == 0.5 ~ 1,
                                       genotype == 1 ~ 2,
                                       .default = genotype))
  # add metadata
  test = test %>%
    dplyr::left_join(metadata[,c("line","treatment","pool","line_pool",
                                 paste0("genotypePC",1:5),paste0("kinshipPC",1:5),
                                 paste0("PC",1:pcnumber),
                                 "scaled_phenotype","one_over_ncells")]) %>%
    dplyr::filter(!is.na(treatment))
    
  # add scaled_phenotype###
  # subset to significant gene from tested variant
  sign_gene = significant_tensor_eQTL %>%
    dplyr::filter(variant_id == var)
  
  for(gen in unique(sign_gene$gene_name)){
  
    
    subset_scaled =  scaled %>%
      dplyr::filter(gene_name == gen) %>%
      tidyr::pivot_longer(!c("chr","start","end","gene_id","gene_name"),
                          names_to = "line_pool",values_to = "scaled_expression")
    

    
    data[[paste(gen,var,sep = "-")]] = test %>%
      dplyr::left_join(subset_scaled) 
    
    # if any category < 5 counts, remove
    # then if there's only one category, omit
    category_to_omit = names(table(data[[paste(gen,var,sep = "-")]]$genotype)[table(data[[paste(gen,var,sep = "-")]]$genotype)<5])
    data[[paste(gen,var,sep = "-")]] = data[[paste(gen,var,sep = "-")]] %>%
      dplyr::filter(!genotype %in% category_to_omit)
    if(length(unique(data[[paste(gen,var,sep = "-")]]$genotype))>1){
      # fit
      if(nrow(subset_scaled)==0){
        message("Not found in scaled gene scaled_phenotype - old eQTL results?")
        next()
        
      }else{
        my_fit[[paste(gen,var,sep = "-")]] = lmerTest::lmer(scaled_expression ~  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 +
                                                              genotype*scaled_phenotype + 
                                                              (1 | kinshipPC1)  + (1 | kinshipPC2)  + (1 | kinshipPC3)  + (1 | kinshipPC4)  + (1 | kinshipPC5),
                                                            data = data[[paste(gen,var,sep = "-")]],
                                                            control = lmerControl(optCtrl=list(maxfun=5000) ),
                                                            REML = FALSE)
       
      }
    }else{
      message("Skipping because of insufficient genotype categories")
      
      next()
    }
    

    
  }
 
}


saveRDS(my_fit,paste0(output_dir,"/",phenotype,"_",treat,"_",prolif,"_LMM_noreps.rds")) 
saveRDS(data,paste0(output_dir,"/",phenotype,"_",treat,"_",prolif,"_data_noreps.rds"))

my_fit = readRDS(paste0(output_dir,"/",phenotype,"_",treat,"_",prolif,"_LMM_noreps.rds"))
data = readRDS(paste0(output_dir,"/",phenotype,"_",treat,"_",prolif,"_data_noreps.rds"))
to_sort = lapply(my_fit,jtools::summ)
pvals = unlist(lapply(to_sort, function(x) x$coeftable["genotype:scaled_phenotype", "p"]))
beta = unlist(lapply(to_sort, function(x) x$coeftable["genotype:scaled_phenotype", "Est."]))

qvals = qvalue::qvalue(pvals, pi0.method = "bootstrap")$qvalues
qvals = sort(qvals)

examine_qvals = qvalue::qvalue(pvals, pi0.method = "bootstrap")

toplot = data.frame(pi0.lambda = examine_qvals$pi0.lambda,
                    lambda = examine_qvals$lambda)

# plot(examine_qvals)


# pdf(paste0(output_dir,"/",treat,"_",phenotype,"qvals_QC_norep.pdf"), height = 7, width = 7)
# p1 = ggplot(toplot,aes(y=pi0.lambda,lambda)) + 
#   geom_point() + 
#   theme_minimal() + 
#   geom_hline(yintercept = examine_qvals$pi0, linetype = "dashed") + 
#   ggtitle(paste0("pi0 =" ,round(examine_qvals$pi0,3)) )
# 
# p2 = hist(examine_qvals)
# p3 = qvals %>%
#   as_tibble() %>%
#   ggplot(.,aes(value)) + geom_histogram() + theme_minimal() + 
#   ggtitle("qvalues")
# 
# (p1 + p2) / p3
# 
# dev.off()

# 
# pdf(paste0(output_dir,"/",treat,"_",phenotype,"pvals_hist_norep.pdf"), height = 8, width = 8)
# p1 = pvals %>%
#   as_tibble() %>%
#   ggplot(.,aes(value)) + geom_histogram() + theme_minimal() + 
#   ggtitle("pvalues")
# p2 = qvals %>%
#   as_tibble() %>%
#   ggplot(.,aes(value)) + geom_histogram() + theme_minimal() + 
#   ggtitle("qvalues")
# 
# p1 / p2
# 
# dev.off()



hist(pvals)
hist(qvals)
table(pvals<0.10) # 422 untreated phagocytosis, 312 LPS phagocytosis
table(pvals<0.05) # 264 untreated phagocytosis, 15 untreated migration, 5 IFN phagocytosis, 6 IFN migration, 171 phagocytosis LPS, 8 migration LPS
table(qvals<0.20) # 117 untreated phagocytosis, 12 untreated migration, 5 IFN phagocytosis, 1 IFN migration, 30 phagocytosis LPS,2 migration LPS
table(qvals<0.10) # 39 untreated phagocytosis, 4 LPS phagocytosis
table(qvals<0.05) # 16 untreated phagocytosis, 3 LPS phagocytosis



# plotting all significant results  
pdf(paste0(output_dir,"/all_q0.2_",phenotype,"_",treat,"_",prolif,"_interactions_expr_x_pheno_noreps.pdf"),
    width = 7,height = 7
           )
qvasl = sort(qvals)
  sign = names(qvals)[qvals<0.20]
for(gene_var in sign){
  
 
  genotype_vals = unique(data[[gene_var]]$genotype)
  
  # with partial residuals
  
  p1 = interact_plot(my_fit[[gene_var]], 
                     pred = scaled_phenotype, 
                     modx = genotype, plot.points = TRUE,
                     modx.values = genotype_vals,
                     partial.residuals = TRUE,
                     interval = FALSE,
                     colors = "Qual1") + 
    ggtitle(paste0(gene_var) ) +
    ylab("Gene expression") + 
    xlab("Phenotype")
  # with partial residuals
  
 
  plot(p1)
}
dev.off()

}
######################
#####################
#################### 

### now with replicates ###########

mean_pheno = readr::read_csv(paste0("../../../OTAR2065_phenotypic_QTLs/data/results/",
                                    phenotype,
                                    "/1.check_line_proportions/line_prop_changes_per_well.csv")) %>%
  dplyr::filter(!log_fraction_mean %in% c(NA,Inf, -Inf)) %>%
  dplyr::filter(treatment  ==  treat) %>%
  # dplyr::group_by(line, pool) %>%
  # dplyr::summarise(scaled_phenotype = mean(scaled_log_fraction)) %>%
  dplyr::rename(scaled_phenotype = scaled_log_fraction) %>%
  dplyr::mutate(line_pool = paste(line,pool,sep = "_"),
                line_pool_rep = paste(line,pool,replicate,sep = "_"),
                pool_rep = paste(pool,replicate,sep = "_"))
# covariates
metadata = readr::read_tsv("../../../OTAR2065_differentiation_efficiency/data/results/5.prepare_differential_expression_pseudobulks/metadata_treatmentxprolifxidxpool_reps.txt") %>%
  dplyr::mutate(line_pool = paste(donor_id,pool,sep = "_"),
                pool_rep = paste(pool,replicate,sep = "_"),
                line_pool_rep =  paste(donor_id,pool,replicate,sep = "_")) %>%
  dplyr::filter(treatment == treat & proliferation_status == prolif & line_pool_rep %in% mean_pheno$line_pool_rep) %>%
  dplyr::rename(line=donor_id)
# scaled_phenotype
pseudobulk = readr::read_tsv("../../../OTAR2065_differentiation_efficiency/data/results/5.prepare_differential_expression_pseudobulks/combined_pseudobulk_raw_counts_treatmentxprolifxidxpool_reps.txt") %>%
  dplyr::select(c("gene",metadata$cols_names)) %>%
  column_to_rownames("gene")


## genotype PCs

full_genotype_pcs = read.table("../../data/genotype/plink_genotypes/all_pools.genotype.MAF05.eigenvec")
rownames(full_genotype_pcs) = full_genotype_pcs$V2 # should be line names, or donor names if line was not available


full_genotype_pcs = full_genotype_pcs[c(-1,-2)]
colnames(full_genotype_pcs) = paste0("genotypePC",1:ncol(full_genotype_pcs))
genotype_pcs = full_genotype_pcs %>%
  tibble::rownames_to_column(var = "line") %>%
  dplyr::select(c("line",paste0("genotypePC",1:5)))

metadata = metadata %>%
  dplyr::left_join(genotype_pcs)
# join scaled_phenotype metadata and phenotype
metadata = metadata %>%
  dplyr::left_join(mean_pheno[,c("line_pool","line_pool_rep","pool_rep","well","scaled_phenotype")],by = "line_pool_rep")

## filters as in regular tensorQTL

message("Filtering by number of cells: remove donors with fewer than 50 cells")
retain = metadata$count >=50
metadata = metadata[retain,]
pseudobulk = pseudobulk[,metadata$cols_names]

# transform counts
metadata = metadata %>%
  dplyr::mutate(one_over_ncells = 1/count)


message("Excluding donors from the blacklist")

retain = !(metadata$line  %in% donor_blacklist)
metadata = metadata[retain,]
pseudobulk = pseudobulk[,metadata$cols_names]

# rename columns of pseudobulk by line_pool_rep
colnames(pseudobulk) = metadata$line_pool_rep


#log2 the summed counts applying size factors
message("Normalising to lib size, calculating log2(counts + 1)")

size_factors = colSums(pseudobulk)/mean(colSums(pseudobulk))
log_sum_counts = as.data.frame(log2(t(t(pseudobulk)/size_factors) + 1))
cpm =  apply(pseudobulk,2, function(x) (x/sum(x))*1000000)

# retain genes whose mean counts across donors are > 1CPM
# and retain genes that are present in at least 30% of the donors with at least 1CPM (rounding down)
filtered = cpm %>%
  as.data.frame() %>%
  dplyr::filter(rowMeans(.) >= 1 & (floor((rowSums(. > 0)/ncol(.))*100))>=30) 
#anyDuplicated(colnames(cpm))
#cpm[,duplicated(colnames(cpm))]
log_sum_counts = log_sum_counts[rownames(filtered),]
message("There are ",nrow(log_sum_counts)," genes left after filtering")


message("Mean-SD plots")

# mean-sd plots

p1 =  vsn::meanSdPlot(as.matrix(log_sum_counts), plot = FALSE)$gg +theme_bw() +
  ggtitle("Libsize-scaled + log2 sum of raw counts (ranks)")

pdf(file = paste0(output_dir,"/mean_sd_sum_pseudobulk.pdf"),
    width = 6, height = 5)
plot(p1)
dev.off()

message("Annotating ensembl ids and positions")

ensembl = read_csv("../../../resources/ENSEMBL_human_v111_2023/Homo_sapiens.GRCh38.111.genes.csv") %>%
  dplyr::select(seqname,start,end,gene_id,gene_name)

log_sum_counts = log_sum_counts %>%
  tibble::rownames_to_column(var = "gene_name") %>%
  dplyr::left_join(.,ensembl,by="gene_name") %>%
  dplyr::filter(!is.na(gene_id)) %>% # removing genes that are not found, around 900 from 12.5k
  dplyr::filter(seqname %in% c(1:22)) %>%  # Remove chromosomes not in 1:22
  dplyr::relocate(seqname,start,end,gene_id,gene_name) %>%
  dplyr::rename(chr = seqname)
message("There are ",nrow(log_sum_counts)," genes left after removing unmatched ensembl ids")


fixed_cols = log_sum_counts %>%
  dplyr::select("chr","start","end","gene_id","gene_name")

numeric_df = log_sum_counts %>%
  dplyr::select(!c(chr,start,end,gene_id, gene_name)) %>%
  t() %>% # transposing because scales per column (and we want to scale per gene)
  as.data.frame() %>%
  mutate_all(as.numeric)
scaled = base::scale(numeric_df) %>%
  t() %>%
  as.data.frame() %>%
  mutate_all(as.numeric) 


scaled = cbind(fixed_cols,scaled)
colnames(scaled) = c("chr","start","end","gene_id","gene_name",colnames(pseudobulk))


####### calculating PCS #####
message("Calculating expression PCs")

numeric_df = scaled %>%
  dplyr::select(!c(chr,start,end,gene_id,gene_name)) %>%
  as.data.frame() %>%
  mutate_all(as.numeric)

pcs=prcomp(numeric_df,center = FALSE) #already centered

if(pcnumber > ncol(pcs$rotation)){
  warning("pcnumber is larger than the number of columns in pcs$rotation.")
  
  break
}else{
  pcs = pcs$rotation[,1:pcnumber] %>%
    as.data.frame() %>%
    rownames_to_column(var = "line_pool_rep")
  metadata = metadata %>%
    dplyr::left_join(pcs, by = "line_pool_rep")
}




message("Read in genotype file with significant eQTL results")


significant_tensor_eQTL = read_csv("../../data/results/4.Inspect_eQTL_results/tensorQTL_variant_gene_60PCs.csv") %>%
  dplyr::filter(treatment == treat & cluster == prolif) %>%
  dplyr::mutate(qval = qvalue::qvalue(pval_beta)$qvalues) %>% # qvalue should already be there?
  dplyr::filter(qval < 0.01) 

message("Read in chunked genotype with allele dosages")
genotype = list()
for(n in 1:48){
  genotype[[n]] = readr::read_csv(paste0("../../../OTAR2065_phenotypic_QTLs/data/genotypes/full_genotype/genotype_minor_allele_dosage_",
                                         n,
                                         ".csv")) %>%
    dplyr::relocate(rn) %>%
    dplyr::rename(variant_id = rn) %>%
    dplyr::mutate(variant_id = str_replace(variant_id,pattern="chr",replacement = "")) %>%
    dplyr::filter(variant_id %in% significant_tensor_eQTL$variant_id)   # filter to shared
  
}
genotype = do.call("rbind",genotype)
saveRDS(genotype,paste0(output_dir,"/pheno_genotype_",treat,"_",prolif,"_",phenotype,".rds")) # to save a bit of time
genotype = readRDS(paste0(output_dir,"/pheno_genotype_",treat,"_",prolif,"_",phenotype,".rds"))

## variance partition #######
var ="2_119355944_A_G"
test = genotype %>%
  dplyr::filter(variant_id == var) %>%
  tidyr::pivot_longer(!variant_id,names_to = "line", values_to = "genotype") %>%
  dplyr::mutate(genotype = case_when(genotype == 0.5 ~ 1,
                                     genotype == 1 ~ 2,
                                     .default = genotype))
# add metadata
test = test %>%
  dplyr::left_join(metadata[,c("line","treatment","pool","line_pool","pool_rep","line_pool_rep",paste0("genotypePC",1:5),
                               paste0("PC",1:pcnumber),
                               "scaled_phenotype","one_over_ncells")]) %>%
  dplyr::filter(!is.na(treatment))

# add scaled_phenotype###
# subset to significant gene from tested variant
sign_gene = significant_tensor_eQTL %>%
  dplyr::filter(variant_id == var)
gen = unique(sign_gene$gene_name)

subset_scaled = scaled %>%
  dplyr::filter(gene_name == gen) %>%
  tidyr::pivot_longer(!c("chr","start","end","gene_id","gene_name"),
                      names_to = "line_pool_rep",values_to = "scaled_expression")

test = test %>%
  dplyr::left_join(subset_scaled)

form = scaled_expression ~  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 +
  genotype*scaled_phenotype + line + (1 | genotypePC1)  + (1 | genotypePC2)  + (1 | genotypePC3)  + (1 | genotypePC4)  + (1 | genotypePC5)

# canonical correlation analysis
canon =  canCorPairs(form, test)
# Plot correlation matrix
# between all pairs of variables
pdf(paste0(output_dir,"/canon_correlation_LMM_fit_DBI_2_119355944_A_G",treat,"_",phenotype,".pdf"), width = 6, height = 6)
plotCorrMatrix(canon,sort = FALSE)
dev.off()

# genotype PCs are correlated
# one over ncells and pool are correlated with scaled_phenotype PCs
# are genotype PCs uncorrelated when there are no repeated donors??
form = ~ genotypePC1 + genotypePC2 + genotypePC3 + genotypePC4 + genotypePC5

canon =  canCorPairs(form, full_genotype_pcs)
# Plot correlation matrix
# between all pairs of variables
pdf(paste0(output_dir,"/canon_correlation_orig_genotype_PCs.pdf"), width = 6, height = 6)
plotCorrMatrix(canon,sort = FALSE)
dev.off()
# they are uncorrelated
####### LMM #########

metadata = metadata %>%
  dplyr::rename(line_pool = line_pool.x,
                pool_rep = pool_rep.x) 
my_fit_reps = list()
test2_reps = list()
for(var in genotype$variant_id){
  
  # subset and fix genotype
  test = genotype %>%
    dplyr::filter(variant_id == var) %>%
    tidyr::pivot_longer(!variant_id,names_to = "line", values_to = "genotype") %>%
    dplyr::mutate(genotype = case_when(genotype == 0.5 ~ 1,
                                       genotype == 1 ~ 2,
                                       .default = genotype))
  # add metadata
  test = test %>%
    dplyr::left_join(metadata[,c("line","treatment","pool","line_pool","line_pool_rep","pool_rep",paste0("genotypePC",1:5),
                                 paste0("PC",1:pcnumber),
                                 "scaled_phenotype","one_over_ncells")]) %>%
    dplyr::filter(!is.na(treatment))
  
  # add scaled_phenotype###
  # subset to significant gene from tested variant
  sign_gene = significant_tensor_eQTL %>%
    dplyr::filter(variant_id == var)
  
  for(gen in unique(sign_gene$gene_name)){
    message("Working on gene", gen)
    subset_scaled = scaled %>%
      dplyr::filter(gene_name == gen) %>%
      tidyr::pivot_longer(!c("chr","start","end","gene_id","gene_name"),
                          names_to = "line_pool_rep",values_to = "scaled_expression")
    
    test2_reps[[paste(gen,var,sep = "-")]] = test %>%
      dplyr::left_join(subset_scaled)
    
    # fit
    # if any category < 5 lines contributing, remove
    # then if there's only one category, omit
    category_to_omit = test2_reps[[paste(gen,var,sep = "-")]] %>%
      dplyr::group_by(genotype) %>%
      dplyr::summarise(count = n_distinct(line)) %>%
      dplyr::mutate(to_omit = count<5) %>%
      dplyr::filter(to_omit==TRUE) %>%
      dplyr::select(genotype)
    test2_reps[[paste(gen,var,sep = "-")]] = test2_reps[[paste(gen,var,sep = "-")]] %>%
      dplyr::filter(!genotype %in% category_to_omit$genotype)
    
    if(length(unique(test2_reps[[paste(gen,var,sep = "-")]]$genotype))>1){
      
      if(nrow(subset_scaled)==0){
        message("Not found in scaled gene scaled_phenotype - old eQTL results?")
        next()
        
      }else{
        message("Fitting LMM for ",paste(gen,var,sep = "-") )
        
        my_fit_reps[[paste(gen,var,sep = "-")]] = lmerTest::lmer(scaled_expression ~   PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 +
                                                                   genotype*scaled_phenotype + (1 | genotypePC1)  + (1 | genotypePC2)  + (1 | genotypePC3)  + (1 | genotypePC4)  + (1 | genotypePC5) ,
                                                                 data = test2_reps[[paste(gen,var,sep = "-")]],
                                                                 control = lmerControl(optCtrl=list(maxfun=5000) ))
        
        
      }
    }else{
      message("Skipping because of insufficient genotype categories")
      
      next()
    }
    
  }
  
}

saveRDS(my_fit_reps,paste0(output_dir,"/LMM_reps.rds"))
to_sort_reps = lapply(my_fit_reps,jtools::summ)
pvals_reps = unlist(lapply(to_sort_reps, function(x) x$coeftable["genotype:scaled_phenotype", "p"]))
beta_reps = unlist(lapply(to_sort_reps, function(x) x$coeftable["genotype:scaled_phenotype", "Est."]))
qvals_reps = qvalue::qvalue(pvals_reps,pi0.method = "bootstrap")$qvalues
qvals_reps = sort(qvals_reps)

hist(pvals_reps)
hist(qvals_reps)
table(pvals_reps<0.05) # 322 untreated phagocytosis,
table(qvals_reps<0.20) # 211 untreated phagocytosis,
table(qvals_reps<0.05) # 26 untreated phagocytosis,


examine_qvals = qvalue::qvalue(pvals_reps, pi0.method = "bootstrap")

toplot = data.frame(pi0.lambda = examine_qvals$pi0.lambda,
                    lambda = examine_qvals$lambda)

plot(examine_qvals)



pdf(paste0(output_dir,"/",treat,"_",phenotype,"qvals_QC_reps.pdf"), height = 7, width = 7)
p1 = ggplot(toplot,aes(y=pi0.lambda,lambda)) + 
  geom_point() + 
  theme_minimal() + 
  geom_hline(yintercept = examine_qvals$pi0, linetype = "dashed") + 
  ggtitle(paste0("pi0 =" ,round(examine_qvals$pi0,3)) )

p2 = hist(examine_qvals)
p3 = qvals_reps %>%
  as_tibble() %>%
  ggplot(.,aes(value)) + geom_histogram() + theme_minimal() + 
  ggtitle("qvalues")

(p1 + p2) / p3

dev.off()


pdf(paste0(output_dir,"/",treat,"_",phenotype,"pvals_hist_reps.pdf"), height = 8, width = 8)
p1 = pvals_reps %>%
  as_tibble() %>%
  ggplot(.,aes(value)) + geom_histogram() + theme_minimal() + 
  ggtitle("pvalues")
p2 = qvals_reps %>%
  as_tibble() %>%
  ggplot(.,aes(value)) + geom_histogram() + theme_minimal() + 
  ggtitle("qvalues")

p1 / p2

dev.off()

pdf(paste0(output_dir,"/top20_",phenotype,"_",treat,"_",prolif,"_interactions_expr_x_pheno_reps.pdf"),width = 10,height = 5
)

for(gene_var in names(qvals_reps)[1:20]){
  
  
  genotype_vals = unique(test2_reps[[gene_var]]$genotype)
  
  # with partial residuals
  
  p1 = interact_plot(my_fit_reps[[gene_var]], 
                     pred = scaled_phenotype, 
                     modx = genotype, plot.points = TRUE,
                     modx.values = genotype_vals,
                     partial.residuals = TRUE,
                     interval = FALSE,
                     colors = "Qual1") + 
    ggtitle(paste0(gene_var) ) +
    ylab("Gene expression") + 
    xlab("Phenotype")
  
  
  
  plot(p1)
}
dev.off()


### pvalue correlations

toplot = data.frame(pvals = pvals,pvals_reps = pvals_reps[names(pvals)],gene_var = names(pvals), beta = beta[names(pvals)],
                    beta_reps = beta_reps[names(pvals)], beta_diff = beta_reps[names(pvals)] - beta[names(pvals)] )

p1 = ggplot(toplot,aes(x = -log10(pvals), y = -log10(pvals_reps))) +
  geom_point() + 
  theme_minimal() + 
  geom_abline() + 
  ggtitle("Interact. pvalue correlation")
p1

# beta correlations

p2 = ggplot(toplot,aes(x = beta, y = beta_reps)) +
  geom_point() + 
  theme_minimal() + 
  geom_abline() + 
  ggtitle("Interact. beta correlation")
p2

# do per significance bin
p3 = ggplot(toplot[toplot$pvals_reps<0.05,],aes(x = beta, y = beta_reps)) +
  geom_point() + 
  theme_minimal() + 
  geom_abline() + 
  ggtitle("Interact. pvals_reps < 0.05")
p3

p4 = toplot %>%
  ggplot(aes(x = sort(-log10(pvals_reps)), y = beta_diff)) +
  geom_point() + 
  theme_minimal() +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
  ggtitle("Interact. beta diff at sorted p-vals (reps) ")
p4

pdf(paste0(output_dir,"/reps_noreps_correlations_pval_betas.pdf"), height = 10,
    width = 10)
plot((p1 + p2) / (p3 + p4))
dev.off()
# check correlation of expression with phago for each genotype category
# permutations

###################
###################
###################

### now check positive, negative or ambiguous expr x phenotype correlations ###########


lines_crossing_xy = list()
slope_sign = list()
for(n in names(my_fit)){
  
  genotype_vals = unique(data[[n]]$genotype)
  
  slopes = jtools::sim_slopes(my_fit[[n]],
                     pred = scaled_phenotype, 
                     modx = genotype,
                     modx.values = genotype_vals,
                     cond.int = TRUE )
  # calculate where genotype regression lines cross
  lines_crossing_xy[[n]] = intersection(slop = slopes$slopes,inter = slopes$ints)
  slope_sign[[n]] = c("is_50percent_slope_sign" = (sum(slopes$slopes$p<0.05) / length(slopes$slopes$p))>0.50)
}

lines_crossing_xy = do.call("rbind",lines_crossing_xy)
lines_crossing_xy = as.data.frame(lines_crossing_xy)
lines_crossing_xy$gene_var = rownames(lines_crossing_xy)
lines_crossing_xy = lines_crossing_xy[names(pvals),]
lines_crossing_xy = na.omit(lines_crossing_xy)
slope_sign = do.call("rbind",slope_sign) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_var")


sum(lines_crossing_xy$x < -1)
sum(lines_crossing_xy$x > 1)
sum(lines_crossing_xy$x > -1 & lines_crossing_xy$x < 1)

sum(lines_crossing_xy$x < -2)
sum(lines_crossing_xy$x > 2)
sum(lines_crossing_xy$x > -2 & lines_crossing_xy$x < 2)
# plot sorted q values vs value of x crossing
lines_crossing_xy$pvals = pvals[lines_crossing_xy$gene_var]
p1 = lines_crossing_xy %>%
  dplyr::arrange(pvals) %>%
  dplyr::filter(x>-100 & x<2000) %>% # remove two extreme outliers
  ggplot(., aes(x = -log10(pvals), y = x)) + 
  geom_point() + 
  geom_hline(yintercept = -1) +
  geom_hline(yintercept = 1) +
  theme_minimal() + 
  ylab("log10(phenotype value) where genotype regression lines intersect")

p1  

# where the intersections are very far away from 0, this means the slopes are minimal (not very different from 0)
p1 = lines_crossing_xy %>%
  dplyr::arrange(pvals) %>%
  dplyr::filter(x>-100 & x<2000) %>% # remove two extreme outliers
  ggplot(., aes(x = -log10(pvals), y = x)) + 
  geom_point() + 
  geom_hline(yintercept = -1) +
  geom_hline(yintercept = 1) +
  theme_minimal() + 
  ylab("log10(phenotype value) where genotype regression lines intersect")

p1  


slopes = sim_slopes(my_fit$`TDRD3-13_60392368_G_C`,
                    pred = scaled_phenotype, 
                    modx = genotype,
                    modx.values = genotype_vals,
                    cond.int = TRUE)
intersection(slop = slopes$slopes,inter = slopes$ints)

GOLT1B-12_21505670_G_T

to_sort$`TDRD3-13_60392368_G_C`

genotype_vals = unique(data$`TDRD3-13_60392368_G_C`$genotype)

interactions::sim_slopes(my_fit$`TDRD3-13_60392368_G_C`,
                         pred = scaled_phenotype, 
                         modx = genotype,
                         modx.values = genotype_vals,
                         johnson_neyman = TRUE)
probe_interaction(my_fit$`TDRD3-13_60392368_G_C`,
                   pred = scaled_phenotype, 
                   modx = genotype,
                  modx.values = genotype_vals,
                  cond.int = TRUE,
                  interval = TRUE,  jnplot = TRUE)
genotype_vals = unique(data$`SLC35F6-2_26928687_G_GTA`$genotype)

probe_interaction(my_fit$`SLC35F6-2_26928687_G_GTA`,
                  pred = scaled_phenotype, 
                  modx = genotype,
                  modx.values = genotype_vals,
                  cond.int = TRUE,
                  interval = TRUE,  jnplot = TRUE)

genotype_vals = unique(data$`CNRIP1-2_68353338_A_G`$genotype)



coef_names = c("genotype" = "genotype", "scaled phenotype" = "scaled_phenotype"  , "interaction" = "genotype:scaled_phenotype")
jtools::plot_summs(my_fit$`TDRD3-13_60392368_G_C`, coefs = coef_names)
jtools::plot_coefs(my_fit$`TDRD3-13_60392368_G_C`, coefs = coef_names)


 probe_interaction(my_fit$`GOLT1B-12_21505670_G_T`,
                           pred = scaled_phenotype, 
                           modx = genotype,
                           modx.values = genotype_vals,
                           cond.int = TRUE,
                           interval = TRUE,  jnplot = TRUE)
 gene_var = "GOLT1B-12_21505670_G_T"
 interact_plot(my_fit$`GOLT1B-12_21505670_G_T`, 
               pred = scaled_phenotype, 
               modx = genotype, plot.points = TRUE,
               modx.values = genotype_vals,
               partial.residuals = TRUE,
               interval = FALSE,
               colors = "Qual1") + 
   ggtitle(paste0(gene_var) ) +
   ylab("Gene expression") + 
   xlab("Phenotype")
