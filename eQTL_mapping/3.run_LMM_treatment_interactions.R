# interaction QTLs for treatments with LMM
library(lme4)
library(lmerTest)
library(interactions)
library(tidyverse)
# try first with significant gene-variant pairs from chromosome 1
## LMM - lm() or something more resistant to outliers
# interaction plots with partial residuals
output_dir = "../../data/results/3.run_LMM_treatment_interactions"
dir.create(output_dir, recursive = TRUE)


prolif = "Not_proliferating"
donor_blacklist="letw_5,lizq_3,zaie_1,romx_2,seru_7,qonc_2,iukl_1,curn_3,boqx_2,garx_2,sojd_3,yoch_6"
pcnumber = 15

donor_blacklist = stringr::str_split(donor_blacklist,pattern = ",")[[1]]

####

# load files


# covariates
metadata = readr::read_tsv("../../../OTAR2065_differentiation_efficiency/data/results/5.prepare_differential_expression_pseudobulks/metadata_treatmentxprolifxidxpool.txt") %>%
  dplyr::mutate(line_pool = paste(donor_id,pool,sep = "_"),
                line_pool_treatment = paste(donor_id,pool,treatment,sep = "_")) %>%
  dplyr::filter(proliferation_status == prolif ) %>%
  dplyr::rename(line=donor_id)
# expression
pseudobulk = readr::read_tsv("../../../OTAR2065_differentiation_efficiency/data/results/5.prepare_differential_expression_pseudobulks/combined_pseudobulk_raw_counts_treatmentxprolifxidxpool.txt") %>%
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


## filters as in regular tensorQTL
if (!identical(colnames(pseudobulk), metadata$cols_names)) {
  stop("Column names of pseudobulk are not identical to metadata$cols_names.")
} else {
  print("Expression and metadata column names match.")
}

message("Filtering by number of cells: remove donors with fewer than 100 cells")
retain = metadata$count >=100
metadata = metadata[retain,]
pseudobulk = pseudobulk[,retain]

# transform counts
metadata = metadata %>%
  dplyr::mutate(one_over_ncells = 1/count)


message("Excluding donors from the blacklist")

retain = !(metadata$line  %in% donor_blacklist)
metadata = metadata[retain,]
pseudobulk = pseudobulk[,retain]

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
colnames(scaled) = c("chr","start","end","gene_id","gene_name",metadata$line_pool_treatment)


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
    rownames_to_column(var = "line_pool_treatment")
  metadata = metadata %>%
    dplyr::left_join(pcs)
}




message("Read in genotype file with significant eQTL results - mashr")

mashr = readRDS("../../data/results/8.4.Shared_eqtls_mashr/mashr_results.rds")
mashr_sign = mashr$result$lfsr[rowSums(mashr$result$lfsr<0.05)==1,]
vars = rownames(mashr_sign)
mashr_sign = mashr_sign %>%
  as.data.frame() %>%
  dplyr::mutate(gene_id_var = vars,
                shared_category = case_when(`60_untreated_Not_proliferating` < 0.05 & `60_IFN_Not_proliferating` > 0.05 & `60_LPS_Not_proliferating` > 0.05 ~ "untreated_only",
                                            `60_untreated_Not_proliferating` > 0.05 & `60_IFN_Not_proliferating` < 0.05 & `60_LPS_Not_proliferating` > 0.05 ~ "IFN_only",
                                            `60_untreated_Not_proliferating` > 0.05 & `60_IFN_Not_proliferating` > 0.05 & `60_LPS_Not_proliferating` < 0.05 ~ "LPS_only",
                                           .default = "shared")) %>%
   dplyr::mutate(gene_id = str_split_i(gene_id_var,pattern = "-",i = 1),
                 variant_id = str_split_i(gene_id_var,pattern = "-",i = -1)) %>%
  dplyr::left_join(ensembl)

  saveRDS(mashr_sign,paste0(output_dir,"/mashr_specific_tensorQTL.rds"))


message("Read in chunked genotype with allele dosages")
genotype = list()
for(n in 1:476){
  genotype[[n]] = readr::read_csv(paste0("../../../OTAR2065_phenotypic_QTLs/data/genotypes/full_genotype/genotype_minor_allele_dosage_",
                                         n,
                                         ".csv")) %>%
    dplyr::relocate(rn) %>%
    dplyr::rename(variant_id = rn) %>%
    dplyr::mutate(variant_id = str_replace(variant_id,pattern="chr",replacement = "")) %>%
    dplyr::filter(variant_id %in% mashr_sign$variant_id)   # filter to shared
  
}
genotype = do.call("rbind",genotype)
saveRDS(genotype,paste0(output_dir,"/genotype_sign_tensorQTL_int_treatments_mashr_",prolif,".rds")) # to save a bit of time

my_fit = list()
test2 = list()
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
    dplyr::left_join(metadata[,c("line","treatment","pool","line_pool","line_pool_treatment",paste0("genotypePC",1:5),
                                 paste0("PC",1:pcnumber),
                                 "one_over_ncells")]) %>%
    dplyr::filter(!is.na(treatment))
  
  # add expression###
  # subset to significant gene from tested variant
  sign_gene = mashr_sign %>%
    dplyr::filter(variant_id == var)
  
  for(gen in unique(sign_gene$gene_name)){
    
    subset_scaled = scaled %>%
      dplyr::filter(gene_name == gen) %>%
      tidyr::pivot_longer(!c("chr","start","end","gene_id","gene_name"),
                          names_to = "line_pool_treatment",values_to = "scaled_expression")
    
    test2[[paste(gen,var,sep = "-")]] = test %>%
      dplyr::left_join(subset_scaled)
    test2[[paste(gen,var,sep = "-")]]$genotype = factor( test2[[paste(gen,var,sep = "-")]]$genotype,
                                                            levels = unique( test2[[paste(gen,var,sep = "-")]]$genotype),
                                                            ordered=TRUE)

    # fit
    if(nrow(subset_scaled)==0){
      message("Not found in scaled gene expression - old eQTL results?")
      next()
      
    }else{
      my_fit[[paste(gen,var,sep = "-")]] = lmerTest::lmer(scaled_expression ~  genotypePC1 + genotypePC2 + genotypePC3 + genotypePC4 + genotypePC5 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 +
                                                            genotype*treatment +  (1|pool) ,
                                                          data = test2[[paste(gen,var,sep = "-")]],
                                                          control = lmerControl(optCtrl=list(maxfun=5000) ))
      
    }
    
    
  }
  
}

to_sort = lapply(my_fit,jtools::summ)
pvals = unlist(lapply(to_sort, function(x) x$coeftable["genotype.L:treatmentuntreated", "p"])) # need to extract every comparison to ref
qvals = qvalue::qvalue(pvals)$qvalues
qvals = sort(qvals)

hist(pvals)
hist(qvals)


gene_var = "CDC37-19_11236817_T_C "


# with partial residuals

interact_plot(my_fit$`CDC37-19_11236817_T_C`, 
              pred = scaled_fraction_mean, 
              modx = genotype, plot.points = TRUE,
              modx.values = unique(test2$`CDC37-19_11236817_T_C`$genotype),
              partial.residuals = TRUE,
              interval = FALSE) + ggtitle(paste0(gene_var),
                                          subtitle = "Accounting for covariates" )
p1 = interact_plot(my_fit$`CDC37-19_11236817_T_C`, 
                   pred = scaled_fraction_mean, 
                   modx = genotype, plot.points = TRUE,
                   partial.residuals = TRUE,modx.values = unique(test2$`CDC37-19_11236817_T_C`$genotype),
                   interval = FALSE, linearity.check = TRUE) + ggtitle(paste0(gene_var),
                                                                       subtitle = "Accounting for covariates" )
p1

ss = sim_slopes(my_fit$`CDC37-19_11236817_T_C`, pred = scaled_fraction_mean, modx = genotype, 
                modx.values =unique(test2$`CDC37-19_11236817_T_C`$genotype))
plot(ss)
# with partial residuals

head(sort(pvals))
head(sort(qvals))


# plotting top 20 most significant results  
pdf(paste0(output_dir,"/top20_IFN_interactions_expr_tensorQTL_interactions.pdf"),width = 10,height = 5
)

for(gene_var in names(qvals)[1:20]){
  
  genotype_vals = unique(test2[[gene_var]]$genotype)
  shared_cat = mashr_sign %>% filter(variant_id ==  unique(test2[[gene_var]]$variant_id)) %>% select(shared_category) %>% unlist()
  p1 = cat_plot(my_fit[[gene_var]], pred = genotype, modx = treatment, plot.points = TRUE,partial.residuals = FALSE) +
    ggtitle(paste0(gene_var), subtitle = shared_cat)
  p2 = cat_plot(my_fit[[gene_var]], pred = genotype, modx = treatment, plot.points = TRUE,partial.residuals = TRUE) +
    ggtitle(label = "",
             subtitle = paste0("Accounting for covariates. ",shared_cat) )
 
  plot(patchwork::wrap_plots(p1, p2,guides = "collect"))
}
dev.off()
