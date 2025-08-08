# interactions of expr, phenotype and PRS in LMM
library(lme4)
library(lmerTest)
library(interactions)
library(tidyverse)
library(variancePartition)
library(ggrepel)
source("./helpers.R")

output_dir = "../../data/results/3.LMM_expr_phenotype_PRS_interactions"
dir.create(output_dir, recursive = TRUE)

phenotype = "phagocytosis"
treat = "untreated"
prolif = "Not_proliferating"

####
for(phenotype in c("phagocytosis","migration")){
 for(treat in c("untreated", "IFN", "LPS")){
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
  
  # load PRS classification
  prs = read_tsv("../../../hipsci_genotype_processing/data/prs_ad_bellenguez/AD_polygenic_hazard.hipsci_ipmar_donors.txt") %>%
    dplyr::mutate(donor = case_when(donor=="Arlene-003" ~ "Arlene_3",
                                    donor=="Cindy-005" ~ "Cindy_5",
                                    donor=="Dexter-006" ~ "Dexter_6",
                                    donor=="Fiona-010" ~ "Fiona_10",
                                    donor=="Gilma-009" ~ "Gilma_9",
                                    donor=="Hector-011" ~ "Hector_11",
                                    donor=="Imani-012" ~ "Imani_12",
                                    donor=="Javier-013" ~ "Javier_13",
                                    donor=="Keoni-014" ~ "Keoni_14",
                                    donor=="Olaf-018" ~ "Olaf_18",
                                    donor=="Bertha-004" ~ "Bertha_4",
                                    donor=="Mindy-016" ~ "Mindy_16",
                                    .default = donor)) %>%
    dplyr::rename(line=donor) %>%
    dplyr::mutate(donor =  str_split_i(line,pattern = "_",i=1))
  old_prs = read_tsv("../../../hipsci_genotype_processing/data/old_prs_info/AD_polygenic_hazard.cell_lines.recoded.txt") %>%
    dplyr::mutate(line = str_split_i(sampleID,pattern = "-",i=2),
                  overall_HR_quartile_old = case_when(overallHR_quantile <= 0.25 ~ "Q1",
                                                      overallHR_quantile <= 0.50 & overallHR_quantile > 0.25 ~ "Q2",
                                                      overallHR_quantile <= 0.75 & overallHR_quantile > 0.50 ~ "Q3",
                                                      overallHR_quantile > 0.75 ~ "Q4"),
                  polygenicHR_quartile_old = case_when(polygenicHR_quantile <= 0.25 ~ "Q1",
                                                       polygenicHR_quantile <= 0.50 & polygenicHR_quantile > 0.25 ~ "Q2",
                                                       polygenicHR_quantile <= 0.75 & polygenicHR_quantile > 0.50 ~ "Q3",
                                                       polygenicHR_quantile > 0.75 ~ "Q4")) %>%
    dplyr::select(!sampleID) %>%
    distinct(donor,.keep_all = TRUE) %>%
    dplyr::select(donor,overall_HR_quartile_old,polygenicHR_quartile_old)
  
  prs = prs %>%
    dplyr::left_join(.,old_prs)
  
  # check if the new and old PRS are in same quartile
  prs = prs %>%
    dplyr::mutate(same_q = case_when(overallHR_quartile == overall_HR_quartile_old ~ "yes",
                                     .default = "no")) %>%
    dplyr::mutate(keep = case_when(source == "IPMAR" ~ "yes",
                                   source=="HipSci" & same_q=="yes" ~ "yes",
                                   .default = "no")) %>%
    dplyr::filter(keep == "yes")
  
  ### expression
  pseudobulk =readr::read_tsv("../../data/results/5.prepare_differential_expression_pseudobulks/combined_pseudobulk_raw_counts_treatmentxprolifxidxpool.txt") %>%
    tibble::column_to_rownames(var = "gene")
  metadata = readr::read_tsv("../../data/results/5.prepare_differential_expression_pseudobulks/metadata_treatmentxprolifxidxpool.txt") %>%
    tibble::column_to_rownames(var = "cols_names") %>%
    dplyr::filter(treatment == treat) %>%
    dplyr::mutate(line_pool = paste(donor_id,pool,sep = "_")) %>%
    dplyr::rename(ncells = count,donor= donor_id) %>%
    dplyr::mutate(log10_ncells = log10(ncells)) %>%
    dplyr::mutate(pseudobulk_colnames=paste(treatment,proliferation_status,donor,pool,sep = "_"),
                  donor = case_when(donor=="Arlene-003" ~ "Arlene_3",
                                    donor=="Cindy-005" ~ "Cindy_5",
                                    donor=="Dexter-006" ~ "Dexter_6",
                                    donor=="Fiona-010" ~ "Fiona_10",
                                    donor=="Gilma-009" ~ "Gilma_9",
                                    donor=="Hector-011" ~ "Hector_11",
                                    donor=="Imani-012" ~ "Imani_12",
                                    donor=="Javier-013" ~ "Javier_13",
                                    donor=="Keoni-014" ~ "Keoni_14",
                                    donor=="Olaf-018" ~ "Olaf_18",
                                    donor=="Bertha-004" ~ "Bertha_4",
                                    donor=="Mindy-016" ~ "Mindy_16",
                                    .default = donor))
  
  # subset to non-proliferating only - do not subset further yet to correct for pool effects with full information
  metadata = metadata %>%
    dplyr::filter(proliferation_status == "Not_proliferating") %>%
    dplyr::rename(line = donor) 
  
  # add phenotype
  metadata = metadata %>%
    dplyr::left_join(mean_pheno[,c("line_pool","scaled_phenotype")],by = "line_pool")
  
  # add PRS information and filter to Q1 and Q4
  metadata = metadata %>%
    dplyr::left_join(.,prs) %>%
    dplyr::filter(!is.na(source)) %>%
    dplyr::filter(overallHR_quartile %in% c("Q1","Q4")) 
  
  # remove rows without phenotype data
  metadata = metadata %>%
    dplyr::filter(!is.na(scaled_phenotype))
  
  length(unique(metadata$line))
  
  pseudobulk = pseudobulk %>%
    dplyr::select(match(metadata$pseudobulk_colnames,colnames(.)))
  # same order
  identical(colnames(pseudobulk),metadata$pseudobulk_colnames)
  
  
  discarded = metadata$ncells < 100
  table(discarded)
  
  pseudobulk = pseudobulk[,!discarded]
  metadata = metadata[!discarded,]
  length(unique(metadata$line)) # 41/250, 16% (per treatment)
  
  
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
  
  
  my_fit = list()
  data = list()
  for(gene in unique(scaled$gene_name)){
    
    # subset gene
    subset_scaled = scaled %>%
      dplyr::filter(gene_name == gene) %>%
      tidyr::pivot_longer(!c("chr","start","end","gene_id","gene_name"),
                          names_to = "line_pool",values_to = "scaled_expression")
    # add metadata
    data[[gene]]  = subset_scaled %>%
      dplyr::left_join(metadata[,c("line","treatment","pool","line_pool",
                                   "scaled_phenotype","log10_ncells", "overallHR_quartile")]) %>%
      dplyr::filter(!is.na(treatment))

      
      # if any category < 5 quartile replicates, remove
      # then if there's only one category, omit
      category_to_omit = names(table(data[[gene]]$overallHR_quartile)[table(data[[gene]]$overallHR_quartile)<5])
      data[[gene]] = data[[gene]] %>%
        dplyr::filter(!overallHR_quartile %in% category_to_omit)
      if(length(unique(data[[gene]]$overallHR_quartile))>1){
        # fit
        if(nrow(subset_scaled)==0){
          message("Error: no rows in scaled expression ")
          next()
          
        }else{
          my_fit[[gene]] = lmerTest::lmer(scaled_expression ~   overallHR_quartile*scaled_phenotype + log10_ncells + (1 | pool),
                                                              data = data[[gene]],
                                                              control = lmerControl(optCtrl=list(maxfun=5000) ),  REML = FALSE)
          
        }
      }else{
        message("Skipping because of insufficient PRS replicates")
        
        next()
      }
      
      
      
    }
  
  saveRDS(my_fit,paste0(output_dir,"/",phenotype,"_",treat,"_",prolif,"_LMM_noreps.rds")) 
  saveRDS(data,paste0(output_dir,"/",phenotype,"_",treat,"_",prolif,"_data_noreps.rds"))
  
 }
}

for(phenotype in c("phagocytosis","migration")){
  for(treat in c("untreated", "IFN", "LPS")){
    
my_fit = readRDS(paste0(output_dir,"/",phenotype,"_",treat,"_",prolif,"_LMM_noreps.rds"))
data = readRDS(paste0(output_dir,"/",phenotype,"_",treat,"_",prolif,"_data_noreps.rds"))
to_sort = lapply(my_fit,jtools::summ,r.squared = FALSE)
pvals = unlist(lapply(to_sort, function(x) x$coeftable["overallHR_quartileQ4:scaled_phenotype", "p"]))
beta = unlist(lapply(to_sort, function(x) x$coeftable["overallHR_quartileQ4:scaled_phenotype", "Est."]))

qvals = qvalue::qvalue(pvals, pi0.method = "bootstrap")$qvalues
qvals = sort(qvals)

hist(pvals)
hist(qvals)
table(pvals<0.10) # 1534 IFN phagocytosis
table(pvals<0.05) # 874 IFN phagocytosis, x untreated migration,  x IFN migration, x phagocytosis LPS, x migration LPS
table(qvals<0.20) # 7 IFN phagocytosis, x untreated migration,  x IFN migration, x phagocytosis LPS,x migration LPS
table(qvals<0.10) # 5 IFN phagocytosis,
table(qvals<0.05) # 2 IFN phagocytosis,

qvals = sort(qvals)
sign = names(qvals)[qvals<0.20]

tosave = data.frame(pvals = pvals[names(qvals)],qvals = qvals,
                    gene = names(qvals), beta = beta[names(qvals)])

write_csv(tosave,
          paste0(output_dir,"/all_",phenotype,"_",treat,"_",prolif,"_interactions_expr_x_pheno_PRS_noreps.csv"))

pdf(paste0(output_dir,"/all_",phenotype,"_",treat,"_",prolif,"_interactions_expr_x_pheno_PRS_noreps.pdf"),
    width = 7,height = 7
)

for(gene in sign){
  
  

  # with partial residuals
  
  p1 = interact_plot(my_fit[[gene]], 
                     pred = scaled_phenotype, 
                     modx = overallHR_quartile, plot.points = TRUE,
                     modx.values = c("Q4","Q1"),
                     partial.residuals = TRUE,
                     interval = FALSE,
                     colors = "Qual1") + 
    ggtitle(paste0(gene) ) +
    ylab("Gene expression") + 
    xlab("Phenotype")
  # with partial residuals
  
  
  plot(p1)
}
dev.off()
  }}
