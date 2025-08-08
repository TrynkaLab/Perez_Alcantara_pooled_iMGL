# Diff expr high vs low PRS
# Differential expression with limma+voom

library(limma)
library(edgeR)
library(tidyverse) # For ggplot2 and easy manipulation of data
library(patchwork) # To combine plots
set.seed(123)
source("./helpers.R")
outdir="../../data/results/5.1.2.Diff_expr_limma_high_low_PRS"
dir.create(outdir, recursive = T)


#### pool as random effect
# load PRS classification
#################
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

### treatment with prs interaction
pseudobulk =readr::read_tsv("../../data/results/5.prepare_differential_expression_pseudobulks/combined_pseudobulk_raw_counts_treatmentxprolifxidxpool.txt") %>%
  tibble::column_to_rownames(var = "gene")
metadata = readr::read_tsv("../../data/results/5.prepare_differential_expression_pseudobulks/metadata_treatmentxprolifxidxpool.txt") %>%
  tibble::column_to_rownames(var = "cols_names") %>%
  dplyr::rename(ncells = count,donor= donor_id) %>%
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
  

length(unique(metadata$line))

pseudobulk = pseudobulk %>%
  dplyr::select(match(metadata$pseudobulk_colnames,colnames(.)))
# same order
identical(colnames(pseudobulk),metadata$pseudobulk_colnames)


discarded = metadata$ncells < 100
table(discarded)

pseudobulk = pseudobulk[,!discarded]
metadata = metadata[!discarded,]
length(unique(metadata$line)) # 194/250, 78%

# remove samples from pool 9
# discarded = metadata$pool %in% c("pool9")
# pseudobulk = pseudobulk[,!discarded]
# metadata = metadata[!discarded,]
# length(unique(metadata$line)) # 167 donors, 189/248 = 67%

### creating DGElist object



dge = edgeR::DGEList(counts=pseudobulk)

# ~    (1/ncells) + pool + treatment)


# remove genes not expressed in 30% samples (min 1 CPM)
isexpr = rowSums(edgeR::cpm(dge) > 1) >= floor(0.3*ncol(dge))
message(paste0("There are ", sum(isexpr), " genes with over 1 count per million in 30% of samples"))
dge = dge[isexpr,,keep.lib.sizes=FALSE]  
dim(dge)


## how many have 0 counts
table(rowSums(dge$counts)==0)
# none

# TMM normalization:
dge = edgeR::calcNormFactors(dge)


### create design matrix
metadata$treatment = factor(metadata$treatment)
metadata$treatment = relevel(metadata$treatment,ref = "untreated")
metadata$line = factor(metadata$line)
metadata$line = relevel(metadata$line,ref = "hegp_3")
mat = model.matrix(~ 0 + treatment +  line + log10(ncells),metadata)
mat = as.data.frame(mat)
qr(mat)$rank
ncol(mat)
is.fullrank(mat)
colnames(mat)
mat = mat[,colSums(mat)>0]
colnames(mat)
is.fullrank(mat)

# drop columns we're not interested in
#mat = mat[,c(1:2,8,169:493)] # leaving shared donor aowh_2 in to 
# estimate mean donor effects (limited to the two donors shared across all pools)
is.fullrank(mat)
colnames(mat) = make.names(colnames(mat))


v = voom(dge, design = mat, plot=TRUE)

cor = duplicateCorrelation(v, design=mat, block = metadata$pool)
cor$consensus # within-pool correlation
# voom weights may have changed:
v2 = voom(dge, design = mat, plot=TRUE, block = metadata$pool, correlation = cor$consensus)
cor = duplicateCorrelation(v2, design = mat, block =  metadata$pool) # extract duplicates again
cor$consensus

# add PRS information and filter to Q1 and Q4
metadata = metadata %>%
  dplyr::left_join(.,prs) %>%
  dplyr::filter(!is.na(source)) %>%
  dplyr::filter(overallHR_quartile %in% c("Q1","Q4")) 

metadata$overallHR_quartile = factor(metadata$overallHR_quartile)
metadata$overallHR_quartile = relevel(metadata$overallHR_quartile,ref = "Q1")

# create metadata matrix again
# if I fit line, it's not full rank, because line is confounded with PRS (subcategory)
# metadata$line = factor(metadata$line, levels = unique(metadata$line))
# newref = unique(metadata$line)[1]
# metadata$line = relevel(metadata$line,ref = newref)
mat = model.matrix(~ 0 + overallHR_quartile*treatment + log10(ncells),metadata)
mat = as.data.frame(mat)
qr(mat)$rank
ncol(mat)
is.fullrank(mat)
colnames(mat)
mat = mat[,colSums(mat)>0]
colnames(mat)
is.fullrank(mat)

# drop columns we're not interested in
#mat = mat[,c(1:2,8,169:493)] # leaving shared donor aowh_2 in to 
# estimate mean donor effects (limited to the two donors shared across all pools)
is.fullrank(mat)
colnames(mat) = make.names(colnames(mat))

# subset v to the remaining elements of the metadata
v2=v2[,metadata$pseudobulk_colnames]

# fit 
fit = lmFit(v2, design=mat,block = metadata$pool, correlation  = cor$consensus)

cont = c( "Q4vQ1" = "overallHR_quartileQ4-overallHR_quartileQ1",
          "InteractionIFN" = "overallHR_quartileQ4.treatmentIFN", # interaction (IFNQ4 - IFNQ1) - (UntreatedQ4-UntreatedQ1)
          "InteractionLPS" =  "overallHR_quartileQ4.treatmentLPS") # interaction (LPSQ4 - LPSQ1) - (UntreatedQ4-UntreatedQ1)


contrast.matrix = makeContrasts(contrasts = cont,
                                levels=mat)
colnames(contrast.matrix) = c("Q4vQ1","InteractionIFN","InteractionLPS")

fit2 = contrasts.fit(fit, contrast.matrix)
fit2 = eBayes(fit2)

summary_additive = summary(decideTests(fit2,lfc = 1, 
                                       adjust.method = "fdr",p.value=0.05))
summary_additive

res_additive = list()
for(cont in colnames(contrast.matrix)){
  res_additive[[cont]] = topTable(fit2, coef=cont,
                                  number = length(fit2$coefficients)) %>%
    dplyr::arrange(adj.P.Val) %>%
    tibble::rownames_to_column(var = "symbol")
  
  write.table(res_additive[[cont]],paste0(outdir,
                                          "/DiffExpr_",cont,"all_genes_negative_lower_in_former.txt"),
              quote = F, sep = "\t",row.names = F, col.names = T)
}
saveRDS(res_additive,paste0(outdir,"/additive_res.rds"))

### remove pool batch effects for plotting
# similar to partial residuals
v=v[,metadata$pseudobulk_colnames]
pseudobulk_corrected = limma::removeBatchEffect(v$E, batch = metadata$pool, design = mat)

for(cont in colnames(contrast.matrix)){
  g = res_additive[[cont]] %>%
    dplyr::slice_head(n=30) 
  pdf(paste0(outdir,"/DiffExpr_",cont,"top30_abs_log2FC.pdf"),width = 7,height = 7)
  for(i in 1:30){
    message("Plotting top 30 genes sorted by abs(log2FC)")
    expr = as.data.frame(pseudobulk_corrected)
    expr$gene = rownames(expr)
    
    gene_to_plot = expr %>%
      dplyr::filter(gene==g$symbol[i]) %>%
      tidyr::pivot_longer(cols = !gene, 
                          values_to = "voom_count",
                          names_to = "pseudobulk_colnames") %>%
      dplyr::left_join(.,metadata)
    
    
    p = ggplot(gene_to_plot, aes(x=treatment, y=voom_count,col=overallHR_quartile)) +
      geom_point(position=position_dodge(width=.1), size=3,alpha=0.6) + 
      theme_bw() + 
      ggtitle( g$symbol[i]) + 
      ylab("Normalised log2(counts+0.5)")
    plot(p)
  }
  
  dev.off()
}

# separate additive effects per treatment #############
#######################################################

for(treat in c("untreated","IFN","LPS")){
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
  
  ##############
  ### treatment with prs interaction ##########
  ############################
  
  pseudobulk =readr::read_tsv("../../data/results/5.prepare_differential_expression_pseudobulks/combined_pseudobulk_raw_counts_treatmentxprolifxidxpool.txt") %>%
    tibble::column_to_rownames(var = "gene")
  metadata = readr::read_tsv("../../data/results/5.prepare_differential_expression_pseudobulks/metadata_treatmentxprolifxidxpool.txt") %>%
    tibble::column_to_rownames(var = "cols_names") %>%
    dplyr::rename(ncells = count,donor= donor_id) %>%
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
    dplyr::filter(proliferation_status == "Not_proliferating" & treatment == treat) %>%
    dplyr::rename(line = donor) 
  
  
  length(unique(metadata$line))
  
  pseudobulk = pseudobulk %>%
    dplyr::select(match(metadata$pseudobulk_colnames,colnames(.)))
  # same order
  identical(colnames(pseudobulk),metadata$pseudobulk_colnames)
  
  
  discarded = metadata$ncells < 100
  table(discarded)
  
  pseudobulk = pseudobulk[,!discarded]
  metadata = metadata[!discarded,]
  length(unique(metadata$line)) # 190/250 untreated
  
  # remove samples from pool 9
  # discarded = metadata$pool %in% c("pool9")
  # pseudobulk = pseudobulk[,!discarded]
  # metadata = metadata[!discarded,]
  # length(unique(metadata$line)) # 167 donors, 189/248 = 67%
  
  ### creating DGElist object
  
  dge = edgeR::DGEList(counts=pseudobulk)
  
  # ~    (1/ncells) + pool + treatment)
  
  
  # remove genes not expressed in 30% samples (min 1 CPM)
  isexpr = rowSums(edgeR::cpm(dge) > 1) >= floor(0.3*ncol(dge))
  message(paste0("There are ", sum(isexpr), " genes with over 1 count per million in 30% of samples"))
  dge = dge[isexpr,,keep.lib.sizes=FALSE]  
  dim(dge)
  
  
  ## how many have 0 counts
  table(rowSums(dge$counts)==0)
  # none
  
  # TMM normalization:
  dge = edgeR::calcNormFactors(dge)
  
  
  ### create design matrix
  metadata$treatment = factor(metadata$treatment)
  metadata$treatment = relevel(metadata$treatment,ref =treat)
  metadata$line = factor(metadata$line)
  metadata$line = relevel(metadata$line,ref = "hegp_3")
  mat = model.matrix(~ 0 +  line + log10(ncells),metadata)
  mat = as.data.frame(mat)
  qr(mat)$rank
  ncol(mat)
  is.fullrank(mat)
  colnames(mat)
  mat = mat[,colSums(mat)>0]
  colnames(mat)
  is.fullrank(mat)
  
  # drop columns we're not interested in
  #mat = mat[,c(1:2,8,169:493)] # leaving shared donor aowh_2 in to 
  # estimate mean donor effects (limited to the two donors shared across all pools)
  is.fullrank(mat)
  colnames(mat) = make.names(colnames(mat))
  
  
  v = voom(dge, design = mat, plot=TRUE)
  
  cor = duplicateCorrelation(v, design=mat, block = metadata$pool)
  cor$consensus # within-pool correlation
  # voom weights may have changed:
  v2 = voom(dge, design = mat, plot=TRUE, block = metadata$pool, correlation = cor$consensus)
  cor = duplicateCorrelation(v2, design = mat, block =  metadata$pool) # extract duplicates again
  cor$consensus
  
  # add PRS information and filter to Q1 and Q4
  metadata = metadata %>%
    dplyr::left_join(.,prs) %>%
    dplyr::filter(!is.na(source)) %>%
    dplyr::filter(overallHR_quartile %in% c("Q1","Q4")) 
  
  metadata$overallHR_quartile = factor(metadata$overallHR_quartile)
  metadata$overallHR_quartile = relevel(metadata$overallHR_quartile,ref = "Q1")
  
  # create metadata matrix again
  # if I fit line, it's not full rank, because line is confounded with PRS (subcategory)
  # metadata$line = factor(metadata$line, levels = unique(metadata$line))
  # newref = unique(metadata$line)[1]
  # metadata$line = relevel(metadata$line,ref = newref)
  mat = model.matrix(~ 0 + overallHR_quartile + log10(ncells),metadata)
  mat = as.data.frame(mat)
  qr(mat)$rank
  ncol(mat)
  is.fullrank(mat)
  colnames(mat)
  mat = mat[,colSums(mat)>0]
  colnames(mat)
  is.fullrank(mat)
  
  # drop columns we're not interested in
  #mat = mat[,c(1:2,8,169:493)] # leaving shared donor aowh_2 in to 
  # estimate mean donor effects (limited to the two donors shared across all pools)
  is.fullrank(mat)
  colnames(mat) = make.names(colnames(mat))
  
  # subset v to the remaining elements of the metadata
  v2=v2[,metadata$pseudobulk_colnames]
  
  # fit 
  fit = lmFit(v2, design=mat,block = metadata$pool, correlation  = cor$consensus)
  
  cont = c( "Q4vQ1" = "overallHR_quartileQ4-overallHR_quartileQ1") 
  
  contrast.matrix = makeContrasts(contrasts = cont,
                                  levels=mat)
  colnames(contrast.matrix) = c("Q4vQ1")
  
  fit2 = contrasts.fit(fit, contrast.matrix)
  fit2 = eBayes(fit2)
  
  summary_additive = summary(decideTests(fit2,lfc = 0, 
                                         adjust.method = "fdr",p.value=0.05))
  summary_additive
  
  res_additive = list()
  for(cont in colnames(contrast.matrix)){
    res_additive[[cont]] = topTable(fit2, coef=cont,
                                    number = length(fit2$coefficients)) %>%
      dplyr::arrange(adj.P.Val) %>%
      tibble::rownames_to_column(var = "symbol")
    
    write.table(res_additive[[cont]],paste0(outdir,
                                            "/DiffExpr_",cont,"_",treat,"all_genes_negative_lower_in_former.txt"),
                quote = F, sep = "\t",row.names = F, col.names = T)
  }
  saveRDS(res_additive,paste0(outdir,"/",treat,"_additive_res.rds"))
  
  ### remove pool batch effects for plotting
  # similar to partial residuals
  v=v[,metadata$pseudobulk_colnames]
  pseudobulk_corrected = limma::removeBatchEffect(v$E, batch = metadata$pool, design = mat)
  
  for(cont in colnames(contrast.matrix)){
    g = res_additive[[cont]] %>%
      dplyr::slice_head(n=30) 
    pdf(paste0(outdir,"/DiffExpr_",cont,"_",treat,"_top30_abs_log2FC.pdf"),width = 7,height = 7)
    for(i in 1:30){
      message("Plotting top 30 genes sorted by abs(log2FC)")
      expr = as.data.frame(pseudobulk_corrected)
      expr$gene = rownames(expr)
      
      gene_to_plot = expr %>%
        dplyr::filter(gene==g$symbol[i]) %>%
        tidyr::pivot_longer(cols = !gene, 
                            values_to = "voom_count",
                            names_to = "pseudobulk_colnames") %>%
        dplyr::left_join(.,metadata)
      
      
      p = ggplot(gene_to_plot, aes(x=treatment, y=voom_count,col=overallHR_quartile)) +
        geom_point(position=position_dodge(width=.1), size=3,alpha=0.6) + 
        theme_bw() + 
        ggtitle( g$symbol[i]) + 
        ylab("Normalised log2(counts+0.5)")
      plot(p)
    }
    
    dev.off()
  }
}

# APOE4/ APOE4 vs APOE2 / APOE2 #############
#######################################################

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
  dplyr::mutate(donor =  str_split_i(line,pattern = "_",i=1),
                apoe_status = paste(apoeAllele1,apoeAllele2,sep = "/"))

old_prs = read_tsv("../../../hipsci_genotype_processing/data/old_prs_info/AD_polygenic_hazard.cell_lines.recoded.txt") %>%
  dplyr::mutate(line = str_split_i(sampleID,pattern = "-",i=2),
                overall_HR_quartile_old = case_when(overallHR_quantile <= 0.25 ~ "Q1",
                                                    overallHR_quantile <= 0.50 & overallHR_quantile > 0.25 ~ "Q2",
                                                    overallHR_quantile <= 0.75 & overallHR_quantile > 0.50 ~ "Q3",
                                                    overallHR_quantile > 0.75 ~ "Q4"),
                polygenicHR_quartile_old = case_when(polygenicHR_quantile <= 0.25 ~ "Q1",
                                                     polygenicHR_quantile <= 0.50 & polygenicHR_quantile > 0.25 ~ "Q2",
                                                     polygenicHR_quantile <= 0.75 & polygenicHR_quantile > 0.50 ~ "Q3",
                                                     polygenicHR_quantile > 0.75 ~ "Q4"),
                apoe_status_old = paste(apoeAllele1,apoeAllele2,sep = "/")) %>%
  dplyr::select(!sampleID) %>%
  distinct(donor,.keep_all = TRUE) %>%
  dplyr::select(donor,overall_HR_quartile_old,polygenicHR_quartile_old,apoe_status_old)

prs = prs %>%
  dplyr::left_join(.,old_prs)

# check if the new and old PRS are in same quartile
prs = prs %>%
  dplyr::mutate(same_apoe = case_when(apoe_status == apoe_status_old ~ "yes",
                                   .default = "no")) %>%
  dplyr::mutate(keep = case_when(source == "IPMAR" ~ "yes",
                                 source=="HipSci" & same_apoe=="yes" ~ "yes",
                                 .default = "no")) %>%
  dplyr::filter(keep == "yes")

### treatment with prs interaction
pseudobulk =readr::read_tsv("../../data/results/5.prepare_differential_expression_pseudobulks/combined_pseudobulk_raw_counts_treatmentxprolifxidxpool.txt") %>%
  tibble::column_to_rownames(var = "gene")
metadata = readr::read_tsv("../../data/results/5.prepare_differential_expression_pseudobulks/metadata_treatmentxprolifxidxpool.txt") %>%
  tibble::column_to_rownames(var = "cols_names") %>%
  dplyr::rename(ncells = count,donor= donor_id) %>%
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


length(unique(metadata$line))

pseudobulk = pseudobulk %>%
  dplyr::select(match(metadata$pseudobulk_colnames,colnames(.)))
# same order
identical(colnames(pseudobulk),metadata$pseudobulk_colnames)


discarded = metadata$ncells < 100
table(discarded)

pseudobulk = pseudobulk[,!discarded]
metadata = metadata[!discarded,]
length(unique(metadata$line)) # 194/250, 78%

# remove samples from pool 9
# discarded = metadata$pool %in% c("pool9")
# pseudobulk = pseudobulk[,!discarded]
# metadata = metadata[!discarded,]
# length(unique(metadata$line)) # 167 donors, 189/248 = 67%

### creating DGElist object

dge = edgeR::DGEList(counts=pseudobulk)

# ~    (1/ncells) + pool + treatment)


# remove genes not expressed in 30% samples (min 1 CPM)
isexpr = rowSums(edgeR::cpm(dge) > 1) >= floor(0.3*ncol(dge))
message(paste0("There are ", sum(isexpr), " genes with over 1 count per million in 30% of samples"))
dge = dge[isexpr,,keep.lib.sizes=FALSE]  
dim(dge)


## how many have 0 counts
table(rowSums(dge$counts)==0)
# none

# TMM normalization:
dge = edgeR::calcNormFactors(dge)


### create design matrix
metadata$treatment = factor(metadata$treatment)
metadata$treatment = relevel(metadata$treatment,ref = "untreated")
metadata$line = factor(metadata$line)
metadata$line = relevel(metadata$line,ref = "hegp_3")
mat = model.matrix(~ 0 + treatment +  line + log10(ncells),metadata)
mat = as.data.frame(mat)
qr(mat)$rank
ncol(mat)
is.fullrank(mat)
colnames(mat)
mat = mat[,colSums(mat)>0]
colnames(mat)
is.fullrank(mat)

# drop columns we're not interested in
#mat = mat[,c(1:2,8,169:493)] # leaving shared donor aowh_2 in to 
# estimate mean donor effects (limited to the two donors shared across all pools)
is.fullrank(mat)
colnames(mat) = make.names(colnames(mat))


v = voom(dge, design = mat, plot=TRUE)

cor = duplicateCorrelation(v, design=mat, block = metadata$pool)
cor$consensus # within-pool correlation
# voom weights may have changed:
v2 = voom(dge, design = mat, plot=TRUE, block = metadata$pool, correlation = cor$consensus)
cor = duplicateCorrelation(v2, design = mat, block =  metadata$pool) # extract duplicates again
cor$consensus

# add PRS information and filter to double apoe4 and apoe2/2 & 2/3
metadata = metadata %>%
  dplyr::left_join(.,prs) %>%
  dplyr::filter(!is.na(source)) %>%
  dplyr::filter(apoe_status %in% c("e4/e4","e2/e2","e2/e3","e3/e2")) 

metadata = metadata %>%
  dplyr::mutate(apoe_status = case_when(apoe_status=="e3/e2" ~ "e2",
                                        apoe_status=="e2/e3" ~ "e2",
                                        apoe_status=="e2/e2" ~ "e2",
                                        apoe_status=="e4/e4" ~ "e4"))
metadata$apoe_status = factor(metadata$apoe_status)
metadata$apoe_status = relevel(metadata$apoe_status,ref = "e2")

# create metadata matrix again
# if I fit line, it's not full rank, because line is confounded with PRS (subcategory)
# metadata$line = factor(metadata$line, levels = unique(metadata$line))
# newref = unique(metadata$line)[1]
# metadata$line = relevel(metadata$line,ref = newref)
mat = model.matrix(~ 0 + apoe_status*treatment + log10(ncells),metadata)
mat = as.data.frame(mat)
qr(mat)$rank
ncol(mat)
is.fullrank(mat)
colnames(mat)
mat = mat[,colSums(mat)>0]
colnames(mat)
is.fullrank(mat)

# drop columns we're not interested in
#mat = mat[,c(1:2,8,169:493)] # leaving shared donor aowh_2 in to 
# estimate mean donor effects (limited to the two donors shared across all pools)
is.fullrank(mat)
colnames(mat) = make.names(colnames(mat))

# subset v to the remaining elements of the metadata
v2=v2[,metadata$pseudobulk_colnames]

# fit 
fit = lmFit(v2, design=mat,block = metadata$pool, correlation  = cor$consensus)

cont = c( "e4ve2" = "apoe_statuse4-apoe_statuse2",
          "InteractionIFN_apoe" = "apoe_statuse4.treatmentIFN", # interaction (IFNe4 - IFNe2) - (Untreatede4-Untreatede2)
          "InteractionLPS_apoe" =  "apoe_statuse4.treatmentLPS") # interaction (LPSe4 - LPSe2) - (Untreatede4-Untreatede2)


contrast.matrix = makeContrasts(contrasts = cont,
                                levels=mat)
colnames(contrast.matrix) = c("e4ve2","InteractionIFN_apoe","InteractionLPS_apoe")

fit2 = contrasts.fit(fit, contrast.matrix)
fit2 = eBayes(fit2)

summary_additive = summary(decideTests(fit2,adjust.method = "fdr",p.value=0.05))
summary_additive

res_additive = list()
for(cont in colnames(contrast.matrix)){
  res_additive[[cont]] = topTable(fit2, coef=cont,
                                  number = length(fit2$coefficients)) %>%
    dplyr::arrange(adj.P.Val) %>%
    tibble::rownames_to_column(var = "symbol")
  
  write.table(res_additive[[cont]],paste0(outdir,
                                          "/DiffExpr_",cont,"all_genes_negative_lower_in_former.txt"),
              quote = F, sep = "\t",row.names = F, col.names = T)
}
saveRDS(res_additive,paste0(outdir,"/apoe_res.rds"))

### remove pool batch effects for plotting
# similar to partial residuals
v=v[,metadata$pseudobulk_colnames]
pseudobulk_corrected = limma::removeBatchEffect(v$E, batch = metadata$pool, design = mat)

for(cont in colnames(contrast.matrix)){
  g = res_additive[[cont]] %>%
    dplyr::slice_head(n=30) 
  pdf(paste0(outdir,"/DiffExpr_",cont,"top30_abs_log2FC.pdf"),width = 7,height = 7)
  for(i in 1:30){
    message("Plotting top 30 genes sorted by abs(log2FC)")
    expr = as.data.frame(pseudobulk_corrected)
    expr$gene = rownames(expr)
    
    gene_to_plot = expr %>%
      dplyr::filter(gene==g$symbol[i]) %>%
      tidyr::pivot_longer(cols = !gene, 
                          values_to = "voom_count",
                          names_to = "pseudobulk_colnames") %>%
      dplyr::left_join(.,metadata)
    
    
    p = ggplot(gene_to_plot, aes(x=treatment, y=voom_count,col=apoe_status)) +
      geom_point(position=position_dodge(width=.1), size=3,alpha=0.6) + 
      theme_bw() + 
      ggtitle( g$symbol[i]) + 
      ylab("Normalised log2(counts+0.5)")
    plot(p)
  }
  
  dev.off()
}

#################
### with continuous APOE HR ############
####################

