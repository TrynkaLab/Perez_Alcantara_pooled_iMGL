# differential expression with PRS as continuous variable
# Diff expr high vs low PRS
# Differential expression with limma+voom

library(limma)
library(edgeR)
library(tidyverse) # For ggplot2 and easy manipulation of data
library(patchwork) # To combine plots
library(PCAtools)
library(variancePartition)
set.seed(123)
source("./helpers.R")
outdir="../../data/results/5.1.2.Diff_expr_limma_high_low_PRS"
dir.create(outdir, recursive = T)

gene_info = read_csv("/lustre/scratch123/hgi/teams/trynka/resources/biomart/Homo_sapiens.GRCh38.111.genes.csv") %>%
  dplyr::rename(gene = gene_name) %>%
  dplyr::filter(gene_biotype == "protein_coding")
### treatment with prs interaction
pseudobulk =readr::read_tsv("../../data/results/5.prepare_differential_expression_pseudobulks/combined_pseudobulk_raw_counts_treatmentxprolifxidxpool.txt") %>%
  dplyr::filter(gene %in% gene_info$gene) %>%
  tibble::column_to_rownames(var = "gene") 

metadata = readr::read_tsv("../../data/results/5.prepare_differential_expression_pseudobulks/metadata_treatmentxprolifxidxpool.txt") %>%
  tibble::column_to_rownames(var = "cols_names") %>%
  dplyr::rename(ncells = count,donor= donor_id) %>%
  dplyr::mutate(pseudobulk_colnames=paste(treatment,proliferation_status,donor,pool,sep = "_"),
                donor = case_when(donor=="Arlene-003" ~ "Arlene",
                                  donor=="Cindy-005" ~ "Cindy",
                                  donor=="Dexter-006" ~ "Dexter",
                                  donor=="Fiona-010" ~ "Fiona",
                                  donor=="Gilma-009" ~ "Gilma",
                                  donor=="Hector-011" ~ "Hector",
                                  donor=="Imani-012" ~ "Imani",
                                  donor=="Javier-013" ~ "Javier",
                                  donor=="Keoni-014" ~ "Keoni",
                                  donor=="Olaf-018" ~ "Olaf",
                                  donor=="Bertha-004" ~ "Bertha",
                                  donor=="Mindy-016" ~ "Mindy",
                                  donor=="Qiana-022" ~ "Qiana",
                                  donor=="Nestor-017" ~ "Nestor",
                                  .default = donor))

# subset to non-proliferating only
metadata = metadata %>%
  dplyr::filter(proliferation_status == "Not_proliferating") %>%
  dplyr::rename(line = donor) 


length(unique(metadata$line)) # 250



# load sex and age information
add_metadata = read.csv("../../data/metadata_info_hipsci_IPMAR.csv") %>%
  dplyr::mutate(line = case_when(donor=="Arlene-003" ~ "Arlene",
                                 donor=="Cindy-005" ~ "Cindy",
                                 donor=="Dexter-006" ~ "Dexter",
                                 donor=="Fiona-010" ~ "Fiona",
                                 donor=="Gilma-009" ~ "Gilma",
                                 donor=="Hector-011" ~ "Hector",
                                 donor=="Imani-012" ~ "Imani",
                                 donor=="Javier-013" ~ "Javier",
                                 donor=="Keoni-014" ~ "Keoni",
                                 donor=="Olaf-018" ~ "Olaf",
                                 donor=="Bertha-004" ~ "Bertha",
                                 donor=="Mindy-016" ~ "Mindy",
                                 donor=="Qiana-022" ~ "Qiana",
                                 donor=="Nestor-017" ~ "Nestor",
                                 
                                 
                                 .default = line)) %>%
  dplyr::select(line,Sex,Age)

# converting NA to unknown to avoid problems with missing values
metadata = metadata %>%
  dplyr::left_join(.,add_metadata) %>%
  dplyr::mutate(Age = case_when(Age == "" ~ 0,  # relevel to numeric
                                is.na(Age) ~ 0,
                                Age == "25-29" ~ 1,
                                Age == "30-34" ~ 2,
                                Age =="35-39" ~ 3,
                                Age =="40-44" ~ 4,
                                Age =="45-49" ~ 5,
                                Age =="50-54" ~ 6,
                                Age =="55-59" ~ 7,
                                Age =="60-64" ~ 8,
                                Age =="65-69" ~ 9,
                                Age =="70-74" ~ 10,
                                Age =="75-79" ~ 11,
                                Age == "80"  ~ 12,
                                Age == "89"  ~ 13,
                                Age == "95"  ~ 15,
                                .default = 0),
                Sex = case_when(Sex == "" ~ "unkown",
                                is.na(Sex) ~ "unknown",
                                .default = Sex)) %>%
  distinct()
# Age will be partly confounded with PRS and sex most likely

pseudobulk = pseudobulk %>%
  dplyr::select(match(metadata$pseudobulk_colnames,colnames(.)))
# same order
identical(colnames(pseudobulk),metadata$pseudobulk_colnames)
# load PRS classification
#################
prs = read_tsv("../../../hipsci_genotype_processing/data/prs_ad_bellenguez/PRSice/hipsci_polygenic_score_AD_Bellenguez_withAPOE.tsv") %>%
  dplyr::select(-FID)

metadata = metadata %>%
  left_join(prs)
##### filter by number of cells


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
metadata$Sex = factor(metadata$Sex)
metadata$Sex = relevel(metadata$Sex,ref = "Female")


mat = model.matrix(~ 0 + treatment   +  log10(ncells),metadata)
mat = as.data.frame(mat)
qr(mat)$rank
ncol(mat)
is.fullrank(mat)
colnames(mat)
mat = mat[,colSums(mat)>0]
colnames(mat)
is.fullrank(mat)

is.fullrank(mat)
colnames(mat) = make.names(colnames(mat))


v = voom(dge, design = mat, plot=TRUE)

cor = duplicateCorrelation(v, design=mat, block = metadata$pool)
cor$consensus # within-pool correlation
# voom weights may have changed:
v2 = voom(dge, design = mat, plot=TRUE, block = metadata$pool, correlation = cor$consensus)
cor = duplicateCorrelation(v2, design = mat, block =  metadata$pool) # extract duplicates again
cor$consensus


# create metadata matrix again
# if I fit line, it's not full rank, because line is confounded with PRS (subcategory)
mat = model.matrix(~ 0 + treatment*full_PRS + Sex + log10(ncells),data = metadata)
mat = as.data.frame(mat)
qr(mat)$rank
ncol(mat)
is.fullrank(mat)
colnames(mat)
mat = mat[,colSums(mat)!=0]
colnames(mat)
is.fullrank(mat)

colnames(mat) = make.names(colnames(mat))

# subset v to the remaining elements of the metadata
v2=v2[,metadata$pseudobulk_colnames]

# fit 
fit = lmFit(v2, design=mat,block = metadata$pool, correlation  = cor$consensus)

cont = c( "PRS" = "full_PRS", # changes in gene expression per unit increase in PRS
          "InteractionIFN" = "treatmentIFN.full_PRS", # interaction (changes in gene expression in IFN per unit increase in PRS) - (changes in gene expression in untreated per unit increase in PRS)
          "InteractionLPS" =  "treatmentLPS.full_PRS") # interaction for LPS


contrast.matrix = makeContrasts(contrasts = cont,
                                levels=mat)
colnames(contrast.matrix) = c("PRS","InteractionIFN","InteractionLPS")

fit2 = contrasts.fit(fit, contrast.matrix)
fit2 = eBayes(fit2)

summary_additive = summary(decideTests(fit2,lfc = 0, 
                                       adjust.method = "fdr",p.value=0.05))
summary_additive
# 0 interactions wether I fit Age and sex or not

res_additive = list()
for(cont in colnames(contrast.matrix)){
  res_additive[[cont]] = topTable(fit2, coef=cont,
                                  number = length(fit2$coefficients)) %>%
    dplyr::arrange(adj.P.Val) %>%
    tibble::rownames_to_column(var = "symbol")
  
  write.csv(res_additive[[cont]],paste0(outdir,
                                        "/DiffExpr_",cont,".csv"),
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
  pdf(paste0(outdir,"/DiffExpr_",cont,"top30_sign.pdf"),width = 7,height = 7)
  message("Plotting top 30 genes sorted by abs(log2FC)")
  
  for(i in 1:30){
    expr = as.data.frame(pseudobulk_corrected)
    expr$gene = rownames(expr)
    
    gene_to_plot = expr %>%
      dplyr::filter(gene==g$symbol[i]) %>%
      tidyr::pivot_longer(cols = !gene, 
                          values_to = "voom_count",
                          names_to = "pseudobulk_colnames") %>%
      dplyr::left_join(.,metadata)
    
    
    p = ggplot(gene_to_plot, aes(x=full_PRS, y=voom_count,col=treatment)) +
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
  message("Working on ", treat)
  pseudobulk =readr::read_tsv("../../data/results/5.prepare_differential_expression_pseudobulks/combined_pseudobulk_raw_counts_treatmentxprolifxidxpool.txt") %>%
    dplyr::filter(gene %in% gene_info$gene) %>%
    tibble::column_to_rownames(var = "gene") 
  
  metadata = readr::read_tsv("../../data/results/5.prepare_differential_expression_pseudobulks/metadata_treatmentxprolifxidxpool.txt") %>%
    tibble::column_to_rownames(var = "cols_names") %>%
    dplyr::rename(ncells = count,donor= donor_id) %>%
    dplyr::mutate(pseudobulk_colnames=paste(treatment,proliferation_status,donor,pool,sep = "_"),
                  donor = case_when(donor=="Arlene-003" ~ "Arlene",
                                    donor=="Cindy-005" ~ "Cindy",
                                    donor=="Dexter-006" ~ "Dexter",
                                    donor=="Fiona-010" ~ "Fiona",
                                    donor=="Gilma-009" ~ "Gilma",
                                    donor=="Hector-011" ~ "Hector",
                                    donor=="Imani-012" ~ "Imani",
                                    donor=="Javier-013" ~ "Javier",
                                    donor=="Keoni-014" ~ "Keoni",
                                    donor=="Olaf-018" ~ "Olaf",
                                    donor=="Bertha-004" ~ "Bertha",
                                    donor=="Mindy-016" ~ "Mindy",
                                    donor=="Qiana-022" ~ "Qiana",
                                    donor=="Nestor-017" ~ "Nestor",
                                    
                                    .default = donor))
  
  # subset to treatment and non-proliferating only
  metadata = metadata %>%
    dplyr::filter(treatment == treat) %>%
    dplyr::filter(proliferation_status == "Not_proliferating") %>%
    dplyr::rename(line = donor) 
  
  message("Number of donors in ", treat, " ...")
  print(length(unique(metadata$line))) # 241 untreated
  
  
  
  # load sex and age information
  add_metadata = read.csv("../../data/metadata_info_hipsci_IPMAR.csv") %>%
    dplyr::mutate(line = case_when(donor=="Arlene-003" ~ "Arlene",
                                   donor=="Cindy-005" ~ "Cindy",
                                   donor=="Dexter-006" ~ "Dexter",
                                   donor=="Fiona-010" ~ "Fiona",
                                   donor=="Gilma-009" ~ "Gilma",
                                   donor=="Hector-011" ~ "Hector",
                                   donor=="Imani-012" ~ "Imani",
                                   donor=="Javier-013" ~ "Javier",
                                   donor=="Keoni-014" ~ "Keoni",
                                   donor=="Olaf-018" ~ "Olaf",
                                   donor=="Bertha-004" ~ "Bertha",
                                   donor=="Mindy-016" ~ "Mindy",
                                   donor=="Qiana-022" ~ "Qiana",
                                   donor=="Nestor-017" ~ "Nestor",
                                   
                                   .default = line)) %>%
    dplyr::select(line,Sex,Age)
  
  # converting NA to unknown to avoid problems with missing values
  metadata = metadata %>%
    dplyr::left_join(.,add_metadata) %>%
    dplyr::mutate(Age = case_when(Age == "" ~ 0,  # relevel to numeric
                                  is.na(Age) ~ 0,
                                  Age == "25-29" ~ 1,
                                  Age == "30-34" ~ 2,
                                  Age =="35-39" ~ 3,
                                  Age =="40-44" ~ 4,
                                  Age =="45-49" ~ 5,
                                  Age =="50-54" ~ 6,
                                  Age =="55-59" ~ 7,
                                  Age =="60-64" ~ 8,
                                  Age =="65-69" ~ 9,
                                  Age =="70-74" ~ 10,
                                  Age =="75-79" ~ 11,
                                  Age == "80"  ~ 12,
                                  Age == "89"  ~ 13,
                                  Age == "95"  ~ 15,
                                  .default = 0),
                  Sex = case_when(Sex == "" ~ "unkown",
                                  is.na(Sex) ~ "unknown",
                                  .default = Sex)) %>%
    distinct()
  
  pseudobulk = pseudobulk %>%
    dplyr::select(match(metadata$pseudobulk_colnames,colnames(.)))
  # same order
  identical(colnames(pseudobulk),metadata$pseudobulk_colnames)
  # load PRS classification
  #################
  prs = read_tsv("../../../hipsci_genotype_processing/data/prs_ad_bellenguez/PRSice/hipsci_polygenic_score_AD_Bellenguez_withAPOE.tsv") %>%
    dplyr::select(-FID)
  
  metadata = metadata %>%
    left_join(prs)
  
  ##### filter by number of cells
  
  discarded = metadata$ncells < 100
  table(discarded)
  
  pseudobulk = pseudobulk[,!discarded]
  metadata = metadata[!discarded,]
  length(unique(metadata$line)) # 190/250, 76% for untreated
  
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
  
  metadata$Sex = factor(metadata$Sex)
  metadata$Sex = relevel(metadata$Sex,ref = "Female")
  
  
  # subset v to the remaining elements of the metadata
  dge=dge[,metadata$pseudobulk_colnames]
  
  
  mat = model.matrix(~ full_PRS + Sex + log10(ncells),metadata)
  mat = as.data.frame(mat)
  qr(mat)$rank
  ncol(mat)
  is.fullrank(mat)
  colnames(mat)
  mat = mat[,colSums(mat)!=0]
  colnames(mat)
  is.fullrank(mat)
  
  colnames(mat) = make.names(colnames(mat))
  
  
  v = voom(dge, design = mat, plot=TRUE)
  
  cor = duplicateCorrelation(v, design=mat, block = metadata$pool)
  cor$consensus # within-pool correlation
  # voom weights may have changed:
  v2 = voom(dge, design = mat, plot=TRUE, block = metadata$pool, correlation = cor$consensus)
  cor = duplicateCorrelation(v2, design = mat, block =  metadata$pool) # extract duplicates again
  cor$consensus
  
  # fit 
  message("Fitting")
  fit = lmFit(v2, design=mat,block = metadata$pool, correlation  = cor$consensus)
  
  cont = c( "PRS" = "full_PRS") # changes in gene expression per unit increase in PRS, per treatment 
  
  
  contrast.matrix = makeContrasts(contrasts = cont,
                                  levels=mat)
  colnames(contrast.matrix) = c("PRS")
  
  fit2 = contrasts.fit(fit, contrast.matrix)
  fit2 = eBayes(fit2)
  
  summary_additive = summary(decideTests(fit2,lfc = 0, # effects seem to be small
                                         adjust.method = "fdr",p.value=0.05))
  print(summary_additive)
  message("Extracting results")
  
  res_additive = list()
  for(cont in colnames(contrast.matrix)){
    res_additive[[cont]] = topTable(fit2, coef=cont,
                                    number = length(fit2$coefficients)) %>%
      dplyr::arrange(adj.P.Val) %>%
      tibble::rownames_to_column(var = "symbol")
    
    write.csv(res_additive[[cont]],paste0(outdir,
                                          "/DiffExpr_",cont,"_",treat,".csv"),
              quote = F, sep = "\t",row.names = F, col.names = T)
  }
  saveRDS(res_additive,paste0(outdir,"/additive_res_",treat,".rds"))
  
  ### remove pool batch effects for plotting
  # similar to partial residuals
  pseudobulk_corrected = limma::removeBatchEffect(v$E, batch = metadata$pool, design = mat)
  
  for(cont in colnames(contrast.matrix)){
    g = res_additive[[cont]] %>%
      dplyr::slice_head(n=30) 
    pdf(paste0(outdir,"/DiffExpr_",cont,"_",treat,"_top30_sign.pdf"),width = 7,height = 7)
    message("Plotting top 30 genes sorted by abs(log2FC)")
    
    for(i in 1:30){
      expr = as.data.frame(pseudobulk_corrected)
      expr$gene = rownames(expr)
      
      gene_to_plot = expr %>%
        dplyr::filter(gene==g$symbol[i]) %>%
        tidyr::pivot_longer(cols = !gene, 
                            values_to = "voom_count",
                            names_to = "pseudobulk_colnames") %>%
        dplyr::left_join(.,metadata)
      
      # slope of limma fit don't match corrected counts, doing fit on the fly
      
      p = ggplot(gene_to_plot, aes(x=full_PRS, y=voom_count)) +
        geom_point( size=3,alpha=0.6) + 
        geom_smooth(method = "lm") +
        ggpubr::stat_cor(method="pearson") +
        theme_bw() + 
        ggtitle( g$symbol[i]) + 
        ylab("Normalised log2(counts+0.5)")
      plot(p)
    }
    
    dev.off()
  }
  
}

# APOE as numeric (APOE PRS) #############
#######################################################


for(treat in c("untreated","IFN","LPS")){
  pseudobulk =readr::read_tsv("../../data/results/5.prepare_differential_expression_pseudobulks/combined_pseudobulk_raw_counts_treatmentxprolifxidxpool.txt") %>%
    dplyr::filter(gene %in% gene_info$gene) %>%
    tibble::column_to_rownames(var = "gene") 
  metadata = readr::read_tsv("../../data/results/5.prepare_differential_expression_pseudobulks/metadata_treatmentxprolifxidxpool.txt") %>%
    tibble::column_to_rownames(var = "cols_names") %>%
    dplyr::rename(ncells = count,donor= donor_id) %>%
    dplyr::mutate(pseudobulk_colnames=paste(treatment,proliferation_status,donor,pool,sep = "_"),
                  donor = case_when(donor=="Arlene-003" ~ "Arlene",
                                    donor=="Cindy-005" ~ "Cindy",
                                    donor=="Dexter-006" ~ "Dexter",
                                    donor=="Fiona-010" ~ "Fiona",
                                    donor=="Gilma-009" ~ "Gilma",
                                    donor=="Hector-011" ~ "Hector",
                                    donor=="Imani-012" ~ "Imani",
                                    donor=="Javier-013" ~ "Javier",
                                    donor=="Keoni-014" ~ "Keoni",
                                    donor=="Olaf-018" ~ "Olaf",
                                    donor=="Bertha-004" ~ "Bertha",
                                    donor=="Mindy-016" ~ "Mindy",
                                    donor=="Qiana-022" ~ "Qiana",
                                    donor=="Nestor-017" ~ "Nestor",
                                    
                                    .default = donor))
  
  # subset to treatment and non-proliferating only
  metadata = metadata %>%
    dplyr::filter(treatment == treat) %>%
    dplyr::filter(proliferation_status == "Not_proliferating") %>%
    dplyr::rename(line = donor) 
  
  
  length(unique(metadata$line))
  
  
  
  # load sex and age information
  add_metadata = read.csv("../../data/metadata_info_hipsci_IPMAR.csv") %>%
    dplyr::mutate(line = case_when(donor=="Arlene-003" ~ "Arlene",
                                   donor=="Cindy-005" ~ "Cindy",
                                   donor=="Dexter-006" ~ "Dexter",
                                   donor=="Fiona-010" ~ "Fiona",
                                   donor=="Gilma-009" ~ "Gilma",
                                   donor=="Hector-011" ~ "Hector",
                                   donor=="Imani-012" ~ "Imani",
                                   donor=="Javier-013" ~ "Javier",
                                   donor=="Keoni-014" ~ "Keoni",
                                   donor=="Olaf-018" ~ "Olaf",
                                   donor=="Bertha-004" ~ "Bertha",
                                   donor=="Mindy-016" ~ "Mindy",
                                   donor=="Qiana-022" ~ "Qiana",
                                   donor=="Nestor-017" ~ "Nestor",
                                   .default = line)) %>%
    dplyr::select(line,Sex,Age)
  
  # converting NA to unknown to avoid problems with missing values
  metadata = metadata %>%
    dplyr::left_join(.,add_metadata) %>%
    dplyr::mutate(Age = case_when(Age == "" ~ 0,  # relevel to numeric
                                  is.na(Age) ~ 0,
                                  Age == "25-29" ~ 1,
                                  Age == "30-34" ~ 2,
                                  Age =="35-39" ~ 3,
                                  Age =="40-44" ~ 4,
                                  Age =="45-49" ~ 5,
                                  Age =="50-54" ~ 6,
                                  Age =="55-59" ~ 7,
                                  Age =="60-64" ~ 8,
                                  Age =="65-69" ~ 9,
                                  Age =="70-74" ~ 10,
                                  Age =="75-79" ~ 11,
                                  Age == "80"  ~ 12,
                                  Age == "89"  ~ 13,
                                  Age == "95"  ~ 15,
                                  .default = 0),
                  Sex = case_when(Sex == "" ~ "unkown",
                                  is.na(Sex) ~ "unknown",
                                  .default = Sex)) %>%
    distinct()
  
  pseudobulk = pseudobulk %>%
    dplyr::select(match(metadata$pseudobulk_colnames,colnames(.)))
  # same order
  identical(colnames(pseudobulk),metadata$pseudobulk_colnames)
  # load PRS classification
  #################
  prs = read_tsv("../../../hipsci_genotype_processing/data/prs_ad_bellenguez/PRSice/hipsci_polygenic_score_AD_Bellenguez_withAPOE.tsv") %>%
    dplyr::select(-FID)
  
  metadata = metadata %>%
    left_join(prs)
  
  ##### filter by number of cells
  
  discarded = metadata$ncells < 100
  table(discarded)
  
  pseudobulk = pseudobulk[,!discarded]
  metadata = metadata[!discarded,]
  length(unique(metadata$line)) # 190/250, 76% for untreated
  
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
  
  metadata$Sex = factor(metadata$Sex)
  metadata$Sex = relevel(metadata$Sex,ref = "Female")
  
  
  
  # subset v to the remaining elements of the metadata
  dge=dge[,metadata$pseudobulk_colnames]
  
  
  mat = model.matrix(~ APOE_sum_scaled + Sex +  log10(ncells),metadata)
  mat = as.data.frame(mat)
  qr(mat)$rank
  ncol(mat)
  is.fullrank(mat)
  colnames(mat)
  mat = mat[,colSums(mat)!=0]
  colnames(mat)
  is.fullrank(mat)
  
  colnames(mat) = make.names(colnames(mat))
  
  
  v = voom(dge, design = mat, plot=TRUE)
  
  cor = duplicateCorrelation(v, design=mat, block = metadata$pool)
  cor$consensus # within-pool correlation
  # voom weights may have changed:
  v2 = voom(dge, design = mat, plot=TRUE, block = metadata$pool, correlation = cor$consensus)
  cor = duplicateCorrelation(v2, design = mat, block =  metadata$pool) # extract duplicates again
  cor$consensus
  
  # fit 
  fit = lmFit(v2, design=mat,block = metadata$pool, correlation  = cor$consensus)
  
  cont = c( "PRS_APOE" = "APOE_sum_scaled") # changes in gene expression per unit increase in PRS, per treatment 
  
  
  contrast.matrix = makeContrasts(contrasts = cont,
                                  levels=mat)
  colnames(contrast.matrix) = c("PRS_APOE")
  
  fit2 = contrasts.fit(fit, contrast.matrix)
  fit2 = eBayes(fit2)
  
  summary_additive = summary(decideTests(fit2,lfc = 0, # effects seem to be small
                                         adjust.method = "fdr",p.value=0.05))
  print(summary_additive)
  
  res_additive = list()
  for(cont in colnames(contrast.matrix)){
    res_additive[[cont]] = topTable(fit2, coef=cont,
                                    number = length(fit2$coefficients)) %>%
      dplyr::arrange(adj.P.Val) %>%
      tibble::rownames_to_column(var = "symbol")
    
    write.csv(res_additive[[cont]],paste0(outdir,
                                          "/DiffExpr_APOE_sum_scaled_",cont,"_",treat,".csv"),
              quote = F, sep = "\t",row.names = F, col.names = T)
  }
  saveRDS(res_additive,paste0(outdir,"/additive_res_APOE_sum_scaled_",treat,".rds"))
  
  ### remove pool batch effects for plotting
  # similar to partial residuals
  pseudobulk_corrected = limma::removeBatchEffect(v$E, batch = metadata$pool, design = mat)
  
  for(cont in colnames(contrast.matrix)){
    g = res_additive[[cont]] %>%
      dplyr::slice_head(n=30) 
    pdf(paste0(outdir,"/DiffExpr_APOE_sum_scaled_",cont,"_",treat,"_top30_sign.pdf"),width = 7,height = 7)
    message("Plotting top 30 genes sorted by abs(log2FC)")
    
    for(i in 1:30){
      expr = as.data.frame(pseudobulk_corrected)
      expr$gene = rownames(expr)
      
      gene_to_plot = expr %>%
        dplyr::filter(gene==g$symbol[i]) %>%
        tidyr::pivot_longer(cols = !gene, 
                            values_to = "voom_count",
                            names_to = "pseudobulk_colnames") %>%
        dplyr::left_join(.,metadata)
      
      # slope of limma fit don't match corrected counts, doing fit on the fly
      
      p = ggplot(gene_to_plot, aes(x=APOE_sum_scaled, y=voom_count)) +
        geom_point( size=3,alpha=0.6) + 
        geom_smooth(method = "lm") +
        ggpubr::stat_cor(method="pearson") +
        theme_bw() + 
        ggtitle( g$symbol[i]) + 
        ylab("Normalised log2(counts+0.5)")
      plot(p)
    }
    
    dev.off()
  }
  #not really continuous, so intermediate risk categories don't seem to be useful in lm context
}

# Polygenic component of PRS only (no APOE) #############
#######################################################


for(treat in c("untreated","IFN","LPS")){
  pseudobulk =readr::read_tsv("../../data/results/5.prepare_differential_expression_pseudobulks/combined_pseudobulk_raw_counts_treatmentxprolifxidxpool.txt") %>%
    dplyr::filter(gene %in% gene_info$gene) %>%
    tibble::column_to_rownames(var = "gene") 
  metadata = readr::read_tsv("../../data/results/5.prepare_differential_expression_pseudobulks/metadata_treatmentxprolifxidxpool.txt") %>%
    tibble::column_to_rownames(var = "cols_names") %>%
    dplyr::rename(ncells = count,donor= donor_id) %>%
    dplyr::mutate(pseudobulk_colnames=paste(treatment,proliferation_status,donor,pool,sep = "_"),
                  donor = case_when(donor=="Arlene-003" ~ "Arlene",
                                    donor=="Cindy-005" ~ "Cindy",
                                    donor=="Dexter-006" ~ "Dexter",
                                    donor=="Fiona-010" ~ "Fiona",
                                    donor=="Gilma-009" ~ "Gilma",
                                    donor=="Hector-011" ~ "Hector",
                                    donor=="Imani-012" ~ "Imani",
                                    donor=="Javier-013" ~ "Javier",
                                    donor=="Keoni-014" ~ "Keoni",
                                    donor=="Olaf-018" ~ "Olaf",
                                    donor=="Bertha-004" ~ "Bertha",
                                    donor=="Mindy-016" ~ "Mindy",
                                    donor=="Qiana-022" ~ "Qiana",
                                    donor=="Nestor-017" ~ "Nestor",
                                    
                                    .default = donor))
  
  # subset to treatment and non-proliferating only
  metadata = metadata %>%
    dplyr::filter(treatment == treat) %>%
    dplyr::filter(proliferation_status == "Not_proliferating") %>%
    dplyr::rename(line = donor) 
  
  
  length(unique(metadata$line))
  
  
  
  # load sex and age information
  add_metadata = read.csv("../../data/metadata_info_hipsci_IPMAR.csv") %>%
    dplyr::mutate(line = case_when(donor=="Arlene-003" ~ "Arlene",
                                   donor=="Cindy-005" ~ "Cindy",
                                   donor=="Dexter-006" ~ "Dexter",
                                   donor=="Fiona-010" ~ "Fiona",
                                   donor=="Gilma-009" ~ "Gilma",
                                   donor=="Hector-011" ~ "Hector",
                                   donor=="Imani-012" ~ "Imani",
                                   donor=="Javier-013" ~ "Javier",
                                   donor=="Keoni-014" ~ "Keoni",
                                   donor=="Olaf-018" ~ "Olaf",
                                   donor=="Bertha-004" ~ "Bertha",
                                   donor=="Mindy-016" ~ "Mindy",
                                   donor=="Qiana-022" ~ "Qiana",
                                   donor=="Nestor-017" ~ "Nestor",
                                   .default = line)) %>%
    dplyr::select(line,Sex,Age)
  
  # converting NA to unknown to avoid problems with missing values
  metadata = metadata %>%
    dplyr::left_join(.,add_metadata) %>%
    dplyr::mutate(Age = case_when(Age == "" ~ 0,  # relevel to numeric
                                  is.na(Age) ~ 0,
                                  Age == "25-29" ~ 1,
                                  Age == "30-34" ~ 2,
                                  Age =="35-39" ~ 3,
                                  Age =="40-44" ~ 4,
                                  Age =="45-49" ~ 5,
                                  Age =="50-54" ~ 6,
                                  Age =="55-59" ~ 7,
                                  Age =="60-64" ~ 8,
                                  Age =="65-69" ~ 9,
                                  Age =="70-74" ~ 10,
                                  Age =="75-79" ~ 11,
                                  Age == "80"  ~ 12,
                                  Age == "89"  ~ 13,
                                  Age == "95"  ~ 15,
                                  .default = 0),
                  Sex = case_when(Sex == "" ~ "unkown",
                                  is.na(Sex) ~ "unknown",
                                  .default = Sex)) %>%
    distinct()
  
  pseudobulk = pseudobulk %>%
    dplyr::select(match(metadata$pseudobulk_colnames,colnames(.)))
  # same order
  identical(colnames(pseudobulk),metadata$pseudobulk_colnames)
  # load PRS classification
  #################
  prs = read_tsv("../../../hipsci_genotype_processing/data/prs_ad_bellenguez/PRSice/hipsci_polygenic_score_AD_Bellenguez_withAPOE.tsv") %>%
    dplyr::select(-FID)
  
  metadata = metadata %>%
    left_join(prs)
  
  ##### filter by number of cells
  
  discarded = metadata$ncells < 100
  table(discarded)
  
  pseudobulk = pseudobulk[,!discarded]
  metadata = metadata[!discarded,]
  length(unique(metadata$line)) # 190/250, 76% for untreated
  
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
  
  metadata$Sex = factor(metadata$Sex)
  metadata$Sex = relevel(metadata$Sex,ref = "Female")
  
  
  
  # subset v to the remaining elements of the metadata
  dge=dge[,metadata$pseudobulk_colnames]
  
  
  mat = model.matrix(~ prs_scaled + Sex + log10(ncells),metadata)
  mat = as.data.frame(mat)
  qr(mat)$rank
  ncol(mat)
  is.fullrank(mat)
  colnames(mat)
  mat = mat[,colSums(mat)!=0]
  colnames(mat)
  is.fullrank(mat)
  
  colnames(mat) = make.names(colnames(mat))
  
  
  v = voom(dge, design = mat, plot=TRUE)
  
  cor = duplicateCorrelation(v, design=mat, block = metadata$pool)
  cor$consensus # within-pool correlation
  # voom weights may have changed:
  v2 = voom(dge, design = mat, plot=TRUE, block = metadata$pool, correlation = cor$consensus)
  cor = duplicateCorrelation(v2, design = mat, block =  metadata$pool) # extract duplicates again
  cor$consensus
  
  # fit 
  fit = lmFit(v2, design=mat,block = metadata$pool, correlation  = cor$consensus)
  
  cont = c( "PRS_polygenic" = "prs_scaled") # changes in gene expression per unit increase in PRS, per treatment 
  
  
  contrast.matrix = makeContrasts(contrasts = cont,
                                  levels=mat)
  colnames(contrast.matrix) = c("PRS_polygenic")
  
  fit2 = contrasts.fit(fit, contrast.matrix)
  fit2 = eBayes(fit2)
  
  summary_additive = summary(decideTests(fit2,lfc = 0, # effects seem to be small
                                         adjust.method = "fdr",p.value=0.05))
  print(summary_additive)
  
  res_additive = list()
  for(cont in colnames(contrast.matrix)){
    res_additive[[cont]] = topTable(fit2, coef=cont,
                                    number = length(fit2$coefficients)) %>%
      dplyr::arrange(adj.P.Val) %>%
      tibble::rownames_to_column(var = "symbol")
    
    write.csv(res_additive[[cont]],paste0(outdir,
                                          "/DiffExpr_prs_scaled_",cont,"_",treat,".csv"),
              quote = F, sep = "\t",row.names = F, col.names = T)
  }
  saveRDS(res_additive,paste0(outdir,"/additive_res_prs_scaled_",treat,".rds"))
  
  ### remove pool batch effects for plotting
  # similar to partial residuals
  pseudobulk_corrected = limma::removeBatchEffect(v$E, batch = metadata$pool, design = mat)
  
  for(cont in colnames(contrast.matrix)){
    g = res_additive[[cont]] %>%
      dplyr::slice_head(n=30) 
    pdf(paste0(outdir,"/DiffExpr_prs_scaled_",cont,"_",treat,"_top30_sign.pdf"),width = 7,height = 7)
    message("Plotting top 30 genes sorted by abs(log2FC)")
    
    for(i in 1:30){
      expr = as.data.frame(pseudobulk_corrected)
      expr$gene = rownames(expr)
      
      gene_to_plot = expr %>%
        dplyr::filter(gene==g$symbol[i]) %>%
        tidyr::pivot_longer(cols = !gene, 
                            values_to = "voom_count",
                            names_to = "pseudobulk_colnames") %>%
        dplyr::left_join(.,metadata)
      
      # slope of limma fit don't match corrected counts, doing fit on the fly
      
      p = ggplot(gene_to_plot, aes(x=prs_scaled, y=voom_count)) +
        geom_point( size=3,alpha=0.6) + 
        geom_smooth(method = "lm") +
        ggpubr::stat_cor(method="pearson") +
        theme_bw() + 
        ggtitle( g$symbol[i]) + 
        ylab("Normalised log2(counts+0.5)")
      plot(p)
    }
    
    dev.off()
  }
  #not really continuous, so intermediate risk categories don't seem to be useful in lm context
}

########################################
### are sex and age confounded with PRS?
# are genotype PCs or expression PCs correlated with age and sex and/or with PRS?

gene_info = read_csv("/lustre/scratch123/hgi/teams/trynka/resources/biomart/Homo_sapiens.GRCh38.111.genes.csv") %>%
  dplyr::rename(gene = gene_name) %>%
  dplyr::filter(gene_biotype == "protein_coding")
### treatment with prs interaction
pseudobulk =readr::read_tsv("../../data/results/5.prepare_differential_expression_pseudobulks/combined_pseudobulk_raw_counts_treatmentxprolifxidxpool.txt") %>%
  dplyr::filter(gene %in% gene_info$gene) %>%
  tibble::column_to_rownames(var = "gene") 

metadata = readr::read_tsv("../../data/results/5.prepare_differential_expression_pseudobulks/metadata_treatmentxprolifxidxpool.txt") %>%
  tibble::column_to_rownames(var = "cols_names") %>%
  dplyr::rename(ncells = count,donor= donor_id) %>%
  dplyr::mutate(pseudobulk_colnames=paste(treatment,proliferation_status,donor,pool,sep = "_"),
                donor = case_when(donor=="Arlene-003" ~ "Arlene",
                                  donor=="Cindy-005" ~ "Cindy",
                                  donor=="Dexter-006" ~ "Dexter",
                                  donor=="Fiona-010" ~ "Fiona",
                                  donor=="Gilma-009" ~ "Gilma",
                                  donor=="Hector-011" ~ "Hector",
                                  donor=="Imani-012" ~ "Imani",
                                  donor=="Javier-013" ~ "Javier",
                                  donor=="Keoni-014" ~ "Keoni",
                                  donor=="Olaf-018" ~ "Olaf",
                                  donor=="Bertha-004" ~ "Bertha",
                                  donor=="Mindy-016" ~ "Mindy",
                                  donor=="Qiana-022" ~ "Qiana",
                                  donor=="Nestor-017" ~ "Nestor",
                                  .default = donor))

# subset to non-proliferating only
metadata = metadata %>%
  dplyr::filter(proliferation_status == "Not_proliferating") %>%
  dplyr::rename(line = donor) 


length(unique(metadata$line)) # 250



# load sex and age information
add_metadata = read.csv("../../data/metadata_info_hipsci_IPMAR.csv") %>%
  dplyr::mutate(line = case_when(donor=="Arlene-003" ~ "Arlene",
                                 donor=="Cindy-005" ~ "Cindy",
                                 donor=="Dexter-006" ~ "Dexter",
                                 donor=="Fiona-010" ~ "Fiona",
                                 donor=="Gilma-009" ~ "Gilma",
                                 donor=="Hector-011" ~ "Hector",
                                 donor=="Imani-012" ~ "Imani",
                                 donor=="Javier-013" ~ "Javier",
                                 donor=="Keoni-014" ~ "Keoni",
                                 donor=="Olaf-018" ~ "Olaf",
                                 donor=="Bertha-004" ~ "Bertha",
                                 donor=="Mindy-016" ~ "Mindy",
                                 donor=="Qiana-022" ~ "Qiana",
                                 donor=="Nestor-017" ~ "Nestor",
                                 
                                 
                                 .default = line)) %>%
  dplyr::select(line,Sex,Age)

# converting NA to unknown to avoid problems with missing values
metadata = metadata %>%
  dplyr::left_join(.,add_metadata) %>%
  dplyr::mutate(Age = case_when(Age == "" ~ 0,  # relevel to numeric
                                is.na(Age) ~ 0,
                                Age == "25-29" ~ 1,
                                Age == "30-34" ~ 2,
                                Age =="35-39" ~ 3,
                                Age =="40-44" ~ 4,
                                Age =="45-49" ~ 5,
                                Age =="50-54" ~ 6,
                                Age =="55-59" ~ 7,
                                Age =="60-64" ~ 8,
                                Age =="65-69" ~ 9,
                                Age =="70-74" ~ 10,
                                Age =="75-79" ~ 11,
                                Age == "80"  ~ 12,
                                Age == "89"  ~ 13,
                                Age == "95"  ~ 15,
                                .default = 0),
                Sex = case_when(Sex == "" ~ "unkown",
                                is.na(Sex) ~ "unknown",
                                .default = Sex)) %>%
  distinct()
# Age will be partly confounded with PRS and sex most likely

pseudobulk = pseudobulk %>%
  dplyr::select(match(metadata$pseudobulk_colnames,colnames(.)))
# same order
identical(colnames(pseudobulk),metadata$pseudobulk_colnames)
# load PRS classification
#################
prs = read_tsv("../../../hipsci_genotype_processing/data/prs_ad_bellenguez/PRSice/hipsci_polygenic_score_AD_Bellenguez_withAPOE.tsv") %>%
  dplyr::select(-FID)

metadata = metadata %>%
  left_join(prs)
##### filter by number of cells


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
metadata$Sex = factor(metadata$Sex)
metadata$Sex = relevel(metadata$Sex,ref = "Female")


mat = model.matrix(~ 0 + treatment   +  Sex + Age + pool + log10(ncells),metadata)
mat = as.data.frame(mat)
colnames(mat) = make.names(colnames(mat))
v = voom(dge, design = mat, plot=TRUE)

#inspecting with PCAtools
met = metadata
rownames(met) = metadata$pseudobulk_colnames
p = pca(v$E, metadata = met, removeVar = 0.1)

# categoricals need to be factors

eigencorplot(p,metavars = c('ncells',"full_PRS",'prs_scaled','APOE_sum_scaled',"treatment","Age","Sex","PC1","PC2"))


### variance partition
met$log10_ncells = log10(met$ncells)
form = ~  (1 | treatment)   +  (1 | Sex) + Age + (1 | pool) + (1 | line) + full_PRS + log10_ncells

varPart = fitExtractVarPartModel(v, form, met)
vp = sortCols(varPart)
plotVarPart(vp)

res = fitVarPartModel(v[1:1000, ], form, met)

colinearityScore(res[[2]]) # only deals with fixed terms

form = ~  treatment   +  Sex + Age + pool + log10_ncells +  full_PRS + prs_scaled + APOE_sum_scaled + PC1 + PC2

# Compute Canonical Correlation Analysis (CCA)
# between all pairs of variables
# returns absolute correlation value
C = canCorPairs(form, met)

# Plot correlation matrix
# between all pairs of variables
plotCorrMatrix(C)


### check variation for TREM2
i = which(rownames(varPart) == "TREM2")
GE = data.frame(
  Expression = v$E[i, ],
  Line = met$line
)

label = paste("Line:", format(varPart$line[i] * 100,
                              digits = 3
), "%")
main = rownames(v)[i]
plotStratify(Expression ~ Line, GE,
             colorBy = NULL,
             text = label, main = main
)

### per treatment
GE = data.frame(
  Expression = v$E[i, ],
  Treatment = met$treatment
)

label = paste("Treatment:", format(varPart$treatment[i] * 100,
                                   digits = 3
), "%")
main = rownames(v)[i]

plotStratify(Expression ~ Treatment, GE,text = label, main = main)


### check variance partition within treatment to be sure that sex explains more variance than age
# things can change depening how many variables you fit

form = ~    (1 | Sex) + Age + (1 | pool) + (1 | line) + full_PRS + log10_ncells

p = list()
for (treat in c("untreated","IFN","LPS")){
  subset = which(met$treatment == treat)
  vs = v
  vs$targets = v$targets[subset,]
  vs$design = v$design[subset,]
  vs$E = v$E[,subset]
  vs$weights = v$weights[,subset]
  mets = met[subset,]
  varPart = fitExtractVarPartModel(vs, form, mets)
  vp = sortCols(varPart)
  p[[treat]] = plotVarPart(vp,main = treat)
}

pdf(paste0(outdir,"/variance_partition_",treat,".pdf"),width = 5,height = 4)
patchwork::wrap_plots(p)

dev.off()
#### correct for pool then plot again


mat = model.matrix(~ 0 + treatment   +  log10(ncells),metadata)
mat = as.data.frame(mat)

is.fullrank(mat)
colnames(mat) = make.names(colnames(mat))


v = voom(dge, design = mat, plot=TRUE)

cor = duplicateCorrelation(v, design=mat, block = metadata$pool)
cor$consensus # within-pool correlation
# voom weights may have changed:
v2 = voom(dge, design = mat, plot=TRUE, block = metadata$pool, correlation = cor$consensus)
cor = duplicateCorrelation(v2, design = mat, block =  metadata$pool) # extract duplicates again
cor$consensus


mat = model.matrix(~ 1,data = metadata) # fit intercept

fit = lmFit(v2, design=mat,block = metadata$pool, correlation  = cor$consensus)
res = residuals(fit, v2)

# fit model on residuals
form = ~   (1 | treatment) +  (1 | Sex) + Age +  (1 | line) + (1 | pool) + full_PRS + log10_ncells

varPart = fitExtractVarPartModel(res, form, met)
plotVarPart(varPart,main = "After correcting for pool")


form = ~    (1 | Sex) + Age + (1 | pool) + (1 | line) + full_PRS + log10_ncells

p = list()
for (treat in c("untreated","IFN","LPS")){
  subset = which(met$treatment == treat)
  vs = res[,subset]
  mets = met[subset,]
  varPart = fitExtractVarPartModel(vs, form, mets)
  vp = sortCols(varPart)
  p[[treat]] = plotVarPart(vp,main = treat)
}

pdf(paste0(outdir,"/variance_partition_pool_corrected.pdf"),width = 5,height = 4)

patchwork::wrap_plots(p)
dev.off()

### ffect is not very strong
