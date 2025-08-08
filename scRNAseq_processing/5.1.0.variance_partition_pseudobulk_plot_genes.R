# variance partition and plotting specific genes

library(limma)
library(edgeR)
library(tidyverse) # For ggplot2 and easy manipulation of data
library(patchwork) # To combine plots
library(PCAtools)
library(variancePartition)
library(readr)
library(ggrepel)
set.seed(123)
source("./helpers.R")
outdir="../../data/results/5.1.0.variance_partition_analysis"
dir.create(outdir, recursive = T)


gene_info = read_csv("/lustre/scratch123/hgi/teams/trynka/resources/biomart/Homo_sapiens.GRCh38.111.genes.csv") %>%
  dplyr::rename(gene = gene_name) %>%
  dplyr::filter(gene_biotype == "protein_coding")

###########
### read in pseudobulk
pseudobulk =read_tsv("../../data/results/5.prepare_differential_expression_pseudobulks/combined_pseudobulk_raw_counts_treatmentxprolifxidxpool.txt") %>%
  dplyr::filter(gene %in% gene_info$gene) %>%
  tibble::column_to_rownames(var = "gene") 

metadata = read_tsv("../../data/results/5.prepare_differential_expression_pseudobulks/metadata_treatmentxprolifxidxpool.txt") %>%
  tibble::column_to_rownames(var = "cols_names") %>%
  dplyr::rename(ncells = count,line= donor_id) %>%
  dplyr::mutate(pseudobulk_colnames=paste(treatment,proliferation_status,line,pool,sep = "_"),
                donor = case_when(line=="Arlene-003" ~ "Arlene",
                                  line=="Cindy-005" ~ "Cindy",
                                  line=="Dexter-006" ~ "Dexter",
                                  line=="Fiona-010" ~ "Fiona",
                                  line=="Gilma-009" ~ "Gilma",
                                  line=="Hector-011" ~ "Hector",
                                  line=="Imani-012" ~ "Imani",
                                  line=="Javier-013" ~ "Javier",
                                  line=="Keoni-014" ~ "Keoni",
                                  line=="Olaf-018" ~ "Olaf",
                                  line=="Bertha-004" ~ "Bertha",
                                  line=="Mindy-016" ~ "Mindy",
                                  line=="Qiana-022" ~ "Qiana",
                                  line=="Nestor-017" ~ "Nestor",
                                  .default =  str_split_i(line,"_",i=1)))

# subset to non-proliferating only
metadata = metadata %>%
  dplyr::filter(proliferation_status == "Not_proliferating")


length(unique(metadata$line)) # 250



# load sex and age information
add_metadata = read.csv("../../data/donor_metadata_complete_with_imputed_sex.csv") %>%
  dplyr::select(donor,Sex,Age,Disease.Status)

summary(add_metadata)


metadata = metadata %>%
  dplyr::left_join(.,add_metadata) %>%
  dplyr::mutate(Age = case_when(Age == "" ~ NA,  # relevel to numeric
                                is.na(Age) ~ NA,
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
                                .default = NA),
                Sex = case_when(Sex == "" ~ "NA",
                                is.na(Sex) ~ "NA",
                                .default = Sex)) %>%
  distinct()

summary(add_metadata)

# Age will be partly confounded with PRS and sex most likely

pseudobulk = pseudobulk %>%
  dplyr::select(match(metadata$pseudobulk_colnames,colnames(.)))
# same order
identical(colnames(pseudobulk),metadata$pseudobulk_colnames)
# load PRS classification
#################
prs = read_tsv("../../../hipsci_genotype_processing/data/prs_ad_bellenguez/PRSice/hipsci_polygenic_score_AD_Bellenguez_withAPOE.tsv") %>%
  dplyr::select(-FID) %>%
  dplyr::mutate(donor =  str_split_i(line,"_",i=1))


metadata = metadata %>%
  dplyr::left_join(prs)

summary(metadata)

# removing NAs
metadata = metadata %>%
  tidyr::drop_na()
length(unique(metadata$line)) # 185/250

select = metadata$pseudobulk_colnames

pseudobulk = pseudobulk[,select]

########## to-do ##############
## adding scaled proliferation


### adding phagocytosis

########## end to-do ##############


### filters

discarded = metadata$ncells < 100
table(discarded)

pseudobulk = pseudobulk[,!discarded]
metadata = metadata[!discarded,]
length(unique(metadata$line)) # 60%


### variance partition

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

## including covariates with large effects, and also pool as random as in final diff expr
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

### full PRS x 

# Define formula
form = ~ pool +  treatment + Sex +  Age + Disease.Status + log10(ncells) + full_PRS + prs_scaled + APOE_sum_scaled 

C = canCorPairs(form, metadata)

# Plot correlation matrix
# between all pairs of variables
pdf(paste0(outdir,"/variancePart_corrmatrix_incl_treatment.pdf"), width = 10, height = 10)
plotCorrMatrix(C)
dev.off()

# variancePartition seamlessly deals with the result of voom()
# by default, it seamlessly models the precision weights
# This can be turned off with useWeights=FALSE

#form = ~ (1 | pool) + (1 | treatment) + Age + (1 | Disease.Status) + log10(ncells) + full_PRS  # singular fit
form = ~ (1 | pool) + (1 | treatment) +  Age   + full_PRS 

varPart = fitExtractVarPartModel(v2["TREM2",], form, metadata)

vp = sortCols(varPart)

plotVarPart(vp)


# Figure 2b
# plot expression stratified by Tissue
label <- paste("Individual:", format(varPart$Individual[i] * 100,
                                     digits = 3
), "%")
main <- rownames(geneExpr)[i]
plotStratify(Expression ~ Individual, GE,
             colorBy = NULL,
             text = label, main = main)
             
# Figure 1a
# Bar plot of variance fractions for the first 10 genes
plotPercentBars(vp[1, ])


colinearityScore(vp[[1]])

## heatmap ###

# heatmap with voom counts
voom_cor = cor(v2$E)
# Plot heatmap

pools = c(rainbow(16))
names(pools) = paste0("pool",2:17)
ann_color = list("treatment"=c("untreated"="black","LPS"="yellow","IFN"="red"),
                 "pool"=pools)
annot = metadata[, c("treatment", "pool","ncells"), drop=F]
rownames(annot) = metadata$pseudobulk_colnames
png(paste0(outdir,"/cor_heatmap_voom_pool_corrected_allpools.png"),
    width = 16, height = 14, res = 400,units = "in", type = "cairo")
voom_cor %>%
  magrittr::set_rownames(metadata$pseudobulk_colnames) %>%
  pheatmap::pheatmap(., 
                     annotation = annot,
           show_colnames = FALSE,
           show_rownames = FALSE,
           annotation_colors = ann_color)

dev.off()

# subset to outlier correlations (pool9, and pools 3,6 ,13 and 17 for reference)
tomatch = metadata %>%
  dplyr::filter(str_detect(pseudobulk_colnames,"pool3") | 
                  str_detect(pseudobulk_colnames,"pool6") | 
                  str_detect(pseudobulk_colnames,"pool17") | 
                  str_detect(pseudobulk_colnames,"pool13") | 
                  str_detect(pseudobulk_colnames,"pool9") ) %>%
  dplyr::select(pseudobulk_colnames) %>%
  unlist() %>%
  unname(.)
voom_cor_subset = voom_cor[tomatch,tomatch]
voom_cor_subset %>%
  magrittr::set_rownames(tomatch) %>%
  pheatmap::pheatmap(., 
                     annotation = annot,
                     show_colnames = FALSE,
                     show_rownames = FALSE,
                     annotation_colors = ann_color)

###################################
# do within treatment to quantify better variance from other covariates ###########

metadata = metadata %>%
  dplyr::filter(treatment == "untreated")

select = metadata$pseudobulk_colnames

pseudobulk = pseudobulk[,select]


discarded = metadata$ncells < 100
table(discarded)

pseudobulk = pseudobulk[,!discarded]
metadata = metadata[!discarded,]
length(unique(metadata$line)) # 145/250, 60%


### variance partition

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

## including covariates with large effects, and also pool as random as in final diff expr
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

### full PRS x 

# Define formula
form = ~ pool  + Sex +  Age + Disease.Status + log10(ncells) + full_PRS + prs_scaled + APOE_sum_scaled 

C = canCorPairs(form, metadata)

# Plot correlation matrix
# between all pairs of variables
pdf(paste0(outdir,"/variancePart_corrmatrix_untreated_only.pdf"), width = 10, height = 10)

plotCorrMatrix(C)

dev.off()
# variancePartition deals with the result of voom()
# by default, it  models the precision weights
# This can be turned off with useWeights=FALSE

form = ~ (1 | pool) + Age + (1 | Sex) + (1 | Disease.Status) + log10(ncells) + full_PRS 

varPart = fitExtractVarPartModel(v2, form, metadata)

vp = sortCols(varPart)
plotPercentBars(vp[c("TREM2","APH1B","SPG11"), ])
plotVarPart(vp)


# checking line effect, to see if the unexplained variance could come from here
form2 = ~ (1 | pool) + (1 | line) +  Age  + (1 |Sex) + full_PRS + log10(ncells)

varPart2 = fitExtractVarPartModel(v2, form2, metadata)

vp2 = sortCols(varPart2)

pdf(paste0(outdir,"/variancePart_variance_explained_untreated_only.pdf"), width = 10, height = 10)

plotVarPart(vp2)

dev.off()

C = canCorPairs(form2, metadata)
# High colinearity between variables:
#   line and Age
# line and Sex
# line and full_PRS
pdf(paste0(outdir,"/variancePart_corrmatrix_untreated_only_with_line.pdf"), width = 10, height = 10)

plotCorrMatrix(C)

dev.off()

# residuals still explain a lot of variance and line is correlated with so many things we adjust for
# we can't include it anyway

### examining the sequenced samples pseudobulks ##########

pseudobulk =read_tsv("../../data/results/5.prepare_differential_expression_pseudobulks/combined_pseudobulk_raw_counts_treatmentxprolifxidxseqsample.txt") %>%
  dplyr::filter(gene %in% gene_info$gene) %>%
  tibble::column_to_rownames(var = "gene") 

metadata = read_tsv("../../data/results/5.prepare_differential_expression_pseudobulks/metadata_treatmentxprolifxidxseqsample.txt") %>%
  tibble::column_to_rownames(var = "cols_names") %>%
  dplyr::rename(ncells = count,line= donor_id) %>%
  dplyr::mutate(pseudobulk_colnames=paste(treatment,proliferation_status,line,orig.ident,sep = "_"),
                donor = case_when(line=="Arlene-003" ~ "Arlene",
                                  line=="Cindy-005" ~ "Cindy",
                                  line=="Dexter-006" ~ "Dexter",
                                  line=="Fiona-010" ~ "Fiona",
                                  line=="Gilma-009" ~ "Gilma",
                                  line=="Hector-011" ~ "Hector",
                                  line=="Imani-012" ~ "Imani",
                                  line=="Javier-013" ~ "Javier",
                                  line=="Keoni-014" ~ "Keoni",
                                  line=="Olaf-018" ~ "Olaf",
                                  line=="Bertha-004" ~ "Bertha",
                                  line=="Mindy-016" ~ "Mindy",
                                  line=="Qiana-022" ~ "Qiana",
                                  line=="Nestor-017" ~ "Nestor",
                                  .default =  str_split_i(line,"_",i=1)))

# subset to non-proliferating only
metadata = metadata %>%
  dplyr::filter(proliferation_status == "Not_proliferating")


length(unique(metadata$line)) # 250



# load sex and age information
add_metadata = read.csv("../../data/donor_metadata_complete_with_imputed_sex.csv") %>%
  dplyr::select(donor,Sex,Age,Disease.Status)

summary(add_metadata)


metadata = metadata %>%
  dplyr::left_join(.,add_metadata) %>%
  dplyr::mutate(Age = case_when(Age == "" ~ NA,  # relevel to numeric
                                is.na(Age) ~ NA,
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
                                .default = NA),
                Sex = case_when(Sex == "" ~ "NA",
                                is.na(Sex) ~ "NA",
                                .default = Sex)) %>%
  distinct()

summary(add_metadata)

pseudobulk = pseudobulk %>%
  dplyr::select(match(metadata$pseudobulk_colnames,colnames(.)))
# same order
identical(colnames(pseudobulk),metadata$pseudobulk_colnames)

prs = read_tsv("../../../hipsci_genotype_processing/data/prs_ad_bellenguez/PRSice/hipsci_polygenic_score_AD_Bellenguez_withAPOE.tsv") %>%
  dplyr::select(-FID) %>%
  dplyr::mutate(donor =  str_split_i(line,"_",i=1)) %>%
  dplyr::mutate(line = case_when(donor=="Arlene" ~ "Arlene-003",
                                          donor=="Cindy" ~ "Cindy-005",
                                          donor=="Dexter"  ~ "Dexter-006",
                                          donor== "Fiona" ~ "Fiona-010",
                                          donor=="Gilma" ~ "Gilma-009" ,
                                          donor=="Hector" ~ "Hector-011" ,
                                          donor== "Imani" ~ "Imani-012"  ,
                                          donor=="Javier" ~ "Javier-013",
                                          donor== "Keoni" ~"Keoni-014" ,
                                          donor== "Olaf" ~ "Olaf-018" ,
                                          donor== "Bertha" ~ "Bertha-004" ,
                                          donor==  "Mindy" ~ "Mindy-016",
                                          donor=="Qiana" ~ "Qiana-022" ,
                                          donor== "Nestor" ~ "Nestor-017",
                                 .default = line)) %>%
  dplyr::select(-donor)
  

metadata = metadata %>%
  dplyr::left_join(prs,by = "line")

summary(metadata)

# NOT removing NAs
length(unique(metadata$line)) # 250

select = metadata$pseudobulk_colnames

pseudobulk = pseudobulk[,select]

### filters

discarded = metadata$ncells < 100
table(discarded)

pseudobulk = pseudobulk[,!discarded]
metadata = metadata[!discarded,]
length(unique(metadata$line)) # 182/250, 73%


### variance partition

dge = edgeR::DGEList(counts=pseudobulk)



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

## including covariates with large effects, and also pool as random as in final diff expr
mat = model.matrix(~ 0 + treatment   +  log10(ncells) + Sex,metadata)
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

# heatmap with voom counts
# heatmap with voom counts
voom_cor = cor(v2$E)
# Plot heatmap

pools = c(rainbow(16))
names(pools) = paste0("pool",2:17)
ann_color = list("treatment"=c("untreated"="black","LPS"="yellow","IFN"="red"),
                 "pool"=pools)
annot = metadata[, c("treatment", "pool","ncells"), drop=F]
rownames(annot) = metadata$pseudobulk_colnames
png(paste0(outdir,"/cor_heatmap_voom_pool_corrected_allpools_seqsamples.png"),
    width = 16, height = 14, res = 400,units = "in", type = "cairo")
voom_cor %>%
  magrittr::set_rownames(metadata$pseudobulk_colnames) %>%
  pheatmap::pheatmap(., 
                     annotation = annot,
                     show_colnames = FALSE,
                     show_rownames = FALSE,
                     annotation_colors = ann_color)

dev.off()

# subset to outlier correlations (pool9, and pools 3,6 ,13 and 17 for reference)
tomatch = metadata %>%
  dplyr::filter(str_detect(pseudobulk_colnames,"pool3") | 
                  str_detect(pseudobulk_colnames,"pool6") | 
                  str_detect(pseudobulk_colnames,"pool17") | 
                  str_detect(pseudobulk_colnames,"pool13") | 
                  str_detect(pseudobulk_colnames,"pool9") ) %>%
  dplyr::select(pseudobulk_colnames) %>%
  unlist() %>%
  unname(.)
voom_cor_subset = voom_cor[tomatch,tomatch]
voom_cor_subset %>%
  magrittr::set_rownames(tomatch) %>%
  pheatmap::pheatmap(., 
                     annotation = annot,
                     show_colnames = FALSE,
                     show_rownames = FALSE,
                     annotation_colors = ann_color)
# pool9 as a whole seems to be an outlier

# exploring correlations
hclustering = hclust(dist(voom_cor,method = "euclidean"))
ctree = cutree(hclustering, k = 6) # number of second level clusters I see in heatmap
ctree = data.frame("cluster" = ctree,"pseudobulk_colnames" = names(ctree))

metadata_ctree = metadata %>%
  dplyr::left_join(ctree) %>%
  dplyr::select(line, treatment,pool,ncells,cluster, orig.ident, pseudobulk_colnames)


table(metadata_ctree$treatment, metadata_ctree$cluster) 
#             1   2   3   4   5   6
# untreated 433  82   0   0   0   0
# IFN         0   0   0   0 337  68
# LPS         0   0 333  72   0   0
table(metadata_ctree$pool, metadata_ctree$cluster)

#         1  2  3  4  5  6
# pool10 94  0 78  0 63  0
# pool11 29  0  9  0 15  0
# pool12  4  0  2  0  2  0
# pool13 37  0 19  0 19  0
# pool14 23  0 30  0 31  0
# pool15 29  0 30  0 30  0
# pool16 46  0 46  0 46  0
# pool17 59  0 57  0 61  0
# pool2  55  0  0  0  7  0
# pool3   6  0  5  0  7  0
# pool4   8  0  9  0  9  0
# pool5   9  0  9  0  9  0
# pool6   9  0  9  0  9  0
# pool7  12  0 14  0 14  0
# pool8  13  0 16  0 15  0
# pool9   0 82  0 72  0 68
  
# all donors in pool 9 are clustering away from the rest
metadata_ctree %>%
  write_csv(paste0(outdir,"/clustering_outliers_metadata_seqsample_results.csv"))


 p1 = metadata_ctree %>%
   dplyr::mutate(pool=factor(pool,levels = paste0("pool",2:17), ordered = TRUE)) %>%
   ggplot( aes(x=cluster,fill = treatment)) +
     geom_bar() +
     geom_text(stat='count', aes(label=..count..), vjust=-1) +
     theme_bw()+
    ggtitle("Number of sequencing samples (line x orig.ident) per cluster") + 
     scale_x_continuous(breaks = 1:6)+
     facet_wrap(vars(pool))

 p2 = metadata_ctree %>%
   dplyr::mutate(pool=factor(pool,levels = paste0("pool",2:17), ordered = TRUE)) %>%
   ggplot( aes(x=cluster, y = log(ncells),col = treatment)) +
   geom_point() +
   theme_bw()+
   scale_x_continuous(breaks = 1:6)+
   facet_wrap(vars(pool))
 
 png(paste0(outdir,"/barplot_cordist_allpools_seqsamples.png"),
           width = 10, height = 9, res = 400,units = "in", type = "cairo")
 plot(p1)
 dev.off()

 png(paste0(outdir,"/dotplot_ncells_cordist_allpools_seqsamples.png"),
          width = 15, height = 14, res = 400,units = "in", type = "cairo")
 plot(p2)
 dev.off()
 
 ### same after explicitly removing batch effects
 rbe = limma::removeBatchEffect(pseudobulk,
                                batch = metadata$pool,
                                group = metadata$treatment) 
 
 rbe_cor = cor(rbe)
 pools = c(rainbow(16))
 names(pools) = paste0("pool",2:17)
 ann_color = list("treatment"=c("untreated"="black","LPS"="yellow","IFN"="red"),
                  "pool"=pools)
 annot = metadata[,
                   c("treatment", "pool","ncells"), drop=F]
 rownames(annot) = metadata$pseudobulk_colnames
 png(paste0(outdir,"/cor_heatmap_remBatchEff_pool_corrected_allpools_seqsamples.png"),
     width = 16, height = 14, res = 400,units = "in", type = "cairo")
 rbe_cor %>%
   magrittr::set_rownames(metadata$pseudobulk_colnames) %>%
   pheatmap::pheatmap(., 
                      annotation = annot,
                      show_colnames = FALSE,
                      show_rownames = FALSE,
                      annotation_colors = ann_color)
 
 dev.off()
 
 hclustering = hclust(dist(rbe_cor,method = "euclidean"))
 ctree = cutree(hclustering, k = 2) # number of second level clusters I see in heatmap
 ctree = data.frame("cluster" = ctree,"pseudobulk_colnames" = names(ctree))
 
 metadata_ctree = metadata %>%
   dplyr::left_join(ctree) %>%
   dplyr::select(line, treatment,pool,ncells,cluster, orig.ident, pseudobulk_colnames)
 
 
 table(metadata_ctree$pool, metadata_ctree$cluster)
 # some specific donors from pool 12, 5 and 6 cluster away
 # still I don't get why it's not separating more clearly by treatment
 metadata_ctree %>%
   write_csv(paste0(outdir,"/clustering_outliers_metadata_seqsample_remBatchEff_results.csv"))
 
 # subset by treatment because even after grouping it doesn't respect treatment
 for(treatment in c("untreated","LPS","IFN")){
   pseudobulk2 = pseudobulk %>%
     dplyr::select(contains(treatment)) 
   metadata2 = metadata %>%
     dplyr::filter(pseudobulk_colnames %in% colnames(pseudobulk2)) 
   
   rbe = limma::removeBatchEffect(pseudobulk2,
                                  batch = metadata2$pool) 
   
   rbe_cor = cor(rbe)
   pools = c(rainbow(16))
   names(pools) = paste0("pool",2:17)
   ann_color = list("treatment"=c("untreated"="black","LPS"="yellow","IFN"="red"),
                    "pool"=pools)
   annot = metadata2[,
                    c("treatment", "pool","ncells"), drop=F]
   rownames(annot) = metadata2$pseudobulk_colnames
   png(paste0(outdir,"/",treatment,"_cor_heatmap_remBatchEff_pool_corrected_allpools_seqsamples.png"),
       width = 16, height = 14, res = 400,units = "in", type = "cairo")
   rbe_cor %>%
     magrittr::set_rownames(metadata2$pseudobulk_colnames) %>%
     pheatmap::pheatmap(., 
                        annotation = annot,
                        show_colnames = FALSE,
                        show_rownames = FALSE,
                        annotation_colors = ann_color)
   
   dev.off()
   
   
   hclustering = hclust(dist(rbe_cor,method = "euclidean"))
   ctree = cutree(hclustering, k = 2) # number of second level clusters I see in heatmap
   ctree = data.frame("cluster" = ctree,"pseudobulk_colnames" = names(ctree))
   
   metadata_ctree = metadata2 %>%
     dplyr::left_join(ctree) %>%
     dplyr::select(line, treatment,pool,ncells,cluster, orig.ident, pseudobulk_colnames)
   
   
   table(metadata_ctree$pool, metadata_ctree$cluster)
   # pool12, 5 and 6, incl lizq_1, zaie_1, qonc_1
 }

 
############## plot specific gene expression


# plot TREM2 expression by pool
# and other genes
 
transposed = v2$E %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% c("BIN1","TREM2","APH1B","TYROBP","DAP12")) %>%
  tidyr::pivot_longer(cols  = -gene,
                      names_to = "pseudobulk_colnames",
                      values_to = "voom_counts") %>%
  dplyr::left_join(metadata) %>%
  dplyr::mutate(pool = factor(pool,levels = paste0("pool",2:17)))
 
 png(paste0(outdir,"/TREM2_BIN1_APH1B_TYROBP_voom_pool_corrected_allpools_seqsamples.png"),
     width = 8, height = 6, res = 400,units = "in", type = "cairo")
transposed %>%
  ggplot(aes(x = pool,y = voom_counts, colour = treatment)) +
  geom_point() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5)) + 
  facet_wrap(vars(gene))
dev.off()

### scatterplot correlation TREM2 and TYROBP
png(paste0(outdir,"/TREM2_TYROBP_voom_line_scatterplot.png"),
    width = 8, height = 6, res = 400,units = "in", type = "cairo")
transposed %>%
  dplyr::filter(gene %in% c("TREM2","TYROBP") & proliferation_status == "Not_proliferating") %>%
  dplyr::group_by(gene, proliferation_status, line) %>%
  dplyr::summarise(voom_counts = mean(voom_counts)) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(id_cols = line,names_from = gene, values_from = voom_counts) %>%
  ggplot(aes(x = TREM2,y = TYROBP,label = line)) +
  geom_point() + 
  geom_label_repel() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5))

dev.off()

### batch corrected explicitly
transposed = rbe %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% c("BIN1","TREM2","APH1B")) %>%
  tidyr::pivot_longer(cols  = -gene,
                      names_to = "pseudobulk_colnames",
                      values_to = "rbe_counts") %>%
  dplyr::left_join(metadata) %>%
  dplyr::mutate(pool = factor(pool,levels = paste0("pool",2:17)))

png(paste0(outdir,"/TREM2_BIN1_APH1B_TYROBP_remBatchEff_pool_corrected_allpools_seqsamples.png"),
    width = 8, height = 6, res = 400,units = "in", type = "cairo")

transposed %>%
  ggplot(aes(x = pool,y = rbe_counts, colour = treatment)) +
  geom_point() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5)) + 
  facet_wrap(vars(gene))

dev.off()





