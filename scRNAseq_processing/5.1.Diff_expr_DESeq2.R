# Differential expression with DEseq2
.libPaths(c('/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/',
            "/software/teamtrynka/ma23/R4.1/libs",
            "/software/R-4.1.0/lib/R/library"
)) # path to libraries
library(DESeq2)
library(edgeR) # for some helper functions
library(tidyverse) # For ggplot2 and easy manipulation of data
library(patchwork) # To combine plots
library(vsn) # Some visualizations
library(AnnotationDbi) # gene annotation
library(org.Hs.eg.db) # gene annotation
library(biomaRt)
library(ggrepel)
library(pheatmap)
library(sva)
source("./helpers.R")
outdir="../../data/results/5.1.Diff_expr_DESeq2"
dir.create(outdir, recursive = T)

## check IFN-activated genes
# enrichment in all 3 treatment conditions
# check if cell cycle genes are very upregulated in certain treatments

######### main ####

# load data
pseudobulk =readr::read_tsv("../../data/results/5.prepare_differential_expression_pseudobulks/combined_pseudobulk_raw_counts_treatmentxprolifxidxpool.txt") %>%
  tibble::column_to_rownames(var = "gene")
metadata = readr::read_tsv("../../data/results/5.prepare_differential_expression_pseudobulks/metadata_treatmentxprolifxidxpool.txt") %>%
  tibble::column_to_rownames(var = "cols_names") %>%
  dplyr::rename(ncells = count)
metadata_rows = rownames(metadata)

## gene object
genes = data.frame(ensembl_gene_ids = rownames( pseudobulk),
                   external_gene_names = rownames( pseudobulk))
genes$ensembl_gene_ids = as.character(genes$ensembl_gene_ids)
genes$external_gene_names = as.character(genes$external_gene_names)
gene_names = symbol_to_ensembl(genes$external_gene_names)
gene_names = data.frame(external_gene_names=names(gene_names), 
                        ensembl_gene_ids=gene_names, row.names=NULL)

genes$ensembl_gene_ids = gene_names$ensembl_gene_ids
rm(gene_names)
rownames(genes) = genes$external_gene_names

# set untreated as baseline
metadata$treatment = factor(metadata$treatment)
metadata$treatment = relevel(metadata$treatment, "untreated")
metadata$pool = factor(metadata$pool)
metadata$pool = relevel(metadata$pool, "pool17")



# ensure they are sorted identically
order_df1 = match(rownames(metadata), colnames(pseudobulk))

# Sort both data frames based on the order of row names
metadata = metadata[order_df1, , drop = FALSE]
pseudobulk = pseudobulk[, order_df1, drop = FALSE]

# creating new grouping
metadata$prolif_x_treatment=paste0(metadata$proliferation_status,"_",metadata$treatment)
metadata$log10_ncells=log10(metadata$ncells)
dds = DESeqDataSetFromMatrix(countData = pseudobulk,
                             colData =  metadata,
                             #design = ~    donor_id + treatment)
                             design = ~    pool + treatment)
length(unique(dds$donor_id)) # 248 donors 
dds

hist(log10(dds$ncells),breaks = 50)
hist(dds$ncells,breaks = 50)
summary(dds$ncells)
quantile(dds$ncells,probs = seq(0,1,0.05))
## remove samples with very low number of cells (min 50 - around the 45th quantile) 
discarded = dds$ncells < 50
dds = dds[,!discarded]
summary(discarded)
length(unique(dds$donor_id)) # 199 donors
# 194/248 is ~80% donors retained

#########
# remove genes not expressed in 30% samples (min 1 CPM)
isexpr = rowSums(edgeR::cpm(dds) > 1) >= floor(0.3*ncol(dds))
message(paste0("There are ", sum(isexpr), " genes with over 1 count per million in 30% of samples"))
dds = dds[isexpr,]  
dim(dds)
genes = genes[genes$external_gene_names %in% rownames(counts(dds)),]
# Transform counts for data visualization
vsttr = vst(dds, blind=TRUE)

# Plot PCA

DESeq2::plotPCA(vsttr, intgroup = "prolif_x_treatment")
# Proliferating IFN is very similar to Not proliferating untreated??
DESeq2::plotPCA(vsttr, intgroup = "log10_ncells")

vst_mat = assay(vsttr)
vst_cor = cor(vst_mat)
# Plot heatmap
ann_color = list("prolif_x_treatment" = c("Not_proliferating_LPS" = "yellow","Proliferating_LPS" = "gold", 
                                          "Not_proliferating_IFN" = "coral","Proliferating_IFN" = "red",
                                          "Not_proliferating_untreated" = "grey","Proliferating_untreated" = "black"),
                 "treatment"=c("untreated"="black","LPS"="yellow","IFN"="red"),
                 "proliferation_status"=c("Not_proliferating"="lightblue","Proliferating"="darkblue"))

png(paste0(outdir,"/cor_heatmap_vst_blindTRUE_allclusters.png"),
    width = 16, height = 14, res = 400,units = "in", type = "cairo")
vst_cor %>%
  magrittr::set_rownames(colData(vsttr)$donor_id) %>%
  pheatmap(., annotation = metadata[, c("log10_ncells","proliferation_status","treatment",
                                        "prolif_x_treatment"), drop=F],
           show_colnames = FALSE,
           annotation_colors = ann_color)

dev.off()

#######################
# subset to non-proliferating only
metadata = metadata %>%
  dplyr::filter(proliferation_status == "Not_proliferating")

pseudobulk = pseudobulk %>%
  dplyr::select(match(rownames(metadata),colnames(.)))

dds = DESeqDataSetFromMatrix(countData = pseudobulk,
                             colData =  metadata,
                             design = ~    (1/ncells) + pool + treatment)
                             #design = ~    donor_id + treatment)
length(unique(dds$donor_id)) # 248 donors 
dds

hist(log10(dds$ncells),breaks = 50)
hist(dds$ncells,breaks = 50)
summary(dds$ncells)
quantile(dds$ncells,probs = seq(0,1,0.05))

## remove samples with very low number of cells (min 100- around 47th quantile now that we exclude proliferating) 
discarded = dds$ncells < 100
dds = dds[,!discarded]
summary(discarded)
length(unique(dds$donor_id)) # 189 donors
# 194/248 is ~76% donors retained

####################
# remove samples from pools 9 adb 10 ########
# discarded = dds$pool %in% c("pool9", "pool10")
# dds = dds[,!discarded]
# summary(discarded)
# length(unique(dds$donor_id)) # 154 donors

#######
# if fitting donors
# Remove donors that are present only in 2 or fewer levels to avoid rank problems
# doesn't seem to be needed
#donors_to_retain = colSums(table(dds$treatment,dds$donor_id))==3
# dds = dds[,dds$donor_id %in% names(donors_to_retain)[donors_to_retain]]
# dds$donor_id = relevel(dds$donor_id,ref = "hegp_3")
# dds$donor_id = droplevels(dds$donor_id)
length(unique(dds$donor_id)) # 143 donors, 58% if dropping, 194 if not dropping
table(dds$treatment,dds$donor_id)
##########


# remove genes not expressed in 30% samples (min 1 CPM)
isexpr = rowSums(edgeR::cpm(dds) > 1) >= floor(0.3*ncol(dds))
message(paste0("There are ", sum(isexpr), " genes with over 1 count per million in 30% of samples"))
dds = dds[isexpr,]  
dim(dds)
genes = genes[genes$external_gene_names %in% rownames(counts(dds)),]
colData(dds)$pool = droplevels(colData(dds)$pool)

# Transform counts for data visualization
# calculates a variance stabilizing transformation (VST) 
# from the fitted dispersion-mean relation(s) and then transforms 
# the count data (normalized by division by the size factors or 
# normalization factors), yielding a matrix of 
# values which are now approximately homoskedastic (having constant variance 
#  along the range of mean values). The transformation also normalizes with 
# respect to library size. 

##### now blind to design (expr ~ donor_id + treatment)
vsttr = vst(dds, blind=TRUE)

# Plot PCA

DESeq2::plotPCA(vsttr, intgroup = "treatment")
# Proliferating IFN is very similar to Not proliferating untreated??
DESeq2::plotPCA(vsttr, intgroup = "log10_ncells")

vst_mat = assay(vsttr)
vst_cor = cor(vst_mat, method = "pearson")
# Plot heatmap
ann_color = list("prolif_x_treatment" = c("Not_proliferating_LPS" = "yellow","Proliferating_LPS" = "gold", 
                                          "Not_proliferating_IFN" = "coral","Proliferating_IFN" = "red",
                                          "Not_proliferating_untreated" = "grey","Proliferating_untreated" = "black"),
                 "treatment"=c("untreated"="black","LPS"="yellow","IFN"="red"),
                 "proliferation_status"=c("Not_proliferating"="lightblue","Proliferating"="darkblue"))

png(paste0(outdir,"/cor_heatmap_vst_blindTRUE.png"),
    width = 16, height = 14, res = 400,units = "in", type = "cairo")
vst_cor %>%
  magrittr::set_rownames(colData(vsttr)$donor_id) %>%
  pheatmap(., annotation = metadata[, c("log10_ncells","pool","treatment"), drop=F],
           show_colnames = FALSE,show_rownames = FALSE,
           annotation_colors = ann_color)

dev.off()
# 
# mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'https://www.ensembl.org')
# genes_bio = biomaRt::getBM(attributes = c("ensembl_gene_id", "transcript_biotype"), filters = c("ensembl_gene_id"),
#                            values = list(genes$ensembl_gene_ids), mart = mart)
# genes_bio = genes_bio[genes_bio$transcript_biotype == "protein_coding",]
# genes = genes[genes$ensembl_gene_ids %in% genes_bio$ensembl_gene_id,]
# dds = dds[rownames(dds) %in% genes$external_gene_names, ]  
# dim(dds)



png(paste0(outdir,"/sc_RNA_allclusters_library_size_after_filter.png"),
    width = 6, height = 4, res = 400,units = "in", type = "cairo")
barplot(colSums(assay(dds))*1e-6, names = dds$orig.ident, 
        ylab="Library size (millions)", las = 2, cex.names = 0.4)
dev.off()



## PCA

plots_pcas = plot_pca(vst_mat, colData(dds))

png(paste0(outdir,"/PCA_treatment.png"),
    width = 9, height = 7, res = 400,units = "in", type = "cairo")
plots_pcas[[1]]
dev.off()

png(paste0(outdir,"/PCA_count.png"),
    width = 9, height = 7, res = 400,units = "in", type = "cairo")
plots_pcas[[2]]
dev.off()

png(paste0(outdir,"/PCA_count_PC2_3.png"),
    width = 9, height = 7, res = 400,units = "in", type = "cairo")
plots_pcas[[3]]
dev.off()

png(paste0(outdir,"/PCA_pool_PC2_3.png"),
    width = 9, height = 7, res = 400,units = "in", type = "cairo")
plots_pcas[[4]]
dev.off()

# pcatool = PCAtools::pca(vst_mat, metadata = colData(dds), removeVar = 0.1)
# 
# PCAtools::eigencorplot(pcatool)
getLoadings = function(dds){
  
  df = assay(dds)
  pca = prcomp(t(df), retx = TRUE)
  
  return(pca$rotation)
}

loadings_vst = getLoadings(vsttr) %>% as.data.frame()

# show the top 10 genes from PC3
loadings_PC3_top100 = loadings_vst %>% 
  dplyr::mutate(symbol = rownames(.)) %>%
  # select only the PCs we are interested in
  dplyr::select(symbol, PC3) %>%
  # convert to "long" format
  pivot_longer(cols = "PC3", names_to = "PC3", values_to = "loadings") %>% 
  # for PC2
  group_by(PC3) %>% 
  # arrange by descending order
  arrange(desc(abs(loadings))) %>%
  slice(1:100)

loadings_PC3_top100 # many pseudogenes? do GO/pathways

## get PC3 samples

pca1 = prcomp(t(vst_mat), retx = TRUE)

plot(pca1, type = "l") #variance vs first 10 components
summary(pca1)          #importance of each component (important line is "proportion of variance")

percentVar = (pca1$sdev)^2 / sum(pca1$sdev^2)
percentVar = round(100 * percentVar)
pcs = as.data.frame(pca1$x)
pcs = cbind(pcs, colData(dds))

donors_pc3_over0 = pcs %>%
  dplyr::filter(PC3>0)
unique(donors_pc3_over0$donor_id)

donors_pc3_over0 = donors_pc3_over0 %>%
  dplyr::left_join(donor_info)

table(donors_pc3_over0$pool)
#### Diffexpr
message("Performing differential expression for all levels")

dds = estimateSizeFactors(dds)

colData(dds)$pool = droplevels(colData(dds)$pool)
colData(dds)$donor_id = droplevels(colData(dds)$donor_id)

check = colData(dds)
table(check$donor_id, check$treatment)
mat = model.matrix(~ donor_id + treatment, colData(dds))
is.fullrank(mat)
mat = model.matrix(~ pool + treatment, colData(dds))
is.fullrank(mat)

BiocParallel::register(BiocParallel::MulticoreParam(workers = 32))

dds = DESeq(dds,parallel=TRUE)


png(paste0(outdir,"/disp_estimates_allpools_pool_treatment_design.png"),
    width = 9, height = 7, res = 400,units = "in", type = "cairo")
plotDispEsts(dds) # weird dispersion, should decrease with more counts and follow the mean
# https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html
# a bit better after removing pools 9 and 10, but still a straight line
# also mostly straight if not excluding those pools
dev.off()
resultsNames(dds)

for(treatment in c("IFN","LPS")){
  message("Getting contrasts for ", treatment," treatment")
  
  res_treatment = results(dds, 
                          name = paste0("treatment_",treatment,"_vs_untreated"),
                          lfcThreshold = 1)
  message("Number of padj<0.01 and abs(log2FC>1)")
  message(table(res_treatment$padj<0.01)) # 676
  identical(rownames(res_treatment),rownames(assay(dds))) # check res_treatment and dds genes are in same order
  # BiocParallel::register(BiocParallel::MulticoreParam(workers = 32))
  # 
  # res_treatment_shrink = lfcShrink(dds, 
  #                                  coef = "treatment_IFN_vs_untreated", 
  #                                  res=res_treatment,
  #                                  lfcThreshold = 1,
  #                                  parallel = TRUE)
  # takes long time and memory for so many donors
  
  message("Of those with abs(log2FC)>1, how many genes are significant at 5% FDR level?")
  message(sum(res_treatment$padj<0.05 & abs(res_treatment$log2FoldChange) > 1, na.rm = TRUE))
  # 692
  
  res_treatment = res_treatment[order(abs(res_treatment$log2FoldChange),decreasing = TRUE),]
  res_treatment
  
  # add ensembl ids
  res_treatment$ensembl = mapIds(org.Hs.eg.db,
                                 keys=row.names(res_treatment),
                                 column="ENSEMBL",
                                 keytype= "SYMBOL",
                                 multiVals="first")
  res_treatment$symbol = row.names(res_treatment)
  
  DESeq2::plotMA(res_treatment, ylim=c(-6,6))
  
  write.table(res_treatment, file = paste0(outdir,"/DiffExpr_",
                                           treatment,
                                           "_all_genes_negative_lower_in_",treatment,".txt"),
              quote = F, sep = "\t",row.names = F, col.names = T)
  res_sign_treatment = res_treatment %>%  
    as.data.frame() %>%
    filter(padj<0.05 & abs(log2FoldChange) > 1)
  
  write.table(res_sign_treatment,
              file = paste0(outdir,"/DiffExpr_",treatment,
                            "_padj05_log2FC1_negative_lower_in_",
                            treatment,".txt"), quote = F, sep = "\t",row.names = F, col.names = T)
  
  
  ### visualization
  genes = rownames(vst_mat)
  names = rownames(colData(vsttr))
  meta = as_tibble(colData(vsttr)) %>%
    dplyr::mutate(name=names)
  
  pdf(paste0(outdir,"/DiffExpr_",treatment,"vsUntr_top30_abs_log2FC.pdf"),width = 7,height = 7)
  for(i in 1:30){
    message("Plotting top 30 genes sorted by abs(log2FC)")
    
    gene_to_plot = vst_mat %>%
      as_tibble() %>%
      dplyr::mutate(gene = genes) %>%
      dplyr::filter(gene==res_sign_treatment$symbol[i]) %>%
      tidyr::pivot_longer(cols = !gene, values_to = "vst_count") %>%
      dplyr::left_join(.,meta)
    
    
    p = ggplot(gene_to_plot, aes(x=treatment, y=vst_count,col=pool)) +
      scale_y_log10() +
      geom_point(position=position_jitter(width=.1,height=0), size=3,alpha=0.6) + 
      theme_bw() + 
      ggtitle( res_sign_treatment$symbol[i]) + 
      ylab("Counts after VST transformation")
    plot(p)
  }
  
  dev.off()
  
}

# are the counts correlated with the number of cells before correction for donor effects

genes = rownames(vst_mat)
names = rownames(colData(vsttr))
meta = as_tibble(colData(vsttr)) %>%
  dplyr::mutate(name=names)
# check top 10 genes
gene_to_plot = vst_mat %>%
  as_tibble() %>%
  dplyr::mutate(gene = genes) %>%   
  dplyr::filter(gene %in% res_sign_treatment$symbol[1:10]) %>%
  tidyr::pivot_longer(cols = !gene, values_to = "vst_count") %>%
  dplyr::left_join(.,meta)

png(paste0(outdir,"/cor_expr_log10ncells_vst_blindTRUE.png"),
    width = 16, height = 14, res = 400,units = "in", type = "cairo")
ggplot(gene_to_plot,aes(x=log10_ncells,y=vst_count, col=treatment)) +
  geom_point(alpha=0.6) + 
  theme_minimal() + 
  facet_wrap(vars(gene))
dev.off()

# are within-treatment differences between donors driven by ncells?
# account for pool with nested design?

# SAM's gene of interest

png(paste0(outdir,"/SAMHD1_normexpr_results.png"),
    width = 5, height = 4, res = 400,units = "in", type = "cairo")

gene_to_plot = vst_mat %>%
  as_tibble() %>%
  dplyr::mutate(gene = genes) %>%
  dplyr::filter(gene=="SAMHD1") %>%
  tidyr::pivot_longer(cols = !gene, values_to = "vst_count") %>%
  dplyr::left_join(.,meta)


p = ggplot(gene_to_plot, aes(x=treatment, y=vst_count,col=pool)) +
  scale_y_log10() +
  geom_point(position=position_jitter(width=.1,height=0), size=3,alpha=0.6) + 
  theme_bw() + 
  ggtitle( "SAMHD1") + 
  ylab("Counts after VST transformation")
plot(p)

dev.off()

# IFN vs LPS. ############

message("Getting contrasts for IFN vs LPS")

res_treatment = results(dds, 
                        name = paste0("treatment_IFN_vs_LPS"),
                        lfcThreshold = 1)
message("Number of padj<0.01 and abs(log2FC>1)")
message(table(res_treatment$padj<0.01)) # 676
identical(rownames(res_treatment),rownames(assay(dds))) # check res_treatment and dds genes are in same order
# BiocParallel::register(BiocParallel::MulticoreParam(workers = 32))
# 
# res_treatment_shrink = lfcShrink(dds, 
#                                  coef = "treatment_IFN_vs_untreated", 
#                                  res=res_treatment,
#                                  lfcThreshold = 1,
#                                  parallel = TRUE)
# takes long time and memory for so many donors

message("Of those with abs(log2FC)>1, how many genes are significant at 5% FDR level?")
message(sum(res_treatment$padj<0.05 & abs(res_treatment$log2FoldChange) > 1, na.rm = TRUE))
# 692

res_treatment = res_treatment[order(abs(res_treatment$log2FoldChange),decreasing = TRUE),]
res_treatment

# add ensembl ids
res_treatment$ensembl = mapIds(org.Hs.eg.db,
                               keys=row.names(res_treatment),
                               column="ENSEMBL",
                               keytype= "SYMBOL",
                               multiVals="first")
res_treatment$symbol = row.names(res_treatment)

DESeq2::plotMA(res_treatment, ylim=c(-6,6))

write.table(res_treatment, file = paste0(outdir,"/DiffExpr_",
                                         treatment,
                                         "_all_genes_negative_lower_in_",treatment,".txt"),
            quote = F, sep = "\t",row.names = F, col.names = T)
res_sign_treatment = res_treatment %>%  
  as.data.frame() %>%
  filter(padj<0.05 & abs(log2FoldChange) > 1)

write.table(res_sign_treatment,
            file = paste0(outdir,"/DiffExpr_",treatment,
                          "_padj05_log2FC1_negative_lower_in_",
                          treatment,".txt"), quote = F, sep = "\t",row.names = F, col.names = T)


### visualization
genes = rownames(vst_mat)
names = rownames(colData(vsttr))
meta = as_tibble(colData(vsttr)) %>%
  dplyr::mutate(name=names)

pdf(paste0(outdir,"/DiffExpr_",treatment,"vsUntr_top30_abs_log2FC.pdf"),width = 7,height = 7)
for(i in 1:30){
  message("Plotting top 30 genes sorted by abs(log2FC)")
  
  gene_to_plot = vst_mat %>%
    as_tibble() %>%
    dplyr::mutate(gene = genes) %>%
    dplyr::filter(gene==res_sign_treatment$symbol[1]) %>%
    tidyr::pivot_longer(cols = !gene, values_to = "vst_count") %>%
    dplyr::left_join(.,meta)
  
  
  p = ggplot(gene_to_plot, aes(x=treatment, y=vst_count,col=pool)) +
    scale_y_log10() +
    geom_point(position=position_jitter(width=.1,height=0), size=3,alpha=0.6) + 
    theme_bw() + 
    ggtitle( res_sign_treatment$symbol[1]) + 
    ylab("Counts after VST transformation")
  plot(p)
}

dev.off()


######## interactions ########
############ treatment effect accounting for donors and pools ########
### pre-correcting for pool and then fitting donor-treatment interactions

dds$donor_id = as.factor(dds$donor_id)
dds$donor_id = droplevels(dds$donor_id)

# adjusted = sva::ComBat_seq(counts(dds), batch=dds$pool, group=NULL)
adjusted = sva::ComBat_seq(counts(dds), batch=dds$pool, 
                           covar_mod=cbind(dds$ncells,dds$treatment,dds$donor_id) # signals to keep as much as possible for DESeq2 step
                           )

mat = model.matrix(~ (1/ncells) + donor_id*treatment,colData(dds))
mat
is.fullrank(mat)
colnames(mat)
mat = mat[,colSums(mat)>0]
colnames(mat)
is.fullrank(mat)

# drop columns we're not interested in
mat = mat[,c(1,190:556)]
is.fullrank(mat)

dds = DESeqDataSetFromMatrix(countData = adjusted,
                             colData =  colData(dds),
                             design = ~ treatment) # dummy here, fit matrix later

dds = estimateSizeFactors(dds)

BiocParallel::register(BiocParallel::MulticoreParam(workers = 32))
dds = DESeq(dds,parallel=TRUE,
            full=mat,
            betaPrior = FALSE) # for Wald tests

write_rds(dds,paste0(outdir,"/dds_interactions_pool_corrected.rds"))

png(paste0(outdir,"/disp_estimates_interactions_pool_corrected.png"),
    width = 9, height = 7, res = 400,units = "in", type = "cairo")
plotDispEsts(dds) # weird dispersion, should decrease with more counts and follow the mean
# https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html
# a bit better after removing pools 9 and 10, but still a straight line
# also mostly straight if not excluding those pools
dev.off()

resultsNames(dds)

for(treatment in c("IFN","LPS")){
  message("Getting contrasts for ", treatment," treatment")
  
  res_treatment = results(dds, 
                          name = paste0("treatment_",treatment,"_vs_untreated"),
                          lfcThreshold = 1)
  message("Number of padj<0.01 and abs(log2FC>1)")
  message(table(res_treatment$padj<0.01)) # 676
  identical(rownames(res_treatment),rownames(assay(dds))) # check res_treatment and dds genes are in same order
  # BiocParallel::register(BiocParallel::MulticoreParam(workers = 32))
  # 
  # res_treatment_shrink = lfcShrink(dds, 
  #                                  coef = "treatment_IFN_vs_untreated", 
  #                                  res=res_treatment,
  #                                  lfcThreshold = 1,
  #                                  parallel = TRUE)
  # takes long time and memory for so many donors
  
  message("Of those with abs(log2FC)>1, how many genes are significant at 5% FDR level?")
  message(sum(res_treatment$padj<0.05 & abs(res_treatment$log2FoldChange) > 1, na.rm = TRUE))
  # 692
  
  res_treatment = res_treatment[order(abs(res_treatment$log2FoldChange),decreasing = TRUE),]
  res_treatment
  
  # add ensembl ids
  res_treatment$ensembl = mapIds(org.Hs.eg.db,
                                 keys=row.names(res_treatment),
                                 column="ENSEMBL",
                                 keytype= "SYMBOL",
                                 multiVals="first")
  res_treatment$symbol = row.names(res_treatment)
  
  DESeq2::plotMA(res_treatment, ylim=c(-6,6))
  
  write.table(res_treatment, file = paste0(outdir,"/DiffExpr_",
                                           treatment,
                                           "_all_genes_negative_lower_in_",treatment,".txt"),
              quote = F, sep = "\t",row.names = F, col.names = T)
  res_sign_treatment = res_treatment %>%  
    as.data.frame() %>%
    filter(padj<0.05 & abs(log2FoldChange) > 1)
  
  write.table(res_sign_treatment,
              file = paste0(outdir,"/DiffExpr_",treatment,
                            "_padj05_log2FC1_negative_lower_in_",
                            treatment,".txt"), quote = F, sep = "\t",row.names = F, col.names = T)
  
  
  ### visualization
  genes = rownames(vst_mat)
  names = rownames(colData(vsttr))
  meta = as_tibble(colData(vsttr)) %>%
    dplyr::mutate(name=names)
  
  pdf(paste0(outdir,"/DiffExpr_",treatment,"vsUntr_top30_abs_log2FC.pdf"),width = 7,height = 7)
  for(i in 1:30){
    message("Plotting top 30 genes sorted by abs(log2FC)")
    
    gene_to_plot = vst_mat %>%
      as_tibble() %>%
      dplyr::mutate(gene = genes) %>%
      dplyr::filter(gene==res_sign_treatment$symbol[i]) %>%
      tidyr::pivot_longer(cols = !gene, values_to = "vst_count") %>%
      dplyr::left_join(.,meta)
    
    
    p = ggplot(gene_to_plot, aes(x=treatment, y=vst_count,col=pool)) +
      scale_y_log10() +
      geom_point(position=position_jitter(width=.1,height=0), size=3,alpha=0.6) + 
      theme_bw() + 
      ggtitle( res_sign_treatment$symbol[i]) + 
      ylab("Counts after VST transformation")
    plot(p)
  }
  
  dev.off()
  
}

