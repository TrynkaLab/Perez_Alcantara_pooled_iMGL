# Differential expression with limma+voom

library(limma)
library(edgeR)
library(tidyverse) # For ggplot2 and easy manipulation of data
library(patchwork) # To combine plots
set.seed(123)
source("./helpers.R")
outdir="../../data/results/5.1.Diff_expr_limma"
dir.create(outdir, recursive = T)


#### pool as random effect
### first the effects of treatment without interaction
######### reading in again to avoid confusion
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
metadata_rows = rownames(metadata)

# subset to non-proliferating only
metadata = metadata %>%
  dplyr::filter(proliferation_status == "Not_proliferating")

pseudobulk = pseudobulk %>%
  dplyr::select(match(metadata$pseudobulk_colnames,colnames(.)))
# same order
identical(colnames(pseudobulk),metadata$pseudobulk_colnames)


discarded = metadata$ncells < 100
table(discarded)

pseudobulk = pseudobulk[,!discarded]
metadata = metadata[!discarded,]
length(unique(metadata$donor)) # 194 donors, 194/250 = 78%
# remove samples from pool 9
# 
# discarded = metadata$pool %in% c("pool9")
# pseudobulk = pseudobulk[,!discarded]
# metadata = metadata[!discarded,]
# length(unique(metadata$donor)) # 172 donors, 172/250 = 69%

### creating DGElist object

### remove pool batch effects

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
metadata$donor = factor(metadata$donor)
metadata$donor = relevel(metadata$donor,ref = "hegp_3")
mat = model.matrix(~ 0 + treatment +  donor + log10(ncells),metadata)
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
v = voom(dge, design = mat, plot=TRUE, block = metadata$pool, correlation = cor$consensus)
cor = duplicateCorrelation(v, design = mat, block =  metadata$pool) # extract duplicates again
cor$consensus


# heatmap with voom-transformed counts
# pseudobulk_corrected = limma::removeBatchEffect(pseudobulk, metadata$pool)
# need to correct per treatment
# vst_cor = cor(pseudobulk_corrected)
# # Plot heatmap
# 
# pools = c(rainbow(16))
# names(pools) = paste0("pool",2:17)
# ann_color = list("treatment"=c("untreated"="black","LPS"="yellow","IFN"="red"),
#                  "proliferation_status"=c("Not_proliferating"="lightblue","Proliferating"="darkblue"),
#                  "pool"=pools)
# 
# png(paste0(outdir,"/cor_heatmap_voom_pool_corrected_allpools.png"),
#     width = 16, height = 14, res = 400,units = "in", type = "cairo")
# vst_cor %>%
#   magrittr::set_rownames(metadata$donor) %>%
#   pheatmap::pheatmap(., annotation = metadata[, c("ncells","proliferation_status","treatment",
#                                         "pool"), drop=F],
#            show_colnames = FALSE,
#            show_rownames = FALSE,
#            annotation_colors = ann_color)
# 
# dev.off()

# fit 
fit = lmFit(v, design=mat,block = metadata$pool, correlation  = cor$consensus)

cont = c("treatmentIFN-treatmentuntreated",
         "treatmentLPS-treatmentuntreated",
         "treatmentIFN-treatmentLPS")

contrast.matrix = makeContrasts(contrasts = cont,
                                levels=mat)
colnames(contrast.matrix) = c("IFNvsUntreated","LPSvsUntreated","IFNvsLPS")

fit2 = contrasts.fit(fit, contrast.matrix)
fit2 = eBayes(fit2)

summary_additive = summary(decideTests(fit2,lfc = 1, 
                                       adjust.method = "fdr",p.value=0.01))

res_additive = list()
for(cont in colnames(contrast.matrix)){
  res_additive[[cont]] = topTable(fit2, coef=cont,
                                  number = length(fit2$coefficients)) %>%
    dplyr::arrange(desc(abs(logFC))) %>%
    tibble::rownames_to_column(var = "symbol")
  
  write.table(res_additive[[cont]],paste0(outdir,
                                          "/DiffExpr_",cont,"all_genes_negative_lower_in_latter.txt"),
              quote = F, sep = "\t",row.names = F, col.names = T)
}
saveRDS(res_additive,paste0(outdir,"/additive_res.rds"))

for(cont in colnames(contrast.matrix)){
  g = res_additive[[cont]] %>%
    dplyr::slice_head(n=30)
  pdf(paste0(outdir,"/DiffExpr_",cont,"top30_abs_log2FC.pdf"),width = 7,height = 7)
  for(i in 1:30){
    message("Plotting top 30 genes sorted by abs(log2FC)")
    expr = as.data.frame(v$E)
    expr$gene = rownames(expr)
    
    gene_to_plot = expr %>%
      dplyr::filter(gene==g$symbol[i]) %>%
      tidyr::pivot_longer(cols = !gene, 
                          values_to = "voom_count",
                          names_to = "pseudobulk_colnames") %>%
      dplyr::left_join(.,metadata)
    
    
    p = ggplot(gene_to_plot, aes(x=treatment, y=voom_count,col=pool)) +
      geom_point(position=position_jitter(width=.1,height=0), size=3,alpha=0.6) + 
      theme_bw() + 
      ggtitle( g$symbol[i]) + 
      ylab("Normalised log2(counts+0.5)")
    plot(p)
  }
  
  dev.off()
}

# 
# # DESeq2 vs limma
# deseq_ifn=read.table("../../data/results/5.1.Diff_expr_DESeq2/fit_pool_all_pools/DiffExpr_IFN_all_genes_negative_lower_in_IFN.txt",
#                      header = TRUE) %>%
#   dplyr::rename(logFC=log2FoldChange) %>%
#   dplyr::select(symbol,logFC) %>%
#   dplyr::mutate(method="DESeq2")
# 
# 
# p1 = res_additive$IFNvsUntreated %>%
#   dplyr::select(symbol,logFC) %>%
#   dplyr::mutate(method="limma") %>%
#   dplyr::bind_rows(deseq_ifn) %>%
#   tidyr::pivot_wider(id_cols = symbol,values_from = logFC,names_from = method) %>%
#   ggplot(.,aes(x=DESeq2,y=limma)) + 
#   geom_point(alpha=0.6) + 
#   theme_minimal() + 
#   ggtitle("Correlation DESeq2-limma log2FC in IFN vs untreated")
# 
# pdf(paste0(outdir,"/log2FC_IFNvsUntreated_limma_DESeq2.pdf"),width = 7,height = 7)
# 
# p1
# 
# dev.off()
# 
# deseq_ifn=read.table("../../data/results/5.1.Diff_expr_DESeq2/fit_pool_all_pools/DiffExpr_IFN_all_genes_negative_lower_in_IFN.txt",
#                      header = TRUE) %>%
#   dplyr::select(symbol,padj) %>%
#   dplyr::mutate(method="DESeq2")
# 
# p2 = res_additive$IFNvsUntreated %>%
#   dplyr::select(symbol,adj.P.Val) %>%
#   dplyr::rename(padj = adj.P.Val) %>%
#   dplyr::mutate(method="limma") %>%
#   dplyr::bind_rows(deseq_ifn) %>%
#   dplyr::mutate(padj=-log10(padj)) %>%
#   dplyr::mutate(padj=case_when(padj=="Inf"~350,
#                                .default = padj)) %>% # substitute for very extreme pvalue, at the limit of the precision
#   
#   tidyr::pivot_wider(id_cols = symbol,values_from = padj,names_from = method) %>%
#   ggplot(.,aes(x=DESeq2,y=limma)) + 
#   scale_x_continuous(trans=scales::pseudo_log_trans(base = 10)) +
#   scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
#   geom_point(alpha=0.6) + 
#   theme_minimal() + 
#   ggtitle("Correlation DESeq2-limma p-adj in IFN vs untreated")
# 
# pdf(paste0(outdir,"/pvals_IFNvsUntreated_limma_DESeq2.pdf"),width = 7,height = 7)
# 
# p2
# 
# dev.off()
# 
# 
# ##### pool as random effects with interactions ###############
# #############
# pseudobulk =readr::read_tsv("../../data/results/5.prepare_differential_expression_pseudobulks/combined_pseudobulk_raw_counts_treatmentxprolifxidxpool.txt") %>%
#   tibble::column_to_rownames(var = "gene")
# metadata = readr::read_tsv("../../data/results/5.prepare_differential_expression_pseudobulks/metadata_treatmentxprolifxidxpool.txt") %>%
#   tibble::column_to_rownames(var = "cols_names") %>%
#   dplyr::rename(ncells = count,donor= donor_id) %>%
#   dplyr::mutate(pseudobulk_colnames=paste(treatment,proliferation_status,donor,pool,sep = "_"),
#                 donor = case_when(donor=="Arlene-003" ~ "Arlene_3",
#                                   donor=="Cindy-005" ~ "Cindy_5",
#                                   donor=="Dexter-006" ~ "Dexter_6",
#                                   donor=="Fiona-010" ~ "Fiona_10",
#                                   donor=="Gilma-009" ~ "Gilma_9",
#                                   donor=="Hector-011" ~ "Hector_11",
#                                   donor=="Imani-012" ~ "Imani_12",
#                                   donor=="Javier-013" ~ "Javier_13",
#                                   donor=="Keoni-014" ~ "Keoni_14",
#                                   donor=="Olaf-018" ~ "Olaf_18",
#                                   donor=="Bertha-004" ~ "Bertha_4",
#                                   donor=="Mindy-016" ~ "Mindy_16",
#                                   .default = donor))
# metadata_rows = rownames(metadata)
# 
# # subset to non-proliferating only
# metadata = metadata %>%
#   dplyr::filter(proliferation_status == "Not_proliferating")
# 
# pseudobulk = pseudobulk %>%
#   dplyr::select(match(metadata$pseudobulk_colnames,colnames(.)))
# # same order
# identical(colnames(pseudobulk),metadata$pseudobulk_colnames)
# 
# 
# discarded = metadata$ncells < 100
# table(discarded)
# 
# pseudobulk = pseudobulk[,!discarded]
# metadata = metadata[!discarded,]
# length(unique(metadata$donor)) # 189 donors, 189/248 = 76%
# # remove samples from pool 9
# 
# discarded = metadata$pool %in% c("pool9")
# pseudobulk = pseudobulk[,!discarded]
# metadata = metadata[!discarded,]
# length(unique(metadata$donor)) # 167 donors, 189/248 = 67%
# 
# ### creating DGElist object
# 
# ### remove pool batch effects
# #pseudobulk_corrected = limma::removeBatchEffect(pseudobulk, metadata$pool)
# 
# dge = edgeR::DGEList(counts=pseudobulk)
# 
# # ~    (1/ncells) + pool + treatment)
# 
# 
# # remove genes not expressed in 30% samples (min 1 CPM)
# isexpr = rowSums(edgeR::cpm(dge) > 1) >= floor(0.3*ncol(dge))
# message(paste0("There are ", sum(isexpr), " genes with over 1 count per million in 30% of samples"))
# dge = dge[isexpr,,keep.lib.sizes=FALSE]  
# dim(dge)
# 
# 
# ## how many have 0 counts
# table(rowSums(dge$counts)==0)
# # none
# 
# # TMM normalization:
# dge = edgeR::calcNormFactors(dge)
# 
# 
# ### create design matrix
# metadata$treatment = factor(metadata$treatment)
# metadata$treatment = relevel(metadata$treatment,ref = "LPS")
# metadata$donor = factor(metadata$donor)
# metadata$donor = relevel(metadata$donor,ref = "hegp_3")
# mat = model.matrix(~ 0 + donor*treatment + log10(ncells),metadata)
# mat = as.data.frame(mat)
# qr(mat)$rank
# ncol(mat)
# is.fullrank(mat)
# colnames(mat)
# mat = mat[,colSums(mat)>0]
# colnames(mat)
# is.fullrank(mat)
# 
# # drop columns we're not interested in
# mat = mat[,c(1,7,168:493)] # leaving shared donor hegp_3 and aowh_2 in to 
# # estimate mean donor effects (limited to the two donors shared across all pools)
# is.fullrank(mat)
# colnames(mat) = make.names(colnames(mat))
# 
# 
# v = voom(dge, design = mat, plot=TRUE)
# 
# cor = duplicateCorrelation(v, design=mat, block = metadata$pool)
# 
# cor$consensus # within-pool correlation
# # voom weights may have changed:
# v = voom(dge, design = mat, plot=TRUE, block = metadata$pool, correlation = cor$consensus)
# cor = duplicateCorrelation(v, design = mat, block =  metadata$pool) # extract duplicates again
# cor$consensus
# 
# # fit 
# fit = lmFit(v, design=mat,block = metadata$pool, correlation  = cor$consensus)
# 
# donor_contrasts = stringr::str_split(colnames(mat),pattern = "\\.",n=2,simplify = TRUE)[,1]
# donor_contrasts = stringr::str_replace(donor_contrasts,
#                                        pattern = "donor",replacement = "")
# donor_contrasts = unique(donor_contrasts[6:328])
# # exclude some donors that don't have levels in certain treatments
# donor_contrasts = donor_contrasts[!donor_contrasts %in% c("kuco_1","qayj_4","seru_1","vils_1","wuye_3", # not in untreated
#                                                           "nurh","peoj_1","puhk_2","tolg_4")] # not in IFN
# 
# length(donor_contrasts) # 
# cont = c(paste0("(donor",donor_contrasts,".treatmentIFN-donor",
#                 donor_contrasts,".treatmentuntreated)-(treatmentIFN-treatmentuntreated)-((donorhegp_3 + donoraowh_2)/2) "))
# ## fitting interaction contrasts to test the null hypothesis that 
# # donor specific effect of treatment - AVERAGE across all donors effect of treatment =0
# # I've done it this way because there is no clear baseline for the donor effects
# 
# contrast.matrix = makeContrasts(contrasts = cont,
#                                 levels=mat)
# colnames(contrast.matrix) = c(paste0(donor_contrasts,"_IFN_interaction"))
# 
# fit2 = contrasts.fit(fit, contrast.matrix)
# fit2 = eBayes(fit2)
# 
# summary(decideTests(fit2,lfc=1))[,1:3]
# summary_interactions5 = summary(decideTests(fit2,lfc = 1, 
#                                             adjust.method = "fdr",p.value=0.01))
# 
# # all results for interactions and IFN treatment effect
# ifn_interactions = list()
# # hegp_3 and some others to check
# for(cont in  colnames(contrast.matrix)){
#   ifn_interactions[[cont]] = topTable(fit2, coef=cont,
#                                       number = length(fit2$coefficients)) %>%
#     dplyr::arrange(desc(abs(logFC)))
# }
# 
# saveRDS(ifn_interactions,paste0(outdir,"/IFN_interactions.rds"))
# ## plot interactions
# 
# pdf(paste0(outdir,"/DiffExpr_IFN_vs_untreated_interactions_top30_abs_log2FC.pdf"),width = 15,height = 18)
# 
# for(cont in colnames(contrast.matrix)){
#   message("Plotting top 30 genes sorted by abs(log2FC)")
#   
#   g = ifn_interactions[[cont]] %>%
#     dplyr::slice_head(n=30) %>%
#     tibble::rownames_to_column(var = "symbol")
#   plist = list()
#   for(i in 1:30){
#     gen=g$symbol[i]
#     expr = as.data.frame(v$E)
#     expr$gene = rownames(expr)
#     d = stringr::str_replace(cont,pattern = "_IFN_interaction",replacement = "")
#     gene_to_plot = expr %>%
#       dplyr::filter(gene==gen) %>%
#       tidyr::pivot_longer(cols = !gene, 
#                           values_to = "voom_count",
#                           names_to = "pseudobulk_colnames") %>%
#       dplyr::left_join(.,metadata) %>%
#       dplyr::filter(donor == d)
#     average_baseline_to_plot =    expr %>%
#       dplyr::filter(gene==gen) %>%
#       tidyr::pivot_longer(cols = !gene, 
#                           values_to = "voom_count",
#                           names_to = "pseudobulk_colnames") %>%
#       dplyr::left_join(.,metadata) %>%
#       dplyr::filter(donor %in% c("aowh_2","hegp_3")) %>%
#       dplyr::group_by(treatment,pool,gene) %>%
#       dplyr::summarise(voom_count = mean(voom_count)) %>%
#       dplyr::mutate(donor = "average") %>%
#       dplyr::ungroup()
#     
#     average_gene_to_plot = gene_to_plot %>%
#       dplyr::group_by(treatment,gene,donor) %>%
#       dplyr::summarise(voom_count = mean(voom_count)) %>%
#       dplyr::ungroup()
#     averages = average_baseline_to_plot %>%
#       dplyr::group_by(treatment,gene,donor) %>%
#       dplyr::summarise(voom_count = mean(voom_count)) %>%
#       dplyr::bind_rows(.,average_gene_to_plot) %>%
#       dplyr::filter(treatment !="LPS") 
#     
#     colors = c("#17BEBB","#EB5160")
#     names(colors) = c("average",d)
#     plist[[i]] = gene_to_plot %>%
#       dplyr::select(gene,voom_count,donor,treatment,pool) %>%
#       dplyr::bind_rows(average_baseline_to_plot) %>%
#       dplyr::filter(treatment !="LPS") %>%
#       
#       ggplot(., aes(x=treatment, y=voom_count,col=donor)) +
#       geom_point(position=position_jitter(width=.1,height=0), size=3,alpha=0.6) + 
#       geom_line(data = averages, aes(x = treatment, y = voom_count, color = donor,group = donor), linewidth = 1) +
#       scale_color_manual(values = colors)+
#       theme_bw() + 
#       ggtitle( gen) + 
#       ylab("Normalised log2(counts+0.5)")
#   }
#   
#   plot(patchwork::wrap_plots(plist,ncol = 5,nrow = 6) + patchwork::plot_annotation(title = cont))
# }
# dev.off()
# 
# ### for LPS vs untreated
# 
# pseudobulk =readr::read_tsv("../../data/results/5.prepare_differential_expression_pseudobulks/combined_pseudobulk_raw_counts_treatmentxprolifxidxpool.txt") %>%
#   tibble::column_to_rownames(var = "gene")
# metadata = readr::read_tsv("../../data/results/5.prepare_differential_expression_pseudobulks/metadata_treatmentxprolifxidxpool.txt") %>%
#   tibble::column_to_rownames(var = "cols_names") %>%
#   dplyr::rename(ncells = count,donor= donor_id) %>%
#   dplyr::mutate(pseudobulk_colnames=paste(treatment,proliferation_status,donor,pool,sep = "_"),
#                 donor = case_when(donor=="Arlene-003" ~ "Arlene_3",
#                                   donor=="Cindy-005" ~ "Cindy_5",
#                                   donor=="Dexter-006" ~ "Dexter_6",
#                                   donor=="Fiona-010" ~ "Fiona_10",
#                                   donor=="Gilma-009" ~ "Gilma_9",
#                                   donor=="Hector-011" ~ "Hector_11",
#                                   donor=="Imani-012" ~ "Imani_12",
#                                   donor=="Javier-013" ~ "Javier_13",
#                                   donor=="Keoni-014" ~ "Keoni_14",
#                                   donor=="Olaf-018" ~ "Olaf_18",
#                                   donor=="Bertha-004" ~ "Bertha_4",
#                                   donor=="Mindy-016" ~ "Mindy_16",
#                                   .default = donor))
# metadata_rows = rownames(metadata)
# 
# # subset to non-proliferating only
# metadata = metadata %>%
#   dplyr::filter(proliferation_status == "Not_proliferating")
# 
# pseudobulk = pseudobulk %>%
#   dplyr::select(match(metadata$pseudobulk_colnames,colnames(.)))
# # same order
# identical(colnames(pseudobulk),metadata$pseudobulk_colnames)
# 
# 
# discarded = metadata$ncells < 100
# table(discarded)
# 
# pseudobulk = pseudobulk[,!discarded]
# metadata = metadata[!discarded,]
# length(unique(metadata$donor)) # 189 donors, 189/248 = 76%
# # remove samples from pool 9
# 
# discarded = metadata$pool %in% c("pool9")
# pseudobulk = pseudobulk[,!discarded]
# metadata = metadata[!discarded,]
# length(unique(metadata$donor)) # 167 donors, 189/248 = 67%
# 
# ### creating DGElist object
# 
# ### remove pool batch effects
# #pseudobulk_corrected = limma::removeBatchEffect(pseudobulk, metadata$pool)
# 
# dge = edgeR::DGEList(counts=pseudobulk)
# 
# # ~    (1/ncells) + pool + treatment)
# 
# 
# # remove genes not expressed in 30% samples (min 1 CPM)
# isexpr = rowSums(edgeR::cpm(dge) > 1) >= floor(0.3*ncol(dge))
# message(paste0("There are ", sum(isexpr), " genes with over 1 count per million in 30% of samples"))
# dge = dge[isexpr,,keep.lib.sizes=FALSE]  
# dim(dge)
# 
# 
# ## how many have 0 counts
# table(rowSums(dge$counts)==0)
# # none
# 
# # TMM normalization:
# dge = edgeR::calcNormFactors(dge)
# 
# 
# ### create design matrix
# metadata$treatment = factor(metadata$treatment)
# metadata$treatment = relevel(metadata$treatment,ref = "IFN")
# metadata$donor = factor(metadata$donor)
# metadata$donor = relevel(metadata$donor,ref = "hegp_3")
# mat = model.matrix(~ 0 + donor*treatment + log10(ncells),metadata)
# mat = as.data.frame(mat)
# qr(mat)$rank
# ncol(mat)
# is.fullrank(mat)
# colnames(mat)
# mat = mat[,colSums(mat)>0]
# colnames(mat)
# is.fullrank(mat)
# 
# # drop columns we're not interested in
# mat = mat[,c(1,7,168:487)] # leaving shared donor hegp_3 and aowh_2 in to 
# # estimate mean donor effects (limited to the two donors shared across all pools)
# is.fullrank(mat)
# colnames(mat) = make.names(colnames(mat))
# 
# 
# v = voom(dge, design = mat, plot=TRUE)
# 
# cor = duplicateCorrelation(v, design=mat, block = metadata$pool)
# 
# cor$consensus # within-pool correlation
# # voom weights may have changed:
# v = voom(dge, design = mat, plot=TRUE, block = metadata$pool, correlation = cor$consensus)
# cor = duplicateCorrelation(v, design = mat, block =  metadata$pool) # extract duplicates again
# cor$consensus
# 
# # fit 
# fit = lmFit(v, design=mat,block = metadata$pool, correlation  = cor$consensus)
# 
# donor_contrasts = stringr::str_split(colnames(mat),pattern = "\\.",n=2,simplify = TRUE)[,1]
# donor_contrasts = stringr::str_replace(donor_contrasts,
#                                        pattern = "donor",replacement = "")
# donor_contrasts = unique(donor_contrasts[6:322])
# # exclude some donors that don't have levels in certain treatments
# donor_contrasts = donor_contrasts[!donor_contrasts %in% c("kuco_1","qayj_4","seru_1","vils_1","wuye_3", # not in untreated
#                                                           "coxy_33","cuhk_2","letw_1","nurh","peoj_1","qolg_1","tolg_4","veku_2")] # not in LPS
# 
# 
# length(donor_contrasts) # 
# cont = c(paste0("(donor",donor_contrasts,".treatmentLPS-donor",
#                 donor_contrasts,".treatmentuntreated)-(treatmentLPS-treatmentuntreated)-((donorhegp_3 + donoraowh_2)/2) "))
# ## fitting interaction contrasts to test the null hypothesis that 
# # donor specific effect of treatment - AVERAGE across all donors effect of treatment =0
# # I've done it this way because there is no clear baseline for the donor effects
# 
# contrast.matrix = makeContrasts(contrasts = cont,
#                                 levels=mat)
# colnames(contrast.matrix) = c(paste0(donor_contrasts,"_LPS_interaction"))
# 
# fit2 = contrasts.fit(fit, contrast.matrix)
# fit2 = eBayes(fit2)
# 
# summary(decideTests(fit2,lfc=1))[,1:3]
# summary_interactions6 = summary(decideTests(fit2,lfc = 1, 
#                                             adjust.method = "fdr",p.value=0.01))
# 
# # all results for interactions and IFN treatment effect
# lps_interactions = list()
# # hegp_3 and some others to check
# for(cont in  colnames(contrast.matrix)){
#   lps_interactions[[cont]] = topTable(fit2, coef=cont,
#                                       number = length(fit2$coefficients)) %>%
#     dplyr::arrange(desc(abs(logFC)))
# }
# 
# saveRDS(lps_interactions,paste0(outdir,"/LPS_interactions.rds"))
# ## plot interactions
# 
# pdf(paste0(outdir,"/DiffExpr_LPS_vs_untreated_interactions_top30_abs_log2FC.pdf"),width = 15,height = 18)
# 
# for(cont in colnames(contrast.matrix)){
#   message("Plotting top 30 genes sorted by abs(log2FC)")
#   
#   g = lps_interactions[[cont]] %>%
#     dplyr::slice_head(n=30) %>%
#     tibble::rownames_to_column(var = "symbol")
#   plist = list()
#   for(i in 1:30){
#     gen=g$symbol[i]
#     expr = as.data.frame(v$E)
#     expr$gene = rownames(expr)
#     d = stringr::str_replace(cont,pattern = "_LPS_interaction",replacement = "")
#     gene_to_plot = expr %>%
#       dplyr::filter(gene==gen) %>%
#       tidyr::pivot_longer(cols = !gene, 
#                           values_to = "voom_count",
#                           names_to = "pseudobulk_colnames") %>%
#       dplyr::left_join(.,metadata) %>%
#       dplyr::filter(donor == d)
#     average_baseline_to_plot =    expr %>%
#       dplyr::filter(gene==gen) %>%
#       tidyr::pivot_longer(cols = !gene, 
#                           values_to = "voom_count",
#                           names_to = "pseudobulk_colnames") %>%
#       dplyr::left_join(.,metadata) %>%
#       dplyr::filter(donor %in% c("aowh_2","hegp_3")) %>%
#       dplyr::group_by(treatment,pool,gene) %>%
#       dplyr::summarise(voom_count = mean(voom_count)) %>%
#       dplyr::mutate(donor = "average") %>%
#       dplyr::ungroup()
#     
#     average_gene_to_plot = gene_to_plot %>%
#       dplyr::group_by(treatment,gene,donor) %>%
#       dplyr::summarise(voom_count = mean(voom_count)) %>%
#       dplyr::ungroup()
#     averages = average_baseline_to_plot %>%
#       dplyr::group_by(treatment,gene,donor) %>%
#       dplyr::summarise(voom_count = mean(voom_count)) %>%
#       dplyr::bind_rows(.,average_gene_to_plot) %>%
#       dplyr::filter(treatment !="IFN") 
#     colors = c("#17BEBB","#EB5160")
#     names(colors) = c("average",d)
#     plist[[i]] = gene_to_plot %>%
#       dplyr::select(gene,voom_count,donor,treatment,pool) %>%
#       dplyr::bind_rows(average_baseline_to_plot) %>%
#       dplyr::filter(treatment !="IFN") %>%
#       
#       ggplot(., aes(x=treatment, y=voom_count,col=donor)) +
#       geom_point(position=position_jitter(width=.1,height=0), size=3,alpha=0.6) + 
#       geom_line(data = averages, aes(x = treatment, y = voom_count, color = donor,group = donor), linewidth = 1) +
#       scale_color_manual(values = colors)+
#       theme_bw() + 
#       ggtitle( gen) + 
#       ylab("Normalised log2(counts+0.5)")
#   }
#   
#   plot(patchwork::wrap_plots(plist,ncol = 5,nrow = 6) + patchwork::plot_annotation(title = cont))
# }
# dev.off()
# 
# 
# # just check I'm picking up the interactions correctly per gene
# # the coefficients are not easy to interpret