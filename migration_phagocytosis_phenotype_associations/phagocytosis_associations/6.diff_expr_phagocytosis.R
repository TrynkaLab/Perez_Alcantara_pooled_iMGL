# Diff expr vs phagocytosis

# differential expression with PRS as continuous variable
# Diff expr phagocytosis
# Differential expression with limma+voom

library(limma)
library(edgeR)
library(tidyverse) # For ggplot2 and easy manipulation of data
library(patchwork) # To combine plots
library(PCAtools)
library(variancePartition)
set.seed(123)
source("./functions.R")
outdir="../../../data/results/phagocytosis/6.Diff_expr_limma_phagocytosis"
dir.create(outdir, recursive = T)

line_prop_changes_path = "../../../data/results/phagocytosis/1.check_line_proportions/line_prop_changes_per_well.csv"
genotype_pc_path = "../../../../OTAR2065_sc_eQTL/data/genotype/plink_genotypes/all_pools.all_donors.genotype.MAF05.eigenvec"

gene_info = read_csv("/lustre/scratch123/hgi/teams/trynka/resources/biomart/Homo_sapiens.GRCh38.111.genes.csv") %>%
  dplyr::rename(gene = gene_name) %>%
  dplyr::filter(gene_biotype == "protein_coding")
### treatment with prs interaction
pseudobulk =readr::read_tsv("../../../../OTAR2065_differentiation_efficiency/data/results/5.prepare_differential_expression_pseudobulks/combined_pseudobulk_raw_counts_treatmentxprolifxidxpool.txt") %>%
  dplyr::filter(gene %in% gene_info$gene) %>%
  tibble::column_to_rownames(var = "gene") 


full_genotype_pcs = read.table(genotype_pc_path)
rownames(full_genotype_pcs) = ifelse(
  full_genotype_pcs$V1 == full_genotype_pcs$V2,
  yes = full_genotype_pcs$V1,
  no = paste(full_genotype_pcs$V1, full_genotype_pcs$V2, sep = "_")
)
full_genotype_pcs = full_genotype_pcs[c(2:7)] %>%
  dplyr::rename(
    line = V2,
    genotypePC1 = V3,
    genotypePC2 = V4,
    genotypePC3 = V5,
    genotypePC4 = V6,
    genotypePC5 = V7
  )

metadata = readr::read_tsv("../../../../OTAR2065_differentiation_efficiency/data/results/5.prepare_differential_expression_pseudobulks/metadata_treatmentxprolifxidxpool.txt") %>%
  tibble::column_to_rownames(var = "cols_names") %>%
  dplyr::rename(ncells = count,line= donor_id) %>%
  dplyr::mutate(pseudobulk_colnames=paste(treatment,proliferation_status,line,pool,sep = "_"))

# subset to non-proliferating only
metadata = metadata %>%
  dplyr::filter(proliferation_status == "Not_proliferating") 

pseudobulk = pseudobulk %>%
  dplyr::select(metadata$pseudobulk_colnames)
length(unique(metadata$line)) # 250

line_prop_changes = read.csv(line_prop_changes_path)

message("There are ", length(unique(line_prop_changes$line)), " lines in the analysis") # 202

line_prop_changes_summary = line_prop_changes %>%
  dplyr::filter(prop_unadjusted_min_value > 0.005) %>%
  dplyr::filter(!log_fraction_mean %in% c(NA,Inf,-Inf)) %>%
  dplyr::filter(!is.na(log_fraction_mean)) %>%
  dplyr::select(log_fraction_mean,replicate,line,condition,treatment,pool,sex,prop_unadjusted_min_value) %>%
  dplyr::group_by(line,treatment,pool,sex) %>%
  dplyr::summarise(log_fraction_mean_per_pool = mean(log_fraction_mean,na.rm = TRUE), prop_unadjusted_min_value_mean = mean(prop_unadjusted_min_value,na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  distinct() 

message("There are ", length(unique(line_prop_changes_summary$line)), " lines in the analysis") # 159

metadata = metadata %>%
  dplyr::left_join(line_prop_changes_summary) %>%
  dplyr::left_join(full_genotype_pcs) %>%
  distinct()
# Age will be partly confounded with PRS and sex most likely

# remove those that have NA values for the phagocytosis
metadata = metadata %>%
  dplyr::filter(!is.na(log_fraction_mean_per_pool))
message("There are ", length(unique(metadata$line)), " lines in the analysis") # 158

pseudobulk = pseudobulk %>%
  dplyr::select(match(metadata$pseudobulk_colnames,colnames(.)))
# same order
identical(colnames(pseudobulk),metadata$pseudobulk_colnames)

##### filter by number of cells


discarded = metadata$ncells < 100
table(discarded)

pseudobulk = pseudobulk[,!discarded]
metadata = metadata[!discarded,]
length(unique(metadata$line)) # 154/261, 59%

# analyse separately per treatment ########

for(treat in c("untreated","IFN","LPS")){
  
  keep = metadata$treatment == treat
  table(keep)
### creating DGElist object
  
  pseudobulk2 = pseudobulk[,keep]
  metadata2 = metadata[keep,]

dge = edgeR::DGEList(counts=pseudobulk2)



# remove genes not expressed in 30% samples (min 0.1 CPM) - checking very small effects
isexpr = rowSums(edgeR::cpm(dge) > 0.1) >= floor(0.3*ncol(dge))
message(paste0("There are ", sum(isexpr), " genes with over 0.1 count per million in 30% of samples"))
# ~13000 genes
dge = dge[isexpr,,keep.lib.sizes=FALSE]  
dim(dge)

## how many have 0 counts
table(rowSums(dge$counts)==0)
# none

# TMM normalization:
dge = edgeR::calcNormFactors(dge)

  ### create design matrix
  
  metadata2 = metadata2 %>%
    dplyr::filter(!is.na(genotypePC1)) # clones were not included in the genotype PC calculations
  
  metadata2$sex = as.factor(metadata2$sex)
  metadata2$sex = relevel(metadata2$sex,ref = "Female")
  
  # subset v to the remaining elements of the metadata
  dge=dge[,metadata2$pseudobulk_colnames]
  
  
  # ncells is correlated with prop_unadjusted_min_value_mean
  mat = model.matrix(~ log_fraction_mean_per_pool + sex + prop_unadjusted_min_value_mean + genotypePC1 + genotypePC2,metadata2)
  mat = as.data.frame(mat)
  rownames(mat) = metadata2$pseudobulk_colnames
  qr(mat)$rank
  ncol(mat)
  is.fullrank(mat)
  colnames(mat)
  mat = mat[,colSums(mat)!=0]
  colnames(mat)
  is.fullrank(mat)
  
  colnames(mat) = make.names(colnames(mat))
  
  v = voom(dge, design = mat, plot=TRUE)
  
  cor = duplicateCorrelation(v, design=mat, block = metadata2$pool)
  cor$consensus # within-pool correlation
  # voom weights may have changed:
  v2 = voom(dge, design = mat, plot=TRUE, block = metadata2$pool, correlation = cor$consensus)
  cor = duplicateCorrelation(v2, design = mat, block =  metadata2$pool) # extract duplicates again
  cor$consensus
  
  # fit 
  message("Fitting")
  fit = lmFit(v2, design=mat,block = metadata2$pool, correlation  = cor$consensus)
  
  cont = c( "log_fraction_phago" = "log_fraction_mean_per_pool") # changes in gene expression per unit increase in phagocytosis, per treatment 
  
  
  contrast.matrix = makeContrasts(contrasts = cont,
                                  levels=mat)
  colnames(contrast.matrix) = c("log_fraction_phago")
  
  fit2 = contrasts.fit(fit, contrast.matrix)
  fit2 = eBayes(fit2)
  
  summary_additive = summary(decideTests(fit2,lfc = 0.5, # effects seem to be small
                                         adjust.method = "fdr",p.value=0.05))
  print(summary_additive)
  message("Extracting results")
  
  res_additive = list()
  for(cont in colnames(contrast.matrix)){
    res_additive[[cont]] = topTable(fit2, coef=cont,
                                    number = length(fit2$coefficients)) %>%
      dplyr::arrange(adj.P.Val) %>%
      tibble::rownames_to_column(var = "symbol")
    
    write_csv(res_additive[[cont]],paste0(outdir,
                                          "/DiffExpr_",cont,"_",treat,".csv"))
  }
  saveRDS(res_additive,paste0(outdir,"/additive_res_",treat,".rds"))
  
  ### remove pool batch effects for plotting
  # similar to partial residuals
  pseudobulk_corrected = limma::removeBatchEffect(v$E, batch = metadata2$pool, design = mat)
  saveRDS(pseudobulk_corrected,paste0(outdir,"/pseudobulk_corrected_",treat,".rds"))
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
        dplyr::left_join(.,metadata2)
      
      # slope of limma fit don't match corrected counts, doing fit on the fly
      
      p = ggplot(gene_to_plot, aes(x=log_fraction_mean_per_pool, y=voom_count)) +
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
