# Inspecting eQTL results.

library(stringr)
library(tidyverse)
# library(UpSetR) # install
library(patchwork)
library(magrittr)
library(vcfR)
library(data.table)
library(ggrepel)
library(PCAtools)
library(grid)
library(qvalue)
# for interactive plots
library(plotly)
library(htmlwidgets)
library(htmltools)

source("./functions.R")
treatment_cols =  c(untreated = "#8D918B", IFN = "#3A5683", LPS = "#F8766D")

outDir = "../../data/results/4.Inspect_eQTL_results"
dir.create(outDir)
##### Read in ########
# read tensorQTL files
filenames = list.files("../../data/results/tensorqtl", pattern="*window_tensorQTL.txt", 
                       full.names=TRUE,recursive = TRUE)
tensorqtl = lapply(filenames, read.delim)

## TensorQTL - Tidying up names of tested categories
split_file =  str_split(filenames,pattern = "/")
split_filename1 = lapply(X = split_file,FUN = function(x) x[6])
split_filename2 = lapply(X = split_file,FUN = function(x) str_remove(x[7],"sum_sizefactorsNorm_log2_scaled_centered_" ))
split_filename2 = lapply(X = split_filename2,FUN = function(x) str_remove(x[1],"_common_500kb_window_tensorQTL.txt" ))
names(tensorqtl) = paste0(unlist(split_filename1),"_",unlist(split_filename2))

# read scaled expression
filenames = list.files("../../data/for_tensorQTL", pattern="[expr_sum_sizefactorsNorm_log2_scaled_centered_]*bed.gz", 
                       full.names=TRUE,recursive=FALSE)
scaled_expr = lapply(filenames, read_delim)

split_filename =  str_split(filenames,pattern = "/")
split_filename1 = lapply(X = split_filename,FUN = function(x) x[5])
split_filename2 = lapply(X = split_filename1,FUN = function(x) str_remove(x[1],"expr_sum_sizefactorsNorm_log2_scaled_centered_" ))
split_filename2 = lapply(X = split_filename2,FUN = function(x) str_remove(x[1],".bed.gz" ))
names(scaled_expr) = unlist(split_filename2)

# read logNormalised expression - to check abundance of results
# with low expression 

filenames = paste0(paste("../../data/for_tensorQTL/expr_sum_sizefactorsNorm_log2",rep(c("untreated","IFN","LPS"),each=2),
                         rep(c("Not_proliferating","Proliferating"),3),
                         sep = "_" ),".bed")


expression = lapply(filenames, read_delim)

split_filename =  str_split(filenames,pattern = "/")
split_filename1 = lapply(X = split_filename,FUN = function(x) x[5])
split_filename2 = lapply(X = split_filename1,FUN = function(x) str_remove(x[1],"expr_sum_sizefactorsNorm_log2_" ))
split_filename2 = lapply(X = split_filename2,FUN = function(x) str_remove(x[1],".bed" ))
names(expression) = unlist(split_filename2)


# Read metadata

filenames = list.files("../../data/for_tensorQTL", pattern="tensorQTL_metadata_sum_sizefactorsNorm_log2_*", 
                       full.names=TRUE,recursive = TRUE)
metadata = lapply(filenames, read.delim)
split_filename =  str_split(filenames,pattern = "/")
split_filename1 = lapply(X = split_filename,FUN = function(x) x[5])
split_filename2 = lapply(X = split_filename,FUN = function(x) str_remove(x[6],"tensorQTL_metadata_sum_sizefactorsNorm_log2_" ))
split_filename2 = lapply(X = split_filename2,FUN = function(x) str_remove(x[1],".txt" ))
names(metadata) = paste0(unlist(split_filename1),"_",unlist(split_filename2))
rm(split_filename,split_filename1,split_filename2)
metadata = lapply(metadata,t)
metadata = lapply(metadata,as.data.frame)

# In reality we only need to keep a set with a certain number of PCs
# that are the only thing that changes within each treatment and cluster
metadata = metadata[c("73_IFN_Not_proliferating" ,  "84_LPS_Not_proliferating" , "63_untreated_Not_proliferating",
                      "15_LPS_Proliferating"  ,    "15_untreated_Proliferating"  )]

# reading in metadata without PCs - contains more info
filenames = list.files("../../data/for_tensorQTL", pattern="metadata_noPCs_*", 
                       full.names=TRUE,recursive = TRUE)
metadata_general = lapply(filenames, read.delim)
split_filename =  str_split(filenames,pattern = "/")
split_filename1 = lapply(X = split_filename,FUN = function(x) x[5])
split_filename2 = lapply(X = split_filename1,FUN = function(x) str_remove(x[1],"metadata_noPCs_" ))
split_filename2 = lapply(X = split_filename2,FUN = function(x) str_remove(x[1],".txt" ))
names(metadata_general) = unlist(split_filename2)
rm(split_filename,split_filename1,split_filename2)
metadata_general = lapply(metadata_general,as.data.frame)
metadata_general = do.call("rbind",metadata_general)
metadata_general$name = paste(metadata_general$treatment,metadata_general$proliferation_status,sep = "_")
metadata_general$donor_id = str_replace(metadata_general$donor_id,"\\-","\\.")
for(s in names(metadata)){
  
  metadata[[s]] = metadata[[s]] %>%
    dplyr::mutate(PC1 = as.numeric(PC1), PC2 = as.numeric(PC2),
                  genotypePC1=as.numeric(genotypePC1), genotypePC2 = as.numeric(genotypePC2),
                  
                  name = paste(str_split(s,pattern = "_")[[1]][-1],collapse = "_")) %>%
    tibble::rownames_to_column(var = "donor_id") %>%
    dplyr::left_join(.,
                     metadata_general[,c("count","proliferation_status", "donor_id", "treatment", "name")],
                     by = c("donor_id","name")) %>%
    dplyr::rename(ncells = count) %>%
    dplyr::mutate(name = paste(treatment,proliferation_status,sep="_"))
  
  # ncells = as.numeric(ncells),
  # pool = case_when(is.na(pool)==TRUE ~ "shared",
  #                  .default=pool),
  # highlight_shared = case_when(pool=="shared" ~ "shared",
  #                              pool!="shared" ~ "not_shared"))
}
# merge lists and save

metadata = do.call("rbind",metadata)
anyNA(metadata)
readr::write_csv(metadata,paste0(outDir,"/full_metadata.csv"))
# Double-check IPMAR genotypes are being matched to metadata and there are no problems with "-" vs "."

for(s in names(scaled_expr)){
  donors = colnames(scaled_expr[[s]])[!colnames(scaled_expr[[s]]) %in% c("#chr"   ,    "start"   ,   "end"    ,    "gene_id"  )]
  
  expr_subset = expression[[s]] %>%
    dplyr::select(!c("#chr"   ,    "start"   ,   "end"    ,    "gene_id" ,"gene_name"))
  colnames(expr_subset) = donors
  expression[[s]] = expression[[s]] %>%
    dplyr::select(c("#chr"   ,    "start"   ,   "end"    ,    "gene_id" ,"gene_name")) %>%
    dplyr::bind_cols(.,expr_subset)
  
}
# get each dataframe to long format, because they have different donor numbers after filtering
expression = expression %>%
  purrr::imap(~ .x %>%
                dplyr::mutate(
                  name = .y
                ) %>%
                dplyr::relocate(name) %>%
                dplyr::rename("chr" = "#chr")
  ) %>%
  purrr::map( ~ tidyr::pivot_longer(.x,cols = -c(name:gene_name), names_to = "donor_id", values_to = "expression"))

expression = do.call("rbind",expression) %>%
  dplyr::mutate(proliferation_status = dplyr::case_when(grepl("Not",name) ~ "Not_proliferating",
                                                        .default = "Proliferating"),
                treatment = dplyr::case_when(grepl("IFN",name) ~ "IFN",
                                             grepl("LPS",name) ~ "LPS",
                                             grepl("untreated",name) ~ "untreated"))
expression$donor_id = str_replace(expression$donor_id,"\\-","\\.")

readr::write_csv(expression,paste0(outDir,"/full_expression_logNorm.csv"))

gc()

gene_names = expression %>%
  dplyr::select("gene_name","gene_id") %>%
  dplyr::distinct()
# same with the other expression list
# get each dataframe to long format, because they have different donor numbers after filtering
scaled_expr = scaled_expr %>%
  purrr::imap(~ .x %>%
                dplyr::mutate(
                  name = .y
                ) %>%
                dplyr::relocate(name) %>%
                dplyr::rename("chr" = "#chr")
  ) %>%
  purrr::map( ~ tidyr::pivot_longer(.x,cols = -c(name:gene_id), names_to = "donor_id", values_to = "expression"))
scaled_expr = do.call("rbind",scaled_expr) %>%
  dplyr::mutate(proliferation_status = dplyr::case_when(grepl("Not",name) ~ "Not_proliferating",
                                                        .default = "Proliferating"),
                treatment = dplyr::case_when(grepl("IFN",name) ~ "IFN",
                                             grepl("LPS",name) ~ "LPS",
                                             grepl("untreated",name) ~ "untreated")) %>%
  dplyr::left_join(.,gene_names)
readr::write_csv(scaled_expr,paste0(outDir,"/full_expression_scaled.csv"))

gc()


### plot several metadata elements
# for(s in c("35_Not_proliferating_IFN" ,  "35_Not_proliferating_LPS" , "35_Not_proliferating_untreated",
#                 "35_Proliferating_LPS"  ,    "35_Proliferating_untreated"  )){
#   # PCAtools to recalculate PCs for prettier expression PCA plots
#   simple_expression = scaled_expr %>%
#     dplyr::filter(name == s) %>%
#     dplyr::select(donor,expression,gene_id) %>%
#     dplyr::ungroup()
#   simple_expression = simple_expression %>%
#     tidyr::pivot_wider(names_from = donor,values_from = expression)   %>%
#     column_to_rownames(var = "gene_id")
#   
#   metadata_simple = metadata %>%
#     dplyr::filter(name == s)
#   
#   rownames(metadata_simple) = metadata_simple$donor_id
#   
#     
#   pcs = pca(simple_expression, 
#            metadata = metadata_simple, removeVar = NULL)
#   
#   pdf(file =  paste0(outDir,"/PCAtools_expression_pairsplot_",s,".pdf"),width = 10, height = 10)
#   p = pairsplot(pcs,colby = "pool")
#   plot(p)
#   dev.off()
#   
#   pdf(file =  paste0(outDir,"/PCAtools_expression_corrplot_",sample,".pdf"),
#       width = 20, height = 5)
#   # p = eigencorplot(pcs,
#   #              metavars = c('pool','ncells','PC1','genotypePC1')) # PC1 from other PCA method highly correlates
#    p = PCAtools::eigencorplot(pcs,
#                 metavars = c('pool','ncells','genotypePC1'),
#                 components = getComponents(pcs, seq_len(30)),)
#   plot(p)
#   dev.off()
# 
# 
#   
#   # plot expression PCs
#   p1 = ggplot(metadata_simple,aes(x=PC1,y=PC2,label=donor_id,color=pool))+
#     geom_point() + 
#     geom_text_repel()+
#     theme_bw() + 
#     ggtitle("Expression PCs")
#   
#   p2 = ggplot(metadata_simple,aes(x=genotypePC1,y=genotypePC2,label=donor_id,color=pool))+
#     geom_point() + 
#     geom_text_repel()+
#     theme_bw() + 
#     ggtitle("Genotype PCs")
#   
#   p3 = ggplot(metadata_simple,aes(x=PC1,y=PC2,label=donor_id,color=highlight_shared))+
#     geom_point() + 
#     geom_text_repel()+
#     theme_bw() + 
#     ggtitle("Expression PCs")
#   
#   p4 = ggplot(metadata_simple,aes(x=genotypePC1,y=genotypePC2,label=donor_id,color=highlight_shared))+
#     geom_point() + 
#     geom_text_repel()+
#     theme_bw() + 
#     ggtitle("Genotype PCs")
#   
#   pdf(file =  paste0(outDir,"/PC1_PC2_expression_genotype_",sample,".pdf"),width = 15, height = 15)
#   plot((p1 +p3)/ (p2 + p4))
#   dev.off()
#   
# 
#   
#   # number of cells per donor
#   
#   p = ggplot(metadata_simple,aes(x=PC1,y=PC2,label=donor_id,color=log10(ncells)))+
#     geom_point() + 
#     geom_text_repel()+
#     theme_bw() + 
#     ggtitle("Expression PCs")
#   
#   pdf(file =  paste0(outDir,"/PC1_PC2_ncells_",sample,".pdf"),width = 7, height = 5)
#   
#   plot(p)
#   
#   dev.off()
#   
# }
# 

### fix the following
# 
# # explore genotype PCs with 1k genomes
# eur_metadata=read_csv("../../../resources/1000Genomes_eur_ancestry_metadata.csv") %>%
#   dplyr::rename(donor=id)
# 
# genotype_1k_hipsci=read.table("../../data/genotype/plink_genotypes/all_pools.1kgenomes.genotype.MAF05.LD.eigenvec")
# names = ifelse(genotype_1k_hipsci$V1==genotype_1k_hipsci$V2,
#                yes = genotype_1k_hipsci$V1, 
#                no = paste(genotype_1k_hipsci$V1,genotype_1k_hipsci$V2, sep = "_"))
# names = ifelse(grepl("NA|HG",names),
#                yes =unlist(lapply(str_split(names, "_"), function(x) x[2])),
#                no = names)
# genotype_1k_hipsci = genotype_1k_hipsci[c(-1:-2)]
# colnames(genotype_1k_hipsci) = paste0("genotypePC",1:ncol(genotype_1k_hipsci))
# genotype_1k_hipsci$donor = names
# genotype_1k_hipsci$origin = ifelse(grepl("NA|HG",genotype_1k_hipsci$donor),
#                                    yes = "1000 Genomes", no = "hiPSCi")
# genotype_1k_hipsci = merge(genotype_1k_hipsci,eur_metadata,all.x=TRUE, by="donor") %>%
#   mutate(pop = if_else(is.na(pop),true = "hiPSCi",false=pop))
# genotype_1k_hipsci = genotype_1k_hipsci %>%
#   mutate(hipsci_labels = if_else(pop=="hiPSCi",true = donor,false=NA),
#          pop_labels = if_else(!duplicated(pop),true = pop,false = NA))
# 
# p1 = ggplot(genotype_1k_hipsci,aes(x=genotypePC1,y=genotypePC2,label=hipsci_labels,color=origin))+
#   geom_point() + 
#   geom_text_repel()+
#   theme_bw() + 
#   ggtitle("Genotype PCs with 1000 genomes")
# p2 = ggplot(genotype_1k_hipsci,aes(x=genotypePC1,y=genotypePC2,label=pop_labels,color=pop,fill=pop))+
#   geom_point(alpha=0.7) + 
#   geom_text_repel(max.overlaps = 100, col="black",fontface="bold")+
#   theme_bw() + 
#   ggtitle("Genotype PCs with 1000 genomes")
# pdf(file =  paste0(outDir,"/PC1_PC2_genotype_1kgenomes_hipsci.pdf"),width = 10, height = 10)
# plot(p2)
# dev.off()
# 
# pdf(file =  paste0(outDir,"/PC1_PC2_genotype_1kgenomes_hipsci_highlight.pdf"),width = 10, height = 10,
#     units = "in", res = 200
# )
# plot(p1 )
# dev.off()

### Explanation of tensorQTL cis mapping with permutation results
# phenotype_id = ID of the tested molecular phenotype (in this particular case, the gene ID)
# num_var = Number of variants tested in cis for this phenotype
# beta_shape1 = MLE of the shape1 parameter of the Beta distribution
# beta_shape2 = MLE of the shape2 parameter of the Beta distribution
# true_df = effective degrees of freedom the beta distribution approximation
# pval_true_df = empirical P-value for the beta distribution approximation
# variant_id = ID of the best variant found for this molecular phenotypes (i.e. with the smallest p-value)
# tss_distance = Distance between the molecular phenotype - variant pair
# ma_samples =  number of samples carrying the minor allele
# ma_count = total number of minor alleles across individuals (?)
# maf = Minor allele frequency of the variant
# ref_factor =  flag indicating if the alternative allele is the minor allele in the cohort (1 if AF <= 0.5, -1 if not)
# pval_nominal = The nominal p-value of association that quantifies how significant from 0, the regression coefficient is
# slope = The slope associated with the nominal p-value of association [only in version > v2-184]
# pval_perm = A first permutation p-value directly obtained from the permutations with the direct method. This is basically a corrected version of the nominal p-value that accounts for the fact that multiple variants are tested per molecular phenotype.
# pval_beta = A second permutation p-value obtained via beta approximation. We advice to use this one in any downstream analysis.

## Checking the beta approximated p values and permuted p values agree well: sanity check
for(n in names(tensorqtl)){
  message(n)
  plot(tensorqtl[[n]]$pval_perm, tensorqtl[[n]]$pval_beta, xlab="Direct method", ylab="Beta approximation", main=n)
  abline(0, 1, col="red")
  anyDuplicated(tensorqtl[[n]]$phenotype_id) # No duplicated genes
  anyDuplicated(tensorqtl[[n]]$variant_id) # There are duplicated variants (associated with more than one gene)
  
}

## Number of eQTL SNP-gene pairs detected 

## Multiple testing correction for the number of genes genome-wide with the Storey-Tibshirani method
########
for (nam in names(tensorqtl)){
  message(nam)
  has_NA <- FALSE  # Flag variable to track NA values
  
  tryCatch({
    # Check for NAs in pval_beta
    if (any(is.na(tensorqtl[[nam]]$pval_beta))) {
      has_NA <- TRUE
      warning("pval_beta contains NA values for ", n)
    }
    
    # Calculate qval if no NA values
    if (!has_NA) {
      tensorqtl[[nam]]$qval = qvalue(tensorqtl[[nam]]$pval_beta)$qvalues
    }
  }, warning = function(w) {
    message(w)
    has_NA <- TRUE
  })
  
  # remove NA values
  if (has_NA) {
    check = which(is.na(tensorqtl[[nam]]$pval_beta))
    message("\n There are ",length(check)," genes with NA values in pval_beta. Removing those")
    tensorqtl[[nam]] = tensorqtl[[nam]][which(!is.na(tensorqtl[[nam]]$pval_beta)),]
    tensorqtl[[nam]]$qval = qvalue(tensorqtl[[nam]]$pval_beta)$qvalues
    
  }
  tensorqtl[[nam]]$name = paste(str_split(nam,pattern = "_")[[1]][-1],collapse = "_")
}

###### summaries
tensor_results = lapply(names(tensorqtl), function(nam) tensorQTL_summary(QTL = tensorqtl[[nam]], 
                                                                          metadata = metadata[metadata$name %in%  paste(str_split(nam,pattern = "_")[[1]][-1],collapse = "_"),]))
tensor_results = as.data.frame(do.call("cbind",tensor_results))

colnames(tensor_results) = names(tensorqtl)

tensor_results = tensor_results %>% 
  rownames_to_column(var="variables") %T>%
  write_csv(paste0(outDir,"/tensorQTL_summary.csv"))

gc()
# number of eGenes / eQTLs per fitted PCs

toplot =  tensor_results %>% 
  pivot_longer(cols = !variables,
               names_to = "sample", values_to = "number") %>%
  mutate(nPCs = as.double(str_split_i(sample,pattern = "_",i = 1)),
         cluster = case_when(grepl("Not_proliferating",sample)==TRUE ~ "Not_proliferating",
                             grepl("Proliferating",sample)==TRUE ~ "Proliferating"),
         treatment = case_when(grepl("IFN",sample)==TRUE ~ "IFN",
                               grepl("LPS",sample)==TRUE ~ "LPS",
                               grepl("untreated",sample)==TRUE ~ "untreated")
  ) %>%
  dplyr::filter(variables %in% c("n_genes_05")) %>%
  dplyr::mutate(variables = case_when(variables == "n_genes_05" ~ "Significant eGenes 5%",
                                      .default = variables),
                cluster = case_when(cluster == "Not_proliferating" ~ "Not proliferating",
                                    .default = cluster)) %>%
  dplyr::group_by(treatment,cluster) %>%
  dplyr::mutate(colors = if_else(number == max(number), true = treatment, false = "blank"),
                is_max = if_else(number == max(number), true = "max", false = "not_max")) %>%
  dplyr::ungroup() %>%
  dplyr::filter(nPCs %in% c(seq(5,60,by=5),63,67,73,77,84,seq(90,120,by=5))) ### added this after checking the full table of results so that 
# the plot is not so busy
p1 = toplot %>%
  dplyr::filter(cluster == "Not proliferating") %>%
  ggplot( aes(x=nPCs, y=number, col = treatment, fill = colors)) + 
  geom_point(aes(shape=is_max),size=3,stroke = 1,position = position_jitter(h = 0.3, w = 0))+
  scale_shape_manual(values = c(21,23))+
  theme_bw()+
  ylab("Number of eGenes")+
  xlab("Number of expression PCs fitted")+
  ggtitle("Not proliferating") +
  scale_color_manual(values=treatment_cols) +
  scale_fill_manual(values=c(treatment_cols,"blank" ="white")) +
  guides(fill = "none",shape="none") +
  theme(strip.background = element_rect(fill="white")) +
  geom_vline(xintercept = 63, col = "#8D918B", linetype = "dashed") + 
  geom_vline(xintercept = 73, col = "#3A5683", linetype = "dashed") + 
  geom_vline(xintercept = 84, col = "#F8766D", linetype = "dashed") 

p2 = toplot %>%
  dplyr::filter(cluster == "Proliferating") %>%
  ggplot( aes(x=nPCs, y=number, col = treatment, fill = colors)) + 
  geom_point(aes(shape=is_max),size=3,stroke = 1,position = position_jitter(h = 0.3, w = 0))+
  scale_shape_manual(values = c(21,23))+
  theme_bw()+
  ylab("Number of eGenes")+
  xlab("Number of expression PCs fitted")+
  ggtitle("Proliferating") +
  scale_color_manual(values=treatment_cols) +
  scale_fill_manual(values=c(treatment_cols,"blank" ="white")) +
  guides(fill = "none",shape="none", col = "none") +
  theme(strip.background = element_rect(fill="white")) +
  geom_vline(xintercept = 15, col = "#8D918B", linetype = "dashed") + 
  geom_vline(xintercept = 15, col = "#F8766D", linetype = "dashed") 


pdf(file =  paste0(outDir,"/eQTLs_per_PCs.pdf"),width = 7, height = 5)
plot(p1 / p2 + patchwork::plot_layout(guides='collect') &
       theme(legend.position='bottom'))
dev.off()

# eQTL vs mean/median/CV number of cells per cluster
# check number of donors per cluster and eQTL as well
to_plot = tensor_results %>% 
  pivot_longer(cols = !variables,
               names_to = "sample", values_to = "number") %>%
  mutate(nPCs= factor(str_split_i(sample,pattern = "_",i = 1),levels = c(seq(5,60,by=5),seq(61,94,by=1),seq(95,120,by=5)),ordered = TRUE),
         cluster = case_when(grepl("Not_proliferating",sample)==TRUE ~ "Not_proliferating",
                             grepl("Proliferating",sample)==TRUE ~ "Proliferating"),
         treatment = case_when(grepl("IFN",sample)==TRUE ~ "IFN",
                               grepl("LPS",sample)==TRUE ~ "LPS",
                               grepl("untreated",sample)==TRUE ~ "untreated") ) %>%
  dplyr::filter(nPCs ==65) %>%
  pivot_wider(id_cols = sample,names_from = variables,values_from = number)


pdf(file =  paste0(outDir,"/eGenes_mean_median_cells_per_cluster_vs_eQTL_donors_CV_qval_05.pdf"),width = 12, height = 8)
p1 = ggplot(to_plot, aes(x=mean_n_cells, y=n_genes_05, col = sample)) + 
  geom_point(aes(size=n_donors)) + 
  ylim(min(to_plot$n_genes_05)-100,max(to_plot$n_genes_05)+100) +
  xlab("Mean cells across donors per cluster") +
  ylab("N eGenes (q value < 0.05)") +
  labs(size = "Number of donors")  +
  theme_classic() 
p2 =  ggplot(to_plot, aes(x=mean_n_cells, y=n_genes_05,col = sample)) + 
  geom_point(aes(size=CV_cells)) + 
  ylim(min(to_plot$n_genes_05)-100,max(to_plot$n_genes_05)+100) +
  xlab("Mean cells across donors per cluster") +
  ylab("N eGenes (q value < 0.05)") +
  labs(size = "Coef. of variation\n of N cells")  +
  theme_classic() 

p3 = ggplot(to_plot, aes(x=median_n_cells, y=n_genes_05, col = sample)) + 
  geom_point(aes(size=n_donors)) + 
  ylim(min(to_plot$n_genes_05)-100,max(to_plot$n_genes_05)+100) +
  xlab("Median cells across donors per cluster") +
  ylab("N eGenes (q value < 0.05)") +
  labs(size = "Number of donors") +
  theme_classic() 

p4 = ggplot(to_plot, aes(x=median_n_cells, y=n_genes_05, col = sample)) + 
  geom_point(aes(size=CV_cells)) + 
  ylim(min(to_plot$n_genes_05)-100,max(to_plot$n_genes_05)+100) +
  xlab("Median cells across donors per cluster") +
  ylab("N eGenes (q value < 0.05)") +
  labs(size = "Coef. of variation\n of N cells") +
  theme_classic() 


(p1 + p2) / (p3 + p4)
dev.off()

pdf(file =  paste0(outDir,"/eQTLs_mean_median_cells_per_cluster_vs_eQTL_donors_CV_q_val_05.pdf"),width = 12, height = 8)
p1 = ggplot(to_plot, aes(x=mean_n_cells, y=n_variants_05, col = sample)) + 
  geom_point(aes(size=n_donors)) + 
  ylim(min(to_plot$n_variants_05)-100,max(to_plot$n_variants_05)+100) +
  xlab("Mean cells across donors per cluster") +
  ylab("N eQTLs (q value < 0.05)") +
  labs(size = "Number of donors")  +
  theme_classic() 
p2 =  ggplot(to_plot, aes(x=mean_n_cells, y=n_variants_05, col = sample)) + 
  geom_point(aes(size=CV_cells)) + 
  ylim(min(to_plot$n_variants_05)-100,max(to_plot$n_variants_05)+100) +
  xlab("Mean cells across donors per cluster") +
  ylab("N eQTLs (q value < 0.05)") +
  labs(size = "Coef. of variation\n of N cells") +
  theme_classic() 

p3 = ggplot(to_plot, aes(x=median_n_cells, y=n_variants_05, col = sample)) + 
  geom_point(aes(size=n_donors)) + 
  ylim(min(to_plot$n_variants_05)-100,max(to_plot$n_variants_05)+100) +
  xlab("Median cells across donors per cluster") +
  ylab("N eQTLs (q value < 0.05)") +
  labs(size = "Number of donors") +
  theme_classic() 

p4 = ggplot(to_plot, aes(x=median_n_cells, y=n_variants_05, col = sample)) + 
  geom_point(aes(size=CV_cells)) + 
  ylim(min(to_plot$n_variants_05)-100,max(to_plot$n_variants_05)+100) +
  xlab("Median cells across donors per cluster") +
  ylab("N eQTLs (q value < 0.05)") +
  labs(size = "Coef. of variation\n of N cells") +
  theme_classic() 

(p1 + p2) / (p3 + p4)
dev.off()

## eQTL shared between groups

# Bind all info together for best number of PCs

for(n in names(tensorqtl)){
  tensorqtl[[n]]$group = n # there's a "sample" with different meaning in eQTL results already
}

tensorqtl_all = do.call("rbind",args = tensorqtl)
tensorqtl_all = tensorqtl_all %>% 
  dplyr::mutate(nPCs= as.double(str_split_i(group,pattern = "_",i = 1)),
                cluster = case_when(grepl("Not_proliferating",group)==TRUE ~ "Not_proliferating",
                                    grepl("Proliferating",group)==TRUE ~ "Proliferating"),
                treatment = case_when(grepl("IFN",group)==TRUE ~ "IFN",
                                      grepl("LPS",group)==TRUE ~ "LPS",
                                      grepl("untreated",group)==TRUE ~ "untreated")
  ) %>%
  dplyr::mutate(toselect = paste(nPCs,treatment,cluster,sep = "-")) %>%
  dplyr::filter(toselect %in% c("63-untreated-Not_proliferating",
                                "73-IFN-Not_proliferating",
                                "84-LPS-Not_proliferating",
                                "15-untreated-Proliferating",
                                "15-LPS-Proliferating"))

# map ensembl gene ids to symbol

tensorqtl_all = tensorqtl_all %>%
  dplyr::rename(gene_id = phenotype_id) %>%
  dplyr::left_join(x = .,y = gene_names, by="gene_id") %>%
  dplyr::relocate(gene_id,gene_name)

# how many significant eQTLs across conditions for selected PCs
eGenes = tensorqtl_all %>%
  dplyr::filter( qval<0.05) %>%
  dplyr::select(gene_name) %>%
  distinct()
nrow(eGenes) # 7121

# only non-prolif
eGenes = tensorqtl_all %>%
  dplyr::filter( qval<0.05 & cluster == "Not_proliferating") %>%
  dplyr::select(gene_name) %>%
  distinct()
nrow(eGenes) # 7068

eQTLs = tensorqtl_all %>%
  dplyr::filter( qval<0.05) %>%
  dplyr::select(variant_id) %>%
  distinct()
nrow(eQTLs) # 13205
# not proliferating:
eQTLs = tensorqtl_all %>%
  dplyr::filter( qval<0.05 & cluster == "Not_proliferating") %>%
  dplyr::select(variant_id) %>%
  distinct()
nrow(eQTLs) # 12230

# add position column

split_pos = str_split(tensorqtl_all$variant_id,pattern = "_")
tensorqtl_all$position = unlist(lapply(X = split_pos,FUN = function(x) paste0(x[1],"_",x[2])))

## read in vcf file 
vcf = vcfR::read.vcfR("../../data/genotype/tensor.vcf.gz")
# pre-filtered VCF file to tested eQTL variants so that this step
# and the next don't take so long
ref_alt_donor = make_dosages_ref_alt(vcf) 

ref_alt_donor = as_tibble(ref_alt_donor) %>%
  relocate(c("names" ,"position")) %>%
  mutate(names = gsub("chr","", names),
         position = gsub("chr","", position))
colnames(ref_alt_donor) = str_replace(colnames(ref_alt_donor),"\\-","\\.")


ref_alt_donor = ref_alt_donor %>%
  dplyr::filter(position %in% tensorqtl_all$position)
tensorqtl_all = tensorqtl_all %>%
  dplyr::filter(position %in% ref_alt_donor$position) %>%
  dplyr::arrange(qval)
# save this file
write_rds(ref_alt_donor,paste0(outDir,"/vcf_tibble_sign_tensorQTL_REF_ALT.rds"))

# save other important files
write_rds(scaled_expr,paste0(outDir,"/scaled_expr.rds"))
write_rds(expression,paste0(outDir,"/expr.rds"))
write_rds(metadata,paste0(outDir,"/metadata.rds"))


###### can reload big files from here #####
ref_alt_donor = read_rds(paste0(outDir,"/vcf_tibble_sign_tensorQTL_REF_ALT.rds"))
scaled_expr = read_rds(paste0(outDir,"/scaled_expr.rds"))
expression = read_rds(paste0(outDir,"/expr.rds"))
metadata = read_rds(paste0(outDir,"/metadata.rds"))

table(ref_alt_donor$names %in% tensorqtl_all$variant_id) # FALSE are potential multiallelic variants?
table(tensorqtl_all$variant_id  %in% ref_alt_donor$names) # FALSE are potential ALT REF swaps
# PLINK by default puts the minor allele on dosage 2
# so variants that where the ALT allele in the VCF is more abundant in the population will be swapped
# save on the tensorQTL result file the variants that are swapped with the VCF
tensorqtl_all = tensorqtl_all %>%
  mutate(swapped_tensor_alleles = unlist(map(variant_id,swap_REF_ALT_alelle_names))) %>%
  mutate(swapped_tensor_in_VCF = case_when(swapped_tensor_alleles %in% ref_alt_donor$names == TRUE ~ TRUE,
                                           swapped_tensor_alleles %in% ref_alt_donor$names == FALSE ~ FALSE)) %>%
  mutate(swapped_slope = case_when(swapped_tensor_in_VCF == TRUE ~ -slope,
                                   swapped_tensor_in_VCF == FALSE ~ NA)) %>%
  mutate(gene_variant_pair = paste0(gene_name,"-",variant_id))  %T>% # create gene-variant pairs
  write_csv(paste0(outDir,"/tensorQTL_variant_gene.csv"))

#### can reload this file from here
tensorqtl_all = read_csv(paste0(outDir,"/tensorQTL_variant_gene.csv"))

prop.table(table(tensorqtl_all$swapped_tensor_in_VCF)) * 100 
#3% swapped alleles
# be careful with the slopes when reporting these effects with respect to REF (need to swap signs)
prop.table(table(tensorqtl_all$variant_id %in% ref_alt_donor$names )) * 100 

# check top 9 results per group
for (g in unique(tensorqtl_all$group)){
  plist = list()
  for(n in 1:9){
    subset = tensorqtl_all %>%
      relocate(c("gene_name" ,"position")) %>%
      dplyr::filter(group == g) %>%
      dplyr::arrange(qval) %>%
      dplyr::slice(n)
    
    # Try-catch block to check the number of rows in my_tibble
    tryCatch({
      # Check the number of rows in my_tibble
      if(nrow(subset) > 1) {
        stop("Error: df has more than one row")
      }
    }, error = function(e) {
      # Handle the error by printing the error message
      message(e)
    })
    
    expression_subset = expression %>%
      dplyr::filter(name ==unique(subset$name))
    metadata_subset = metadata %>%
      dplyr::filter(name ==unique(subset$name))
    plist[[n]] = boxplot_eQTL(variant = subset$variant_id,
                              gene = subset$gene_name,
                              group = str_remove(g,pattern = "35_"),
                              genotype_df = ref_alt_donor,
                              expression_df = expression_subset,
                              metadata_df = metadata_subset,
                              color_by_pool=FALSE # keep this false for now
    )
  }
  p = wrap_plots(plist, ncol = 3) + plot_layout(guides = "collect")
  
  
  pdf(file =  paste0(outDir,"/eQTLs_top_results_",g,".pdf"),
      width = 15, height = 15
  )
  plot(p)
  dev.off()
  
}


# upsetplot shared / specific eQTLs - gene pairs

tensorqtl_all = read_csv(paste0(outDir,"/tensorQTL_variant_gene.csv"))
upset_tensor = tensorqtl_all %>% 
  dplyr::filter(qval < 0.05) %>% # filter for significance
  group_by(cluster, treatment) %>% 
  summarize(gene_variant_pair = list(setNames(gene_variant_pair, paste0(cluster, "_", treatment)))) %>% # format for upsetR
  ungroup() %>%
  dplyr::select(gene_variant_pair)%>%
  .$gene_variant_pair

names(upset_tensor) = sapply(upset_tensor, function(x) names(x[1])) # name every element of the list
upset_tensor = lapply(upset_tensor, function(x) { # remove names within lists
  names(x) <- NULL
  x
})

UpSetR::upset(fromList(upset_tensor), order.by = "freq")
# most unique, but likely shared signal because of LD - need to use coloc/LD map to assess
# because there are many more shared genes (but not that many?):
# checking again selecting only genes

upset_tensor = tensorqtl_all %>% 
  dplyr::filter(qval < 0.05) %>% # filter for significance
  group_by(cluster, treatment) %>% 
  summarize(gene_name = list(setNames(gene_name, paste0(cluster, "_", treatment)))) %>% # format for upsetR
  ungroup() %>%
  dplyr::select(gene_name)%>%
  .$gene_name

names(upset_tensor) = sapply(upset_tensor, function(x) names(x[1])) # name every element of the list
upset_tensor = lapply(upset_tensor, function(x) { # remove names within lists
  names(x) <- NULL
  x
})



# Specific boxplots - fix
totest = tensorqtl_all %>%
  dplyr::filter(group ==  "73_IFN_Not_proliferating", 
                gene_name == "TREM2")

res = boxplot_eQTL(variant = totest$variant_id,
                   gene = totest$gene_name,
                   group = "untreated_Not_proliferating",
                   genotype_df = ref_alt_donor,
                   expression_df = scaled_expr,
                   metadata_df = metadata,
                   color_by_pool=FALSE
)
pdf(file =  paste0(outDir,"/TREM2_tophit_73_Not_proliferating_IFN.pdf"),
    width = 6, height = 5
)


plot(res)
dev.off()
