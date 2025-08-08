library(tidyverse)
library(qvalue)
library(pheatmap)
source("./functions.R")
set.seed(123)

external_eqtl_datasets = c("macromap_Panousis_2023_Ctrl_6","macromap_Panousis_2023_Ctrl_24",
                           "macromap_Panousis_2023_IFNG_6","macromap_Panousis_2023_IFNG_24",
                           "macromap_Panousis_2023_sLPS_6","macromap_Panousis_2023_sLPS_24")
for(dataset in external_eqtl_datasets){
outDir = paste0("../../data/results/8.colocalisation_analysis/coloc_results/",dataset)
name = dataset
external_eQTL_folder = as.character(paste0("../../../resources/macromap_Panousis_2023/eQTL/",dataset)) 

external_eQTL_lead_path = list.files(external_eQTL_folder,pattern = "*permuted.tsv.gz",full.names = TRUE) 

if(!is_empty(external_eQTL_lead_path)){
  external_eQTL_variants_path = list.files(external_eQTL_folder,pattern = "*permuted.tsv.gz",full.names = TRUE) 
  
  
  external_eQTL_lead =  readr::read_tsv(external_eQTL_lead_path) %>%
    dplyr::rename(variant_id = variant,chr = chromosome,snp_pos = position) %>%
    # omit genes without significant eQTLs
    dplyr::filter(!is.na(variant_id)) %>%
    dplyr::filter(!is.na(p_beta)) %>%
    # calculating q values
    dplyr::mutate(qval =  qvalue(p_beta)$qvalues) %>%
    dplyr::filter(qval <0.05) %>%
    dplyr::rename(gene_id = ensembl_id)
} else{
  external_eQTL_path = list.files(external_eQTL_folder,pattern = "*.all.tsv.gz",full.names = TRUE) 
  # exclude index, if present
  external_eQTL_path = external_eQTL_path[!grepl("tbi",external_eQTL_path)]
  
  external_eQTL_lead = readr::read_tsv(external_eQTL_path) %>%
    dplyr::rename(standard_error = if_else(condition = "se" %in% colnames(.),true = "se", false = "standard_error")) %>%
    dplyr::rename(ref = REF) %>%
    dplyr::mutate(variant = paste(chromosome,position,ref,alt,sep = "_")) %>%
    dplyr::select(molecular_trait_id,variant, beta,pvalue,alt, chromosome,position,
                  standard_error,maf,significant_by_2step_FDR) %>%
    dplyr::rename(variant_id = variant,
                  pval_nominal = pvalue,
                  effect_allele = alt,
                  gene_id = molecular_trait_id) %>%
    dplyr::mutate(chr = chromosome,
                  pos=position, # ensuring types are correct
                  varbeta = standard_error^2) %>%
    # no need to adjust for N, see https://github.com/chr1swallace/coloc/issues/14
    
    dplyr::select(gene_id,chr, pos,variant_id,effect_allele, beta,varbeta,pval_nominal,significant_by_2step_FDR) %>%
    dplyr::filter(significant_by_2step_FDR == "Yes")
  
  external_eQTL_lead = external_eQTL_lead %>%
    dplyr::group_by(gene_id) %>%
    dplyr::arrange(pval_nominal) %>%
    dplyr::slice_head(n=1)  %>% # take the lead variant per eGene as in tensorQTL
    dplyr::ungroup()
  

}

eqtl = readr::read_csv("../../data/results/4.Inspect_eQTL_results/tensorQTL_variant_gene_60PCs.csv")
# extract significant results
signif = eqtl %>%
  dplyr::filter(group %in% c(paste0("60_",c("untreated","IFN","LPS"),"_Not_proliferating"))) %>%
  dplyr::filter(qval<0.05)



message("There are ",length(unique(external_eQTL_lead$gene_id)), " significant eGenes in ",name)
message("There are ",length(unique(signif$gene_id)), " significant eGenes in our data") # 5127

# measuring proportion of genes in their data that are present in ours
# over the total of all THEIR signf eGenes
external_eQTL_lead = external_eQTL_lead %>%
  dplyr::mutate(shared_eGenes_all = gene_id %in% signif$gene_id,
                shared_eGenes_untreated = gene_id %in% signif[signif$treatment=="untreated",]$gene_id,
                shared_eGenes_IFN = gene_id %in% signif[signif$treatment=="IFN",]$gene_id,
                shared_eGenes_LPS = gene_id %in% signif[signif$treatment=="LPS",]$gene_id)

pdf(paste0(outDir,"/shared_primary_eQTL_ours_",name,".pdf"),width = 6, height = 5)
p = external_eQTL_lead %>%
  dplyr::select(gene_id,shared_eGenes_all,shared_eGenes_untreated,shared_eGenes_IFN,shared_eGenes_LPS) %>%
  dplyr::rename(all =shared_eGenes_all,untreated = shared_eGenes_untreated,
                IFN = shared_eGenes_IFN,LPS = shared_eGenes_LPS) %>%
  tidyr::pivot_longer(cols = -gene_id, names_to = "category",values_to = "shared eGene") %>%
  dplyr::mutate(category = factor(category, levels = c("all","untreated","IFN","LPS"), ordered = TRUE)) %>%
  ggplot(aes(x = category,fill = `shared eGene`)) +
  scale_fill_brewer(palette = "Set1") +
  geom_bar() +
  theme_minimal() + 
  scale_x_discrete(labels = c("all" = paste0("all (",
   round(sum(external_eQTL_lead$shared_eGenes_all) / sum(table(external_eQTL_lead$shared_eGenes_all)),2)*100 ,"%)"), 
  "IFN" = paste0("IFN (",
  round(sum(external_eQTL_lead$shared_eGenes_IFN) / sum(table(external_eQTL_lead$shared_eGenes_IFN)),2)*100 ,"%)"),
  "LPS" = paste0("LPS (",
  round(sum(external_eQTL_lead$shared_eGenes_LPS) / sum(table(external_eQTL_lead$shared_eGenes_LPS)),2)*100 ,"%)"),
  "untreated" = paste0("untreated (",
  round(sum(external_eQTL_lead$shared_eGenes_untreated) / sum(table(external_eQTL_lead$shared_eGenes_untreated)),2)*100 ,"%)"))) +
  xlab("treatment (% shared)") +
  ggtitle("Shared eGenes with ", name)
plot(p)

dev.off()


round(table(external_eQTL_lead$shared_eGenes_all) / sum(table(external_eQTL_lead$shared_eGenes_all)),2) # FALSE  TRUE 
# 0.59  0.41.  Young
# 0.45 0.55 Fujita
round(table(external_eQTL_lead$shared_eGenes_untreated) / sum(table(external_eQTL_lead$shared_eGenes_untreated)),2)# FALSE  TRUE 
# 0.67  0.33 Young
#  0.54  0.46 Fujita
round(table(external_eQTL_lead$shared_eGenes_IFN) / sum(table(external_eQTL_lead$shared_eGenes_IFN)),2) # 0.7. 0.3 Young
# 0.63 0.37 Fujita
round(table(external_eQTL_lead$shared_eGenes_LPS) / sum(table(external_eQTL_lead$shared_eGenes_LPS)),2) # 0.67. 0.33 Young
# 0.59. 0.41 Fujita


### checking the oposite case
# measuring proportion of genes in their data that are present in ours
# over the total of all THEIR signf eGenes
signif = signif %>%
  dplyr::mutate(shared_eGenes = gene_id %in% external_eQTL_lead$gene_id)
table(signif$shared_eGenes)
 

}

#### how many of non-shared eGenes are not expressed in my data? #####
# reading in expression after eQTL filters
my_expression = list()
for(treat in c("untreated","IFN","LPS")){
  my_expression[[treat]] = read_tsv(paste0("../../data/for_tensorQTL/expr_sum_sizefactorsNorm_log2_",treat,
                                    "_Not_proliferating.bed")) %>%
    # pivot to long and average per gene
    tidyr::pivot_longer(cols = c(-`#chr`, -start, -end, -gene_id, -gene_name),
                        names_to = "sample", values_to = "log2_norm") %>%
    dplyr::mutate(treatment = treat) %>%
    dplyr::group_by(gene_id,gene_name,treatment) %>%
    dplyr::reframe(median_log2_norm = median(log2_norm)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(gene_id,.keep_all = TRUE)
  
}

my_expression = do.call("rbind",my_expression)

# pivot wider and fill NAs with 0

my_expression = my_expression %>%
  dplyr::arrange(gene_id) %>%
  tidyr::pivot_wider(names_from = treatment, values_from = median_log2_norm,
                     values_fill = 0)

# check external eGenes not in my dataset
not_there = external_eQTL_lead$gene_id[!external_eQTL_lead$shared_eGenes_all]

toplot = my_expression %>%
  dplyr::filter(gene_id %in% not_there) %>%
  dplyr::select(-gene_id) %>%
  dplyr::distinct(gene_name,.keep_all = TRUE) %>%
  tibble::column_to_rownames("gene_name") %>%
  as.matrix()



pdf(paste0(outDir,"/eGenes_from_",name,"_not_in_ours_log2_median_expr.pdf"),width = 6, height = 27)

pheatmap(toplot,cluster_cols = TRUE,fontsize = 10)
dev.off()
  

