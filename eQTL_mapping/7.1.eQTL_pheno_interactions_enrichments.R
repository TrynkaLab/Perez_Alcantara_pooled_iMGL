# enrichments eQTL-pheno interactions


library(tidyverse)
library(patchwork)
library(fgsea)
library(msigdbr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
source("./functions.R")
options(future.globals.maxSize = 40000 * 1024^2) # 40Gb


output_dir = "../../data/results/3.run_LMM_pheno_interactions"

##### data from 3.run_LMM_pheno_interactions.R
prolif = "Not_proliferating"
all_my_genes = list()
# for(pheno in c("phagocytosis","migration")){
for(pheno in c("phagocytosis")){
  for(treat in c("untreated","IFN","LPS") ){
    message("Working on ",treat)
    
    my_fit = readRDS(paste0(output_dir,"/",pheno,"_",treat,"_",prolif,"_LMM_noreps.rds"))
    data = readRDS(paste0(output_dir,"/",pheno,"_",treat,"_",prolif,"_data_noreps.rds"))
    to_sort = lapply(my_fit,jtools::summ)
    pvals = unlist(lapply(to_sort, function(x) x$coeftable["genotype:scaled_phenotype", "p"]))
    beta = unlist(lapply(to_sort, function(x) x$coeftable["genotype:scaled_phenotype", "Est."]))
    
    qvals = qvalue::qvalue(pvals, pi0.method = "bootstrap")$qvalues
    qvals = sort(qvals)
    
    data = do.call("rbind",data)
    
    genes = data %>%
      dplyr::mutate(gene_var = paste0(gene_name,"-",variant_id)) %>%
      dplyr::select(gene_var,gene_name) %>%
      distinct() %>%
      column_to_rownames("gene_var")
    all_my_genes[[paste0(pheno,"_",treat)]] = data.frame(pvals = pvals[names(qvals)],
                                                         qvals = qvals,
                                                         gene_var = names(qvals), 
                                                         gene = genes[names(qvals),"gene_name"],
                                                         beta = beta[names(qvals)], 
                                                         phenotype = pheno,
                                                         treatment = treat)
    
  }
}

all_my_genes = do.call("rbind",all_my_genes)
write_tsv(all_my_genes,paste0(output_dir,"/",pheno,"_LMM_noreps_all_results.tsv"))
all_my_genes = read_tsv(paste0(output_dir,"/",pheno,"_LMM_noreps_all_results.tsv"))
# filter to 20% FDR
# all_my_genes = all_my_genes %>%
#   dplyr::filter(qvals < 0.20)
#### compare with Sam's phago CRISPR screen results ######
### check phagocytosis genes
# from Sam
sam_neuroid = read_csv("../../../resources/Sam_NeuroID/NeuroID_L20vsT20_Forward_gdata.csv") %>%
  dplyr::rename(gene = id, CRISPR_score = Score, CRISPR_FDR = FDR)
table(all_my_genes$gene %in% sam_neuroid$id)

all_my_genes = all_my_genes %>%
  dplyr::left_join(sam_neuroid)

##### pathway analysis ########


sign_res_untreated_phago = sign_res %>%
  dplyr::filter(condition == "phagocytosis_untreated_Not_proliferating")
sign_res_untreated_migr = sign_res %>%
  dplyr::filter(condition == "migration_untreated_Not_proliferating")

all_untreated_phago = int$phagocytosis_untreated_Not_proliferating %>%
  dplyr::arrange(pval_adj_bh)
all_untreated_migr = int$migration_untreated_Not_proliferating %>%
  dplyr::arrange(pval_adj_bh)


all_untreated_phago_sorted = 1:length(all_untreated_phago$pval_adj_bh)
names(all_untreated_phago_sorted) = symbol_to_entrez(all_untreated_phago$gene_name)

all_untreated_migr_sorted = 1:length(all_untreated_migr$pval_adj_bh)
names(all_untreated_migr_sorted) = symbol_to_entrez(all_untreated_migr$gene_name)

pathways = msigdbr("human", category="H")
pathways = split(as.character(pathways$entrez_gene), pathways$gs_name)

phago_res = list()
migr_res = list()

phago_res[["gsea_msigdb_hallmark"]] = fgsea::fgsea(pathways=pathways, 
                                                 stats=all_untreated_phago_sorted,
                                                 eps=0) %>%
  tidyr::as_tibble() %>%
  dplyr::arrange(desc(NES))

p1 = phago_res[["gsea_msigdb_hallmark"]] %>%
  dplyr::filter(padj < 0.05) %>%
  ggplot(aes(reorder(pathway, NES), NES)) +
  geom_col() +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score") + 
  theme_minimal() + 
  ggtitle("Untreated phagocytosis eGene interactions: GSEA")

migr_res[["gsea_msigdb_hallmark"]] = fgsea::fgsea(pathways=pathways, 
                                                   stats=all_untreated_migr_sorted,
                                                   eps=0) %>%
  tidyr::as_tibble() %>%
  dplyr::arrange(desc(NES))


p2 = migr_res[["gsea_msigdb_hallmark"]] %>%
  dplyr::filter(padj < 0.05) %>%
  ggplot(aes(reorder(pathway, NES), NES)) +
  geom_col() +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score") + 
  theme_minimal() + 
  ggtitle("Untreated migration eGene interactions: GSEA")

png(paste0(outdir,"/msigdb_hallmark_gsea.png"),
    width = 12, height = 6, res = 400,units = "in", type = "cairo")
p1 + p2 + patchwork::plot_annotation(title = "DE genes: Hallmark pathways NES from GSEA")
dev.off()