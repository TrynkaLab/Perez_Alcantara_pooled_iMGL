# inspect pheno interactions

library(lme4)
library(lmerTest)
library(interactions)
library(tidyverse)
library(variancePartition)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(pheatmap)
library(ggrepel)
library(decoupleR)
library(enrichplot)
library(clusterProfiler)
library(msigdbr)
library(qvalue)
source("./functions.R")

treatment_cols =  c(untreated = "#8D918B", IFN = "#3A5683", LPS = "#F8766D")

# interaction plots with partial residuals
output_dir = "../../data/results/3.run_LMM_pheno_interactions"
dir.create(output_dir, recursive = TRUE)

prolif = "Not_proliferating"
donor_blacklist="letw_5,lizq_3,zaie_1,romx_2,seru_7,qonc_2,iukl_1,curn_3,boqx_2,garx_2,sojd_3,yoch_6"
pcnumber = 15

##### data from 3.run_LMM_pheno_interactions.R
# 
for(pheno in c("phagocytosis","migration")){
for(treat in c("untreated","IFN","LPS") ){
my_fit = readRDS(paste0(output_dir,"/",pheno,"_",treat,"_",prolif,"_LMM_noreps.rds"))
data = readRDS(paste0(output_dir,"/",pheno,"_",treat,"_",prolif,"_data_noreps.rds"))
to_sort = lapply(my_fit,jtools::summ)
pvals = unlist(lapply(to_sort, function(x) x$coeftable["genotype:scaled_phenotype", "p"]))
tstats = unlist(lapply(to_sort, function(x) x$coeftable["genotype:scaled_phenotype", "t val."]))
beta = unlist(lapply(to_sort, function(x) x$coeftable["genotype:scaled_phenotype", "Est."]))
se=  unlist(lapply(to_sort, function(x) x$coeftable["genotype:scaled_phenotype", "S.E."]))
qvals = qvalue::qvalue(pvals, pi0.method = "bootstrap")$qvalues
qvals = sort(qvals)

data = do.call("rbind",data)

genes = data %>%
  dplyr::mutate(gene_var = paste0(gene_name,"-",variant_id)) %>%
  dplyr::select(gene_var,gene_name) %>%
  distinct() %>%
  column_to_rownames("gene_var")
all_my_genes = data.frame(pvals = pvals[names(qvals)],qvals = qvals,
                    gene_var = names(qvals), gene = genes[names(qvals),"gene_name"],
                    beta = beta[names(qvals)],se = se,
                    t= tstats,
                    phenotype = pheno,
                    treatment = treat)

write_csv(all_my_genes,paste0(output_dir,"/",pheno,"_",treat,"_",prolif,"_LMM_noreps_formatted.csv"))

}
}
all_my_genes = list()
for(pheno in c("phagocytosis","migration")){
  for(treat in c("untreated","IFN","LPS") ){
    all_my_genes[[paste0(pheno,"_",treat)]]=read_csv(paste0(output_dir,"/",pheno,"_",treat,"_",prolif,"_LMM_noreps_formatted.csv"))
  }
}
all_my_genes = do.call("rbind",all_my_genes) %>%
  dplyr::mutate(FDR_20_sign  = if_else(qvals < 0.20,"SIGN","NOT SIGN"),
                pheno_treat = paste0(phenotype,"_",treatment))

table(all_my_genes$pheno_treat,all_my_genes$FDR_20_sign)

# save
write_csv(all_my_genes,paste0(output_dir,"/all_LMM_noreps_formatted.csv"))

### extracting interesting candidate genes
jeremy_ad_candidates = read_table("../../../resources/Jeremy_medrXiv_AD_candidate_genes.txt")
jeremy_set1_ad_candidates = read_csv("../../../resources/AD_PD_gene_sets_Andrew/Set1_Jeremy_candidates.csv")
jeremy_set3_pd_candidates = read_csv("../../../resources/AD_PD_gene_sets_Andrew/Set3_manually_curated.csv") %>%
  dplyr::filter(str_detect(disease,"PD")) %>%
  pull(gene_name)
# lead phagocytosis sam
sam_crispr = read_tsv("../../../CRISPR/OTAR2065_phagocytosis_CRISPR/data/2024_07_sam_screen_results.tsv") %>%
  dplyr::arrange(FDR) %>%
  slice_head(n = 30) %>%
  pull(id)

# lead phagocytosis kampmann
kampmann_crispr_phago = read_csv("../../../resources/CRISPRbrain_Kampmanm/kampmann_CRISPRi_iTF_microglia_phagocytosis.csv") %>%
  dplyr::mutate(FDR = qvalue(`P Value`)$qvalues) %>% # nothing seems to be significant
  dplyr::arrange(`P Value`) %>%
  slice_head(n = 30) %>%
  pull(Gene)

# lead immune activation kampmann
kampmann_crispr_immune = read_csv("../../../resources/CRISPRbrain_Kampmanm/kampmann_CRISPRi_iTF_microglia_immune_activation.csv") %>%
  dplyr::mutate(FDR = qvalue(`P Value`)$qvalues) %>% # nothing seems to be significant
  dplyr::arrange(`P Value`) %>%
  slice_head(n = 30) %>%
  pull(Gene)


pathways = list()
pathways[["jeremy_ad_candidates"]] = jeremy_ad_candidates$symbol
pathways[["set1_jeremy_ad_candidates"]] = jeremy_set1_ad_candidates$gene_sym # larger
pathways[["set3_jeremy_pd_candidates"]] = jeremy_set3_pd_candidates
pathways[["lead_phagocytosis_sam"]] = sam_crispr 
pathways[["lead_phagocytosis_kampmann"]] = kampmann_crispr_phago
pathways[["lead_immune_kampmann"]] = kampmann_crispr_immune 

# extract all genes from my data
# sort by p-val of interaction
all_my_genes$entrez = symbol_to_entrez(all_my_genes$gene)
all_my_genes = all_my_genes %>%
  dplyr::mutate(minus_log10pval = -log10(pvals))

res_df = list()
ifn_sorted = list()
untreated_sorted = list()
lps_sorted = list()
for(pheno in c("phagocytosis","migration")){
  
  untreated_sorted[[pheno]] =  all_my_genes %>%
    dplyr::filter(treatment == "untreated" & phenotype == pheno) %>%
    dplyr::arrange( desc(beta)) %>%
     pull(beta)
  
  names(untreated_sorted[[pheno]]) = all_my_genes %>%
    dplyr::filter(treatment == "untreated"  & phenotype == pheno) %>%
    dplyr::arrange( desc(beta)) %>%
    pull(gene)
  untreated_sorted[[pheno]] = untreated_sorted[[pheno]][!is.na(names(untreated_sorted))]

  ifn_sorted[[pheno]] =  all_my_genes %>%
    dplyr::filter(treatment == "IFN"  & phenotype == pheno) %>%
    dplyr::arrange( desc(beta)) %>%
    pull(beta)
  
  names(ifn_sorted[[pheno]]) = all_my_genes %>%
    dplyr::filter(treatment == "IFN"  & phenotype == pheno) %>%
    dplyr::arrange( desc(beta)) %>%
    pull(gene)
  ifn_sorted[[pheno]] = ifn_sorted[[pheno]][!is.na(names(ifn_sorted))]

  lps_sorted[[pheno]] =  all_my_genes %>%
    dplyr::filter(treatment == "LPS"  & phenotype == pheno) %>%
    dplyr::arrange( desc(beta)) %>%
     pull(beta)
  
  names(lps_sorted[[pheno]]) = all_my_genes %>%
    dplyr::filter(treatment == "LPS" & phenotype == pheno) %>%
    dplyr::arrange( desc(beta)) %>%
    pull(gene)
  lps_sorted[[pheno]] = lps_sorted[[pheno]][!is.na(names(lps_sorted))]

  ### GSEA
  

  
  untreated_res = fgsea::fgsea(pathways=pathways, 
                                                    stats=untreated_sorted[[pheno]],
                                                    # scoreType="pos",
                                                    eps=0) %>%
    tidyr::as_tibble() %>%
    dplyr::arrange(desc(NES))
  
  lps_res = fgsea::fgsea(pathways=pathways, 
                                              stats=lps_sorted[[pheno]],
                                              # scoreType="pos",
                                              eps=0) %>%
    tidyr::as_tibble() %>%
    dplyr::arrange(desc(NES))
  
  
  ifn_res = fgsea::fgsea(pathways=pathways, 
                                              stats=ifn_sorted[[pheno]],
                                              # scoreType="pos",
                                              eps=0) %>%
    tidyr::as_tibble() %>%
    dplyr::arrange(desc(NES))
  
  
  ifn_res  = ifn_res %>%
    dplyr::mutate(contrast = paste0("IFN ", pheno))
  lps_res  = lps_res %>%
    dplyr::mutate(contrast =  paste0("LPS ", pheno))
  untreated_res = untreated_res %>%
    dplyr::mutate(contrast =  paste0("untreated ", pheno))
  
  res_df[[pheno]] = rbind(ifn_res , 
                 lps_res,
                 untreated_res  )
}

# tried abs(beta), -log10pval - nothing significant

########### human hallmark pathways #################
# with directionality

pathways = msigdbr("human", category="H")
pathways = split(as.character(pathways$entrez_gene), pathways$gs_name)

res_df = list()
ifn_sorted = list()
untreated_sorted = list()
lps_sorted = list()
for(pheno in c("phagocytosis","migration")){
  
  untreated_sorted[[pheno]] =  all_my_genes %>%
    dplyr::filter(treatment == "untreated" & phenotype == pheno) %>%
    dplyr::arrange( desc( minus_log10pval)) %>%
    pull( minus_log10pval)
  
  names(untreated_sorted[[pheno]]) = all_my_genes %>%
    dplyr::filter(treatment == "untreated"  & phenotype == pheno) %>%
    dplyr::arrange( desc( minus_log10pval)) %>%
    pull(entrez)
  untreated_sorted[[pheno]] = untreated_sorted[[pheno]][!is.na(names(untreated_sorted))]
  
  ifn_sorted[[pheno]] =  all_my_genes %>%
    dplyr::filter(treatment == "IFN"  & phenotype == pheno) %>%
    dplyr::arrange( desc( minus_log10pval)) %>%
    pull( minus_log10pval)
  
  names(ifn_sorted[[pheno]]) = all_my_genes %>%
    dplyr::filter(treatment == "IFN"  & phenotype == pheno) %>%
    dplyr::arrange( desc( minus_log10pval)) %>%
    pull(entrez)
  ifn_sorted[[pheno]] = ifn_sorted[[pheno]][!is.na(names(ifn_sorted))]
  
  lps_sorted[[pheno]] =  all_my_genes %>%
    dplyr::filter(treatment == "LPS"  & phenotype == pheno) %>%
    dplyr::arrange( desc( minus_log10pval)) %>%
    pull( minus_log10pval)
  
  names(lps_sorted[[pheno]]) = all_my_genes %>%
    dplyr::filter(treatment == "LPS" & phenotype == pheno) %>%
    dplyr::arrange( desc( minus_log10pval)) %>%
    pull(entrez)
  lps_sorted[[pheno]] = lps_sorted[[pheno]][!is.na(names(lps_sorted))]
  
  ### GSEA
  
  untreated_res = fgsea::fgsea(pathways=pathways, 
                               stats=untreated_sorted[[pheno]],
                               scoreType="pos",
                               eps=0) %>%
    tidyr::as_tibble() %>%
    dplyr::arrange(desc(NES)) %>%
    dplyr::mutate(pathway = str_remove(pathway,"HALLMARK_")) %>%
    dplyr::mutate(pathway = str_replace_all(pathway,"_", " "))  %>%
    dplyr::mutate(pathway = str_to_sentence(pathway))
  
  lps_res = fgsea::fgsea(pathways=pathways, 
                         stats=lps_sorted[[pheno]],
                         scoreType="pos",
                         eps=0) %>%
    tidyr::as_tibble() %>%
    dplyr::arrange(desc(NES)) %>%
    dplyr::mutate(pathway = str_remove(pathway,"HALLMARK_")) %>%
    dplyr::mutate(pathway = str_replace_all(pathway,"_", " "))  %>%
    dplyr::mutate(pathway = str_to_sentence(pathway))
  
  
  ifn_res = fgsea::fgsea(pathways=pathways, 
                         stats=ifn_sorted[[pheno]],
                         scoreType="pos",
                         eps=0)  %>%
    tidyr::as_tibble() %>%
    dplyr::arrange(desc(NES)) %>%
    dplyr::mutate(pathway = str_remove(pathway,"HALLMARK_")) %>%
    dplyr::mutate(pathway = str_replace_all(pathway,"_", " "))  %>%
    dplyr::mutate(pathway = str_to_sentence(pathway))
  
  
  ifn_res  = ifn_res %>%
    dplyr::mutate(contrast = paste0("IFN ", pheno))
  lps_res  = lps_res %>%
    dplyr::mutate(contrast =  paste0("LPS ", pheno))
  untreated_res = untreated_res %>%
    dplyr::mutate(contrast =  paste0("untreated ", pheno))

  # nothing significant
  p1 = untreated_res%>%
    dplyr::filter(padj < 0.05) %>%
    ggplot(aes(reorder(pathway, NES), NES)) +
    geom_col() +
    coord_flip() +
    labs(x="Pathway", y="") + 
    theme_minimal() + 
    ggtitle(paste0("Untreated ", pheno))

p2 = ifn_res%>%
  dplyr::filter(padj < 0.05) %>%
  ggplot(aes(reorder(pathway, NES), NES)) +
  geom_col() +
  coord_flip() +
  labs(x="", y="Normalized Enrichment Score") + 
  theme_minimal() + 
  ggtitle(paste0("IFN ", pheno))

p3 = lps_res %>%
  dplyr::filter(padj < 0.05) %>%
  ggplot(aes(reorder(pathway, NES), NES)) +
  geom_col() +
  coord_flip() +
  labs(x="", y="Normalized Enrichment Score") + 
  theme_minimal() + 
  ggtitle(paste0("LPS ", pheno))



png(paste0(outdir,"/msigdb_hallmark_gsea_", pheno, ".png"),
    width = 12, height = 5, res = 400,units = "in", type = "cairo")
(p1 + p2 + p3 ) + patchwork::plot_annotation(title = "DE genes: Hallmark pathways NES from GSEA")
dev.off()


# which genes in these pathways
ifn_res  = ifn_res  %>% 
  dplyr::filter(padj < 0.05) %>%
  dplyr::select("pathway", "leadingEdge","padj","NES","size") %>% 
  tidyr::unnest(cols = c(leadingEdge)) %>% 
  dplyr::rename(entrez = leadingEdge) %>%
  dplyr::inner_join(all_my_genes[c("gene","entrez")], by="entrez") %>%
  dplyr::group_by(pathway) %>%
  dplyr::reframe(padj = padj, NES = NES,size=size,
                 leading_edge = paste( gene, collapse = ", ")) %>%
  distinct() %>%
  dplyr::arrange(padj) %>%
  dplyr::mutate(treatment = "IFNvsUntreated")

lps_res  = lps_res  %>% 
  dplyr::filter(padj < 0.05) %>%
  dplyr::select("pathway", "leadingEdge","padj","NES","size") %>% 
  tidyr::unnest(cols = c(leadingEdge)) %>% 
  dplyr::rename(entrez = leadingEdge) %>%
  dplyr::inner_join(all_my_genes[c("gene","entrez")], by="entrez") %>%
  dplyr::group_by(pathway) %>%
  dplyr::reframe(padj = padj, NES = NES,size=size,
                 leading_edge = paste( gene, collapse = ", ")) %>%
  distinct() %>%
  dplyr::arrange(padj)  %>%
  dplyr::mutate(treatment = "LPSvsUntreated")

ifn_lps_res  = ifn_lps_res  %>% 
  dplyr::filter(padj < 0.05) %>%
  dplyr::select("pathway", "leadingEdge","padj","NES","size") %>% 
  tidyr::unnest(cols = c(leadingEdge)) %>% 
  dplyr::rename(entrez = leadingEdge) %>%
  dplyr::inner_join(all_my_genes[c("gene","entrez")], by="entrez") %>%
  dplyr::group_by(pathway) %>%
  dplyr::reframe(padj = padj, NES = NES,size=size,
                 leading_edge = paste( gene, collapse = ", ")) %>%
  distinct() %>%
  dplyr::arrange(padj) %>%
  dplyr::mutate(treatment = "IFNvsLPS")

# save results

res_df[[pheno]] = rbind(ifn_res  ,
      lps_res  ,
      ifn_lps_res )


}

write_csv(res_df[[pheno]] ,file=paste0(outdir,"/msigdb_halmmark_gsea_genelists.csv"))
# tried with p-value and signed beta, nothing significant
