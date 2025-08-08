# Inspect limma treatments DEG

library(tidyverse)
library(patchwork)
library(fgsea)
source("./helpers.R")
library(AnnotationDbi)
library(org.Hs.eg.db)
library(pheatmap)
library(ggrepel)
library(decoupleR)
library(enrichplot)
library(clusterProfiler)
library(msigdbr)
library(qvalue)
options(future.globals.maxSize = (120000 * 1024 ^ 2)) #(~120 Gb)


treatment_cols =  c(untreated = "#8D918B", IFN = "#3A5683", LPS = "#F8766D")


outdir = "../../data/results/5.2.Inspect_DEG_limma_treatments"
dir.create(outdir,recursive = TRUE)
res_additive = readRDS("../../data/results/5.1.Diff_expr_limma/additive_res.rds")
for(cont in names(res_additive)){
  res_additive[[cont]] =   res_additive[[cont]]%>%
    dplyr::mutate(diffexp = case_when(logFC > 1 & adj.P.Val < 0.05 ~ "UP",
                                      logFC < -1 & adj.P.Val < 0.05 ~ "DOWN",
                                      .default = "NO"),
                  contrast = cont) 
}


deg = do.call("rbind",res_additive)

# total DEG
table(deg$contrast,deg$diffexp)
# grand total (TRUE == sign DEG)
table(deg$contrast,deg$diffexp!="NO")


# save tables
write_csv(res_additive$IFNvsUntreated,"../../data/results/5.1.Diff_expr_limma/IFNvsUntreated_poslogFC_higher_IFN.csv")
write_csv(res_additive$LPSvsUntreated,"../../data/results/5.1.Diff_expr_limma/LPSvsUntreated_poslogFC_higher_LPS.csv")
write_csv(res_additive$IFNvsLPS,"../../data/results/5.1.Diff_expr_limma/IFNvsLPS_poslogFC_higher_IFN.csv")

### volcano plots
p1 = res_additive$IFNvsUntreated  %>%
  dplyr::filter(diffexp == "NO") %>%   # Filter the rows that meet the conditions
  dplyr::sample_frac(0.05) %>%         # Sample 10% of the filtered data
  dplyr::bind_rows(res_additive$IFNvsUntreated  %>% filter(!(diffexp == "NO"))) %>%  # Combine with the rest of the data
  dplyr::mutate(pval=-log10(adj.P.Val)) %>% 
  
  dplyr::mutate(pval=case_when(pval=="Inf"~350,
                               .default = pval)) %>% # substitute for very extreme pvalue, at the limit of the precision
  
  ggplot(aes(
    x = logFC,
    y = pval,
    col = diffexp,
    label = symbol
  )) +
  geom_point() +
  theme_minimal() +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),breaks = c(0,10,50,100,200,400)) +
  ggrepel::geom_text_repel(color = "black") +
  scale_color_manual(values = c("DOWN"="#3A5683",
                                "NO"= "#A7A9AE",
                                "UP"="#74121D")) +
  geom_vline(xintercept = c(-1, 1),linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 2,linetype = "dashed", color = "grey") +
  ylab("-log10(adj. p-val)") + 
  labs(fill = "Diff. expr.") +
  ggtitle("IFN vs untreated")


p2 = res_additive$LPSvsUntreated %>%
  dplyr::filter(diffexp == "NO") %>%   # Filter the rows that meet the conditions
  dplyr::sample_frac(0.05) %>%         # Sample 10% of the filtered data
  dplyr::bind_rows( res_additive$LPSvsUntreated %>% filter(!(diffexp == "NO"))) %>%  # Combine with the rest of the data
  dplyr::mutate(pval=-log10(adj.P.Val)) %>% 
  
  dplyr::mutate(pval=case_when(pval=="Inf"~350,
                               .default = pval)) %>% # substitute for very extreme pvalue, at the limit of the precision
  
  ggplot(aes(
    x = logFC,
    y = pval,
    col = diffexp,
    label = symbol
  )) +
  geom_point() +
  theme_minimal() +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),breaks = c(0,10,50,100,200,400)) +
  ggrepel::geom_text_repel(color = "black") +
  scale_color_manual(values = c("DOWN"="#3A5683",
                                "NO"= "#A7A9AE",
                                "UP"="#74121D")) +
  
  geom_vline(xintercept = c(-1, 1),linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 2,linetype = "dashed", color = "grey") +
  ylab("-log10(adj. p-val)") + 
  labs(fill = "Diff. expr.") +
  ggtitle("LPS vs untreated")

p3 = res_additive$IFNvsLPS %>%
  dplyr::filter(diffexp == "NO") %>%   # Filter the rows that meet the conditions
  dplyr::sample_frac(0.05) %>%         # Sample 10% of the filtered data
  dplyr::bind_rows(res_additive$IFNvsLPS %>% filter(!(diffexp == "NO"))) %>%  # Combine with the rest of the data
  dplyr::mutate(pval=-log10(adj.P.Val)) %>% 
  
  dplyr::mutate(pval=case_when(pval=="Inf"~350,
                               .default = pval)) %>% # substitute for very extreme pvalue, at the limit of the precision
  
  ggplot( aes(x=logFC, 
              y=pval, 
              col=diffexp, label=symbol)) +
  geom_point() + 
  theme_minimal() +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),breaks = c(0,10,50,100,200,400)) +
  ggrepel::geom_text_repel(color = "black") +
  scale_color_manual(values = c("DOWN"="#3A5683",
                                "NO"= "#A7A9AE",
                                "UP"="#74121D")) +
  
  geom_vline(xintercept = c(-1, 1),linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 2,linetype = "dashed", color = "grey") +
  ylab("-log10(adj. p-val)") + 
  labs(fill = "Diff. expr.") +
  ggtitle("IFN vs LPS")

png(paste0(outdir,"/volcanos_res_additive.png"),
    width = 15, height = 4, res = 400,units = "in", type = "cairo")
patchwork::wrap_plots(p1,p2,p3, guides = "collect")
dev.off()

## jeremy genes
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

# signal significant differentially expressed genes

res_additive$IFNvsUntreated = res_additive$IFNvsUntreated %>%
  dplyr::arrange(desc(abs(t)))

res_additive$LPSvsUntreated = res_additive$LPSvsUntreated %>%
  dplyr::arrange(desc(abs(t)))

res_additive$IFNvsLPS = res_additive$IFNvsLPS %>%
  dplyr::arrange(desc(abs(t)))


ifn_sorted = abs(res_additive$IFNvsUntreated$t)
names(ifn_sorted) = res_additive$IFNvsUntreated$symbol
ifn_sorted = sort(ifn_sorted,decreasing = TRUE)
ifn_sorted = ifn_sorted[!is.na(names(ifn_sorted))]
res_additive$IFNvsUntreated$entrez = symbol_to_entrez(res_additive$IFNvsUntreated$symbol)

lps_sorted = abs(res_additive$LPSvsUntreated$t)
names(lps_sorted) = res_additive$LPSvsUntreated$symbol
lps_sorted = sort(lps_sorted,decreasing = TRUE)
lps_sorted = lps_sorted[!is.na(names(lps_sorted))]
res_additive$LPSvsUntreated$entrez = symbol_to_entrez(res_additive$LPSvsUntreated$symbol)


IFN_LPS_sorted = abs(res_additive$IFNvsLPS$t)
names(IFN_LPS_sorted) = res_additive$IFNvsLPS$symbol
IFN_LPS_sorted = sort(IFN_LPS_sorted,decreasing = TRUE)
IFN_LPS_sorted = IFN_LPS_sorted[!is.na(names(IFN_LPS_sorted))]
res_additive$IFNvsLPS$entrez = symbol_to_entrez(res_additive$IFNvsLPS$symbol)

pathways = list()
pathways[["jeremy_ad_candidates"]] = jeremy_ad_candidates$symbol
pathways[["set1_jeremy_ad_candidates"]] = jeremy_set1_ad_candidates$gene_sym # larger
pathways[["set3_jeremy_pd_candidates"]] = jeremy_set3_pd_candidates
pathways[["lead_phagocytosis_sam"]] = sam_crispr 
pathways[["lead_phagocytosis_kampmann"]] = kampmann_crispr_phago
pathways[["lead_immune_kampmann"]] = kampmann_crispr_immune 

### GSEA

ifn_res = list()
lps_res = list()
ifn_lps_res = list()

ifn_res[["custom_genesets"]] = fgsea::fgsea(pathways=pathways, 
                                                      stats=ifn_sorted,
                                                      scoreType="pos",
                                                      eps=0) %>%
  tidyr::as_tibble() %>%
  dplyr::arrange(desc(NES))

lps_res[["custom_genesets"]] = fgsea::fgsea(pathways=pathways, 
                                                      stats=lps_sorted,
                                                      scoreType="pos",
                                                      eps=0) %>%
  tidyr::as_tibble() %>%
  dplyr::arrange(desc(NES))

ifn_lps_res[["custom_genesets"]] = fgsea::fgsea(pathways=pathways, 
                                                      stats=IFN_LPS_sorted,
                                                      scoreType="pos",
                                                      eps=0) %>%
  tidyr::as_tibble() %>%
  dplyr::arrange(desc(NES))

ifn_res[["custom_genesets"]]  = ifn_res[["custom_genesets"]] %>%
  dplyr::mutate(contrast = "IFN vs untreated")
lps_res[["custom_genesets"]]  = lps_res[["custom_genesets"]] %>%
  dplyr::mutate(contrast = "LPS vs untreated")
lps_res$custom_genesets$leadingEdge
ifn_lps_res[["custom_genesets"]] = ifn_lps_res[["custom_genesets"]] %>%
  dplyr::mutate(contrast = "IFN vs LPS")
plotEnrichment(pathways$set1_jeremy_ad_candidates, lps_sorted)

res_df = rbind(ifn_res[["custom_genesets"]] , lps_res[["custom_genesets"]],ifn_lps_res[["custom_genesets"]]  )


p = res_df %>%
  dplyr::filter(pathway != "jeremy_ad_candidates") %>%
  dplyr::mutate(contrast = factor(contrast,levels = c("IFN vs untreated","LPS vs untreated","IFN vs LPS")),
                pathway = case_when(pathway == "set3_jeremy_pd_candidates" ~ "PD candidates",
                                    pathway == "set1_jeremy_ad_candidates" ~ "AD candidates",
                                    pathway == "lead_phagocytosis_sam" ~ "Lead phagocytosis CRISPR",
                                    pathway == "lead_phagocytosis_kampmann" ~ "Lead phagocytosis CRISPR (Kampmann)",
                                    pathway == "lead_immune_kampmann" ~ "Lead immune CRISPR (Kampmann)"),
                minus_log10_pval = -log10(padj),
                p_val_sig_plot =  c("***", "**", "*", "")[findInterval(padj, c(0.001, 0.01, 0.05)) + 1]) %>%
  dplyr::mutate(pathway = factor(pathway,levels = c( "AD candidates", "PD candidates",
                                                     "Lead phagocytosis CRISPR", "Lead phagocytosis CRISPR (Kampmann)","Lead immune CRISPR (Kampmann)"))) %>%
ggplot( aes(y=pathway, x=contrast, fill = minus_log10_pval)) + 
  geom_raster() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), # lables vertical
        strip.text.y = element_blank()) +  #remove facet bar on y 
  scale_fill_gradient(low = "white", high = "#74121D",name = bquote(-log10(P^{GSEA}))) +
  geom_text(aes(label = p_val_sig_plot), col = "white", size = 5) + 
  ggtitle("GSEA enrichment in candidate gene sets") +
  theme_minimal() + 
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12)) +
  # facet_grid(cols = vars(`gene list` )) +
  theme(legend.position="right")


png(paste0(outdir,"/GSEA_candidate_genes_tileplot.png"), 
    height = 4.5,width =10,units = "in", res = 500)
plot(p)
dev.off()

# with directionality

ifn_sorted = res_additive$IFNvsUntreated %>%
  dplyr::arrange(desc(t)) %>%
  pull(t)
names(ifn_sorted) = res_additive$IFNvsUntreated %>%
  dplyr::arrange(desc(t)) %>%
  pull(symbol)
lps_sorted = res_additive$LPSvsUntreated %>%
  dplyr::arrange(desc(t)) %>%
  pull(t)
names(lps_sorted) = res_additive$LPSvsUntreated %>%
  dplyr::arrange(desc(t)) %>%
  pull(symbol)
IFN_LPS_sorted = res_additive$IFNvsLPS %>%
  dplyr::arrange(desc(t)) %>%
  pull(t)
names(IFN_LPS_sorted) = res_additive$IFNvsLPS %>%
  dplyr::arrange(desc(t)) %>%
  pull(symbol)


ifn_res_signed = list()
lps_res_signed = list()
ifn_lps_res_signed = list()

ifn_res_signed[["custom_genesets"]] = fgsea::fgsea(pathways=pathways, 
                                            stats=ifn_sorted,
                                            eps=0) %>%
  tidyr::as_tibble() %>%
  dplyr::arrange(desc(NES))

lps_res_signed[["custom_genesets"]] = fgsea::fgsea(pathways=pathways, 
                                            stats=lps_sorted,
                                            eps=0) %>%
  tidyr::as_tibble() %>%
  dplyr::arrange(desc(NES))

ifn_lps_res_signed[["custom_genesets"]] = fgsea::fgsea(pathways=pathways, 
                                                stats=IFN_LPS_sorted,
                                                eps=0) %>%
  tidyr::as_tibble() %>%
  dplyr::arrange(desc(NES))

ifn_res_signed[["custom_genesets"]]
lps_res_signed[["custom_genesets"]]
lps_res_signed$custom_genesets$leadingEdge
ifn_lps_res_signed[["custom_genesets"]]
ifn_lps_res_signed$custom_genesets$leadingEdge

plotEnrichment(pathways$set1_jeremy_ad_candidates, lps_sorted)

# they accumulate on both sides

# AD candidate genes are enriched in differentially expressed genes


### GO enrichment - with upregulation / downregulation directionality
ifn_sorted = res_additive$IFNvsUntreated %>%
  dplyr::arrange(desc(t)) %>%
  pull(t)
names(ifn_sorted) = res_additive$IFNvsUntreated %>%
  dplyr::arrange(desc(t)) %>%
  pull(entrez)
lps_sorted = res_additive$LPSvsUntreated %>%
  dplyr::arrange(desc(t)) %>%
  pull(t)
names(lps_sorted) = res_additive$LPSvsUntreated %>%
  dplyr::arrange(desc(t)) %>%
  pull(entrez)
IFN_LPS_sorted = res_additive$IFNvsLPS %>%
  dplyr::arrange(desc(t)) %>%
  pull(t)
names(IFN_LPS_sorted) = res_additive$IFNvsLPS %>%
  dplyr::arrange(desc(t)) %>%
  pull(entrez)

ifn_res[["GO_BP"]] <- gseGO(geneList     = ifn_sorted,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              minGSSize    = 50,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE,
              eps = 0)


lps_res[["GO_BP"]] <- gseGO(geneList     = lps_sorted,
                            OrgDb        = org.Hs.eg.db,
                            ont          = "BP",
                            minGSSize    = 50,
                            maxGSSize    = 500,
                            pvalueCutoff = 0.05,
                            verbose      = FALSE,
                            eps = 0)

ifn_lps_res[["GO_BP"]] <- gseGO(geneList     = IFN_LPS_sorted,
                            OrgDb        = org.Hs.eg.db,
                            ont          = "BP",
                            minGSSize    = 50,
                            maxGSSize    = 500,
                            pvalueCutoff = 0.05,
                            verbose      = FALSE,
                            eps = 0)
# filter
ifn_res[["GO_BP"]] = filter(ifn_res[["GO_BP"]], p.adjust < .05, qvalue < 0.2)
lps_res[["GO_BP"]] = filter(lps_res[["GO_BP"]], p.adjust < .05, qvalue < 0.2)
ifn_lps_res[["GO_BP"]] = filter(ifn_lps_res[["GO_BP"]], p.adjust < .05, qvalue < 0.2)

p1 = dotplot(ifn_res[["GO_BP"]], showCategory=10) + ggtitle("DEGs IFN vs untreated")
p2 = dotplot(lps_res[["GO_BP"]], showCategory=10) + ggtitle("DEGs LPS vs untreated")
p3 = dotplot(ifn_lps_res[["GO_BP"]], showCategory=10) + ggtitle("DEGs IFN vs LPS")


png(paste0(outdir,"/GO_BP_gsea.png"),
    width = 15, height = 7, res = 400,units = "in", type = "cairo")
(p1 + p2 + p3) + patchwork::plot_annotation(title = "GO BP")

dev.off()
# human hallmark pathways #########


pathways = msigdbr("human", category="H")
pathways = split(as.character(pathways$entrez_gene), pathways$gs_name)

ifn_res[["gsea_msigdb_hallmark"]] = fgsea::fgsea(pathways=pathways, 
                                                 stats=ifn_sorted,
                                                 eps=0) %>%
  tidyr::as_tibble() %>%
  dplyr::arrange(desc(NES)) %>%
  dplyr::mutate(pathway = str_remove(pathway,"HALLMARK_")) %>%
  dplyr::mutate(pathway = str_replace_all(pathway,"_", " "))  %>%
  dplyr::mutate(pathway = str_to_sentence(pathway))
  

lps_res[["gsea_msigdb_hallmark"]] = fgsea::fgsea(pathways=pathways, 
                                                 stats=lps_sorted,
                                                 eps=0) %>%
  tidyr::as_tibble() %>%
  dplyr::arrange(desc(NES)) %>%
  dplyr::mutate(pathway = str_remove(pathway,"HALLMARK_")) %>%
  dplyr::mutate(pathway = str_replace(pathway,"_", " "))  %>%
  dplyr::mutate(pathway = str_to_sentence(pathway))

ifn_lps_res[["gsea_msigdb_hallmark"]] = fgsea::fgsea(pathways=pathways, 
                                                 stats=IFN_LPS_sorted,
                                                 eps=0) %>%
  tidyr::as_tibble() %>%
  dplyr::arrange(desc(NES)) %>%
  dplyr::mutate(pathway = str_remove(pathway,"HALLMARK_")) %>%
  dplyr::mutate(pathway = str_replace(pathway,"_", " "))  %>%
  dplyr::mutate(pathway = str_to_sentence(pathway))

p1 = ifn_res[["gsea_msigdb_hallmark"]] %>%
  dplyr::filter(padj < 0.05) %>%
  ggplot(aes(reorder(pathway, NES), NES)) +
  geom_col() +
  coord_flip() +
  labs(x="Pathway", y="") + 
  theme_minimal() + 
  ggtitle("IFN vs untreated")

p2 = lps_res[["gsea_msigdb_hallmark"]] %>%
  dplyr::filter(padj < 0.05) %>%
  ggplot(aes(reorder(pathway, NES), NES)) +
  geom_col() +
  coord_flip() +
  labs(x="", y="Normalized Enrichment Score") + 
  theme_minimal() + 
  ggtitle("LPS vs untreated")

p3 = ifn_lps_res[["gsea_msigdb_hallmark"]] %>%
  dplyr::filter(padj < 0.05) %>%
  ggplot(aes(reorder(pathway, NES), NES)) +
  geom_col() +
  coord_flip() +
  labs(x="", y="") + 
  theme_minimal() + 
  ggtitle("IFN vs LPS")

png(paste0(outdir,"/msigdb_hallmark_gsea.png"),
    width = 12, height = 5, res = 400,units = "in", type = "cairo")
(p1 + p2 + p3 ) + patchwork::plot_annotation(title = "DE genes: Hallmark pathways NES from GSEA")
dev.off()


# which genes in these pathways
ifn_res[["gsea_msigdb_hallmark"]] = ifn_res[["gsea_msigdb_hallmark"]] %>% 
  dplyr::filter(padj < 0.05) %>%
  dplyr::select("pathway", "leadingEdge","padj","NES","size") %>% 
  tidyr::unnest(cols = c(leadingEdge)) %>% 
  dplyr::rename(entrez = leadingEdge) %>%
  dplyr::inner_join(res_additive$IFNvsUntreated[c("symbol","entrez")], by="entrez") %>%
  dplyr::group_by(pathway) %>%
  dplyr::reframe(padj = padj, NES = NES,size=size,
                 leading_edge = paste(symbol, collapse = ", ")) %>%
  distinct() %>%
  dplyr::arrange(padj) %>%
  dplyr::mutate(contrast = "IFNvsUntreated")

lps_res[["gsea_msigdb_hallmark"]] = lps_res[["gsea_msigdb_hallmark"]] %>% 
  dplyr::filter(padj < 0.05) %>%
  dplyr::select("pathway", "leadingEdge","padj","NES","size") %>% 
  tidyr::unnest(cols = c(leadingEdge)) %>% 
  dplyr::rename(entrez = leadingEdge) %>%
  dplyr::inner_join(res_additive$LPSvsUntreated[c("symbol","entrez")], by="entrez") %>%
  dplyr::group_by(pathway) %>%
  dplyr::reframe(padj = padj, NES = NES,size=size,
                 leading_edge = paste(symbol, collapse = ", ")) %>%
  distinct() %>%
  dplyr::arrange(padj)  %>%
  dplyr::mutate(contrast = "LPSvsUntreated")

ifn_lps_res[["gsea_msigdb_hallmark"]] = ifn_lps_res[["gsea_msigdb_hallmark"]] %>% 
  dplyr::filter(padj < 0.05) %>%
  dplyr::select("pathway", "leadingEdge","padj","NES","size") %>% 
  tidyr::unnest(cols = c(leadingEdge)) %>% 
  dplyr::rename(entrez = leadingEdge) %>%
  dplyr::inner_join(res_additive$IFNvsLPS[c("symbol","entrez")], by="entrez") %>%
  dplyr::group_by(pathway) %>%
  dplyr::reframe(padj = padj, NES = NES,size=size,
                 leading_edge = paste(symbol, collapse = ", ")) %>%
  distinct() %>%
  dplyr::arrange(padj) %>%
  dplyr::mutate(contrast = "IFNvsLPS")

# save results

rbind(ifn_res[["gsea_msigdb_hallmark"]] ,
      lps_res[["gsea_msigdb_hallmark"]] ,
      ifn_lps_res[["gsea_msigdb_hallmark"]]) %>%
  write_csv(.,file=paste0(outdir,"/msigdb_halmmark_gsea_genelists.csv"))


######### TF activity inference ######
net = get_collectri(organism='human', split_complexes=FALSE)
net

# check TF function in microglia
# https://www.jci.org/articles/view/90604

# from t-statistic
# Run ulm


deg_t_filtered = deg %>%
  dplyr::select(symbol, t, contrast) %>%
  dplyr::arrange(symbol) %>%
  tidyr::pivot_wider(id_cols = c(symbol), 
                     names_from = contrast, values_from = t, values_fill =0) %>%
  # remove rows with NA
  tidyr::drop_na() %>%
  column_to_rownames(var = "symbol") %>% 
  as.matrix() 

TF_scores_t = run_ulm(mat=deg_t_filtered, net=net, .source='source', 
                      .target='target',
                               .mor='mor', minsize = 5)

table(TF_scores_t$p_value < 0.05, TF_scores_t$condition) # 393 in total

n_tfs = 30
top_sign_tfs = TF_scores_t %>%
  dplyr::filter(p_value < 0.05) %>%
  distinct(source) %>%
  head(n_tfs) %>%
  pull(source)  

toplot = TF_scores_t %>%
  dplyr::filter(source %in% top_sign_tfs) %>%
  tidyr::pivot_wider(id_cols = 'condition', names_from = 'source',
                     values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix() %>%
  t()

palette_length = 100
my_color = colorRampPalette(c("#3A5683", "white","#74121D"))(palette_length)

my_breaks = c(seq(min(toplot), 0, length.out=ceiling(palette_length/2) + 1),
              seq(0.05, max(toplot), length.out=floor(palette_length/2)))

p_t_unscaled = ggplotify::as.ggplot(pheatmap(toplot, border_color = NA, 
                                             color=my_color, breaks = my_breaks,
                                             cluster_cols = FALSE)) +
  ggtitle("TF target scores - unscaled")

pdf(paste0(outdir,"/TF_activity_plots_full_DEG_t_statistic_heatmap.pdf"), 
    height = 8, width = 3)
p_t_unscaled

dev.off()



###### pathway inference with Progeny ########
net = get_progeny(organism = 'human', top = 500)
net

# Run mlm
progeny_scores_t = run_mlm(mat=deg_t_filtered, 
                      net=net, .source='source', .target='target',
                       .mor='weight', minsize = 5)



table(progeny_scores_t$p_value < 0.05, progeny_scores_t$condition) # 19 in total

n_pathways = 30
top_sign_tfs = progeny_scores_t %>%
  dplyr::filter(p_value < 0.05) %>%
  distinct(source) %>%
  head(n_pathways) %>%
  pull(source)  

toplot = progeny_scores_t %>%
  dplyr::filter(source %in% top_sign_tfs) %>%
  tidyr::pivot_wider(id_cols = 'condition', names_from = 'source',
                     values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix() %>%
  t()

palette_length = 100
my_color = colorRampPalette(c("#3A5683", "white","#74121D"))(palette_length)

my_breaks = c(seq(min(toplot), 0, length.out=ceiling(palette_length/2) + 1),
              seq(0.05, max(toplot), length.out=floor(palette_length/2)))

progeny_unscaled = ggplotify::as.ggplot(pheatmap(toplot, border_color = NA, 
                                             color=my_color, breaks = my_breaks,
                                             cluster_cols = FALSE)) +
  ggtitle("Progeny pathway scores - unscaled") +

pdf(paste0(outdir,"/progeny_activity_plots_full_DEG_t_statistic_heatmap.pdf"), 
    height = 5, width = 3.5)
progeny_unscaled

dev.off()

# positive: pathway is upregulated
# negative: pathway is downregulated
