# GO/pathway analyses with GSEA

library(tidyverse)
library(patchwork)
library(fgsea)
library(msigdbr)
source("./helpers.R")
library(AnnotationDbi)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)

outdir = "../../data/results/5.2.GSEA_GO_pathway"
dir.create(outdir,recursive = TRUE)
ifn_vs_untreated = read_table("../../data/results/5.1.Diff_expr_DESeq2/fit_pool_all_pools/DiffExpr_IFN_all_genes_negative_lower_in_IFN.txt")
lps_vs_untreated = read_table("../../data/results/5.1.Diff_expr_DESeq2/fit_pool_all_pools/DiffExpr_LPS_all_genes_negative_lower_in_LPS.txt")

ifn_sorted = ifn_vs_untreated$log2FoldChange
names(ifn_sorted) = symbol_to_entrez(ifn_vs_untreated$symbol)
ifn_sorted = sort(ifn_sorted,decreasing = TRUE)
ifn_sorted = ifn_sorted[!is.na(names(ifn_sorted))]
ifn_vs_untreated$entrez = symbol_to_entrez(ifn_vs_untreated$symbol)

lps_sorted = lps_vs_untreated$log2FoldChange
names(lps_sorted) = symbol_to_entrez(lps_vs_untreated$symbol)
lps_sorted = sort(lps_sorted,decreasing = TRUE)
lps_sorted = lps_sorted[!is.na(names(lps_sorted))]
lps_vs_untreated$entrez = symbol_to_entrez(lps_vs_untreated$symbol)


# human hallmark pathways
pathways = msigdbr("human", category="H")
pathways = split(as.character(pathways$entrez_gene), pathways$gs_name)

ifn_res = list()
lps_res = list()

 ifn_res[["gsea_msigdb_hallmark"]] = fgsea::fgsea(pathways=pathways, 
                                                  stats=ifn_sorted,
                                                  eps=0) %>%
   tidyr::as_tibble() %>%
   dplyr::arrange(desc(NES))
 
 p1 = ifn_res[["gsea_msigdb_hallmark"]] %>%
   dplyr::filter(padj < 0.05) %>%
   ggplot(aes(reorder(pathway, NES), NES)) +
   geom_col() +
   coord_flip() +
   labs(x="Pathway", y="Normalized Enrichment Score") + 
   theme_minimal() + 
   ggtitle("IFN vs untreated")

 
 lps_res[["gsea_msigdb_hallmark"]] = fgsea::fgsea(pathways=pathways, 
                                                  stats=lps_sorted,
                                                  eps=0) %>%
   tidyr::as_tibble() %>%
   dplyr::arrange(desc(NES))
 
 p2 = lps_res[["gsea_msigdb_hallmark"]] %>%
   dplyr::filter(padj < 0.05) %>%
   ggplot(aes(reorder(pathway, NES), NES)) +
   geom_col() +
   coord_flip() +
   labs(x="Pathway", y="Normalized Enrichment Score") + 
   theme_minimal() + 
   ggtitle("LPS vs untreated")
 
 png(paste0(outdir,"/msigdb_hallmark_gsea.png"),
     width = 12, height = 6, res = 400,units = "in", type = "cairo")
 p1 + p2 + patchwork::plot_annotation(title = "DE genes: Hallmark pathways NES from GSEA")
 dev.off()
 
 
 # which genes in these pathways
 ifn_res[["gsea_msigdb_hallmark"]] = ifn_res[["gsea_msigdb_hallmark"]] %>% 
   dplyr::filter(padj < 0.05) %>%
   dplyr::select("pathway", "leadingEdge","padj","NES","size") %>% 
   tidyr::unnest(cols = c(leadingEdge)) %>% 
   dplyr::rename(entrez = leadingEdge) %>%
   dplyr::inner_join(ifn_vs_untreated[c("symbol","entrez")], by="entrez") %>%
   dplyr::group_by(pathway) %>%
   dplyr::reframe(padj = padj, NES = NES,size=size,
                  leading_edge = paste(symbol, collapse = ", ")) %>%
   distinct() %>%
   dplyr::arrange(padj)
 
 # other enrichments
 
 ifn_res[["gsea_GO"]] = clusterProfiler::gseGO(geneList     = ifn_sorted,
               OrgDb        = org.Hs.eg.db,
               ont          = "BP",
               minGSSize    = 100,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               verbose      = TRUE,
               eps=0)
 
 lps_res[["gsea_GO"]] = clusterProfiler::gseGO(geneList     = lps_sorted,
                                               OrgDb        = org.Hs.eg.db,
                                               ont          = "BP",
                                               minGSSize    = 100,
                                               maxGSSize    = 500,
                                               pvalueCutoff = 0.05,
                                               verbose      = TRUE,
                                               eps=0)
 ##### miscellaneous visualizations #####
 # https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
 edo2 = gseDO(ifn_sorted)
 edox = setReadable(edo2, 'org.Hs.eg.db', 'ENTREZID')


 
 t2g = tibble::enframe(pathways,name = "gs_name",value = "entrez_gene") %>%
   tidyr::unnest(cols = c(entrez_gene))
   
 ifn = clusterProfiler::GSEA(ifn_sorted, TERM2GENE = t2g,
                             eps=0) 
 lps = clusterProfiler::GSEA(lps_sorted, TERM2GENE = t2g,
                             eps=0)
 upsetplot(ifn) 
 ridgeplot( lps_res[["gsea_GO"]])
 
 p1 = enrichplot::heatplot(ifn, foldChange=ifn_sorted, showCategory=5)
 p2 =  enrichplot::heatplot(lps, foldChange=lps_sorted, showCategory=5)
 
 png(paste0(outdir,"/msigdb_hallmark_gsea_heatmap.png"),
     width = 16, height = 6, res = 400,units = "in", type = "cairo")
p1 / p2

dev.off()

### as clusters ########
# hypergeometric test enrichments
ifn_sign_pos = ifn_vs_untreated %>%
  dplyr::filter(log2FoldChange > 1 & padj<0.05) %>%
ifn_sign_neg = ifn_vs_untreated %>%
  dplyr::filter(log2FoldChange < -1 & padj<0.05)
lps_sign_pos = lps_vs_untreated %>%
  dplyr::filter(log2FoldChange > 1 & padj<0.05)
lps_sign_neg = lps_vs_untreated %>%
  dplyr::filter(log2FoldChange < -1 & padj<0.05)

go = clusterProfiler::compareCluster(geneCluster = list("IFNvsUntreated+" = ifn_sign_pos$entrez,
                                       "IFNvsUntreated-" = ifn_sign_neg$entrez,
                                       "LPSvsUntreated+" = lps_sign_pos$entrez,
                                       "LPSvsUntreated-" = lps_sign_neg$entrez), 
                    fun = enrichGO, 
                    OrgDb='org.Hs.eg.db',
                    ont="BP", 
                    pvalueCutoff=0.01,
                    universe = ifn_vs_untreated$symbol,
                    minGSSize    = 100,
                    maxGSSize    = 500)

go = clusterProfiler::setReadable(go, 
                                  OrgDb = org.Hs.eg.db, 
                                  keyType="ENTREZID")
head(go) 


png(paste0(outdir,"/go_hypergeomtric_dotplot.png"),
    width = 10, height = 6, res = 400,units = "in", type = "cairo")
clusterProfiler::dotplot(go)
dev.off()
# pathway and KEGG don't work

 