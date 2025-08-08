# checking certain genes
# familial AD, AD candidate genes enrichments

library(tidyverse)
library(patchwork)
library(fgsea)
source("./helpers.R")
library(AnnotationDbi)
library(org.Hs.eg.db)
library(pheatmap)
library(ggrepel)

outdir = "../../data/results/5.3.enrichment_candidate_genes"
dir.create(outdir,recursive = TRUE)
res_additive = readRDS("../../data/results/5.1.Diff_expr_limma/additive_res.rds")
for(cont in names(res_additive)){
  res_additive[[cont]] =   res_additive[[cont]]%>%
    dplyr::mutate(diffexp = case_when(logFC > 1 & adj.P.Val < 0.05 ~ "UP",
                                      logFC < -1 & adj.P.Val < 0.05 ~ "DOWN",
                                      .default = "NO")) 
}


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
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),breaks = c(0,10,50,100,200,300,400)) +
  ggrepel::geom_text_repel(color = "darkgrey") +
  scale_color_manual(values = c("DOWN"="#53B3CB",
                                "NO"= "#A7A9AE",
                                "UP"="#D56062")) +

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
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),breaks = c(0,10,50,100,200,300,400)) +
  ggrepel::geom_text_repel(color = "darkgrey") +
  scale_color_manual(values = c("DOWN"="#53B3CB",
                                "NO"= "#A7A9AE",
                                "UP"="#D56062")) +
  
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
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),breaks = c(0,10,50,100,200,300,400)) +
  ggrepel::geom_text_repel(color = "darkgrey") +
  scale_color_manual(values = c("DOWN"="#53B3CB",
                                "NO"= "#A7A9AE",
                                "UP"="#D56062")) +
  
  geom_vline(xintercept = c(-1, 1),linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 2,linetype = "dashed", color = "grey") +
  ylab("-log10(adj. p-val)") + 
  labs(fill = "Diff. expr.") +
  ggtitle("IFN vs LPS")

png(paste0(outdir,"/volcanos_res_additive.png"),
    width = 20, height = 5, res = 400,units = "in", type = "cairo")
patchwork::wrap_plots(p1,p2,p3, guides = "collect")
dev.off()

## jeremy genes
jeremy_ad_candidates = read_table("../../../resources/Jeremy_medrXiv_AD_candidate_genes.txt")

## heatmap of AD candiate genes
# signal significant differentially expressed genes

pheatmap::pheatmap()

ifn_sorted = ifn_vs_untreated$logFC
names(ifn_sorted) = ifn_vs_untreated$symbol
ifn_sorted = sort(ifn_sorted,decreasing = TRUE)
ifn_sorted = ifn_sorted[!is.na(names(ifn_sorted))]
ifn_vs_untreated$entrez = symbol_to_entrez(ifn_vs_untreated$symbol)

lps_sorted = lps_vs_untreated$logFC
names(lps_sorted) = lps_vs_untreated$symbol
lps_sorted = sort(lps_sorted,decreasing = TRUE)
lps_sorted = lps_sorted[!is.na(names(lps_sorted))]
lps_vs_untreated$entrez = symbol_to_entrez(lps_vs_untreated$symbol)

pathways = list()
pathways[["jeremy_ad_candidates"]] = jeremy_ad_candidates$symbol

### GSEA

ifn_res = list()
lps_res = list()

ifn_res[["gsea_jeremy_ad_candidates"]] = fgsea::fgsea(pathways=pathways, 
                                                      stats=ifn_sorted,
                                                      eps=0) %>%
  tidyr::as_tibble() %>%
  dplyr::arrange(desc(NES))

lps_res[["gsea_jeremy_ad_candidates"]] = fgsea::fgsea(pathways=pathways, 
                                                      stats=lps_sorted,
                                                      eps=0) %>%
  tidyr::as_tibble() %>%
  dplyr::arrange(desc(NES))

# check hypergeomtric
