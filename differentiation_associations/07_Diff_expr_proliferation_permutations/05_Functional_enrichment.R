
library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(fgsea)
library(sessioninfo)

#-------------------------------------------------------------------------------
#         7.5 Functional enrichment analysis for proliferation DEGs
#-------------------------------------------------------------------------------
#  Code to perform enrichment analysis for AD/PD candidate genes, and 
#  macrophage survival/ essential genes, among DEGs for cell line proliferation 
#  from macrophage precursor -> microglia in each treatment.
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

## Working dir
setwd("/lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_Daianna/")

## Output data and plot dirs 
outdir = paste(getwd(), "output_data", "07_Diff_expr_proliferation_permutations", "05_Functional_enrichment", sep = "/")
plotdir = paste(getwd(), "plots", "07_Diff_expr_proliferation_permutations", "05_Functional_enrichment", sep = "/")
dir.create(outdir, recursive = T)
dir.create(plotdir, recursive = T)

## Subdirs
sapply(paste0(outdir, c("/01_ORA", "/02_GSEA")), function(p){dir.create(p, recursive = T)})
sapply(paste0(plotdir, c("/01_ORA", "/02_GSEA")), function(p){dir.create(p, recursive = T)})

## Input dir
input_dir <- paste(getwd(), "output_data", "07_Diff_expr_proliferation_permutations", sep = "/")
#input_dir01 <- paste(getwd(), "output_data", "01_Diff_expr_PRS_permutations", "03_Functional_enrichment", sep = "/")
input_dir04 <- paste(getwd(), "output_data", "04_Burden_test_proliferation", "02_Functional_enrichment", sep = "/")

## DEGs per treatment(s)
DEGs_prolif <- read_csv(paste0(input_dir, "/04_DEGs_overlaps_explorations/DEGs_in_treatments_only_all_phenotypes.csv")) %>% 
  filter(phenotype == "proliferation")

## DGE results per treatment
load(paste0(input_dir, "/01_run_proliferation_DGE/top_genes_proliferation_IFN.Rdata"), verbose = T)
load(paste0(input_dir, "/01_run_proliferation_DGE/top_genes_proliferation_LPS.Rdata"), verbose = T)
load(paste0(input_dir, "/01_run_proliferation_DGE/top_genes_proliferation_untreated.Rdata"), verbose = T)


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
##                   7.5.1 Overrepresentation analysis (ORA)
##- - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  
## Gene universe = all expressed genes across all treatments 
geneUniverse <- unique(c(top_genes_proliferation_IFN$gene, 
                         top_genes_proliferation_LPS$gene, 
                         top_genes_proliferation_untreated$gene))

## DEG groups
DEGs_list <- list("DEGs in IFN and LPS only" = subset(DEGs_prolif, DGE_across_treatments == "DEG_in_IFN_and_LPS_only")$gene, 
                  "DEGs in IFN only" = subset(DEGs_prolif, DGE_across_treatments == "DEG_in_IFN_only")$gene,
                  "DEGs in LPS only" = subset(DEGs_prolif, DGE_across_treatments == "DEG_in_LPS_only")$gene)
names(DEGs_list) <- paste0(names(DEGs_list), " (n = ", unlist(lapply(DEGs_list, length)), ")")

## Gene sets:
# * AD/PD candidate genes  
load(paste0(input_dir01, "/AD_PD_leadCRISPR_candidate_gene_sets.Rdata"), verbose = T)
# candidate_gene_sets
candidate_gene_sets <- lapply(candidate_gene_sets, function(.){.[,"external_gene_name"]})
  
# * Essential/Viability macrophage genes  
load(paste0(input_dir04, "/Essential_and_macrophageSurvival_gene_sets.Rdata"), verbose = T)
# gene_sets
essential_surv_genes <- gene_sets
names(essential_surv_genes) <- sapply(names(essential_surv_genes), function(.){gsub("_", " ", .)})
gene_sets <- c(candidate_gene_sets, essential_surv_genes)
  

fisher_ps <- matrix(data = NA, nrow = length(gene_sets), ncol = length(DEGs_list))
fisher_OR <- matrix(data = NA, nrow = length(gene_sets), ncol = length(DEGs_list))
intersect_size <-  matrix(data = NA, nrow = length(gene_sets), ncol = length(DEGs_list))
rownames(fisher_ps) <- rownames(fisher_OR) <- rownames(intersect_size) <- names(gene_sets)
colnames(fisher_ps) <- colnames(fisher_OR) <- colnames(intersect_size) <- names(DEGs_list)

for(j in 1:length(DEGs_list)){
  
  DEGs <- DEGs_list[[j]]
  DEGs_name <- names(DEGs_list)[j]
  
  for(i in 1:length(gene_sets)){
    
    gene_set <- gene_sets[[i]]
    gene_set_name <- names(gene_sets)[i]
    
    ## Subset to genes in set that are in universe
    gene_set <- gene_set[gene_set %in% geneUniverse] %>% unique
    
    ## Define sets:
    non_DEGs <- geneUniverse[!geneUniverse %in% DEGs]
    DEGs_in_set <- DEGs[DEGs %in% gene_set]
    non_DEGs_in_set <- non_DEGs[non_DEGs %in% gene_set]
    
    DEGs_not_in_set <- DEGs[!DEGs %in% gene_set]
    non_DEGs_not_in_set <- non_DEGs[!non_DEGs %in% gene_set]
    
    ## Confirm length of universe, signif group, non-signif group, and gene set
    stopifnot(length(DEGs) + length(non_DEGs) == length(geneUniverse))
    stopifnot(length(DEGs_in_set) + length(DEGs_not_in_set) == length(DEGs))
    stopifnot(length(non_DEGs_in_set) + length(non_DEGs_not_in_set) == length(non_DEGs))
    stopifnot(length(DEGs_in_set) + length(non_DEGs_in_set) == length(gene_set))
    
    m <- matrix(data = c(length(DEGs_in_set), 
                         length(DEGs_not_in_set), 
                         length(non_DEGs_in_set), 
                         length(non_DEGs_not_in_set)), ncol = 2)
    
    f <- fisher.test(m, alternative = "greater")
    
    fisher_ps[gene_set_name, DEGs_name] <- f$p.value
    fisher_OR[gene_set_name, DEGs_name] <- f$estimate
    intersect_size[gene_set_name, DEGs_name] <- paste0(length(DEGs_in_set), "/",  length(gene_set))
    
  }
}

## Heatmap
h_row_title <- "Gene sets"
h <- ComplexHeatmap::Heatmap(-log10(as.matrix(fisher_ps)), 
                             row_title = h_row_title, 
                             row_title_gp = gpar(fontsize = 10),
                             name = "-log10(pval)", 
                             col = colorRamp2(c(0, 0.5, 1, 1.5, 2), c("white", "lavenderblush2", "lightpink", "salmon", "salmon4")), 
                             border = T, 
                             rect_gp = gpar(col = "black", lwd = 1),
                             row_names_gp = gpar(fontsize = 8.5), 
                             column_names_gp = gpar(fontsize = 6.5), 
                             row_names_max_width = unit(10, "cm"),
                             cluster_columns = F, 
                             cluster_rows = F, 
                             row_order = names(gene_sets),
                             row_labels = names(gene_sets), 
                             cell_fun = function(j, i, x, y, width, height, fill) {
                               if(log10(as.matrix(fisher_ps))[i,j] < log10(0.05)){
                                 grid.text("*", x, y, gp = gpar(fontsize = 12.5, col = "yellow"), just = c(0.5,-1.2))
                                 grid.text(paste0(signif(fisher_OR[i, j], digits = 3), "\n", "(", intersect_size[i,j], ")"), 
                                           x, y, gp = gpar(fontsize = 8.5)) 
                               }
                               else{
                                 grid.text(paste0(signif(fisher_OR[i, j], digits = 3), "\n", "(", intersect_size[i,j], ")"), x, y, gp = gpar(fontsize = 8.5))
                               }
                             })

pdf(file = paste0(plotdir, "/01_ORA", "/ORA_prolif_DEGs_in_candidate_gene_sets", ".pdf"), height = 7, width = 6)
h
dev.off()



## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
##                   7.5.2 Gene Set Enrichment Analysis (GSEA)
##- - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# - - - - - - - - - - - - - - - - IFN GSEA  - - - - - - - - - - - - - - - - - - 

## Rank genes by t-stats
ranked_genes_metrics_IFN <- top_genes_proliferation_IFN[order(top_genes_proliferation_IFN$t, decreasing = T), ]
ranked_ts_IFN <- ranked_genes_metrics_IFN$t
names(ranked_ts_IFN) <- ranked_genes_metrics_IFN$gene

fgsea_candidates_IFN <- fgsea(pathways = gene_sets, 
                              stats = ranked_ts_IFN,
                              minSize = 1,  
                              maxSize = Inf,  
                              eps = 0) %>%
  tidyr::as_tibble() %>%
  dplyr::arrange(desc(NES)) %>% as.data.frame()

rownames(fgsea_candidates_IFN) <- fgsea_candidates_IFN$pathway
save(fgsea_candidates_IFN, file = paste0(outdir, "/02_GSEA", "/Candidates_GSEA_results_", "IFN", ".Rdata"))

## Results
fgsea_candidates_IFN 
#                                   pathway         pval         padj    log2err         ES        NES size  leadingEdge
# 1       Covarrubias MAGeCK macro survival 7.247311e-16 1.630645e-15 1.02766987  0.4066371  2.2806862  346 PSMB4, C....
# 2                  DEPMAP essential genes 1.495912e-36 1.346321e-35 1.57637362  0.3626277  2.2410569 1174 ADSL, ER....
# 3              Covarrubias macro survival 3.170539e-32 9.511616e-32 1.47456250  0.2793122  1.8069286 2841 ADSL, PD....
# 4 Covarrubias Mann Whitney macro survival 2.478937e-32 9.511616e-32 1.47456250  0.2773574  1.7916075 2818 ADSL, PD....
# 5           Lead immune CRISPR (Kampmann) 3.414097e-01 5.036101e-01 0.10672988  0.3008232  1.0840265   27 DCPS, EN....
# 6     Lead phagocytosis CRISPR (Kampmann) 3.572985e-01 5.036101e-01 0.10319747  0.2980894  1.0642248   26 DCPS, CS....
# 7                           PD candidates 5.277778e-01 5.937500e-01 0.07977059  0.2829595  0.9617895   22 RAB29, V....
# 8                           AD candidates 9.309091e-01 9.309091e-01 0.04459616 -0.1763750 -0.7184821   60 APOE, RA....
# 9                Lead phagocytosis CRISPR 3.916968e-01 5.036101e-01 0.08679498 -0.2993862 -1.0483298   28 PRKN, AT....

plotGseaTable(pathways = gene_sets, 
              stats = ranked_ts_IFN, 
              fgseaRes = fgsea_candidates_IFN)
ggsave(filename = paste0(plotdir, "/02_GSEA", "/Results_fgsea_", "candidate_gene_sets_IFN", ".pdf"), height = 4, width = 10)


# - - - - - - - - - - - - - - - - LPS GSEA  - - - - - - - - - - - - - - - - - - 

## Rank genes by t-stats
ranked_genes_metrics_LPS <- top_genes_proliferation_LPS[order(top_genes_proliferation_LPS$t, decreasing = T), ]
ranked_ts_LPS <- ranked_genes_metrics_LPS$t
names(ranked_ts_LPS) <- ranked_genes_metrics_LPS$gene

fgsea_candidates_LPS <- fgsea(pathways = gene_sets, 
                              stats = ranked_ts_LPS,
                              minSize = 1,  
                              maxSize = Inf,  
                              eps = 0) %>%
  tidyr::as_tibble() %>%
  dplyr::arrange(desc(NES)) %>% as.data.frame()

rownames(fgsea_candidates_LPS) <- fgsea_candidates_LPS$pathway
save(fgsea_candidates_LPS, file = paste0(outdir, "/02_GSEA", "/Candidates_GSEA_results_", "LPS", ".Rdata"))

## Results
fgsea_candidates_LPS 
#                                   pathway         pval         padj    log2err         ES        NES size  leadingEdge
# 1                  DEPMAP essential genes 2.377324e-67 2.139592e-66 2.14409993  0.4454157  2.7029087 1180 MCM2, VP....
# 2       Covarrubias MAGeCK macro survival 2.166703e-21 4.875081e-21 1.19534448  0.4476775  2.4544827  349 SMC3, MB....
# 3              Covarrubias macro survival 1.154959e-46 5.197316e-46 1.77997643  0.3069551  1.9197405 2851 MCM2, CD....
# 4 Covarrubias Mann Whitney macro survival 6.053354e-46 1.816006e-45 1.76830433  0.3057414  1.9141314 2827 MCM2, CD....
# 5           Lead immune CRISPR (Kampmann) 9.638554e-01 9.638554e-01 0.04049348 -0.1854827 -0.6108491   27 MKNK1, L....
# 6                Lead phagocytosis CRISPR 7.744755e-01 8.712850e-01 0.05121843 -0.2356992 -0.7846941   28 CASS4, T....
# 7     Lead phagocytosis CRISPR (Kampmann) 6.898955e-01 8.712850e-01 0.05641184 -0.2643135 -0.8533282   26 MKNK1, P....
# 8                           AD candidates 4.491115e-01 6.736672e-01 0.07362127 -0.2579697 -1.0145388   61 APOE, CA....
# 9                           PD candidates 1.880492e-01 3.384886e-01 0.13214726 -0.3997854 -1.2418949   22 ATP13A2,....

plotGseaTable(pathways = gene_sets, 
              stats = ranked_ts_LPS, 
              fgseaRes = fgsea_candidates_LPS)
ggsave(filename = paste0(plotdir, "/02_GSEA", "/Results_fgsea_", "candidate_gene_sets_LPS", ".pdf"), height = 4, width = 10)


# - - - - - - - - - - - - - - - - untreated GSEA  - - - - - - - - - - - - - - - - - - 

## Rank genes by t-stats
ranked_genes_metrics_untreated <- top_genes_proliferation_untreated[order(top_genes_proliferation_untreated$t, decreasing = T), ]
ranked_ts_untreated <- ranked_genes_metrics_untreated$t
names(ranked_ts_untreated) <- ranked_genes_metrics_untreated$gene

fgsea_candidates_untreated <- fgsea(pathways = gene_sets, 
                              stats = ranked_ts_untreated,
                              minSize = 1,  
                              maxSize = Inf,  
                              eps = 0) %>%
  tidyr::as_tibble() %>%
  dplyr::arrange(desc(NES)) %>% as.data.frame()

rownames(fgsea_candidates_untreated) <- fgsea_candidates_untreated$pathway
save(fgsea_candidates_untreated, file = paste0(outdir, "/02_GSEA", "/Candidates_GSEA_results_", "untreated", ".Rdata"))

## Results
fgsea_candidates_untreated 
#                                   pathway         pval         padj    log2err         ES        NES size  leadingEdge
# 1                  DEPMAP essential genes 3.703093e-56 1.110928e-55 1.95728650  0.4284668  2.5809382 1180 HNRNPM, ....
# 2       Covarrubias MAGeCK macro survival 1.389625e-19 3.126655e-19 1.14219120  0.4399751  2.4087408  350 SRSF1, A....
# 3              Covarrubias macro survival 8.895482e-60 5.431110e-59 2.01976875  0.3395181  2.1420101 2865 HNRNPM, ....
# 4 Covarrubias Mann Whitney macro survival 1.206913e-59 5.431110e-59 2.01463591  0.3397141  2.1390197 2841 HNRNPM, ....
# 5     Lead phagocytosis CRISPR (Kampmann) 6.848960e-03 1.232813e-02 0.40701792  0.5079770  1.7144559   26 DCPS, DD....
# 6           Lead immune CRISPR (Kampmann) 4.015595e-01 5.412196e-01 0.08971047  0.3016618  1.0233000   27 DCPS, HP....
# 7                           PD candidates 4.209486e-01 5.412196e-01 0.08783126  0.3142913  1.0104887   22 RAB29, P....
# 8                Lead phagocytosis CRISPR 8.895706e-01 8.895706e-01 0.05205700 -0.2008931 -0.6965698   28 DCAF7, I....
# 9                           AD candidates 7.290837e-01 8.202191e-01 0.06011861 -0.2049794 -0.8594300   62 APOE, AP....

plotGseaTable(pathways = gene_sets, 
              stats = ranked_ts_untreated, 
              fgseaRes = fgsea_candidates_untreated)
ggsave(filename = paste0(plotdir, "/02_GSEA", "/Results_fgsea_", "candidate_gene_sets_untreated", ".pdf"), height = 4, width = 10)







## Reproducibility info
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.1 (2023-06-16)
# os       Ubuntu 22.04.4 LTS
# system   x86_64, linux-gnu
# ui       RStudio
# language (EN)
# collate  en_GB.UTF-8
# ctype    en_GB.UTF-8
# tz       Europe/Belfast
# date     2025-02-04
# rstudio  2024.04.0+735 Chocolate Cosmos (server)
# pandoc   3.1.12.3 @ /opt/view/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version     date (UTC) lib source
# abind                  1.4-5       2016-07-21 [1] CRAN (R 4.3.1)
# AnnotationDbi        * 1.64.1      2023-11-03 [1] Bioconductor
# Biobase              * 2.62.0      2023-10-24 [1] Bioconductor
# BiocFileCache          2.10.2      2024-03-27 [1] Bioconductor 3.18 (R 4.3.1)
# BiocGenerics         * 0.48.1      2023-11-01 [1] Bioconductor
# BiocParallel           1.36.0      2023-10-24 [1] Bioconductor
# biomaRt              * 2.58.2      2024-01-30 [1] Bioconductor 3.18 (R 4.3.1)
# Biostrings             2.70.3      2024-03-13 [1] Bioconductor 3.18 (R 4.3.1)
# bit                    4.0.5       2022-11-15 [1] CRAN (R 4.3.1)
# bit64                  4.0.5       2020-08-30 [1] CRAN (R 4.3.1)
# bitops                 1.0-7       2021-04-24 [1] CRAN (R 4.3.1)
# blob                   1.2.4       2023-03-17 [1] CRAN (R 4.3.1)
# cachem                 1.0.8       2023-05-01 [1] CRAN (R 4.3.1)
# Cairo                  1.6-2       2023-11-28 [1] CRAN (R 4.3.1)
# circlize             * 0.4.16      2024-02-20 [1] CRAN (R 4.3.1)
# cli                    3.6.2       2023-12-11 [1] CRAN (R 4.3.1)
# clue                   0.3-65      2023-09-23 [1] CRAN (R 4.3.1)
# cluster                2.1.6       2023-12-01 [1] CRAN (R 4.3.1)
# codetools              0.2-20      2024-03-31 [1] CRAN (R 4.3.1)
# colorspace             2.1-0       2023-01-23 [1] CRAN (R 4.3.1)
# ComplexHeatmap       * 2.18.0      2023-10-24 [1] Bioconductor
# cowplot              * 1.1.3       2024-01-22 [1] CRAN (R 4.3.1)
# crayon                 1.5.2       2022-09-29 [1] CRAN (R 4.3.1)
# curl                   5.2.1       2024-03-01 [1] CRAN (R 4.3.1)
# data.table             1.15.4      2024-03-30 [1] CRAN (R 4.3.1)
# DBI                    1.2.2       2024-02-16 [1] CRAN (R 4.3.1)
# dbplyr                 2.5.0       2024-03-19 [1] CRAN (R 4.3.1)
# DelayedArray           0.28.0      2023-10-24 [1] Bioconductor
# digest                 0.6.35      2024-03-11 [1] CRAN (R 4.3.1)
# doParallel             1.0.17      2022-02-07 [1] CRAN (R 4.3.1)
# dplyr                * 1.1.4       2023-11-17 [1] CRAN (R 4.3.1)
# fansi                  1.0.6       2023-12-08 [1] CRAN (R 4.3.1)
# farver                 2.1.1       2022-07-06 [1] CRAN (R 4.3.1)
# fastmap                1.1.1       2023-02-24 [1] CRAN (R 4.3.1)
# fastmatch              1.1-4       2023-08-18 [1] CRAN (R 4.3.1)
# fgsea                * 1.28.0      2023-10-24 [1] Bioconductor
# filelock               1.0.3       2023-12-11 [1] CRAN (R 4.3.1)
# forcats              * 1.0.0       2023-01-29 [1] CRAN (R 4.3.1)
# foreach                1.5.2       2022-02-02 [1] CRAN (R 4.3.1)
# generics               0.1.3       2022-07-05 [1] CRAN (R 4.3.1)
# GenomeInfoDb         * 1.38.8      2024-03-15 [1] Bioconductor 3.18 (R 4.3.1)
# GenomeInfoDbData       1.2.11      2024-11-06 [1] Bioconductor
# GenomicRanges        * 1.54.1      2023-10-29 [1] Bioconductor
# GetoptLong             1.0.5       2020-12-15 [1] CRAN (R 4.3.1)
# ggplot2              * 3.5.1       2024-04-23 [1] CRAN (R 4.3.1)
# ggrepel                0.9.5       2024-01-10 [1] CRAN (R 4.3.1)
# GlobalOptions          0.1.2       2020-06-10 [1] CRAN (R 4.3.1)
# glue                   1.7.0       2024-01-09 [1] CRAN (R 4.3.1)
# gtable                 0.3.4       2023-08-21 [1] CRAN (R 4.3.1)
# hms                    1.1.3       2023-03-21 [1] CRAN (R 4.3.1)
# httr                   1.4.7       2023-08-15 [1] CRAN (R 4.3.1)
# IRanges              * 2.36.0      2023-10-24 [1] Bioconductor
# iterators              1.0.14      2022-02-05 [1] CRAN (R 4.3.1)
# jtools                 2.2.2       2023-07-11 [1] CRAN (R 4.3.1)
# KEGGREST               1.42.0      2023-10-24 [1] Bioconductor
# labeling               0.4.3       2023-08-29 [1] CRAN (R 4.3.1)
# lattice                0.22-6      2024-03-20 [1] CRAN (R 4.3.1)
# lifecycle              1.0.4       2023-11-07 [1] CRAN (R 4.3.1)
# lubridate            * 1.9.3       2023-09-27 [1] CRAN (R 4.3.1)
# magick                 2.8.3       2024-02-18 [1] CRAN (R 4.3.1)
# magrittr               2.0.3       2022-03-30 [1] CRAN (R 4.3.1)
# Matrix                 1.6-5       2024-01-11 [1] CRAN (R 4.3.1)
# MatrixGenerics       * 1.14.0      2023-10-24 [1] Bioconductor
# matrixStats          * 1.2.0       2023-12-11 [1] CRAN (R 4.3.1)
# memoise                2.0.1       2021-11-26 [1] CRAN (R 4.3.1)
# munsell                0.5.1       2024-04-01 [1] CRAN (R 4.3.1)
# org.Hs.eg.db         * 3.18.0      2024-11-06 [1] Bioconductor
# pander                 0.6.5       2022-03-18 [1] CRAN (R 4.3.1)
# pillar                 1.9.0       2023-03-22 [1] CRAN (R 4.3.1)
# pkgconfig              2.0.3       2019-09-22 [1] CRAN (R 4.3.1)
# png                    0.1-8       2022-11-29 [1] CRAN (R 4.3.1)
# prettyunits            1.2.0       2023-09-24 [1] CRAN (R 4.3.1)
# progress               1.2.3       2023-12-06 [1] CRAN (R 4.3.1)
# purrr                * 1.0.2       2023-08-10 [1] CRAN (R 4.3.1)
# R6                     2.5.1       2021-08-19 [1] CRAN (R 4.3.1)
# ragg                   1.3.0       2024-03-13 [1] CRAN (R 4.3.1)
# rappdirs               0.3.3       2021-01-31 [1] CRAN (R 4.3.1)
# RColorBrewer           1.1-3       2022-04-03 [1] CRAN (R 4.3.1)
# Rcpp                   1.0.12      2024-01-09 [1] CRAN (R 4.3.1)
# RCurl                  1.98-1.14   2024-01-09 [1] CRAN (R 4.3.1)
# readr                * 2.1.5       2024-01-10 [1] CRAN (R 4.3.1)
# rjson                  0.2.21      2022-01-09 [1] CRAN (R 4.3.1)
# rlang                * 1.1.3       2024-01-10 [1] CRAN (R 4.3.1)
# RSQLite                2.3.6       2024-03-31 [1] CRAN (R 4.3.1)
# rstudioapi             0.16.0      2024-03-24 [1] CRAN (R 4.3.1)
# S4Arrays               1.2.1       2024-03-04 [1] Bioconductor 3.18 (R 4.3.1)
# S4Vectors            * 0.40.2      2023-11-23 [1] Bioconductor 3.18 (R 4.3.1)
# scales                 1.3.0       2023-11-28 [1] CRAN (R 4.3.1)
# sessioninfo          * 1.2.2       2021-12-06 [1] CRAN (R 4.3.1)
# shape                  1.4.6.1     2024-02-23 [1] CRAN (R 4.3.1)
# SparseArray            1.2.4       2024-02-11 [1] Bioconductor 3.18 (R 4.3.1)
# stringi                1.8.3       2023-12-11 [1] CRAN (R 4.3.1)
# stringr              * 1.5.1       2023-11-14 [1] CRAN (R 4.3.1)
# SummarizedExperiment * 1.32.0      2023-10-24 [1] Bioconductor
# systemfonts            1.0.6       2024-03-07 [1] CRAN (R 4.3.1)
# textshaping            0.3.7       2023-10-09 [1] CRAN (R 4.3.1)
# tibble               * 3.2.1       2023-03-20 [1] CRAN (R 4.3.1)
# tidyr                * 1.3.1       2024-01-24 [1] CRAN (R 4.3.1)
# tidyselect             1.2.1       2024-03-11 [1] CRAN (R 4.3.1)
# tidyverse            * 2.0.0       2023-02-22 [1] CRAN (R 4.3.1)
# timechange             0.3.0       2024-01-18 [1] CRAN (R 4.3.1)
# tzdb                   0.4.0       2023-05-12 [1] CRAN (R 4.3.1)
# utf8                   1.2.4       2023-10-22 [1] CRAN (R 4.3.1)
# vctrs                  0.6.5       2023-12-01 [1] CRAN (R 4.3.1)
# vroom                  1.6.5       2023-12-05 [1] CRAN (R 4.3.1)
# withr                  3.0.0       2024-01-16 [1] CRAN (R 4.3.1)
# XML                    3.99-0.16.1 2024-01-22 [1] CRAN (R 4.3.1)
# xml2                   1.3.6       2023-12-04 [1] CRAN (R 4.3.1)
# XVector                0.42.0      2023-10-24 [1] Bioconductor
# zlibbioc               1.48.2      2024-03-13 [1] Bioconductor 3.18 (R 4.3.1)
# 
# [1] /opt/view/rlib/R/library
# [2] /software/hgi/installs/softpack/rstudio/.spack/linux-ubuntu22.04-x86_64_v3/gcc-11.4.0/r-4.3.1-bfwldrk76z6f52upk47zepliekn7ayqz/rlib/R/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────


