
library(tidyverse)
library(dplyr)
library(biomaRt)
library(readr)
library(rlang)
library(ggrepel)
library(AnnotationHub)
library(locuszoomr)
library(sessioninfo)


#-------------------------------------------------------------------------------
#                    9.3 Check proliferation TWMR results 
#-------------------------------------------------------------------------------
#  Code to examine and plot the gene-wise TWMR results for proliferation in 
#  each treatment.  
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

## Define dirs
input_dir = paste("..", "..", "input_data", "09_TWMR_proliferation", "03_check_results", sep = "/")
plotdir = paste("..", "..", "plots", "09_TWMR_proliferation", "03_check_results", sep = "/")

dir.create(input_dir, recursive = T, showWarnings = F)
dir.create(plotdir, recursive = T, showWarnings = F)

## Input dirs
input_dir0701 = paste("..", "..", "output_data", "07_Diff_expr_proliferation_permutations", "01_run_proliferation_DGE", sep = "/")
input_dir0802 = paste("..", "..", "output_data", "08_GWAS_proliferation", "02_Check_association_results", sep = "/")
input_dir0803 = paste(getwd(), "output_data", "08_GWAS_proliferation", "03_Closest_genes_explorations", sep = "/")
input_dir01 = paste("..", "..", "input_data", "09_TWMR_proliferation", "01_IVs_and_exposures_selection", sep = "/")

## Load and process results x treatment 
for(treatment in c("IFN", "LPS", "untreated")){
  
  results_dir <- paste("..", "..", "output_data", "09_TWMR_proliferation", "02_run_TWMR", treatment, sep = "/") 
  filenames <- list.files(results_dir, full.names = T)
  genes <- str_split_fixed(list.files(results_dir), pattern = "_", n = 2)[,1]
  results <- lapply(filenames, function(gene_path){as.data.frame(get(load(gene_path)))}) %>% do.call(rbind, .)
  results <- cbind("gene_symbol" = genes, results)
  
  ## Add Bonferroni corrected pval
  results$p_adj <- sapply(results$P, function(P){min(1, P * dim(results)[1])})
  
  ## Add gene ensembl IDs
  mart <- useEnsembl(biomart="ENSEMBL_MART_ENSEMBL",
                     dataset="hsapiens_gene_ensembl",
                     GRCh = "GRCh38")
  
  ids <- getBM(values = unique(results$gene_symbol),
               mart = mart,
               attributes = c("external_gene_name", "ensembl_gene_id"),
               filters = "external_gene_name") %>% 
    rename(gene_symbol = external_gene_name,
           gene_id = ensembl_gene_id)
  
  results = results %>% left_join(ids, multiple = "first") 
  
  assign(paste0(treatment, "_results_all_genes"), results)
  
}

save(IFN_results_all_genes, file = paste0(outdir, "/IFN_results_all_genes.Rdata"))
save(LPS_results_all_genes, file = paste0(outdir, "/LPS_results_all_genes.Rdata"))
save(untreated_results_all_genes, file = paste0(outdir, "/untreated_results_all_genes.Rdata"))


## Subset only to genes with at least one signif cis-eQTL (q<0.05)
tensorQTL_variant_gene <- read_csv("~/OTAR2065_sc_eQTL/data/results/8.5.eQTL_MR/TWMR/input/tensorQTL_variant_gene.csv")
signif_genes_IFN <- tensorQTL_variant_gene %>% filter(treatment == "IFN" & qval < 0.05) %>% select(gene_name) %>% unlist %>% unique
length(signif_genes_IFN)
# [1] 4450
signif_genes_LPS <- tensorQTL_variant_gene %>% filter(treatment == "LPS" & qval < 0.05) %>% select(gene_name) %>% unlist %>% unique
length(signif_genes_LPS)
# [1] 4521
signif_genes_untreated <- tensorQTL_variant_gene %>% filter(treatment == "untreated" & qval < 0.05) %>% select(gene_name) %>% unlist %>% unique
length(signif_genes_untreated)
# [1] 5182

## Subset and re-compute Bonf-corrected pval with new num of tested genes
IFN_results_genes_with_signif_eQTLs <- IFN_results_all_genes %>% filter(gene_symbol %in% signif_genes_IFN) %>% 
  mutate(p_adj = pmin(1, P * dim(.)[1]))
LPS_results_genes_with_signif_eQTLs <- LPS_results_all_genes %>% filter(gene_symbol %in% signif_genes_LPS) %>% 
  mutate(p_adj = pmin(1, P * dim(.)[1]))
untreated_results_genes_with_signif_eQTLs <- untreated_results_all_genes %>% filter(gene_symbol %in% signif_genes_untreated) %>% 
  mutate(p_adj = pmin(1, P * dim(.)[1]))

save(IFN_results_genes_with_signif_eQTLs, file = paste0(outdir, "/IFN_results_genes_with_signif_eQTLs.Rdata"))
save(LPS_results_genes_with_signif_eQTLs, file = paste0(outdir, "/LPS_results_genes_with_signif_eQTLs.Rdata"))
save(untreated_results_genes_with_signif_eQTLs, file = paste0(outdir, "/untreated_results_genes_with_signif_eQTLs.Rdata"))



## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#                       9.3.1 Explore gene TWMR models
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
##    All genes:
###################################  IFN  ######################################  
dim(IFN_results_all_genes)
# [1] 10548     9

pdf(file = paste0(plotdir, "/results_TWMR_", "IFN", "_all_genes.pdf"), height = 9, width = 7)
par(mfrow=c(3,2))

summary(IFN_results_all_genes$Nsnps)
hist(IFN_results_all_genes$Nsnps)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   1.000   4.031   1.000 243.000 

summary(IFN_results_all_genes$Ngene)
hist(IFN_results_all_genes$Ngene)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   2.000   3.000   3.178   4.000  13.000 

summary(IFN_results_all_genes$P)
hist(IFN_results_all_genes$P)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000056 0.4360296 0.6703471 0.6211438 0.8362741 1.0000000 

summary(IFN_results_all_genes$p_adj)
hist(IFN_results_all_genes$p_adj)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.05925 1.00000 1.00000 0.99962 1.00000 1.00000 

summary(IFN_results_all_genes$Z)
hist(IFN_results_all_genes$Z)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -4.540307 -0.426814  0.000000  0.001848  0.424368  4.357610 

dev.off()

###################################  LPS  ######################################  
dim(LPS_results_all_genes)
# [1] 10600     9

pdf(file = paste0(plotdir, "/results_TWMR_", "LPS", "_all_genes.pdf"), height = 9, width = 7)
par(mfrow=c(3,2))

summary(LPS_results_all_genes$Nsnps)
hist(LPS_results_all_genes$Nsnps)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   1.000   3.496   1.000 273.000 

summary(LPS_results_all_genes$Ngene)
hist(LPS_results_all_genes$Ngene)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   2.000   3.000   3.202   4.000  14.000 

summary(LPS_results_all_genes$P)
hist(LPS_results_all_genes$P)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000019 0.4450898 0.6710373 0.6244169 0.8408080 1.0000000 

summary(LPS_results_all_genes$p_adj)
hist(LPS_results_all_genes$p_adj)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.02033 1.00000 1.00000 0.99977 1.00000 1.00000 

summary(LPS_results_all_genes$Z)
hist(LPS_results_all_genes$Z)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -4.761863 -0.418258  0.000000  0.002837  0.431223  4.012199 

dev.off()

################################  untreated  ###################################  
dim(untreated_results_all_genes)
# [1] 10895     9

pdf(file = paste0(plotdir, "/results_TWMR_", "untreated", "_all_genes.pdf"), height = 9, width = 7)
par(mfrow=c(3,2))

summary(untreated_results_all_genes$Nsnps)
hist(untreated_results_all_genes$Nsnps)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   1.000   4.577   2.000 280.000 

summary(untreated_results_all_genes$Ngene)
hist(untreated_results_all_genes$Ngene)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   2.000   3.000   3.223   4.000  12.000 

summary(untreated_results_all_genes$P)
hist(untreated_results_all_genes$P)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.4213  0.6692  0.6163  0.8450  1.0000 

summary(untreated_results_all_genes$p_adj)
hist(untreated_results_all_genes$p_adj)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000355 1.000000 1.000000 0.999053 1.000000 1.000000 

summary(untreated_results_all_genes$Z)
hist(untreated_results_all_genes$Z)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -5.526974 -0.424792  0.003628 -0.001065  0.428930  4.509917 

dev.off()


##    Genes with signif. cis-eQTLs: 
###################################  IFN  ######################################  
dim(IFN_results_genes_with_signif_eQTLs)
# [1] 4446     9

pdf(file = paste0(plotdir, "/results_TWMR_", "IFN", "_genes_with_signif_eQTLs.pdf"), height = 9, width = 7)
par(mfrow=c(3,2))

summary(IFN_results_genes_with_signif_eQTLs$Nsnps)
hist(IFN_results_genes_with_signif_eQTLs$Nsnps)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   1.000   4.351   1.000 170.000 

summary(IFN_results_genes_with_signif_eQTLs$Ngene)
hist(IFN_results_genes_with_signif_eQTLs$Ngene)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   2.000   3.000   3.048   4.000  13.000 

summary(IFN_results_genes_with_signif_eQTLs$P)
hist(IFN_results_genes_with_signif_eQTLs$P)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000056 0.3892760 0.6304431 0.5907788 0.8092184 1.0000000 

summary(IFN_results_genes_with_signif_eQTLs$p_adj)
hist(IFN_results_genes_with_signif_eQTLs$p_adj)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.02497 1.00000 1.00000 0.99914 1.00000 1.00000 

summary(IFN_results_genes_with_signif_eQTLs$Z)
hist(IFN_results_genes_with_signif_eQTLs$Z)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -4.540307 -0.490773  0.000000 -0.006466  0.466374  4.046308 

dev.off()

###################################  LPS  ######################################  
dim(LPS_results_genes_with_signif_eQTLs)
# [1] 4512      9

pdf(file = paste0(plotdir, "/results_TWMR_", "LPS", "_genes_with_signif_eQTLs.pdf"), height = 9, width = 7)
par(mfrow=c(3,2))

summary(LPS_results_genes_with_signif_eQTLs$Nsnps)
hist(LPS_results_genes_with_signif_eQTLs$Nsnps)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   1.000   3.548   1.000 146.000 

summary(LPS_results_genes_with_signif_eQTLs$Ngene)
hist(LPS_results_genes_with_signif_eQTLs$Ngene)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   2.000   3.000   2.995   4.000  12.000 

summary(LPS_results_genes_with_signif_eQTLs$P)
hist(LPS_results_genes_with_signif_eQTLs$P)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000333 0.3947785 0.6311792 0.5904866 0.8092899 1.0000000 

summary(LPS_results_genes_with_signif_eQTLs$p_adj)
hist(LPS_results_genes_with_signif_eQTLs$p_adj)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1502  1.0000  1.0000  0.9993  1.0000  1.0000 

summary(LPS_results_genes_with_signif_eQTLs$Z)
hist(LPS_results_genes_with_signif_eQTLs$Z)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -4.149714 -0.481415  0.000000 -0.000971  0.479975  4.012199  

dev.off()

################################  untreated  ################################### 
dim(untreated_results_genes_with_signif_eQTLs)
# [1] 5170     9

pdf(file = paste0(plotdir, "/results_TWMR_", "untreated", "_genes_with_signif_eQTLs.pdf"), height = 9, width = 7)
par(mfrow=c(3,2))

summary(untreated_results_genes_with_signif_eQTLs$Nsnps)
hist(untreated_results_genes_with_signif_eQTLs$Nsnps)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   1.000   4.629   2.000 105.000 

summary(untreated_results_genes_with_signif_eQTLs$Ngene)
hist(untreated_results_genes_with_signif_eQTLs$Ngene)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   2.000   3.000   3.074   4.000  12.000  

summary(untreated_results_genes_with_signif_eQTLs$P)
hist(untreated_results_genes_with_signif_eQTLs$P)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000008 0.3713357 0.6225835 0.5838084 0.8223386 1.0000000 

summary(untreated_results_genes_with_signif_eQTLs$p_adj)
hist(untreated_results_genes_with_signif_eQTLs$p_adj)
#     Min.  1st Qu. Median    Mean  3rd Qu.   Max. 
#  0.00433 1.00000 1.00000 0.99875 1.00000 1.00000  

summary(untreated_results_genes_with_signif_eQTLs$Z)
hist(untreated_results_genes_with_signif_eQTLs$Z)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -4.92641 -0.50954  0.00000 -0.01843  0.47200  4.31182 

dev.off()



## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#                  9.3.2 Explore genome-wide results
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

## Function for combined plot of GWAS + TWMR results   
miami_plotting <- function(treatment, TWMR_subset){
  
  ## GWAS summary stats
  GWAS_stats <- get(load(paste0(input_dir0802, "/summary_stats_microglia_vs_premac_", treatment, ".Rdata")))
  
  ## TWMR results
  if(TWMR_subset == "all genes"){
    TWMR_stats <- eval(parse_expr(paste0(treatment, "_results_all_genes")))
  } else{
    TWMR_stats <- eval(parse_expr(paste0(treatment, "_results_genes_with_signif_eQTLs")))
  }
    
  
  #x x x x x x x x x x x x x x x x   GWAS plot   x x x x x x x x x x x x x x x x 
  
  ## Add row num for filtering for plotting
  GWAS_stats = GWAS_stats %>% 
    rename(SNP = variant_id) %>% 
    arrange(P) %>% 
    mutate(row_num = 1:nrow(GWAS_stats))
  
  ## Hide n.s. points for plotting:
  ## Keep all variants with P<=0.01 and only 50% variants with P>0.01
  GWAS_stats$keep_row = apply(GWAS_stats, 1, 
                                 function(x){if(as.numeric(x["P"]) > 0.01 & as.numeric(x["row_num"]) %% 2 == 1){FALSE} else{TRUE}})
  GWAS_stats_subset <- GWAS_stats %>% filter(keep_row == TRUE) %>% mutate(row_num = 1:nrow(.))
  
  ## Keep all variants with P<=0.05 and only 10% of variants with P>0.05
  GWAS_stats_subset$keep_row2 = apply(GWAS_stats_subset, 1, 
                                         function(x){if(as.numeric(x["P"]) > 0.05 & as.numeric(x["row_num"]) %% 10 != 0){FALSE} else{TRUE}})
  GWAS_stats_subset2 <- GWAS_stats_subset %>% filter(keep_row2 == TRUE) %>% mutate(row_num = 1:nrow(.))
  
  ## Keep all variants with P<=0.3 and only 0.3% of variants with P>0.3
  GWAS_stats_subset2$keep_row3 = apply(GWAS_stats_subset2, 1, 
                                          function(x){if(as.numeric(x["P"]) > 0.3 & as.numeric(x["row_num"]) %% 30 != 0){FALSE} else{TRUE}})
  GWAS_stats_subset3 <- GWAS_stats_subset2 %>% filter(keep_row3 == TRUE) 
  
  ## Order variants by BP per CHR 
  GWAS_stats_subset3 = GWAS_stats_subset3 %>% 
    arrange(as.numeric(CHR), as.numeric(BP)) %>% 
    mutate(plot_pos = 1:nrow(GWAS_stats_subset3))
  
  ## Signif line with Bonferroni correction
  genomewideline <- 0.05 / (10^6)  ## (6x10^6 tests actually)  
  
  #x x x x x x x x x x x x x x x x   TWMR plot   x x x x x x x x x x x x x x x x 
    
  ## Add gene positions
  TWMR_stats_genes_pos = readr::read_csv("/lustre/scratch123/hgi/teams/trynka/resources/biomart/Homo_sapiens.GRCh38.111.genes.csv") %>% 
    dplyr::rename(CHR = seqname, BP = start, gene_symbol = gene_name) %>%
    dplyr::right_join(TWMR_stats) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(CHR = as.numeric(CHR),
                  BP = as.numeric(BP)) %>% 
    arrange(CHR, BP)
  
  dim(TWMR_stats_genes_pos) # we lose some genes without position info
  # All genes plotted in IFN: 9986
  # All genes plotted in LPS: 10034 
  # All genes plotted in untreated: 10311
  
  # Genes with signif eQTLs plotted in IFN: 4157
  # Genes with signif eQTLs plotted in LPS: 4246
  # Genes with signif eQTLs plotted in untreated: 4864
  
  ## Place genes in the plot position of the closest variant 
  for(i in 1:length(TWMR_stats_genes_pos$gene_symbol)){
    gene = TWMR_stats_genes_pos$gene_symbol[i]
    gene_chr <- TWMR_stats_genes_pos %>% filter(gene_symbol == gene) %>% select(CHR) %>% as.numeric()
    gene_pos <- TWMR_stats_genes_pos %>% filter(gene_symbol == gene) %>% select(BP) %>% as.numeric()
    chr_variants <- GWAS_stats_subset3 %>% filter(CHR == gene_chr)
    closest_pos_in_plot <- chr_variants[which.min(abs(chr_variants$BP - gene_pos)), "plot_pos"]
    TWMR_stats_genes_pos[i, "plot_pos"] = closest_pos_in_plot
  }
  
  genomewideline_twas = 0.05/dim(TWMR_stats)[1]
  
  ## Bind GWAS and TWMR results
  GWAS_stats_subset3$id <- GWAS_stats_subset3$label
  TWMR_stats_genes_pos$id <- TWMR_stats_genes_pos$gene_symbol
  GWAS_stats_subset3$type = "variant"
  TWMR_stats_genes_pos$type = "gene"
  
  data <- rbind(GWAS_stats_subset3[, c("CHR", "BP", "P", "plot_pos", "id", "type")],
                TWMR_stats_genes_pos[, c("CHR", "BP", "P", "plot_pos", "id", "type")])
  
  ## Add -log10(p)
  data$logP = -log10(data$P)
  ## Label genes above Bonf-corrected signif threshold
  data$label = data$type =="gene" & data$P < genomewideline_twas
  
  ## Params
  chr_colors = rep(c("grey70", "grey20"), 11)
  names(chr_colors) <- paste(1:22)
  data$CHR <- factor(data$CHR, levels = paste(1:22))
  
  ## Position of x-axis CHR labels 
  x = tapply(data$plot_pos, data$CHR, function(pos){floor((max(pos) - min(pos) + 1)/2)}) 
  x = tapply(data$plot_pos, data$CHR, min)  +  x
  
  ## Reverse log-pvals for variants to show them at the bottom plot
  data[data$type == "variant", "logP"] = -data[data$type == "variant", "logP"]
  
  add_y <- ifelse(treatment == "untreated" & TWMR_subset == "all genes", 0.8, 0.2)
  
  miami_plot <- ggplot(data, aes(plot_pos, logP)) +
    geom_point(data = subset(data, label == FALSE), 
               aes(colour = as.factor(CHR)), size = 1, alpha = 0.6) + 
    scale_color_manual(values = chr_colors) + 
    geom_point(data = subset(data, label == TRUE),
               color = "darkred", size = 1, alpha = 0.6) +
    labs(title = paste(treatment,"TWMR (above) vs GWAS (below) p-values"),
         x = "Chromosome", 
         y = expression(-log[10](italic(p)))) +
    scale_y_continuous(expand = c(0, 0), limit = c(-(-log10(genomewideline) + 0.2), (-log10(genomewideline) + add_y)),
                       breaks = seq(from = -floor(-log10(genomewideline)), to = floor(-log10(genomewideline)), by = 1), 
                       labels = c(7:0, 1:7)) +
    scale_x_continuous(breaks = x, labels = 1:22) +
    geom_hline(yintercept = -(-log10(genomewideline) + 0.2), linetype = "solid",col = "grey10", linewidth = 0.7) +
    geom_hline(yintercept = 0.0001,linetype = "solid",col = "grey10", linewidth = 1) +
    geom_hline(yintercept = -log10(genomewideline_twas), color = "darkred", linetype = "dotted") + 
    geom_hline(yintercept =  log10(genomewideline), color = "darkred", linetype = "dotted") + 
    geom_text(data = subset(data, label == TRUE), aes(label = id),
                     size = 3, color = "grey10", vjust = -0.5, hjust = 0.5) +
    theme_bw() +
    theme(legend.position = "none",
          strip.background = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(1,1,3,1), units = "lines"))
  
  return(miami_plot)
  
}

miami_plot_IFN_all_genes <- miami_plotting("IFN", "all genes")
ggsave(miami_plot_IFN_all_genes, filename = paste0(plotdir, "/miami_plot_IFN_all_genes.pdf"), width = 12, height = 6)
miami_plot_LPS_all_genes <- miami_plotting("LPS", "all genes")
ggsave(miami_plot_LPS_all_genes, filename = paste0(plotdir, "/miami_plot_LPS_all_genes.pdf"), width = 12, height = 6)
miami_plot_untreated_all_genes <- miami_plotting("untreated", "all genes")
ggsave(miami_plot_untreated_all_genes, filename = paste0(plotdir, "/miami_plot_untreated_all_genes.pdf"), width = 12, height = 6)

miami_plot_IFN_genes_with_signif_eQTLs <- miami_plotting("IFN", "genes with significant cis-eQTLs")
ggsave(miami_plot_IFN_genes_with_signif_eQTLs, filename = paste0(plotdir, "/miami_plot_IFN_genes_with_signif_eQTLs.pdf"), width = 12, height = 6)
miami_plot_LPS_genes_with_signif_eQTLs <- miami_plotting("LPS", "genes with significant cis-eQTLs")
ggsave(miami_plot_LPS_genes_with_signif_eQTLs, filename = paste0(plotdir, "/miami_plot_LPS_genes_with_signif_eQTLs.pdf"), width = 12, height = 6)
miami_plot_untreated_genes_with_signif_eQTLs <- miami_plotting("untreated", "genes with significant cis-eQTLs")
ggsave(miami_plot_untreated_genes_with_signif_eQTLs, filename = paste0(plotdir, "/miami_plot_untreated_genes_with_signif_eQTLs.pdf"), width = 12, height = 6)



## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#              9.3.3 Explore significant (putative causal) genes*
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# *genes with at least one significant cis-eQTL 

###################################  IFN  ###################################### 
## Nominally significant genes
nominal_signif_genes_IFN <- IFN_results_genes_with_signif_eQTLs %>% filter(P <0.05)
dim(nominal_signif_genes_IFN)
# [1] 122   9

## Signif genes after Bonf correction
bonf_signif_genes_IFN <- IFN_results_genes_with_signif_eQTLs %>% filter(p_adj <0.05) 
#  gene_symbol      alpha         Z         SE            P Nsnps Ngene      p_adj         gene_id
#      NDUFA10 -0.4081067 -4.540307 0.08988527 5.617235e-06     1     1 0.02497423 ENSG00000130414

###################################  LPS  ######################################  
nominal_signif_genes_LPS <- LPS_results_genes_with_signif_eQTLs %>% filter(P <0.05)
dim(nominal_signif_genes_LPS)
# [1] 132   9

bonf_signif_genes_LPS <- LPS_results_genes_with_signif_eQTLs %>% filter(p_adj <0.05) # zero

################################  untreated  ###################################  
nominal_signif_genes_untreated <- untreated_results_genes_with_signif_eQTLs %>% filter(P <0.05)
dim(nominal_signif_genes_untreated)
# [1] 196   9

bonf_signif_genes_untreated <- untreated_results_genes_with_signif_eQTLs %>% filter(p_adj <0.05)
#  gene_symbol      alpha         Z         SE            P Nsnps Ngene       p_adj         gene_id
#        HDHD2 -0.3872739 -4.543805 0.08523117 5.524770e-06     5     2 0.028563062 ENSG00000167220
#         LY86 -0.3465846 -4.926408 0.07035240 8.375496e-07     1     1 0.004330132 ENSG00000112799



## Function for regional plot around genes of interest
ah <- AnnotationHub()
ensDb_v106 <- ah[["AH100643"]] # GRCh38 genome build

regional_plot_around_gene <- function(gene_name, treatment){

  ## GWAS summary stats
  GWAS_stats <- get(load(paste0(input_dir0802, "/summary_stats_microglia_vs_premac_", treatment, ".Rdata")))
  
  ## Lead GWS variants 
  lead_variants <- lapply(get(paste0("gene_tracks_across_chrs_microglia_vs_premac_", treatment)), names) %>% unlist
  
  ## Zoom to variants around gene
  loc = locuszoomr::locus(GWAS_stats, gene = gene_name, 
                          ens_db = ensDb_v106, p = "P", labs = "label", pos = "BP",
                          flank = 2.5e5) 

  ## cis-eQTL effects for the same variants around gene, on that gene
  eQTL_stats_full <- read_csv(paste0(input_dir01, "/nominal_eqtl_results_", treatment, ".csv"))
  eQTL_stats = eQTL_stats_full %>% dplyr::filter(gene_symbol == gene_name & variant_id %in% loc$data$variant_id) %>% 
    dplyr::rename(eqtl_p = pval_nominal)
  
  ## Bind eQTL and GWAS variant stats
  GWAS_eQTL_data <- inner_join(loc$data, eQTL_stats, by = "variant_id") %>% 
    dplyr::rename(gwas_p = P) %>% dplyr::select("variant_id", "gwas_p", "eqtl_p", "BP", "CHR", "label")
  
  GWAS_eQTL_data_melted <- melt(data = GWAS_eQTL_data, id.vars = c("variant_id", "BP", "CHR", "label")) %>% 
    dplyr::rename(P = value, study = variable)
  
  ## Label variants above GWAS suggestive line and GWS lead variants
  genomewideline = 0.05 / (10^6)
  suggestiveline <- genomewideline / 0.005
  GWAS_eQTL_data_melted$gwas_suggestive_signif <- factor(GWAS_eQTL_data_melted$P <= suggestiveline & GWAS_eQTL_data_melted$study == "gwas_p") %>% as.character()
  GWAS_eQTL_data_melted$gwas_lead <- as.character(GWAS_eQTL_data_melted$label %in% lead_variants & GWAS_eQTL_data_melted$study == "gwas_p")
  
  ## Label signif eQTLs 
  eqtl_signif_line = 0.05 / dim(eQTL_stats_full)[1]
  GWAS_eQTL_data_melted$eqtl_signif <- factor(GWAS_eQTL_data_melted$P <= eqtl_signif_line & GWAS_eQTL_data_melted$study == "eqtl_p") %>% as.character()

  ## Reverse eQTL log-pvals to show them at the bottom of x-axis
  GWAS_eQTL_data_melted$logP = -log10(GWAS_eQTL_data_melted$P)
  GWAS_eQTL_data_melted[which(GWAS_eQTL_data_melted$study == "eqtl_p"), "logP"] <-  - GWAS_eQTL_data_melted[which(GWAS_eQTL_data_melted$study == "eqtl_p"), "logP"]

  
  #x x x x x x x x x x x x x x x x   GWAS and eQTL plot   x x x x x x x x x x x x x x x x 
  ## Signif colors 
  colors <- c("FALSE" = "gray", "TRUE"="purple")

  ## Label top most signif GWS SNP 
  top_SNP_gwas <- subset(GWAS_eQTL_data_melted, study == "gwas_p" & gwas_suggestive_signif == TRUE) %>% arrange(P) %>% .[1,]
  ## same as lead variant?
  if(!is.na(top_SNP_gwas$variant_id)){
    top_SNP_gwas$variant_id == GWAS_eQTL_data_melted[which(GWAS_eQTL_data_melted$gwas_lead == TRUE), "variant_id"]
  }
  ## Label top most signif eQTL 
  top_SNP_eqtl <- subset(GWAS_eQTL_data_melted, study == "eqtl_p" & eqtl_signif == TRUE) %>% arrange(P) %>% .[1,]
  
  gwas_eqtl_plot <- ggplot(data = GWAS_eQTL_data_melted, aes(x = BP/1e+06, y = logP)) +
    geom_point(data = subset(GWAS_eQTL_data_melted, study == "gwas_p"), aes(fill = gwas_suggestive_signif),
               shape = 21, color = "grey20", size = 1.65) + 
    geom_point(data = subset(GWAS_eQTL_data_melted, study == "eqtl_p"), aes(fill = eqtl_signif),
               shape = 21, color = "grey20", size = 1.65) + 
    scale_fill_manual(values = colors) + 
    geom_text_repel(data = top_SNP_gwas,
                    aes(label = label),
                    max.overlaps = Inf, min.segment.length = 0,
                    size = 2.5, segment.size = 0.3,
                    nudge_y = 0.45, fontface = "bold") +
    geom_text_repel(data = top_SNP_eqtl, 
                    aes(label = label),
                    max.overlaps = Inf, min.segment.length = 0, 
                    size = 2.5, segment.size = 0.6, nudge_x = -0.05,
                    nudge_y = -0.5, fontface = "bold") +
    geom_hline(yintercept = -log10(suggestiveline), color = "gray70", linetype = "dotted") +
    geom_hline(yintercept = -log10(genomewideline), color = "darkred", linetype = "dotted") +
    geom_hline(yintercept = log10(eqtl_signif_line), color = "darkred", linetype = "dotted") +
    geom_hline(yintercept = 0, linetype = "solid", col = "grey10", linewidth = 1) +
    scale_x_continuous(limit = range(GWAS_eQTL_data_melted$BP/1e+06)) +
    scale_y_continuous(expand = c(0, 0), limit = c((min(GWAS_eQTL_data_melted$logP) - 0.8), (-log10(genomewideline) + 0.4)),
                       breaks = seq(from = floor((min(GWAS_eQTL_data_melted$logP) - 0.8)), to = floor(-log10(genomewideline)), by = 1), 
                       labels = abs(seq(from = floor((min(GWAS_eQTL_data_melted$logP) - 0.8)), to = floor(-log10(genomewideline)), by = 1))) +
    labs(title = paste0("GWAS, eQTL, and TWMR results for ", gene_name, " in ", treatment), 
         y = expression(-log[10](italic(p)))) +
    theme_pubr() +
    theme(legend.position = "none",      
          plot.title = element_text(hjust = (0.5), size =  10),
          axis.title.y = element_text(size = 9),
          axis.text.y = element_text(size = 8),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(), 
          axis.ticks.length.x = unit(0, "cm"))
  
  #x x x x x x x x x x x x x x x x   TWMR plot   x x x x x x x x x x x x x x x x 
  
  ## TWAS stats
  TWMR_stats <- eval(parse_expr(paste0(treatment, "_results_genes_with_signif_eQTLs")))
  gene_p_nominal <- signif(TWMR_stats[which(TWMR_stats$gene_symbol == gene_name), "P"], digits = 2)
  gene_padj <- signif(TWMR_stats[which(TWMR_stats$gene_symbol == gene_name), "p_adj"], digits = 2)
    
  ## Add gene track 
  gene_plot <- gg_genetracks(loc, filter_gene_name = gene_name) +  
    scale_x_continuous(limit = range(GWAS_eQTL_data_melted$BP/1e+06)) +
    geom_text(x = max(GWAS_eQTL_data_melted$BP/1e+06)-0.02, y = 0.77, label = paste0("p: ", gene_p_nominal, "\n",
                                                                                     "padj: ", gene_padj), size = 2.8, fontface = "bold") +
    theme(axis.title.x = element_text(size = 10),
          axis.text = element_text(size = 9))
  
  plot = plot_grid(plotlist = list(gwas_eqtl_plot, gene_plot), ncol = 1, align = "v", rel_heights = c(1, 0.23))

  ggsave(plot, filename = paste0(plotdir, "/GWAS_eQTL_TWMR_results_for_", gene_name, "_in_", treatment, ".pdf"), width = 6, height = 5)
  
}

## Regional plots for signif TWMR genes after Bonf correction
regional_plot_around_gene("NDUFA10", "IFN")
regional_plot_around_gene("LY86", "untreated")
regional_plot_around_gene("HDHD2", "untreated")


# ______________________________________________________________________________
#                9.3.3.1 Overlaps with proliferation DEGs* 
# ______________________________________________________________________________
# *for proliferation from preMac -> microglia 

load(paste0(input_dir0701, "/de_genes_proliferation_IFN.Rdata"), verbose = T)
load(paste0(input_dir0701, "/de_genes_proliferation_LPS.Rdata"), verbose = T)
load(paste0(input_dir0701, "/de_genes_proliferation_untreated.Rdata"), verbose = T)

de_genes_proliferation_IFN <- de_genes_proliferation_IFN$gene_id
de_genes_proliferation_LPS <- de_genes_proliferation_LPS$gene_id
de_genes_proliferation_untreated <- de_genes_proliferation_untreated$gene_id

#~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ 
#                           Nominal causal genes in IFN
#~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~
paste0("Of ", length(de_genes_proliferation_IFN), " DEGs in IFN and ", dim(nominal_signif_genes_IFN)[1], 
       " TWMR nominally signif genes in IFN, ", length(intersect(nominal_signif_genes_IFN$gene_id, de_genes_proliferation_IFN)),
       " are both.")
# [1] "Of 6579 DEGs in IFN and 122 TWMR nominally signif genes in IFN, 60 are both."
paste0("Of ", length(de_genes_proliferation_LPS), " DEGs in LPS and ", dim(nominal_signif_genes_IFN)[1], 
       " TWMR nominally signif genes in IFN, ", length(intersect(nominal_signif_genes_IFN$gene_id, de_genes_proliferation_LPS)),
       " are both.")
# [1] "Of 3035 DEGs in LPS and 122 TWMR nominally signif genes in IFN, 34 are both."
paste0("Of ", length(de_genes_proliferation_untreated), " DEGs in untreated and ", dim(nominal_signif_genes_IFN)[1], 
       " TWMR nominally signif genes in IFN, ", length(intersect(nominal_signif_genes_IFN$gene_id, de_genes_proliferation_untreated)),
       " are both.")
# [1] "Of 7454 DEGs in untreated and 122 TWMR nominally signif genes in IFN, 70 are both."

#~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ 
#                            Nominal causal genes in LPS
#~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ 
paste0("Of ", length(de_genes_proliferation_IFN), " DEGs in IFN and ", dim(nominal_signif_genes_LPS)[1], 
       " TWMR nominally signif genes in LPS, ", length(intersect(nominal_signif_genes_LPS$gene_id, de_genes_proliferation_IFN)),
       " are both.")
# [1] "Of 6579 DEGs in IFN and 132 TWMR nominally signif genes in LPS, 58 are both."
paste0("Of ", length(de_genes_proliferation_LPS), " DEGs in LPS and ", dim(nominal_signif_genes_LPS)[1], 
       " TWMR nominally signif genes in LPS, ", length(intersect(nominal_signif_genes_LPS$gene_id, de_genes_proliferation_LPS)),
       " are both.")
# [1] "Of 3035 DEGs in LPS and 132 TWMR nominally signif genes in LPS, 28 are both."
paste0("Of ", length(de_genes_proliferation_untreated), " DEGs in untreated and ", dim(nominal_signif_genes_LPS)[1], 
       " TWMR nominally signif genes in LPS, ", length(intersect(nominal_signif_genes_LPS$gene_id, de_genes_proliferation_untreated)),
       " are both.")
# [1] "Of 7454 DEGs in untreated and 132 TWMR nominally signif genes in LPS, 63 are both."

#~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ 
#                       Nominal causal genes in untreated
#~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ 
paste0("Of ", length(de_genes_proliferation_IFN), " DEGs in IFN and ", dim(nominal_signif_genes_untreated)[1], 
       " TWMR nominally signif genes in untreated, ", length(intersect(nominal_signif_genes_untreated$gene_id, de_genes_proliferation_IFN)),
       " are both.")
# [1] "Of 6579 DEGs in IFN and 196 TWMR nominally signif genes in untreated, 96 are both."
paste0("Of ", length(de_genes_proliferation_LPS), " DEGs in LPS and ", dim(nominal_signif_genes_untreated)[1], 
       " TWMR nominally signif genes in untreated, ", length(intersect(nominal_signif_genes_untreated$gene_id, de_genes_proliferation_LPS)),
       " are both.")
# [1] "Of 3035 DEGs in LPS and 196 TWMR nominally signif genes in untreated, 47 are both."
paste0("Of ", length(de_genes_proliferation_untreated), " DEGs in untreated and ", dim(nominal_signif_genes_untreated)[1], 
       " TWMR nominally signif genes in untreated, ", length(intersect(nominal_signif_genes_untreated$gene_id, de_genes_proliferation_untreated)),
       " are both.")
# [1] "Of 7454 DEGs in untreated and 196 TWMR nominally signif genes in untreated, 99 are both."


# ______________________________________________________________________________
#   9.3.3.2 Overlaps with genes with proliferation-related deleterious burdens
# ______________________________________________________________________________

del_burden_results <- as.data.frame(read_csv("/lustre/scratch123/hgi/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/data/results/2.2.rare_vars_vs_prolif/del_burden_scaled_centered_prop_pvals.csv"))
deleterious_Burden_results_signif <- subset(del_burden_results, p_Bonf < 0.05)

## Divide by prolif comparison
deleterious_Burden_signif_premac_vs_iPSC <- subset(deleterious_Burden_results_signif, comparison == "line_prop_changes_premac_iPSC")$gene_name
deleterious_Burden_signif_old_vs_young_premac <- subset(deleterious_Burden_results_signif, comparison == "line_prop_changes_old_vs_young_premac")$gene_name
deleterious_Burden_signif_microglia_vs_premac <- subset(deleterious_Burden_results_signif, comparison == "line_prop_changes_microglia_premac")$gene_name

## Add ensemble IDs
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
deleterious_Burden_signif_premac_vs_iPSC_ids <- getBM(values = deleterious_Burden_signif_premac_vs_iPSC,
                                                      mart = ensembl,
                                                      attributes = c("ensembl_gene_id"),
                                                      filters = "external_gene_name") %>% unname() %>% unlist()

deleterious_Burden_signif_old_vs_young_premac_ids <- getBM(values = deleterious_Burden_signif_old_vs_young_premac,
                                                           mart = ensembl,
                                                           attributes = c("ensembl_gene_id"),
                                                           filters = "external_gene_name") %>% unname() %>% unlist()

deleterious_Burden_signif_microglia_vs_premac_ids <- getBM(values = deleterious_Burden_signif_microglia_vs_premac,
                                                           mart = ensembl,
                                                           attributes = c("ensembl_gene_id"),
                                                           filters = "external_gene_name") %>% unname() %>% unlist()

#~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ 
#                           Nominal causal genes in IFN
#~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~
intersect(nominal_signif_genes_IFN$gene_id, deleterious_Burden_signif_premac_vs_iPSC_ids) # zero
intersect(nominal_signif_genes_IFN$gene_id, deleterious_Burden_signif_old_vs_young_premac_ids) # zero
intersect(nominal_signif_genes_IFN$gene_id, deleterious_Burden_signif_microglia_vs_premac_ids) # zero
#~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ 
#                           Nominal causal genes in LPS
#~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~
intersect(nominal_signif_genes_LPS$gene_id, deleterious_Burden_signif_premac_vs_iPSC_ids) # zero
intersect(nominal_signif_genes_LPS$gene_id, deleterious_Burden_signif_old_vs_young_premac_ids) # zero
intersect(nominal_signif_genes_LPS$gene_id, deleterious_Burden_signif_microglia_vs_premac_ids) 
# [1] "ENSG00000081014" (AP4E1)
## Regional plot for overlapping gene
regional_plot_around_gene("AP4E1", "LPS")
#~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ 
#                           Nominal causal genes in untreated
#~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~
intersect(nominal_signif_genes_untreated$gene_id, deleterious_Burden_signif_premac_vs_iPSC_ids) # zero
intersect(nominal_signif_genes_untreated$gene_id, deleterious_Burden_signif_old_vs_young_premac_ids) # zero
intersect(nominal_signif_genes_untreated$gene_id, deleterious_Burden_signif_microglia_vs_premac_ids) # zero


# ______________________________________________________________________________
#            9.3.3.3 Overlaps with nearest genes to GWAS lead variants
# ______________________________________________________________________________

## Closest genes: 
load(paste0(input_dir0802, "/closest_genes_leadSNPs_microglia_vs_premac_IFN.Rdata"), verbose = T)
#   gene_tracks_across_chrs_microglia_vs_premac_IFN
closest_genes_microglia_vs_premac_IFN <- lapply(gene_tracks_across_chrs_microglia_vs_premac_IFN, 
                                                function(chr)lapply(chr, function(leadSNP){leadSNP[,"gene_id"]})) %>% unlist %>% unique
length(closest_genes_microglia_vs_premac_IFN)
# 106

load(paste0(input_dir0802, "/closest_genes_leadSNPs_microglia_vs_premac_LPS.Rdata"), verbose = T)
#   gene_tracks_across_chrs_microglia_vs_premac_LPS
closest_genes_microglia_vs_premac_LPS <- lapply(gene_tracks_across_chrs_microglia_vs_premac_LPS, 
                                                function(chr)lapply(chr, function(leadSNP){leadSNP[,"gene_id"]})) %>% unlist %>% unique
length(closest_genes_microglia_vs_premac_LPS)
# 123

load(paste0(input_dir0802, "/closest_genes_leadSNPs_microglia_vs_premac_untreated.Rdata"), verbose = T)
#   gene_tracks_across_chrs_microglia_vs_premac_untreated
closest_genes_microglia_vs_premac_untreated <- lapply(gene_tracks_across_chrs_microglia_vs_premac_untreated, 
                                                      function(chr)lapply(chr, function(leadSNP){leadSNP[,"gene_id"]})) %>% unlist %>% unique
length(closest_genes_microglia_vs_premac_untreated)
# 51

#~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ 
#                           Nominal causal genes in IFN
#~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~
## Overlap with closest genes to lead SNPs in IFN:
closest_genes_IFN_nominal_signif_IFN <- intersect(closest_genes_microglia_vs_premac_IFN, nominal_signif_genes_IFN$gene_id) 
subset(nominal_signif_genes_IFN, gene_id %in% closest_genes_IFN_nominal_signif_IFN)
#   gene_symbol     alpha        Z        SE          P Nsnps Ngene p_adj         gene_id
#          CYCS 0.3161764 1.965718 0.1608452 0.04933116     1     1     1 ENSG00000172115
## Regional plot 
regional_plot_around_gene("CYCS", "IFN")

## Overlap with closest genes to lead SNPs in LPS:
intersect(closest_genes_microglia_vs_premac_LPS, nominal_signif_genes_IFN$gene_id) # zero

## Overlap with closest genes to lead SNPs in untreated:
intersect(closest_genes_microglia_vs_premac_untreated, nominal_signif_genes_IFN$gene_id) # zero

#~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ 
#                           Nominal causal genes in LPS
#~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~
## Overlap with closest genes to lead SNPs in IFN:
intersect(closest_genes_microglia_vs_premac_IFN, nominal_signif_genes_LPS$gene_id) # zero

## Overlap with closest genes to lead SNPs in LPS:
closest_genes_LPS_nominal_signif_LPS <- intersect(closest_genes_microglia_vs_premac_LPS, nominal_signif_genes_LPS$gene_id) 
subset(nominal_signif_genes_LPS, gene_id %in% closest_genes_LPS_nominal_signif_LPS)
# gene_symbol     alpha        Z        SE           P Nsnps Ngene p_adj         gene_id
#        ELP3 -1.035331 -2.93437 0.3528291 0.003342259     1     1     1 ENSG00000134014
regional_plot_around_gene("ELP3", "LPS")

## Overlap with closest genes to lead SNPs in untreated:
intersect(closest_genes_microglia_vs_premac_untreated, nominal_signif_genes_LPS$gene_id) # zero

#~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ 
#                       Nominal causal genes in untreated
#~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~ ~.~
## Overlap with closest genes to lead SNPs in IFN:
closest_genes_IFN_nominal_signif_untreated <- intersect(closest_genes_microglia_vs_premac_IFN, nominal_signif_genes_untreated$gene_id)
subset(nominal_signif_genes_untreated, gene_id %in% closest_genes_IFN_nominal_signif_untreated)
#  gene_symbol      alpha         Z        SE          P Nsnps Ngene p_adj         gene_id
#       OCIAD1 -0.5893622 -2.095749 0.2812179 0.03610446     1     1     1 ENSG00000109180

## Overlap with closest genes to lead SNPs in LPS:
closest_genes_LPS_nominal_signif_untreated <- intersect(closest_genes_microglia_vs_premac_LPS, nominal_signif_genes_untreated$gene_id)
subset(nominal_signif_genes_untreated, gene_id %in% closest_genes_LPS_nominal_signif_untreated)
# gene_symbol      alpha         Z        SE          P Nsnps Ngene p_adj         gene_id
#      OCIAD1 -0.5893622 -2.095749 0.2812179 0.03610446     1     1     1 ENSG00000109180
#       RAP1A  0.5152072  2.259997 0.2279681 0.02382143     4     2     1 ENSG00000116473

## Overlap with closest genes to lead SNPs in untreated:
intersect(closest_genes_microglia_vs_premac_untreated, nominal_signif_genes_untreated$gene_id) # zero


# ______________________________________________________________________________
#        9.3.3.4 Comparison of putative causal genes across treatments 
# ______________________________________________________________________________

nominal_signif_genes_IFN_LPS <- intersect(nominal_signif_genes_IFN$gene_id, nominal_signif_genes_LPS$gene_id) # 18
nominal_signif_genes_IFN_untreated <- intersect(nominal_signif_genes_IFN$gene_id, nominal_signif_genes_untreated$gene_id) # 14
nominal_signif_genes_LPS_untreated <- intersect(nominal_signif_genes_LPS$gene_id, nominal_signif_genes_untreated$gene_id) # 20
nominal_signif_genes_IFN_LPS_untreated <- intersect(nominal_signif_genes_IFN_LPS, nominal_signif_genes_LPS_untreated) # 7







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
# date     2025-04-01
# rstudio  2024.04.0+735 Chocolate Cosmos (server)
# pandoc   3.1.12.3 @ /opt/view/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package                * version     date (UTC) lib source
# abind                    1.4-5       2016-07-21 [1] CRAN (R 4.3.1)
# AnnotationDbi          * 1.64.1      2023-11-03 [1] Bioconductor
# AnnotationFilter       * 1.26.0      2023-10-24 [1] Bioconductor
# AnnotationHub          * 3.10.0      2023-10-24 [1] Bioconductor
# backports                1.4.1       2021-12-13 [1] CRAN (R 4.3.1)
# Biobase                * 2.62.0      2023-10-24 [1] Bioconductor
# BiocFileCache          * 2.10.2      2024-03-27 [1] Bioconductor 3.18 (R 4.3.1)
# BiocGenerics           * 0.48.1      2023-11-01 [1] Bioconductor
# BiocIO                   1.12.0      2023-10-24 [1] Bioconductor
# BiocManager              1.30.22     2023-08-08 [1] CRAN (R 4.3.1)
# BiocParallel             1.36.0      2023-10-24 [1] Bioconductor
# BiocVersion              3.18.1      2023-11-15 [1] Bioconductor
# biomaRt                * 2.58.2      2024-01-30 [1] Bioconductor 3.18 (R 4.3.1)
# Biostrings               2.70.3      2024-03-13 [1] Bioconductor 3.18 (R 4.3.1)
# bit                      4.0.5       2022-11-15 [1] CRAN (R 4.3.1)
# bit64                    4.0.5       2020-08-30 [1] CRAN (R 4.3.1)
# bitops                   1.0-7       2021-04-24 [1] CRAN (R 4.3.1)
# blob                     1.2.4       2023-03-17 [1] CRAN (R 4.3.1)
# broom                    1.0.5       2023-06-09 [1] CRAN (R 4.3.1)
# cachem                   1.0.8       2023-05-01 [1] CRAN (R 4.3.1)
# car                      3.1-2       2023-03-30 [1] CRAN (R 4.3.1)
# carData                  3.0-5       2022-01-06 [1] CRAN (R 4.3.1)
# cli                      3.6.2       2023-12-11 [1] CRAN (R 4.3.1)
# codetools                0.2-20      2024-03-31 [1] CRAN (R 4.3.1)
# colorspace               2.1-0       2023-01-23 [1] CRAN (R 4.3.1)
# cowplot                * 1.1.3       2024-01-22 [1] CRAN (R 4.3.1)
# crayon                   1.5.2       2022-09-29 [1] CRAN (R 4.3.1)
# curl                     5.2.1       2024-03-01 [1] CRAN (R 4.3.1)
# data.table               1.15.4      2024-03-30 [1] CRAN (R 4.3.1)
# DBI                      1.2.2       2024-02-16 [1] CRAN (R 4.3.1)
# dbplyr                 * 2.5.0       2024-03-19 [1] CRAN (R 4.3.1)
# DelayedArray             0.28.0      2023-10-24 [1] Bioconductor
# digest                   0.6.35      2024-03-11 [1] CRAN (R 4.3.1)
# dplyr                  * 1.1.4       2023-11-17 [1] CRAN (R 4.3.1)
# edgeR                    4.0.16      2024-02-18 [1] Bioconductor 3.18 (R 4.3.1)
# ensembldb              * 2.26.0      2023-10-24 [1] Bioconductor
# fansi                    1.0.6       2023-12-08 [1] CRAN (R 4.3.1)
# farver                   2.1.1       2022-07-06 [1] CRAN (R 4.3.1)
# fastmap                  1.1.1       2023-02-24 [1] CRAN (R 4.3.1)
# filelock                 1.0.3       2023-12-11 [1] CRAN (R 4.3.1)
# forcats                * 1.0.0       2023-01-29 [1] CRAN (R 4.3.1)
# generics                 0.1.3       2022-07-05 [1] CRAN (R 4.3.1)
# GenomeInfoDb           * 1.38.8      2024-03-15 [1] Bioconductor 3.18 (R 4.3.1)
# GenomeInfoDbData         1.2.11      2025-02-07 [1] Bioconductor
# GenomicAlignments        1.38.2      2024-01-16 [1] Bioconductor 3.18 (R 4.3.1)
# GenomicFeatures        * 1.54.4      2024-03-13 [1] Bioconductor 3.18 (R 4.3.1)
# GenomicRanges          * 1.54.1      2023-10-29 [1] Bioconductor
# gggrid                   0.2-0       2022-01-11 [1] CRAN (R 4.3.1)
# ggplot2                * 3.5.1       2024-04-23 [1] CRAN (R 4.3.1)
# ggpubr                 * 0.6.0       2023-02-10 [1] CRAN (R 4.3.1)
# ggrepel                * 0.9.5       2024-01-10 [1] CRAN (R 4.3.1)
# ggsignif                 0.6.4       2022-10-13 [1] CRAN (R 4.3.1)
# glue                     1.7.0       2024-01-09 [1] CRAN (R 4.3.1)
# gtable                   0.3.4       2023-08-21 [1] CRAN (R 4.3.1)
# hms                      1.1.3       2023-03-21 [1] CRAN (R 4.3.1)
# htmltools                0.5.8       2024-03-25 [1] CRAN (R 4.3.1)
# htmlwidgets              1.6.4       2023-12-06 [1] CRAN (R 4.3.1)
# httpuv                   1.6.15      2024-03-26 [1] CRAN (R 4.3.1)
# httr                   * 1.4.7       2023-08-15 [1] CRAN (R 4.3.1)
# interactiveDisplayBase   1.40.0      2023-10-24 [1] Bioconductor
# IRanges                * 2.36.0      2023-10-24 [1] Bioconductor
# jsonlite                 1.8.8       2023-12-04 [1] CRAN (R 4.3.1)
# KEGGREST                 1.42.0      2023-10-24 [1] Bioconductor
# labeling                 0.4.3       2023-08-29 [1] CRAN (R 4.3.1)
# later                    1.3.2       2023-12-06 [1] CRAN (R 4.3.1)
# lattice                  0.22-6      2024-03-20 [1] CRAN (R 4.3.1)
# lazyeval                 0.2.2       2019-03-15 [1] CRAN (R 4.3.1)
# LDlinkR                  1.4.0       2024-04-10 [1] CRAN (R 4.3.1)
# lifecycle                1.0.4       2023-11-07 [1] CRAN (R 4.3.1)
# limma                    3.58.1      2023-10-31 [1] Bioconductor
# locfit                   1.5-9.9     2024-03-01 [1] CRAN (R 4.3.1)
# locuszoomr             * 0.2.1       2024-02-17 [1] CRAN (R 4.3.1)
# lubridate              * 1.9.3       2023-09-27 [1] CRAN (R 4.3.1)
# magrittr                 2.0.3       2022-03-30 [1] CRAN (R 4.3.1)
# Matrix                   1.6-5       2024-01-11 [1] CRAN (R 4.3.1)
# MatrixGenerics           1.14.0      2023-10-24 [1] Bioconductor
# matrixStats              1.2.0       2023-12-11 [1] CRAN (R 4.3.1)
# memoise                  2.0.1       2021-11-26 [1] CRAN (R 4.3.1)
# mime                     0.12        2021-09-28 [1] CRAN (R 4.3.1)
# munsell                  0.5.1       2024-04-01 [1] CRAN (R 4.3.1)
# patchwork                1.2.0       2024-01-08 [1] CRAN (R 4.3.1)
# pillar                   1.9.0       2023-03-22 [1] CRAN (R 4.3.1)
# pkgconfig                2.0.3       2019-09-22 [1] CRAN (R 4.3.1)
# pkgload                  1.3.4       2024-01-16 [1] CRAN (R 4.3.1)
# plotly                   4.10.4      2024-01-13 [1] CRAN (R 4.3.1)
# plyr                     1.8.9       2023-10-02 [1] CRAN (R 4.3.1)
# png                      0.1-8       2022-11-29 [1] CRAN (R 4.3.1)
# prettyunits              1.2.0       2023-09-24 [1] CRAN (R 4.3.1)
# progress                 1.2.3       2023-12-06 [1] CRAN (R 4.3.1)
# promises                 1.2.1       2023-08-10 [1] CRAN (R 4.3.1)
# ProtGenerics             1.34.0      2023-10-24 [1] Bioconductor
# purrr                  * 1.0.2       2023-08-10 [1] CRAN (R 4.3.1)
# R6                       2.5.1       2021-08-19 [1] CRAN (R 4.3.1)
# ragg                     1.3.0       2024-03-13 [1] CRAN (R 4.3.1)
# rappdirs                 0.3.3       2021-01-31 [1] CRAN (R 4.3.1)
# Rcpp                     1.0.12      2024-01-09 [1] CRAN (R 4.3.1)
# RCurl                    1.98-1.14   2024-01-09 [1] CRAN (R 4.3.1)
# readr                  * 2.1.5       2024-01-10 [1] CRAN (R 4.3.1)
# reshape2               * 1.4.4       2020-04-09 [1] CRAN (R 4.3.1)
# restfulr                 0.0.15      2022-06-16 [1] CRAN (R 4.3.1)
# rjson                    0.2.21      2022-01-09 [1] CRAN (R 4.3.1)
# rlang                  * 1.1.3       2024-01-10 [1] CRAN (R 4.3.1)
# Rsamtools                2.18.0      2023-10-24 [1] Bioconductor
# RSQLite                  2.3.6       2024-03-31 [1] CRAN (R 4.3.1)
# rstatix                  0.7.2       2023-02-01 [1] CRAN (R 4.3.1)
# rstudioapi               0.16.0      2024-03-24 [1] CRAN (R 4.3.1)
# rtracklayer              1.62.0      2023-10-24 [1] Bioconductor
# S4Arrays                 1.2.1       2024-03-04 [1] Bioconductor 3.18 (R 4.3.1)
# S4Vectors              * 0.40.2      2023-11-23 [1] Bioconductor 3.18 (R 4.3.1)
# scales                   1.3.0       2023-11-28 [1] CRAN (R 4.3.1)
# sessioninfo            * 1.2.2       2021-12-06 [1] CRAN (R 4.3.1)
# shiny                    1.8.1.1     2024-04-02 [1] CRAN (R 4.3.1)
# SparseArray              1.2.4       2024-02-11 [1] Bioconductor 3.18 (R 4.3.1)
# statmod                  1.5.0       2023-01-06 [1] CRAN (R 4.3.1)
# stringi                  1.8.3       2023-12-11 [1] CRAN (R 4.3.1)
# stringr                * 1.5.1       2023-11-14 [1] CRAN (R 4.3.1)
# SummarizedExperiment     1.32.0      2023-10-24 [1] Bioconductor
# systemfonts              1.0.6       2024-03-07 [1] CRAN (R 4.3.1)
# textshaping              0.3.7       2023-10-09 [1] CRAN (R 4.3.1)
# tibble                 * 3.2.1       2023-03-20 [1] CRAN (R 4.3.1)
# tidyr                  * 1.3.1       2024-01-24 [1] CRAN (R 4.3.1)
# tidyselect               1.2.1       2024-03-11 [1] CRAN (R 4.3.1)
# tidyverse              * 2.0.0       2023-02-22 [1] CRAN (R 4.3.1)
# timechange               0.3.0       2024-01-18 [1] CRAN (R 4.3.1)
# tzdb                     0.4.0       2023-05-12 [1] CRAN (R 4.3.1)
# utf8                     1.2.4       2023-10-22 [1] CRAN (R 4.3.1)
# vctrs                    0.6.5       2023-12-01 [1] CRAN (R 4.3.1)
# viridisLite              0.4.2       2023-05-02 [1] CRAN (R 4.3.1)
# vroom                    1.6.5       2023-12-05 [1] CRAN (R 4.3.1)
# withr                    3.0.0       2024-01-16 [1] CRAN (R 4.3.1)
# XML                      3.99-0.16.1 2024-01-22 [1] CRAN (R 4.3.1)
# xml2                     1.3.6       2023-12-04 [1] CRAN (R 4.3.1)
# xtable                   1.8-4       2019-04-21 [1] CRAN (R 4.3.1)
# XVector                  0.42.0      2023-10-24 [1] Bioconductor
# yaml                     2.3.8       2023-12-11 [1] CRAN (R 4.3.1)
# zlibbioc                 1.48.2      2024-03-13 [1] Bioconductor 3.18 (R 4.3.1)
# zoo                      1.8-12      2023-04-13 [1] CRAN (R 4.3.1)
# 
# [1] /opt/view/rlib/R/library
# [2] /software/hgi/installs/softpack/rstudio/.spack/linux-ubuntu22.04-x86_64_v3/gcc-11.4.0/r-4.3.1-bfwldrk76z6f52upk47zepliekn7ayqz/rlib/R/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────






