
library(tidyverse)
library(dplyr)
library(rlang)
library(cowplot)
library(ggrepel)
library(AnnotationHub)
library(locuszoomr)
library(sessioninfo)

#-------------------------------------------------------------------------------
#                8.2 Checking results of proliferation GWASes 
#-------------------------------------------------------------------------------
#  Code to explore and plot results from each proliferation GWAS across all 
#  variant chunks.  
#  * Note: code based on Marta's code and analysis pipeline.
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

## Working dir
setwd("/lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_Daianna/")

## Output data and plot dirs 
outdir = paste(getwd(), "output_data", "08_GWAS_proliferation", "02_Check_association_results", sep = "/")
input_dir = paste(getwd(), "input_data", "08_GWAS_proliferation", "02_Check_association_results", sep = "/")
plotdir = paste(getwd(), "plots", "08_GWAS_proliferation", "02_Check_association_results", sep = "/")
dir.create(outdir, recursive = T)
dir.create(input_dir, recursive = T)
dir.create(plotdir, recursive = T)

## Subdirs for locus zoom plots
phenotypes <- c("premac_vs_iPSC", 
                "old_vs_young_premac", 
                "microglia_vs_premac_IFN", 
                "microglia_vs_premac_LPS", 
                "microglia_vs_premac_untreated")
plot_subdirs = paste(plotdir, "regional_plots", phenotypes, sep = "/")
for(subdir in plot_subdirs){dir.create(subdir, recursive = T)}

## Input dir
input_dir_01 <- paste(getwd(), "output_data", "08_GWAS_proliferation", "01_Association_testing", sep = "/")


## Create summary stats files
create_summary_stats <- function(phenotype){
  
  ## Sample sizes
  num_samples <- c("premac_vs_iPSC" = 228, 
                   "old_vs_young_premac" = 210,
                   "microglia_vs_premac_IFN" = 230,
                   "microglia_vs_premac_LPS" = 229,
                   "microglia_vs_premac_untreated" = 209)
  
  ## Bind results across all variant chunks
  summary_stats <- vector()
  
  for(chunk in 1:582){
    load(paste0(input_dir_01, "/", phenotype, "_results/GWAS_", phenotype, "_results_", chunk, ".Rdata"), verbose = T)
    # Loading objects:
    #   results
    summary_stats <-  rbind(summary_stats,
                            lapply(results, function(snp){if(length(snp$residuals) == num_samples[phenotype]){ 
                              snp$coefficients["allele_dosage",] } else{print("Error")}}) %>% do.call(rbind, .))
    
  }
  
  ## Add variant ID and sample size
  summary_stats <- as.data.frame(summary_stats)
  summary_stats$variant_id <- rownames(summary_stats)
  
  ## Confirm number of variants
  dim(summary_stats)
  # iPSC to preMac: 5811307       5
  # young to old preMac: 5811307       5
  # preMac to microglia in IFN: 5811307      5
  # preMac to microglia in LPS: 5811307      5
  # preMac to microglia in untreated: 5811307     5

  length(unique(summary_stats$variant_id))
  # iPSC to preMac: 5811307
  # young to old preMac: 5811307
  # preMac to microglia in IFN: 5811307
  # preMac to microglia in LPS: 5811307
  # preMac to microglia in untreated: 5811307
  
  ## No NAs
  which(is.na(summary_stats))
  # integer(0)
  
  ## Add sample size
  summary_stats$N <- num_samples[phenotype]
  
  ## Bonferroni correction = p * num tests
  summary_stats$p_bonf <- sapply(summary_stats$`Pr(>|t|)` * dim(summary_stats)[1], function(x){min(1, x)})
  ## Same with p.adjust()
  # summary_stats$p_bonf2 <-  p.adjust(summary_stats$`Pr(>|t|)`, method = "bonferroni")
  
  assign(paste0("summary_stats_", phenotype), summary_stats)
  
}

for (phenotype in phenotypes){
  create_summary_stats(phenotype)  
}


## Create input variant file to search their rsIDs  
summary_stats_premac_vs_iPSC$CHR <- str_split_fixed(summary_stats_premac_vs_iPSC$variant_id, "_", 4)[,1] 
summary_stats_premac_vs_iPSC$BP  <- str_split_fixed(summary_stats_premac_vs_iPSC$variant_id, "_", 4)[,2] 

write.table(summary_stats_premac_vs_iPSC[, c("CHR", "BP")], file = paste0(outdir, "/variants_for_GWAS.txt"), col.names = F, row.names = F, quote = F, sep = "\t")

## Read output from find_rsIDs.sh
variant_IDs = read_table(paste0(input_dir, "/dbsnp.156.variants_for_GWAS.subset.tsv")) %>% as.data.frame()
variant_IDs = rbind(colnames(variant_IDs), variant_IDs)
colnames(variant_IDs) <- c("CHR", "BP", "rsID", "Ref", "Alt")
variant_IDs$SNP = paste(variant_IDs$CHR, variant_IDs$BP, variant_IDs$Ref, variant_IDs$Alt, sep = "_")
save(variant_IDs, file = paste0(outdir, "/variants_for_GWAS_with_rsIDs.Rdata"))


## Add rsIDs to variants (CHR-BP-REF-ALT if not)
for(phenotype in phenotypes){
  
  summary_stats <- eval(parse_expr(paste0("summary_stats_", phenotype)))
  
  summary_stats$rsID <- variant_IDs[match(summary_stats$variant_id, variant_IDs$SNP, nomatch = NA), "rsID"]
  summary_stats = summary_stats %>% mutate(label = if_else(is.na(rsID), variant_id, rsID))
  
  summary_stats$CHR <- str_split_fixed(summary_stats$variant_id, "_", 4)[,1] %>% as.numeric()
  summary_stats$BP <- str_split_fixed(summary_stats$variant_id, "_", 4)[,2] %>% as.numeric()
  summary_stats$P <- summary_stats$`Pr(>|t|)`
  
  assign(paste0("summary_stats_", phenotype), summary_stats)
  save(summary_stats, file = paste0(outdir, "/summary_stats_", phenotype, ".Rdata"))
  write_csv(x = summary_stats, file = paste0(outdir, "/summary_stats_", phenotype, ".csv"), col_names = T)
  
}


## Q-Q plots
Q_Q_plot <- function(phenotype){
   
   summary_stats <- eval(parse_expr(paste0("summary_stats_", phenotype)))
   
   ## Discard variants in LD
   not_LD = read_tsv("../OTAR2065_sc_eQTL/data/genotype/plink_genotypes/all_pools.all_donors.genotype.MAF05.LDpruned.prune.in",
                     col_names = FALSE) %>%
     dplyr::mutate(variant_id = str_remove(X1,"chr"))
   
   summary_stats_notLD = summary_stats[which(summary_stats$variant_id %in% not_LD$variant_id), ]
   dim(summary_stats_notLD)
   # [1] 97404     7
   
   ## Generate expected p-vals as (1:n) - a / (n + 1 - 2a) with a = 0.375
   exp_pvals <- ppoints(n = length(summary_stats_notLD$`Pr(>|t|)`))
   ## Exp pvals are uniformly distibuted
   hist(exp_pvals)

   ## Obs pvals
   obs_pvals <- summary_stats_notLD$`Pr(>|t|)`
   
   ## Q-Q plot
   df = data.frame("expected_logp"= sort(-log10(exp_pvals)), 
                   "observed_logp" = sort(-log10(obs_pvals)))
   
   phenotypes <- c("premac_vs_iPSC", "old_vs_young_premac", 
                   "microglia_vs_premac_IFN", 
                   "microglia_vs_premac_LPS", 
                   "microglia_vs_premac_untreated")
   titles <- c("iPSC to preMac", "Young to old preMac", 
               "preMac to microglia in IFN", 
               "preMac to microglia in LPS",
               "preMac to microglia in untreated")
   names(titles) <- phenotypes
   
   ## Deflated pvals in -log10 scale
   p <- ggplot(data = df,
           aes(x = expected_logp, y = observed_logp)) +
           geom_point(color = "gray", size = 0.6, alpha = 0.5) + 
           geom_line(aes(x = expected_logp, y = expected_logp), color = "red", linewidth = 0.4) +
           labs(x = "Expected -log10(P)", y = "Observed -log10(P)", title = titles[phenotype]) +
           theme_bw() +
           theme(plot.title = element_text(face="bold"))
   
    return(p)
}


## Manhattan plots
Manhattan_plot <- function(phenotype){
  
  summary_stats <- eval(parse_expr(paste0("summary_stats_", phenotype)))
  
  ## Add row num for filtering for plotting
  summary_stats = summary_stats %>% 
    rename(SNP = variant_id) %>% 
    arrange(P) %>% 
    mutate(row_num = 1:nrow(summary_stats))

  ## Keep all variants with P<=0.01
  ## Hide half variants with P>0.01
  summary_stats$keep_row = apply(summary_stats, 1, 
                                 function(x){if(as.numeric(x["P"]) > 0.01 & as.numeric(x["row_num"]) %% 2 == 1){FALSE} else{TRUE}})
  summary_stats_subset <- summary_stats %>% filter(keep_row == TRUE) %>% mutate(row_num = 1:nrow(.))
  
  ## Keep all variants with P<=0.05
  ## Only keep 10% of variants with P>0.05
  summary_stats_subset$keep_row2 = apply(summary_stats_subset, 1, 
                                  function(x){if(as.numeric(x["P"]) > 0.05 & as.numeric(x["row_num"]) %% 10 != 0){FALSE} else{TRUE}})
  summary_stats_subset2 <- summary_stats_subset %>% filter(keep_row2 == TRUE) %>% mutate(row_num = 1:nrow(.))

  ## Keep all variants with P<=0.3
  ## Only keep 0.3% of variants with P>0.3
  summary_stats_subset2$keep_row3 = apply(summary_stats_subset2, 1, 
                                         function(x){if(as.numeric(x["P"]) > 0.3 & as.numeric(x["row_num"]) %% 30 != 0){FALSE} else{TRUE}})
  summary_stats_subset3 <- summary_stats_subset2 %>% filter(keep_row3 == TRUE) 
  dim(summary_stats_subset3)
  # iPSC -> preMac: 273974     16
  # young -> old preMac: 279760     16
  # preMac -> microglia IFN: 276862     14
  # preMac -> microglia LPS: 273595    15
  # preMac -> microglia untreated: 276736     14
  
  pdf(paste0(plotdir, "/Raw_pvals_filtered_variants_for_Manhattan_", phenotype, ".pdf"), 
      width = 5, height = 4)
  hist(summary_stats_subset3$P)
  dev.off()
  
  ## Order variants by BP per CHR 
  summary_stats_subset3 = summary_stats_subset3 %>% 
    arrange(as.numeric(CHR), as.numeric(BP)) %>% 
    mutate(plot_pos = 1:nrow(summary_stats_subset3))
  
  df <- summary_stats_subset3
  
  ## Params
  chr_colors = rep(c("grey70", "grey20"), 11)
  names(chr_colors) <- paste(1:22)
  df$CHR <- factor(df$CHR, levels = paste(1:22))
  
  ## Signif line with Bonferroni correction
  genomewideline <- 0.05 / (10^6)  ## (6x10^6 tests actually)  
  ## Suggestive line
  suggestiveline <- genomewideline / 0.005

  ## Variants to highlight
  SNPs_to_highlight <- as.character(df$SNP[df$P < suggestiveline])
  
  ## rsID label for top 10 most signif variants per CHR 
  top_SNPs <- df %>% 
    filter(P < suggestiveline) %>% 
    arrange(as.numeric(CHR), P) %>% group_by(CHR) %>% slice_head(n = 5) %>% 
    ungroup() %>%  select(SNP) %>% as.vector()
    
  df$label <- sapply(df$SNP, function(snp){if(snp %in% top_SNPs$SNP){ if(snp %in% variant_IDs$SNP){variant_IDs[match(snp, variant_IDs$SNP), "rsID"]} else{snp} } else{NA}})
   
  ## Top SNP x chr above suggestive line 
  lead_SNPs <- df %>% 
    filter(P < suggestiveline) %>% 
    arrange(as.numeric(CHR), P) %>% group_by(CHR) %>% slice_head(n = 1) %>% 
    ungroup() %>%  select(SNP)  %>% unlist()
  lead_SNPs_labs <- sub("_", ">", sub("_", " ", sub("_", ":", sub("", "chr", lead_SNPs))))
  names(lead_SNPs_labs) <- lead_SNPs
  
  ## Position of x-axis CHR labels 
  x = tapply(df$plot_pos, df$CHR, function(pos){floor((max(pos) - min(pos) + 1)/2)}) 
  x = tapply(df$plot_pos, df$CHR, min)  +  x
  
  g1 <- ggplot(data = df, aes(plot_pos, -log10(P))) +
          geom_point(data = df[! df$SNP %in% SNPs_to_highlight, ], 
                     aes(colour = as.factor(CHR)), size = 1) + 
          scale_color_manual(values = chr_colors) + 
          labs(title = paste("GWAS for proliferation from", titles[phenotype]), 
               x = "Chromosome", 
               y = expression(-log[10](italic(p)))) +
          geom_point(data = df[df$SNP %in% SNPs_to_highlight, ], 
                     shape = 21, fill = "orange", color = "black", size = 1) +
          scale_y_continuous(expand = c(0, 0), limit = c(0, -log10(genomewideline) + 0.2),
                             breaks = seq(from = 0, to = floor(-log10(genomewideline)), by = 1)) +
          scale_x_continuous(breaks = x, labels = 1:22) +
          geom_hline(yintercept = -log10(suggestiveline), color = "gray70", linetype = "dashed") + 
          geom_hline(yintercept = -log10(genomewideline), color = "red", linetype = "dashed") + 
          geom_text_repel(aes(label = label),
                          max.overlaps = Inf, min.segment.length = 0, 
                          size = 2.7, segment.size = 0.3, nudge_y = -0.08) +
          theme(legend.position = "none",
                axis.line.y = element_line(linewidth = 0.4, color = "black"),
                axis.line.x = element_line(linewidth = 0.4, color = "black"),
                panel.background = element_blank(), 
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                plot.title = element_text(hjust = (0.5), size =  12, face = "bold"),
                axis.title.y = element_text(size = 10),
                axis.title.x = element_text(size = 10),
                axis.text = element_text(size = 9))

  ## Simpler plots for manuscript
  g2 <- ggplot(data = df, aes(plot_pos, -log10(P))) +
    geom_point(data = df[! df$SNP %in% SNPs_to_highlight, ], 
               aes(colour = as.factor(CHR)), size = 1, alpha = 0.6) + 
    scale_color_manual(values = chr_colors) + 
    labs(title = paste("GWAS for proliferation from", titles[phenotype]), 
         x = "Chromosome", 
         y = expression(-log[10](italic(p)))) +
    geom_point(data = df[df$SNP %in% SNPs_to_highlight, ], 
               color = "orange", size = 1, alpha = 0.6) +
    scale_y_continuous(expand = c(0, 0), limit = c(0, -log10(genomewideline) + 0.2),
                       breaks = seq(from = 0, to = floor(-log10(genomewideline)), by = 1)) +
    scale_x_continuous(breaks = x, labels = 1:22) +
    geom_hline(yintercept = -log10(suggestiveline), color = "gray70", linetype = "dotted") + 
    geom_hline(yintercept = -log10(genomewideline), color = "red", linetype = "dotted") + 
    geom_text_repel(data = df[df$SNP %in% lead_SNPs, ], label = lead_SNPs_labs, 
                    max.overlaps = Inf, min.segment.length = 0, 
                    size = 2.7, segment.size = 0.3, nudge_y = 0.08) +
    theme(legend.position = "none",
          axis.line.y = element_line(linewidth = 0.4, color = "black"),
          axis.line.x = element_line(linewidth = 0.4, color = "black"),
          panel.background = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = (0.5), size =  12, face = "bold"),
          axis.title.y = element_text(size = 10),
          axis.title.x = element_text(size = 10),
          axis.text = element_text(size = 9))
  
  return(list(g1, g2))  
}


## Regional plots x chr

## Load Ensembl data base
ah <- AnnotationHub()
query(ah, c("EnsDb", "Homo sapiens"))
## Fetch Ensembl 106 EnsDb for Homo sapiens
ensDb_v106 <- ah[["AH100643"]]

## Genomic zoom plots encompassing all variants above suggestive line in chr
locuszoom_plot <- function(chr, phenotype){
  
  summary_stats <- eval(parse_expr(paste0("summary_stats_", phenotype)))
  
  ## All significant variants in chr
  lead_variants <- summary_stats %>% dplyr::filter(CHR == chr & P<=1e-05) 
  
  if(dim(lead_variants)[1] == 0){
    return(NULL)
  }

  ## Regional plots around each lead SNP until all variants above suggestive line are covered
  signif_variants_left <- lead_variants
  num_signif_variants_left <- length(signif_variants_left$label)
  covered_variants <- ""
  
  plots <- list()
  gene_tracks <- list()
  
  while(num_signif_variants_left>0){
    
    ## Lead variant
    lead_var <- signif_variants_left %>% arrange(P) %>% .[1,"label"]
    
    ## Exclude already covered variants
    summary_stats <- summary_stats[! summary_stats$label %in% covered_variants, ]
    ## Subset to genomic range
    loc <- locus(data = summary_stats, 
                 chrom = "CHR", pos = "BP", p = "P", 
                 index_snp = lead_var,
                 labs = "label", 
                 fix_window = 2.5e5,
                 ens_db = ensDb_v106)
    
    genomewideline = 5e-08
    suggestiveline <- genomewideline / 0.005
    
    loc$data$suggestive_signif <- factor(loc$data$P<=suggestiveline) %>% as.character()
    loc$data$lead <- as.character(loc$data$label == lead_var & loc$data$suggestive_signif == "TRUE")
    loc$data = loc$data %>% arrange(P)
    
    ## Variants covered
    covered_variants <- loc$data$label
    ## Signif variants not covered yet
    signif_variants_left <- signif_variants_left[! signif_variants_left$label %in% covered_variants, ]
    num_signif_variants_left <- dim(signif_variants_left)[1]
    
    ## Signif colors 
    colors <- c("FALSE" = "gray", "TRUE"="yellow2")
    
    p <- ggplot(data = loc$data, aes(x = BP/1e+06, y = -log10(P), fill = suggestive_signif)) +
      geom_point(shape = 21, color = "grey50", size = 1.65) + 
      scale_fill_manual(values = colors) + 
      geom_point(data = loc$data[which.min(loc$data$P), ], 
                 shape = 21, color = "black", fill = "orangered", 
                 size = 1.65, stroke = 0.5) +
      geom_text_repel(data = subset(loc$data, suggestive_signif == "TRUE" & lead == "FALSE")[1:9,], aes(label = label),
                      max.overlaps = Inf, min.segment.length = 0, 
                      size = 2.5, segment.size = 0.3, 
                      nudge_y = -0.08, fontface = "bold") +
      geom_text_repel(data = subset(loc$data, lead == "TRUE"), aes(label = label),
                      max.overlaps = 0, min.segment.length = 0, 
                      size = 2.5, segment.size = 0.6, 
                      nudge_y = 0.45, fontface = "bold") +
      geom_hline(yintercept = -log10(suggestiveline), color = "gray70", linetype = "dashed") +
      geom_hline(yintercept = -log10(genomewideline), color = "indianred2", linetype = "dashed") +
      scale_y_continuous(expand = c(0, 0), limit = c(0, -log10(genomewideline) + 0.4),
                         breaks = seq(from = 0, to = floor(-log10(genomewideline)), by = 1)) +
      labs(title = paste0("Variants around lead SNP ", lead_var), 
           x = paste("Chromosome", chr, "(Mb)"), 
           y = expression(-log[10](italic(p)))) +
      theme_classic() +
      theme(legend.position = "none",      
            plot.title = element_text(hjust = (0.5), size =  9, face = "bold"),
            axis.title.y = element_text(size = 10),
            axis.title.x = element_text(size = 10),
            axis.text = element_text(size = 9))
    
    ## Add gene tracks in 2.5e5 window around lead variant
    g <- gg_genetracks(loc) +  
      theme(axis.title.y = element_text(size = 10),
            axis.title.x = element_text(size = 10),
            axis.text = element_text(size = 9))
    
    plots[[lead_var]] = plot_grid(plotlist = list(p, g), ncol = 1, align = "v", 
                                  rel_heights = c(1, 0.33))
    gene_tracks[[lead_var]] <- loc$TX
    
  }
  
  multiple_plots <- plot_grid(plotlist = plots, align = "vh", ncol = 3)
  ggsave(filename = paste0(plotdir, "/regional_plots/", phenotype, "/locus_zoom_chr_", chr, ".pdf"), width = 15, height = 5.2*ceiling(length(plots) / 3))
  
  return_plots <- plot_grid(plotlist = plots, align = "vh", ncol = 3, nrow = 2)
  return(list(return_plots, gene_tracks))
}



################################################################################
##                    iPSC -> young preMac proliferation 
################################################################################

phenotype <- "premac_vs_iPSC"
qqplots <- list()
qqplots[[1]] <- Q_Q_plot(phenotype)
manhattans  <- Manhattan_plot(phenotype)
ggsave(manhattans[[1]], filename = paste0(plotdir, "/Manhattan_plot_", phenotype, ".pdf"), width = 15, height = 6)
ggsave(manhattans[[2]], filename = paste0(plotdir, "/Manhattan_plot_", phenotype, "_for_publication.pdf"), width = 15, height = 6)


## Locus zoom plots 
gene_tracks_across_chrs_premac_vs_iPSC <- list()
plots <- list()

for(chr in 1:22){
  results <- locuszoom_plot(chr, phenotype)
  plots[[chr]] <- results[[1]]
  gene_tracks_across_chrs_premac_vs_iPSC[[chr]] <- results[[2]]
}

save(gene_tracks_across_chrs_premac_vs_iPSC, file = paste0(outdir, "/closest_genes_leadSNPs_premac_vs_iPSC.Rdata"))
pdf(paste0(plotdir, "/regional_plots/", phenotype, "/locus_zoom_plots.pdf"), width = 15, height = 11)
plots[which(! lapply(plots, is.null) %>% unlist())]
dev.off()


################################################################################
##                 young preMac -> old preMac proliferation
################################################################################

phenotype <- "old_vs_young_premac"
qqplots[[2]] <- Q_Q_plot(phenotype)
manhattans  <- Manhattan_plot(phenotype)
ggsave(manhattans[[1]], filename = paste0(plotdir, "/Manhattan_plot_", phenotype, ".pdf"), width = 15, height = 6)
ggsave(manhattans[[2]], filename = paste0(plotdir, "/Manhattan_plot_", phenotype, "_for_publication.pdf"), width = 15, height = 6)

gene_tracks_across_chrs_old_vs_young_premac <- list()
plots <- list()

for(chr in 1:22){
  results <- locuszoom_plot(chr, phenotype)
  plots[[chr]] <- results[[1]]
  gene_tracks_across_chrs_old_vs_young_premac[[chr]] <- results[[2]]
}

save(gene_tracks_across_chrs_old_vs_young_premac, file = paste0(outdir, "/closest_genes_leadSNPs_old_vs_young_premac.Rdata"))
pdf(paste0(plotdir, "/regional_plots/", phenotype, "/locus_zoom_plots.pdf"), width = 15, height = 11)
plots[which(! lapply(plots, is.null) %>% unlist())]
dev.off()


################################################################################
##           preMac -> microglia (in IFN/LPS/untreated) proliferation
################################################################################

#############################   Microglia in IFN   #############################
phenotype <- "microglia_vs_premac_IFN"
qqplots[[3]] <- Q_Q_plot(phenotype)
manhattans  <- Manhattan_plot(phenotype)
ggsave(manhattans[[1]], filename = paste0(plotdir, "/Manhattan_plot_", phenotype, ".pdf"), width = 15, height = 6)
ggsave(manhattans[[2]], filename = paste0(plotdir, "/Manhattan_plot_", phenotype, "_for_publication.pdf"), width = 15, height = 6)

## Regional plots and closest genes to lead SNPs in all chrs
gene_tracks_across_chrs_microglia_vs_premac_IFN <- list()
plots <- list()

for(chr in 1:22){
  results <- locuszoom_plot(chr, phenotype)
  plots[[chr]] <- results[[1]]
  gene_tracks_across_chrs_microglia_vs_premac_IFN[[chr]] <- results[[2]]
}

save(gene_tracks_across_chrs_microglia_vs_premac_IFN, file = paste0(outdir, "/closest_genes_leadSNPs_microglia_vs_premac_IFN.Rdata"))
pdf(paste0(plotdir, "/regional_plots/", phenotype, "/locus_zoom_plots.pdf"), width = 15, height = 11)
plots[which(! lapply(plots, is.null) %>% unlist())]
dev.off()

#############################   Microglia in LPS   #############################
phenotype <- "microglia_vs_premac_LPS"
qqplots[[4]] <- Q_Q_plot(phenotype)
manhattans  <- Manhattan_plot(phenotype)
ggsave(manhattans[[1]], filename = paste0(plotdir, "/Manhattan_plot_", phenotype, ".pdf"), width = 15, height =6)
ggsave(manhattans[[2]], filename = paste0(plotdir, "/Manhattan_plot_", phenotype, "_for_publication.pdf"), width = 15, height =6)

gene_tracks_across_chrs_microglia_vs_premac_LPS <- list()
plots <- list()

for(chr in 1:22){
  results <- locuszoom_plot(chr, phenotype)
  plots[[chr]] <- results[[1]]
  gene_tracks_across_chrs_microglia_vs_premac_LPS[[chr]] <- results[[2]]
}

save(gene_tracks_across_chrs_microglia_vs_premac_LPS, file = paste0(outdir, "/closest_genes_leadSNPs_microglia_vs_premac_LPS.Rdata"))
pdf(paste0(plotdir, "/regional_plots/", phenotype, "/locus_zoom_plots.pdf"), width = 15, height = 11)
plots[which(! lapply(plots, is.null) %>% unlist())]
dev.off()

##########################   Microglia in untreated   ##########################
phenotype <- "microglia_vs_premac_untreated"
qqplots[[5]] <- Q_Q_plot(phenotype)
plot_grid(plotlist = qqplots, ncol = 5)
ggsave(filename = paste(plotdir, "Q_Q_plots.pdf", sep = "/"), width = 20, height = 4)

manhattans  <- Manhattan_plot(phenotype)
ggsave(manhattans[[1]], filename = paste0(plotdir, "/Manhattan_plot_", phenotype, ".pdf"), width = 15, height = 6)
ggsave(manhattans[[2]], filename = paste0(plotdir, "/Manhattan_plot_", phenotype, "_for_publication.pdf"), width = 15, height = 6)

gene_tracks_across_chrs_microglia_vs_premac_untreated <- list()
plots <- list()

for(chr in 1:22){
  results <- locuszoom_plot(chr, phenotype)
  plots[[chr]] <- results[[1]]
  gene_tracks_across_chrs_microglia_vs_premac_untreated[[chr]] <- results[[2]]
}

save(gene_tracks_across_chrs_microglia_vs_premac_untreated, file = paste0(outdir, "/closest_genes_leadSNPs_microglia_vs_premac_untreated.Rdata"))
pdf(paste0(plotdir, "/regional_plots/", phenotype, "/locus_zoom_plots.pdf"), width = 15, height = 11)
plots[which(! lapply(plots, is.null) %>% unlist())]
dev.off()


## Raw pval histograms
hists <- list()
for(i in 1:length(phenotypes)){
  
  summary_stats <- eval(parse_expr(paste0("summary_stats_", phenotypes[i])))
  
  hists[[i]] <- ggplot(summary_stats, aes(x = `Pr(>|t|)`)) +
    geom_histogram(bins = 100, colour = "black", fill = "gray90") +
    labs(x = "Raw pval", Y = "Count", title = titles[phenotypes[i]]) +
    theme_bw()
}

plot_grid(plotlist = hists, ncol = 5, nrow = 1)
ggsave(filename = paste(plotdir, "Raw_pvals_histograms.pdf", sep = "/"))


## Corrected pvals
summary(summary_stats_premac_vs_iPSC$p_bonf)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.4288  1.0000  1.0000  1.0000  1.0000  1.0000 

summary(summary_stats_old_vs_young_premac$p_bonf)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#     1       1       1       1       1       1 

summary(summary_stats_microglia_vs_premac_IFN$p_bonf)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3904  1.0000  1.0000  1.0000  1.0000  1.0000 

summary(summary_stats_microglia_vs_premac_LPS$p_bonf)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    1       1       1       1       1       1 

summary(summary_stats_microglia_vs_premac_untreated$p_bonf)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    1       1       1       1       1       1 


## Raw pvals
summary(summary_stats_premac_vs_iPSC$`Pr(>|t|)`)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000001 0.2507765 0.5008256 0.5005265 0.7503970 0.9999995 

summary(summary_stats_old_vs_young_premac$`Pr(>|t|)`)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000011 0.2494559 0.5018708 0.5003565 0.7508763 0.9999993 

summary(summary_stats_microglia_vs_premac_IFN$`Pr(>|t|)`)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000001 0.2514116 0.5013586 0.5004907 0.7491411 0.9999994 

summary(summary_stats_microglia_vs_premac_LPS$`Pr(>|t|)`)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000002 0.2523750 0.5015611 0.5007665 0.7494996 0.9999991 

summary(summary_stats_microglia_vs_premac_untreated$`Pr(>|t|)`)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000003 0.2500746 0.4994495 0.4994690 0.7489923 0.9999996 







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
# date     2025-02-17
# rstudio  2024.04.0+735 Chocolate Cosmos (server)
# pandoc   3.1.12.3 @ /opt/view/bin/ (via rmarkdown)
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version     date (UTC) lib source
# abind                  1.4-5       2016-07-21 [1] CRAN (R 4.3.1)
# AnnotationDbi          1.64.1      2023-11-03 [1] Bioconductor
# AnnotationFilter       1.26.0      2023-10-24 [1] Bioconductor
# ape                    5.7-1       2023-03-13 [1] CRAN (R 4.3.1)
# backports              1.4.1       2021-12-13 [1] CRAN (R 4.3.1)
# beachmat               2.18.1      2024-02-14 [1] Bioconductor 3.18 (R 4.3.1)
# Biobase                2.62.0      2023-10-24 [1] Bioconductor
# BiocFileCache          2.10.2      2024-03-27 [1] Bioconductor 3.18 (R 4.3.1)
# BiocGenerics           0.48.1      2023-11-01 [1] Bioconductor
# BiocIO                 1.12.0      2023-10-24 [1] Bioconductor
# BiocParallel           1.36.0      2023-10-24 [1] Bioconductor
# BiocSingular           1.18.0      2023-10-24 [1] Bioconductor
# biomaRt                2.58.2      2024-01-30 [1] Bioconductor 3.18 (R 4.3.1)
# Biostrings             2.70.3      2024-03-13 [1] Bioconductor 3.18 (R 4.3.1)
# bit                    4.0.5       2022-11-15 [1] CRAN (R 4.3.1)
# bit64                  4.0.5       2020-08-30 [1] CRAN (R 4.3.1)
# bitops                 1.0-7       2021-04-24 [1] CRAN (R 4.3.1)
# blob                   1.2.4       2023-03-17 [1] CRAN (R 4.3.1)
# boot                   1.3-30      2024-02-26 [1] CRAN (R 4.3.1)
# broom                  1.0.5       2023-06-09 [1] CRAN (R 4.3.1)
# cachem                 1.0.8       2023-05-01 [1] CRAN (R 4.3.1)
# car                    3.1-2       2023-03-30 [1] CRAN (R 4.3.1)
# carData                3.0-5       2022-01-06 [1] CRAN (R 4.3.1)
# circlize               0.4.16      2024-02-20 [1] CRAN (R 4.3.1)
# cli                    3.6.2       2023-12-11 [1] CRAN (R 4.3.1)
# cluster                2.1.6       2023-12-01 [1] CRAN (R 4.3.1)
# codetools              0.2-20      2024-03-31 [1] CRAN (R 4.3.1)
# colorspace             2.1-0       2023-01-23 [1] CRAN (R 4.3.1)
# cowplot              * 1.1.3       2024-01-22 [1] CRAN (R 4.3.1)
# crayon                 1.5.2       2022-09-29 [1] CRAN (R 4.3.1)
# curl                   5.2.1       2024-03-01 [1] CRAN (R 4.3.1)
# data.table             1.15.4      2024-03-30 [1] CRAN (R 4.3.1)
# DBI                    1.2.2       2024-02-16 [1] CRAN (R 4.3.1)
# dbplyr                 2.5.0       2024-03-19 [1] CRAN (R 4.3.1)
# DelayedArray           0.28.0      2023-10-24 [1] Bioconductor
# digest                 0.6.35      2024-03-11 [1] CRAN (R 4.3.1)
# dplyr                * 1.1.4       2023-11-17 [1] CRAN (R 4.3.1)
# edgeR                  4.0.16      2024-02-18 [1] Bioconductor 3.18 (R 4.3.1)
# ensembldb              2.26.0      2023-10-24 [1] Bioconductor
# evaluate               0.23        2023-11-01 [1] CRAN (R 4.3.1)
# fansi                  1.0.6       2023-12-08 [1] CRAN (R 4.3.1)
# farver                 2.1.1       2022-07-06 [1] CRAN (R 4.3.1)
# fastmap                1.1.1       2023-02-24 [1] CRAN (R 4.3.1)
# filelock               1.0.3       2023-12-11 [1] CRAN (R 4.3.1)
# forcats              * 1.0.0       2023-01-29 [1] CRAN (R 4.3.1)
# generics               0.1.3       2022-07-05 [1] CRAN (R 4.3.1)
# GenomeInfoDb           1.38.8      2024-03-15 [1] Bioconductor 3.18 (R 4.3.1)
# GenomeInfoDbData       1.2.11      2025-02-07 [1] Bioconductor
# GenomicAlignments      1.38.2      2024-01-16 [1] Bioconductor 3.18 (R 4.3.1)
# GenomicFeatures        1.54.4      2024-03-13 [1] Bioconductor 3.18 (R 4.3.1)
# GenomicRanges          1.54.1      2023-10-29 [1] Bioconductor
# gggrid                 0.2-0       2022-01-11 [1] CRAN (R 4.3.1)
# ggplot2              * 3.5.1       2024-04-23 [1] CRAN (R 4.3.1)
# ggpubr                 0.6.0       2023-02-10 [1] CRAN (R 4.3.1)
# ggrepel              * 0.9.5       2024-01-10 [1] CRAN (R 4.3.1)
# ggsignif               0.6.4       2022-10-13 [1] CRAN (R 4.3.1)
# GlobalOptions          0.1.2       2020-06-10 [1] CRAN (R 4.3.1)
# glue                   1.7.0       2024-01-09 [1] CRAN (R 4.3.1)
# gtable                 0.3.4       2023-08-21 [1] CRAN (R 4.3.1)
# hms                    1.1.3       2023-03-21 [1] CRAN (R 4.3.1)
# htmltools              0.5.8       2024-03-25 [1] CRAN (R 4.3.1)
# htmlwidgets            1.6.4       2023-12-06 [1] CRAN (R 4.3.1)
# httr                   1.4.7       2023-08-15 [1] CRAN (R 4.3.1)
# IRanges                2.36.0      2023-10-24 [1] Bioconductor
# irlba                  2.3.5.1     2022-10-03 [1] CRAN (R 4.3.1)
# jsonlite               1.8.8       2023-12-04 [1] CRAN (R 4.3.1)
# jtools                 2.2.2       2023-07-11 [1] CRAN (R 4.3.1)
# KEGGREST               1.42.0      2023-10-24 [1] Bioconductor
# knitr                  1.45        2023-10-30 [1] CRAN (R 4.3.1)
# labeling               0.4.3       2023-08-29 [1] CRAN (R 4.3.1)
# lattice                0.22-6      2024-03-20 [1] CRAN (R 4.3.1)
# lazyeval               0.2.2       2019-03-15 [1] CRAN (R 4.3.1)
# LDlinkR                1.4.0       2024-04-10 [1] CRAN (R 4.3.1)
# lifecycle              1.0.4       2023-11-07 [1] CRAN (R 4.3.1)
# limma                  3.58.1      2023-10-31 [1] Bioconductor
# lme4                   1.1-35.2    2024-03-28 [1] CRAN (R 4.3.1)
# lmerTest               3.1-3       2020-10-23 [1] CRAN (R 4.3.1)
# locfit                 1.5-9.9     2024-03-01 [1] CRAN (R 4.3.1)
# locuszoomr             0.2.1       2024-02-17 [1] CRAN (R 4.3.1)
# lubridate            * 1.9.3       2023-09-27 [1] CRAN (R 4.3.1)
# magrittr               2.0.3       2022-03-30 [1] CRAN (R 4.3.1)
# MASS                   7.3-60.0.1  2024-01-13 [1] CRAN (R 4.3.1)
# Matrix                 1.6-5       2024-01-11 [1] CRAN (R 4.3.1)
# MatrixGenerics         1.14.0      2023-10-24 [1] Bioconductor
# matrixStats            1.2.0       2023-12-11 [1] CRAN (R 4.3.1)
# memoise                2.0.1       2021-11-26 [1] CRAN (R 4.3.1)
# mgcv                   1.9-1       2023-12-21 [1] CRAN (R 4.3.1)
# minqa                  1.2.6       2023-09-11 [1] CRAN (R 4.3.1)
# munsell                0.5.1       2024-04-01 [1] CRAN (R 4.3.1)
# nlme                   3.1-164     2023-11-27 [1] CRAN (R 4.3.1)
# nloptr                 2.0.3       2022-05-26 [1] CRAN (R 4.3.1)
# numDeriv               2016.8-1.1  2019-06-06 [1] CRAN (R 4.3.1)
# pander                 0.6.5       2022-03-18 [1] CRAN (R 4.3.1)
# patchwork              1.2.0       2024-01-08 [1] CRAN (R 4.3.1)
# permute                0.9-7       2022-01-27 [1] CRAN (R 4.3.1)
# pillar                 1.9.0       2023-03-22 [1] CRAN (R 4.3.1)
# pinfsc50               1.3.0       2023-12-05 [1] CRAN (R 4.3.1)
# pkgconfig              2.0.3       2019-09-22 [1] CRAN (R 4.3.1)
# pkgload                1.3.4       2024-01-16 [1] CRAN (R 4.3.1)
# plotly                 4.10.4      2024-01-13 [1] CRAN (R 4.3.1)
# png                    0.1-8       2022-11-29 [1] CRAN (R 4.3.1)
# prettyunits            1.2.0       2023-09-24 [1] CRAN (R 4.3.1)
# progress               1.2.3       2023-12-06 [1] CRAN (R 4.3.1)
# ProtGenerics           1.34.0      2023-10-24 [1] Bioconductor
# purrr                * 1.0.2       2023-08-10 [1] CRAN (R 4.3.1)
# R6                     2.5.1       2021-08-19 [1] CRAN (R 4.3.1)
# ragg                   1.3.0       2024-03-13 [1] CRAN (R 4.3.1)
# rappdirs               0.3.3       2021-01-31 [1] CRAN (R 4.3.1)
# Rcpp                   1.0.12      2024-01-09 [1] CRAN (R 4.3.1)
# RCurl                  1.98-1.14   2024-01-09 [1] CRAN (R 4.3.1)
# readr                * 2.1.5       2024-01-10 [1] CRAN (R 4.3.1)
# restfulr               0.0.15      2022-06-16 [1] CRAN (R 4.3.1)
# rjson                  0.2.21      2022-01-09 [1] CRAN (R 4.3.1)
# rlang                * 1.1.3       2024-01-10 [1] CRAN (R 4.3.1)
# rmarkdown              2.26        2024-03-05 [1] CRAN (R 4.3.1)
# Rsamtools              2.18.0      2023-10-24 [1] Bioconductor
# RSQLite                2.3.6       2024-03-31 [1] CRAN (R 4.3.1)
# rstatix                0.7.2       2023-02-01 [1] CRAN (R 4.3.1)
# rstudioapi             0.16.0      2024-03-24 [1] CRAN (R 4.3.1)
# rsvd                   1.0.5       2021-04-16 [1] CRAN (R 4.3.1)
# rtracklayer            1.62.0      2023-10-24 [1] Bioconductor
# S4Arrays               1.2.1       2024-03-04 [1] Bioconductor 3.18 (R 4.3.1)
# S4Vectors              0.40.2      2023-11-23 [1] Bioconductor 3.18 (R 4.3.1)
# ScaledMatrix           1.10.0      2023-10-24 [1] Bioconductor
# scales                 1.3.0       2023-11-28 [1] CRAN (R 4.3.1)
# sessioninfo          * 1.2.2       2021-12-06 [1] CRAN (R 4.3.1)
# shape                  1.4.6.1     2024-02-23 [1] CRAN (R 4.3.1)
# SparseArray            1.2.4       2024-02-11 [1] Bioconductor 3.18 (R 4.3.1)
# statmod                1.5.0       2023-01-06 [1] CRAN (R 4.3.1)
# stringi                1.8.3       2023-12-11 [1] CRAN (R 4.3.1)
# stringr              * 1.5.1       2023-11-14 [1] CRAN (R 4.3.1)
# SummarizedExperiment   1.32.0      2023-10-24 [1] Bioconductor
# systemfonts            1.0.6       2024-03-07 [1] CRAN (R 4.3.1)
# textshaping            0.3.7       2023-10-09 [1] CRAN (R 4.3.1)
# tibble               * 3.2.1       2023-03-20 [1] CRAN (R 4.3.1)
# tidyr                * 1.3.1       2024-01-24 [1] CRAN (R 4.3.1)
# tidyselect             1.2.1       2024-03-11 [1] CRAN (R 4.3.1)
# tidyverse            * 2.0.0       2023-02-22 [1] CRAN (R 4.3.1)
# timechange             0.3.0       2024-01-18 [1] CRAN (R 4.3.1)
# tzdb                   0.4.0       2023-05-12 [1] CRAN (R 4.3.1)
# utf8                   1.2.4       2023-10-22 [1] CRAN (R 4.3.1)
# vcfR                   1.15.0      2023-12-08 [1] CRAN (R 4.3.1)
# vctrs                  0.6.5       2023-12-01 [1] CRAN (R 4.3.1)
# vegan                  2.6-4       2022-10-11 [1] CRAN (R 4.3.1)
# viridisLite            0.4.2       2023-05-02 [1] CRAN (R 4.3.1)
# vroom                  1.6.5       2023-12-05 [1] CRAN (R 4.3.1)
# withr                  3.0.0       2024-01-16 [1] CRAN (R 4.3.1)
# xfun                   0.43        2024-03-25 [1] CRAN (R 4.3.1)
# XML                    3.99-0.16.1 2024-01-22 [1] CRAN (R 4.3.1)
# xml2                   1.3.6       2023-12-04 [1] CRAN (R 4.3.1)
# XVector                0.42.0      2023-10-24 [1] Bioconductor
# yaml                   2.3.8       2023-12-11 [1] CRAN (R 4.3.1)
# zlibbioc               1.48.2      2024-03-13 [1] Bioconductor 3.18 (R 4.3.1)
# zoo                    1.8-12      2023-04-13 [1] CRAN (R 4.3.1)
# 
# [1] /opt/view/rlib/R/library
# [2] /software/hgi/installs/softpack/rstudio/.spack/linux-ubuntu22.04-x86_64_v3/gcc-11.4.0/r-4.3.1-bfwldrk76z6f52upk47zepliekn7ayqz/rlib/R/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

