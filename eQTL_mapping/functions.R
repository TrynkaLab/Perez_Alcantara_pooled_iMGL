#functions
extract_unique_gene_name<- function(df) {
  df %>%
    dplyr::distinct(gene_name) %>%
    pull(gene_name)
}


intersection <- function(slop, inter) {
  # m1, m2 = slopes; b1, b2 = intercepts                       / 
  # m1x - y = -b1
  # m2x - y = -b2
  A =  rbind(c(slop$Est.[1], -1),  
             c(slop$Est.[2], -1))
  b =  c(-inter$Est[1],
         -inter$Est[2])
  res = solve(A,b)
  names(res) = c("x","y")
  res
}

tensorQTL_summary= function(QTL, metadata){
  metadata$ncells = as.numeric(as.character(metadata$ncells))
  result = data.frame(n_genes_10=1,n_genes_05=1,n_genes_01=1,n_variants_10=1,n_variants_05=1, n_variants_01=1, 
                      mean_n_cells = 1, median_n_cells = 1, CV_cells = 1,
                      sd_cells = 1, total_genes_tested=1)
  result$n_genes_10 = sum(QTL$qval < 0.10)
  result$n_genes_05 = sum(QTL$qval < 0.05)
  result$n_genes_01 = sum(QTL$qval < 0.01)
  result$n_variants_10 = QTL %>% 
    dplyr::filter(qval<0.10) %>%
    dplyr::filter(!duplicated(variant_id)) %>% nrow()
  result$n_variants_05 = QTL %>% 
    dplyr::filter(qval<0.05) %>%
    dplyr::filter(!duplicated(variant_id)) %>% nrow()
  result$n_variants_01 = QTL %>% 
    dplyr::filter(qval<0.01) %>%
    dplyr::filter(!duplicated(variant_id)) %>% nrow()
  result$n_donors = nrow(metadata)
  result$mean_n_cells = mean(metadata$ncells)
  result$median_n_cells = median(metadata$ncells)
  result$CV_cells = sd(metadata$ncells) / mean(metadata$ncells)
  result$sd_cells = sd(metadata$ncells) 
  result$total_genes_tested = length(unique(QTL$phenotype_id))
  return(unlist(result))
}
ensembl_to_gene_name = function(IDs,IDFrom,IDTo){
  require(AnnotationHub)
  # check ensembl releases with:
  hub = AnnotationHub()
  query(hub, c("Homo sapiens","EnsDb"))
  hub = hub[["AH83216"]]
  idmap=mapIds(x = hub,
               keys = IDs,column = IDTo,keytype = IDFrom,multiVals = "first")
  na_vec=names(idmap[is.na(idmap)==T])
  idmap=idmap[is.na(idmap)==F]
  idmap_df=data.frame("From"=names(idmap),"To"=unlist(unname(idmap)),stringsAsFactors = F)
  return(list(map=idmap_df,noMap=na_vec))
}

make_dosages_ref_alt = function(vcf){
  ## Extracting GT
  gt = as.data.table(extract.gt(vcf,return.alleles=TRUE))
  # first make new names from chromosome, position, ref and alt (so there are no duplicates)
  names = paste0(vcf@fix[,"CHROM"],"_",vcf@fix[,"POS"],"_",vcf@fix[,"REF"],"_",vcf@fix[,"ALT"])
  position = paste0(vcf@fix[,"CHROM"],"_",vcf@fix[,"POS"])
  # There are rows for same SNP but different minor alleles: those will be removed with the previous line (does not include "ALT")
  # there are still some duplicates so remove those
  gt=gt[!duplicated(names),]
  position=position[!duplicated(names)]
  
  message("Removing ",sum(duplicated(names)), " duplicated rows from vcf")
  
  ## removing those from the names
  names=names[!duplicated(names)]
  gt[, ("names") := names] 
  gt[, ("position") := position] 
  
  # setkey(gt,rn) # this sorts the data.table, be careful
  object.size(gt)
  
  return(gt)
}

# swap the two last elements in an underscore-separated string
# in this case to swap the REF and ALT alleles
# in format chr_pos_REF_ALT
swap_REF_ALT_alelle_names = function(variant_name){
  
  # Split the input string by the underscore character
  string_parts = strsplit(variant_name, "_")[[1]]
  
  # Extract the last element and swap the order of the last two letters
  last_element = string_parts[length(string_parts)]
  new_last_element =string_parts[length(string_parts)-1]
  
  # Combine the parts back together with the modified last element
  new_variant_name = paste0(paste0(string_parts[1:2], collapse = "_"), "_", last_element,"_", new_last_element)
  
  # Print the new string
  return(new_variant_name)
  
}

# plot eQTL result boxplots per variant
boxplot_eQTL = function(variant,gene,group,genotype_df,expression_df,metadata_df, color_by_pool=FALSE){
  message("Plotting variant ",variant)
  # Try-catch block to check for ref / alt swaps
  mod_variant = tryCatch({
    if (!variant %in% genotype_df$names) {
      warning(paste0("Variant ", variant, " not found, trying swapping REF and ALT"))
      
      
    }
    variant
  }, warning = function(w) {
    message("Warning: ", w)
    # for some reason it can't find the externally-defined swap_REF_ALT_alelle_names() though it's loaded
    
    # Split the input string by the underscore character
    string_parts = strsplit(variant, "_")[[1]]
    
    # Extract the last element and swap the order of the last two letters
    last_element = string_parts[length(string_parts)]
    new_last_element =string_parts[length(string_parts)-1]
    
    # Combine the parts back together with the modified last element
    new_variant_name = paste0(paste0(string_parts[1:2], collapse = "_"), "_", last_element,"_", new_last_element)
    
    # Print the new string
    new_variant_name
    
    
  }, error = function(e) {
    message("Error: ", e)
    NULL
  })
  
  # check if modified variant is found in the genotype file, stop function if not
  if(mod_variant != variant){
    if (!mod_variant %in% genotype_df$names) {
      stop(paste0("Modified variant ", mod_variant, " not found"))
    }
    else{
      variant = mod_variant
    }
  }
  
  genotype = genotype_df %>%
    dplyr::filter(names == variant) 
  
  
  ref = unlist(str_split(genotype$names,"_"))[3]
  alt = unlist(str_split(genotype$names,"_"))[4]
  hom_ref = paste0(ref,"|",ref)
  het = paste0(ref,"|",alt)
  other_het =  paste0(alt,"|",ref) # some are swapped, shouldn't matter
  hom_alt = paste0(alt,"|",alt)
  
  
  genotype = genotype %>%
    pivot_longer(cols = !c(names,position),
                 names_to = "line",
                 values_to = "genotype") %>%
    mutate(genotype = if_else(genotype==other_het,het,genotype)) %>% # renaming for levels in plot
    mutate(genotype = factor(genotype,
                             levels = c(hom_ref,het,hom_alt),
                             ordered = TRUE))
  
  metadata_df = metadata_df %>%
    dplyr::filter(name == paste0(str_split_1(group,"_")[2:length(str_split_1(group,"_"))],collapse = "_")) %>%
    dplyr::rename(line=donor_id)
  expression = expression_df %>%
    dplyr::filter(gene_name == gene) %>%
    dplyr::select(!c(name:gene_name,proliferation_status,treatment)) %>%
    dplyr::rename(line=donor_id) %>%
    dplyr::left_join(metadata_df, by = "line") %>%
    # dplyr::mutate(pool = case_when(is.na(pool) == TRUE ~ "shared",
    #                         is.na(pool) == FALSE ~ pool)) %>% 
    dplyr::mutate(ncells = as.numeric(trimws(ncells))) %>%
    left_join(genotype, by = "line")
  
  
  # Create the plot with box plot and points
  
  if(color_by_pool == TRUE){
    p = expression %>% 
      ggplot(aes(x = genotype, y = expression)) +
      geom_boxplot() + # Create the box plot
      geom_jitter(aes(size = ncells, col = pool), alpha = 0.6,width = 0.2) + # Create points with size based on 'size_var'
      theme_classic() +
      ggtitle(paste0(unique(expression$names)," - ", unique(expression$gene_name), " gene pair")) + 
      ylab(paste0("log2-normalised pseudobulk counts for ", unique(expression$gene_name)))
    
  } else{
    p = expression %>% 
      ggplot(aes(x = genotype, y = expression)) +
      geom_boxplot() + # Create the box plot
      geom_jitter(aes(size = ncells), alpha = 0.6,width = 0.2) + # Create points with size based on 'size_var'
      theme_classic() +
      ggtitle(paste0(unique(expression$names)," - ", gene, " gene pair")) + 
      ylab(paste0("log2-normalised pseudobulk counts for ",gene))
    plot(p)
  }
  
  
  return(p)
  
}

# same eQTL boxplot, interactive

interactive_boxplot_eQTL = function(variant,gene,group,genotype_df,expression_df,metadata_df, color_by_pool=FALSE){
  message("Plotting variant ",variant)
  # Try-catch block to check for ref / alt swaps
  mod_variant = tryCatch({
    if (!variant %in% genotype_df$names) {
      warning(paste0("Variant ", variant, " not found, trying swapping REF and ALT"))
      
      
    }
    variant
  }, warning = function(w) {
    message("Warning: ", w)
    # for some reason it can't find the externally-defined swap_REF_ALT_alelle_names() though it's loaded
    
    # Split the input string by the underscore character
    string_parts = strsplit(variant, "_")[[1]]
    
    # Extract the last element and swap the order of the last two letters
    last_element = string_parts[length(string_parts)]
    new_last_element =string_parts[length(string_parts)-1]
    
    # Combine the parts back together with the modified last element
    new_variant_name = paste0(paste0(string_parts[1:2], collapse = "_"), "_", last_element,"_", new_last_element)
    
    # Print the new string
    new_variant_name
    
    
  }, error = function(e) {
    message("Error: ", e)
    NULL
  })
  
  # check if modified variant is found in the genotype file, stop function if not
  if(mod_variant != variant){
    if (!mod_variant %in% genotype_df$names) {
      stop(paste0("Modified variant ", mod_variant, " not found"))
    }
    else{
      variant = mod_variant
    }
  }
  
  genotype = genotype_df %>%
    dplyr::filter(names == variant) 
  
  
  ref = unlist(str_split(genotype$names,"_"))[3]
  alt = unlist(str_split(genotype$names,"_"))[4]
  hom_ref = paste0(ref,"|",ref)
  het = paste0(ref,"|",alt)
  other_het =  paste0(alt,"|",ref) # some are swapped, shouldn't matter
  hom_alt = paste0(alt,"|",alt)
  
  
  genotype = genotype %>%
    pivot_longer(cols = !c(names,position),
                 names_to = "line",
                 values_to = "genotype") %>%
    mutate(genotype = if_else(genotype==other_het,het,genotype)) %>% # renaming for levels in plot
    mutate(genotype = factor(genotype,
                             levels = c(hom_ref,het,hom_alt),
                             ordered = TRUE))
  
  metadata_df = metadata_df %>%
    mutate(line = donor_id)
  expression = expression_df %>%
    dplyr::filter(gene_name == gene) %>%
    dplyr::select(!c(chr,  start ,   end)) %>%
    pivot_longer(cols = !c(gene_name,gene_name),
                 names_to = "line",
                 values_to = "logNorm_expr") %>%
    left_join(metadata_df, by = "line") %>%
    mutate(pool = case_when(is.na(pool) == TRUE ~ "shared",
                            is.na(pool) == FALSE ~ pool),
           ncells = as.numeric(trimws(ncells))) %>%
    left_join(genotype, by = "line")
  
  
  # Create the plot with box plot and points
  
  if(color_by_pool == TRUE){
    
    
    # Convert genotype to numeric and add jitter
    
    p = plot_ly(expression, x = ~genotype, y = ~logNorm_expr,
                text = ~paste("Line: ", line, '<br>ncells:', ncells), 
                size = ~ncells, color = ~pool,
                type = 'scatter', mode = 'markers',
                #Choosing the range of the bubbles' sizes:
                sizes = c(5, 50),
                marker = list(opacity = 0.5, sizemode = 'diameter')) %>%
      layout(title = expression$gene_variant_pair)
    
    # Combine scatter plot and layout
    
    
    
  }
  
  # Convert ggplot to plotly object
  
  return(p)
  
}

########################################
####### colocalization functions #######
######################################


runColoc = function(region, p1, p2, p12,gwas_summary_stats,eqtl_summary_stats) {
  #' Perform colocalization analysis between eQTLs and GWAS loci
  #'
  #' This function performs colocalization analysis between eQTLs and GWAS loci within a specified region.
  #'
  #' @param region A data frame containing information about the region of interest.
  #' @param p1 The prior probability of association for dataset 1 (eQTLs).
  #' @param p2 The prior probability of association for dataset 2 (GWAS loci).
  #' @param p12 The prior probability of colocalization.
  #' @param gwas_summary_stats GWAS summary stats (dense format). Preferably harmonised.
  #' @param eqtl_summary_stats eQTL summary stats (dense format).  Preferably harmonised, make sure it's the same genome build 
  #' as the GWAS summary stats.
  #' @return A numeric vector containing the colocalization posterior probabilities.
  #'
  #' @export
  #' 
  gene = as.character(region["gene_name"])
  variant = as.character(region["variant_id"])
  chr = as.integer(region["chr"])
  start = as.integer(region["start"])
  end = as.integer(region["end"])
  
  # dplyr::filtering GWAS variants within the window
  gwas_variants_in_window = gwas_summary_stats %>%
    dplyr::filter(chr == chr, pos >= start, pos <= end) %>%
    dplyr::filter(!is.na(beta)) %>%
    distinct(variant_id)
  
  # dplyr::filtering eQTL variants within the window
  eqtl_variants_in_window = eqtl_summary_stats %>%
    dplyr::filter(gene_name == gene, chr == chr, pos >= start, pos <= end) %>%
    distinct(variant_id)
  
  # Counting the number of shared variants
  number_of_shared_variants = sum(eqtl_variants_in_window$variant_id %in% gwas_variants_in_window$variant_id)
  
  if (number_of_shared_variants == 0) {
    print("No shared variants identified between eQTL and GWAS at this locus")
    coloc_posteriors = c(gene_name = gene, variant_id = variant, rep(NA, 13))
  } else {
    coloc_results = coloc::coloc.signals(
      dataset1 = list(
        pvalues = eqtl_variants_in_window$p_nominal,
        N = eqtl_sample_size,
        MAF = eqtl_variants_in_window$maf,
        beta = eqtl_variants_in_window$beta,
        varbeta = eqtl_variants_in_window$var_beta,
        type = "quant",
        snp = eqtl_variants_in_window$variant_id
      ),
      dataset2 = list(
        beta = gwas_variants_in_window$beta,
        varbeta = gwas_variants_in_window$var_beta^2,
        type = "cc",
        N = gwas_sample_size,
        s = gwas_case_control_ratio,
        snp = gwas_variants_in_window$variant_id
      ),
      method = "single",
      p1 = p1,
      p2 = p2,
      p12 = p12
    )
    
    coloc_posteriors = c(gene_name = gene, variant_id = variant, unlist(coloc_results$summary))
  }
  
  return(coloc_posteriors)
}

myManhattan = function(df, graph.title = "", highlight = NULL, highlight.col = "green",
                       col = c("lightblue", "navy"), even.facet = FALSE, chrom.lab = NULL,
                       suggestiveline = 1e-05, suggestivecolor = "blue",
                       genomewideline = 5e-08, genomewidecolor = "red",
                       font.size = 12, axis.size = 0.5, significance = NULL, report = FALSE,
                       inf.corr = 0.95, y.step = 2, point.size = 1){
  # modified from https://github.com/alfonsosaera/myManhattan
  myMin <- min(df$P[df$P != 0]) * inf.corr
  df$P[df$P == 0] <- myMin
  require(ggplot2)
  require(stats)
  y.title <- expression(-log[10](italic(p)))
  if (length(col) > length(unique(df$CHR))){
    chrom.col <- col[1:length(unique(df$CHR))]
  } else if (!(length(col) > length(unique(df$CHR)))){
    chrom.col <- rep(col, length(unique(df$CHR))/length(col))
    if (length(chrom.col) < length(unique(df$CHR))){
      dif <- length(unique(df$CHR)) - length(chrom.col)
      chrom.col <- c(chrom.col, col[1:dif])
    }
  }
  y.max <- floor(max(-log10(df$P))) + 1
  if (y.max %% 2 != 0){
    y.max <- y.max + 1
  }
  if (!is.null(chrom.lab)){
    if (length(unique(df$CHR)) != length(chrom.lab)){
      warning("Number of chrom.lab different of number of chromosomes in dataset, argument ignored.")
    } else {
      df$CHR <- factor(df$CHR, levels = unique(df$CHR), labels=chrom.lab)
    }
  }
  g <- ggplot(df) +
    geom_point(aes(BP, -log10(P), colour = as.factor(CHR)), size = point.size)
  if (!is.null(significance)){
    if (is.numeric(significance)){
      genomewideline <- significance
      suggestiveline <- genomewideline / 0.005
    } else if (significance == "Bonferroni"){
      BFlevel <- 0.05 / length(df$SNP)
      cat("Bonferroni correction significance level:", BFlevel, "\n")
      genomewideline <- BFlevel
      suggestiveline <- BFlevel / 0.005
    } else if (significance == "FDR"){
      df$fdr <- p.adjust(df$P, "fdr")
      genomewideline <- 0.05
      suggestiveline <- FALSE
      y.title <- expression(-log[10](italic(q)))
      g <- ggplot(df) +
        geom_point(aes(BP, -log10(fdr), colour = as.factor(CHR)), size = point.size)
      if (!is.null(highlight)) {
        if (is.numeric(highlight)){
          highlight <- as.character(df$SNP[df$P < highlight])
        }
        if (any(!(highlight %in% df$SNP))){
          warning("Cannot highlight SNPs not present in the dataset. Argument is ignored.")
        } else {
          g <- g + geom_point(data = df[which(df$SNP %in% highlight), ],
                              aes(BP, -log10(fdr), group=SNP, colour=SNP),
                              color = highlight.col, size = point.size)
          highlight <- NULL
          y.max <- floor(max(-log10(df$fdr))) + 1
          if (y.max %% 2 != 0){
            y.max <- y.max + 1
          }
        }
      }
    }
  }
  if (even.facet){
    g <- g + facet_grid(.~CHR, scale = "free_x", switch = "x")
  } else {
    g <- g + facet_grid(.~CHR, scale = "free_x", space = "free_x", switch = "x")
  }
  g <- g + scale_colour_manual(values = chrom.col) +
    scale_y_continuous(expand = c(0, 0), limit = c(0, y.max),
                       breaks = seq(from = 0, to = y.max, by = y.step)) +
    scale_x_continuous() +
    ggpubr::theme_pubr()+
    theme(strip.background = element_blank(), legend.position = "none",
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.spacing.x=unit(0.1, "lines"),
          axis.line.y = element_line(size = axis.size, color = "black"),
          axis.ticks.y = element_line(size = axis.size, color = "black"),
          axis.ticks.length = unit(axis.size * 10, "points"),
          plot.title = element_text(hjust = (0.5), size = font.size + 8),
          axis.title.y = element_text(size = font.size + 5),
          axis.title.x = element_text(size = font.size + 5),
          axis.text = element_text(size = font.size),
          strip.text.x = element_text(size = font.size))+
    labs(title = graph.title, x = "Chromosome", y = y.title)
  if (!is.null(highlight)) {
    if (is.numeric(highlight)){
      highlight <- as.character(df$SNP[df$P < highlight])
    }
    if (any(!(highlight %in% df$SNP))){
      warning("Cannot highlight SNPs not present in the dataset. Argument is ignored.")
    } else {
      g <- g + geom_point(data = df[which(df$SNP %in% highlight), ],
                          aes(BP, -log10(P), group=SNP, colour=SNP),
                          color = highlight.col, size = point.size)
    }
  }
  if (suggestiveline){
    g <- g + geom_hline(yintercept = -log10(suggestiveline), color = suggestivecolor,
                        linetype = "dotted")
  }
  if (genomewideline){
    g <- g + geom_hline(yintercept = -log10(genomewideline), color = genomewidecolor,
                        linetype = "dotted" ) + 
      ylim(c(0,-log10(genomewideline) + 0.5))
  }
  if (report){
    if (significance == "FDR"){
      rep <- df[df$fdr < 0.05, ]
    } else if (significance == "Bonferroni"){
      rep <- df[df$P < BFlevel, ]
    } else if (is.numeric(significance)){
      rep <- df[df$P < significance, ]
    } else {
      cat("using default significance level, 5e-8")
      rep <- df[df$P < 5e-8, ]
    }
    print(rep)
  }
  return(g)
}

## gene symbol to Emsemblr ID ############
symbol_to_ensembl <- function(x){
  ensembl = mapIds(org.Hs.eg.db, keys=x, columns=c("ENSEMBL","SYMBOL"), keytype="ENSEMBL")
  return(ensembl)
}
symbol_to_entrez <- function(x){
  ensembl = mapIds(org.Hs.eg.db, x, "ENTREZID","SYMBOL")
  return(ensembl)
}


# Hypergeometric test

# hypergeometric test with permutations

# Total balls in the urn: background dataset (all RNA-seq genes)
# x --> number of white balls drawn without replacement from an urn which contains both black and white balls (GWAS genes in my tested dataset)
# m --> White balls in the urn: characteristic I'm testing for enrichment (GWAS genes)
# n --> Black balls in the urn: total-white balls (non-GWAS genes)
# k --> Number of balls drawn from the urn: depends on the size of the dataset tested 

# I want the p-value of getting (number of GWAS genes in my DE dataset) or more white balls in a sample 
# of size (number of genes in my DE dataset) from an urn with (number of GWAS genes) white balls and 
# (total genes RNAseq-GWAS genes in that set) black balls.

# 1-phyper(q, m, n, k)

# q in this case is x-1 (diff exp dataset -1), the probability of x or more. It is a quantile defined as the smallest value x such that F(x) ≥ p, 
# where F is the distribution function.
# 


# comparison of probabilities

phyperRandom <- function(myGeneList, myGeneSet, genome){
  myRandomGS <- sample( genome,size=length(myGeneSet) )  # take random sample of genes from whole genome, same size as set tested
  myX <- length(which(myGeneList %in% myRandomGS))  # number of GWAS in random set
  myM <- length(which(myGeneList %in% genome)) #  GWAS genes in my pool
  myN <- length(genome) - length(which(myGeneList %in% genome))  # genes in genome not GWAS 
  myK <- length(myGeneSet) # size of set tested
  return(list(prob=1-phyper(q=myX-1,
                            m=myM, n=myN, k=myK), n.overlap=myX))
}

# selected locuszoom plots for coloc
plot_selected_coloc_locuszoom = function(condition,locus_name,eGene,genes,dbsnp,nudge_y = 0){
  
  
  eQTL = coloc_res[[condition]][[locus_name]][[eGene]]$results %>%
    dplyr::select(snp,V.df2,z.df2, r.df2 ,lABF.df2, SNP.PP.H4) %>%
    dplyr::rename(variant_id = snp) %>%
    tidyr::separate_wider_delim(cols = variant_id,delim = "_",names = c("chrom","pos","other_allele","effect_allele"),cols_remove = FALSE) %>%
    dplyr::mutate(Position_GRCh38 = paste0("chr",chrom,":",pos),
                  chrom = as.numeric(chrom),
                  pos = as.numeric(pos),
                  gene_name = eGene,
                  pval_nominal = 10 ^(( pnorm(-abs(z.df2),log.p=TRUE) + log(2) ) / log(10))) %>% # transforming this bc coloc uses p-vals derived from z score if betas + se from eQTL / GWAS are available
    dplyr::left_join(genes) %>%
    dplyr::left_join(dbsnp, by = "variant_id") 
  
  top_coloc_var = eQTL %>%
    dplyr::filter(gene_name == eGene) %>% # filter to eGene, beware this removes variants without rsId which I hope are not the lead (make a warning for that)
    dplyr::slice_min(pval_nominal) %>% # if there are ties by pval, sort by rsid
    dplyr::slice_min(rsid,na_rm = FALSE) %>% # if there are NAs in rsid, but not in all rows, take minimum rsid (not natural order)
    dplyr::distinct(rsid,.keep_all = TRUE)
  
  tryCatch({
    if (is.na(top_coloc_var$rsid)) stop("Warning: top coloc var has no rsID - skipping ...")
  }, error=function(e){})
  
  eQTL_subset = eQTL %>%
    dplyr::filter(gene_name == eGene & !is.na(rsid)) %>% # filter to eGene, beware this removes variants without rsId which I hope are not the lead (make a warning for that)
    data.frame()
  
  
  loc_eqtl = locuszoomr::locus(eQTL_subset, gene = eGene, flank = 2.5e5,
                               ens_db = ensDb_v106, p = "pval_nominal", labs = "rsid", pos = "pos")
  
  loc_eqtl = locuszoomr::link_LD(loc_eqtl, token = Sys.getenv("LDLINK_TOKEN"), pop="CEU",r2d = "r2",genome_build = "grch38_high_coverage")
  p_eqtl = locuszoomr::gg_scatter(loc_eqtl,labels = "index", nudge_y = nudge_y) + ggtitle("eQTL")
  
  
  GWAS = coloc_res[[condition]][[locus_name]][[eGene]][["results"]] %>%
    dplyr::select(snp,V.df1,z.df1, r.df1 ,lABF.df1, SNP.PP.H4) %>%
    dplyr::rename(variant_id = snp) %>%
    tidyr::separate_wider_delim(cols = variant_id,delim = "_",names = c("chrom","pos","other_allele","effect_allele"),cols_remove = FALSE) %>%
    dplyr::mutate(Position_GRCh38 = paste0("chr",chrom,":",pos),
                  chrom = as.numeric(chrom),
                  pos = as.numeric(pos),
                  gene_name = eGene,
                  pval_nominal = 10 ^(( pnorm(-abs(z.df1),log.p=TRUE) + log(2) ) / log(10))) %>%
    dplyr::left_join(genes) %>%
    dplyr::left_join(dbsnp, by = "variant_id") 
  
  GWAS_subset = GWAS %>%
    dplyr::filter(gene_name == eGene & !is.na(rsid)) %>% # filter to eGene, beware this removes variants without rsId which I hope are not the lead (make a warning for that)
    data.frame()
  
  
  loc = locuszoomr::locus(GWAS_subset, gene = eGene, flank = 2.5e5,
                          ens_db = ensDb_v106, p = "pval_nominal", labs = "rsid", pos = "pos",
                          index_snp = loc_eqtl$index_snp) # specifying the index SNP from the eQTL data!
  
  loc = locuszoomr::link_LD(loc, token = Sys.getenv("LDLINK_TOKEN"), pop="CEU",
                            r2d = "r2",genome_build = "grch38_high_coverage")
  p_gwas = locuszoomr::gg_scatter(loc,
                                  pcutoff = 5e-100, # setting to NULL doesn't work so I give high limit
                                  legend_pos = NULL,
                                  labels = "index", nudge_y=nudge_y) + ggtitle("GWAS")
  g = locuszoomr::gg_genetracks(loc,gene_col = "grey40",exon_col = "grey40",exon_border=NA,
                                filter_gene_biotype = "protein_coding" )
  p = (p_gwas / p_eqtl / g) + patchwork::plot_annotation(title = paste0(locus_name, " locus , ", eGene, " eGene - chr",
                                                                        top_coloc_var$chrom, " ",top_coloc_var$pos, " ",
                                                                        top_coloc_var$other_allele," \u2192 ", # unicode arrow
                                                                        top_coloc_var$effect_allele ),
                                                         subtitle =  paste0("PP for colocalization = ", 
                                                                            floor(coloc_res[[condition]][[locus_name]][[eGene]]$summary[[6]] * 100), "%")) + 
    patchwork::plot_layout(ncol = 1,nrow = 3,heights = c(2,2,0.5))
  
  return(list("p" = p, "eQTL"= eQTL,"GWAS"=GWAS))    
  
}

plot_selected_GWAS_locuszoom = function(gwas_dataset,gene,locus_name,eqtl_index_snp,nudge_y = 0){
  
  gwas_file = list.files(paste0("/lustre/scratch123/hgi/mdt1/projects/healthy_imm_expr/resources/summary_statistics/public/",gwas_dataset,"/harmonised"),
                         pattern="*.h.tsv.gz",full.names = TRUE)
  chrom = stringr::str_split_i(locus_name,pattern="_",i=1)
  GWAS = readr::read_tsv(gwas_file) %>%
    dplyr::filter(hm_chrom==chrom) %>% # to save memory
    dplyr::select(hm_variant_id,hm_rsid,hm_chrom, hm_pos ,hm_other_allele, hm_effect_allele,p_value) %>%
    dplyr::rename(variant_id = hm_variant_id,
                  rsid = hm_rsid,
                  chrom=hm_chrom,
                  pos=hm_pos, 
                  other_allele=hm_other_allele,
                  effect_allele=hm_effect_allele,
                  pval_nominal=p_value) %>%
    dplyr::filter(!is.na(rsid)) %>%
    data.frame()
  
  loc = locuszoomr::locus(GWAS, gene = gene, flank = 2.5e5,
                          ens_db = ensDb_v106, p = "pval_nominal", labs = "rsid", pos = "pos",
                          index_snp = eqtl_index_snp) # specifying the index SNP from the eQTL data!
  
  loc = locuszoomr::link_LD(loc, token = Sys.getenv("LDLINK_TOKEN"), pop="CEU",
                            r2d = "r2",genome_build = "grch38_high_coverage")
  p_gwas = locuszoomr::gg_scatter(loc,
                                  showLD = ifelse(sum(is.na(loc$data$ld)) == length(loc$data$ld), # don't plot LD if all NA
                                                  FALSE,TRUE),
                                  scheme = c("grey","darkgrey","purple"),
                                  labels = "index", nudge_y=nudge_y)
  g = locuszoomr::gg_genetracks(loc,gene_col = "grey40",exon_col = "grey40",exon_border=NA,
                                filter_gene_biotype = "protein_coding" )
  p = (p_gwas / g) + patchwork::plot_annotation(title = paste0(gwas_dataset, " GWAS ,", gene, " locus")) + 
    patchwork::plot_layout(ncol = 1,nrow = 2,heights = c(4,1))
  
  return(list("p" = p,"locus" = loc))    
  
}

plot_genetracks_with_y = function (loc, 
                                   p1_limits = p1_limits,
                                   cex.axis = 1, cex.lab = 1, cex.text = 0.7, 
                                   gene_col = "blue4", 
                                   text_pos = "top", xticks = TRUE, xlab = NULL) 
{
  if (!inherits(loc, "locus")) 
    stop("Object of class 'locus' required")
  
  xlab <- paste("Chromosome", loc$seqname, "(Mb)")
  
  # manually adjusting segments to be within limits
  p1_min <- p1_limits[1]
  p1_max <- p1_limits[2]
  
  loc$TX <- loc$TX %>%
    dplyr::filter(gene_biotype == "protein_coding") %>%
    group_by(gene_id) %>%
    dplyr::mutate(
      start_clipped = max(start/1000000, p1_min), # in Mb
      end_clipped = min(end/1000000, p1_max)
    ) %>%
    ungroup()
  
  # Plot
  p = ggplot(loc$TX, aes(x = start_clipped, xend = end_clipped, y = logP, yend = logP,label = symbol)) +
    coord_cartesian(xlim = p1_limits ) + # Apply limits from gwas scatter while clipping lines
    geom_segment(size = 0.5, color = "grey20") +  # Draws the lines
    geom_text_repel() + 
    geom_hline(yintercept = -log10(0.05/20000),linetype = "dashed",col = "royalblue") +
    labs(x = xlab, y =expression(-log[10]~P)) +
    theme_classic()
  
  return(p)
  
}

### modified ldsc
# to return the SE of the covariance matrix I
modified_ldsc_rg = function (munged_sumstats, ancestry, sample_prev = NA, population_prev = NA, 
                             ld, wld, n_blocks = 200, chisq_max = NA, chr_filter = seq(1, 
                                                                                       22, 1)) 
{
  if (missing(ancestry)) {
    cli::cli_progress_step("No ancestry specified, checking for user-specified `ld` and `wld`")
    checkmate::assert_directory_exists(ld)
    checkmate::assert_directory_exists(wld)
  }
  else {
    checkmate::assert_choice(ancestry, c("AFR", "AMR", "CSA", 
                                         "EAS", "EUR", "MID"), null.ok = FALSE)
    cli::cli_progress_step("Using {ancestry} reference from Pan-UKB")
  }
  checkmate::assert_list(munged_sumstats)
  if (missing(sample_prev)) {
    cli::cli_alert_info("No sample prevalence data provided. Estimating heritabilities on the observed scale.")
  }
  checkmate::assert_number(n_blocks)
  n.blocks <- n_blocks
  n.traits <- length(munged_sumstats)
  n.V <- n.traits * (n.traits + 1)/2
  if (n.traits > 18) {
    n.blocks <- (((n.traits + 1) * (n.traits + 2))/2) + 
      1
    cli::cli_alert_info("Setting the number of blocks used to perform the block jacknife used to estimate the sampling covariance matrix (V) to {n.blocks}")
    if (n.blocks > 1000) {
      cli_alert_warning("The number of blocks needed to estimate V is > 1000, which may result in sampling dependencies across the blocks used to estimate standard errors and can bias results.")
    }
  }
  cov <- matrix(NA, nrow = n.traits, ncol = n.traits)
  V.hold <- matrix(NA, nrow = n.blocks, ncol = n.V)
  N.vec <- matrix(NA, nrow = 1, ncol = n.V)
  Liab.S <- rep(1, n.traits)
  I <- matrix(NA, nrow = n.traits, ncol = n.traits)
  h2_res <- tibble()
  cli::cli_progress_step("Reading LD Scores")
  x <- ldscr:::read_ld(ancestry, ld)
  x$CM <- x$MAF <- NULL
  cli::cli_progress_step("Reading weights")
  w <- ldscr:::read_wld(ancestry, wld)
  w$CM <- w$MAF <- NULL
  colnames(w)[ncol(w)] <- "wLD"
  cli::cli_progress_step("Reading M")
  m <- ldscr:::read_m(ancestry, ld)
  M.tot <- sum(m)
  m <- M.tot
  cli::cli_progress_step("Reading summary statistics")
  all_y <- purrr::imap(munged_sumstats, ~{
    if (is.character(.x)) {
      cli::cli_progress_step("Reading summary statistics for '{.y}' from {.x}")
      sumstats_df <- vroom::vroom(.x, col_types = vroom::cols())
    }
    else {
      cli::cli_progress_step("Reading summary statistics for '{.y}' from dataframe")
      sumstats_df <- .x
    }
    cli::cli_progress_step("Merging '{.y}' with LD-score files")
    merged <- ldscr:::merge_sumstats(sumstats_df, w, x, chr_filter)
    cli::cli_alert_info(glue::glue("{nrow(merged)}/{nrow(sumstats_df)} SNPs remain after merging '{.y}' with LD-score files"))
    if (is.na(chisq_max)) {
      chisq_max <- max(0.001 * max(merged$N), 80)
    }
    rm <- (merged$Z^2 > chisq_max)
    merged <- merged[!rm, ]
    cli::cli_alert_info(glue::glue("Removed {sum(rm)} SNPs with Chi^2 > {chisq_max} from '{.y}'; {nrow(merged)} SNPs remain"))
    return(merged)
  })
  s <- 1
  for (j in 1:n.traits) {
    y1 <- all_y[[j]]
    y1$chi1 <- y1$Z^2
    for (k in j:n.traits) {
      if (j == k) {
        trait <- names(munged_sumstats[j])
        cli::cli_progress_step("Estimating heritability for '{trait}'")
        samp.prev <- sample_prev[j]
        pop.prev <- population_prev[j]
        merged <- y1
        n.snps <- nrow(merged)
        merged$intercept <- 1
        merged$x.tot <- merged$L2
        merged$x.tot.intercept <- 1
        initial.w <- ldscr:::make_weights(chi1 = merged$chi1, 
                                          L2 = merged$L2, wLD = merged$wLD, N = merged$N, 
                                          M.tot)
        merged$weights <- initial.w/sum(initial.w)
        N.bar <- mean(merged$N)
        weighted.LD <- as.matrix(cbind(merged$L2, merged$intercept) * 
                                   merged$weights)
        weighted.chi <- as.matrix(merged$chi1 * merged$weights)
        analysis_res <- ldscr:::perform_analysis(n.blocks, n.snps, 
                                                 weighted.LD, weighted.chi, N.bar, m)
        V.hold[, s] <- analysis_res$pseudo.values
        N.vec[1, s] <- analysis_res$N.bar
        lambda.gc <- median(merged$chi1)/qchisq(0.5, 
                                                df = 1)
        mean.Chi <- mean(merged$chi1)
        ratio <- (analysis_res$intercept - 1)/(mean.Chi - 
                                                 1)
        ratio.se <- analysis_res$intercept.se/(mean.Chi - 
                                                 1)
        if (is.na(population_prev) == F & is.na(sample_prev) == 
            F) {
          h2_lia <- ldscr:::h2_liability(h2 = analysis_res$reg.tot, 
                                         sample_prev, population_prev)
          h2_res <- h2_res %>% bind_rows(tibble(trait = trait, 
                                                mean_chisq = mean.Chi, lambda_gc = lambda.gc, 
                                                intercept = analysis_res$intercept, intercept_se = analysis_res$intercept.se, 
                                                ratio = ratio, ratio_se = ratio.se, h2_observed = analysis_res$reg.tot, 
                                                h2_observed_se = analysis_res$tot.se, h2_Z = analysis_res$reg.tot/analysis_res$tot.se, 
                                                h2_p = 2 * pnorm(abs(h2_Z), lower.tail = FALSE), 
                                                h2_liability = h2_lia, h2_liability_se = h2_lia/h2_Z))
        }
        else {
          h2_res <- h2_res %>% bind_rows(tibble(trait = trait, 
                                                mean_chisq = mean.Chi, lambda_gc = lambda.gc, 
                                                intercept = analysis_res$intercept, intercept_se = analysis_res$intercept.se, 
                                                ratio = ratio, ratio_se = ratio.se, h2_observed = analysis_res$reg.tot, 
                                                h2_observed_se = analysis_res$tot.se, h2_Z = analysis_res$reg.tot/analysis_res$tot.se, 
                                                h2_p = 2 * pnorm(abs(h2_Z), lower.tail = FALSE)))
        }
        cov[j, j] <- analysis_res$reg.tot
        I[j, j] <- analysis_res$intercept
      }
      if (j != k) {
        trait1 <- names(munged_sumstats[j])
        trait2 <- names(munged_sumstats[k])
        cli::cli_progress_step("Estimating genetic covariance for for '{trait1}' and '{trait2}'")
        y2 <- all_y[[k]]
        y <- merge(y1, y2[, c("SNP", "N", "Z", "A1")], 
                   by = "SNP", sort = FALSE)
        y$Z.x <- ifelse(y$A1.y == y$A1.x, y$Z.x, -y$Z.x)
        y$ZZ <- y$Z.y * y$Z.x
        y$chi2 <- y$Z.y^2
        merged <- na.omit(y)
        n.snps <- nrow(merged)
        merged$intercept <- 1
        merged$x.tot <- merged$L2
        merged$x.tot.intercept <- 1
        initial.w <- ldscr:::make_weights(chi1 = merged$chi1, 
                                          L2 = merged$L2, wLD = merged$wLD, N = merged$N.x, 
                                          M.tot)
        initial.w2 <- ldscr:::make_weights(chi1 = merged$chi2, 
                                           L2 = merged$L2, wLD = merged$wLD, N = merged$N.y, 
                                           M.tot)
        merged$weights_cov <- (initial.w + initial.w2)/sum(initial.w + 
                                                             initial.w2)
        N.bar <- sqrt(mean(merged$N.x) * mean(merged$N.y))
        weighted.LD <- as.matrix(cbind(merged$L2, merged$intercept) * 
                                   merged$weights)
        weighted.chi <- as.matrix(merged$ZZ * merged$weights_cov)
        covariance_res <- ldscr:::perform_analysis(n.blocks, 
                                                   n.snps, weighted.LD, weighted.chi, N.bar, 
                                                   m)
        V.hold[, s] <- covariance_res$pseudo.values
        N.vec[1, s] <- covariance_res$N.bar
        cov[k, j] <- cov[j, k] <- covariance_res$reg.tot
        I[k, j] <- I[j, k] <- covariance_res$intercept
      }
      s <- s + 1
    }
  }
  v.out <- cov(V.hold)/crossprod(N.vec * (sqrt(n.blocks)/m))
  ratio <- tcrossprod(sqrt(Liab.S))
  S <- cov * ratio
  scaleO <- gdata::lowerTriangle(ratio, diag = TRUE)
  V <- v.out * tcrossprod(scaleO)
  colnames(S) <- if (is.null(names(munged_sumstats))) 
    paste0("V", 1:ncol(S))
  else names(munged_sumstats)
  rownames(S) <- if (is.null(names(munged_sumstats))) 
    paste0("V", 1:ncol(S))
  else names(munged_sumstats)
  if (mean(Liab.S) != 1) {
    r <- nrow(S)
    SE <- matrix(0, r, r)
    SE[lower.tri(SE, diag = TRUE)] <- sqrt(diag(V))
    colnames(SE) <- colnames(S)
    rownames(SE) <- rownames(S)
  }
  if (all(diag(S) > 0)) {
    ratio <- tcrossprod(1/sqrt(diag(S)))
    S_Stand <- S * ratio
    scaleO <- gdata::lowerTriangle(ratio, diag = TRUE)
    V_Stand <- V * tcrossprod(scaleO)
    r <- nrow(S)
    SE_Stand <- matrix(0, r, r)
    SE_Stand[lower.tri(SE_Stand, diag = TRUE)] <- sqrt(diag(V_Stand))
    colnames(SE_Stand) <- colnames(S)
    rownames(SE_Stand) <- rownames(S)
  }
  else {
    cli::cli_alert_warning("Your genetic covariance matrix includes traits estimated to have a negative heritability.")
  }
  ind <- which(lower.tri(S, diag = F), arr.ind = TRUE)
  rg_res <- tibble(trait1 = dimnames(S_Stand)[[2]][ind[, 2]], 
                   trait2 = dimnames(S_Stand)[[1]][ind[, 1]], rg = S_Stand[ind], 
                   rg_se = SE_Stand[ind], rg_p = 2 * pnorm(abs(rg/rg_se), 
                                                           lower.tail = FALSE))
  output <- list(h2 = h2_res, rg = rg_res, raw = list(V = V, 
                                                      S = S, I = I, N = N.vec, m = m, V_Stand = V_Stand, S_Stand = S_Stand, 
                                                      SE_Stand = SE_Stand),
                 covariance_full_result = covariance_res)
  class(output) <- c("ldscr_list", "list")
  return(output)
}

calculate_F_statistic_MR = function(pval,samplesize){
  F_stat = qf(pval, df1=1, df2=samplesize-1, lower.tail=FALSE)
  return(F_stat)
}

# miami plot GWAS vs TWMR eQTL-weighted gene p-vals
# color in significant TWMR means directionality for alpha (red = positive, blue = negative)
miami_TWMR_GWAS = function(gwas,twmr){
  
  manhattan = gwas %>%
    dplyr::mutate(CHR=as.numeric(stringr::str_split_fixed(string = variant_id,pattern = "_",n = 4)[,1]),
                  BP=as.numeric(stringr::str_split_fixed(string = variant_id,pattern = "_",n = 4)[,2])) %>%
    dplyr::rename(P=pval_nominal,SNP=variant_id) %>%
    dplyr::mutate(row_id = row_number()) %>% # Add a row id for tracking
    dplyr::group_by(P >0.05) %>%               # Group by the condition
    dplyr::mutate(
      keep_row = if_else(P >0.30, row_number() %% 30 == 1, # Keep every 30th row if P >0.30
                         if_else(P >0.05, row_number() %% 5 == 1,  # Keep every 5th row if P >0.05
                                 if_else(P >0.01, row_number() %% 2 == 1, TRUE))) # Keep every other row if P >0.01
    ) %>%
    dplyr::filter(keep_row) %>%              # Filter rows to keep
    dplyr::ungroup() %>%
    dplyr::select(SNP,P,CHR,BP)                
  
  genes_to_location = readr::read_csv("/lustre/scratch123/hgi/teams/trynka/resources/biomart/Homo_sapiens.GRCh38.111.genes.csv") %>%
    dplyr::rename(CHR = seqname,BP = start, gene = gene_name) %>%
    dplyr::right_join(twmr) %>%
    tidyr::drop_na() %>%
    dplyr::rename(P = p) %>%
    dplyr::select(gene,CHR,BP,P) %>%
    dplyr::mutate(CHR = as.numeric(CHR),
                  BP = as.numeric(BP)) %>%
    dplyr::mutate(logP_TWMR = -log10(P)) %>%
    dplyr::select(CHR,BP,logP_TWMR,gene) 
  
  manhattan = manhattan %>%
    dplyr::mutate(logP_GWAS = -log10(P)) %>%
    dplyr::select(CHR,BP,logP_GWAS) 
  
  
  # Custom function to format y-axis labels (remove "-" from negatives)
  
  format_y_labels = function(x) {
    abs(x)
  }
  
  # joining GWAS and TWMR tables so that the TWMR positions match the closest GWAS positions
  # for plotting in same region
  find_closest = function(value, reference_values) {
    reference_values[which.min(abs(reference_values - value))]
  }
  
  toplot = list()
  for(chr in 1:22){
    toplot[[chr]] = genes_to_location %>%
      dplyr::filter(CHR == chr) %>%
      dplyr::mutate(BP = sapply(BP, find_closest, 
                                reference_values = manhattan %>% 
                                  dplyr::filter(CHR == chr) %>%
                                  dplyr::pull(BP))) 
    
  }
  toplot = do.call(rbind,toplot) %>%
    dplyr::right_join(manhattan,by = join_by(CHR,BP)) %>%
    dplyr::mutate(logP_TWMR = case_when(logP_TWMR > 15 ~ 15, # fixing very small p-value to 15 - examine what's happening there
                                        .default = logP_TWMR)) %>%
    dplyr::arrange(CHR,BP) %>%
    dplyr::mutate(new_pos = 1:nrow(.),
                  col_GWAS = case_when(CHR %% 2 == 1 ~ "grey70",
                                       CHR %% 2 != 1 ~ "grey20"),
                  
                  col_TWMR = case_when(CHR %% 2 == 1 ~ "grey70",
                                       CHR %% 2 != 1 ~ "grey20"
                  )) %>%
    dplyr::mutate(col_GWAS = case_when(logP_GWAS > -log10(5*1e-8) ~ "darkred",
                                       .default = col_GWAS),
                  col_TWMR = case_when(logP_TWMR > -log10(0.05/8971)~ "darkred",
                                       .default = col_TWMR)) %>%
    dplyr::group_by(CHR) %>%
    dplyr::mutate(midpoints= mean(new_pos)) %>% # creating midpoints for label placement later
    dplyr::ungroup()
  # Create the Miami plot
  miami_plot = toplot %>% 
    
    ggplot() +
    # Plot Study 1 (above the x-axis)
    geom_point(aes(x = interaction(new_pos, CHR), y = logP_TWMR, color = col_TWMR),
               alpha = 0.6) +
    geom_text(
      aes(x = new_pos, y = logP_TWMR, label = ifelse(col_TWMR == "darkred", gene, NA)),
      size = 3, color = "grey10",
      vjust = -0.5, hjust = 0.5,
    ) +
    # Plot Study 2 (below the x-axis, as negative values)
    geom_point(aes(x = new_pos, y = -logP_GWAS, color = col_GWAS), alpha = 0.6) +
    # Facet by chromosome
    # facet_grid(~ CHR,  scales = "free_x",
    #            #nrow = 1,strip.position = "bottom",
    #             space='free') +
    # Customize colors
    scale_color_identity() +
    # Customize y-axis to show negative values but with positive labels
    scale_y_continuous(
      labels = format_y_labels,  # Use custom function to format labels
      breaks = seq(-ceiling(max(toplot$logP_TWMR,na.rm = TRUE)), ceiling(max(toplot$logP_TWMR,na.rm = TRUE)), by = 3),  # Adjust breaks as needed
      limits = c(-max(toplot$logP_TWMR,na.rm = TRUE) - 2, max(toplot$logP_TWMR,na.rm = TRUE) + 0.5)
    ) +
    # Add labels and title
    labs(
      x = "Chromosome",
      y = expression(-log[10]~P),
      title = paste(unique(twmr$phenotype),
                    unique(twmr$treatment),"TWMR (above) vs GWAS (below) p-values",sep = " "),
      color = "Study"
    ) +
    # Use a minimal theme
    theme_bw() +
    geom_hline(yintercept = -log10(0.05/20000),linetype = "dotted",col = "darkred") +
    geom_hline(yintercept = log10( 5*1e-8),linetype = "dotted",col = "darkred") +
    
    # Customize theme
    theme(
      strip.background = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.text.x = element_blank(),
      plot.margin = unit(c(1,1,3,1), units = "lines")
    )  +
    
    coord_cartesian(ylim = c(-16,16),clip = "off") +
    
    annotate(geom = "text", x = unique(toplot$midpoints),
             y = -16, 
             label = unique(toplot$CHR), size = 3) +
    
    # 0 and x border line
    geom_hline(yintercept = 0.0001,linetype = "solid",col = "grey10", linewidth = 1) +
    geom_hline(yintercept = -15.5,linetype = "solid",col = "grey10", linewidth = 0.7)
  
  # Display the plot
  return(miami_plot)
}
