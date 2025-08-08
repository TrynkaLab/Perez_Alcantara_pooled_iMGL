# plot colocs from GWAS, eQTL and pQTL (our phenotype GWAS)
library(tidyverse)
library(GenomicRanges)
library(ensembldb)
library(Gviz)
library(EnsDb.Hsapiens.v86)

pqtl_dir = "../../../data/results/phagocytosis/2.check_association_results/lm_1pct_filtered_deflated/coloc_results/"
eqtl_dir = "../../../../OTAR2065_sc_eQTL/data/results/8.colocalisation_analysis/coloc_results/"
gwas = "GCST90027158"
GWAS_locus = "15_63277703_APH1B"
treatment = "LPS"
outDir = paste0("../../../data/results/phagocytosis/2.check_association_results/lm_1pct_filtered_deflated/coloc_results/",
                gwas)
pqtl_coloc = readRDS(paste0(pqtl_dir,gwas,"/coloc_results_single_causal_variant_500kb.rds"))
eqtl_coloc = readRDS(paste0(eqtl_dir,gwas,"/coloc_results_single_causal_variant_500kb.rds"))

pqtl_res = pqtl_coloc[[treatment]][[GWAS_locus]]
eqtl_res = eqtl_coloc[[paste("35",treatment,
                             "Not_proliferating",sep = "_")]][[GWAS_locus]]

scheme = Gviz::getScheme()
scheme$GeneRegionTrack$fill = "darkgrey"
scheme$GeneRegionTrack$col = NULL
scheme$GeneRegionTrack$transcriptAnnotation = "symbol"
Gviz::addScheme(scheme, "myScheme")
options(Gviz.scheme = "myScheme")

# genome build GRCh38
edb = EnsDb.Hsapiens.v86
edb

pdf(file = paste0(outDir,"/",treatment,"_",GWAS_locus,"_eqtl_GWAS_pqtl_coloc_plots.pdf"),width = 8,height = 8)

for(locus in names(eqtl_res)){
  
  message("Working on ", locus)
  # color palette
  Pal = colorRampPalette(c('white','blue'))
  
  eqtl_toplot= eqtl_res[[locus]]$results %>%
    dplyr::mutate(GWAS =  - ( pnorm(-abs(z.df1),log.p=TRUE) + log(2) ) / log(10),
                  eqtl =  - ( pnorm(-abs(z.df2),log.p=TRUE) + log(2) ) / log(10)) %>%
    # double-check in 8.coloc... this works as intended - 
    # it does for GWAS, recovers correct nominal pvals
    # but seem inflated for the eQTL? investigate
    dplyr::select(snp,position,GWAS,eqtl,SNP.PP.H4) %>%
    tidyr::pivot_longer(cols=-c(snp,position,SNP.PP.H4),
                        names_to = "dataset",
                        values_to = "minus_logp") %>%
    dplyr::mutate(  color_pp4 = Pal(100)[ceiling(100*(SNP.PP.H4+0.000001))])
  
  ### pqtl to plot
  
  pqtl_toplot= pqtl_res$results %>%
    dplyr::mutate(GWAS =  - ( pnorm(-abs(z.df1),log.p=TRUE) + log(2) ) / log(10),
                  pqtl =  - ( pnorm(-abs(z.df2),log.p=TRUE) + log(2) ) / log(10)) %>%
    # double-check in 8.coloc... this works as intended - 
    # it does for GWAS, recovers correct nominal pvals
    # but seem inflated for the eQTL? investigate
    dplyr::select(snp,position,GWAS,pqtl,SNP.PP.H4) %>%
    tidyr::pivot_longer(cols=-c(snp,position,SNP.PP.H4),
                        names_to = "dataset",
                        values_to = "minus_logp") %>%
    dplyr::mutate(  color_pp4 = Pal(100)[ceiling(100*(SNP.PP.H4+0.000001))])
  
  toplot = pqtl_toplot %>%
    dplyr::bind_rows(eqtl_toplot)
  
  ##This adds a column of color values
  ## based on the y values
  top_coloc_var_eqtl= unique(eqtl_toplot[eqtl_toplot$SNP.PP.H4 == max(eqtl_toplot$SNP.PP.H4),"snp"]) %>%
    tidyr::separate_wider_delim(
      cols =  snp,
      names = c("chr", "pos", "minor_allele", "major_allele"),
      delim = "_",
      cols_remove = FALSE
    ) %>%
    dplyr::slice(1)
  
  top_coloc_var_pqtl= unique(pqtl_toplot[pqtl_toplot$SNP.PP.H4 == max(pqtl_toplot$SNP.PP.H4),"snp"]) %>%
    tidyr::separate_wider_delim(
      cols =  snp,
      names = c("chr", "pos", "minor_allele", "major_allele"),
      delim = "_",
      cols_remove = FALSE
    ) %>%
    dplyr::slice(1)
  
  manhattan_plot = ggplot(toplot,aes(x=position,y = minus_logp)) + 
    geom_point(shape = 21, color = "black", fill = toplot$color_pp4) +
    facet_wrap(nrow = 3,facets =vars(dataset), scales = "free_y",
               labeller = labeller(dataset = 
                                     c("eqtl" = paste0("eQTL: ",locus, " locus, ",
                                                       top_coloc_var_eqtl$chr,
                                                       ":",top_coloc_var_eqtl$pos, " ",
                                                       top_coloc_var_eqtl$minor_allele," \u2192 ", # unicode arrow
                                                       top_coloc_var_eqtl$major_allele ,
                                                       " PP for colocalization = ", 
                                      floor(eqtl_res[[locus]]$summary[[6]] * 100), "%"),
                                       "GWAS" = "GWAS",
                                      "pqtl" = paste0("pQTL: ",
                                                      top_coloc_var_pqtl$chr,
                                                      ":",top_coloc_var_pqtl$pos, " ",
                                                      top_coloc_var_pqtl$minor_allele," \u2192 ", # unicode arrow
                                                      top_coloc_var_pqtl$major_allele ,
                                                      " PP for colocalization = ", 
                                                      floor(pqtl_res$summary[[6]] * 100), "%"))
               )) + 
    # indicating locus (GWAS), eGene (eQTL) and top coloc variant (might be different to lead GWAS or eQTL variant)
    ggtitle(paste0(locus, " locus ",GWAS_locus)) + 
    ylab("-log10(p-value)") +
    xlab(paste0("chr",str_split(eqtl_toplot$snp[1],pattern = "_")[[1]][1])) + 
    theme_minimal()

 
  
  
  flt = AnnotationFilter(~ seq_name == top_coloc_var_pqtl$chr &
                           gene_start > min(toplot$position) &
                           gene_end < max(toplot$position) &
                           gene_biotype == "protein_coding" )
  Tx = ensembldb::genes(edb, filter = flt)
  
  if(length(Tx$gene_id)!=0){
    eTrack = Gviz::GeneRegionTrack(Tx,name = "Gene")
    
    track = grid::grid.grabExpr(Gviz::plotTracks(eTrack,collapseTranscripts = TRUE, shape = "arrow"))
    
    # define layout
    lay = matrix(c(rep(1,7),2), nrow = 8, ncol = 1)
    
    
    gridExtra::grid.arrange(manhattan_plot, track, layout_matrix = lay)
  } else{
    
    
    text = grid.text("No genes in the vicinity")
    
    # define layout
    lay = matrix(c(rep(1,7),2), nrow = 8, ncol = 1)
    
    
    gridExtra::grid.arrange(manhattan_plot,text, layout_matrix = lay)
  }
  
  
}

dev.off()
