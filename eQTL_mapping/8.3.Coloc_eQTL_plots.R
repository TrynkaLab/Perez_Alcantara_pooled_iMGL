# coloc eQTL plots


library(tidyverse)
library(ensembldb)
library(patchwork)
library(AnnotationHub)
source("functions.R")
options(future.globals.maxSize = 90000 * 1024^2) # 90Gb
options(ucscChromosomeNames=FALSE) # to allow for arbitrary chromosome identifiers.
readRenviron(path = "/nfs/users/nfs_m/ma23/.Renviron") # to use LDlinkR token
options(scipen = 999) # prevents from showing very small numbers as 0
dirs = list.dirs("../../data/results/8.colocalisation_analysis/coloc_results")[-1]

for(outDir in dirs){
  
  primary_eqtl =  c("macromap", "primary_eQTL", "QTD000559")
  
  if(any(str_detect(outDir, primary_eqtl))){
    type = "primary eQTL" # "primary eQTL" or "GWAS"
    
  }else{
    type = "GWAS"
  }
  scheme = Gviz::getScheme()
  scheme$GeneRegionTrack$fill <- "darkgrey"
  scheme$GeneRegionTrack$col <- NULL
  scheme$GeneRegionTrack$transcriptAnnotation <- "symbol"
  Gviz::addScheme(scheme, "myScheme")
  options(Gviz.scheme = "myScheme")
  
  # genome build GRCh38
  edb = EnsDb.Hsapiens.v86
  edb
  
  message("Working on ", outDir)
  
  ### load coloc results
  coloc_res = read_rds(paste0(outDir,"/coloc_results_single_causal_variant_",500000/1000,"kb.rds"))
  
  for(nms in names(coloc_res)){
    # sort in natural order
    coloc_res[[nms]] = coloc_res[[nms]][gtools::mixedsort(names(coloc_res[[nms]]))]
    
    pdf(file = paste0(outDir,"/",nms,"_coloc_plots.pdf"),width = 8,height = 8)
    
    for(locus in names(coloc_res[[nms]])){
      
      for(eGene in names(coloc_res[[nms]][[locus]])){
        
        message("Working on ", nms, " ", locus, " ",eGene)
        # color palette
        Pal <- colorRampPalette(c('white','blue'))
        
        toplot= coloc_res[[nms]][[locus]][[eGene]]$results %>%
          dplyr::mutate(GWAS =  - ( pnorm(-abs(z.df1),log.p=TRUE) + log(2) ) / log(10),
                        eQTL =  - ( pnorm(-abs(z.df2),log.p=TRUE) + log(2) ) / log(10)) %>%
          # double-check in 8.coloc... this works as intended - 
          # it does for GWAS, recovers correct nominal pvals
          # but seem inflated for the eQTL? investigate
          dplyr::select(snp,position,GWAS,eQTL,SNP.PP.H4) %>%
          tidyr::pivot_longer(cols=-c(snp,position,SNP.PP.H4),
                              names_to = "dataset",
                              values_to = "minus_logp") %>%
          dplyr::mutate(  color_pp4 = Pal(100)[ceiling(100*(SNP.PP.H4+0.000001))])
        
        # rename dataset values
        
        toplot = toplot %>%
          dplyr::mutate(dataset = case_when(dataset == "GWAS" ~ if_else(type == "GWAS","GWAS","primary eQTL"),
                                            .default = dataset))
        ##This adds a column of color values
        ## based on the y values
        top_coloc_var= unique(toplot[toplot$SNP.PP.H4 == max(toplot$SNP.PP.H4),"snp"]) %>%
          tidyr::separate_wider_delim(
            cols =  snp,
            names = c("chr", "pos", "minor_allele", "major_allele"),
            delim = "_",
            cols_remove = FALSE
          ) %>%
          dplyr::slice(1)
        
        manhattan_plot = ggplot(toplot,aes(x=position,y = minus_logp)) + 
          geom_point(shape = 21, color = "black", fill = toplot$color_pp4) +
          facet_wrap(nrow = 2,facets =vars(dataset), scales = "free_y") + 
          # indicating locus (GWAS), eGene (eQTL) and top coloc variant (might be different to lead GWAS or eQTL variant)
          ggtitle(paste0(locus, " locus , ", eGene, " eGene - chr",
                         top_coloc_var$chr, " ",top_coloc_var$pos, " ",
                         top_coloc_var$minor_allele," \u2192 ", # unicode arrow
                         top_coloc_var$major_allele ),
                  subtitle = paste0("PP for colocalization = ", 
                                    floor(coloc_res[[nms]][[locus]][[eGene]]$summary[[6]] * 100), "%")) + 
          ylab("-log10(p-value)")+
          xlab(paste0("chr",str_split(toplot$snp[1],pattern = "_")[[1]][1]))
        
        
        flt = AnnotationFilter(~ seq_name == top_coloc_var$chr &
                                 gene_start > min(toplot$position) &
                                 gene_end < max(toplot$position) &
                                 gene_biotype == "protein_coding" )
        Tx <- ensembldb::genes(edb, filter = flt)
        
        if(length(Tx$gene_id)!=0){
          eTrack <- Gviz::GeneRegionTrack(Tx,name = "Gene")
          
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
    }
    dev.off()
  }
}

