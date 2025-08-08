# process TWMR results
library(patchwork)
library(tidyverse)
library(locuszoomr)
library(AnnotationHub)
source("./functions.R")
library(ggrepel)

readRenviron(path = "/nfs/users/nfs_m/ma23/.Renviron") # Use LDlinkR token



########

twmr_res = list()
twmr_p = list()

for (phenotype in c("migration","phagocytosis")){
for(treat in c("untreated","LPS","IFN")){
  dir = paste0("../../data/results/8.5.eQTL_MR/TWMR/output/",treat,"/", phenotype,"/")
  
  
  files = list.files(dir,pattern = ".rds")

  for(gene in stringr::str_remove(files,".rds")){
    twmr_res[[paste(phenotype,treat,gene,sep = "_")]] = readr::read_rds(paste0(dir,gene,".rds"))
    
    if(anyNA(twmr_res[[paste(phenotype,treat,gene,sep = "_")]])){
      # if not tested for whatever reason, fill in with NA
      twmr_res[[paste(phenotype,treat,gene,sep = "_")]]$alpha = NA
      twmr_res[[paste(phenotype,treat,gene,sep = "_")]]$SE = NA
      twmr_res[[paste(phenotype,treat,gene,sep = "_")]]$P = NA
      twmr_res[[paste(phenotype,treat,gene,sep = "_")]]$Nsnps = NA
      twmr_res[[paste(phenotype,treat,gene,sep = "_")]]$Ngene = NA
      twmr_res[[paste(phenotype,treat,gene,sep = "_")]]$LD_matrix = NA
      twmr_res[[paste(phenotype,treat,gene,sep = "_")]]$eQTL_GWAS_effects = NA
      
    }
    twmr_res[[paste(phenotype,treat,gene,sep = "_")]]$gene = gene
    twmr_p[[paste(phenotype,treat,gene,sep = "_")]] = data.frame("p" = twmr_res[[paste(phenotype,treat,gene,sep = "_")]]$P,
                                                       "alpha" = twmr_res[[paste(phenotype,treat,gene,sep = "_")]]$alpha,
                                                       "alpha_SE" = twmr_res[[paste(phenotype,treat,gene,sep = "_")]]$SE,
                                                       "Nsnps" = twmr_res[[paste(phenotype,treat,gene,sep = "_")]]$Nsnps,
                                                       "Ngene" = twmr_res[[paste(phenotype,treat,gene,sep = "_")]]$Ngene,
                                "gene" = gene,
                                "treatment" = treat,
                                "phenotype" = phenotype)
  }
  
}
}
twmr_p = do.call("rbind",twmr_p) 

# TWMR paper PTWMR<3×10−6=0.05/16,000, number if genes tested per phenotype
# our unique genes 
twmr_p %>%
  tidyr::drop_na() %>%
  dplyr::group_by(phenotype,treatment) %>%
  dplyr::summarise(ngenes = n())
# max number of tests: 8971
# 5.5 x 10-6
twmr_p = twmr_p %>%
  as_tibble() %>%
  dplyr::group_by(phenotype,treatment) %>%
  dplyr::mutate(logP = -log10(p)) %>%
  dplyr::arrange(desc(logP)) %>%
  dplyr::mutate(p_Bonf = p.adjust(p,method = "bonferroni")) %>%
  dplyr::ungroup()


# bonferroni per gene, taking into account the number of snps tested 
readr::write_rds(twmr_res,
                 paste0("../../data/results/8.5.eQTL_MR/TWMR/output/phagocytosis_migration_twmr_res.rds"))

readr::write_csv(twmr_p,
                 paste0("../../data/results/8.5.eQTL_MR/TWMR/output/phagocytosis_migration_twmr_p.csv"))


twmr_p = readr::read_csv("../../data/results/8.5.eQTL_MR/TWMR/output/phagocytosis_migration_twmr_p.csv")
##### plotting specific gene regions ######


ah <- AnnotationHub()
ensDb_v106 <- ah[["AH100643"]] # GRCh38 genome build

treat = "LPS"
gene = "LRRK2"
gwas_file = "../../../OTAR2065_phenotypic_QTLs/data/results/phagocytosis/2.check_association_results/lm_1pct_filtered_deflated/all_res_pqtl_phagocytosis.csv"
GWAS = readr::read_csv(gwas_file) %>%
  dplyr::select(contains(treat),snp) %>%
  tidyr::separate_wider_delim(snp,"_",names = c("chrom", "pos","other_allele","effect_allele"),cols_remove = FALSE) %>%
  dplyr::rename(variant_id = snp,
                pval_nominal=paste("p",treat,sep = "_"),
                beta = paste("coef",treat,sep = "_")) %>%
  dplyr::mutate(chrom = as.numeric(chrom),
                pos = as.numeric(pos),
                pval_nominal = as.numeric(pval_nominal),
                beta = as.numeric(beta)) 

loc = locuszoomr::locus(as.data.frame(GWAS), gene = gene, 
                        ens_db = ensDb_v106, p = "pval_nominal", labs = "variant_id", pos = "pos",
                        flank = 1e5 ) 


p_gwas = locuszoomr::gg_scatter(loc,
                                showLD = FALSE,
                                scheme = c("grey","darkgrey","purple"),
                                labels = NULL,
                                nudge_y=1) 
# Extract x-axis limits from p1
p1_limits <- ggplot_build(p_gwas)$layout$panel_params[[1]]$x.range

loc$TX <-  loc$TX %>%
  dplyr::left_join(twmr_p %>% dplyr::filter(treatment == treat & phenotype == "phagocytosis"),by = join_by("symbol" == "gene")) %>%
  dplyr::mutate(logP = tidyr::replace_na(logP,0)) %>%
  dplyr::filter(gene_biotype == "protein_coding") %>%
  group_by(gene_id) %>%
  dplyr::mutate(
    start_clipped = max(start/1000000, p1_limits[1]), # in Mb
    end_clipped = min(end/1000000, p1_limits[2])
  ) %>%
  ungroup()


# Plot
p = ggplot() +
  geom_segment(data = loc$TX,
               aes(x = start_clipped, xend = end_clipped, y = logP, yend = logP),size = 0.5, color = "grey20",
               inherit.aes = FALSE) +  # Draws the lines
  geom_label_repel(data = loc$TX,
                  aes(x = end_clipped, y = logP,label =symbol), inherit.aes = FALSE) + 
  theme_classic() 

# adding plots together and coordinates for prettier results
p = p_gwas + 
  p$layers[[2]] + 
  geom_hline(yintercept = -log10(0.05/8971),linetype = "dashed",col = "royalblue") +
  geom_hline(yintercept = -log10(5*1e-8),linetype = "dashed",col = "darkred") +
  coord_cartesian(clip = "off" )+ # Apply limits from gwas scatter while clipping lines
  p$layers[[1]] + 
  scale_x_continuous(expand = expansion(add=0),
                     limits = c(min(loc$TX$start_clipped),p1_limits[2]))
  


p = p + patchwork::plot_annotation(title = paste0("Phagocytosis GWAS + TWMR: ", treat)) 
# 
# p2 = (p_gwas / pgg) + patchwork::plot_annotation(title = paste0("Phagocytosis GWAS: ", treatment)) + 
#   patchwork::plot_layout(ncol = 1,nrow = 2,heights = c(4,1))

outputDir = paste0("../../data/results/8.5.eQTL_MR/TWMR/output/")
plot(p)
ggsave(filename = paste0(outputDir,gene,"_",treat,"_phagocytosis_TWMR_plot.png"),
       device = "png",width = 7,height = 5.5)

# check how many genes are genome-wide sign
genome_wide_sign_threshold = -log10(0.05/0.05/8971) # genome-wide line 0.05/ ngenes (8971 tested max) ~ 5.5 x 10e-6 

# genes that modify the p-val from 0 in darkgrey / lightgrey depending on wether they 
# are uplifted from logP = 0

jeremy = read_delim("../../../resources/Jeremy_medrXiv_AD_candidate_genes.txt",col_names = TRUE)
sam = read_csv("../../../resources/Sam_NeuroID/NeuroID_L20vsT20_Forward_gdata.csv")

twmr_p %>%
  dplyr::filter(gene %in% unique(sam$id))
# p  alpha alpha_SE Nsnps Ngene gene   treatment phenotype     logP   p_Bonf
# <dbl>  <dbl>    <dbl> <dbl> <dbl> <chr>  <chr>     <chr>        <dbl>    <dbl>
#   1 4.30e-14  0.675   0.0894     2     1 LRRK2  LPS       phagocytosis 13.4  3.67e-10
# 2 9.18e- 9 -0.538   0.0936     1     1 MTHFSD LPS       phagocytosis  8.04 7.84e- 5
# 3 1.92e- 6 -1.22    0.256      2     2 SORL1  IFN       migration     5.72 1.65e- 2
# 4 8.96e- 6  0.332   0.0748     1     1 MTHFSD untreated migration     5.05 7.99e- 2
twmr_p %>%
  dplyr::filter(gene %in% unique(jeremy$symbol))
# p  alpha alpha_SE Nsnps Ngene gene    treatment phenotype     logP p_Bonf
# <dbl>  <dbl>    <dbl> <dbl> <dbl> <chr>   <chr>     <chr>        <dbl>  <dbl>
#   1 0.00000192 -1.22     0.256     2     2 SORL1   IFN       migration     5.72 0.0165
# 2 0.00617    -0.304    0.111     1     1 HS3ST1  untreated migration     2.21 1     
# do the same for genes a range of kbs away from AD loci

GWAS = list()
  
for (phenotype in c("migration","phagocytosis")){
  for(treat in c("untreated","LPS","IFN")){
    gwas_file = paste0("../../../OTAR2065_phenotypic_QTLs/data/results/",
    phenotype,
    "/2.check_association_results/lm_1pct_filtered_deflated/all_res_pqtl_",phenotype,".csv")
    
    GWAS[[paste(phenotype,treat,sep = "_")]] = readr::read_csv(gwas_file) %>%
      dplyr::select(contains(treat),snp) %>%
      tidyr::separate_wider_delim(snp,"_",names = c("chrom", "pos","other_allele","effect_allele"),cols_remove = FALSE) %>%
      dplyr::rename(variant_id = snp,
                    pval_nominal=paste("p",treat,sep = "_"),
                    beta = paste("coef",treat,sep = "_")) %>%
      dplyr::mutate(chrom = as.numeric(chrom),
                    pos = as.numeric(pos),
                    pval_nominal = as.numeric(pval_nominal),
                    beta = as.numeric(beta)) 
  }
  
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

outputDir = "../../data/results/8.5.eQTL_MR/TWMR/output/"
pdf(paste0(outputDir,"miami_plots_microglial_phenotypes.pdf"),
      width = 10,height = 5)
for (pheno in c("migration","phagocytosis")){
  for(treat in c("untreated","LPS","IFN")){
    twmr = twmr_p %>%
      dplyr::filter(treatment==treat & phenotype == pheno)
p = miami_TWMR_GWAS( GWAS[[paste(pheno,treat,sep = "_")]] ,
                       twmr)

plot(p)
  }
}

dev.off()
