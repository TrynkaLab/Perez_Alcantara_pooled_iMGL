# process TWMR results
library(patchwork)
library(tidyverse)
library(locuszoomr)
library(AnnotationHub)
source("./functions.R")
library(ggrepel)

# readRenviron(path = "/nfs/users/nfs_m/ma23/.Renviron") # Use LDlinkR token



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

# TWMR paper PTWMR<3×10−6=0.05/16,000, number of genes tested per phenotype
# our unique genes 
twmr_p %>%
  tidyr::drop_na() %>%
  dplyr::group_by(phenotype,treatment) %>%
  dplyr::summarise(ngenes = n())
# max number of tests: 1430
# 5.5 x 10-6
twmr_p = twmr_p %>%
  as_tibble() %>%
  dplyr::group_by(phenotype,treatment) %>%
  dplyr::mutate(logP = -log10(p)) %>%
  dplyr::arrange(desc(logP)) %>%
  dplyr::mutate(p_Bonf = p.adjust(p,method = "bonferroni")) %>%
  dplyr::ungroup()

# number of significant results
twmr_p %>%
  dplyr::filter(p_Bonf<0.05) %>%
  dplyr::group_by(treatment,phenotype) %>%
  dplyr::summarise(n = n())

# treatment phenotype        n
# <chr>     <chr>        <int>
#   1 LPS       phagocytosis     2
# 2 untreated phagocytosis     2
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
