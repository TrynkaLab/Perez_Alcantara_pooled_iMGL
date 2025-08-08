# prepare lead variant file from coloc
# at the moment it's approximate - should do clumping based on LD to identify independent signals
library(tidyverse)
outdir = "../../../data/results/migration/2.check_association_results/lm_1pct_filtered_deflated/"
gwas = read_csv(paste0(outdir,"nominally_significant_results.csv"))
all_res = read_csv(paste0(outdir,"all_res_pqtl_migration.csv"))
selected = list()
for(treat in c("untreated","IFN","LPS")){
  
  selected[[treat]] = gwas %>% 
    dplyr::select(snp, ends_with(treat)) %>%
    dplyr::mutate(chr = as.numeric(str_split_i(snp,pattern = "_",i = 1)),
                  snp_pos = as.numeric(str_split_i(snp,pattern = "_",i = 2)),
                  REF = str_split_i(snp,pattern = "_",i = 3),
                  ALT = str_split_i(snp,pattern = "_",i = 4),
                  locus_name = snp) %>%
    dplyr::rename(coef = paste("coef",treat,sep = "_"),
                  p_value = paste("p",treat,sep = "_"),
                  se = paste("se",treat,sep = "_"),
                  variant_id = snp) %>%
    dplyr::select(variant_id,locus_name,coef,p_value,se,chr,snp_pos,REF,ALT) %>%
    dplyr::arrange(chr,snp_pos) 
  dir.create(paste0(outdir,treat),recursive = TRUE)
  # save nominal results
  all_res %>%
    dplyr::select(snp, ends_with(treat)) %>%
    dplyr::mutate(chr = as.numeric(str_split_i(snp,pattern = "_",i = 1)),
                  snp_pos = as.numeric(str_split_i(snp,pattern = "_",i = 2)),
                  REF = str_split_i(snp,pattern = "_",i = 3),
                  ALT = str_split_i(snp,pattern = "_",i = 4),
                  locus_name = snp) %>%
    dplyr::rename(coef = paste("coef",treat,sep = "_"),
                  p_value = paste("p",treat,sep = "_"),
                  se = paste("se",treat,sep = "_"),
                  variant_id = snp) %>%
    dplyr::select(variant_id,locus_name,coef,p_value,se,chr,snp_pos,REF,ALT) %>%
    dplyr::arrange(chr,snp_pos)  %>%
  write_tsv(.,paste0(outdir,treat,"/",treat,"_GWAS_GRCh38.tsv.gz"))
  
  selected[[treat]] = selected[[treat]] %>%
    dplyr::filter(p_value < 1e-5) # nominal significance threshold
  
  # create it manually so far
  # should have a genomic range solution to automate? 
  # probably easier to just do clumpiong with PLINK
  
  
  if(treat == "untreated"){
    selected[[treat]] = selected[[treat]] %>% # all further than 500kb away from each other
      dplyr::filter(variant_id %in% c("1_47874252_A_ATTTTT",
                               "8_117213529_C_T",
                               "14_36407696_A_G",
                               "3_344466_A_AACAC", 
                               "6_25031556_GAT_G")) 
    
    selected[[treat]] %>%
    write_tsv(.,paste0(outdir,treat,"/",treat,"_lead_vars_GRCh38_loci.txt"))
  }
  if(treat == "IFN"){
    selected[[treat]] = selected[[treat]] %>% # all further than 500kb away from each other
      dplyr::filter(variant_id %in% c("6_45594223_A_G",
                                      "8_71531011_T_G",
                                      "1_55291846_T_C",
                                      "16_63320491_C_T",
                                      "16_62181394_G_T",
                                      "2_225564230_A_G")) 
    
    selected[[treat]] %>%
      write_tsv(.,paste0(outdir,treat,"/",treat,"_lead_vars_GRCh38_loci.txt"))
  } 
  if(treat == "LPS"){
    selected[[treat]] = selected[[treat]]  %>% # all further than 500kb away from each other
      dplyr::filter(variant_id %in% c("3_149768114_T_TG",
                                      "18_22330647_T_G",
                                      "22_49026451_A_G",
                                      "22_37185382_G_C",
                                      "14_77901362_CT_C",
                                      "17_51625399_T_G")) 
    selected[[treat]] %>%
      write_tsv(.,paste0(outdir,treat,"/",treat,"_lead_vars_GRCh38_loci.txt"))
  } 
}

selected = list()
for(treat in c("untreated","IFN","LPS")){
  selected[[treat]] = read_tsv(paste0(outdir,treat,"/",treat,"_lead_vars_GRCh38_loci.txt"))
}
# find all genes at a distance from the locus

genes = read_csv("../../../../../../../teams/trynka/resources/biomart/Homo_sapiens.GRCh38.111.genes.csv")

genes = genes %>%
  dplyr::filter(gene_biotype == "protein_coding") %>%
  dplyr::select(seqname,start,end,strand,gene_name) %>%
  distinct()

# get gene within 500kb of lead
dist = 500*1000

genes_at_dist = list()

for(treat in c("untreated","IFN","LPS")){
  genes_at_dist[[treat]] = list()
  for(row in 1:nrow(selected[[treat]])){
    sel = selected[[treat]][row,]
    
    genes_at_dist[[treat]][[sel$locus_name]] = genes %>%
      dplyr::filter(seqname == sel$chr & start > (sel$snp_pos - dist) & end < (sel$snp_pos + dist)) %>%
      dplyr::mutate(locus = sel$locus_name,
                    tratment = treat)
 
  }
 
  genes_at_dist[[treat]] = do.call("rbind",genes_at_dist[[treat]])
}

genes_at_dist = do.call("rbind",genes_at_dist)


# sam's candidates
sam = read_tsv("../../../../CRISPR/OTAR2065_phagocytosis_CRISPR/data/2024_07_sam_screen_results.tsv") %>%
  dplyr::arrange(FDR) 

table(genes_at_dist$gene_name %in% sam$id)
# none