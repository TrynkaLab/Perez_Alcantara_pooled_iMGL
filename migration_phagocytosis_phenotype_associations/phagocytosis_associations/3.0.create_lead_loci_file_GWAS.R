# prepare lead variant file from coloc
# at the moment it's approximate - should do clumping based on LD to identify independent signals
library(tidyverse)
outdir = "../../../data/results/phagocytosis/2.check_association_results/lm_1pct_filtered_deflated/"
gwas = read_csv(paste0(outdir,"nominally_significant_results.csv"))
all_res = read_csv(paste0(outdir,"all_res_pqtl_phagocytosis.csv"))
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
  selected
  
  if(treat == "untreated"){
    dir.create(paste0(outdir,treat), recursive = TRUE)
    selected[[treat]] = selected[[treat]] %>% # all further than 500kb away from each other
      dplyr::filter(variant_id %in% c("2_1442973_C_T",
                               "2_14403297_G_C",
                               "2_28853835_C_CA",
                               "2_230589088_C_CAA", # more or less in the middle, still lowest p-value
                               "11_126116412_C_G")) 
    
    selected[[treat]] %>%
    write_tsv(.,paste0(outdir,treat,"/",treat,"_lead_vars_GRCh38_loci.txt"))
  }
  if(treat == "IFN"){
    dir.create(paste0(outdir,treat), recursive = TRUE)
    selected[[treat]] = selected[[treat]] %>% # all further than 500kb away from each other
      dplyr::filter(variant_id %in% c("2_48532706_C_T",
                                      "6_102537150_T_A",
                                      "11_128639769_G_GA",
                                      "12_53650704_G_A",
                                      "13_28388529_T_C")) 
    
    selected[[treat]] %>%
      write_tsv(.,paste0(outdir,treat,"/",treat,"_lead_vars_GRCh38_loci.txt"))
  } 
  if(treat == "LPS"){
    dir.create(paste0(outdir,treat), recursive = TRUE)
    selected[[treat]] = selected %>% # all further than 500kb away from each other
      dplyr::filter(variant_id %in% c("2_230591860_C_T",
                                      "3_64871716_C_CTA",
                                      "9_35638781_C_A",
                                      "11_62248820_C_CT",
                                      "11_132104609_G_GTAA",
                                      "13_28388529_T_C",
                                      "16_70165041_A_G")) 
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
