# inspect diffexpr

library(tidyverse)
library(ggrepel)
outdir="../../data/results/5.1.2.Diff_expr_limma_high_low_PRS"

sam_neuroid = read_csv("../../../resources/Sam_NeuroID/NeuroID_L20vsT20_Forward_gdata.csv")
res = list()
for(comparison in c("Q4vQ1","InteractionIFN","InteractionLPS")){
res[[comparison]] = read_tsv(paste0("../../data/results/5.1.2.Diff_expr_limma_high_low_PRS/DiffExpr_",
         comparison,
         "all_genes_negative_lower_in_latter.txt")) %>%
  dplyr::filter(adj.P.Val < 0.05) 

}

lapply(res,nrow)

table(res$Q4vQ1$symbol %in% sam_neuroid$id)
table(res$InteractionIFN$symbol %in% sam_neuroid$id)
table(res$InteractionLPS$symbol %in% sam_neuroid$id)

res$InteractionLPS$symbol[res$InteractionLPS$symbol %in% sam_neuroid$id]
# "PICALM"   "SIGLEC11" "CD55"     "GTPBP2"   "RAB29"   

res_pheno = list()
for(phenotype in c("phagocytosis","migration")){
for(treatment in c("untreated","IFN","LPS")){
  res_pheno[[phenotype]][[treatment]] = read_csv(paste0("../../data/results/3.LMM_expr_phenotype_PRS_interactions/all_",
                                      phenotype,
                                      "_", treatment,
                                      "_Not_proliferating_interactions_expr_x_pheno_PRS_noreps.csv")) %>%
    dplyr::filter(qvals < 0.1) 
  

}
}

table(res_pheno$phagocytosis$untreated$qvals<0.2)
table(res_pheno$phagocytosis$IFN$qvals<0.2)
table(res_pheno$phagocytosis$LPS$qvals<0.2)

table(res_pheno$migration$untreated$qvals<0.2)
table(res_pheno$migration$IFN$qvals<0.2)
table(res_pheno$migration$LPS$qvals<0.2)

table(res_pheno$phagocytosis$untreated$gene %in% sam_neuroid$id)
table(res_pheno$phagocytosis$IFN$gene %in% sam_neuroid$id)
table(res_pheno$phagocytosis$LPS$gene %in% sam_neuroid$id)

table(res_pheno$migration$untreated$gene %in% sam_neuroid$id)
table(res_pheno$migration$IFN$gene %in% sam_neuroid$id)
res_pheno$migration$IFN$gene[res_pheno$migration$IFN$gene %in% sam_neuroid$id]
# "APP"  "LRP1"
table(res_pheno$migration$LPS$gene %in% sam_neuroid$id)

