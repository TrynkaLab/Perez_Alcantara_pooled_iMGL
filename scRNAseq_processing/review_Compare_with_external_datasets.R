# compare our data to external datasets
# data not available on GEO was downloaded from stemformatics
library(Seurat)
library(future)
library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(biomaRt)
library(limma)
library(edgeR)
library(ggrepel)
library(ggpubr)
options(stringsAsFactors = FALSE)

outdir = "../../data/results/1.5.compare_with_external_datasets"
dir.create(file.path(outdir))


source("./helpers.R")
# set this option when analyzing large datasets
options(future.globals.maxSize = 110000 * 1024^2) # 110Gb
treatment_cols =  c(untreated = "#8D918B", IFN = "#3A5683", LPS = "#F8766D")

# for new seurat assay slots
options(Seurat.object.assay.version = "v5")

directory = "../../data/results/1.4.markers"

# read in annotated seurat objects
seu_subset = list()
for(treatment in c("untreated", "IFN", "LPS")) {
  
  seu_subset[[treatment]] = LoadSeuratRds(
    file = paste0(directory,"/",treatment,"_filtered_harmony/",treatment,"_filtered_harmony.Rds"))
  
  
}


## functions ### ######
#to ensembl gene ids
geneName_to_ensembl <- function(expr_df) {
  
  library(biomaRt)
  
  ensembl38 <- useMart(
    "ensembl",
    dataset = "hsapiens_gene_ensembl"
  )
  
  # map gene symbols → Ensembl IDs
  gene_map <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters = "external_gene_name",
    values = rownames(expr_df),
    mart = ensembl38,
    useCache = FALSE
  )
  
  # remove empty mappings
  gene_map <- gene_map[gene_map$ensembl_gene_id != "", ]
  
  # keep only the first mapping per gene symbol
  gene_map <- gene_map[!duplicated(gene_map$external_gene_name), ]
  
  # map symbols → Ensembl
  ensembl_ids <- gene_map$ensembl_gene_id[
    match(rownames(expr_df), gene_map$external_gene_name)
  ]
  
  # drop unmapped genes
  keep <- !is.na(ensembl_ids)
  expr_df <- expr_df[keep, , drop = FALSE]
  ensembl_ids <- ensembl_ids[keep]
  
  # assign Ensembl IDs (temporarily)
  rownames(expr_df) <- ensembl_ids
  
  # collapse duplicate Ensembl IDs by summing
  expr_df_collapsed <- aggregate(
    expr_df,
    by = list(ensembl_id = rownames(expr_df)),
    FUN = sum
  )
  
  # set final row names
  rownames(expr_df_collapsed) <- expr_df_collapsed$ensembl_id
  expr_df_collapsed$ensembl_id <- NULL
  
  return(expr_df_collapsed)
}

# 
# ####main ######
# # https://www.stemformatics.org/datasets/view?id=7268
# Abud = read.table("../../../resources/Abud_2017_microglia_bulk/stemformatics_dataset_7268.raw.tsv")
# #metadata from https://www.stemformatics.org/datasets/view?id=7268
# colnames(Abud) = c("Abud_iMG_1","Abud_iMG_2","Abud_iMG_3","Abud_iMG_4","Abud_iMG_5","Abud_iMG_6",
#                    "Abud_iHPC_1","Abud_iHPC_2","Abud_iHPC_3", # hematopoietic precursor
#                    "Abud_iPSC_1","Abud_iPSC_2","Abud_iPSC_3","Abud_iPSC_4",
#                    "Abud_monocyte_CD14_1","Abud_monocyte_CD14_2","Abud_monocyte_CD14_3","Abud_monocyte_CD14_4","Abud_monocyte_CD14_5",
#                    "Abud_monocyte_CD14_CD16_1","Abud_monocyte_CD14_CD16_2","Abud_monocyte_CD14_CD16_3","Abud_monocyte_CD14_CD16_4",
#                    "Abud_pFMGL_1","Abud_pFMGL_2","Abud_pFMGL_3", # foetal primary microglia
#                    "Abud_pAMGL_1","Abud_pAMGL_2","Abud_pAMGL_3", # adult primary microglia
#                    "Abud_DC_1","Abud_DC_2","Abud_DC_3", # dendritic cells
#                    "Abud_iMG_neuron_1", "Abud_iMG_neuron_2","Abud_iMG_neuron_3","Abud_iMG_neuron_4","Abud_iMG_neuron_5", "Abud_iMG_neuron_6", # iPSC-derived microglia cocultured with neuron
#                    "Abud_iMG_noCD200_noCX3CL1_1", "Abud_iMG_noCD200_noCX3CL1_2", "Abud_iMG_noCD200_noCX3CL1_3",  # microglial cell differentiated from induced pluripotent stem cell with CD200 and CX3CL1 withdrawl
#                    "Abud_iMG_noTGFb_1","Abud_iMG_noTGFb_2","Abud_iMG_noTGFb_3") # microglial cell differentiated from induced pluripotent stem cell with TGFβ withdrawal
# 
# # ex vivo FACS-isolated microglia from fresh postmortem (Galatro_2017) or surgery-resected human brain (Gosselin_2017)
# # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99074
# Galatro = read.table("../../../resources/Galatro_2017/GSE99074_HumanMicrogliaBrainCounts.txt") # Grubman says these are raw counts but it's not true (there are decimal points)?
# 
# colnames(Galatro) = c( paste0("Galatro_exMGL_",1:39),
#                        paste0("Galatro_exBrain_",1:16),
#                        "S969" ,     "S947"  ,    "S861"   ,   "S810"  ,    "S805"  ,    "S755"    ,  "S291"  ,    "S262"    ,  "S243"    ,  "S147"   )
# 
# # Gosselin data is mouse, so not included
# 
# # Douvaras 2017
# #https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97744 for metadata
# # data from https://www.stemformatics.org/datasets/view?id=7193
# Douvaras = read.table("../../../resources/Douvaras_2017/stemformatics_dataset_7193.raw.tsv")
# 
# colnames(Douvaras) = c(paste0("Douvaras_pMAC_",1:8), # primary macrophages
#                        paste0("Douvaras_pAMGL_serum_",1:2), # primary microglia + serum
#                        paste0("Douvaras_iMG_",1:10), # iPSC-derived microglia
#                        paste0("Douvaras_pAMGL_",1:2), # primary microglia no serum
#                        paste0("Douvaras_PMAC_hepatic_",1:2) # primary hepatic macrophages 
# )
# 
# 
# # Muffat 2016
# # https://www.stemformatics.org/datasets/view?id=7062
# Muffat = read.table("../../../resources/Muffat_2016/stemformatics_dataset_7062.raw.tsv")
# colnames(Muffat) = c("Muffat_pFMGL_1","Muffat_pFMGL_2","Muffat_pFMGL_3",
#                      paste0("Muffat_iMG_",1:5),
#                      "Muffat_eMG",
#                      paste0("Muffat_iMG_neuron_",1:3), # iPSC-microglia cultured in neural conditioned medium
#                      paste0("Muffat_iNPC_",1:4) # iPSC-derived neuronal progenitor
# )
# 
# 
# # Mancuso data ex vivo - 2019
# # ex_vivo_mancuso_2019=readRDS("../../../resources/Mancuso_Nat_Neurosci_2019/ex_vivo_mancuso_filtered_annotated_allSamples_renamed_clusters.rds")
# # 
# # mg.sce = as.SingleCellExperiment(ex_vivo_mancuso)
# # 
# # # pseudobulk
# # ex_vivo_mancuso = aggregateAcrossCells(mg.sce, use_exprs_values = "counts",
# #                                        ids=colData(mg.sce)[,c( "renamed_clusters")]) 
# # 
# # ex_vivo_mancuso@colData
# # head(counts(ex_vivo_mancuso))
# # 
# # 
# # 
# # ex_vivo_mancuso = geneName_to_ensembl(counts(ex_vivo_mancuso))
# 
# ## Mancuso 2024
mancuso_2024 = readRDS("../../../resources/Mancuso_2024/data/Mancuso2024_filtered_with_pvM_renamed.rds")
# mg.sce = as.SingleCellExperiment(mancuso_2024)
# 
# # pseudobulk
# mancuso_2024 = aggregateAcrossCells(mg.sce, use_exprs_values = "counts",
#                                     ids=colData(mg.sce)[,c( "cell.state_v1")]) 
# 
# mancuso_2024@colData
# head(counts(mancuso_2024))
# 
# mancuso_2024 = geneName_to_ensembl(counts(mancuso_2024))
# 
# ## data used in Dolan et al 2023 fig 2 - they processed it from Gazestani et al 2023 ###
# # I annotated the samples
# Gazestani = readRDS("../../../resources/Gazestani_2023_early_late_AD/processed_by_Dolan/Gazestani_Dolan_annotated_seurat_2024.rds")
# # Vector of cell types to exclude
# exclude_types <- c("Myeloid_CD300E", "Macrophage_CD200R1")
# # Subset Seurat object to exclude these
# Gazestani <- subset(Gazestani, subset = !cell_type %in% exclude_types)
# # exclude non-microglia / pvM
# mg.sce = as.SingleCellExperiment(Gazestani)
# 
# # pseudobulk
# Gazestani_by_cell_type = aggregateAcrossCells(mg.sce, use_exprs_values = "counts",
#                                     ids=colData(mg.sce)[,c( "cell_type")]) 
# Gazestani_by_disease = aggregateAcrossCells(mg.sce, use_exprs_values = "counts",
#                                               ids=colData(mg.sce)[,c( "Status")]) 
# Gazestani_by_cell_type@colData
# Gazestani_by_disease@colData
# 
# head(counts(Gazestani_by_cell_type))
# Gazestani_by_cell_type = geneName_to_ensembl(counts(Gazestani_by_cell_type))
# Gazestani_by_disease = geneName_to_ensembl(counts(Gazestani_by_disease))
# 
# # our data - make pseudobulks per annotated cluster #### 
# 
# untreated_pseudo = AggregateExpression(seu_subset$untreated, 
#                                        assays = "RNA",
#                                        return.seurat = FALSE,
#                                        group.by = c("merged_clusters"))
# untreated_pseudo = as.data.frame(untreated_pseudo)
# colnames(untreated_pseudo) = str_remove_all(colnames(untreated_pseudo),"RNA.")
# 
# LPS_pseudo = AggregateExpression(seu_subset$LPS, 
#                                  assays = "RNA",
#                                  return.seurat = FALSE,
#                                  group.by = c("merged_clusters"))
# LPS_pseudo = as.data.frame(LPS_pseudo)
# colnames(LPS_pseudo) = str_remove_all(colnames(LPS_pseudo),"RNA.")
# 
# IFNg_pseudo = AggregateExpression(seu_subset$IFN, 
#                                   assays = "RNA",
#                                   return.seurat = FALSE,
#                                   group.by = c("merged_clusters"))
# IFNg_pseudo = as.data.frame(IFNg_pseudo)
# colnames(IFNg_pseudo) = str_remove_all(colnames(IFNg_pseudo),"RNA.")
# 
# untreated_pseudo = geneName_to_ensembl(untreated_pseudo)
# LPS_pseudo = geneName_to_ensembl(LPS_pseudo)
# IFNg_pseudo = geneName_to_ensembl(IFNg_pseudo)
# 
# ### saving pseudobulks just in case #####
# write.csv(untreated_pseudo,"../../data/results/1.5.compare_with_external_datasets/untreated.csv")
# write.csv(LPS_pseudo,"../../data/results/1.5.compare_with_external_datasets/LPS.csv")
# write.csv(IFNg_pseudo,"../../data/results/1.5.compare_with_external_datasets/IFNg.csv")
# 
# write.csv(Muffat,"../../data/results/1.5.compare_with_external_datasets/Muffat.csv")
# write.csv(mancuso_2024,"../../data/results/1.5.compare_with_external_datasets/mancuso_2024.csv")
# write.csv(Gazestani_by_cell_type,"../../data/results/1.5.compare_with_external_datasets/Gazestani_by_cell_type.csv")
# write.csv(Gazestani_by_disease,"../../data/results/1.5.compare_with_external_datasets/Gazestani_by_disease.csv")
# write.csv(Abud,"../../data/results/1.5.compare_with_external_datasets/Abud.csv")
# write.csv(Douvaras,"../../data/results/1.5.compare_with_external_datasets/Douvaras.csv")
# 
# ###################### merging and processing with limma
# mancuso_2024 = read.csv("../../data/results/1.5.compare_with_external_datasets/mancuso_2024.csv",row.names = "X")
# untreated_pseudo = read.csv("../../data/results/1.5.compare_with_external_datasets/untreated.csv",row.names = "X")
# LPS_pseudo = read.csv("../../data/results/1.5.compare_with_external_datasets/LPS.csv",row.names = "X")
# IFNg_pseudo = read.csv("../../data/results/1.5.compare_with_external_datasets/IFNg.csv",row.names = "X")
# 
# message("Merging")
# 
# colnames(mancuso_2024) =  paste0(colnames(mancuso_2024),".mancuso_2024")
# colnames(untreated_pseudo) = paste0(colnames(untreated_pseudo),".current_study_untreated")
# colnames(LPS_pseudo) = paste0(colnames(LPS_pseudo),".current_study_LPS")
# colnames(IFNg_pseudo) =paste0(colnames(IFNg_pseudo),".current_study_IFNg")
# 
# combined_commongenes = merge(mancuso_2024, untreated_pseudo, by = 0, all = TRUE)  # merge by row names (by=0 or by="row.names")
# rownames(combined_commongenes) = combined_commongenes$Row.names  # give row names again
# combined_commongenes = merge(combined_commongenes, LPS_pseudo, by = 0, all = TRUE)  # merge by row names (by=0 or by="row.names")
# rownames(combined_commongenes) = combined_commongenes$Row.names  # give row names again
# combined_commongenes = merge(combined_commongenes, IFNg_pseudo, by = 0, all = TRUE)  # merge by row names (by=0 or by="row.names")
# rownames(combined_commongenes) = combined_commongenes$Row.names  # give row names again
# 
# 
# combined_commongenes = combined_commongenes[,c(4:ncol(combined_commongenes))] # remove first column
# combined_commongenes = na.omit(combined_commongenes)   # remove rows that contain NA values (reduces table to size of smallest table)
# 
# message(paste0("There are ",nrow(combined_commongenes)," genes remaining in the combined dataframe")) # 17026
# 
# 
# combined_commongenes = DGEList(combined_commongenes, 
#                                genes = rownames(combined_commongenes),
#                                samples = colnames(combined_commongenes))
# keep = edgeR::filterByExpr(combined_commongenes)
# 
# combined_commongenes = combined_commongenes[keep,]
# nrow(combined_commongenes) # 11755
# 
# combined_commongenes = calcNormFactors(combined_commongenes)    # Calculate normalization factors. TMM by default
# 
# 
# samples = c(colnames(mancuso_2024), colnames(untreated_pseudo), 
#             colnames(LPS_pseudo),  colnames(IFNg_pseudo))   # create the design matrix
# study = c(rep("Mancuso_2024",ncol(mancuso_2024)),
#           rep("current_study_untreated",ncol(untreated_pseudo)),
#           rep("current_study_LPS",ncol(LPS_pseudo)),
#           rep("current_study_IFNg",ncol(IFNg_pseudo))
#           )
# 
# cell_type = c(str_remove(colnames(mancuso_2024),".mancuso_2024"),
#               str_remove(colnames(untreated_pseudo),".current_study_untreated"),
#               str_remove(colnames(LPS_pseudo),".current_study_LPS"),
#               str_remove(colnames(IFNg_pseudo),".current_study_IFNg")
#               
# )
# cell_type[cell_type == "Prolif"] = "Proliferating"
# cell_type[cell_type == "Interferon\n response (IRM)" ] = "IFN"
# 
# metadata = as.data.frame(cbind(samples,study,cell_type))
# metadata$samples = as.factor(metadata$samples)
# metadata$study = as.factor(metadata$study)
# metadata$cell_type = as.factor(metadata$cell_type)
# 
# 
# design = model.matrix(~study + cell_type)  
# 
# v = voom(combined_commongenes,design,plot = F) # voom normalize the read counts
# 
# v$E = removeBatchEffect(v$E,metadata$study)
# 
# plotPca <- function(expr_mat, metadata., exclude_samples = NULL) {
#   # Step 1: optionally remove samples
#   if(!is.null(exclude_samples)){
#     keep_samples <- setdiff(colnames(expr_mat), exclude_samples)
#   } else {
#     keep_samples <- colnames(expr_mat)
#   }
#   
#   # Step 2: match metadata
#   metadata_sub <- metadata.[match(keep_samples, metadata.$samples), ]
#   
#   # Sanity check
#   if(!all(keep_samples == metadata_sub$samples)){
#     stop("Column names of expr_mat and metadata$samples do not match after subsetting")
#   }
#   
#   # Step 3: compute PCA on full expr_mat
#   pca_res <- prcomp(t(expr_mat), retx = TRUE)
#   
#   # Step 4: extract PCA coordinates for the kept samples
#   pcs <- as.data.frame(pca_res$x[keep_samples, , drop = FALSE])
#   
#   # Step 5: add metadata
#   pcs <- cbind(
#     pcs,
#     Samples = metadata_sub$cell_type,
#     Study = metadata_sub$study,
#     original_name = metadata_sub$samples
#   )
#   
#   # Step 6: factor and colors
#   pcs$Study <- relevel(pcs$Study, ref = "current_study_untreated")
#   pcs$Samples <- factor(
#     pcs$Samples, ordered = TRUE,
#     levels = c("HM", "Disease\n associated (DAM)", "IFN", "DAM.IFN",
#                "CRM", "HLA","Interferon\n response (IRM)" , "Proliferating",
#                "pvM", "RM", "tCRM",
#                "ChoPlexEpi", "ChoPlexMHC")
#   )
#   
#   colors <- c(
#     "HM" = "#208AAE", "Disease\n associated (DAM)" = "#EBBAB9", "IFN" = "#FFDD4A", "DAM.IFN" = "#A44200",
#     "CRM" = "#000000", "HLA" = "#C4D88C", "Interferon\n response (IRM)"  = "#4A7C59", "Proliferating" = "red",
#     "pvM" = "darkred", "RM" = "cyan", "tCRM" = "grey20",
#     "ChoPlexEpi" = "grey40", "ChoPlexMHC" = "grey70"
#   )
#   
#   # Step 7: variance explained
#   percentVar <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)))
#   
#   # Step 8: plot
#   p <- ggplot(pcs, aes(x = PC1, y = PC2, color = Samples, fill = Samples,
#                        label = Samples, shape = Study)) +
#     geom_point(size = 5, alpha = 0.6, stroke = 1) +
#     geom_text_repel(max.overlaps = Inf, box.padding = 0.5) +
#     geom_point() +
#     scale_color_manual(values = colors) +
#     xlab(paste0("PC1: ", percentVar[1], "% variance")) +
#     ylab(paste0("PC2: ", percentVar[2], "% variance")) +
#     theme_pubr()
#   
#   return(p)
# }
# 
# p1 = plotPca(v$E, metadata. = metadata)
# p1
# 
# # zooming in 
# # Remove some samples while keeping PCs the same
# p1_zoom <- plotPca(
#   expr_mat = v$E,
#   metadata. = metadata,
#   exclude_samples = c("ChoPlexEpi.mancuso_2024","ChoPlexMHC.mancuso_2024",
#                       "Prolif.mancuso_2024","Proliferating.current_study_untreated",
#                       "Proliferating.current_study_LPS","Proliferating.current_study_IFNg",
#                       "HLA.current_study_IFNg")
# )
# p1_zoom
# 
# png("../../data/results/1.5.compare_with_external_datasets/PCA_external_datasets.png", 
#     width = 13, height = 9, units = "in", res = 400)
# plot(p1)
# dev.off()
# 
# png("../../data/results/1.5.compare_with_external_datasets/PCA_external_datasets_zoom.png", 
#     width = 13, height = 9, units = "in", res = 400)
# plot(p1_zoom)
# dev.off()
# 
# ### not very good clustering from pseudobulks
# # trying liger
# 
# library(rliger)
# library(Matrix)
# library(Seurat)
# options(ligerDotSize = 0.5)
# 
# mancuso_2024_sc = readRDS("../../../resources/Mancuso_2024/data/Mancuso2024_filtered_with_pvM_renamed.rds")
# 
# # https://welch-lab.github.io/liger/articles/liger_with_seurat.html#joint-clustering
# 
# # Your two Seurat objects
# obj1 = Seurat::UpdateSeuratObject(mancuso_2024_sc)
# obj2 = seu_subset$untreated
# 
# obj1@meta.data$origin = "mancuso_2024"
# obj2@meta.data$origin = "current_study_untreated"
# 
# obj1@meta.data$merged_clusters = obj1@meta.data$cell.state_v1
# obj1@meta.data$merged_clusters = as.character(obj1@meta.data$merged_clusters)
# obj1@meta.data$merged_clusters[obj1@meta.data$merged_clusters == "Prolif"] = "Proliferating"
# 
# obj1_sub <- subset(x = obj1, downsample = 10000) 
# obj2_sub <- subset(x = obj2, downsample = 10000) 
# 
# rm(obj1,obj2)
# gc()
# 
# mat1 <- GetAssayData(obj1_sub, layer = "counts")
# mat2 <- GetAssayData(obj2_sub, layer = "counts")
# 
# 
# seurat.list <- list(mancuso_2024 = mat1, current_study_untreated = mat2)
# ligerObj <- createLiger(seurat.list)
# gc()
# 
# merged_obj <- ligerObj %>%
#   normalize() %>%
#   selectGenes() %>%
#   scaleNotCenter()
# merged_obj

### not working well either

## with seurat only ######

# Load and prepare the Mancuso reference 
mancuso_2024_sc <- readRDS("../../../resources/Mancuso_2024/data/Mancuso2024_filtered_with_pvM_renamed.rds")
mancuso_2024_sc <- Seurat::UpdateSeuratObject(mancuso_2024_sc)
#mancuso <- subset(mancuso_2024_sc, downsample = 10000)
mancuso = mancuso_2024_sc
DefaultAssay(mancuso) <- "RNA"
mancuso$cell_type <- Idents(mancuso)

# Run UMAP on reference
mancuso <- RunUMAP(mancuso, dims = 1:30, 
                       reduction = "pca_integrated_SCT", return.model = TRUE)

# map a query Seurat object to Mancuso reference
  downsample_n = 80000
  query_mapped = list()
  for(condition in names(seu_subset)){
    mancuso_ref = mancuso
    query_seurat = seu_subset[[condition]]
  query_down <- subset(query_seurat, downsample = downsample_n)
  
  query_prepped <- query_down
  DefaultAssay(query_prepped) <- "RNA"
  query_prepped <- NormalizeData(query_prepped)
  
  shared.features <- intersect(
    rownames(mancuso_ref[["RNA"]]),
    rownames(query_prepped[["RNA"]])
  )
  
  anchors <- FindTransferAnchors(
    reference = mancuso_ref,
    query = query_prepped,
    dims = 1:30,
    k.anchor = 30,
    reference.assay = "RNA",
    normalization.method = "LogNormalize",
    reference.reduction = "pca_integrated_SCT",
    features = shared.features
  )
  
  predictions <- TransferData(
    anchorset = anchors,
    refdata = mancuso_ref$cell_type,
    dims = 1:30,
    weight.reduction = "pcaproject",
  )
  
  query_prepped <- AddMetaData(query_prepped, metadata = predictions)
  
  query_mapped[condition] <- MapQuery(
    anchorset = anchors,
    reference = mancuso_ref,
    query = query_prepped,
    refdata = list(cell_type = "cell_type"),
    reference.reduction = "pca_integrated_SCT",
    reduction.model = "umap"
  )
 
     #save_file <- paste0("../../data/results/1.5.compare_with_external_datasets/integrated_mancuso_2024_", condition, ".rds")
     
     #saveRDS(query_mapped, save_path) # can't be saved like this due to recursion issues
  }
  #query_mapped = readRDS(paste0("../../data/results/1.5.compare_with_external_datasets/integrated_mancuso_2024_", condition, ".rds"))

    
    
    for(condition in names(seu_subset)){
    ### confusion matrix for original cluster numbers
    conf_mat <- table(query_mapped[[condition]]$predicted.cell_type,
                      query_mapped[[condition]]$cluster_full)
    
    # Convert to data frame for ggplot (no melt)
    df_conf <- as.data.frame(conf_mat)
    colnames(df_conf) <- c("Predicted", "Cluster", "Count")
    
    df_conf <- df_conf %>%
      group_by(Cluster) %>%
      mutate(Percent = round(Count / sum(Count) * 100,2))
    
    
    p = ggplot(df_conf, aes(x = Cluster, y = Predicted, fill = Percent)) +
      geom_tile(color = "grey70") +
      geom_text(aes(label = Percent), color = "black") +
      scale_fill_gradient(low = "white", high = "red") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Confusion Matrix", x = "Cluster", y = "Mancuso Cell Type", 
           fill = "Percent")
    
    png( paste0("../../data/results/1.5.compare_with_external_datasets/integrated_mancuso_2024_", 
                condition, 
                "_confusion_original_clusters.png"),width = 5,height = 4,res=400,units = "in")
    plot(p)
    dev.off()
    

      if(condition == "untreated"){
        query_mapped[[condition]]@meta.data = query_mapped[[condition]]@meta.data %>%
          dplyr::mutate(cluster_full = as.numeric(cluster_full)) %>%
          # ordered by number of cells checked in the barplot below
          dplyr::mutate(merged_clusters = case_when(cluster_full %in% c(0,1,2,5,8,9) ~ "Disease\n associated (DAM)",
                                                    cluster_full %in% c(3,4) ~ "Proliferating",
                                                    cluster_full %in% c(6,7) ~ "CNS-associated\n macrophages (CAM)" ,
                                                    cluster_full == 10 ~ "Interferon\n response (IRM)" 
          ))  %>%
          dplyr::mutate(merged_clusters = factor(merged_clusters,
                                                 levels = c("Proliferating","CNS-associated\n macrophages (CAM)" ,
                                                            "Disease\n associated (DAM)","Interferon\n response (IRM)" ),
                                                 ordered = TRUE)) %>%
          as.data.frame()
        
        seu_subset[[condition]]@meta.data = seu_subset[[condition]]@meta.data %>%
          dplyr::mutate(cluster_full = as.numeric(cluster_full)) %>%
          # ordered by number of cells checked in the barplot below
          dplyr::mutate(merged_clusters = case_when(cluster_full %in% c(0,1,2,5,8,9) ~ "Disease\n associated (DAM)",
                                                    cluster_full %in% c(3,4) ~ "Proliferating",
                                                    cluster_full %in% c(6,7) ~ "CNS-associated\n macrophages (CAM)" ,
                                                    cluster_full == 10 ~ "Interferon\n response (IRM)" 
          ))  %>%
          dplyr::mutate(merged_clusters = factor(merged_clusters,
                                                 levels = c("Proliferating","CNS-associated\n macrophages (CAM)" ,
                                                            "Disease\n associated (DAM)",
                                                            "Interferon\n response (IRM)" ),
                                                 ordered = TRUE)) %>%
          as.data.frame()
      }
      if(condition == "IFN"){
        query_mapped[[condition]]@meta.data = query_mapped[[condition]]@meta.data %>%
          dplyr::mutate(cluster_full = as.numeric(cluster_full)) %>%
          dplyr::mutate(merged_clusters = case_when(cluster_full %in% c(1,3,5) ~ "Disease\n associated (DAM)",
                                                    cluster_full == 2   ~ "CNS-associated\n macrophages (CAM)" ,
                                                    cluster_full == 4   ~ "Proliferating",
                                                    cluster_full == 0 ~ "Interferon\n response (IRM)" ))  %>%
          dplyr::mutate(merged_clusters = factor(merged_clusters,
                                                 levels = c("Proliferating","Disease\n associated (DAM)",
                                                            "CNS-associated\n macrophages (CAM)",
                                                            "Interferon\n response (IRM)" ),
                                                 ordered = TRUE)) %>%
          as.data.frame()
        
        seu_subset[[condition]]@meta.data = seu_subset[[condition]]@meta.data %>%
          dplyr::mutate(cluster_full = as.numeric(cluster_full)) %>%
          dplyr::mutate(merged_clusters = case_when(cluster_full %in% c(1,3,5) ~ "Disease\n associated (DAM)",
                                                    cluster_full == 2   ~ "CNS-associated\n macrophages (CAM)" ,
                                                    cluster_full == 4   ~ "Proliferating",
                                                    cluster_full == 0 ~ "Interferon\n response (IRM)" ))  %>%
          dplyr::mutate(merged_clusters = factor(merged_clusters,
                                                 levels = c("Proliferating","Disease\n associated (DAM)",
                                                            "CNS-associated\n macrophages (CAM)",
                                                            "Interferon\n response (IRM)" ),
                                                 ordered = TRUE)) %>%
          as.data.frame()
      }
      if(condition == "LPS"){
        query_mapped[[condition]]@meta.data = query_mapped[[condition]]@meta.data %>%
          dplyr::mutate(cluster_full = as.numeric(cluster_full)) %>%
          dplyr::mutate(merged_clusters = case_when(cluster_full %in% c(0,1,6) ~ "Disease\n associated (DAM)",
                                                    cluster_full %in% c(2,3,8) ~ "CNS-associated\n macrophages (CAM)" ,
                                                    cluster_full %in% c(4,5,7) ~ "Proliferating"
          )) %>%
          dplyr::mutate(merged_clusters = factor(merged_clusters,
                                                 levels = c("Proliferating","Disease\n associated (DAM)",
                                                            "CNS-associated\n macrophages (CAM)" ),
                                                 ordered = TRUE)) %>%
          as.data.frame()
        
        seu_subset[[condition]]@meta.data = seu_subset[[condition]]@meta.data %>%
          dplyr::mutate(cluster_full = as.numeric(cluster_full)) %>%
          dplyr::mutate(merged_clusters = case_when(cluster_full %in% c(0,1,6) ~ "Disease\n associated (DAM)",
                                                    cluster_full %in% c(2,3,8) ~ "CNS-associated\n macrophages (CAM)" ,
                                                    cluster_full %in% c(4,5,7) ~ "Proliferating"
          )) %>%
          dplyr::mutate(merged_clusters = factor(merged_clusters,
                                                 levels = c("Proliferating","Disease\n associated (DAM)",
                                                            "CNS-associated\n macrophages (CAM)" ),
                                                 ordered = TRUE)) %>%
          as.data.frame()
        
      }
      Idents(seu_subset[[treatment]]) = seu_subset[[treatment]]$merged_clusters
    
    
    
    ### final annotated clusters
    conf_mat <- table(query_mapped[[condition]]$predicted.cell_type,
                      query_mapped[[condition]]$merged_clusters)
    
    # Convert to data frame for ggplot (no melt)
    df_conf <- as.data.frame(conf_mat)
    colnames(df_conf) <- c("Predicted", "Cluster", "Count")
    
    df_conf <- df_conf %>%
      group_by(Cluster) %>%
      mutate(Percent = round(Count / sum(Count) * 100,2))
    
    
    p = ggplot(df_conf, aes(x = Cluster, y = Predicted, fill = Percent)) +
      geom_tile(color = "grey70") +
      geom_text(aes(label = Percent), color = "black") +
      scale_fill_gradient(low = "white", high = "red") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Confusion Matrix", x = "Cluster", y = "Mancuso Cell Type", 
           fill = "Percent")
    
    png( paste0("../../data/results/1.5.compare_with_external_datasets/integrated_mancuso_2024_", 
                condition, 
                "_confusion.png"),width = 5,height = 4,res=400,units = "in")
    plot(p)
    dev.off()
    
    
    #### plotting integrated
    cats1 <- unique(seu_subset[[condition]]$merged_clusters)
    cats2 <- unique(query_mapped[[condition]]$predicted.cell_type)
    
    all_cats <- sort(unique(c(cats1, cats2)))
    # merged colors
    library(scales)
    
    categories = c( "Antigen-presenting\n response (HLA)", "CNS-associated\n macrophages (CAM)" ,
                    "Disease\n associated (DAM)",  "Homeostatic (HM)"   , "Interferon\n response (IRM)" ,       
                     "Proliferating",  "Transitioning CRM" ,"Translational\n response (TRM)" )
    
    colors <- hue_pal()(length(categories))
    
    names(colors) <- categories
    
    library(patchwork)
    
    p1 <- DimPlot(
      seu_subset[[condition]],
      reduction = "umap",
      group.by = "merged_clusters",
      order = rev(c(
        "Transitioning CRM",
        "Homeostatic (HM)",
        "Translational\n response (TRM)",
        "Antigen-presenting\n response (HLA)",
        "CNS-associated\n macrophages (CAM)",
        "Disease\n associated (DAM)",
        "Proliferating",
        "Interferon\n response (IRM)"
      )),
      cols = colors,
      #raster = TRUE,
      #pt.size = 1.5
    ) + ggtitle("") +
      theme(legend.position = "none")
    
    p2 <- DimPlot(
      query_mapped[[condition]],
      reduction = "umap",
      group.by = "predicted.cell_type",
      order = rev(c(
        "Transitioning CRM",
        "Homeostatic (HM)",
        "Translational\n response (TRM)",
        "Antigen-presenting\n response (HLA)",
        "Interferon\n response (IRM)",
        "CNS-associated\n macrophages (CAM)",
        "Disease\n associated (DAM)",
        "Proliferating"
      )),
      cols = colors,
      #raster = TRUE,
      #pt.size = 1.5
    ) + ggtitle("") + 
      theme(legend.position = "bottom")
    
    # Remove legend unless IFN
    if (condition != "IFN") {
      p2 <- p2 + theme(legend.position = "none")
    }
    
    # Combine plots
    combined <- p1 + p2 + patchwork::plot_layout(ncol = 1)
    

    
    # Save
    png(
      paste0("../../data/results/1.5.compare_with_external_datasets/integrated_mancuso_2024_", 
             condition, ".png"),
      width = 8,
      height = 9,
      res = 400,
      units = "in"
    )
    
    print(combined)
    dev.off()
}

