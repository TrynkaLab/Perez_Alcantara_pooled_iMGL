# Annotation with seurat
.libPaths(c("/software/R-4.1.0/lib/R/library","/software/teamtrynka/ma23/R4.1/libs"))
library(patchwork)
library(tidyverse)
.libPaths(c("/software/teamtrynka/ma23/R4.1/libs","/software/R-4.1.0/lib/R/library"))
library(Seurat)
library(SCpubr)
source("./helpers.R")
options(future.globals.maxSize= 140097152000) #140 Gb

directory = "../../data/results/3.annotation_seurat/"
dir.create(directory, recursive = T)

## reference data -  Mancuso ex vivo (in mice) matured iPSC derived microglia human
ex_vivo_mancuso=readRDS("../../../resources/Mancuso_Nat_Neurosci_2019/ex_vivo_mancuso_filtered_annotated_allSamples_renamed_clusters.rds")

# Reduce automated annotation to microglia or perivascular macrophages
ex_vivo_mancuso$reduced_renamed_clusters = "Microglia"
ex_vivo_mancuso@meta.data[ex_vivo_mancuso@meta.data$seurat_clusters %in% c(1,4),"reduced_renamed_clusters"] = "Perivascular_macrophages"
ex_vivo_mancuso = subset(x = ex_vivo_mancuso, idents = c("1", "0"))

markers = FindMarkers(ex_vivo_mancuso, ident.1 = "1",ident.2="0", 
                      only.pos = FALSE, logfc.threshold = 0.1)
print(head(markers)) #markers in right direction
# loop through the samples and annotate, then save annotations

input_merged = readRDS("../../data/results/1.QC/filtered_media_merged_noRegression.rds")

# input_merged = readRDS("../../data/results/1.QC/filtered_media_merged_noRegression_50k.rds") # for tests
# remove pool3.h3_diff2_LPS from test object, which throws errors later
# input_merged = subset(input_merged, subset = orig.ident != "pool3.h3_diff2_LPS")

seurat_list = SplitObject(input_merged, split.by = "orig.ident")
rm(input_merged)
gc()




for(n in names(seurat_list)){
  # for(n in c("pool3.h3_diff2_IFNg",      "pool3.h3_diff2_LPS" ,      "pool3.h3_diff2_untreated",
  #      "pool4.h4_diff2_IFNg",  "pool4.h4_diff2_LPS",   "pool4.h4_diff2_untreated")){
  #   
  
  message("Working on sample ",n)
  seurat_list[[n]] = RenameCells(object = seurat_list[[n]], 
                                 new.names = paste(seurat_list[[n]]@meta.data$orig.ident,
                                                   seurat_list[[n]]@meta.data$orig_cell,
                                                   sep = "_"))
  
  # cell cycle scoring
  seurat_list[[n]]  = NormalizeData(seurat_list[[n]] )
  seurat_list[[n]]  = FindVariableFeatures(seurat_list[[n]] , verbose = FALSE)
  seurat_list[[n]]  = ScaleData(seurat_list[[n]] , verbose = FALSE)
  seurat_list[[n]]  = CellCycleScoring(seurat_list[[n]] ,  
                                       s.features=cc.genes.updated.2019$s.genes, 
                                       g2m.features=cc.genes.updated.2019$g2m.genes)
  seurat_list[[n]]$CC.Difference = seurat_list[[n]]$S.Score - seurat_list[[n]]$G2M.Score 
  
  # SCTransform
  seurat_list[[n]] = SCTransform(
    seurat_list[[n]],
    method = "glmGamPoi", # faster
    verbose = TRUE,
    vars.to.regress = c("CC.Difference","percent.mt") # group G2M and S cells 
  )
  seurat_list[[n]] = calculate_PCA_UMAP_neighbors_clusters_merged(seurat_list[[n]])
  
  query.anchors = FindTransferAnchors(reference = ex_vivo_mancuso, 
                                      query = seurat_list[[n]],
                                      dims = 1:30, reference.reduction = "pca",
                                      normalization.method = "SCT")
  message("Running TransferData()")
  predictions = TransferData(anchorset = query.anchors, 
                             refdata = ex_vivo_mancuso$reduced_renamed_clusters,
                             dims = 1:30)
  seurat_list[[n]]  = AddMetaData(seurat_list[[n]], metadata = predictions)
  
  message("Plotting UMAP")
  
  p = UMAPS_per_lib(seurat_list[[n]])
  
  
  png(paste0(directory,n,"_seurat_prediction_UMAP.png"),
      width = 10, height = 15, units = "in", res = 400)
  plot(p )
  dev.off()
  
  p1 = VlnPlot( seurat_list[[n]],
                features = c("CD163","COLEC12","LYVE1","F13A1"),
                assay = "SCT",
                slot = "data",
                group.by = "predicted.id") 
  
  p2 = VlnPlot(seurat_list[[n]],
               features = c("CX3CR1","TREM2","P2RY12","OLFML3",
                            #"CXCL8",# not found in mancuso diff expr
                            "TGFB1","GPR34","IL1B","CSF1R","GAS6",
                            #"PROS1", # not found in mancuso diff expr
                            "C1QA"),
               assay = "SCT",
               slot = "data",   
               group.by = "predicted.id") 
  
  png(paste0(directory,n,"_seurat_prediction_pvM_violin.png"),
      width = 10, height = 15, units = "in", res = 400)
  plot(p1)
  dev.off()
  
  png(paste0(directory,n,"_seurat_prediction_microglia_violin.png"),
      width = 10, height = 15, units = "in", res = 400)
  plot(p2)
  dev.off()
  
  p = SCpubr::do_FeaturePlot(seurat_list[[n]],
                             features = c("CD163","COLEC12","LYVE1","F13A1",
                                          "CX3CR1","TREM2","P2RY12","OLFML3",
                                          #"CXCL8",# not found in mancuso diff expr
                                          "TGFB1","GPR34","IL1B","CSF1R","GAS6",
                                          #"PROS1", # not found in mancuso diff expr
                                          "C1QA"))
  png(paste0(directory,n,"_seurat_pvM_microglia_UMAP.png"),
      width = 10, height = 15, units = "in", res = 400)
  plot(p )
  dev.off()
  
}

merged = merge(seurat_list[[1]], 
               seurat_list[2:length(seurat_list)])
rm(seurat_list)
gc()
saveRDS(merged,paste0(directory,"filtered_media_merged_CCdiff_mitPercent_regressed_annotated.rds")) 
# saveRDS(merged,"../../data/results/3.annotation_seurat/filtered_media_merged_CCdiff_mitPercent_regressed_50k_annotated.rds") # for tests
