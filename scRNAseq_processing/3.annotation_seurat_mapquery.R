# Annotation with seurat - then MapQuery()
.libPaths(c("/software/R-4.1.0/lib/R/library","/software/teamtrynka/ma23/R4.1/libs"))
library(patchwork)
library(tidyverse)
.libPaths(c("/software/teamtrynka/ma23/R4.1/libs","/software/R-4.1.0/lib/R/library"))
library(Seurat)
library(SCpubr)
source("./helpers.R")
options(future.globals.maxSize= 290097152000) #290 Gb
# options(future.globals.maxSize= 14097152000) #14 Gb

directory = "../../data/results/3.annotation_seurat_mapquery/"
dir.create(directory, recursive = T)

###### start
# loop through the samples and annotate, then save annotations

input_merged = readRDS("../../data/results/1.QC/filtered_media_merged_noRegression.rds")

# input_merged = readRDS("../../data/results/1.QC/filtered_media_merged_noRegression_50k.rds") # for tests

# test - integrate some samples from pool 2, pool 4 and pool 8 and do automatic annotation on these
# then transfer the rest
seurat_list = SplitObject(input_merged, split.by = "orig.ident")
rm(input_merged)
gc()

# integrate and annotate with reference

# first normalise everything

for(n in names(seurat_list)){
  message("Working on sample ",n)
  seurat_list[[n]] = RenameCells(object = seurat_list[[n]], 
                                 new.names = paste(seurat_list[[n]]@meta.data$orig.ident,
                                                   seurat_list[[n]]@meta.data$orig_cell,
                                                   sep = "_"))
  
  # cell cycle scoring
  seurat_list[[n]]  = NormalizeData(seurat_list[[n]] )
  seurat_list[[n]]  = FindVariableFeatures(seurat_list[[n]] , verbose = TRUE)
  seurat_list[[n]]  = ScaleData(seurat_list[[n]] , verbose = TRUE)
  seurat_list[[n]]  = CellCycleScoring(seurat_list[[n]] ,  
                                       s.features=cc.genes.updated.2019$s.genes, 
                                       g2m.features=cc.genes.updated.2019$g2m.genes)
  seurat_list[[n]]$CC.Difference = seurat_list[[n]]$S.Score - seurat_list[[n]]$G2M.Score 
  
  # SCTransform
  seurat_list[[n]] = SCTransform(
    seurat_list[[n]],
    verbose = FALSE,
    method = "glmGamPoi", # faster
    vars.to.regress = "CC.Difference" # group G2M and S cells 
  )
  
  
}

save(seurat_list,file = paste0(directory,"checkpoint1.RData"))
load(paste0(directory,"checkpoint1.RData"))

# see if annotation works well with untreated only, then extend to others
seurat_list = seurat_list[grep("plst|untreated|AG|YC|ITMG",x = names(seurat_list),value = TRUE)]

reference.list = list()
for(sample in names(seurat_list)){
  ncells = ifelse(ncol(seurat_list[[sample]]) * 0.05 < 1000, 
                  yes = ncol(seurat_list[[sample]]), 
                  no = floor(0.05* ncol(seurat_list[[sample]])))
  reference.list[[sample]] = seurat_list[[sample]][, sample(colnames(seurat_list[[sample]]),
                                                            size = ncells, # sampling 5% of each sample, or 1000 cells
                                                            replace=F)]
}


rm(seurat_list)
gc()
# integrate subset
features = SelectIntegrationFeatures(object.list = reference.list, nfeatures = 3000)
reference.list = PrepSCTIntegration(object.list = reference.list, anchor.features = features)
for (i in 1:length(reference.list)) {
  reference.list[[i]] = RunPCA(object = reference.list[[i]], verbose = TRUE, features = features)
}

anchors = FindIntegrationAnchors(object.list = reference.list, normalization.method = "SCT",
                                 anchor.features = features)

rm(reference.list)
gc()
save(anchors,file = paste0(directory,"anchors.RData"))
load(file = paste0(directory,"anchors.RData"))

integrated = IntegrateData(anchorset = anchors, 
                           normalization.method = "SCT",
                           k.weight = 80) # have to decrease from 100 because of small dataset, see https://github.com/satijalab/seurat/issues/3930

DefaultAssay(integrated) = "integrated"
rm(anchors)
gc()
# Run the standard workflow for visualization and clustering
integrated = calculate_PCA_UMAP_neighbors_clusters_integrated(integrated)
save(integrated,file = paste0(directory,"integrated.RData"))
load(file = paste0(directory,"integrated.RData"))

# annotate

## reference data -  Mancuso ex vivo (in mice) matured iPSC derived microglia human
ex_vivo_mancuso=readRDS("../../../resources/Mancuso_Nat_Neurosci_2019/ex_vivo_mancuso_filtered_annotated_allSamples_renamed_clusters.rds")

# Reduce automated annotation to microglia or perivascular macrophages
ex_vivo_mancuso$reduced_renamed_clusters = "Microglia"
ex_vivo_mancuso@meta.data[ex_vivo_mancuso@meta.data$seurat_clusters %in% c(1,4),"reduced_renamed_clusters"] = "Perivascular_macrophages"

query.anchors = FindTransferAnchors(reference = ex_vivo_mancuso, 
                                    query = integrated,
                                    dims = 1:30, 
                                    reference.reduction = "pca",
                                    normalization.method = "SCT")
predictions = TransferData(anchorset = query.anchors, 
                           refdata = ex_vivo_mancuso$reduced_renamed_clusters,
                           dims = 1:30)
integrated = AddMetaData(integrated, metadata = predictions)


p = UMAPS_all_libs(integrated)

p1 = SCpubr::do_DimPlot(integrated,
                        label = F,  group.by = "predicted.id",
                        legend.position = "right")  +
  ggtitle("Seurat predictions")

p2 = SCpubr::do_FeaturePlot(integrated,
                            features = "prediction.score.Microglia"
) + 
  ggtitle("Microglia score")

png(paste0(directory,"media_integrated_UMAP_reference_subset.png"),
    width = 10, height = 18, units = "in", res = 400)
p / (p1 + p2)
dev.off()

save(integrated,file = paste0(directory,"integrated_annotated.RData"))
## check here and if umap not satisfactory annotate each treatment set separately

rm(ex_vivo_mancuso,p,p1,predictions,query.anchors)
load(paste0(directory,"integrated_annotated.RData"))
load(paste0(directory,"checkpoint1.RData"))

# map the query - remove cells that are already present in the integrated dataset
to_exclude = colnames(integrated)

query.list = list()
for(sample in names(seurat_list)){
  seurat_list[[sample]] =  subset( seurat_list[[sample]] ,
                                   subset = !colnames(seurat_list[[sample]]) %in% to_exclude)
  
}

query.list = seurat_list[setdiff(names(seurat_list),c("pool2.YC_10","pool2.AG_20","pool2.AG_40","pool2.AG_60", "pool2.ITMG","pool4.h4_diff2_IFNg","pool4.h4_diff2_LPS",
                                                      "pool4.h4_diff2_untreated"))]

rm(seurat_list)
gc()

for(n in names(query.list)){
  message("Working on sample ",n)
  
  query.list[[n]] = calculate_PCA_UMAP_neighbors_clusters_merged(query.list[[n]])
  
  query.anchors = FindTransferAnchors(reference = integrated, 
                                      query = query.list[[n]],
                                      dims = 1:30, reference.reduction = "pca",
                                      normalization.method = "SCT")
  predictions = TransferData(anchorset = query.anchors, 
                             refdata = integrated$predicted.id,
                             dims = 1:30)
  
  
  seurat_list[[n]]  = AddMetaData(seurat_list[[n]], metadata = predictions)
  
  p = UMAPS_per_lib(seurat_list[[n]])
  
  png(paste0(directory,n,"_seurat_prediction_UMAP.png"),
      width = 10, height = 15, units = "in", res = 400)
  plot(p )
  dev.off()
  
  # additional steps from mapQuery after TransferData
  query.list[[n]]= IntegrateEmbeddings(anchorset = query.anchors, 
                                       reference = integrated,
                                       query = query.list[[n]], 
                                       new.reduction.name = "ref.pca")
}

merged_query = merge(query.list[[1]], 
                     query.list[2:length(query.list)],
                     merge.data = TRUE, merge.dr = )

# don't forget to add prediction data to the subset used in the integration
AddMetaData(seurat_list[[n]], metadata = integrated$predicted.id)
rm(seurat_list)
gc()


pancreas.query = ProjectUMAP(query = pancreas.query, query.reduction = "ref.pca", reference = pancreas.integrated,
                             reference.reduction = "pca", reduction.model = "umap")


saveRDS(merged,paste0(directory,"filtered_media_merged_CCdiff_regressed_annotated.rds")) # for tests
# saveRDS(merged,"../../data/results/3.annotation_seurat/filtered_media_merged_CCdiff_regressed_50k_annotated.rds") # for tests
