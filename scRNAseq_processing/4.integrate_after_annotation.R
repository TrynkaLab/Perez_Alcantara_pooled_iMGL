# Integrating data after annotation
.libPaths(c("/software/teamtrynka/ma23/R4.1/libs","/software/R-4.1.0/lib/R/library"))
library(patchwork)
library(tidyverse)
library(Seurat)
library(future)
library(SCpubr)
source("./helpers.R")


plan("multisession", workers = 2) 
plan()
options(future.globals.maxSize= 140097152000) #140 Gb
# options(future.globals.maxSize= 10097152000) #10 Gb
directory = "../../data/results/4.integrate_after_annotation/"
dir.create(directory, recursive = T)

merged = readRDS("../../data/results/3.annotation_seurat/filtered_media_merged_CCdiff_regressed_50k_annotated.rds") # for tests

# Integration ###########
# Integrate data for joint analysis - might obfuscate some of the subtler biological differences between the experiments
# with atomic sketch https://satijalab.org/seurat/articles/atomic_integration.html

seurat_list = SplitObject(input_merged, split.by = "orig.ident")

rm(input_merged)
gc()
# remove objects that have fewer than 2200 cells
#$pool3.h3_diff2_LPS - too few cells. What to do with this? Remove for now
# [1] 1601
lapply(seurat_list,ncol)
seurat_list = Filter(function(x) ncol(x) > 2200, seurat_list)


# Rename cells to avoid future conflicts
for(n in names(seurat_list)){
  seurat_list[[n]] = RenameCells(object = seurat_list[[n]], 
                                 new.names = paste(seurat_list[[n]]@meta.data$orig.ident,
                                                   seurat_list[[n]]@meta.data$orig_cell,
                                                   sep = "_"))
  
}


atoms.list = list()
for (n in names(seurat_list)) {
  subset_seurat = seurat_list[[n]]
  DefaultAssay(subset_seurat) = "RNA"
  
  # basic preprocessing
  subset_seurat = NormalizeData(subset_seurat)
  subset_seurat = FindVariableFeatures(subset_seurat)
  
  # calculate leverage score and sample 5000 cells based on leverage score
  # atoms.i = LeverageScoreSampling(object = subset_seurat, num.cells = 5000)
  atoms.i = LeverageScoreSampling(object = subset_seurat, num.cells = floor(0.05 * nrow(subset_seurat@meta.data)))
  atoms.list[[n]] = atoms.i
  
  # delete full object from memory note that this is optional, if you can store the full
  # datasets in memory, you dont have to reload them later
  rm(subset_seurat)
}

## Next we perform integrative analysis on the ‘atoms’ from each of the dataset
for (n in names(seurat_list)) {
  atoms.list[[n]] = SCTransform(atoms.list[[n]], verbose = TRUE)
}

# perform integration
features = SelectIntegrationFeatures(object.list = atoms.list)
atoms.merge = FastRPCAIntegration(object.list = atoms.list, dims = 1:30, 
                                  normalization.method = "SCT",
                                  anchor.features = features)

# we can generate a 2D visualization representing the integrated atoms
atom.reduction = "integrated_dr"
atoms.merge = RunUMAP(atoms.merge, 
                      reduction = atom.reduction, 
                      dims = 1:30, return.model = TRUE)
SCpubr::do_DimPlot(atoms.merge, group.by = "orig.ident")
SCpubr::do_DimPlot(atoms.merge, group.by = "pool")

## we can now place other cells from each dataset in this space as well.
# note: all datasets need to have a subset placed in the atoms. Otherwise it won't work

integrated_objects = list()
for (n in names(seurat_list)) {
  
  subset_seurat = seurat_list[[n]]
  DefaultAssay(subset_seurat) = "RNA"
  
  subset_seurat = NormalizeData(subset_seurat)
  subset_seurat = FindVariableFeatures(subset_seurat, verbose = TRUE)
  subset_seurat = ScaleData(subset_seurat, verbose = TRUE)
  subset_seurat = CellCycleScoring(subset_seurat,  
                                   s.features=cc.genes.updated.2019$s.genes, 
                                   g2m.features=cc.genes.updated.2019$g2m.genes)
  # Integrate all cells into the same space as the atoms
  subset_seurat = IntegrateSketchEmbeddings(object = subset_seurat, 
                                            atom.sketch.object = atoms.merge, 
                                            atom.sketch.reduction = atom.reduction,
                                            features = features)
  
  
  integrated_objects[[n]] = subset_seurat
  rm(subset_seurat)
}

media_integrated = merge(integrated_objects[[1]], 
                         integrated_objects[2:length(integrated_objects)],
                         merge.dr = c("integrated_dr"))
length(media_integrated@assays$RNA@var.features) # variable features are not present

media_integrated = FindNeighbors(media_integrated, 
                                 dims = 1:30, 
                                 reduction = "integrated_dr")
media_integrated = FindClusters(media_integrated, 
                                resolution = 0.3, # double-check this
                                verbose = TRUE)
media_integrated = RunUMAP(media_integrated, 
                           reduction = "integrated_dr",
                           dims = 1:30)

p = UMAPS_all_libs(media_integrated)

p1 = SCpubr::do_DimPlot(media_integrated,
                        label = F,  group.by = "predicted.id")  +
  ggtitle("Seurat predictions")

png(paste0(directory,"media_integrated_UMAP.png"),
    width = 10, height = 18, units = "in", res = 400)
p / (p1 + plot_spacer())
dev.off()

# UMAP split for predicted id

p1 = SCpubr::do_DimPlot(media_integrated,
                        label = F,  group.by = "predicted.id", split.by = "predicted.id")  +
  ggtitle("Seurat predictions")

png(paste0(directory,"media_integrated_UMAP.png"),
    width = 18, height = 10, units = "in", res = 400)
p1
dev.off()


png(paste0(directory,"media_integrated_UMAP_clusters.png"),
    width = 5, height = 6, units = "in", res = 400)
SCpubr::do_DimPlot(media_integrated, 
                   reduction = "umap", 
                   group.by = "seurat_clusters")
dev.off()


SCpubr::do_FeaturePlot(media_integrated, 
                       features  = "APOE", assay = "RNA",slot = "data", order = TRUE) 
SCpubr::do_FeaturePlot(media_integrated, 
                       features  = "APOE", assay = "RNA",slot = "counts", order = TRUE) 
SCpubr::do_FeaturePlot(media_integrated, 
                       features  = "APOE", assay = "RNA",slot = "counts",
                       split.by = "seurat_clusters") 
SCpubr::do_FeaturePlot(media_integrated, 
                       features  = "APOE", assay = "RNA",slot = "data",
                       split.by = "seurat_clusters") 
SCpubr::do_FeaturePlot(media_integrated, 
                       features  = "APOE", assay = "RNA",slot = "scale.data",
                       split.by = "seurat_clusters") # this does

saveRDS(media_integrated, paste0(directory,"media_integrated_50k_annotated.rds"))
