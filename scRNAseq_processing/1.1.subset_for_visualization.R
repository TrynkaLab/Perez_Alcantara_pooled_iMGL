# subset to 100k cells for visualization and then save in old seurat format

library(Seurat)
library(BPCells)
library(SeuratObject)
source("./helpers_seuratv5.R")
# set this option when analyzing large datasets
options(future.globals.maxSize = 3e+09)
# for new seurat assay slots
options(Seurat.object.assay.version = "v5")

directory = "../../data/results/1.1.subset_for_visualization/"
dir.create(directory, recursive = T)

# get downsampled treatments
 for(treatment in c("untreated", "IFN", "LPS")) {
   rm(input_merged, input_subset)
   gc()
  options(Seurat.object.assay.version = "v5")
  input_merged = readRDS(paste0("../../data/results/1.QC_v5/",treatment,"_filtered_harmony/",treatment,"_filtered_harmony.Rds"))
  DefaultAssay(input_merged) = "RNA"
  Idents(input_merged) = "treatment"
  input_subset = subset(x = input_merged, downsample = 50000) # downsample to 50k cells
  input_subset = JoinLayers(input_subset)
  # seuratv3 = as(object = input_seurat[["RNA"]], Class = "Assay")
  # options(seurat.object.assay.version = "v3")
  # 
  # seuratv3 = CreateSeuratObject(counts = seuratv3@counts,
  #                               meta.data =input_seurat@meta.data )
  
  # DimPlot(input_subset,group.by = "pool") # check the dimensional reductions are there
  # DimPlot(input_subset,group.by = "pool",reduction = "umap.full")
  write_matrix_dir(mat = input_subset[["RNA"]]$counts, 
                   dir = paste0(directory,"/",treatment,"_filtered_harmony_subset"),
                   overwrite = TRUE)
  input_subset[["RNA"]]$counts = open_matrix_dir(dir = paste0(directory,"/",treatment,"_filtered_harmony_subset"))
  
  saveRDS(
    object = input_subset,
    file = paste0(directory,"/",treatment,"_filtered_harmony_subset/",treatment,"_filtered_harmony_50k_subset_seuratv5.Rds"))
  

}
