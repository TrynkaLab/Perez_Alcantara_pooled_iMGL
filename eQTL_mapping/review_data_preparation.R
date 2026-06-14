# data preparation for saigeQTL in QTLight
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(hdf5r)
# library(reticulate)
# library(scCustomize)
options(future.globals.maxSize = 180000 * 1024 ^ 2) # 180Gb
# py_require("anndata")
# py_module_available("anndata")

# for new seurat assay slots
options(Seurat.object.assay.version = "v5")


directory = "../../../OTAR2065_differentiation_efficiency/data/results/1.4.markers"

for (treatment in c("untreated",
                    "IFN", 
                    "LPS"))
{
  # read in annotated seurat objects
  seu_subset = readRDS(
    paste0(
      directory,
      "/",
      treatment,
      "_harmony_addedMancusoModules/",
      treatment,
      "__harmony_addedMancusoModules.Rds"
    )
  )
  
  obj <- JoinLayers(seu_subset, assay = "RNA")
  obj@reductions <- list()
  obj <- DietSeurat(
    obj,
    assays = "RNA",
    counts = TRUE,
    data = FALSE,
    scale.data = FALSE,
    dimreducs = NULL,
    graphs = NULL,
    misc = NULL
  )
  obj[["RNA"]]$data <- NULL
  
  #class(GetAssayData(obj, layer = "counts"))
  #object.size(GetAssayData(obj, layer = "counts"))
  
  #str(obj@commands, max.level = 1)
  
  BPCells::write_matrix_hdf5(
    GetAssayData(obj, layer = "counts"),
    path = paste0("../../data/for_saigeQTL/", treatment, ".h5"),
    group = "counts"
  )
  
  # gene and cell names
  file <- H5File$new(paste0("../../data/for_saigeQTL/", treatment, ".h5"),
                     mode = "a")
  
  meta <- file$create_group("meta")
  
  meta[["row_names"]] <- rownames(obj[["RNA"]])
  meta[["col_names"]] <- colnames(obj)
  
  file$close_all()
  
  
  write.csv(obj@meta.data,
            paste0("../../data/for_saigeQTL/", treatment, "_metadata.csv"))
  rm(obj)
  rm(seu_subset)
  gc()
  
}


### another one


library(Seurat)
library(SeuratData)
library(SeuratDisk)

options(future.globals.maxSize = 180000 * 1024 ^ 2) # 180Gb
options(Seurat.object.assay.version = "v5")

directory <- "../../../OTAR2065_differentiation_efficiency/data/results/1.4.markers"

for (treatment in c("untreated", "IFN", "LPS")) {
  
  # read in annotated seurat objects
  seu_subset <- readRDS(
    paste0(
      directory, "/",
      treatment,
      "_harmony_addedMancusoModules/",
      treatment,
      "__harmony_addedMancusoModules.Rds"
    )
  )
  
  obj <- JoinLayers(seu_subset, assay = "RNA")
  
  obj@reductions <- list()
  
  obj <- DietSeurat(
    obj,
    assays = "RNA",
    counts = TRUE,
    data = FALSE,
    scale.data = FALSE,
    dimreducs = NULL,
    graphs = NULL,
    misc = NULL
  )
  obj[["RNA"]]$data <- NULL

  
  write_matrix_10x_hdf5(mat = obj[["RNA"]]$counts,
                        path = paste0("../../data/for_saigeQTL/", treatment, ".h5"))
  
  # gene / cell names
 # writeLines(rownames(counts), paste0(out_dir, "genes.txt"))
  #writeLines(colnames(counts), paste0(out_dir, "cells.txt"))
  
  # -----------------------------
  # metadata (cell-level)
  # -----------------------------
  write.csv(
    obj@meta.data,
    paste0("../../data/for_saigeQTL/", treatment, "_metadata.csv")
  )
  
  rm(obj)
  rm(seu_subset)
  gc()
}



