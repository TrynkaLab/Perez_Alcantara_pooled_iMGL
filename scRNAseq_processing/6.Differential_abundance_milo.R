# differential abundance analysis with miloR
library(miloR)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(patchwork)

options(future.globals.maxSize = 200000 * 1024^2) # 200Gb
treatment_cols =  c(untreated = "#8D918B", IFN = "#3A5683", LPS = "#F8766D")
options(Seurat.object.assay.version = "v5")

directory = "../../data/results/6.Differential_abundance_milo"
dir.create(directory, recursive = T)


seu_subset = list()
for(treatment in c("untreated", "IFN", "LPS")) {
  message("Working on treatment ", treatment)
  seu_subset[[treatment]] = readRDS(paste0("../../data/results/1.QC_v5/",treatment,
                                           "_filtered_harmony/",treatment,"_filtered_harmony.Rds"))
  format(object.size(seu_subset[[treatment]]), units = "Mb")
  metadata =  readr::read_csv(paste0( "../../data/results/1.4.markers/",treatment,"_module_annotated_metadata.csv"))
  seu_subset[[treatment]]@meta.data = metadata
  seu_subset[[treatment]][["RNA"]] = JoinLayers(seu_subset[[treatment]][["RNA"]])
  seu_subset[[treatment]] = as.SingleCellExperiment(seu_subset[[treatment]])
}