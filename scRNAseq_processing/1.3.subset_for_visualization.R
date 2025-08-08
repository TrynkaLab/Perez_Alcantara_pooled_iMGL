# 1.4. marker visualization

#.libPaths(c("/software/teamtrynka/conda/seuratv5/lib/R/library")) # path to libraries
library(Seurat)
library(tidyverse)
library(SCpubr)
library(SingleR)
library(patchwork)
library()
source("./helpers.R")
# set this option when analyzing large datasets
options(future.globals.maxSize = 60000 * 1024^2) # 60Gb
# for new seurat assay slots
options(Seurat.object.assay.version = "v5")

directory = "../../data/results/1.4.markers"
dir.create(directory, recursive = T)

# markers from microglial subtypes
subtypes_Sun_full = readr::read_tsv("../../../resources/Sun_2023/Sun_2023_microglial_states_in_AD_markers_renamed.txt") %>%
  group_by(microglial_state) %>%
  dplyr::slice_head(n=30) %>%
  summarise(gene = list(gene)) %>%
  ungroup()
subtypes_Sun_full = as.list(set_names(subtypes_Sun_full$gene, subtypes_Sun_full$microglial_state))

# get subsets of treatments
seu_subset = list()
for(treatment in c("untreated", "IFN", "LPS")) {
  message("Working on treatment ", treatment)
  seu_subset[[treatment]] = readRDS(paste0("../../data/results/1.3.subset_for_visualization/",treatment,"_filtered_harmony_subset/",
                                           treatment,"_filtered_harmony_50k_subset_seuratv5.Rds"))
  
  DefaultLayer(seu_subset[[treatment]][["RNA"]]) = 'data'
}
seuratobj = merge(seu_subset$untreated, y = c(seu_subset$IFN, seu_subset$LPS),
                  add.cell.ids = c("untreated", "IFN", "LPS"), project = "visualization",
                  merge.data = TRUE)



# # score cells
seuratobj = AddModuleScore(seuratobj,
                           features = subtypes_Sun_full,
                           assay = "RNA",slot = "data",
                           name = names(subtypes_Sun_full))

# saving
saveRDS(seuratobj, file = paste0(directory,"/merged_filtered_harmony_subset_addedModule.Rds"))


seuratobj = readRDS(paste0(
  directory,
  "/merged_filtered_harmony_subset_addedModule.Rds"
))


# subset of genes used for initial annotation in Sun's paper
selected_gene_expr = list("Pluripotency" = c("POU5F1","SOX2"),
                          "General_microglia" = c("CSF1R","CD74","C3"),
                          "Homeostatic"=c("P2RY12","CX3CR1"),
                          "Brain_macrophages_phagocytic"=c("MRC1","CD163","LYVE1","F13A1"),
                          "Neuronal_surveillance" = c("FRMD4A",
                                                      # "GRCK3", # not found 
                                                      # "SLC6A3", # 0
                                                      "ASMTL"),
                          "Inflammatory_I" = c("INO80D","TMEM163","CPEB4"),
                          "Ribosome_biogenesis" = c("FTL","FTH1"),
                          "Lipid_processing" = c("MYO1E","PTPRG"),
                          "Stress" = c("HSP90AA1","HSPH1"),
                          "Glycolytic" = c("NAMPT","SLC2A3"),
                          "Inflammatory_II" = c("SPON1","LRRK2"),
                          "Inflammatory_III" = c("CCL3","IL1B"),
                          "Antiviral" = c("IFI44L1","MX1"),
                          "Cycling" = c("EZH2","BRIP1")
                          
)

pdf(paste0(directory,"/subtypes_Sun_top_",treatment,"_gene_expr.pdf"),width = 10,height = 10)

for(nam in  paste0(names(selected_gene_expr))){
  message("Working on microglia type ",nam)
  # do_DimPlot doesn't work on continuous measures
  p =  FeaturePlot(subset(x = seu_subset[[treatment]], downsample = 10000), # need to downsample otherwise takes forever
                   reduction = "umap.full",
                   features  = selected_gene_expr[[nam]],
                   raster = TRUE, 
                   order=TRUE,
                   pt.size = 3,
                   # legend.title = "Normalised gene expression",
                   # plot.title = 
                   #  assay = "RNA",
                   slot = "data") 
  
  plot(p)
}

dev.off()

pdf(paste0(directory,"/subtypes_Sun_top30_",treatment,"_module_scores.pdf"),width = 10,height = 10)

for(columns in  paste0(names(subtypes_Sun_full),c(1:12))){
  message("Working on column ",columns)
  # do_DimPlot doesn't work on continuous measures
  p =  SCpubr::do_FeaturePlot(seu_subset[[treatment]],
                              reduction = "umap.full",
                              features  = columns,
                              raster = TRUE, 
                              pt.size = 1,
                              order = TRUE) + 
    ggtitle(paste0(columns, " module score"))
  
  plot(p)
}


dev.off()


#

