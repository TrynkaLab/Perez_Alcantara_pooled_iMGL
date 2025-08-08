# Analysing differentiation efficiency
# pools 2, 3, 4, 5, 6, 7, 8, 9, 10
.libPaths(c("/software/R-4.1.0/lib/R/library","/software/teamtrynka/ma23/R4.1/libs"))
library(patchwork)
library(tidyverse)
.libPaths(c("/software/teamtrynka/ma23/R4.1/libs","/software/R-4.1.0/lib/R/library"))
library(Seurat)
library(SCpubr)
library(RColorBrewer)
library(janitor)
source("./helpers.R")

plan("multisession", workers = 2) 
plan()
options(future.globals.maxSize= 100097152000) #100 Gb
# options(future.globals.maxSize= 10097152000) #10 Gb
directory = "../../data/results/1.QC/"
dir.create(directory, recursive = T)

### start of analysis ---------
# list of microglia input data from different pools and conditions: 

# Loads a list of gene expression (feature) counts per cell
# cellranger_data = list()
# for(pool in c("pool2","pool3","pool4","pool5","pool6", "pool7", "pool8","pool9","pool10")){
#   message(pool)
#   files = list.dirs( paste0("../../../OTAR2065_scRNA_data_processing/data/cellranger/",
#                             pool),recursive = FALSE, full.names = FALSE)
#   for(file in files){
#     cellranger_data[[paste0(pool,".",file)]] =  Read10X(data.dir = paste0("../../../OTAR2065_scRNA_data_processing/data/cellranger/",
#                                                                           pool,"/",file,
#                                                                           "/outs/filtered_feature_bc_matrix"),
#                                                         strip.suffix = T)  ## remove trailing -1 from cell barcodes
#   }
# }
# str(cellranger_data)
# 
# # Initialize the Seurat object with the raw (non-normalized data). No filtering yet
# 
# seurat_list = cellranger_to_seurat(cellranger_data)
# rm(cellranger_data)
# gc()
# 
# for(n in names(seurat_list)){
# 
#   seurat_list[[n]][["percent.mt"]] = PercentageFeatureSet(seurat_list[[n]], pattern = "^MT-")
# 
# }
# # 
# # 
# # # Add Vireo line ID info to each cell
# seurat_list  = incorporate_vireo_info(seurat_list)
# gc()
# ##### Check count metrics before filtering ####
# 
# ncells = sapply(seurat_list,ncol)
# sum(ncells) # ~1.4 million cells pools 2-10
# 
# # per pool and per treatment so the object is not impossibly large
# 
# for(treatment in c("plst|untreated|AG|YC|ITMG", "IFN", "LPS")) {
#   message("Working on ",treatment," ...")
#   message("Subsetting...")
#   subset_list = seurat_list[grep(treatment, names(seurat_list))]
#   message("Merging...")
#   input_merged =  merge(subset_list[[1]],
#                         subset_list[2:length(subset_list)],
#                         add.cell.ids = names(subset_list))
#   input_merged@active.ident = factor(x = input_merged@active.ident, 
#                                      levels = names(subset_list), 
#                                      ordered = T)
#   rm(subset_list)
#   gc()
#   
#   # Visualize QC metrics 
#   message("Plotting QC metrics...")
#   
#   p1 = plot_QC_histograms(input_merged, plot.legend = FALSE)
#   
#   # check correlation of assignment probability and UMI counts, and also SNPs covered
#   p2 = plot_QC_vireo_dotplots(input_merged)
#   
#   message("Saving...")
#   if(treatment == "plst|untreated|AG|YC|ITMG"){
#     saveRDS(input_merged, paste0(directory,"untreated_all_pools_unfiltered_media_merged_noRegression.rds"))
#     input_merged@meta.data %>% 
#       group_by(pool) %>%
#       summarise(ncell = n()) %>%
#       write.csv(., paste0(directory,"/untreated_all_pools_unfiltered_media_merged_cellN.csv"),
#                 quote = FALSE, row.names = FALSE)
#     png(paste0(directory,"/untreated_all_pools_media_merged_metrics_cell_density_perUMI_perGene_perMitPercent.png"),
#         width = 12, height = 3, units = "in", res = 400)
#     plot(p1)
#     dev.off()
#     
#     png(paste0(directory,"/untreated_all_pools_media_merged_metrics_QC_vireo_dotplots.png"),
#         width = 12, height = 5, units = "in", res = 400)
#     # This info doesn't include doublets 
#     plot(p2)
#     dev.off()
#   }else{
#     saveRDS(input_merged, paste0(directory,treatment,"_all_pools_unfiltered_media_merged_noRegression.rds"))
#     input_merged@meta.data %>% 
#       group_by(pool) %>%
#       summarise(ncell = n()) %>%
#       write.csv(., paste0(directory,treatment,"_all_pools_unfiltered_media_merged_cellN.csv"),
#                 quote = FALSE, row.names = FALSE)
#     
#     png(paste0(directory,treatment,"_all_pools_media_merged_metrics_cell_density_perUMI_perGene_perMitPercent.png"),
#         width = 12, height = 3, units = "in", res = 400)
#     plot(p1)
#     dev.off()
#     
#     png(paste0(directory,treatment,"_all_pools_media_merged_metrics_QC_vireo_dotplots.png"),
#         width = 12, height = 5, units = "in", res = 400)
#     # This info doesn't include doublets 
#     plot(p2)
#     dev.off()
#     
#   }
#   
# }
# # per pool  
# for(pl in c("pool2","pool3","pool4","pool5","pool6", "pool7", "pool8","pool9","pool10")){
#   # There's one UMI count and gene count outlier (pool 3 LPS). Check by plotting per pool:
#   
#   message("Working on ",pl," ...")
#   message("Subsetting...")
#   subset_list = seurat_list[grep(pl, names(seurat_list))]
#   message("Merging...")
#   input_merged =  merge(subset_list[[1]],
#                         subset_list[2:length(subset_list)],
#                         add.cell.ids = names(subset_list))
#   input_merged@active.ident = factor(x = input_merged@active.ident, 
#                                      levels = names(subset_list), 
#                                      ordered = T)
#   rm(subset_list)
#   gc()
#   
#   # Visualize QC metrics 
#   message("Plotting QC metrics...")
#   
#   p = plot_QC_histograms(seurat_merged_object = subset(x = input_merged, subset = pool == pl),
#                          plot.legend = TRUE)
#   
#   png(paste0(directory,"all_treatment_media_merged_metrics_cell_density_perUMI_perGene_perMitPercent_" ,
#              pl,
#              ".png"),width = 12, height = 3, units = "in", res = 400)
#   # Cells have been filtered for around 500 UMIs per cell min in each cell when creating the Seurat object, or does cellranger do that?
#   plot(p)
#   dev.off()
#   
#   input_merged@meta.data %>%
#     group_by(treatment) %>%
#     summarise(ncell = n()) %>%
#     write.csv(., paste0(directory,pl,"_all_treatment_unfiltered_media_merged_cellN.csv"),
#               quote = FALSE, row.names = FALSE)
#   
# }
# 
# 
# ######### End of count check
# rm ( input_merged,p)
# gc()
# 
# ######### Filtering #########
# 
# ### RE-DO FILTER SO THAT IT IS DONE FOR BOTTOM 5% OF FEATURE COUNTS
# ### We filter cells that have > 10% mitochondrial counts (same filter as before, the usual one for human cells [it depends on cell type and tissue], although many here seem above)
# ## and we remove doublets and unassigned cells
# # splits object for filtering
# 
# seurat_list = filter_x_UMIcount_bottompercent_per_library(seurat_list, bottom_percent = 5)
# ncells = sapply(seurat_list,ncol)
# sum(ncells) # ~1.4 million cells pools 2-10
# # filter doublets, unassigned, and those with higher than 10% mitochondrial percentage
# for (sample in names(seurat_list)){
#   message("working on ",sample)
#   seurat_list[[sample]] = subset(seurat_list[[sample]], subset = percent.mt < 10)
#   seurat_list[[sample]] = subset(seurat_list[[sample]], subset = donor_id!="unassigned")
#   seurat_list[[sample]] = subset(seurat_list[[sample]], subset = donor_id!="doublet")
# 
# }
# 
# ncells = sapply(seurat_list,ncol)
# sum(ncells) # 926554 million cells remaining after filter pools 2-10
# 
# 
# saveRDS(seurat_list, paste0(directory,"filtered_media_noRegression_list.rds"))

##############
## normalise per treatment otherwise throws memory error #########
###############
gc()
# for(treat in c("untreated","LPS","IFN") ){
for(treat in c("LPS","IFN") ){
  
  seurat_list = readRDS(paste0(directory,"filtered_media_noRegression_list.rds"))
  gc()
  ncells = sapply(seurat_list,ncol)
  message("Number of cells in full filtered object:", 
          sum(ncells)) # 926554 million cells remaining after filter pools 2-10
  
  if(treat =="untreated"){
    seurat_list = seurat_list[grep("plst|untreated|AG|YC|ITMG", names(seurat_list))] # subset per treatment and merge
    
  }else{
    seurat_list = seurat_list[grep(treat, names(seurat_list))] # subset per treatment and merge
    
  }
  
  message("Merging...")
  input_merged =  merge(seurat_list[[1]],
                        seurat_list[2:length(seurat_list)],
                        add.cell.ids = names(seurat_list))
  input_merged@active.ident = factor(x = input_merged@active.ident, 
                                     levels = names(seurat_list), 
                                     ordered = T)
  rm(seurat_list)
  gc()
  message("Number of cells in full filtered object for treatment " ,treat,": ", 
          ncol(input_merged)) 
  # untreated: 217177 cells
  #### merge per treatment and normalise
  message("Normalising...")
  
  input_merged= NormalizeData(input_merged, verbose = TRUE)
  message("Finding variable features...")
  
  input_merged = FindVariableFeatures(input_merged, verbose = TRUE)
  message("Scaling...")
  
  input_merged = ScaleData(input_merged, verbose = TRUE)
  gc()
  ## Scoring genes by cell cycle gene expression - from Tirosh et al., 2016 ----------
  message("Scoring cell cycle...")
  
  input_merged = CellCycleScoring(input_merged,  
                                  s.features=cc.genes.updated.2019$s.genes, 
                                  g2m.features=cc.genes.updated.2019$g2m.genes)
  message("Running PCA...")
  
  input_merged = RunPCA(input_merged, verbose = TRUE)
  message("Finding kNN...")
  
  input_merged = FindNeighbors(input_merged, dims = 1:10)
  message("Clustering...")
  
  input_merged = FindClusters(input_merged, resolution = 0.5)
  message("Running UMAP...")
  
  input_merged = RunUMAP(input_merged, dims = 1:50)
  
  ### Perform dimensionality reduction by PCA and UMAP embedding
  
  input_merged$orig.ident = factor(input_merged$orig.ident, 
                                   levels = unique(input_merged$orig.ident), 
                                   ordered = T)
  message("Saving object...")
  
  saveRDS(input_merged, paste0("../../data/results/1.QC/filtered_media_merged_normalised_noRegression_",treat,".rds"))
  
  
  rm(input_merged)
  gc()
}


input_merged = readRDS("../../data/results/1.QC/filtered_media_merged_normalised_noRegression_untreated.rds")

message("Downsampling...")

input_merged_downsampled = input_merged[, sample(colnames(input_merged),
                                                 size =200000, replace=F)]
message("Saving smaller object...")

saveRDS(input_merged_downsampled, paste0(directory,"filtered_media_merged_normalised_noRegression_200k_untreated.rds"))

p = UMAPS_all_libs(input_merged_downsampled)

png(paste0("../../data/results/1.QC/filtered_media_merged_UMAP_nothingRegressed_200k_untreated.png"),
    width = 11, height = 14, units = "in", res = 400)
plot(p)
dev.off()

input_merged = readRDS("../../data/results/1.QC/filtered_media_merged_normalised_noRegression_LPS.rds")


message("Downsampling...")

input_merged_downsampled = input_merged[, sample(colnames(input_merged),
                                                 size =200000, replace=F)]
message("Saving smaller object...")

saveRDS(input_merged_downsampled, paste0(directory,"filtered_media_merged_normalised_noRegression_200k_LPS.rds"))

p = UMAPS_all_libs(input_merged_downsampled)

png(paste0("../../data/results/1.QC/filtered_media_merged_UMAP_nothingRegressed_200k_LPS.png"),
    width = 11, height = 14, units = "in", res = 400)
plot(p)
dev.off()

input_merged = readRDS("../../data/results/1.QC/filtered_media_merged_normalised_noRegression_IFN.rds")
message("Downsampling...")

input_merged_downsampled = input_merged[, sample(colnames(input_merged),
                                                 size =200000, replace=F)]
message("Saving smaller object...")

saveRDS(input_merged_downsampled, paste0(directory,"filtered_media_merged_normalised_noRegression_200k_IFN.rds"))

p = UMAPS_all_libs(input_merged_downsampled)

png(paste0("../../data/results/1.QC/filtered_media_merged_UMAP_nothingRegressed_200k_IFN.png"),
    width = 11, height = 14, units = "in", res = 400)
plot(p)
dev.off()
# thinking about using harmony for pool regression
# check donor abundance per cluster (barplots)
