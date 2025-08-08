# integration on untreated pools
# pools 2 - 17
# using seurat v5 

library(Seurat)
library(BPCells)
library(SeuratObject)
library(Signac)
# library(patchwork)
library(tidyverse)
library(harmony)
# library(future)
#library(SCpubr)
source("./helpers_seuratv5.R")
# set this option when analyzing large datasets
options(future.globals.maxSize = (50000 * 1024 ^ 2)) #(~50 Gb)
# for new seurat assay slots
options(Seurat.object.assay.version = "v5")

directory = "../../data/results/1.QC_v5"
dir.create(directory, recursive = T)

### start of analysis ---------

cellranger_data = list()
for(pool in c(paste0("pool",2:17))){
  message(pool)
  files = list.dirs( paste0("../../../OTAR2065_scRNA_data_processing/data/cellranger/",
                            pool),recursive = FALSE, full.names = FALSE)
  for(file in files){
    cellranger_data[[paste0(pool,".",file)]] =  Read10X(data.dir = paste0("../../../OTAR2065_scRNA_data_processing/data/cellranger/",
                                                                          pool,"/",file,
                                                                          "/outs/filtered_feature_bc_matrix"),
                                                        strip.suffix = T)  ## remove trailing -1 from cell barcodes
  }
}
str(cellranger_data)
seuratobj = cellranger_to_seurat(cellranger_data )
class(seuratobj[["RNA"]])
# NULL
class(seuratobj$pool3.h3_diff2_IFNg[["RNA"]])
# [1] "Assay5"
# attr(,"package")
# [1] "SeuratObject"
for(n in names(seuratobj)){
  
  SaveSeuratRds(
    object = seuratobj[[n]],
    file = paste0(directory,"/",n,".Rds"))
  
  message("BPCells count matrix size: ", format(object.size(seuratobj[[n]]), units = "Mb"))
  
}

### In case I need to load
names = list.files(directory,pattern="pool.*.Rds$",)
names = str_replace(names,".Rds","")
names = names[!str_detect(names,"filtered")] # remove this matches that will appear later in this script
seuratobj = list()
for(n in names){
  seuratobj[[n]] = LoadSeuratRds(
    file = paste0(directory,"/",n,".Rds"))
}

sapply(seuratobj,function(x) summary(x@meta.data$nCount_RNA) ) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("values") %>%
  tidyr::pivot_longer(cols = -values,names_to = "object",values_to = "nCount_RNA") %>%
  group_by(values) %>%
  dplyr::summarise(min = min(nCount_RNA),
                   max = max(nCount_RNA),
                   mean = mean(nCount_RNA),
                   median = median(nCount_RNA))
#####

for(n in names(seuratobj)){
  message(paste0("Working on ",n))
  seuratobj[[n]][["percent.mt"]] = PercentageFeatureSet(seuratobj[[n]], pattern = "^MT-")
  seuratobj[[n]][["percent.ribo"]] = PercentageFeatureSet(seuratobj[[n]], pattern = "^RPS|^RPL")
  
}
head(seuratobj$pool2.AG_10)
#
# # # Add Vireo line ID info to each cell
seuratobj  = incorporate_vireo_info(seuratobj,pools = paste0("pool",2:17))

#### general stats  pre-filter ######
ncells = sapply(seuratobj,ncol)
sum(ncells) # 3,404,020 cells pools 2-17
sapply(seuratobj,function(x) summary(x@meta.data$nCount_RNA) ) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("values") %>%
  tidyr::pivot_longer(cols = -values,names_to = "object",values_to = "nCount_RNA") %>%
  group_by(values) %>%
  dplyr::summarise(min = min(nCount_RNA),
                   max = max(nCount_RNA),
                   mean = mean(nCount_RNA),
                   median = median(nCount_RNA))
# values     min     max   mean median
# <chr>    <dbl>   <dbl>  <dbl>  <dbl>
#   1 1st Qu.  1490   33053.  6473.  6185 
# 2 3rd Qu.  5376   62621  13857. 12924.
# 3 Max.    26220  147810  57924. 55281 
# 4 Mean     3827.  48273. 10700. 10199.
# 5 Median   2934   47064   9640.  9005 
# 6 Min.      500    2490   1155.  1039 

sapply(seuratobj,function(x) summary(x@meta.data$percent.mt) ) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("values") %>%
  tidyr::pivot_longer(cols = -values,names_to = "object",values_to = "percent.mt") %>%
  group_by(values) %>%
  dplyr::summarise(min = min(percent.mt),
                   max = max(percent.mt),
                   mean = mean(percent.mt),
                   median = median(percent.mt))
# values     min     max     mean median
# <chr>    <dbl>   <dbl>    <dbl>  <dbl>
#   1 1st Qu.  0.643  6.01    3.24      3.44
# 2 3rd Qu.  2.30  11.5     5.88      6.51
# 3 Max.    24.9   93.6    64.5      69.4 
# 4 Mean     2.36   9.10    5.59      5.71
# 5 Median   1.56   7.39    4.23      4.58
# 6 Min.     0      0.0847  0.00805   0   


#######
head(seuratobj$pool2.AG_10)


# now filter per original cellranger sample,
seuratobj = filter_x_UMIcount_bottompercent_per_library(seuratobj, bottom_percent = 5)

ncells = sapply(seuratobj,ncol)
sum(ncells) # 3,233,612 cells pools 2-17
head(seuratobj$pool2.AG_10)

for(sample in names(seuratobj)){
  message(paste0("Working on ",sample))
  seuratobj[[sample]] = subset(seuratobj[[sample]], subset = percent.mt < 10)
  seuratobj[[sample]] = subset(seuratobj[[sample]], subset = donor_id!="unassigned")
  seuratobj[[sample]] = subset(seuratobj[[sample]], subset = donor_id!="doublet")
  
}

ncells = sapply(seuratobj,ncol)
sum(ncells) #2,118,070 pools 2-17
names(seuratobj)

######## REMOVING MALAT1 and ribosomal genes
for(sample in names(seuratobj)){
  message(paste0("Working on ",sample))
  counts = GetAssayData(seuratobj[[sample]] , assay = "RNA", layer="counts")
  ribo_genes = rownames(seuratobj[[sample]] )[grep("^RPS|^RPL", rownames(seuratobj[[sample]] ))]
  counts = counts[-(which(rownames(counts) %in% c("MALAT1", ribo_genes))),]
  seuratobj[[sample]]  = subset(seuratobj[[sample]] , features = rownames(counts))
}

# for(n in names(seuratobj)){
#   
#   SaveSeuratRds(
#     object = seuratobj[[n]],
#     file = paste0(directory,"/",n,"_filtered.Rds"))
#   
#   message("BPCells count matrix size: ", format(object.size(seuratobj[[n]]), units = "Mb"))
#   
# }
# 
# ### In case I need to load
# 
# names = list.files(directory,pattern="*_filtered.Rds")
# names = str_replace(names,"_filtered.Rds","")
# seuratobj = list()
# for(n in names){
#   
#   seuratobj[[n]] = LoadSeuratRds(
# 
#     file = paste0(directory,"/",n,"_filtered.Rds"))
#   
# }

########## problematic bits depending on CPU architecture ########

for(n in names(seuratobj)){
  
  # with BPcells
  write_matrix_dir(mat = seuratobj[[n]][["RNA"]]$counts,
                   dir = paste0(directory,"/",n,"_filtered"),
                   overwrite = TRUE,
                   compress=TRUE)
  seuratobj[[n]][["RNA"]]$counts = open_matrix_dir(dir = paste0(directory,"/",n,"_filtered"))
  
  # SaveSeuratRds(
  #   object = seuratobj[[n]],
  #   file = paste0(directory,"/",n,"_filtered/",n,"_filtered.Rds"))
  
  saveRDS(
    object = seuratobj[[n]],
    file = paste0(directory,"/",n,"_filtered/",n,"_filtered.Rds"))
  
  message("BPCells count matrix size: ", format(object.size(seuratobj[[n]]), units = "Mb"))
  
}



#######
### In case I need to load
# names = list.files(directory,pattern="*_filtered.Rds")
# names = str_replace(names,"_filtered.Rds","")
# seuratobj = list()
# for(n in names){
#   seuratobj[[n]] = LoadSeuratRds(
#    file =paste0(directory,"/",n,"_filtered/",n,"_filtered.Rds"))
# 
# }

for(treatment in c("IFN", "LPS","plst|untreated|Untreated|AG|YC|ITMG")) {
  # am I missing INFg? -No
  
  
  ##########
  message("Working on ",treatment," ...")
  message("Subsetting...")
  print(names(seuratobj))
  subset_list = seuratobj[grep(treatment, names(seuratobj))]
  message("Merging...")
  input_merged =  merge(subset_list[[1]],
                        subset_list[2:length(subset_list)],
                        add.cell.ids = names(subset_list))
  input_merged@active.ident = factor(x = input_merged@active.ident,
                                     levels = names(subset_list),
                                     ordered = T)
  rm(subset_list)
  gc()
  
  # join layers
  input_merged[["RNA"]] = JoinLayers(input_merged[["RNA"]])
  Layers(input_merged[["RNA"]])
  
  
  
  
  message("BPCells count matrix size: ", format(object.size(input_merged), units = "Mb"))
  
  
  message("Saving...")
  if(treatment == "plst|untreated|Untreated|AG|YC|ITMG"){
    write_matrix_dir(mat = input_merged[["RNA"]]$counts,
                     dir = paste0(directory,"/untreated_filtered"),
                     overwrite = TRUE)
    input_merged[["RNA"]]$counts = open_matrix_dir(dir = paste0(directory,"/untreated_filtered"))
    
    saveRDS(
      object = input_merged,
      file = paste0(directory,"/untreated_filtered/untreated_filtered.Rds"))
    message("BPCells count matrix size for one-layer treatment: untreated: ", format(object.size(input_merged), units = "Mb"))
    print(unique(input_merged$pool))
    
  }else{
    write_matrix_dir(mat = input_merged[["RNA"]]$counts,
                     dir = paste0(directory,"/",treatment,"_filtered"),
                     overwrite = TRUE)
    input_merged[["RNA"]]$counts = open_matrix_dir(dir = paste0(directory,"/",treatment,"_filtered"))
    
    saveRDS(
      object = input_merged,
      file = paste0(directory,"/",treatment,"_filtered/",treatment,"_filtered.Rds"))
    message("BPCells count matrix size for one-layer treatment: ",treatment,": ", format(object.size(input_merged), units = "Mb"))
    print(unique(input_merged$pool))
    
  }
  
}


# sketch shared donors and integrate the rest
for(treatment in c("IFN", "LPS","untreated")) {

  rm(input_merged)
  gc()
  input_merged = readRDS(paste0(directory,"/",treatment,"_filtered/",treatment,"_filtered.Rds"))
  message("BPCells count matrix size for one-layer treatment: ",
          treatment,": ", format(object.size(input_merged), units = "Mb"))
  print(unique(input_merged$pool))
  input_merged 
  
  # without any integration or correction
  # split layers
  input_merged[["RNA"]] = split(input_merged[["RNA"]], f = input_merged$pool)
  Layers(input_merged[["RNA"]])
  input_merged = NormalizeData(input_merged)
  input_merged = FindVariableFeatures(input_merged)
  
  input_merged = SketchData(object = input_merged, ncells = 5000, 
                            method = "LeverageScore", sketched.assay = "sketch")
  input_merged
  DefaultAssay(input_merged) = "sketch"
  input_merged = FindVariableFeatures(input_merged)
  input_merged = ScaleData(input_merged)
  input_merged = RunPCA(input_merged,verbose = TRUE) 
  input_merged = FindNeighbors(input_merged, dims = 1:50,verbose = TRUE)
  input_merged = FindClusters(input_merged, resolution = 2,verbose = TRUE)
  input_merged = RunUMAP(input_merged, dims = 1:50, return.model = T,verbose = TRUE)
  
  
  message("BPCells count matrix size after normalisation per poolt: ",
          treatment,": ", format(object.size(input_merged), units = "Mb"))
  # 4750.1 Mb
  
  tosave = input_merged
  tosave[["RNA"]] = JoinLayers(tosave[["RNA"]])
  
  
  write_matrix_dir(mat = tosave[["RNA"]]$counts,
                   dir = paste0(directory,"/",treatment,"_filtered_normalised"),
                   overwrite = TRUE)
  tosave[["RNA"]]$counts = open_matrix_dir(dir = paste0(directory,"/",treatment,"_filtered_normalised"))
  
  saveRDS(
    object = tosave,
    file = paste0(directory,"/",treatment,"_filtered_normalised/",treatment,"_filtered_normalised.Rds"))
  
  # 
  # p =  DimPlot(input_merged, 
  #              reduction = "umap", 
  #              group.by = "pool",
  #              #split.by = "pool", # throws error
  #              alpha = 0.8) 
  # 
  # # p = SCpubr::do_DimPlot(input_merged, 
  # #                        group.by = "pool",
  # #                        split.by = "pool")
  # 
  # png(paste0(directory,"/",treatment,"_sketch_UMAP_pools_nointegration.png"),
  #     width = 14, height = 14, units = "in", res = 400)
  # plot(p)
  # dev.off()
  # 
  # p =  DimPlot(input_merged, 
  #              reduction = "umap", 
  #              group.by = "orig.ident",
  #              #split.by = "pool",
  #              alpha = 0.8) 
  # 
  # # p = SCpubr::do_DimPlot(input_merged, 
  # #                        group.by = "orig.ident",
  # #                        split.by = "orig.ident")
  # 
  # png(paste0(directory,"/",treatment,"_sketch_UMAP_orig.ident_nointegration.png"),
  #     width = 14, height = 14, units = "in", res = 400)
  # plot(p)
  # dev.off()
  
  
  # integrate sketch
  
  # integrate the datasets
  input_merged = IntegrateLayers(input_merged, method = HarmonyIntegration, orig = "pca",
                                 new.reduction = "harmony",
                                 dims = 1:30, k.anchor = 20,
                                 verbose = TRUE)
  # cluster the integrated data
  input_merged = FindNeighbors(input_merged, reduction = "harmony", dims = 1:30)
  input_merged = FindClusters(input_merged, resolution = 0.3)
  input_merged = RunUMAP(input_merged, reduction = "harmony", dims = 1:30, return.model = T, verbose = F)
  
  ##### important ######
  # There is no need to regress out sample-specific attributes like percent.mito 
  # after integration, since batch correction takes care of this at the feature level.
  # see https://github.com/satijalab/seurat/issues/3579
  #####
  
  
  
  ########
  # p =  DimPlot(input_merged, 
  #              reduction = "umap", 
  #              group.by = "pool",
  #              #split.by = "pool", # throws error
  #              alpha = 0.8) 
  
  # p = SCpubr::do_DimPlot(input_merged, 
  #                        group.by = "pool",
  #                        split.by = "pool")
  
  # png(paste0(directory,"/",treatment,"_sketch_UMAP_pools_harmony_integration.png"),
  #     width = 14, height = 14, units = "in", res = 400)
  # plot(p)
  # dev.off()
  
  # integrate full object
  # beware the sketch Layers need to stay split
  input_merged[["sketch"]]
  
  input_merged = ProjectIntegration(object = input_merged, sketched.assay = "sketch", assay = "RNA", reduction = "harmony")
  
  
  input_merged = ProjectData(object = input_merged, sketched.assay = "sketch",
                             assay = "RNA", sketched.reduction = "harmony",
                             full.reduction = "harmony.full", dims = 1:30,
                             refdata = list(cluster_full = "seurat_clusters"))
  input_merged = RunUMAP(input_merged, reduction = "harmony.full", dims = 1:30, reduction.name = "umap.full",
                         reduction.key = "UMAPfull_")
  
  # p =  DimPlot(input_merged, 
  #              reduction = "umap.full", 
  #              #group.by = "pool",
  #              split.by = "pool", # throws error
  #              alpha = 0.8,
  #              shuffle = TRUE) 
  
  # p = SCpubr::do_DimPlot(input_merged, 
  #                        group.by = "pool",
  #                        split.by = "pool")
  
  # png(paste0(directory,"/",treatment,"_full_UMAP_pools_harmony_integration.png"),
  #     width = 14, height = 14, units = "in", res = 400)
  # plot(p)
  # dev.off()
  
  # join layers
  input_merged[["RNA"]] = JoinLayers(input_merged[["RNA"]])
  Layers(input_merged[["RNA"]])
  DefaultAssay(input_merged) = "RNA"
  
  
  input_merged = CellCycleScoring(input_merged,  
                                  s.features=cc.genes.updated.2019$s.genes, 
                                  g2m.features=cc.genes.updated.2019$g2m.genes)
  
  message("BPCells count matrix size: ", format(object.size(input_merged), units = "Mb"))
  write_matrix_dir(mat = input_merged[["RNA"]]$counts, 
                   dir = paste0(directory,"/",treatment,"_filtered_harmony"),
                   overwrite = TRUE)
  input_merged[["RNA"]]$counts = open_matrix_dir(dir = paste0(directory,"/",treatment,"_filtered_harmony"))
  
  saveRDS(
    object = input_merged,
    file = paste0(directory,"/",treatment,"_filtered_harmony/",treatment,"_filtered_harmony.Rds"))
  message("BPCells count matrix size for one-layer treatment: ",treatment,": ", format(object.size(input_merged), units = "Mb"))
  print(unique(input_merged$pool))
  
}


