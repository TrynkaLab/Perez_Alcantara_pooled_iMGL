#Inspect integration

#.libPaths(c("/software/teamtrynka/conda/seuratv5/lib/R/library")) # path to libraries
library(Seurat)
library(BPCells)
library(SeuratObject)
library(Signac)
library(patchwork)
library(tidyverse)
library(harmony)
# library(future)
library(SCpubr)
source("./helpers_seuratv5.R")
# set this option when analyzing large datasets
options(future.globals.maxSize = (70000 * 1024 ^ 2)) #(~70 Gb)
# for new seurat assay slots
options(Seurat.object.assay.version = "v5")

directory = "../../data/results/1.2.Inspect_integration"
dir.create(directory, recursive = T)

metadata = list()
for(treatment in c("untreated", "IFN", "LPS")) {
  input_merged = readRDS(paste0("../../data/results/1.QC_v5/",treatment,"_filtered_harmony/",treatment,"_filtered_harmony.Rds"))
  Layers(input_merged[["RNA"]])
  DefaultAssay(input_merged) = "RNA"
  Idents(input_merged) = "cluster_full"
  
  message("BPCells count matrix size for one-layer treatment: ",
          treatment,": ", format(object.size(input_merged), units = "Mb"))
  print(unique(input_merged$pool))
  print(unique(input_merged$treatment))
  input_merged # 862,274 cells pools 2-17 for untreated, 1,064,807 for IFN
  
  # check how pools distribute among clusters
  p1 =  SCpubr::do_BarPlot(input_merged, 
                           group.by = "pool", 
                           split.by = "cluster_full",
                           plot.title = "Number of cells per pool in each cluster",
                           position = "fill")
  
  p2 =  SCpubr::do_BarPlot(input_merged, 
                           group.by = "pool", 
                           split.by = "cluster_full",
                           plot.title = "Number of cells per pool in each cluster",
                           position = "stack")
  
  png(paste0(directory,"/",treatment,"_barplot_cells_per_pool_x_cluster.png"),
      width = 10, height = 16, units = "in", res = 400)
  plot((p1 / p2 ) +patchwork::plot_layout(guides="collect") & ggplot2::theme(legend.position = 'bottom'))
  dev.off()
  
  write.table(table(input_merged@meta.data$pool, input_merged@meta.data$cluster_full),
              paste0(directory,"/",treatment,"_table_poolxcluster_ncells.txt"))
  
  # barplot cluster x cell cycle phase
  p1 = SCpubr::do_BarPlot(input_merged,
                          group.by = "Phase",
                          split.by = "cluster_full",
                          plot.title = "Proportion of cells per phase in each cluster",
                          position = "fill")
  
  png(paste0(directory,"/",treatment,"_barplot_prop_per_phase_x_cluster.png"),
      width = 10, height = 16, units = "in", res = 400)
  p1
  dev.off()
  # untreated: cycling cells in cluster 3 and 4
  # LPS = 4, 5, 7
  # IFN = 4
  
  input_merged@meta.data$proliferation_status = "Not_proliferating"
  if(treatment == "untreated"){
    input_merged@meta.data[input_merged@meta.data$cluster_full %in%  c(3,4) , "proliferation_status"] =  "Proliferating"
    
  }
  if(treatment == "LPS"){
    input_merged@meta.data[input_merged@meta.data$cluster_full %in%  c(4,5,7) , "proliferation_status"] =  "Proliferating"
    
  }
  if(treatment == "IFN"){  # much smaller proliferating clusters on IFN treated cells
    input_merged@meta.data[input_merged@meta.data$cluster_full %in%  c(4) , "proliferation_status"] =  "Proliferating"
    
  }
  
  
  metadata[[treatment]] = input_merged@meta.data
  
}


metadata = do.call("rbind",metadata)
nrow(metadata) # 2,118,070
write_csv(metadata,paste0(directory,"/metadata.csv"))


for(treatment in c("untreated", "IFN", "LPS")){
# plot 50k subset
  input_subset = readRDS(
    file = paste0("../../data/results/1.1.subset_for_visualization/",treatment,"_filtered_harmony_subset/",treatment,"_filtered_harmony_50k_subset_seuratv5.Rds"))
  
  p = SCpubr::do_DimPlot(input_subset,
                         group.by = "pool",
                         split.by = "pool",
                         reduction="umap.full", # full or sketch (second takes less, but first shows some missing cells that group)
                         raster=TRUE,# helps with plotting many cells
                         pt.size = 3, # so the points are not so sparse
                         ncol = 4) 
  
  png(paste0(directory,"/",treatment,"_full_UMAP_pools_harmony_integration.png"),
      width = 16, height = 16, units = "in", res = 400)
  plot(p)
  dev.off()
  
  p = SCpubr::do_DimPlot(input_subset,
                         group.by = "Phase",
                         split.by = "cluster_full",
                         reduction="umap.full", # full or sketch (second takes less, but first shows some missing cells that group)
                         raster=TRUE,# helps with plotting many cells
                         pt.size = 3) # so the points are not so sparse
  
  png(paste0(directory,"/",treatment,"_full_UMAP_clusters_phase_harmony_integration.png"),
      width = 16, height = 16, units = "in", res = 400)
  plot(p)
  dev.off()
  
  p = UMAPS_all_libs(input_subset,raster.status = TRUE,pt.size.status = 2)
  
  png(paste0(directory,"/",treatment,"_full_UMAP_various_info_harmony_integration.png"),
      width = 16, height = 16, units = "in", res = 400)
  plot(p)
  dev.off()
  
  
  # check if lines from same pool cluster always near each other:
  lines_in_pools = read.table(paste0("../../../OTAR2065_scRNA_data_processing/data/","lines_in_pools.txt"), 
                              header = TRUE)
  for(p in unique(input_merged$pool)){
    lines = unlist(strsplit(lines_in_pools[lines_in_pools$Pool == p,"Lines"],split = ";"))
    seurat_subset = subset(input_merged, 
                           subset = donor_id %in% lines) 
    seurat_subset = subset(seurat_subset, 
                           subset = pool == p) 
    # plotting only categories over 10 cells, otherwise it throws error 
    seurat_subset = seurat_subset %>%
      subset(donor_id %in% names(which(table(seurat_subset$donor_id) >10)))
    seurat_subset@assays$sketch <- NULL # remove the sketch assay
    message("plotting plot 1")
    p1 = SCpubr::do_DimPlot(seurat_subset, 
                            group.by = "donor_id", 
                            split.by = "donor_id",
                            plot.title = paste0("Donor ids over 10 cells:" ,p),
                            reduction="umap.full")
    
    png(paste0(directory,"/",treatment,"_UMAP_",p,"_all_donors.png"),
        width = 32, height = 16, units = "in", res = 400)
    plot(p1)
    dev.off()
    
  }
  
  # find shared donors and plot distribution 
  
  sc_ncells = load_sc_donors(doublets_unassigned_in_proportion = FALSE,
                             # if TRUE, consider doublet and unassigned in final proportion 
                             # (final proportions for individual lines will be underestimated, particularly for donors with high proportions?)
                             # fix so doublets are assigned to their individual donors
                             pools = c(paste0("pool",2:17)),
                             directory =  "../../../OTAR2065_scRNA_data_processing/data/")
  
  
  # find shared donors 
  shared_donors = sc_ncells %>% 
    select(!c(set,Frequency,sample,n_lines,starting_prop,final_prop)) %>%
    unique() %>%
    group_by(Line) %>% 
    summarise(common=n()) %>% 
    filter(common>1) %>% 
    .$Line 
  
  # remove donors without cells
  
  shared_donors = shared_donors[shared_donors %in% unique(input_merged$donor_id)]
  
  subset_seurat = subset(input_merged, subset = donor_id %in% shared_donors)
  p1 =  SCpubr::do_BarPlot(subset_seurat, 
                           group.by = "donor_id", 
                           split.by = "cluster_full",
                           plot.title = "Proportion of cells per shared donor in each cluster",
                           position = "fill")
  
  p2 =  SCpubr::do_BarPlot(subset_seurat, 
                           group.by = "donor_id", 
                           split.by = "cluster_full",
                           plot.title = "Number of cells per shared donor in each cluster",
                           position = "stack")
  
  png(paste0(directory,"/",treatment,"_barplot_cells_per_pool_x_shared_donor.png"),
      width = 10, height = 16, units = "in", res = 400)
  plot((p1 / p2 ) +patchwork::plot_layout(guides="collect") & theme(legend.position = 'bottom'))
  dev.off()
  
  
  
  
  sc_ncells %>% 
    select(!c(set,Frequency,sample,n_lines,starting_prop,final_prop)) %>%
    unique() %>%
    filter(Line %in% c("civh_1", # pool9,10
                       "kuul_1")) # pool7, 9
  
  # frequency of donor doesn't seem to correlate with frequency of pool per cluster, but would probably need a proper correlation test
  # check if clones are together (more often in same cluster - should be)
  
  sc_ncells %>% 
    select(Line) %>% 
    unique() %>% 
    arrange(Line) %>% .$Line
  
  # clones so far are: zaie_1 and zaie_5, letw_1 and  letw_5 , lizq_1 and lizq_3, romx_2 and romx_1, seru_7 and seru_1, qonc_2 and qonc_1, 
  
  
  subset_seurat = subset(input_merged, subset = donor_id %in% c("zaie_1","zaie_5",
                                                                "letw_1",  "letw_5" ,
                                                                "lizq_1" , "lizq_3",
                                                                "romx_2", "romx_1",
                                                                "seru_7", "seru_1",
                                                                "qonc_2", "qonc_1",
                                                                "sebn_3","sebn_4"))
  p1 =  SCpubr::do_BarPlot(subset_seurat, 
                           group.by = "donor_id", 
                           split.by = "cluster_full",
                           plot.title = "Proportion of cells per clone in each cluster",
                           position = "fill")
  
  p2 =  SCpubr::do_BarPlot(subset_seurat, 
                           group.by = "donor_id", 
                           split.by = "cluster_full",
                           plot.title = "Number of cells per clone in each cluster",
                           position = "stack")
  
  png(paste0(directory,"/",treatment,"_barplot_cells_per_pool_x_donor_clones.png"),
      width = 10, height = 16, units = "in", res = 400)
  plot((p1 / p2 ) +patchwork::plot_layout(guides="collect") & theme(legend.position = 'bottom'))
  dev.off()
  
  
  
  # pseudobulk and get CPM
  # seurat's native AggregateExpression creates infinite counts using the normalised slot
  input_merged = JoinLayers(input_merged)
  pseudobulk = AggregateExpression(input_merged, 
                                   assays = "RNA", 
                                   slot = "count", # raw counts
                                   group.by = c("treatment"))
  pseudobulk = as.data.frame(pseudobulk)
  pseudobulk %>%
    mutate(CPM = ( RNA / sum(RNA)) * 10^6) %>%
    rownames_to_column( var = "gene") %>%
    dplyr::rename(raw_counts=RNA) %>%
    arrange(desc(raw_counts)) %>%
    write.table(paste0(directory,"/",treatment,"_pseudobulk_aggregated_raw_counts_CPMs.txt"),
                sep = "\t",
                row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  pseudobulk = AggregateExpression(input_merged, 
                                   assays = "RNA", 
                                   slot = "count", # raw counts
                                   group.by = c("treatment","donor_id"))
  pseudobulk = as.data.frame(pseudobulk)
  pseudobulk %>%
    mutate_all(~ (. / sum(.)) * 10^6) %>%
    rownames_to_column( var = "gene") %>%
    write.table(paste0(directory,"/",treatment,"_pseudobulk_aggregated_raw_counts_treatmentxdonor_CPMs.txt"),
                sep = "\t",
                row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  input_merged@meta.data$proliferation_status = "Not_proliferating"
  if(treatment == "untreated"){
    input_merged@meta.data[input_merged@meta.data$cluster_full %in%  c(5,6) , "proliferation_status"] =  "Proliferating"
    
  }
  if(treatment == "LPS"){
    input_merged@meta.data[input_merged@meta.data$cluster_full %in%  c(4,5,7) , "proliferation_status"] =  "Proliferating"
    
  }
  if(treatment == "IFN"){  # much smaller proliferating clusters on IFN treated cells
    input_merged@meta.data[input_merged@meta.data$cluster_full %in%  c(4) , "proliferation_status"] =  "Proliferating"
    
  }
  pseudobulk = AggregateExpression(input_merged, 
                                   assays = "RNA", 
                                   slot = "count", # raw counts
                                   group.by = c("treatment","proliferation_status","donor_id"))
  pseudobulk = as.data.frame(pseudobulk)
  pseudobulk %>%
    mutate_all(~ (. / sum(.)) * 10^6) %>%
    rownames_to_column( var = "gene") %>%
    write.table(paste0(directory,"/",treatment,"_pseudobulk_aggregated_raw_counts_treatmentxprolifxdonor_CPMs.txt"),
                sep = "\t",
                row.names = FALSE, col.names = TRUE, quote = FALSE)


# raw counts pseudobulk
  
  
  pseudobulk = AggregateExpression(input_merged, 
                                   assays = "RNA", 
                                   slot = "count", # raw counts
                                   group.by = c("treatment","donor_id"))
  pseudobulk = as.data.frame(pseudobulk)
  pseudobulk %>%
    rownames_to_column( var = "gene") %>%
    write.table(paste0(directory,"/",treatment,"_pseudobulk_aggregated_raw_counts_treatmentxdonor.txt"),
                sep = "\t",
                row.names = FALSE, col.names = TRUE, quote = FALSE)
  
 
  # pseudobulk = AggregateExpression(input_merged, 
  #                                  assays = "RNA", 
  #                                  slot = "count", # raw counts
  #                                  group.by = c("treatment","proliferation_status","donor_id"))
  # pseudobulk = as.data.frame(pseudobulk)
  # pseudobulk %>%
  #   rownames_to_column( var = "gene") %>%
  #   write.table(paste0(directory,"/",treatment,"_pseudobulk_aggregated_raw_counts_treatmentxprolifxdonor.txt"),
  #               sep = "\t",
  #               row.names = FALSE, col.names = TRUE, quote = FALSE)
  
    
  
}
