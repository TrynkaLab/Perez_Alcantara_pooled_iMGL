# 1.4. marker visualization

library(Seurat)
library(tidyverse)
library(SCpubr)
library(BPCells)
library(scCustomize)
library(SingleR)
library(patchwork)

source("./helpers.R")
# set this option when analyzing large datasets
options(future.globals.maxSize = 200000 * 1024^2) # 200Gb
treatment_cols =  c(untreated = "#8D918B", IFN = "#3A5683", LPS = "#F8766D")

# for new seurat assay slots
options(Seurat.object.assay.version = "v5")

directory = "../../data/results/1.4.markers"
dir.create(directory, recursive = T)


seu_subset = list()
for(treatment in c("untreated", "IFN", "LPS")) {
  message("Working on treatment ", treatment)
  seu_subset[[treatment]] = readRDS(paste0("../../data/results/1.QC_v5/",treatment,
                                           "_filtered_harmony/",treatment,"_filtered_harmony.Rds"))
  format(object.size(seu_subset[[treatment]]), units = "Mb")
  
  gc()
}

# merging clusters from the aggregate information provided later in the code
for(treatment in c("untreated", "IFN", "LPS")) {
  
  if(treatment == "untreated"){
    seu_subset[[treatment]]@meta.data = seu_subset[[treatment]]@meta.data %>%
      dplyr::mutate(cluster_full = as.numeric(cluster_full)) %>%
      # ordered by number of cells checked in the barplot below
      dplyr::mutate(merged_clusters = case_when(cluster_full %in% c(0) ~ "HM",
                                                cluster_full %in% c(1,2,5,6,8,9,10) ~ "DAM",
                                                cluster_full %in% c(3,4) ~ "Proliferating",
                                                cluster_full == 7 ~ "IFN"
      ))  %>%
      dplyr::mutate(merged_clusters = factor(merged_clusters,
                                             levels = c("Proliferating","IFN","DAM","HM"),
                                             ordered = TRUE)) %>%
      as.data.frame()
  }
  if(treatment == "IFN"){
    seu_subset[[treatment]]@meta.data = seu_subset[[treatment]]@meta.data %>%
      dplyr::mutate(cluster_full = as.numeric(cluster_full)) %>%
      dplyr::mutate(merged_clusters = case_when(cluster_full %in% c(0,1,2,3) ~ "DAM-IFN",
                                                cluster_full == 4 ~ "Proliferating",
                                                cluster_full == 5 ~ "HLA"))  %>%
      dplyr::mutate(merged_clusters = factor(merged_clusters,
                                             levels = c("Proliferating","HLA","DAM-IFN"),
                                             ordered = TRUE)) %>%
      as.data.frame()
  }
  if(treatment == "LPS"){
    seu_subset[[treatment]]@meta.data = seu_subset[[treatment]]@meta.data %>%
      dplyr::mutate(cluster_full = as.numeric(cluster_full)) %>%
      dplyr::mutate(merged_clusters = case_when(cluster_full %in% c(0,3) ~ "CRM",
                                                cluster_full == 1 ~ "DAM",
                                                cluster_full == 2 ~ "IFN",
                                                cluster_full %in% c(4,5,7) ~ "Proliferating",
                                                cluster_full %in% c(6,8) ~ "HLA"
      )) %>%
      dplyr::mutate(merged_clusters = factor(merged_clusters,
                                             levels = c("Proliferating","HLA","IFN","DAM","CRM"),
                                             ordered = TRUE)) %>%
      as.data.frame()
  }
  Idents(seu_subset[[treatment]]) = seu_subset[[treatment]]$merged_clusters
}

# markers from microglial subtypes
subtypes_Mancuso_full = readr::read_csv("../../../resources/Mancuso_2024/data/QC/top10_markers_with_pvM.csv") %>%
  group_by(cluster) %>%
  dplyr::slice_head(n=30) %>%
  summarise(gene = list(gene)) %>%
  ungroup()
subtypes_Mancuso_no_TRM = readr::read_csv("../../../resources/Mancuso_2024/data/QC/top10_markers_with_pvM.csv") %>%
  dplyr::filter(cluster!="Translational\n response (TRM)") %>%
  group_by(cluster) %>%
  dplyr::slice_head(n=30) %>%
  summarise(gene = list(gene)) %>%
  ungroup()

subtypes_Mancuso_no_TRM = subtypes_Mancuso_no_TRM$gene
names(subtypes_Mancuso_no_TRM) = c("HLA","CAM-CRM","PvM","CRM","DAM","HM","IRM","Neutrophils","Proliferating","Transit. CRM") # shorter for module scoring

subtypes_Mancuso_full = as.list(set_names(subtypes_Mancuso_full$gene, subtypes_Mancuso_full$cluster))


for(treatment in c("untreated", "IFN", "LPS")) {
  
  seu_subset[[treatment]] = AddModuleScore( seu_subset[[treatment]] ,
                                            features = subtypes_Mancuso_no_TRM,
                                            assay = "RNA",slot = "data",
                                            name = names(subtypes_Mancuso_no_TRM))
  
  
  write_matrix_dir(mat = seu_subset[[treatment]][["RNA"]]$counts,
                   dir = paste0(directory,"/",treatment,"_harmony_addedMancusoModules"),
                   overwrite = TRUE)
  seu_subset[[treatment]][["RNA"]]$counts = open_matrix_dir(dir = paste0(directory,"/",treatment,"_harmony_addedMancusoModules"))
  
  saveRDS(
    object = seu_subset[[treatment]],
    file = paste0(directory,"/",treatment,"_harmony_addedMancusoModules/",treatment,"_harmony_addedMancusoModules.Rds"))
  
  write_csv(seu_subset[[treatment]]@meta.data,paste0(directory,"/",treatment,"_module_annotated_metadata.csv"))
  
}

# read in to save time in case I need to replot
seu_subset = list()
for(treatment in c("untreated", "IFN", "LPS")) {
  
  seu_subset[[treatment]] = readRDS(
    file = paste0(directory,"/",treatment,"_harmony_addedMancusoModules/",treatment,"_harmony_addedMancusoModules.Rds"))
  
}


for(treatment in c("untreated", "IFN", "LPS")) {
  
  # merging 50k subsets doesn't work (I get error "object of type 'S4' is not subsettable")
  # maybe problem with subset object?
  # so need to do the analysis per treatment
  # p1 = scCustomize::DimPlot_scCustom(seu_subset[[treatment]] ,group.by = "pool",split.by = "pool")
  # 
  # pdf(paste0(directory,"/",treatment,"_umap_pools.pdf"),width = 15,height = 15)
  # plot(p1)
  # 
  # dev.off()
  
  p1 = scCustomize::DimPlot_scCustom(seu_subset[[treatment]] ,group.by = "cluster_full",split.by = "cluster_full")
  
  pdf(paste0(directory,"/",treatment,"_umap_cluster.pdf"),width = 15,height = 15)
  plot(p1)
  
  dev.off()
}

# with selected mancuso cluster genes

for(treatment in c("untreated", "IFN", "LPS")) {
  
  selected_gene_expr = list("General microglia" = c("CSF1R","CD74","C3"),
                            "Homeostatic"=c("P2RY12","CX3CR1"),
                            "Proliferating"=c("MKI67"  ,  "ASPM"  ,   "UBE2C" ,   "RRM2" ,    "DLGAP5"),
                            "CRM_merged" = c( "CCL2" ,  "IL1B", "CCL3" ,  "CCL4"),
                            "HLA" = c( "HLA-DRA" , "HLA-DRB5"),
                            "PvM" = c("F13A1", "CD209" ,"LYVE1"),
                            "DAM" = c("LGALS1"   , "CD9"  ,"PLA2G7" , "LGALS3" , "CXCR4"),
                            "IFN" = c("ISG15" , "IFIT1" , "IFIT3"),
                            "Neutrophils" = c("CCDC173"  , "FCER1A"  , "PTGS2"),
                            "Transitioning CRM" = c( "DUSP1" , "KLF2" ,  "JUN"  , "FOS"))
  
  seu_subset[[treatment]] = AddModuleScore( seu_subset[[treatment]] ,
                                            features = selected_gene_expr,
                                            assay = "RNA",slot = "data",
                                            name = names(selected_gene_expr))
  
  # rename
  seu_subset[[treatment]]@meta.data = seu_subset[[treatment]]@meta.data %>%
    dplyr::rename("General microglia" = "General microglia1",                   
                  "Homeostatic"= "Homeostatic2",
                  "Proliferating" = "Proliferating3",
                  "CRM" = "CRM_merged4",
                  "HLA" = "HLA5",
                  "PvM" = "PvM6",
                  "DAM" = "DAM7",         
                  "IFN" = "IFN8",
                  "Neutrophils" = "Neutrophils9",
                  "Transitioning CRM" = "Transitioning CRM10"
    ) %>%
    as.data.frame()
  
  for(treatment in c("untreated", "IFN", "LPS")) {
    
    # add proliferation info
    if(treatment == "untreated"){
      seu_subset[[treatment]]@meta.data[seu_subset[[treatment]]@meta.data$cluster_full %in%  c(3,4) , "proliferation_status"] =  "Proliferating"
      
    }
    if(treatment == "LPS"){
      seu_subset[[treatment]]@meta.data[seu_subset[[treatment]]@meta.data$cluster_full %in%  c(4,5,7) , "proliferation_status"] =  "Proliferating"
      
    }
    if(treatment == "IFN"){  # much smaller proliferating clusters on IFN treated cells
      seu_subset[[treatment]]@meta.data[seu_subset[[treatment]]@meta.data$cluster_full %in%  c(4) , "proliferation_status"] =  "Proliferating"
      
    }
  }
  
  ## fixing cluster 0 from untreated: majority of HM signal is very weak and concentrated at bottom of cluster 0
  seu_subset[["untreated"]]@meta.data = seu_subset[["untreated"]]@meta.data  %>%
    dplyr::mutate(merged_clusters = case_when(cluster_full == 0 & Homeostatic > 1 ~ "HM", 
                                              cluster_full == 0 & Homeostatic < 1 ~ "DAM",
                                              .default = merged_clusters)) %>%
    as.data.frame()
  
  # some HLA as well
  seu_subset[["untreated"]]@meta.data = seu_subset[["untreated"]]@meta.data  %>%
    dplyr::mutate(merged_clusters = case_when(cluster_full == 0 & HLA > 2 ~ "HLA", 
                                              .default = merged_clusters)) %>%
    as.data.frame()
  
  for(treatment in c("untreated","IFN","LPS")){
    if(treatment == "untreated"){
      seu_subset[[treatment]]@meta.data = seu_subset[[treatment]]@meta.data %>%
        dplyr::mutate(merged_clusters = factor(merged_clusters,
                                               levels = c("Proliferating","HLA","IFN","DAM","HM"),
                                               ordered = TRUE)) %>%
        as.data.frame()
    }
    if(treatment == "IFN"){
      seu_subset[[treatment]]@meta.data = seu_subset[[treatment]]@meta.data %>%
        dplyr::mutate(merged_clusters = factor(merged_clusters,
                                               levels = c("Proliferating","HLA","DAM-IFN"),
                                               ordered = TRUE)) %>%
        as.data.frame()
    }
    if(treatment == "LPS"){
      seu_subset[[treatment]]@meta.data = seu_subset[[treatment]]@meta.data %>%
        dplyr::mutate(merged_clusters = factor(merged_clusters,
                                               levels = c("Proliferating","HLA","IFN","DAM","CRM"),
                                               ordered = TRUE)) %>%
        as.data.frame()
    }
    Idents(seu_subset[[treatment]]) = seu_subset[[treatment]]$merged_clusters
    
  }
  
  ### pools by cluster
  for(treatment in c("untreated", "IFN", "LPS")) {
    
    p1 = SCpubr::do_BarPlot(seu_subset[[treatment]] ,
                            group.by = "pool",split.by = "merged_clusters",
                            position = "fill",
                            flip = FALSE,
                            add.n = TRUE,
                            add.n.size = 4,)
    
    pdf(paste0(directory,"/",treatment,"_barplot_pools_by_merged_cluster.pdf"),width = 7,height = 7)
    plot(p1)
    
    dev.off()
    
  }
  
  ### proportion of cells in each merged category
  ## venn diagram
  metadata_df = list()
  for(treatment in c("untreated", "IFN", "LPS")) {
    
    metadata_df[[treatment]] = seu_subset[[treatment]]@meta.data %>%
      dplyr::select(cell,merged_clusters) %>%
      dplyr::mutate(treatment = treatment)
  }
  metadata_df = do.call("rbind",metadata_df)
  
  # Summarize proportions
  df_summary = metadata_df %>%
    dplyr::group_by(treatment) %>%
    dplyr::count(merged_clusters, name = "count") %>%
    dplyr::mutate(proportion = count / sum(count)) %>%
    dplyr::ungroup()
  
  df_summary %>%
    readr::write_csv(paste0(directory,"/pie_proportions.csv"))
  # treatment merged_clusters  count proportion
  # <chr>     <ord>            <int>      <dbl>
  #   1 IFN       Proliferating     4517   0.00735 
  # 2 IFN       HLA                 64   0.000104
  # 3 IFN       DAM-IFN         610391   0.993   
  # 4 LPS       Proliferating    69604   0.109   
  # 5 LPS       HLA               7487   0.0117  
  # 6 LPS       IFN              57592   0.0899  
  # 7 LPS       DAM              64264   0.100   
  # 8 LPS       CRM             441877   0.690   
  # 9 untreated Proliferating    93933   0.109   
  # 10 untreated HLA               7011   0.00813 
  # 11 untreated IFN               6437   0.00747 
  # 12 untreated DAM             734796   0.852   
  # 13 untreated HM               20097   0.0233 
  # Plot pie chart
  p = ggplot(df_summary, aes(x = "", y = proportion, 
                             fill = merged_clusters, label = merged_clusters)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    theme_void() +
    facet_wrap(vars(treatment)) + 
    theme(legend.position = "top")
  
  pdf(paste0(directory,"/pie_merged_cluster.pdf"),
      width = 7,height = 3)
  plot(p)
  
  dev.off()
  
  for(treatment in c("untreated", "IFN", "LPS")) {
    
    #  The following feature lists do not have enough features present in the object: Translational response (TRM)
    # that's why I removed it above
    
    # scCustomize::Iterate_FeaturePlot_scCustom(seurat_object = seu_subset[[treatment]] ,
    #                                          features = paste0(names(subtypes_Mancuso_no_TRM),1:10),
    #                                          file_name =  paste0(directory,"/",treatment,"_umap_mancuso_genes"),
    #                                          single_pdf = T,features_per_page = 10,num_columns = 4,na_cutoff = NA,
    #                                          raster = TRUE,order = TRUE
    #)
    
    
    my_features = paste0(names(subtypes_Mancuso_no_TRM),1:10)
    
    # Retrieve your scores.
    scores = seu_subset[[treatment]] @meta.data[, my_features]
    
    # Turn it into a matrix (metadata columns become the features - rows, cells remain as columns).
    scores = t(as.matrix(scores))
    
    # Create an assay object.
    scores_assay = Seurat::CreateAssayObject(counts = scores)
    
    # Add it to a Seurat object (make copy so it doesn't complain on saving later)
    copy_seurat = seu_subset[[treatment]]
    copy_seurat@assays$mancuso = scores_assay
    
    # Plot the features.
    p = SCpubr::do_ExpressionHeatmap(sample = copy_seurat ,
                                     features = my_features,
                                     assay = "mancuso",cluster = TRUE)
    
    pdf(paste0(directory,"/",treatment,"_heatmap_mancuso.pdf"),width = 5,height = 5)
    
    plot(p)
    
    dev.off()
    
    
    medians = Median_Stats(seurat_object = seu_subset[[treatment]], group_by_var = "merged_clusters",
                           median_var  =  paste0(names(subtypes_Mancuso_no_TRM),1:10))
    write_csv(medians,paste0(directory,"/",treatment,"_medians_per_merged_cluster_mancuso_genes.csv"))
    medians = Median_Stats(seurat_object = seu_subset[[treatment]], group_by_var = "pool",
                           median_var  =  paste0(names(subtypes_Mancuso_no_TRM),1:10))
    write_csv(medians,paste0(directory,"/",treatment,"_medians_per_pool_mancuso_genes.csv"))
    
  }
  
  
  
  for(treatment in c("untreated", "IFN", "LPS")) {
    
    
    write_matrix_dir(mat = seu_subset[[treatment]][["RNA"]]$counts,
                     dir = paste0(directory,"/",treatment,"_harmony_annotated"),
                     overwrite = TRUE)
    seu_subset[[treatment]][["RNA"]]$counts = open_matrix_dir(dir = paste0(directory,"/",treatment,"_harmony_annotated"))
    
    
    saveRDS(
      object = seu_subset[[treatment]],
      file = paste0(directory,"/",treatment,"_harmony_annotated/",treatment,"_harmony_annotated.Rds"))
    
    write_csv(seu_subset[[treatment]]@meta.data,paste0(directory,"/",treatment,"_harmony_annotated_metadata.csv"))
  }
  
  # load in
  
  for(treatment in c("untreated", "IFN", "LPS")) {
    seu_subset[[treatment]] = readRDS(paste0(directory,"/",treatment,"_harmony_annotated/",treatment,"_harmony_annotated.Rds"))
  }
  
  for(treatment in c("untreated", "IFN", "LPS")) {
    
    # Retrieve your scores.
    scores = seu_subset[[treatment]] @meta.data[, c("General microglia", "Homeostatic", "Proliferating" ,
                                                    "CRM" , "HLA","PvM" ,"DAM","IFN","Neutrophils",  "Transitioning CRM"  )]
    
    # Turn it into a matrix (metadata columns become the features - rows, cells remain as columns).
    scores = t(as.matrix(scores))
    
    # Create an assay object.
    scores_assay = Seurat::CreateAssayObject(counts = scores,
                                             key = "mancuso_selected")
    
    # Add it to a Seurat object (make copy so it doesn't complain on saving later)
    copy_seurat = seu_subset[[treatment]]
    copy_seurat@assays$mancuso_selected = scores_assay
    
    # Plot the features.
    p = SCpubr::do_ExpressionHeatmap(sample = copy_seurat ,
                                     features =   c("General microglia", "Homeostatic", "Proliferating" ,
                                                    "CRM" , "HLA","PvM" ,"DAM","IFN","Neutrophils",  "Transitioning CRM"  ),
                                     assay = "mancuso_selected",cluster = TRUE)
    
    pdf(paste0(directory,"/",treatment,"_heatmap_mancuso_selected_genes.pdf"),width = 5,height = 5)
    
    plot(p)
    
    dev.off()
    
    
    
  }
  
  for(treatment in c("untreated", "IFN", "LPS")) {
    
    # p1 =  scCustomize::DotPlot_scCustom(seurat_object = seu_subset[[treatment]],
    #                                     features = c( "CSF1R","CD74","C3", # General microglia
    #                                                          "P2RY12" ,  "CX3CR1" , # homeostatic microglia
    #                                                          "MKI67"  ,  "ASPM"  ,   "UBE2C" ,   "RRM2" ,    "DLGAP5" , # Proliferating
    #                                                          "CCL2" ,  "IL1B", "CCL3" ,  "CCL4", # CRM
    #                                                          "HLA-DRA" , "HLA-DRB5" , # antigen-presenting response
    #                                                          "F13A1", "CD209" ,"LYVE1", # pvM / CNS-associated macrophage
    #                                                          "LGALS1"   , "CD9"  ,"PLA2G7" , "LGALS3" , "CXCR4", # disease-associated microglia
    #                                                          "ISG15" , "IFIT1" , "IFIT3", # interferon response
    #                                                          "CCDC173"  , "FCER1A"  , "PTGS2", # Neutrophils
    #                                                          "DUSP1" , "KLF2" ,  "JUN"  , "FOS", # Transitioning CRM
    #                                                          "RPS4X", "RPL3" , "RPL29" # Translational response
    #                                                          ),
    #               cluster.idents = TRUE) +
    #   theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust = 1))
    
    
    
    p1 =  scCustomize::DotPlot_scCustom(seurat_object = seu_subset[[treatment]],
                                        features = c( "TREM2","ITGAX","INPP5D", # Sam phagocytosis
                                                      "AXL" ,  "SIRPA" , # phago lively and schlichter
                                                      
                                                      "NOX1", "NOX4" , # ROS lively
                                                      "P2RX7" , "P2RY2" , "P2RY6", "P2RY12"   # migration lively
                                                      
                                        ),
                                        cluster.idents = TRUE) +
      theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust = 1))
    
    pdf(paste0(directory,"/",treatment,"_dotplot_phenotype_genes.pdf"),width = 7,height = 4)
    plot(p1)
    
    dev.off()
    
    # plotting features from lively and schlichter
    scCustomize::Iterate_FeaturePlot_scCustom(seurat_object = seu_subset[[treatment]],
                                              features =  c( "TREM2","ITGAX","INPP5D", # Sam phagocytosis
                                                             "AXL" ,  "SIRPA" , # phago lively and schlichter
                                                             
                                                             "NOX1", "NOX4" , # ROS lively
                                                             "P2RX7" , "P2RY2" , "P2RY6", "P2RY12"   # migration lively
                                                             
                                              ),file_name =  paste0(directory,"/",treatment,"_umap_phenotype_genes"),
                                              single_pdf = T,features_per_page = 11,num_columns = 4,
                                              raster = TRUE,order = TRUE
    )
    
    
    p1 = scCustomize::FeaturePlot_scCustom(seurat_object = seu_subset[[treatment]], features = "TREM2", split.by = "pool",
                                           num_columns = 4,order = TRUE,raster = TRUE,layer = "data",combine = TRUE)
    pdf(paste0(directory,"/",treatment,"_umap_TREM2_pools.pdf"),width = 15,height = 4)
    
    plot(p1)
    dev.off()
    
    
    # seu_subset[[treatment]] = NA
    gc()
    
    #### add LPS/Ifng response genes
    # try merging before adding module score
    
  }
  
  # subset of genes used for initial annotation in Mancuso's paper - some merged categories
  
  
  selected_gene_expr = list("General microglia" = c("CSF1R","CD74","C3"),
                            "Homeostatic"=c("P2RY12","CX3CR1"),
                            "Proliferating"=c("MKI67"  ,  "ASPM"  ,   "UBE2C" ,   "RRM2" ,    "DLGAP5"),
                            "CRM" = c( "CCL2" ,  "IL1B", "CCL3" ,  "CCL4"),
                            "HLA" = c( "HLA-DRA" , "HLA-DRB5"),
                            "PvM" = c("F13A1", "CD209" ,"LYVE1"),
                            "DAM" = c("LGALS1"   , "CD9"  ,"PLA2G7" , "LGALS3" , "CXCR4"),
                            "IFN" = c("ISG15" , "IFIT1" , "IFIT3"),
                            "Neutrophils" = c("CCDC173"  , "FCER1A"  , "PTGS2"),
                            "Transitioning CRM" = c( "DUSP1" , "KLF2" ,  "JUN"  , "FOS"))
  
  p = list()
  for(treatment in c("untreated","IFN","LPS")){
    # Idents(seu_subset[[treatment]]) = "merged_clusters"
    
    p[[treatment]] = scCustomize::DotPlot_scCustom(seurat_object = seu_subset[[treatment]],
                                                   group.by = "merged_clusters",
                                                   features = c( "CSF1R","CD74","C3", # General microglia
                                                                 "P2RY12" ,  "CX3CR1" , # homeostatic microglia
                                                                 "LGALS1"   , "CD9"  ,"PLA2G7" , "LGALS3" , "CXCR4", # disease-associated microglia
                                                                 
                                                                 "CCL2" ,  "IL1B", "CCL3" ,  "CCL4", # CRM
                                                                 "DUSP1" , "KLF2" ,  "JUN"  , "FOS", # Transitioning CRM
                                                                 "HLA-DRA" , "HLA-DRB5" , # antigen-presenting response
                                                                 "ISG15" , "IFIT1" , "IFIT3", # interferon response
                                                                 "F13A1", "CD209" ,"LYVE1", # pvM / CNS-associated macrophage
                                                                 
                                                                 "CCDC173"  , "FCER1A"  , "PTGS2", # Neutrophils
                                                                 "MKI67"  ,  "ASPM"  ,   "UBE2C" ,   "RRM2" ,    "DLGAP5"  # Proliferating
                                                                 
                                                   ),
                                                   cluster.idents = FALSE, assay = "sketch",
                                                   # colors_use = viridis_inferno_dark_high,
                                                   dot.scale =3.5,
                                                   col.min = -1, col.max = 1) +
      scale_color_gradient2(low = "#047098",
                            mid = "white",
                            high = "#B20505") +
      ggpubr::theme_pubr() +
      ylab(treatment) 
    
    
  }
  
  pdf(paste0(directory,"/all_treatments_dotplot_Mancuso_categories.pdf"),width = 5.3,height = 4,pointsize = 5)
  patchwork::wrap_plots(p,ncol = 1,heights = c(4,3,5)) +
    patchwork::plot_layout(axes = "collect", guides = "collect") &
    theme(legend.position = "top",
          legend.key.size =  unit(1,"line"),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 7),
          
          axis.title.y = element_text(size = 7),
          axis.title.x = element_blank(),
          
          axis.text.x = element_text(size = 7,
                                     angle = 90,
                                     hjust = 1,
                                     vjust = 0.5),
          axis.text.y = element_text(size = 6))                                                                    
  
  dev.off()
  
  
  ###### create merged seurat object
  
  full_seurat = merge(seu_subset$untreated,y=list(seu_subset$IFN,seu_subset$LPS),
                      merge.dr = FALSE, merge.data = TRUE)
  
  # working on sketch
  
  DefaultAssay(full_seurat) = "sketch"
  full_seurat = FindVariableFeatures(full_seurat)
  full_seurat = ScaleData(full_seurat)
  full_seurat = RunPCA(full_seurat,verbose = TRUE) 
  full_seurat = FindNeighbors(full_seurat, dims = 1:50,verbose = TRUE)
  full_seurat = FindClusters(full_seurat, resolution = 2,verbose = TRUE)
  full_seurat = RunUMAP(full_seurat, dims = 1:50, return.model = T,verbose = TRUE)
  
  colnames(full_seurat@meta.data)
  
  full_seurat[["sketch"]] = JoinLayers(full_seurat[["sketch"]])
  Layers(full_seurat[["sketch"]])
  DefaultAssay(full_seurat) = "sketch"
  # 
  # full_seurat
  
  # full_seurat = AddModuleScore( full_seurat ,
  #                                           features = subtypes_Mancuso_no_TRM,
  #                                           assay = "RNA",slot = "data",
  #                                           name = names(subtypes_Mancuso_no_TRM))
  
  message("BPCells count matrix size: ", format(object.size(full_seurat), units = "Mb"))
  write_matrix_dir(mat = full_seurat[["RNA"]]$counts,
                   dir = paste0(directory,"/merged_sketch"),
                   overwrite = TRUE)
  full_seurat[["RNA"]]$counts = open_matrix_dir(dir = paste0(directory,"/merged_sketch"))
  
  saveRDS(
    object = full_seurat,
    file = paste0(directory,"/merged_sketch/merged_sketch.Rds"))
  
  
  #### dimplots
  full_seurat = readRDS("../../data/results/1.4.markers/merged_sketch/merged_sketch.Rds")
  
  seurat_merged_subset = subset(x = full_seurat, downsample = 10000)
  
  # p = SCpubr::do_DimPlot(seurat_merged_subset,
  #                        label = TRUE,
  #                        reduction = "umap",
  #                        group.by = "merged_clusters",
  #                        raster = TRUE)
  # plot(p)
  # 
  # p = SCpubr::do_DimPlot(seurat_merged_subset,
  #                        label = TRUE,
  #                        reduction = "umap",
  #                        group.by = "merged_clusters",
  #                        split.by = "treatment",
  #                        
  #                        raster = TRUE)
  
  FeaturePlot(
    object = full_seurat,
    features = "DAM"
  )
  
  p= SCpubr::do_FeaturePlot(full_seurat,
                            features = "DAM",
                            reduction = "umap",
                            raster = TRUE,
                            plot.title = "DAM",
                            min.cutoff = -1,
                            max.cutoff = 1,
                            use_viridis = TRUE,
                            viridis.direction = -1,
                            viridis.palette = "viridis")
  p
  p2 = SCpubr::do_FeaturePlot(full_seurat,
                         features = "CRM",
                         reduction = "umap",
                         raster = TRUE,
                         plot.title = "CRM",
                         min.cutoff = -1,
                         max.cutoff = 1,
                         use_viridis = TRUE,
                         viridis.direction = -1,
                         viridis.palette = "viridis")
  p3 = SCpubr::do_FeaturePlot(full_seurat,
                              features = "IFN",
                              reduction = "umap",
                              raster = TRUE,
                              plot.title = "IFN",
                              min.cutoff = -1,
                              max.cutoff = 1,
                              use_viridis = TRUE,
                              viridis.direction = -1,
                              viridis.palette = "viridis")
  p4 = SCpubr::do_FeaturePlot(full_seurat,
                         features = "HLA",
                         reduction = "umap",
                         raster = TRUE,
                         plot.title = "HLA",
                         min.cutoff = -1,
                         max.cutoff = 1,
                         use_viridis = TRUE,
                         viridis.direction = -1,
                         viridis.palette = "viridis")
  
(p + p2 ) / (p3 + p4) + patchwork::plot_layout(guides = "collect")
  
  
  
  
  # subset of genes used for initial annotation in Mancuso's paper
  selected_gene_expr = list("Pluripotency" = c("POU5F1","SOX2"),
                            "General_microglia" = c("CSF1R","CD74","C3"),
                            "Homeostatic"=c("P2RY12","CX3CR1"),
                            "Brain_macrophages_phagocytic"=c("MRC1","CD163","LYVE1","F13A1"),
                            "Neuronal_surveillance" = c("FRMD4A",
                                                        # "GRCK3", # not found
                                                        # "SLC6A3", # 0
                                                        "ASMTL"),
                            "Ribosome_biogenesis" = c("FTL","FTH1"), # both ferritin chains?
                            "Lipid_processing" = c("MYO1E","PTPRG"),
                            "Stress" = c("HSP90AA1","HSPH1"),
                            "Glycolytic" = c("NAMPT","SLC2A3"),
                            "Inflammatory_I" = c("INO80D","TMEM163","CPEB4","IL4R"),
                            "Inflammatory_II" = c("SPON1","LRRK2","FOXP1"),
                            "Inflammatory_III" = c("CCL3","IL1B","RELB"),
                            "Antiviral" = c("IFI44L1","MX1"),
                            "Cycling" = c("EZH2","BRIP1")
                            
  )
  pdf(paste0(directory,"/subtypes_Mancuso_top_",treatment,"_gene_expr.pdf"),width = 10,height = 10)
  
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
    
    p + patchwork::plot_annotation(title = nam)
    
    plot(p)
  }
  
  dev.off()
  
  pdf(paste0(directory,"/subtypes_Mancuso_",treatment,"_module_scores.pdf"),width = 10,height = 10)
  
  for(columns in  paste0(names(subtypes_Mancuso_full),c(1:12))){
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
  
  # scaled module metadata heatmap
  # merged dataset
  scaled_module_score = full_seurat@meta.data %>%
    dplyr::select("cluster_full",
                  "treatment",
                  "Antiviral1",    "Cycling2"  ,   "Glycolytic3"  ,    "Homeostatic4"    ,
                  "Inflammatory_15"  ,   "Inflammatory_26"  ,  "Inflammatory_37",
                  "Lipid_processing8"   ,  "Neuronal_surveillance9" ,    "Phagocytic10"  ,   "Ribosome_biogenesis11",
                  "Stress_signature12"
    ) %>%
    # scaling across the three treatments per set of modules
    # to highlight differences between treatments and clusters
    dplyr::mutate_at(vars(-cluster_full,-treatment), function(x) scale(x,center = TRUE,scale = TRUE)) %>% # scaling all but cluster
    dplyr::as_tibble() %>%
    dplyr::rename( "cluster" = "cluster_full",
                   "Antiviral" = "Antiviral1",
                   "Cycling" ="Cycling2"  ,
                   "Glycolytic"  = "Glycolytic3"  ,
                   "Homeostatic"= "Homeostatic4"    ,
                   "Inflammatory_1" = "Inflammatory_15"  ,
                   "Inflammatory_2"  = "Inflammatory_26"  ,
                   "Inflammatory_3" = "Inflammatory_37",
                   "Lipid_processing"  = "Lipid_processing8"   ,
                   "Neuronal_surveillance" = "Neuronal_surveillance9" ,
                   "Phagocytic" = "Phagocytic10"  ,
                   "Ribosome_biogenesis" = "Ribosome_biogenesis11",
                   "Stress_signature" = "Stress_signature12") %>%
    dplyr::group_by(cluster,treatment) %>%
    dplyr::summarise_at(vars(-group_cols()), colMeans) %>%
    dplyr::mutate(cluster_by_treat = paste(cluster,treatment,sep = "_")) %>%
    column_to_rownames(var = "cluster_by_treat")
  
  
  
  annotation_r = data.frame(
    cluster = scaled_module_score$cluster,
    treatment = scaled_module_score$treatment
  )
  rownames(annotation_r) =rownames(scaled_module_score)
  
  scaled_module_score = scaled_module_score %>%
    dplyr::select(-cluster,-treatment)
  
  ann_color = list("treatment"=c("untreated"="black","LPS"="yellow","IFN"="red"),
                   "cluster" = viridis::viridis(10))
  names(ann_color$cluster) = 0:9
  
  p = pheatmap::pheatmap(scaled_module_score,
                         annotation_row = annotation_r,
                         annotation_colors = ann_color)
  
  pdf(paste0(directory,"/subtypes_Mancuso_top30_merged_module_scores_heatmap_by_cluster_treatment.pdf"),
      width = 10,height = 10)
  
  p
  dev.off()
  
  # per treatment
  
  for(treatment in c("untreated","IFN","LPS")){
    scaled_module_score = seu_subset[[treatment]]@meta.data %>%
      dplyr::select("cluster_full",
                    "treatment",
                    "Antiviral1",    "Cycling2"  ,   "Glycolytic3"  ,    "Homeostatic4"    ,
                    "Inflammatory_15"  ,   "Inflammatory_26"  ,  "Inflammatory_37",
                    "Lipid_processing8"   ,  "Neuronal_surveillance9" ,    "Phagocytic10"  ,   "Ribosome_biogenesis11",
                    "Stress_signature12"
      ) %>%
      # scaling across the three treatments per set of modules
      # to highlight differences between treatments and clusters
      dplyr::mutate_at(vars(-cluster_full,-treatment), function(x) scale(x,center = TRUE,scale = TRUE)) %>% # scaling all but cluster
      dplyr::as_tibble() %>%
      dplyr::rename( "cluster" = "cluster_full",
                     "Antiviral" = "Antiviral1",
                     "Cycling" ="Cycling2"  ,
                     "Glycolytic"  = "Glycolytic3"  ,
                     "Homeostatic"= "Homeostatic4"    ,
                     "Inflammatory_1" = "Inflammatory_15"  ,
                     "Inflammatory_2"  = "Inflammatory_26"  ,
                     "Inflammatory_3" = "Inflammatory_37",
                     "Lipid_processing"  = "Lipid_processing8"   ,
                     "Neuronal_surveillance" = "Neuronal_surveillance9" ,
                     "Phagocytic" = "Phagocytic10"  ,
                     "Ribosome_biogenesis" = "Ribosome_biogenesis11",
                     "Stress_signature" = "Stress_signature12") %>%
      dplyr::group_by(cluster,treatment) %>%
      dplyr::summarise_at(vars(-group_cols()), colMeans) %>%
      dplyr::mutate(cluster_by_treat = paste(cluster,treatment,sep = "_")) %>%
      column_to_rownames(var = "cluster_by_treat")
    
    
    
    annotation_r = data.frame(
      cluster = scaled_module_score$cluster,
      treatment = scaled_module_score$treatment
    )
    rownames(annotation_r) =rownames(scaled_module_score)
    
    scaled_module_score = scaled_module_score %>%
      dplyr::select(-cluster,-treatment)
    
    ann_color = list("treatment"=c("untreated"="black","LPS"="yellow","IFN"="red"),
                     "cluster" = viridis::viridis(10))
    names(ann_color$cluster) = 0:9
    
    p = pheatmap::pheatmap(scaled_module_score,
                           annotation_row = annotation_r,
                           annotation_colors = ann_color)
    
    pdf(paste0(directory,"/subtypes_Mancuso_top30_",treatment,"_module_scores_heatmap_by_cluster_treatment.pdf"),
        width = 10,height = 10)
    
    p
    dev.off()
    
  }
  all.markers = list()
  for(treatment in c("untreated","IFN","LPS")){
    
    Idents(seu_subset[[treatment]]) = seu_subset[[treatment]]$cluster_full
    seurobj = subset(x = seu_subset[[treatment]], downsample = 10000)
    seurobj[["RNA"]]$data = as(object = seurobj[["RNA"]]$data, Class = "dgCMatrix") # workaround FildAllMarkers bug for this version
    all.markers[[treatment]] = tibble::tibble(Seurat::FindAllMarkers(object = seurobj,
                                                                     assay = "RNA",
                                                                     logfc.threshold = 0.6)) # log2FC: 2^0.6 ~ 1.5x
    
    
  }
  rm(seurobj)
  gc()
  
  saveRDS(all.markers,
          file = paste0(directory,"/findAllMarkers_results.Rds"))
  
  all.markers = readRDS(paste0(directory,"/findAllMarkers_results.Rds"))
  for(treatment in c("untreated","IFN","LPS")){
    
    p = SCpubr::do_GroupwiseDEPlot(sample = seu_subset[[treatment]],
                                   de_genes = all.markers[[treatment]],
                                   top_genes = 10)
    
    
    
    pdf(paste0(directory,"/findallmarkers_",treatment,".pdf"),
        width = 10,height = 10)
    
    p
    dev.off()
    
  }
  
  
  # do find markers between clusters, then GO BP and enrichments
  # re-run plots with module scores from full gene sets
  
  # Use AUCell scoring, similar to seurat addmodulescore
  p <- SCpubr::do_EnrichmentHeatmap(sample = sample,
                                    input_gene_list = genes,
                                    flip = TRUE,
                                    cluster_cols = FALSE,
                                    cluster_rows = TRUE,
                                    flavor = "AUCell",
                                    viridis_direction = -1)
  p
  