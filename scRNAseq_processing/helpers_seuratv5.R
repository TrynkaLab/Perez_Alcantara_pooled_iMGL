cellranger_to_seurat = function(cellranger_df){
  seurat_list = vector("list", length(cellranger_df))
  names(seurat_list) = names(cellranger_df)
  for(n in names(cellranger_df)){
    message("Working on sample ", n)
    seurat_list[[n]] = CreateSeuratObject(counts = cellranger_df[[n]] , project = n)
    print(paste0("Median UMI counts are ", median(seurat_list[[n]]$nCount_RNA))) # Median UMI counts - as cellRanger summary
    print(paste0("Median genes per cell are ",median(seurat_list[[n]]$nFeature_RNA))) # Median genes per cell - as cellRanger summary
    seurat_list[[n]][["pool"]] = stringr::str_split(n,pattern = "\\.")[[1]][1]
    seurat_list[[n]][["sample"]] = stringr::str_split(n,pattern = "\\.")[[1]][2]
    if(grepl("plst|untreated|Untreated|AG|YC|ITMG", n)){
      seurat_list[[n]][["treatment"]] = "untreated"
      
    }
    if(grepl("IFN", n)){
      seurat_list[[n]][["treatment"]] = "IFN"
      
    }
    if(grepl("LPS", n)){
      seurat_list[[n]][["treatment"]] = "LPS"
      
    }
    if(grepl("ITMG", n)){
      seurat_list[[n]][["sequencing_batch"]] = "batch1"
      
    }
    if(grepl("AG|YC", n)){
      seurat_list[[n]][["sequencing_batch"]] = "batch2"
    }
    if(grepl("RGM", n)){
      seurat_list[[n]][["sequencing_batch"]] = "batch3"
    }
    if(unique(seurat_list[[n]][["pool"]]) == "pool3" | unique(seurat_list[[n]][["pool"]]) == "pool4"){
      seurat_list[[n]][["sequencing_batch"]] = "batch4"
      
    }
    if(unique(seurat_list[[n]][["pool"]]) == "pool5" | unique(seurat_list[[n]][["pool"]]) == "pool6"){
      seurat_list[[n]][["sequencing_batch"]] = "batch5"
      
    }
    if(unique(seurat_list[[n]][["pool"]]) ==  "pool7" | unique(seurat_list[[n]][["pool"]]) == "pool8" ){
      seurat_list[[n]][["sequencing_batch"]] = "batch6"
      
    }
    # fix batch info for pools 9 and 10
    if(unique(seurat_list[[n]][["pool"]]) ==  "pool9" | unique(seurat_list[[n]][["pool"]]) == "pool10" ){
      seurat_list[[n]][["sequencing_batch"]] = "batch7"
      
    }
    if(unique(seurat_list[[n]][["pool"]]) ==  "pool11"){
      seurat_list[[n]][["sequencing_batch"]] = "batch8"
      
    }
    
    if(unique(seurat_list[[n]][["pool"]]) ==  "pool13"){
      seurat_list[[n]][["sequencing_batch"]] = "batch9"
      
    }
    if(unique(seurat_list[[n]][["pool"]]) ==  "pool12"){
      seurat_list[[n]][["sequencing_batch"]] = "batch10"
      
    }
    if(unique(seurat_list[[n]][["pool"]]) ==  "pool14"){
      seurat_list[[n]][["sequencing_batch"]] = "batch11"
      
    }
    
    if(unique(seurat_list[[n]][["pool"]]) ==  "pool15"){
      seurat_list[[n]][["sequencing_batch"]] = "batch12"
      
    }
      if(unique(seurat_list[[n]][["pool"]]) ==  "pool16"){
      seurat_list[[n]][["sequencing_batch"]] = "batch13"
      
    }
          if(unique(seurat_list[[n]][["pool"]]) ==  "pool17"){
      seurat_list[[n]][["sequencing_batch"]] = "batch14"
      
    }
  }
  return(seurat_list)
}

incorporate_vireo_info = function(seurat_list,
                                  pools = paste0("pool",2:17),
                                  directory =  "../../../OTAR2065_scRNA_data_processing/data/"){
  
  vireo_donor_ID = vector("list", length(seurat_list))
  names(vireo_donor_ID) = names(seurat_list)
  for (pool in pools) {
    message("Working on pool ", pool,"...")
    files = list.dirs(
      paste0(
        directory,"cellranger/",
        pool
      ),
      recursive = FALSE,
      full.names = FALSE
    )
    for (file in files) {
      message("... file ", file,"...")
      
      vireo_donor_ID[[paste0(pool,".",file)]]  = read.table(paste0("../../../OTAR2065_scRNA_data_processing/data/cellranger/",
                                                                   pool,"/",file,"/vireoOutput/donor_ids.tsv"),
                                                            header = TRUE)
      vireo_donor_ID[[paste0(pool,".",file)]] = vireo_donor_ID[[paste0(pool,".",file)]] %>% dplyr::select_if(!names(.) %in% c('doublet_logLikRatio')) # drop column if it exists
      vireo_donor_ID[[paste0(pool,".",file)]]$sample =  unique(seurat_list[[paste0(pool,".",file)]]$sample)
      
    }
  }
  vireo_donor_ID = do.call("rbind",vireo_donor_ID)
  vireo_donor_ID$cell = gsub(pattern = "-1",replacement = "",vireo_donor_ID$cell)
  vireo_donor_ID$cell = paste(vireo_donor_ID$sample,vireo_donor_ID$cell,sep = "_")
  for(n in names(seurat_list)){
    seurat_list[[n]][["orig_cell"]] = rownames(seurat_list[[n]]@meta.data)
    seurat_list[[n]][["cell"]] = paste(seurat_list[[n]]$sample, rownames(seurat_list[[n]]@meta.data),sep = "_")
    seurat_list[[n]]@meta.data = merge(seurat_list[[n]]@meta.data,vireo_donor_ID)
    rownames(seurat_list[[n]]@meta.data) = seurat_list[[n]]$orig_cell
  }
  rm(vireo_donor_ID)
  return(seurat_list)
}

filter_x_UMIcount_bottompercent_per_library = function(seurat_object_list, bottom_percent){
  message("Removing bottom ", bottom_percent, 
          " percent of cells by UMI counts per library.")
  
  
  for (sample in names(seurat_object_list)){
    message("Working on sample ", sample)
    bottom_threshold = quantile(seurat_object_list[[sample]]$nCount_RNA, bottom_percent / 100) 
    message("Bottom count threshold ", bottom_threshold)
    
    
    seurat_object_list[[sample]] = subset(seurat_object_list[[sample]], subset = nCount_RNA > bottom_threshold)
    
  }
  
  return(seurat_object_list)
  
}

UMAPS_all_libs = function(seurat_object,raster.status = TRUE,
                          pt.size.status=3) {
  seurat_object$log10_nCount_RNA = log10( seurat_object$nCount_RNA)
  p1 = SCpubr::do_DimPlot(seurat_object, label = F,
                          reduction="umap.full",
                          group.by = "cluster_full",
                          raster=raster.status, pt.size = pt.size.status) + NoLegend() + ggtitle("Clusters")
  p2 = SCpubr::do_DimPlot(seurat_object,reduction="umap.full", raster=raster.status, pt.size = pt.size.status,
                          label = F, group.by = "orig.ident")  + 
    ggtitle("Libraries") + NoLegend()
  
  p3 = SCpubr::do_FeaturePlot(seurat_object,  features  = "percent.mt",
                              reduction="umap.full",raster=raster.status, pt.size = pt.size.status) + 
    ggtitle("Mitochondrial percentage") 
  
  p4 = SCpubr::do_FeaturePlot(seurat_object,  features  = "log10_nCount_RNA",
                              reduction="umap.full",raster=raster.status, pt.size = pt.size.status) + 
    ggtitle("log10(UMI counts)")
  
  p5 = SCpubr::do_DimPlot(seurat_object,reduction="umap.full",raster=raster.status, pt.size = pt.size.status,
                          label = F,  group.by = "sequencing_batch")  +
    ggtitle("Seq. batch")
  
  p6 = SCpubr::do_DimPlot(seurat_object,   reduction="umap.full",
                          raster=raster.status, pt.size = pt.size.status,
                          group.by = "Phase", shuffle = TRUE)  + # shuffle is TRUE by default
    ggtitle("Cell cycle")
  
 
  p7 = SCpubr::do_DimPlot(seurat_object,reduction="umap.full",
                          raster=raster.status, pt.size = pt.size.status,label = F,  
                          group.by = "pool"
  )  + 
    ggtitle("Pool")
  
 
  p = (p1 | p2) / (p3 | p4) / (p5 | p6) / (p7 | patchwork::plot_spacer())
  
  return(p)
}


# load scRNA-seq donor summaries######


load_sc_donors = function(doublets_unassigned_in_proportion = TRUE,
                          pools = paste0("pool",2:17),
                          directory =  "../../../OTAR2065_scRNA_data_processing/data/"){
  lines_in_pools = read.table(paste0(directory,"lines_in_pools.txt"), 
                              header = TRUE)
  #lines_in_pools has three columns: Pool, N_lines, Lines (last one has lines separated by ;)
  sc_ncells = list()
  for (pool in pools) {
    message("Working on pool ", pool,"...")
    files = list.dirs(
      paste0(
        directory,"cellranger/",
        pool
      ),
      recursive = FALSE,
      full.names = FALSE
    )
    for (file in files) {
      message("... file ", file,"...")
      
      sc_ncells[[paste0(pool,".",file)]]  = read.table(paste0("../../../OTAR2065_scRNA_data_processing/data/cellranger/",
                                                              pool,"/",file,"/vireoOutput/summary.tsv"),
                                                       header = TRUE)
      colnames( sc_ncells[[paste0(pool,".",file)]]) = c("Line","Frequency")
      
      
      # fill in missing donors
      n_lines = lines_in_pools[lines_in_pools$Pool == pool,"N_lines"]
      lines = unlist(strsplit(lines_in_pools[lines_in_pools$Pool == pool,"Lines"],split = ";"))
      sc_ncells[[paste0(pool, ".", file)]] =  sc_ncells[[paste0(pool, ".", file)]] %>%
        rows_insert(tibble(Line = lines[!lines %in% sc_ncells[[paste0(pool, ".", file)]]$Line])) %>%
        replace_na(list(Frequency = 0))
      
      sc_ncells[[paste0(pool,".",file)]]  = sc_ncells[[paste0(pool,".",file)]] %>%
        mutate(pool = pool, sample = file, 
               n_lines = n_lines) %>%
        mutate(starting_prop = 1/n_lines)
      if(doublets_unassigned_in_proportion){  # consider doublet and unassigned in final proportion 
        sc_ncells[[paste0(pool,".",file)]] = sc_ncells[[paste0(pool,".",file)]]%>%
          mutate(final_prop = Frequency / sum(Frequency)) 
        # calculate frequency and final_prop for the sum of doublet and unassigned
        unassigned_plus_doublets = sum(sc_ncells[[paste0(pool,".",file)]][sc_ncells[[paste0(pool,".",file)]]$Line=="unassigned","Frequency"],
                                       sc_ncells[[paste0(pool,".",file)]][sc_ncells[[paste0(pool,".",file)]]$Line=="doublet","Frequency"])
        unassigned_plus_doublets_pcnt = sum(sc_ncells[[paste0(pool,".",file)]][sc_ncells[[paste0(pool,".",file)]]$Line=="unassigned","final_prop"],
                                            sc_ncells[[paste0(pool,".",file)]][sc_ncells[[paste0(pool,".",file)]]$Line=="doublet","final_prop"])
        
        sc_ncells[[paste0(pool,".",file)]] = sc_ncells[[paste0(pool,".",file)]]%>%
          add_row(Line = "unassigned_plus_doublets",
                  Frequency = unassigned_plus_doublets, 
                  final_prop = unassigned_plus_doublets_pcnt,
                  pool = pool,
                  sample = file,
                  n_lines = n_lines
          )
        # sumar doublets and unassigned por separado y anadir como row extra
      }else{ # prefilter doublet and unassigned first 
        sc_ncells[[paste0(pool,".",file)]] = sc_ncells[[paste0(pool,".",file)]]%>%
          filter(!Line %in% c("doublet","unassigned")) %>%
          mutate(final_prop = Frequency / sum(Frequency))
      }
      
      
      
    }
  }
  
  sc_ncells = do.call(rbind,Map(cbind,set = names(sc_ncells),sc_ncells))
  return(sc_ncells)
}
