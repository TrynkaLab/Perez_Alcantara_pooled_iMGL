# helper functions


# load scRNA-seq donor summaries######


load_sc_donors = function(doublets_unassigned_in_proportion = TRUE,
                          pools = paste0("pool", 2:17) ,
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
        dplyr::mutate(pool = pool, sample = file, 
               n_lines = n_lines) %>%
        dplyr::mutate(starting_prop = 1/n_lines)
      if(doublets_unassigned_in_proportion){  # consider doublet and unassigned in final proportion 
        sc_ncells[[paste0(pool,".",file)]] = sc_ncells[[paste0(pool,".",file)]]%>%
          dplyr::mutate(final_prop = Frequency / sum(Frequency)) 
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
         dplyr::filter(!Line %in% c("doublet","unassigned")) %>%
          dplyr::mutate(final_prop = Frequency / sum(Frequency))
      }
      
      
      
    }
  }
  
  sc_ncells = do.call(rbind,Map(cbind,set = names(sc_ncells),sc_ncells))
  return(sc_ncells)
}

load_sc_donors_demuxlet = function(doublets_unassigned_in_proportion = TRUE,
                          pools = c("pool2", "pool3", "pool4", "pool5", "pool6","pool7",
                                    "pool8","pool9","pool10") ,
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
                                                              pool,"/",file,"/demuxlet/demuxlet.best"),
                                                       header = TRUE)
      
      singlets = sc_ncells[[paste0(pool,".",file)]] %>%
       dplyr::filter(DROPLET.TYPE == "SNG") %>%
        group_by(SNG.BEST.GUESS) %>%
        summarise(Frequency = n()) %>%
        rename(Line = SNG.BEST.GUESS)
        
      unassigned = sc_ncells[[paste0(pool,".",file)]] %>%
       dplyr::filter(DROPLET.TYPE == "AMB") %>%
        group_by(DROPLET.TYPE) %>%
        summarise(Frequency = n())  %>%
        rename(Line = DROPLET.TYPE) %>%
        dplyr::mutate(Line = str_replace(Line, "AMB", "unassigned"))

      doublets = sc_ncells[[paste0(pool,".",file)]] %>%
       dplyr::filter(DROPLET.TYPE == "DBL") %>%
        group_by(DROPLET.TYPE) %>%
        summarise(Frequency = n())  %>%
        rename(Line = DROPLET.TYPE) %>%
        dplyr::mutate(Line = str_replace(Line, "DBL", "doublet"))
      
      sc_ncells[[paste0(pool,".",file)]] = rbind(singlets,unassigned,doublets)
      
      
      
      # fill in missing donors
      n_lines = lines_in_pools[lines_in_pools$Pool == pool,"N_lines"]
      lines = unlist(strsplit(lines_in_pools[lines_in_pools$Pool == pool,"Lines"],split = ";"))
      sc_ncells[[paste0(pool, ".", file)]] =  sc_ncells[[paste0(pool, ".", file)]] %>%
        rows_insert(tibble(Line = lines[!lines %in% sc_ncells[[paste0(pool, ".", file)]]$Line])) %>%
        replace_na(list(Frequency = 0))
      
      sc_ncells[[paste0(pool,".",file)]]  = sc_ncells[[paste0(pool,".",file)]] %>%
        dplyr::mutate(pool = pool, sample = file, 
               n_lines = n_lines) %>%
        dplyr::mutate(starting_prop = 1/n_lines)
      if(doublets_unassigned_in_proportion){  # consider doublet and unassigned in final proportion 
        sc_ncells[[paste0(pool,".",file)]] = sc_ncells[[paste0(pool,".",file)]]%>%
          dplyr::mutate(final_prop = Frequency / sum(Frequency)) 
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
         dplyr::filter(!Line %in% c("doublet","unassigned")) %>%
          dplyr::mutate(final_prop = Frequency / sum(Frequency))
      }
      
      
      
    }
  }
  
  sc_ncells = do.call(rbind,Map(cbind,set = names(sc_ncells),sc_ncells))
  return(sc_ncells)
}

# plot histograms of scaled proportions######

plot_prop_hist = function(df, 
                          directory = "../../data/results/2.efficiency/",
                          plot.name = "scaled_prop_histograms.png",
                          nbins = 50,
                          plot.width = 11,
                          plot.height = 7){
  
  p1 = df %>% ggplot( aes(x=final_prop)) + 
    geom_histogram(color="black", fill="white", bins = nbins) + theme_minimal() + ggtitle("Final proportion")
  p2 = df %>% ggplot( aes(x=starting_prop)) + 
    geom_histogram(color="black", fill="white", bins = nbins) + theme_minimal() + ggtitle("Starting proportion")
  p3 = df %>% ggplot( aes(x=centered_prop)) + 
    geom_histogram(color="black", fill="white", bins = nbins) + theme_minimal() + 
    ggtitle("Final_prop - starting_prop")
  p4 = df %>% ggplot( aes(x=log1p_scaled_prop)) + 
    geom_histogram(color="black", fill="white", bins = nbins) + theme_minimal() + 
    ggtitle("log1p(final_prop) / log1p(starting_prop)")
  p5 = df %>% ggplot( aes(x=scaled_centered_prop)) + 
    geom_histogram(color="black", fill="white", bins = nbins) + theme_minimal() + 
    ggtitle("(final_prop - starting_prop) / sd(final_prop)")
  p6 =  df %>% ggplot( aes(x=minmax_centered_prop)) + 
    geom_histogram(color="black", fill="white", bins = nbins) + theme_minimal() + 
    ggtitle("(abs(centered_prop) - min(abs(centered_prop))) / (max(abs(centered_prop)) - min(abs(centered_prop)) ")
  
  p =(p1 + p2) / (p3 + p4) / (p5 + p6)
  png(paste0(directory,plot.name),
      width = plot.width, height = plot.height, units = "in", res = 400)
  plot(p)
  dev.off()
}

# plot scatterplots of scaled proportions######

plot_prop_scatter = function(df, 
                             directory = "../../data/results/2.efficiency/",
                             plot.name = "scaled_prop_scatter.png",
                             plot.width = 11,
                             plot.height = 6){
  
  p1 = df %>%
    ggplot(aes(x = centered_prop, y = scaled_centered_prop, 
               col = ifelse(final_prop == 0, "grey","black"))) + 
    geom_point() +
    scale_color_identity() + 
    theme_bw() +
    ggtitle( "grey: final prop. = 0")
  
  p2 = df %>%
    ggplot(aes(x = centered_prop, y = scaled_centered_prop), col = pool) + 
    geom_point() +
    theme_bw() +
    ggtitle("color by pool")
  
  p3 = df %>%
    ggplot(aes(x = scaled_centered_prop, y = log1p_scaled_prop, 
               col = ifelse(final_prop < starting_prop, "darkred","darkblue"))) + 
    geom_point() +
    scale_color_identity() + 
    theme_bw() +
    ggtitle("red: final prop. < starting prop.")
  
  p4 = df %>%
    ggplot(aes(x = scaled_centered_prop, y = log1p_scaled_prop),col = pool) + 
    geom_point() +
    theme_bw() +
    ggtitle("color by pool")
  
  p5 = df %>%
    ggplot(aes(x = minmax_centered_prop, y = log1p_scaled_prop, 
               col = ifelse(final_prop < starting_prop, "darkred","darkblue"))) + 
    geom_point() +
    scale_color_identity() + 
    theme_bw() +
    ggtitle("red: final prop. < starting prop.")
  
  p6 = df %>%
    ggplot(aes(x = minmax_centered_prop, y = log1p_scaled_prop, 
               col = pool)) + 
    geom_point() +
    theme_bw() +
    ggtitle("color by pool")
  
  
  p =(p1 | p2) / (p3 | p4) / (p5 | p6) 
  png(paste0(directory,plot.name),
      width = plot.width, height = plot.height, units = "in", res = 400)
  plot(p)
  dev.off()
}

# plot density plots of scaled proportions separated by pool #####

plot_prop_density = function(df, 
                             directory = "../../data/results/2.efficiency/",
                             plot.name = "scaled_prop_density_pool.png",
                             col.by = "set",
                             plot.legend = FALSE,
                             plot.width = 11,
                             plot.height = 7){
  p1 = df %>% 
    ggplot(aes(color=!!sym(col.by), x=scaled_centered_prop, fill= !!sym(col.by))) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    ylab("Cell density") +
    xlab("(final_prop - starting_prop) / sd(final_prop)") +
    theme(legend.position=ifelse(plot.legend,"right","none"))
  
  
  p2 = df %>% 
    ggplot(aes(color=!!sym(col.by), x=log1p_scaled_prop, fill= !!sym(col.by))) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    ylab("Cell density") +
    xlab("(final_prop - mean(2 shared donors) ) / sd(final_prop)") +
    theme(legend.position=ifelse(plot.legend,"right","none"))
  
  
  p3 = df %>% 
    ggplot(aes(color=!!sym(col.by), x=minmax_centered_prop, fill= !!sym(col.by))) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    ylab("Cell density") +
    xlab("log1p(final_prop) / log1p(starting_prop) (Puigdevall's prolif rate)") + 
    theme(legend.position=ifelse(plot.legend,"right","none"))
  p4 = df %>% 
    ggplot(aes(color=!!sym(col.by), x=log1p_substract_starting_scaled, fill= !!sym(col.by))) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    ylab("Cell density") +
    xlab("(log1p(final_prop) - log1p(starting_prop)) / log1p(mean(2 shared donors))") +
    theme(legend.position=ifelse(plot.legend,"right","none"))
  
  
  p =(p1 | p2) / (p3 | p4)  +   plot_layout(guides = "collect")
  png(paste0(directory,plot.name),
      width = plot.width, height = plot.height, units = "in", res = 400)
  plot(p)
  dev.off()
}

filter_seurat = function(seurat_object){
  
  message("Performing filter by number of genes, UMI counts, mitochondrial percentage and Vireo doublets.")
  seurat_object = subset(seurat_object, subset = nFeature_RNA > 1000  & percent.mt < 10)
  seurat_object = subset(seurat_object, subset = donor_id != "doublet" ) 
  seurat_object = subset(seurat_object, subset = donor_id != "unassigned") 
  message("Now the object has ", dim(seurat_object)[1], " genes and ", dim(seurat_object)[2], " cells.")
  return(seurat_object)
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
UMAPS_all_libs = function(seurat_object) {
  p1 = SCpubr::do_DimPlot(seurat_object, label = F) + NoLegend() + ggtitle("Clusters")
  p2 = SCpubr::do_DimPlot(seurat_object, label = F, group.by = "orig.ident")  + 
    ggtitle("Libraries") + NoLegend()
  
  p3 = SCpubr::do_FeaturePlot(seurat_object,  features  = "percent.mt") + 
    ggtitle("Mitochondrial percentage") 
  
  p4 = SCpubr::do_FeaturePlot(seurat_object,  features  = "nCount_RNA") + 
    ggtitle("UMI counts")
  
  p5 = SCpubr::do_DimPlot(seurat_object,
                          label = F,  group.by = "sequencing_batch")  +
    ggtitle("Seq. batch")
  
  p6 = SCpubr::do_DimPlot(seurat_object,   
                          group.by = "Phase", shuffle = TRUE)  + # shuffle is TRUE by default
    ggtitle("Cell cycle")
  
  # p7 = SCpubr::do_DimPlot(seurat_object, label = F,  
  #                         shuffle = TRUE,
  #                         group.by = "pool", 
  #                         order = c( "pool2",  "pool3","pool4","pool5","pool6"), # fix colors
  #                         cols = c("#0d3b66","#ee964b","#588b8b","#E34A6F","#9EE493"))  + 
  #   ggtitle("Pool")
  p7 = SCpubr::do_DimPlot(seurat_object, label = F,  
                          group.by = "pool"
  )  + 
    ggtitle("Pool")
  
  p8 = SCpubr::do_DimPlot(seurat_object,  
                          label = F,  group.by = "treatment")  + 
    ggtitle("Treatment")
  
  p = (p1 | p2) / (p3 | p4) / (p5 | p6) / (p7 | p8)
  
  return(p)
}

UMAPS_per_lib = function(seurat_object) {
  p1 = SCpubr::do_DimPlot(seurat_object, label = F) + NoLegend() + ggtitle("Clusters")
  
  p2 = SCpubr::do_FeaturePlot(seurat_object,  features  = "percent.mt") + 
    ggtitle("Mitochondrial percentage") 
  
  p3 = SCpubr::do_FeaturePlot(seurat_object,  features  = "nCount_RNA") + 
    ggtitle("UMI counts")
  
  
  p4 = SCpubr::do_DimPlot(seurat_object,   
                          group.by = "Phase", shuffle = TRUE)  + 
    ggtitle("Cell cycle")
  
  p5 = SCpubr::do_DimPlot(seurat_object,
                          label = F,  group.by = "predicted.id")  +
    ggtitle("Seurat predictions")
  
  p6 = SCpubr::do_FeaturePlot(seurat_object,
                              features = "prediction.score.Microglia"
  ) + 
    ggtitle("Microglia score")
  
  p = (p1 | p2) / (p3 | p4) / (p5 | p6)
  
  return(p)
}

calculate_PCA_UMAP_neighbors_clusters_merged = function(seurat_object){
  seurat_object = RunPCA(seurat_object, verbose = FALSE)
  
  seurat_object = FindNeighbors(seurat_object, dims = 1:30, verbose = FALSE)
  seurat_object = FindClusters(seurat_object, 
                               resolution = 0.3, # Selected in 1.1.Determine_clustering_params - choosing lower loses CRM and pvM info
                               verbose = FALSE)
  
  seurat_object = RunUMAP(seurat_object, dims = 1:30, verbose = FALSE)
  
  
  
  return(seurat_object)
}

calculate_PCA_UMAP_neighbors_clusters_integrated = function(seurat_object){
  DefaultAssay(seurat_object) = "integrated"
  seurat_object = RunPCA(seurat_object, verbose = FALSE)
  
  seurat_object = FindNeighbors(seurat_object, dims = 1:30, verbose = FALSE)
  seurat_object = FindClusters(seurat_object, 
                               resolution = 0.3, # Selected in 1.1.Determine_clustering_params - choosing lower loses CRM and pvM info
                               verbose = FALSE)
  
  seurat_object = RunUMAP(seurat_object, dims = 1:30, verbose = FALSE)
  
  
  
  return(seurat_object)
}

calculate_PCA_UMAP_neighbors_clusters = function(seurat_object){
  seurat_object = RunPCA(seurat_object, verbose = FALSE)
  
  seurat_object = FindNeighbors(seurat_object, dims = 1:30, verbose = FALSE)
  seurat_object = FindClusters(seurat_object,
                               verbose = FALSE)
  
  seurat_object = RunUMAP(seurat_object, dims = 1:30, verbose = FALSE)
  
  
  
  return(seurat_object)
}

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
    if(unique(seurat_list[[n]][["pool"]]) ==  "pool11" ){
      seurat_list[[n]][["sequencing_batch"]] = "batch8"
      
    }
    
    if(unique(seurat_list[[n]][["pool"]]) ==  "pool13" ){
      seurat_list[[n]][["sequencing_batch"]] = "batch9"
      
    }
    if(unique(seurat_list[[n]][["pool"]]) ==  "pool12" ){
      seurat_list[[n]][["sequencing_batch"]] = "batch10"
      
    }
    if(unique(seurat_list[[n]][["pool"]]) ==  "pool14" ){
      seurat_list[[n]][["sequencing_batch"]] = "batch11"
      
    }
    if(unique(seurat_list[[n]][["pool"]]) ==  "pool15" ){
      seurat_list[[n]][["sequencing_batch"]] = "batch12"
      
    }
  }
  return(seurat_list)
}

incorporate_vireo_info = function(seurat_list,
                                  pools = paste0("pool",2:15),
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

plot_QC_histograms = function(seurat_merged_object,
                              plot.legend = TRUE){
  p1 = seurat_merged_object@meta.data %>% 
    ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    geom_vline(xintercept = 40000, color = "red") +
    ylab("Cell density") +
    xlab("UMI counts") + 
    NoLegend()
  p2 = seurat_merged_object@meta.data %>% 
    ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10() + 
    geom_vline(xintercept = 1000) +
    ylab("Cell density") +
    xlab("Number of genes") + 
    NoLegend()
  
  if(plot.legend == FALSE){
    p3 = seurat_merged_object@meta.data %>% 
      ggplot(aes(color=orig.ident, x=percent.mt, fill=orig.ident)) + 
      geom_density(alpha = 0.2) + 
      theme_classic() +
      geom_vline(xintercept = 10) +
      ylab("Cell density") +
      xlab("% of mitochondrial genes") + 
      NoLegend()
  }else{
    p3 = seurat_merged_object@meta.data %>% 
      ggplot(aes(color=orig.ident, x=percent.mt, fill=orig.ident)) + 
      geom_density(alpha = 0.2) + 
      theme_classic() +
      geom_vline(xintercept = 10) +
      ylab("Cell density") +
      xlab("% of mitochondrial genes")
  }
  
  
  return(p1 + p2 + p3)
}

plot_QC_vireo_dotplots = function(seurat_merged_object){
  # the merged object needs to have vireo information
  # Remove doublets first
  seurat_merged_object =  subset(seurat_merged_object, 
                                 subset = donor_id != c("doublet")) 
  
  p1 = seurat_merged_object@meta.data %>%
    ggplot(aes(x=nCount_RNA, y = prob_max, size = n_vars, 
               col = donor_id == "unassigned")) + 
    geom_point() +
    scale_colour_manual(name = 'Vireo ID = unassigned', 
                        values = setNames(c('red','darkblue'),c(TRUE, FALSE))) +
    
    theme_classic() +
    ylab("Vireo singlet assignment probability") +
    xlab("UMI counts") + 
    NoLegend()
  
  p2 = seurat_merged_object@meta.data  %>%
    ggplot(aes(x=percent.mt, y = prob_max, size = n_vars, col = donor_id == "unassigned")) + 
    geom_point() +
    scale_colour_manual(name = 'Vireo ID = unassigned', 
                        values = setNames(c('red','darkblue'),c(TRUE, FALSE))) +
    theme_classic() +
    ylab("Vireo singlet assignment probability") +
    xlab("Mitochondrial percentage")
  
  p3 = seurat_merged_object@meta.data  %>%
    ggplot(aes(x=n_vars, y = prob_max,col = donor_id == "unassigned")) + 
    geom_point() +
    scale_x_log10() +
    scale_colour_manual(name = 'Vireo ID = unassigned', 
                        values = setNames(c('red','darkblue'),c(TRUE, FALSE))) +
    theme_classic() +
    ylab("Vireo singlet assignment probability") +
    xlab("Number of variants used for Vireo")
  
  return(p1 + p2 + p3)
}

# get values from linear model regresion

#### fix below to include the corr test estimate
# cor.test(df$x,df$y)$estimate
lm_eqn <- function(df){
  colnames(df) = c("x","y")
  #m = lm(y ~ x, df) # I think I should just use the Pearson R because there is no dependent / independent variable here
  p =  cor.test(df$x,df$y)$p.value
  r =  as.numeric(cor.test(df$x,df$y)$estimate)
  eq <- substitute(italic(R) == r*","~~italic("p")~"="~p, 
                   list(r =  format(r, digits = 3),
                        p = format(p, digits = 2,  scientific = TRUE)))
  as.character(as.expression(eq));
  
}

cor_test_estimate <- function(df){
  colnames(df) = c("x","y")
  r =  as.numeric(cor.test(df$x,df$y)$estimate)
  return(r)
}

#### correlation plots for proliferation rate
within_pool_within_cond_between_rep = function(df, pool, condition) {
  p = pool
  t = condition
  subset = df %>%
    dplyr::mutate(treatment = case_when(
      grepl("plst|untreated|Untreated|AG|YC|ITMG", sample) ~ "untreated",
      grepl("IFN", sample) ~ "IFN",
      grepl("LPS", sample) ~ "LPS"
    )) %>%
    
    dplyr::filter(pool == p, treatment == t) %>%
    dplyr::select(Line, sample, scaled_centered_prop) %>%
    pivot_wider(names_from = sample, values_from = scaled_centered_prop)
  
  pairwise_comb = combn(colnames(subset[, -1]), 2)
  p = vector(mode = "list", length = ncol(pairwise_comb))
  
  for (n in 1:ncol(pairwise_comb)) {
    subset2 = subset %>%
      dplyr::select(c(pairwise_comb[, n]))
    
    p[[n]] = subset2 %>%
      ggplot(aes_string(x = pairwise_comb[1, n], y = pairwise_comb[2, n])) +
      geom_point() +
      theme_bw() +
      geom_smooth(
        method = lm,
        col = "black",
        alpha = 0.6,
        linetype = "dashed"
      ) +
      geom_text(
        x = 1,
        y = max(subset[, 2]) - 1,
        label = lm_eqn(subset2),
        parse = TRUE
      ) +
      xlab(paste0(
        "Line proliferation: ",
        pairwise_comb[1, n]
      )) +
      ylab(paste0(
        "Line proliferation: ",
        pairwise_comb[2, n]
      ))
    
  }
  p = gridExtra::grid.arrange(grobs = p,
                              ncol = ceiling(sqrt(ncol(pairwise_comb))))
  return(p)
}

between_pools_within_conditions = function(df,condition){
  subset = df %>%
    dplyr::mutate(treatment = case_when(
      grepl("plst|untreated|Untreated|AG|YC|ITMG", sample) ~ "untreated",
      grepl("IFN", sample) ~ "IFN",
      grepl("LPS", sample) ~ "LPS"
    )) %>%
    dplyr::filter(treatment == condition) %>%
    # average between replicates of pool 2, 9 and 10
    dplyr::group_by(Line,pool) %>%
    dplyr::summarise(scaled_centered_prop_median = median(scaled_centered_prop)) %>% # can't work with means because I have negative numbers
    dplyr::select(Line, pool, scaled_centered_prop_median) %>%
    pivot_wider(names_from = pool, values_from = scaled_centered_prop_median) %>%
    dplyr::ungroup() %>%
    dplyr::filter(rowSums(is.na(.))<(ncol(.)-2)) %>% # remove rows with number of NAs equal or higher than number of pools - 1
    # this means the line is only present in one pool, or the line name is missing for some reason
    tidyr::pivot_longer(cols = !Line,names_to = "pool", values_to = "scaled_centered_prop_median") %>%
    na.omit() 
  
  final = list()
  for(line in unique(subset$Line)){
    subset2 = subset %>%
      dplyr::filter(Line == line) %>%
      ungroup()
    if(nrow(subset2) == 2){
      final[[line]] = subset2
    }else{
      pairwise_comb = combn(subset2$pool, 2) # calculate pairwise combinations
      # check how many there are and rename lines with new number
      more_than_2 = list()
      for(n in 1:ncol(pairwise_comb)){
        row1 = subset2 %>%
          dplyr::filter(pool %in% pairwise_comb[1,n]) %>%
          add_row(Line = paste0(.$Line,n), pool = pairwise_comb[1,n],
                  scaled_centered_prop_median = unlist(subset2[subset2$pool == pairwise_comb[1,n],"scaled_centered_prop_median"])) %>%
          dplyr::slice_tail(n = 1)
        
        row2 = subset2 %>%
          dplyr::filter(pool %in% pairwise_comb[2,n]) %>%
          add_row(Line = paste0(.$Line,n), pool = pairwise_comb[2,n],
                  scaled_centered_prop_median = unlist(subset2[subset2$pool == pairwise_comb[2,n],"scaled_centered_prop_median"])) %>%
          dplyr::slice_tail(n = 1)
        
        more_than_2[[n]] = rbind(row1,row2)
      }
      final[[line]] = do.call("rbind",more_than_2)
      
    }
    
  }
  final = do.call("rbind",final)
  
  final2 = final %>%
    # assign pools A or B at random within lines
    group_by(Line) %>% 
    dplyr::mutate(pool=sample(pool)) %>% # random shuffle within line (sample without replacement, so shouldn't repeat pools)
    dplyr::arrange(Line,pool) %>% # re-arrange
    dplyr::mutate(pool_renamed = ifelse(row_number()%%2, "pool A" ,"pool B")) %>%
    dplyr::select(Line,pool_renamed,scaled_centered_prop_median) %>% # to check lines
    pivot_wider(names_from = pool_renamed, values_from = scaled_centered_prop_median) %>%
    ungroup() %>%
    dplyr::select(`pool A`,`pool B`)
  
  
  p = final2 %>%
    ggplot(aes(x = `pool A`, y =`pool B`)) +
    geom_point() +
    theme_bw() +
    geom_smooth(
      method = lm,
      col = "black",
      alpha = 0.6,
      linetype = "dashed"
    ) +
    geom_text(
      x = 1,
      y = max(final2[, 2]) - 1,
      label = lm_eqn(final2),
      parse = TRUE
    ) +
    xlab("Line proliferation: pool A") +
    ylab("Line proliferation: pool B") + 
    ggtitle(paste0("Between pools, within conditions: ",condition))
  
  return(p)
}

between_pools_within_conditions_correlations = function(df,condition){
  subset = df %>%
    dplyr::mutate(treatment = case_when(
      grepl("plst|untreated|Untreated|AG|YC|ITMG", sample) ~ "untreated",
      grepl("IFN", sample) ~ "IFN",
      grepl("LPS", sample) ~ "LPS"
    )) %>%
    dplyr::filter(treatment == condition) %>%
    # average between replicates of pool 2, 9 and 10
    dplyr::group_by(Line,pool) %>%
    dplyr::summarise(scaled_centered_prop_median = median(scaled_centered_prop)) %>% # can't work with means because I have negative numbers
    dplyr::select(Line, pool, scaled_centered_prop_median) %>%
    pivot_wider(names_from = pool, values_from = scaled_centered_prop_median) %>%
    dplyr::ungroup() %>%
    dplyr::filter(rowSums(is.na(.))<(ncol(.)-2)) %>% # remove rows with number of NAs equal or higher than number of pools - 1
    # this means the line is only present in one pool, or the line name is missing for some reason
    tidyr::pivot_longer(cols = !Line,names_to = "pool", values_to = "scaled_centered_prop_median") %>%
    na.omit() 
  
  final = list()
  for(line in unique(subset$Line)){
    subset2 = subset %>%
      dplyr::filter(Line == line) %>%
      ungroup()
    if(nrow(subset2) == 2){
      final[[line]] = subset2
    }else{
      pairwise_comb = combn(subset2$pool, 2) # calculate pairwise combinations
      # check how many there are and rename lines with new number
      more_than_2 = list()
      for(n in 1:ncol(pairwise_comb)){
        row1 = subset2 %>%
          dplyr::filter(pool %in% pairwise_comb[1,n]) %>%
          add_row(Line = paste0(.$Line,n), pool = pairwise_comb[1,n],
                  scaled_centered_prop_median = unlist(subset2[subset2$pool == pairwise_comb[1,n],"scaled_centered_prop_median"])) %>%
          dplyr::slice_tail(n = 1)
        
        row2 = subset2 %>%
          dplyr::filter(pool %in% pairwise_comb[2,n]) %>%
          add_row(Line = paste0(.$Line,n), pool = pairwise_comb[2,n],
                  scaled_centered_prop_median = unlist(subset2[subset2$pool == pairwise_comb[2,n],"scaled_centered_prop_median"])) %>%
          dplyr::slice_tail(n = 1)
        
        more_than_2[[n]] = rbind(row1,row2)
      }
      final[[line]] = do.call("rbind",more_than_2)
      
    }
    
  }
  final = do.call("rbind",final)
  
  cor_sample = list()
  for(shuffle in 1:1000){
    
    final2 = final %>%
      # assign pools A or B at random within lines
      group_by(Line) %>% 
      dplyr::mutate(pool=sample(pool)) %>% # random shuffle within line (sample without replacement, so shouldn't repeat pools)
      dplyr::arrange(Line,pool) %>% # re-arrange
      dplyr::mutate(pool_renamed = ifelse(row_number()%%2, "pool A" ,"pool B")) %>%
      dplyr::select(Line,pool_renamed,scaled_centered_prop_median) %>% # to check lines
      pivot_wider(names_from = pool_renamed, values_from = scaled_centered_prop_median) %>%
      ungroup() %>%
      dplyr::select(`pool A`,`pool B`)
    
    cor_sample[[shuffle]] = cor_test_estimate(final2)
  }
  cor_sample = unlist(cor_sample)
  return(cor_sample)
  
}

## gene symbol to Emsemblr ID ############
symbol_to_ensembl <- function(x){
  ensembl = mapIds(org.Hs.eg.db, x, "ENSEMBL","SYMBOL")
  return(ensembl)
}

# Ensembl ID to symbol
ensembl_to_symbol <- function(x){
  symbol = mapIds(org.Hs.eg.db, x, "SYMBOL","ENSEMBL")
  return(symbol)
}

symbol_to_entrez <- function(x){
  symbol = mapIds(org.Hs.eg.db, x, "ENTREZID","SYMBOL")
  return(symbol)
}

plot_pca = function(x,d ){
  pca1 = prcomp(t(x), retx = TRUE)
  
  plot(pca1, type = "l") #variance vs first 10 components
  summary(pca1)          #importance of each component (important line is "proportion of variance")
  
  percentVar = (pca1$sdev)^2 / sum(pca1$sdev^2)
  percentVar = round(100 * percentVar)
  pcs = as.data.frame(pca1$x)
  pcs = cbind(pcs,d)
  
  
  p_treatment = ggplot(pcs, aes(PC1, PC2, colour = treatment)) + 
    geom_point(size = 4) + xlab(paste0("PC1: " ,percentVar[ 1 ],"% variance")) + 
    ylab(paste0( "PC2: ",percentVar[ 2 ],"% variance" )) 
  p_treatment = p_treatment + theme(	panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                     panel.background = element_blank(),
                                     panel.border = element_rect(fill = NA, colour = "black"), 
                                     legend.key = element_blank(),# legend.position = c(0.5,0.5),
                                     axis.title.y = element_text(face = "bold", angle = 90, size = 16, vjust = 0.2),
                                     axis.title.x = element_text(face = "bold", size = 16, vjust = 0),
                                     axis.text.x = element_text(face = "bold", colour = "black", angle = 90, size = 16, vjust = 0.2, hjust = 1 ),
                                     axis.text.y = element_text(face = "bold", colour = "black",size = 16),
                                     axis.ticks = element_line(colour = "black"),
                                     axis.line = element_line(colour = "black"),
                                     legend.text = element_text(size = 14,face = "bold"),
                                     legend.title = element_text(size = 16,face = "bold"))
  
  
  p_ncells = ggplot(pcs, aes(PC1, PC2, colour = log10_ncells, shape = treatment)) + 
    geom_point(size = 4) + xlab(paste0("PC1: " ,percentVar[ 1 ],"% variance")) + 
    ylab(paste0( "PC2: ",percentVar[ 2 ],"% variance" )) 
  p_ncells = p_ncells + theme(	panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                               panel.background = element_blank(),
                               panel.border = element_rect(fill = NA, colour = "black"), 
                               legend.key = element_blank(),# legend.position = c(0.5,0.5),
                               axis.title.y = element_text(face = "bold", angle = 90, size = 16, vjust = 0.2),
                               axis.title.x = element_text(face = "bold", size = 16, vjust = 0),
                               axis.text.x = element_text(face = "bold", colour = "black", angle = 90, size = 16, vjust = 0.2, hjust = 1 ),
                               axis.text.y = element_text(face = "bold", colour = "black",size = 16),
                               axis.ticks = element_line(colour = "black"),
                               axis.line = element_line(colour = "black"),
                               legend.text = element_text(size = 14,face = "bold"),
                               legend.title = element_text(size = 16,face = "bold"))
  
  p_treatment_PC23 = ggplot(pcs, aes(x=PC2, y=PC3, colour = log10_ncells, shape = treatment)) + 
    geom_point(size = 4) + xlab(paste0("PC2: " ,percentVar[ 2 ],"% variance")) + 
    ylab(paste0( "PC3: ",percentVar[ 3 ],"% variance" )) 
  p_treatment_PC23 = p_treatment_PC23 + theme(	panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                               panel.background = element_blank(),
                                               panel.border = element_rect(fill = NA, colour = "black"), 
                                               legend.key = element_blank(),# legend.position = c(0.5,0.5),
                                               axis.title.y = element_text(face = "bold", angle = 90, size = 16, vjust = 0.2),
                                               axis.title.x = element_text(face = "bold", size = 16, vjust = 0),
                                               axis.text.x = element_text(face = "bold", colour = "black", angle = 90, size = 16, vjust = 0.2, hjust = 1 ),
                                               axis.text.y = element_text(face = "bold", colour = "black",size = 16),
                                               axis.ticks = element_line(colour = "black"),
                                               axis.line = element_line(colour = "black"),
                                               legend.text = element_text(size = 14,face = "bold"),
                                               legend.title = element_text(size = 16,face = "bold"))
  
  p_pool_PC23 = ggplot(pcs, aes(x=PC2, y=PC3, colour = pool)) + 
    geom_point(size = 4) + xlab(paste0("PC2: " ,percentVar[ 2 ],"% variance")) + 
    ylab(paste0( "PC3: ",percentVar[ 3 ],"% variance" )) 
  p_pool_PC23 = p_pool_PC23 + theme(	panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                               panel.background = element_blank(),
                                               panel.border = element_rect(fill = NA, colour = "black"), 
                                               legend.key = element_blank(),# legend.position = c(0.5,0.5),
                                               axis.title.y = element_text(face = "bold", angle = 90, size = 16, vjust = 0.2),
                                               axis.title.x = element_text(face = "bold", size = 16, vjust = 0),
                                               axis.text.x = element_text(face = "bold", colour = "black", angle = 90, size = 16, vjust = 0.2, hjust = 1 ),
                                               axis.text.y = element_text(face = "bold", colour = "black",size = 16),
                                               axis.ticks = element_line(colour = "black"),
                                               axis.line = element_line(colour = "black"),
                                               legend.text = element_text(size = 14,face = "bold"),
                                               legend.title = element_text(size = 16,face = "bold"))
  
  return(list(p_treatment,p_ncells,p_treatment_PC23,p_pool_PC23))
}


fix_burden_matrices <- function(df) {
  df %>%
    tibble::as_tibble(rownames = "gene") %>%
    dplyr::rename_with(~str_extract(., "(?<=-)\\w+"), .cols = contains("H"))
}
extract_burden_results <- function(df) {
  res = data.frame("failed" = df[["failed"]]$coefficients[3,"Pr(>|t|)"], ### careful here
                  "successful" = df[["successful"]]$coefficients[3,"Pr(>|t|)"],
                  "gene" = df[["failed"]]$gene) %>%
 pivot_longer(cols=-gene,names_to = "proliferation_failure",values_to = "pval")
  
}

# subset column names to shared lines, including "gene"
subset_burden_shared_lines = function(x,y) { 
  x %>%
    dplyr::select(c("gene",sort(lubridate::intersect(colnames(.),y)))) 
}

filter_rare_variant_genes = function(tb,fraction) { 
  message("The fraction corresponds to ",floor((ncol(tb)-1)*fraction), " lines")
  tb %>%
    tibble::column_to_rownames(var="gene") %>%
    dplyr::filter(rowSums(. > 0) >= floor(ncol(.)*fraction)) %>%
    dplyr::filter(rowSums(. == 0) >= floor(ncol(.)*fraction)) %>%
    tibble::rownames_to_column(var="gene")
  
}


rare_variant_to_long = function(tbl_list) { 
  # returns a tibble
  for(i in names(tbl_list)){
    tbl_list[[i]] = tbl_list[[i]] %>%
      tidyr::pivot_longer(.,cols = !gene,names_to = "line",values_to = "rare_burden") %>%
      dplyr::mutate(rare_mutation_type = i)
  }
  return(do.call("rbind",tbl_list))
  
}

scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

center_this <- function(x){
  (x - mean(x, na.rm=TRUE)) 
}

fill_columns_df2_with_df1_NA = function(df1,df2){
  # used in MOFA_phenotype_integration.R
  cols_to_bind= df1 %>%
    as_tibble() %>%
    dplyr::slice_head(n=1) %>%
    dplyr::select(setdiff(colnames(df1 ), colnames(df2 ))) %>%
    mutate_all(~ NA)
  
  
  df2 = df2 %>%
    as_tibble() %>%
    dplyr::bind_cols(cols_to_bind)
  
  df2 = df2 %>%
    dplyr::relocate(colnames(df1))
  return(df2)
  
}
