# # pseudobulk 
# per treatment, pool and donor

library(Seurat)
library(BPCells)
library(SeuratObject)
#library(Signac)
library(readr)
library(future)
library(tidyverse)
library(scCustomize)
library(patchwork)
# set this option when analyzing large datasets
options(future.globals.maxSize = 60000 * 1024^2) # 60Gb

# for new input_seurat assay slots
options(seurat.object.assay.version = "v5")

output_dir="../../data/results/5.prepare_differential_expression_pseudobulks/"

################## sampling half of cells of treatmentxprolifxidxpool combinations
pseudobulk = list()
metadata = list()
for(treat in c("untreated","LPS","IFN")){
  message(treat)
  
  message("Reading in data")
  
  # seurat data
  input_seurat = readRDS(paste0("../../data/results/1.QC_v5/",treat,"_filtered_harmony/",treat,"_filtered_harmony.Rds"))
  
  gc()
  
  Layers(input_seurat[["RNA"]])
  DefaultAssay(input_seurat) = "RNA"
  Idents(input_seurat) = "cluster_full"
  # message("Joining layers")
  # #pseudobulk on donor + pool + treatment
  message("Naming proliferation categories")
  
  
  input_seurat@meta.data$proliferation_status = "Not_proliferating"
  if(treat == "untreated"){
    input_seurat@meta.data[input_seurat@meta.data$cluster_full %in%  c(3,4) , "proliferation_status"] =  "Proliferating"
    
  }
  if(treat == "LPS"){
    input_seurat@meta.data[input_seurat@meta.data$cluster_full %in%  c(4,5,7) , "proliferation_status"] =  "Proliferating"
    
  }
  if(treat == "IFN"){  # much smaller proliferating clusters on IFN treated cells
    input_seurat@meta.data[input_seurat@meta.data$cluster_full %in%  c(4) , "proliferation_status"] =  "Proliferating"
    
  }
  message("Aggregating expression: subset to avoid memory issues")
  
  pseudobulk[[treat]] = list()
  metadata[[treat]] = list()
  for(prolif in c("Not_proliferating","Proliferating")){
    subset_seurat = subset(x = input_seurat, subset = proliferation_status == prolif)
    
    subset_seurat@meta.data = subset_seurat@meta.data %>%
      group_by(donor_id, pool) %>%
      mutate(row_num = row_number()) %>%
      mutate(replicate = ifelse(row_num %% 2 == 0, "rep2", "rep1")) %>%
      select(-row_num) %>%
      ungroup() %>%
      as.data.frame()
    rownames(subset_seurat@meta.data) = paste(subset_seurat@meta.data$pool,subset_seurat@meta.data$cell,sep = ".")
    
    ### create new category at random within donor id and pool
    
    
    pseudobulk[[treat]][[prolif]] = AggregateExpression(subset_seurat, 
                                                        assays = "RNA",
                                                        return.seurat = FALSE,
                                                        slot = "counts", # raw counts
                                                        group.by = c("donor_id","pool","replicate"))
    pseudobulk[[treat]][[prolif]]  = as.data.frame(pseudobulk[[treat]][[prolif]] )
    message("Adding metadata")
    
    metadata[[treat]][[prolif]]  = subset_seurat@meta.data %>%
      dplyr::mutate(cols_names = paste(treat,prolif,donor_id,pool,replicate,sep = "_")) %>%
      dplyr::group_by(cols_names) %>%
      dplyr::reframe(count = n(),
                     proliferation_status = prolif,
                     donor_id=donor_id,
                     treatment = treat,
                     pool=pool,
                     replicate=replicate) %>%
      dplyr::distinct() %>%
      dplyr::arrange(tolower(cols_names)) %>% # pseudobulk count columns are in alphabetical order (case insensitive)
      as.data.frame()
    
    colnames(pseudobulk[[treat]][[prolif]] ) =metadata[[treat]][[prolif]] $cols_names
    rm(subset_seurat)
    gc()
  }
  pseudobulk[[treat]] = Reduce(cbind, pseudobulk[[treat]])
  metadata[[treat]] = Reduce(rbind, metadata[[treat]])
  gene = rownames(pseudobulk[[treat]])
  
  rm(input_seurat)
  gc()
}

combined_pseudobulk = Reduce(cbind, pseudobulk)
combined_metadata = Reduce(rbind, metadata)
combined_pseudobulk$gene = gene
length(unique(combined_metadata$donor_id)) # 248
write_tsv(combined_pseudobulk,"../../data/results/5.prepare_differential_expression_pseudobulks/combined_pseudobulk_raw_counts_treatmentxprolifxidxpool_reps.txt")
write_tsv(combined_metadata,"../../data/results/5.prepare_differential_expression_pseudobulks/metadata_treatmentxprolifxidxpool_reps.txt")

##########################treatmentxprolifxidxpool without replicates
pseudobulk = list()
metadata = list()
for(treat in c("untreated","LPS","IFN")){
  message(treat)
  
  message("Reading in data")
  
  # seurat data
  input_seurat = readRDS(paste0("../../data/results/1.QC_v5/",treat,"_filtered_harmony/",treat,"_filtered_harmony.Rds"))
  
  gc()
  
  Layers(input_seurat[["RNA"]])
  DefaultAssay(input_seurat) = "RNA"
  Idents(input_seurat) = "cluster_full"
  # message("Joining layers")
  # #pseudobulk on donor + pool + treatment
  # ## transform to singleCellExperiment class
  # input_seurat = JoinLayers(input_seurat)
  
  pdf(paste0(output_dir,"proliferating_scores_seurat_dimplot.pdf"))
  p1 = scCustomize::DimPlot_scCustom(input_seurat,label = FALSE,group.by = "cluster_full", split.by = "Phase" )
  p1
  pdf(paste0(output_dir,"proliferating_scores_seurat_featureplot.pdf"))
  p2 = scCustomize::FeaturePlot_scCustom(input_seurat,features = c("S.Score","G2M.Score"),split.by = "cluster_full",num_columns = 4)
   
  p2
  dev.off()
  
 
  
  message("Naming proliferation categories")
  
  input_seurat@meta.data$proliferation_status = "Not_proliferating"
  if(treat == "untreated"){
    input_seurat@meta.data[input_seurat@meta.data$cluster_full %in%  c(3,4) , "proliferation_status"] =  "Proliferating"
    
  }
  if(treat == "LPS"){
    input_seurat@meta.data[input_seurat@meta.data$cluster_full %in%  c(4,5,7) , "proliferation_status"] =  "Proliferating"
    
  }
  if(treat == "IFN"){  # much smaller proliferating clusters on IFN treated cells
    input_seurat@meta.data[input_seurat@meta.data$cluster_full %in%  c(4) , "proliferation_status"] =  "Proliferating"
    
  }
  message("Aggregating expression: subset to avoid memory issues")
  
  pseudobulk[[treat]] = list()
  metadata[[treat]] = list()
  for(prolif in c("Not_proliferating","Proliferating")){
    subset_seurat = subset(x = input_seurat, subset = proliferation_status == prolif)
    pseudobulk[[treat]][[prolif]] = AggregateExpression(subset_seurat, 
                                                        assays = "RNA",
                                                        return.seurat = FALSE,
                                                        slot = "counts", # raw counts
                                                        group.by = c("donor_id","pool"))
    pseudobulk[[treat]][[prolif]]  = as.data.frame(pseudobulk[[treat]][[prolif]] )
    message("Adding metadata")
    
    metadata[[treat]][[prolif]]  = subset_seurat@meta.data %>%
      dplyr::mutate(cols_names = paste(treat,prolif,donor_id,pool,sep = "_")) %>%
      dplyr::group_by(cols_names) %>%
      dplyr::reframe(count = n(),
                     proliferation_status = prolif,
                     donor_id=donor_id,
                     treatment = treat,
                     pool=pool) %>%
      dplyr::distinct() %>%
      dplyr::arrange(tolower(cols_names)) %>% # pseudobulk count columns are in alphabetical order (case insensitive)
      as.data.frame()
    
    colnames(pseudobulk[[treat]][[prolif]] ) =metadata[[treat]][[prolif]] $cols_names
    rm(subset_seurat)
    gc()
  }
  pseudobulk[[treat]] = Reduce(cbind, pseudobulk[[treat]])
  metadata[[treat]] = Reduce(rbind, metadata[[treat]])
  gene = rownames(pseudobulk[[treat]])
  
  rm(input_seurat)
  gc()
}

combined_pseudobulk = Reduce(cbind, pseudobulk)
combined_metadata = Reduce(rbind, metadata)
combined_pseudobulk$gene = gene
length(unique(combined_metadata$donor_id)) # 250
write_tsv(combined_pseudobulk,"../../data/results/5.prepare_differential_expression_pseudobulks/combined_pseudobulk_raw_counts_treatmentxprolifxidxpool.txt")
write_tsv(combined_metadata,"../../data/results/5.prepare_differential_expression_pseudobulks/metadata_treatmentxprolifxidxpool.txt")

##########################treatmentxprolifxidxseqsample

pseudobulk = list()
metadata = list()
for(treat in c("untreated","LPS","IFN")){
  message(treat)
  
  message("Reading in data")
  
  # seurat data
  input_seurat = readRDS(paste0("../../data/results/1.QC_v5/",treat,"_filtered_harmony/",treat,"_filtered_harmony.Rds"))
  
  gc()
  
  Layers(input_seurat[["RNA"]])
  DefaultAssay(input_seurat) = "RNA"
  Idents(input_seurat) = "cluster_full"
  # message("Joining layers")
  # #pseudobulk on donor + orig.ident (seq sample) + treatment
  # to check if weird behavior from pool 9 donors is due to a specific sample
  # ## transform to singleCellExperiment class
  # input_seurat = JoinLayers(input_seurat)
  message("Naming proliferation categories")
  
  input_seurat@meta.data$proliferation_status = "Not_proliferating"
  if(treat == "untreated"){
    input_seurat@meta.data[input_seurat@meta.data$cluster_full %in%  c(3,4) , "proliferation_status"] =  "Proliferating"
    
  }
  if(treat == "LPS"){
    input_seurat@meta.data[input_seurat@meta.data$cluster_full %in%  c(4,5,7) , "proliferation_status"] =  "Proliferating"
    
  }
  if(treat == "IFN"){  # much smaller proliferating clusters on IFN treated cells
    input_seurat@meta.data[input_seurat@meta.data$cluster_full %in%  c(4) , "proliferation_status"] =  "Proliferating"
    
  }
  message("Aggregating expression: subset to avoid memory issues")
  
  pseudobulk[[treat]] = list()
  metadata[[treat]] = list()
  for(prolif in c("Not_proliferating","Proliferating")){
    subset_seurat = subset(x = input_seurat, subset = proliferation_status == prolif)
    pseudobulk[[treat]][[prolif]] = AggregateExpression(subset_seurat, 
                                                        assays = "RNA",
                                                        return.seurat = FALSE,
                                                        slot = "counts", # raw counts
                                                        group.by = c("donor_id","orig.ident"))
    pseudobulk[[treat]][[prolif]]  = as.data.frame(pseudobulk[[treat]][[prolif]] )
    message("Adding metadata")
    
    metadata[[treat]][[prolif]]  = subset_seurat@meta.data %>%
      dplyr::mutate(cols_names = paste(treat,prolif,donor_id,orig.ident,sep = "_")) %>%
      dplyr::group_by(cols_names) %>%
      dplyr::reframe(count = n(),
                     proliferation_status = prolif,
                     donor_id=donor_id,
                     treatment = treat,
                     orig.ident=orig.ident,
                     pool=pool) %>%
      dplyr::distinct() %>%
      dplyr::arrange(tolower(cols_names)) %>% # pseudobulk count columns are in alphabetical order (case insensitive)
      as.data.frame()
    
    colnames(pseudobulk[[treat]][[prolif]] ) =metadata[[treat]][[prolif]] $cols_names
    rm(subset_seurat)
    gc()
  }
  pseudobulk[[treat]] = Reduce(cbind, pseudobulk[[treat]])
  metadata[[treat]] = Reduce(rbind, metadata[[treat]])
  gene = rownames(pseudobulk[[treat]])
  
  rm(input_seurat)
  gc()
}

combined_pseudobulk = Reduce(cbind, pseudobulk)
combined_metadata = Reduce(rbind, metadata)
combined_pseudobulk$gene = gene
length(unique(combined_metadata$donor_id)) # 250
write_tsv(combined_pseudobulk,"../../data/results/5.prepare_differential_expression_pseudobulks/combined_pseudobulk_raw_counts_treatmentxprolifxidxseqsample.txt")
write_tsv(combined_metadata,"../../data/results/5.prepare_differential_expression_pseudobulks/metadata_treatmentxprolifxidxseqsample.txt")

