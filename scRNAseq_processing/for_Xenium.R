library(tidyverse)
library(Seurat)
metadata = read.csv("/lustre/scratch123/hgi/projects/otar2065/OTAR2065_differentiation_efficiency/data/results/1.2.Inspect_integration/metadata.csv")
untreated = Read10X_h5("../../../OTAR2065_scRNA_data_processing/data/cellranger/pool17/P17_Diff2_Untreated_A/outs/filtered_feature_bc_matrix.h5", 
                    use.names = TRUE, unique.features = TRUE)

vireo_donor_ID = list()
  files = c( "P17_Diff2_IFN_A", "P17_Diff2_LPS_A" , "P17_Diff2_Untreated_A"   )
  for (file in files) {
    message("... file ", file,"...")
    
    vireo_donor_ID[[paste0(file)]]  = read.table(paste0("../../../OTAR2065_scRNA_data_processing/data/cellranger/",
                                                                 "pool17","/",file,"/vireoOutput/donor_ids.tsv"),
                                                          header = TRUE)
    vireo_donor_ID[[paste0(file)]] = vireo_donor_ID[[paste0(file)]] %>% dplyr::select_if(!names(.) %in% c('doublet_logLikRatio')) # drop column if it exists
    vireo_donor_ID[[paste0(file)]] = vireo_donor_ID[[paste0(file)]] %>%
      dplyr::rename(Barcode = cell,
                    Annotation = donor_id) %>%
      dplyr::select(Barcode,Annotation)
    
  }


metadata_untreated = vireo_donor_ID[["P17_Diff2_Untreated_A"]]  

anyNA(metadata_untreated)
table(colnames(untreated) %in% metadata_untreated$Barcode)
metadata_untreated %>%
  dplyr::mutate(Annotation = "untreated") %>%
  write_csv("../../../OTAR2065_scRNA_data_processing/data/cellranger/pool17/metadata_untreated.csv")

IFN = Read10X_h5("../../../OTAR2065_scRNA_data_processing/data/cellranger/pool17/P17_Diff2_IFN_A/outs/filtered_feature_bc_matrix.h5", 
                       use.names = TRUE, unique.features = TRUE)
metadata_IFN = vireo_donor_ID[["P17_Diff2_IFN_A"]]  
metadata_IFN %>%
  dplyr::mutate(Annotation = "IFN") %>%
  write_csv("../../../OTAR2065_scRNA_data_processing/data/cellranger/pool17/metadata_IFN.csv")

anyNA(metadata_IFN)
table(colnames(IFN) %in% metadata_IFN$Barcode)

LPS = Read10X_h5("../../../OTAR2065_scRNA_data_processing/data/cellranger/pool17/P17_Diff2_LPS_A/outs/filtered_feature_bc_matrix.h5", 
                 use.names = TRUE, unique.features = TRUE)
metadata_LPS = vireo_donor_ID[["P17_Diff2_LPS_A"]]  
metadata_LPS %>%
  dplyr::mutate(Annotation = "LPS") %>%
  write_csv("../../../OTAR2065_scRNA_data_processing/data/cellranger/pool17/metadata_LPS.csv")

anyNA(metadata_LPS)
table(colnames(LPS) %in% metadata_LPS$Barcode)
