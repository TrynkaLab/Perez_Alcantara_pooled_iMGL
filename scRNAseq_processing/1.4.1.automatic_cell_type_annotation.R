# label transfer for automatic cell type annotation
# https://satijalab.org/seurat/articles/integration_mapping
# singleR causing problems with Seuratv5
library(Seurat)
library(tidyverse)
library(SCpubr)
library(patchwork)
library(future)
library(scuttle)
library(scater)
source("./helpers.R")
# set this option when analyzing large datasets
options(future.globals.maxSize = 120000 * 1024^2) # 120Gb
# plan("multiprocess", workers = 14)
# plan()
options(stringsAsFactors = FALSE)

# for new seurat assay slots
options(Seurat.object.assay.version = "v5")

directory = "../../data/results/1.4.markers/1.4.1.automatic_cell_type_annotation"
dir.create(directory, recursive = T)
# 
# library(Seurat)
# library(SingleR)

# library(future)
# library(scran)


### main -------

# Query data: my  data
media_merged= readRDS("../../data/results/1.QC_v5/untreated_filtered_harmony/untreated_filtered_harmony.Rds")
input_subset = subset(x = media_merged, downsample = 50000) # downsample to 50k cells
rm(media_merged)
input_subset[["RNA"]] = JoinLayers(input_subset[["RNA"]])

# Reference data: Mancuso 2024 xenografted iMicroglia #######################
load("../../../resources/Mancuso_2024/data/unfiltered_seurat/MancusoFattorelliMartinez2024_AllMyeloid.RData.gz")

# reference named remind_all.seurat

### seurat automatic annotation via label transfer
anchors = FindTransferAnchors(reference = remind_all.seurat, 
                                        query = input_subset, 
                                        dims = 1:30,
                              reference.assay = "integrated_SCT_final",
                              normalization.method = "SCT",
                                        reference.reduction = "pca_integrated_SCT")

predictions = TransferData(anchorset = anchors, 
                           refdata = remind_all.seurat$cell.state_v1, dims = 1:30)
input_subset = AddMetaData(input_subset, metadata = predictions)
table(input_subset$predicted.id)

## Give always same colors per category

group.colors = c(`Homeostatic (HM)` = "#73937E",  
                 `CNS-associated\n macrophages (CAM)` = "#AF999A", 
                 `Interferon\n response (IRM)` = "#8f04d8", 
                 `Cytokines\n response (CRM)` = "#8f84d8", 
                 `Antigen-presenting\n response (HLA)`  = "#70037E",
                 `Proliferating` = "#EF798A", 
                 `Ribosomal\n response (RM)` = "#EF098A", 
                 `Disease\n associated (DAM)` ="#0FA3B1",
                 `Doublets/LowQ`="#AF000A",
                 `Other Myeloid`= "#AA999A",
                 `Transitioning CRM` = "#AB999A")

input_subset@meta.data %>%
  dplyr::count(predicted.id)  %>% mutate(Percent = n / sum(n)*100) %>%
  ggplot(aes(x=predicted.id)) +
  geom_bar() +
  theme_minimal()

SCpubr::do_DimPlot(sample = input_subset,
                   group.by = "predicted.id",
                   colors.use = group.colors)






long = colData(media_merged) %>%
  as_tibble() %>%
  group_by(orig.ident) %>%
  dplyr::count(xenograft_labels)  %>% mutate(Percent = n / sum(n)*100)


barplot2 = ggplot(data=long, aes(x=orig.ident, y=Percent, fill=xenograft_labels)) +
  geom_bar(stat="identity", col="black") + theme_bw() + 
  xlab("Library") + scale_fill_manual(values=group.colors) + 
  theme(axis.text.x = element_text(face = "bold"), axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"), 
        legend.text = element_text(face = "bold"), legend.title = element_blank(),
        plot.title = element_text(face = "bold")) +
  ggtitle("Label transfer: xenografted iPSC-microglia")

png("../../data/results/3.singleR_labelTransfer/label_transfer_barplots_per_library_Fig5c.png", 
    width = 8, height = 5, units = "in", res = 400)
barplot1 / barplot2
dev.off()
