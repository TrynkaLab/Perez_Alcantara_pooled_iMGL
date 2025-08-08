# automatic annotation with singleR
# run conda activate /software/teamtrynka/conda/otar2065
# then don't run in Rstudio, but in interactive R session
# the latter does not work because it can't find libhdf5_hl.so.100 
# if I try to load it from somewhere else with dyn.load it crashes
# run on windows on subset and compare to seurat transfer, after fixing number of genes used
library(patchwork)
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(SingleR)
library(future)
library(scuttle)
library(scater)
library(scran)
library(SCpubr)
plan("multisession", workers = 2)
plan()
options(future.globals.maxSize = (100000 * 1024 ^ 2)) #(~100 Gb)

directory = "../../data/results/3.annotation_singleR/"
dir.create(directory, recursive = T)

# query data
# see this for why to start from raw counts https://github.com/LTLA/SingleR/issues/98
# input_merged = readRDS("../../data/results/1.1.subset_for_visualization/merged_filtered_harmony_subset/merged_filtered_harmony_50k_subset_seuratv5.Rds")
input_merged = readRDS("../../data/results/1.1.subset_for_visualization/merged_filtered_harmony_subset/merged_filtered_harmony_50k_subset_seuratv5.Rds")
input_merged@meta.data$cell = rownames(input_merged@meta.data)
# needs logNormCounts  - my data has that already
my_metadata = input_merged@meta.data
my_lognormcounts = input_merged[["RNA"]]$data 
## reference data -  Mancuso ex vivo (in mice) matured iPSC derived microglia human
ex_vivo_mancuso=readRDS("../../../resources/Mancuso_Nat_Neurosci_2019/ex_vivo_mancuso_filtered_annotated_allSamples_renamed_clusters.rds")
# Reduce automated annotation to microglia or perivascular macrophages
ex_vivo_mancuso$reduced_renamed_clusters = "Microglia"
ex_vivo_mancuso@meta.data[ex_vivo_mancuso@meta.data$seurat_clusters %in% c(1,4),"reduced_renamed_clusters"] = "Perivascular_macrophages"

ex_vivo_mancuso = as.SingleCellExperiment(ex_vivo_mancuso)

# SingleR() expects reference datasets to be normalized and log-transformed.
ex_vivo_mancuso = logNormCounts(ex_vivo_mancuso)
pred.mancuso = SingleR(test=input_merged, ref=ex_vivo_mancuso, 
                      labels=ex_vivo_mancuso$renamed_clusters,)
table(pred.mancuso$labels)

png("../../data/results/3.singleR_labelTransfer/Mancuso_xenografts_score_heatmap.png", 
    width = 12, height = 3, units = "in", res = 400)
plotScoreHeatmap(pred.mancuso)
dev.off()
png("../../data/results/3.singleR_labelTransfer/Mancuso_xenografts_delta_distribution.png", 
    width = 9, height = 9, units = "in", res = 400)
plotDeltaDistribution(pred.mancuso, ncol = 3)
dev.off()
summary(is.na(pred.mancuso$pruned.labels))
all.markers = metadata(pred.mancuso)$de.genes
media_merged$xenograft_labels = pred.mancuso$labels


## Give always same colors per category

group.colors = c(`Microglia` = "#73937E",  `Perivascular macrophages` = "#218380")


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
