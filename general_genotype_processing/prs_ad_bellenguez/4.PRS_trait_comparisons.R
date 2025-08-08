#!/usr/bin/env Rscript
library(reshape2)
library(ggplot2)
library(plyr)
library(summarytools)
library(ggpubr)
library(ggExtra)

options(stringsAsFactors=F)



AD = read.table("../../data/AD/AD_output/AD_cell_lines_other_info.txt",header = T,sep = "\t")

PD = read.table("../../data/PD/PD_output/PD_cell_lines_other_info.txt",header = T,sep = "\t")


merged = merge(AD[c("Line","overallHR_quantile")],PD[c("Line","overallHR_quantile")], by = "Line")
colnames(merged) = c("Line","overallHR_quantile_AD","overallHR_quantile_PD")


ggscatter(merged, y = "overallHR_quantile_AD", x = "overallHR_quantile_PD",
          size = 3, alpha = 0.6) + border() + geom_abline(intercept = 0, slope = 1)


ggsave(
  filename = "../../data/trait_comparisons/overallHR_quantile_ADvsPD.png",
  device = "png",
  height = 4,
  width = 9,
  units = "in",
  dpi = 400
)

# subset for those differentiated

merged = merge(AD[c("Line","overallHR_quantile","Neuroseq_Differentiated")],PD[c("Line","overallHR_quantile")], by = "Line")
merged = merged[merged$Neuroseq_Differentiated == "yes",]
colnames(merged) = c("Line","overallHR_quantile_AD","Neuroseq_Differentiated","overallHR_quantile_PD")


ggscatter(merged, y = "overallHR_quantile_AD", x = "overallHR_quantile_PD",
          size = 3, alpha = 0.6) + border() + geom_abline(intercept = 0, slope = 1)


ggsave(
  filename = "../../data/trait_comparisons/overallHR_quantile_ADvsPD_differentiated_in_Neuroseq.png",
  device = "png",
  height = 4,
  width = 9,
  units = "in",
  dpi = 400
)

# Save merged table with PD and AD info

merged = merge(AD,PD[c("Line","overallHR","overallHR_quantile","quartile")], by = "Line")

# rename columns
colnames(merged)[c(39,40,41)] = c("overallHR-PD","overallHR_quantile-PD","overallHR_quartile-PD")
colnames(merged)[c(3,6,7,
                   9,10,11)] = c("polygenicHR-AD","apoeHR-AD","overallHR-AD",
                                 "polygenicHR_quantile-AD","overallHR_quantile-AD","overallHR_quartile-AD")

merged = merged[c(1:2,8,4:6,3,9,7,10,11,39,40,41,12:38)]

# saving all info
# retaining only Disease-free, feeder-free, European
merged_all = merged[which(merged$Disease.Status == "Normal"),]
merged_all = merged_all[which(merged_all$Culture == "Feeder-free"),]
length(unique(merged_all$donor))

merged_all = merged_all[which(merged_all$Predicted.Population== "European"),]
length(unique(merged_all$donor))

write.table(merged_all,file = "../../data/trait_comparisons/allinfo_hipsci_PD_AD_PRS_non-disease_feeder-free_European.txt",sep = "\t",quote = F, row.names = F, col.names = T)


# don't save neuroseq information
merged = merged[c(1:34)]
merged = unique(merged)
write.table(merged,file = "../../data/trait_comparisons/hipsci_PD_AD_PRS_all.txt",sep = "\t",quote = F, row.names = F, col.names = T)

# retaining only Disease-free, feeder-free, European
merged = merged[which(merged$Disease.Status == "Normal"),]
merged = merged[which(merged$Culture == "Feeder-free"),]
length(unique(merged$donor))

merged = merged[which(merged$Predicted.Population== "European"),]
length(unique(merged$donor))

write.table(merged,file = "../../data/trait_comparisons/hipsci_PD_AD_PRS_non-disease_feeder-free_European.txt",sep = "\t",quote = F, row.names = F, col.names = T)


write.table(merged,file = "../../data/trait_comparisons/hipsci_PD_AD_PRS_nonDiseased_feeder-free.txt",sep = "\t",quote = F, row.names = F, col.names = T)