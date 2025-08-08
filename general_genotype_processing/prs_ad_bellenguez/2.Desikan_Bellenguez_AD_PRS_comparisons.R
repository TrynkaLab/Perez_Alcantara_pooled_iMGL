# old vs new AD PRS comparisons
#!/usr/bin/env Rscript
.libPaths(c('/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/',
            "/software/teamtrynka/ma23/R4.1/libs",
            "/software/R-4.1.0/lib/R/library"
)) 
library(tidyverse)
library(patchwork)
library(ggrepel)
options(stringsAsFactors=F)
outputPath = "../../data/prs_ad_bellenguez"



AD_new = read.table("../../data/prs_ad_bellenguez/AD_polygenic_hazard.donors.txt",header = T,sep = "\t")

AD_old = read.table("../../data/old_prs_info/AD_polygenic_hazard.cell_lines.recoded.txt",header = T,sep = "\t") %>%
  dplyr::mutate(donor = str_split_i(sampleID,pattern = "-",i=2),
                overall_HR_quartile = case_when(overallHR_quantile <= 0.25 ~ "Q1",
                                                overallHR_quantile <= 0.50 & overallHR_quantile > 0.25 ~ "Q2",
                                                overallHR_quantile <= 0.75 & overallHR_quantile > 0.50 ~ "Q3",
                                                overallHR_quantile > 0.75 ~ "Q4"),
                polygenicHR_quartile = case_when(polygenicHR_quantile <= 0.25 ~ "Q1",
                                                polygenicHR_quantile <= 0.50 & polygenicHR_quantile > 0.25 ~ "Q2",
                                                polygenicHR_quantile <= 0.75 & polygenicHR_quantile > 0.50 ~ "Q3",
                                                polygenicHR_quantile > 0.75 ~ "Q4")) %>%
  dplyr::filter(donor %in% c(AD_new$donor)) %>%
  dplyr::rename(overall_HR_quartile_old = overall_HR_quartile, overall_HR_quantile_old = overallHR_quantile,
                overallHR_old = overallHR,apoeHR_old = apoeHR, polygenicHR_quantile_old = polygenicHR_quantile,
                polygenicHR_old = polygenicHR,polygenic_HR_quartile_old = polygenicHR_quartile)

table(AD_old$overall_HR_quartile_old)
table(AD_new$overallHR_quartile)
table(AD_old$polygenic_HR_quartile_old)
table(AD_new$polygenicHR_quartile)

merged = AD_new %>%
  dplyr::left_join(.,AD_old[c("donor","overall_HR_quartile_old","overall_HR_quantile_old","overallHR_old",
                              "apoeHR_old","polygenicHR_quantile_old","polygenicHR_old","polygenic_HR_quartile_old")]) %>%
  dplyr::mutate(apoediff = abs(apoeHR_old - apoeHR)) %>%
  dplyr::mutate(apoediff = case_when(apoediff !=0 ~ "different",
                                     apoediff ==0 ~ "same")) 
table(merged$apoediff)
# different; 22      same:222 for HRC +1K_10K
# different; 24      same:220 for 1K_10K only
# different; 24      same:220 for HRC only

# the quartiles were calculated from the old one from the overall sample distribution
# so that's why they are not evenly distributed
# but that still doesn't explain why some move from Q1 to Q4 and viceversa

p1 = merged %>%
  ggplot(.,aes(x=overall_HR_quantile_old,y=overallHR_quantile)) +
  theme_minimal() +
  geom_rect(xmin=0, xmax=0.25, ymin=0, ymax=0.25, fill="#CDF5CF", alpha=0.1) + 
  geom_rect(xmin = 0.75, xmax = 1, ymin = 0.75, ymax = 1, fill="#CDF5CF", alpha=0.1) +   
  geom_rect(xmin=0.25, xmax=0.5, ymin=0.25, ymax=0.5, fill="#CDF5CF", alpha=0.1) + 
  geom_rect(xmin = 0.5, xmax = 0.75, ymin = 0.5, ymax = 0.75, fill="#CDF5CF", alpha=0.1) +  
  geom_rect(xmin=0, xmax=0.25, ymin=0.25, ymax=0.5, fill="#FAE0A0", alpha=0.1) + 
  geom_rect(xmin = 0.25, xmax = 0.5, ymin = 0, ymax = 0.25, fill="#FAE0A0", alpha=0.1) +  
  geom_rect(xmin=0.50, xmax=0.75, ymin=0.25, ymax=0.5, fill="#FAE0A0", alpha=0.1) + 
  geom_rect(xmin = 0.25, xmax = 0.5, ymin = 0.50, ymax = 0.75, fill="#FAE0A0", alpha=0.1) +  
  geom_rect(xmin = 0.75, xmax = 1, ymin = 0.5, ymax = 0.75, fill="#FAE0A0", alpha=0.1) +   
  geom_rect(xmin = 0.5, xmax = 0.75, ymin = 0.75, ymax = 1, fill="#FAE0A0", alpha=0.1) + 
  geom_rect(xmin = 0, xmax =0.25, ymin = 0.5, ymax = 0.75, fill="#F7BBA3", alpha=0.1) +   
  geom_rect(xmin = 0.5, xmax = 0.75, ymin = 0, ymax = 0.25, fill="#F7BBA3", alpha=0.1) +  
  geom_rect(xmin = 0.25, xmax =0.5, ymin = 0.75, ymax = 1, fill="#F7BBA3", alpha=0.1) +   
  geom_rect(xmin = 0.75, xmax = 1, ymin = 0.25, ymax = 0.5, fill="#F7BBA3", alpha=0.1) + 
  geom_rect(xmin = 0, xmax =0.25, ymin = 0.75, ymax = 1, fill="#F77979", alpha=0.1) +   
  geom_rect(xmin = 0.75, xmax =1, ymin = 0, ymax = 0.25, fill="#F77979", alpha=0.1) +   
  geom_point() 
  
pdf(paste0(outputPath,"/old_vs_new_PRS_overallHR_quantile.pdf"),width = 5,height = 5)
p1 + ggtitle("Old vs new PRS overall HR quantile")
dev.off()

p2 = merged %>%
  dplyr::mutate(tolabel = case_when(apoediff == "different" ~ donor,
                                     apoediff == "same" ~ NA)) %>%
  ggplot(.,aes(x=apoeHR_old,y=apoeHR,label = tolabel, color = apoediff)) +
  geom_point() +
  geom_text_repel(col = "black",max.overlaps = 25) +
  theme_bw()

pdf(paste0(outputPath,"/old_vs_new_PRS_APOEHR.pdf"),width = 5,height = 5)
p2 + ggtitle("New vs old PRS APOE HR")
dev.off()

p3 = merged %>%

  ggplot(.,aes(x=polygenicHR_quantile_old,y=polygenicHR_quantile)) +
  geom_point() +
  theme_bw()

 p4 = merged %>%
  
  ggplot(.,aes(x=overallHR_old,y=overallHR)) +
  geom_point() +
  theme_bw()
 
 p5 = merged %>%
   
   ggplot(.,aes(x=polygenicHR_old,y=polygenicHR)) +
   geom_point() +
   theme_bw()

table(merged$overallHR_quartile,merged$overall_HR_quartile_old)
table(merged$polygenicHR_quartile,merged$polygenic_HR_quartile_old)
head(merged[,c("polygenicHR_quartile","polygenic_HR_quartile_old")])
table(merged$polygenic_HR_quartile)
table(merged$polygenic_HR_quartile_old)
# 
# ggscatter(merged, y = "overallHR_quantile_AD", x = "overallHR_quantile_PD",
#           size = 3, alpha = 0.6) + border() + geom_abline(intercept = 0, slope = 1)
# 
# 
# ggsave(
#   filename = "../../data/trait_comparisons/overallHR_quantile_ADvsPD.png",
#   device = "png",
#   height = 4,
#   width = 9,
#   units = "in",
#   dpi = 400
# )
# 
# # subset for those differentiated
# 
# merged = merge(AD[c("Line","overallHR_quantile","Neuroseq_Differentiated")],PD[c("Line","overallHR_quantile")], by = "Line")
# merged = merged[merged$Neuroseq_Differentiated == "yes",]
# colnames(merged) = c("Line","overallHR_quantile_AD","Neuroseq_Differentiated","overallHR_quantile_PD")
# 
# 
# ggscatter(merged, y = "overallHR_quantile_AD", x = "overallHR_quantile_PD",
#           size = 3, alpha = 0.6) + border() + geom_abline(intercept = 0, slope = 1)
# 
# 
# ggsave(
#   filename = "../../data/trait_comparisons/overallHR_quantile_ADvsPD_differentiated_in_Neuroseq.png",
#   device = "png",
#   height = 4,
#   width = 9,
#   units = "in",
#   dpi = 400
# )
# 
# # Save merged table with PD and AD info
# 
# merged = merge(AD,PD[c("Line","overallHR","overallHR_quantile","quartile")], by = "Line")
# 
# # rename columns
# colnames(merged)[c(39,40,41)] = c("overallHR-PD","overallHR_quantile-PD","overallHR_quartile-PD")
# colnames(merged)[c(3,6,7,
#                    9,10,11)] = c("polygenicHR-AD","apoeHR-AD","overallHR-AD",
#                                  "polygenicHR_quantile-AD","overallHR_quantile-AD","overallHR_quartile-AD")
# 
# merged = merged[c(1:2,8,4:6,3,9,7,10,11,39,40,41,12:38)]
# 
# # saving all info
# # retaining only Disease-free, feeder-free, European
# merged_all = merged[which(merged$Disease.Status == "Normal"),]
# merged_all = merged_all[which(merged_all$Culture == "Feeder-free"),]
# length(unique(merged_all$donor))
# 
# merged_all = merged_all[which(merged_all$Predicted.Population== "European"),]
# length(unique(merged_all$donor))
# 
# write.table(merged_all,file = "../../data/trait_comparisons/allinfo_hipsci_PD_AD_PRS_non-disease_feeder-free_European.txt",sep = "\t",quote = F, row.names = F, col.names = T)
# 
# 
# # don't save neuroseq information
# merged = merged[c(1:34)]
# merged = unique(merged)
# write.table(merged,file = "../../data/trait_comparisons/hipsci_PD_AD_PRS_all.txt",sep = "\t",quote = F, row.names = F, col.names = T)
# 
# # retaining only Disease-free, feeder-free, European
# merged = merged[which(merged$Disease.Status == "Normal"),]
# merged = merged[which(merged$Culture == "Feeder-free"),]
# length(unique(merged$donor))
# 
# merged = merged[which(merged$Predicted.Population== "European"),]
# length(unique(merged$donor))
# 
# write.table(merged,file = "../../data/trait_comparisons/hipsci_PD_AD_PRS_non-disease_feeder-free_European.txt",sep = "\t",quote = F, row.names = F, col.names = T)
# 
# 
# write.table(merged,file = "../../data/trait_comparisons/hipsci_PD_AD_PRS_nonDiseased_feeder-free.txt",sep = "\t",quote = F, row.names = F, col.names = T)