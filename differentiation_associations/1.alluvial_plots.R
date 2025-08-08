# alluvial plots iPSC - premacs - microglia
# all pools
# Checking donor proportions and make alluvial plots
# with WGS and vireo estimates at differentiation stages
# all pools
library(patchwork)
library(tidyverse)
library(stringr)
library(psych)
library(corrplot)
library(ggalluvial)
library(ggfittext)
library(ggpubr)
library(pheatmap)
library(cowplot)
library(grid)
library(gridExtra)
source("./functions.R")

directory = "../../data/results/1.alluvial_plots/"
dir.create(directory, recursive = T)

source("functions.R")

# make corplots scRNA-seq  - WGS for those coming from same preMacs for deconvolution support

#########
# read in info on donors in pools
# from vireo
# unassigned and doublets are assigned to the best one or two lines respectively
microglia_prop = load_sc_donors(pools = paste0("pool",2:17)) %>%
  dplyr::select(!set) %>%
  dplyr::mutate(type = "scRNA-seq")
# WGS: iPSC and precursors prop, also migration and phagocytosis assays
wgs_prop=list()

for(pool in  paste0("pool",2:17)){
  files = list.files(paste0("../../data/w/",pool))
  for(fil in files){
    name = stringr::str_remove( stringr::str_remove(fil,"_w_estimate.txt"),".txt")
    wgs_prop[[name]] = read.table(paste0("../../data/w/",pool,"/",fil))
    colnames(wgs_prop[[name]]) = c("Line","prop")
    wgs_prop[[name]]$pool = pool
    wgs_prop[[name]]$sample = name
  }
}

wgs_prop = do.call("rbind",wgs_prop)

table(wgs_prop$sample,wgs_prop$pool)

wgs_prop = wgs_prop %>%
  mutate(sample = if_else(sample %in% c("iPSC_pool","iPSC_WGS","iPSC_D0","D0_iPSC"),"day0_iPSC",sample)) %>%
  mutate(prop = if_else(prop < 0,0,prop)) %>%
  mutate(stage = case_when(str_detect(sample,"phago") ~ "microglia",
                           str_detect(sample,"iPSC")~ "iPSC",
                           str_detect(sample,"migration") ~ "microglia",
                            str_detect(sample,pattern = regex("premac",ignore_case = TRUE)) ~ "preMac",
                           str_detect(sample,"chrmX") ~ "microglia",
                           str_detect(sample,"_cell") ~ "microglia",
                           str_detect(sample,pattern = regex("phag",ignore_case = TRUE)) ~ "microglia",
                           str_detect(sample,pattern = regex("mig",ignore_case = TRUE)) ~ "microglia",
                           str_detect(sample,"UNTR") ~ "microglia",
                           str_detect(sample,"LPS")  ~ "microglia",
                           str_detect(sample,"INF")  ~ "microglia",
                           str_detect(sample,"IFN")  ~ "microglia"),
         type = "WGS")
if(anyNA(wgs_prop) == TRUE){message("There are some unclassified samples!")}
merged = rbind(wgs_prop,microglia_prop)
rm(wgs_prop,microglia_prop)

# remove sh5y5y
merged = merged %>%
  dplyr::filter(Line!="sh5y5y")

merged %>%
  dplyr::filter(pool %in% c(paste0("pool",2:11),paste0("pool",13:17))) %>%  ### pool 12 is an outlier
  dplyr::filter(!sample %in% c(paste0("AG_",c(20,40,60)))) %>%  ### remove chromium 10 replicates with more cells loaded because it's technical duplicate
  
  write.table(.,paste0(directory,"pools2-11_13-17_changing_props_iPSC_preMacs_microglia_WGS_sc.txt"),
              row.names = FALSE,col.names = TRUE,quote = FALSE,sep = "\t")

######### alluvial plots

### there are inconsistent names so decided to do this manually as samples came out of sequencing
plist = list()
plist[[1]] = merged %>%
  dplyr::filter(pool == "pool2" & sample %in% c("pool2_day0_iPSC","pool2_D7-post-thaw_iPSC",
                                                "pool2_day29_preMAC_factory1","pool2_day57_preMAC_factory2b","pool2_day109_preMAC_factory2b")) %>%
  dplyr::mutate(sample = factor(sample,
                                levels = c("pool2_day0_iPSC","pool2_D7-post-thaw_iPSC",
                                           "pool2_day29_preMAC_factory1","pool2_day57_preMAC_factory2b","pool2_day109_preMAC_factory2b"), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 2 -  I,II",
                            subtitle = "During iPSC freeze-thawing and preMAC aging") +
  theme(axis.title.x = element_blank(),
        legend.position="none")



plist[[1]]


####### important - rerun this
plist[[2]] = merged %>%
  dplyr::filter(pool == "pool2" & sample %in% c("P2-IV_D0_iPSC_310821" ,
                                                "P2-IV_D28_PreMac_041021",
                                                "P2-IV_D60_PreMac_051121", # untreated
                                                "P2-IV_D84_PreMac_291121"
                                               )) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("P2-IV_D0_iPSC_310821" ,
                                          "P2-IV_D28_PreMac_041021",
                                          "P2-IV_D60_PreMac_051121", # untreated
                                          "P2-IV_D84_PreMac_291121" 
                                ), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 2 -  IV",
                            subtitle = "From iPSC to preMac aging") +
  theme(axis.title.x = element_blank(),
        legend.position="none")



plist[[2]]

########

plist[[3]] = merged %>%
  dplyr::filter(pool == "pool2" & sample %in% c("pool2_day0_iPSC", # I
                                                "pool2_D7-post-thaw_iPSC",  # II
                                                "P2-III_D0_iPSC_300521",  # III
                                                "P2-IV_D0_iPSC_310821"   # IV
                                               
  )) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("pool2_day0_iPSC", # I
                                          "pool2_D7-post-thaw_iPSC",  # II
                                          "P2-III_D0_iPSC_300521",  # III
                                          "P2-IV_D0_iPSC_310821"   # IV
                                          
                                ), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "iPSC comparisons pool 2 -  I, II, III, IV",
                            subtitle = " iPSC comparisons") +
  theme(axis.title.x = element_blank(),
        legend.position="none")



plist[[3]]



####### important - rerun this
p1 = merged %>%
  dplyr::filter(pool == "pool2" & sample %in% c("P2-IV_D0_iPSC_310821" ,
                                                "P2-IV_D28_PreMac_041021",
                                                "h2IVdiff4_chrmX_WGS" # chromium X microglia WGS
  )) %>%
  dplyr::mutate(sample = factor(sample,
                                levels = c("P2-IV_D0_iPSC_310821" ,
                                           "P2-IV_D28_PreMac_041021",
                                           "h2IVdiff4_chrmX_WGS"
                                ), 
                                ordered = TRUE)) %>% 
  alluvial_plot() +ggtitle( "WGS")  +
  theme(axis.title.x = element_blank(),
        legend.position="none")


p2 = merged %>%
  dplyr::filter(pool == "pool2" & sample %in% c("P2-IV_D0_iPSC_310821" ,
                                                "P2-IV_D60_PreMac_051121",
                                                "h2IVdiff9_RGM_plst" # scRNA-seq
  )) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("P2-IV_D0_iPSC_310821" ,
                                          "P2-IV_D60_PreMac_051121",
                                          "h2IVdiff9_RGM_plst"
                                ), 
                                ordered = TRUE)) %>% 
  alluvial_plot() +ggtitle( "scRNA-seq")  +
  theme(axis.title.x = element_blank(),
        legend.position="none")


p3 = merged %>%
  dplyr::filter(pool == "pool2" & sample %in% c("P2-IV_D0_iPSC_310821" ,
                                                "P2-IV_D28_PreMac_041021", # closest premac match
                                                "microglia_migration_s81_untreated", # WGS untreated migration
                                                "microglia_migration_s101_IFN"   # IFN untreated migration
  )) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("P2-IV_D0_iPSC_310821" ,
                                          "P2-IV_D28_PreMac_041021", # closest premac match
                                          "microglia_migration_s81_untreated", # WGS untreated migration
                                          "microglia_migration_s101_IFN"  
                                ), 
                                ordered = TRUE)) %>% 
  alluvial_plot() +ggtitle( "scRNA-seq")  +
  theme(axis.title.x = element_blank(),
        legend.position="none")




plist[[4]] = ((p1 + p3) / (p2 + patchwork::plot_spacer()))   +
  patchwork::plot_annotation(title = "Changing donor proportions for pool 2 -  IV",
                                      subtitle = "From iPSC to microglia")
plist[[4]]

########

  p1 = merged %>%
  dplyr::filter(pool == "pool3" & sample %in% c("P3_iPSC_pool_WGS" ,
                                                "P3_d35_preMAC_WGS",
                                                "h3_diff2_untreated", # untreated scRNA-seq
                                                "h3_diff2_IFNg",  # IFN scRNA-seq
                                                "h3_diff2_LPS") # LPS
  ) %>%
  dplyr::mutate(sample = factor(sample,
                                levels = c("P3_iPSC_pool_WGS" ,
                                           "P3_d35_preMAC_WGS",
                                           "h3_diff2_untreated", # untreated
                                           "h3_diff2_IFNg",  # IFN
                                           "h3_diff2_LPS"), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle( "scRNA-seq")  +
  theme(axis.title.x = element_blank(),
        legend.position="none")

p2 = merged %>%
  dplyr::filter(pool == "pool3" & sample %in% c("P3_iPSC_pool_WGS" ,
                                                "P3_d35_preMAC_WGS",
                                                "microglia_migration_s135_untreated_merged", # untreated WGS
                                                "microglia_migration_s147_LPS_merged") # LPS
  ) %>%
  dplyr::mutate(sample = factor(sample,
                                levels = c("P3_iPSC_pool_WGS" ,
                                           "P3_d35_preMAC_WGS",
                                           "microglia_migration_s135_untreated_merged", # untreated WGS
                                           "microglia_migration_s147_LPS_merged"), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle( "WGS migration")  +
  theme(axis.title.x = element_blank(),
        legend.position="none")

p3 = merged %>%
  dplyr::filter(pool == "pool3" & sample %in% c("P3_iPSC_pool_WGS" ,
                                                "P3_d57_preMAC_WGS",
                                                "microglia_phago_s3D-1_untreated_unsorted" , # untreated phago
                                                "microglia_phago_s3E-1_IFN_unsorted" ,
                                                "microglia_phago_s3F-1_LPS_unsorted" ) # LPS
  ) %>%
  dplyr::mutate(sample = factor(sample,
                                levels = c("P3_iPSC_pool_WGS" ,
                                           "P3_d57_preMAC_WGS",
                                           "microglia_phago_s3D-1_untreated_unsorted" , # untreated phago
                                           "microglia_phago_s3E-1_IFN_unsorted" ,
                                           "microglia_phago_s3F-1_LPS_unsorted" ), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle( "WGS phagocytosis")  +
  theme(axis.title.x = element_blank(),
        legend.position="none")

plist[[5]] = ((p2 + p1) / (p3 + patchwork::plot_spacer()))   +
  patchwork::plot_annotation(title = "Changing donor proportions for pool 3",
                             subtitle = "From iPSC to microglia. Best match of preMAC for migration (day 39)")

plist[[5]]


# pool 4

plist[[6]] = merged %>%
  dplyr::filter(pool == "pool4" & sample %in% c("P4_iPSC_pool_WGS" ,
                                                "P4_d35_preMAC_WGS",
                                                "P4_d57_preMAC_WGS"
                                               )) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("P4_iPSC_pool_WGS" ,
                                          "P4_d35_preMAC_WGS",
                                          "P4_d57_preMAC_WGS"), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 4",
                            subtitle = "iPSC to preMACs")



plist[[6]]


p1 = merged %>%
  dplyr::filter(pool == "pool4" & sample %in% c("P4_iPSC_pool_WGS" ,
                                                "P4_d35_preMAC_WGS",
                                                "h4_diff2_untreated", # untreated scRNA-seq
                                                "h4_diff2_IFNg",  # IFN scRNA-seq
                                                "h4_diff2_LPS") # LPS
  ) %>%
  dplyr::mutate(sample = factor(sample,
                                levels = c("P4_iPSC_pool_WGS" ,
                                           "P4_d35_preMAC_WGS",
                                           "h4_diff2_untreated", # untreated
                                           "h4_diff2_IFNg",  # IFN
                                           "h4_diff2_LPS"), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle( "scRNA-seq")  +
  theme(axis.title.x = element_blank(),
        legend.position="none")

p2 = merged %>%
  dplyr::filter(pool == "pool4" & sample %in% c("P4_iPSC_pool_WGS" ,
                                                "P4_d35_preMAC_WGS",
                                                "microglia_migration_s159_untreated_merged",
                                                "microglia_migration_s171_LPS_merged"  ) # LPS
  ) %>%
  dplyr::mutate(sample = factor(sample,
                                levels = c("P4_iPSC_pool_WGS"  ,
                                           "P4_d35_preMAC_WGS",
                                           "microglia_migration_s159_untreated_merged",
                                           "microglia_migration_s171_LPS_merged"  ), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle( "WGS migration")  +
  theme(axis.title.x = element_blank(),
        legend.position="none")

p3 = merged %>%
  dplyr::filter(pool == "pool4" & sample %in% c("P4_iPSC_pool_WGS" ,
                                                "P4_d57_preMAC_WGS",
                                                "microglia_phago_s4D-1_untreated_unsorted" , # untreated phago
                                                "microglia_phago_s4E-1_IFN_unsorted" ,
                                                "microglia_phago_s4F-1_LPS_unsorted" ) # LPS
  ) %>%
  dplyr::mutate(sample = factor(sample,
                                levels = c("P4_iPSC_pool_WGS" ,
                                           "P4_d57_preMAC_WGS",
                                           "microglia_phago_s4D-1_untreated_unsorted" , # untreated phago
                                           "microglia_phago_s4E-1_IFN_unsorted" ,
                                           "microglia_phago_s4F-1_LPS_unsorted" ), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle( "WGS phagocytosis")  +
  theme(axis.title.x = element_blank(),
        legend.position="none")

plist[[7]] = ((p2 + p1) / (p3 + patchwork::plot_spacer()))   +
  patchwork::plot_annotation(title = "Changing donor proportions for pool 4",
                             subtitle = "From iPSC to microglia. Best match of preMAC for migration (day 42)")

plist[[7]]

### pool 5 ######

plist[[8]] = merged %>%
  dplyr::filter(pool == "pool5" & sample %in% c("P5_iPSC_WGS"  ,
                                                "P5_d35_preMAC_WGS" ,  
                                                "P5_D42_PreMac_160522" ,
                                                "P5_D50_PreMac_240522"  
  )) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("P5_iPSC_WGS"  ,
                                          "P5_d35_preMAC_WGS" ,  
                                          "P5_D42_PreMac_160522" ,
                                          "P5_D50_PreMac_240522"  ), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 5",
                            subtitle = "iPSC to preMACs")



plist[[8]]


p1 = merged %>%
  dplyr::filter(pool == "pool5" & sample %in% c("P5_iPSC_WGS"  ,
                                                "P5_d35_preMAC_WGS" ,  
                                                "P5_diff2_untreated", # untreated scRNA-seq
                                                "P5_diff2_IFN" ,
                                                "P5_diff2_LPS"    ) # LPS
  ) %>%
  dplyr::mutate(sample = factor(sample,
                                levels = c("P5_iPSC_WGS"  ,
                                           "P5_d35_preMAC_WGS" ,  
                                           "P5_diff2_untreated", # untreated scRNA-seq
                                           "P5_diff2_IFN" ,
                                           "P5_diff2_LPS" ), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle( "scRNA-seq")  +
  theme(axis.title.x = element_blank(),
        legend.position="none")

p2 = merged %>%
  dplyr::filter(pool == "pool5" & sample %in% c("P5_iPSC_WGS"  ,
                                                "P5_D42_PreMac_160522" ,
                                                "P5_Mig15_UNT_seeded",
                                                "P5_Mig15_IFN_seeded"  ,
                                                "P5_Mig15_LPS_seeded"   ) # LPS
  ) %>%
  dplyr::mutate(sample = factor(sample,
                                levels = c("P5_iPSC_WGS"   ,
                                           "P5_D42_PreMac_160522" ,
                                           "P5_Mig15_UNT_seeded",
                                           "P5_Mig15_IFN_seeded"  ,
                                           "P5_Mig15_LPS_seeded"     ), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle( "WGS migration")  +
  theme(axis.title.x = element_blank(),
        legend.position="none")

p3 = merged %>%
  dplyr::filter(pool == "pool5" & sample %in% c("P5_iPSC_WGS"  ,
                                                "P5_D50_PreMac_240522"   ,
                                                "microglia_phago_s5A-1_untreated_unsorted" ,
                                                "microglia_phago_s5B-1_IFN_unsorted" ,
                                                "microglia_phago_s5C-1_LPS_unsorted"   ) # LPS
  ) %>%
  dplyr::mutate(sample = factor(sample,
                                levels = c("P5_iPSC_WGS"  ,
                                           "P5_D50_PreMac_240522"   ,
                                           "microglia_phago_s5A-1_untreated_unsorted" ,
                                           "microglia_phago_s5B-1_IFN_unsorted" ,
                                           "microglia_phago_s5C-1_LPS_unsorted"  ), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle( "WGS phagocytosis")  +
  theme(axis.title.x = element_blank(),
        legend.position="none")

plist[[9]] = ((p2 + p1) / (p3 + patchwork::plot_spacer()))   +
  patchwork::plot_annotation(title = "Changing donor proportions for pool 5",
                             subtitle = "From iPSC to microglia.")

plist[[9]]

######### pool 6 #########

plist[[10]] = merged %>%
  dplyr::filter(pool == "pool6" & sample %in% c("P6_iPSC_WGS"  ,
                                                "P6_d35_preMAC_WGS" ,
                                                "P6_D39_PreMac_130522" ,
                                                "P6_D50_PreMac_240522"  
  )) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("P6_iPSC_WGS"  ,
                                          "P6_d35_preMAC_WGS" ,  
                                          "P6_D39_PreMac_130522" ,
                                          "P6_D50_PreMac_240522" ), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 6",
                            subtitle = "iPSC to preMACs")



plist[[10]]


p1 = merged %>%
  dplyr::filter(pool == "pool6" & sample %in% c("P6_iPSC_WGS"  ,
                                                "P6_d35_preMAC_WGS" ,  
                                                "P6_diff2_untreated", # untreated scRNA-seq
                                                "P6_diff2_IFN" ,
                                                "P6_diff2_LPS"    ) 
  ) %>%
  dplyr::mutate(sample = factor(sample,
                                levels = c("P6_iPSC_WGS"  ,
                                           "P6_d35_preMAC_WGS"    ,  
                                           "P6_diff2_untreated", # untreated scRNA-seq
                                           "P6_diff2_IFN" ,
                                           "P6_diff2_LPS" ), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle( "scRNA-seq")  +
  theme(axis.title.x = element_blank(),
        legend.position="none")

p2 = merged %>%
  dplyr::filter(pool == "pool6" & sample %in% c("P6_iPSC_WGS"  ,
                                                "P6_D39_PreMac_130522"      ,
                                                "P6_Mig14_UNT_seeded",
                                                "P6_Mig14_IFN_seeded"  ,
                                                "P6_Mig14_LPS_seeded"   ) 
  ) %>%
  dplyr::mutate(sample = factor(sample,
                                levels = c("P6_iPSC_WGS"   ,
                                           "P6_D39_PreMac_130522"     ,
                                           "P6_Mig14_UNT_seeded",
                                           "P6_Mig14_IFN_seeded"  ,
                                           "P6_Mig14_LPS_seeded"     ), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle( "WGS migration")  +
  theme(axis.title.x = element_blank(),
        legend.position="none")

p3 = merged %>%
  dplyr::filter(pool == "pool6" & sample %in% c("P6_iPSC_WGS"   ,
                                                "P6_D50_PreMac_240522"   ,
                                                "microglia_phago_s6A-1_untreated_unsorted",
                                                "microglia_phago_s6B-1_IFN_unsorted" ,
                                                "microglia_phago_s6C-1_LPS_unsorted"     ) # LPS
  ) %>%
  dplyr::mutate(sample = factor(sample,
                                levels = c("P6_iPSC_WGS"   ,
                                           "P6_D50_PreMac_240522" ,
                                           "microglia_phago_s6A-1_untreated_unsorted",
                                           "microglia_phago_s6B-1_IFN_unsorted" ,
                                           "microglia_phago_s6C-1_LPS_unsorted"  ), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle( "WGS phagocytosis")  +
  theme(axis.title.x = element_blank(),
        legend.position="none")

plist[[11]] = ((p2 + p1) / (p3 + patchwork::plot_spacer()))   +
  patchwork::plot_annotation(title = "Changing donor proportions for pool 6",
                             subtitle = "From iPSC to microglia.")

plist[[11]]


####### pool 7 #########
 plist[[9]] = merged %>%
  dplyr::filter(pool == "pool7" & sample %in% c("day0_iPSC","D36_preMAC",
                                                "P7_diff2_untreated",
                                                "P7_diff2_IFN",
                                                "P7_diff2_LPS")) %>%
  dplyr::mutate(sample = factor(sample,
                                levels = c("day0_iPSC","D36_preMAC",
                                           "P7_diff2_untreated",
                                           "P7_diff2_IFN",
                                           "P7_diff2_LPS"), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 7",
                            subtitle = "scRNA-seq samples")



plist[[9]]

plist[[10]] = merged %>%
  dplyr::filter(pool == "pool7" & sample %in% c("day0_iPSC","D36_preMAC","D47_preMAC",
                                                "Mig17_Untr", "Mig17_LPS")) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("day0_iPSC","D36_preMAC","D47_preMAC",
                                          "Mig17_Untr", "Mig17_LPS"), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 7",
                            subtitle = "WGS migration samples")



plist[[10]]

plist[[11]] = merged %>%
  dplyr::filter(pool == "pool7" & sample %in% c("day0_iPSC","D36_preMAC","D47_preMAC","D54_preMAC",
                                                "phago_s7A-1", "phago_s7B-1", "phago_s7C-1")) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("day0_iPSC","D36_preMAC","D47_preMAC","D54_preMAC",
                                          "phago_s7A-1", "phago_s7B-1", "phago_s7C-1"), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 7",
                            subtitle = "WGS phagocytosis samples")



plist[[11]]



############### fix the following later ############

plist[[12]] = merged %>%
  dplyr::filter(pool == "pool8" & sample %in% c("day0_iPSC","D36_preMAC",
                                                "P8_diff2_untreated",
                                                "P8_diff2_IFN",
                                                "P8_diff2_LPS")) %>%
  dplyr::mutate(sample = factor(sample,
                                levels = c("day0_iPSC","D36_preMAC",
                                           "P8_diff2_untreated",
                                           "P8_diff2_IFN",
                                           "P8_diff2_LPS"), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 8",
                            subtitle = "scRNA-seq samples")



plist[[12]]

plist[[13]] = merged %>%
  dplyr::filter(pool == "pool8" & sample %in% c("day0_iPSC","D36_preMAC","D40_preMAC",
                                                "Mig16_Untr", "Mig16_LPS")) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("day0_iPSC","D36_preMAC","D40_preMAC",
                                          "Mig16_Untr", "Mig16_LPS"), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 8",
                            subtitle = "WGS migration samples")



plist[[13]]

plist[[14]] = merged %>%
  dplyr::filter(pool == "pool8" & sample %in% c("day0_iPSC","D36_preMAC","D40_preMAC","D54_preMAC",
                                                "phago_s8A-1", "phago_s8B-1", "phago_s8C-1")) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("day0_iPSC","D36_preMAC","D40_preMAC","D54_preMAC",
                                          "phago_s8A-1", "phago_s8B-1", "phago_s8C-1"), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 8",
                            subtitle = "WGS phagocytosis samples")



plist[[14]]


plist[[15]] = merged %>%
  dplyr::filter(pool == "pool9" & sample %in% c("day0_iPSC","P3_d35_preMAC_WGS",
                                                "P9_Diff2_Untreated_A", "P9_Diff2_Untreated_B",
                                                "P9_Diff2_Untreated_C", "P9_Diff2_Untreated_D",
                                                "P9_Diff2_IFN_A"  ,  "P9_Diff2_IFN_B",
                                                "P9_Diff2_IFN_C", "P9_Diff2_LPS_A" , "P9_Diff2_LPS_B",  "P9_Diff2_LPS_C")) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("day0_iPSC","P3_d35_preMAC_WGS",
                                          "P9_Diff2_Untreated_A", "P9_Diff2_Untreated_B",
                                          "P9_Diff2_Untreated_C", "P9_Diff2_Untreated_D",
                                          "P9_Diff2_IFN_A"  ,  "P9_Diff2_IFN_B",
                                          "P9_Diff2_IFN_C", "P9_Diff2_LPS_A" , "P9_Diff2_LPS_B",  "P9_Diff2_LPS_C"), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 9",
                            subtitle = "scRNA-seq samples ")



plist[[15]]

plist[[16]] = merged %>%
  dplyr::filter(pool == "pool9" & sample %in% c("day0_iPSC","P3_d35_preMAC_WGS","D39_PreMac",
                                               
                                                "Untr_cells","LPS_cells")) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("day0_iPSC","P3_d35_preMAC_WGS","D39_PreMac",
                                       
                                          "Untr_cells","LPS_cells"), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 9",
                            subtitle = "WGS migration samples")



plist[[16]]

plist[[17]] = merged %>%
  dplyr::filter(pool == "pool9" & sample %in% c("day0_iPSC","P3_d35_preMAC_WGS","D39_PreMac","D50_PreMac",
                                                
                                                "phago_s9A-1",
                                                "phago_s9C-1") # LPS
                ) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("day0_iPSC","P3_d35_preMAC_WGS","D39_PreMac","D50_PreMac",
                                          
                                          "phago_s9A-1",
                                          "phago_s9C-1"), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 9",
                            subtitle = "WGS phagocytosis samples")



plist[[17]]

plist[[18]] = merged %>%
  dplyr::filter(pool == "pool10" & sample %in% c("day0_iPSC","D36_PreMac",
                                                "P10_Diff2_Untreated_A", "P10_Diff2_Untreated_B",
                                                "P10_Diff2_Untreated_C", "P10_Diff2_Untreated_D",
                                                "P10_Diff2_IFN_A"  ,  "P10_Diff2_IFN_B",
                                                "P10_Diff2_IFN_C", 
                                                "P10_Diff2_LPS_A" , "P9_Diff2_LPS_B",  "P9_Diff2_LPS_C")) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("day0_iPSC","D36_PreMac",
                                          "P10_Diff2_Untreated_A", "P10_Diff2_Untreated_B",
                                          "P10_Diff2_Untreated_C", "P10_Diff2_Untreated_D",
                                          "P10_Diff2_IFN_A"  ,  "P10_Diff2_IFN_B",
                                          "P10_Diff2_IFN_C", 
                                          "P10_Diff2_LPS_A" , "P9_Diff2_LPS_B",  "P9_Diff2_LPS_C"), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 10",
                            subtitle = "scRNA-seq samples - preMACs day 36 ")



plist[[18]]

plist[[19]] = merged %>%
  dplyr::filter(pool == "pool10" & sample %in% c("day0_iPSC","D36_PreMac","D43_PreMac","D50_PreMac",
                                                 "D54_PreMac",
                                                 "P10_Diff5_Untreated_A", "P10_Diff5_Untreated_B",
                                                 "P10_Diff5_IFN",  "P10_Diff5_LPS" )) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("day0_iPSC","D36_PreMac","D43_PreMac","D50_PreMac",
                                          "D54_PreMac",
                                          "P10_Diff5_Untreated_A", "P10_Diff5_Untreated_B",
                                          "P10_Diff5_IFN",  "P10_Diff5_LPS" ), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 10",
                            subtitle = "scRNA-seq samples - preMACs day 54 ")



plist[[19]]

plist[[20]] = merged %>%
  dplyr::filter(pool == "pool10" & sample %in% c("day0_iPSC","D36_PreMac","D43_PreMac",
                                                 "Untr_cell","LPS_cells")) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("day0_iPSC","D36_PreMac","D43_PreMac",
                                          "Untr_cell","LPS_cells" ), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 10",
                            subtitle = "WGS migration ")



plist[[20]]

plist[[21]] = merged %>%
  dplyr::filter(pool == "pool10" & sample %in% c("day0_iPSC","D36_PreMac","D43_PreMac","D50_PreMac",
                                                 "phago_s10A-1","phago_s10B-1","phago_s10C-1")) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("day0_iPSC","D36_PreMac","D43_PreMac","D50_PreMac",
                                          "phago_s10A-1","phago_s10B-1","phago_s10C-1" ), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 10",
                            subtitle = "WGS phagocytosis ")



plist[[21]]


plist[[22]] = merged %>%
  dplyr::filter(pool == "pool11" & sample %in% c("day0_iPSC","day35_preMAC",
                                                 
                                                 "P11_Diff2_Untreated_A", "P11_Diff2_Untreated_B",
                                                 "P11_Diff2_IFN",  "P11_Diff2_LPS" )) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("day0_iPSC","day35_preMAC",
                                        
                                          "P11_Diff2_Untreated_A", "P11_Diff2_Untreated_B",
                                          "P11_Diff2_IFN",  "P11_Diff2_LPS" ), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 11",
                            subtitle = "scRNA-seq samples ")



plist[[22]]

plist[[23]] = merged %>%
  dplyr::filter(pool == "pool11" & sample %in% c("day0_iPSC","day35_preMAC","day39_preMAC",
                                                 
                                                 "pool11_UNTR","pool11_LPS")) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("day0_iPSC","day35_preMAC","day39_preMAC",
                                          
                                          "pool11_UNTR","pool11_LPS"), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 11",
                            subtitle = "WGS migration ")



plist[[23]]

plist[[24]] = merged %>%
  dplyr::filter(pool == "pool11" & sample %in% c("day0_iPSC","day35_preMAC","day39_preMAC",
                                                 "day49_preMAC",
                                                 
                                                 "pool11_A1","pool11_B1","pool11_C1")) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("day0_iPSC","day35_preMAC","day39_preMAC",
                                          "day49_preMAC",
                                          
                                          "pool11_A1","pool11_B1","pool11_C1"), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 11",
                            subtitle = "WGS phagocytosis ")



plist[[24]]

plist[[25]] = merged %>%
  dplyr::filter(pool == "pool13" & sample %in% c("day0_iPSC","day36_preMAC",
                                                 
                                                 "P13_Diff2_Untreated_A", "P13_Diff2_Untreated_B",
                                                 "P13_Diff2_IFN",  "P13_Diff2_LPS" )) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("day0_iPSC","day36_preMAC",
                                          
                                          "P13_Diff2_Untreated_A", "P13_Diff2_Untreated_B",
                                          "P13_Diff2_IFN",  "P13_Diff2_LPS" ), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 13",
                            subtitle = "scRNA-seq samples ")



plist[[25]]

plist[[26]] = merged %>%
  dplyr::filter(pool == "pool13" & sample %in% c("day0_iPSC","day36_preMAC","day43_preMAC",
                                                 
                                                 "pool13_UNTR","pool13_IFN","pool13_LPS")) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("day0_iPSC","day36_preMAC","day43_preMAC",
                                          
                                          "pool13_UNTR","pool13_IFN","pool13_LPS"), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 13",
                            subtitle = "WGS migration ")



plist[[26]]

plist[[27]] = merged %>%
  dplyr::filter(pool == "pool13" & sample %in% c("day0_iPSC","day36_preMAC","day43_preMAC","day50_preMAC",
                                                 
                                                 "pool13_A1","pool13_B1","pool13_C1")) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("day0_iPSC","day36_preMAC","day43_preMAC","day50_preMAC",
                                          
                                          "pool13_A1","pool13_B1","pool13_C1"), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 13",
                            subtitle = "WGS phagocytosis")



plist[[27]]

plist[[28]] = merged %>%
  dplyr::filter(pool == "pool14" & sample %in% c("day0_iPSC","D36_PreMac",
                                                 
                                                 "P14_Diff2_Untreated_A", "P14_Diff2_Untreated_B",
                                                 "P14_Diff2_IFN_A",   "P14_Diff2_IFN_B",
                                                 "P14_Diff2_LPS_A","P14_Diff2_LPS_B" )) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("day0_iPSC","D36_PreMac",
                                          
                                          "P14_Diff2_Untreated_A", "P14_Diff2_Untreated_B",
                                          "P14_Diff2_IFN_A",   "P14_Diff2_IFN_B",
                                          "P14_Diff2_LPS_A","P14_Diff2_LPS_B"  ), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 14",
                            subtitle = "scRNA-seq samples ")



plist[[28]]

plist[[29]] = merged %>%
  dplyr::filter(pool == "pool14" & sample %in% c("day0_iPSC","D36_PreMac","D39_PreMac",
                                                 
                                                 "Mig24_Untr","Mig24_LPS")) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("day0_iPSC","D36_PreMac","D39_PreMac",
                                          
                                          "Mig24_Untr","Mig24_LPS"), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 14",
                            subtitle = "WGS migration")



plist[[29]]

plist[[30]] = merged %>%
  dplyr::filter(pool == "pool15" & sample %in% c("day0_iPSC","D36_PreMac",
                                                 
                                                 "P15_Diff2_Untreated_A", "P15_Diff2_Untreated_B",
                                                 "P15_Diff2_IFN_A",   "P15_Diff2_IFN_B",
                                                 "P15_Diff2_LPS_A","P15_Diff2_LPS_B" )) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("day0_iPSC","D36_PreMac",
                                          
                                          "P15_Diff2_Untreated_A", "P15_Diff2_Untreated_B",
                                          "P15_Diff2_IFN_A",   "P15_Diff2_IFN_B",
                                          "P15_Diff2_LPS_A","P15_Diff2_LPS_B" ), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 15",
                            subtitle = "scRNA-seq samples ")



plist[[30]]

plist[[31]] = merged %>%
  dplyr::filter(pool == "pool15" & sample %in% c("day0_iPSC","D36_PreMac","D39_PreMac",
                                                 
                                                 "Mig23_Untr", "Mig23_IFN",
                                                 "Mig23_LPS")) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("day0_iPSC","D36_PreMac","D39_PreMac",
                                          
                                          "Mig23_Untr", "Mig23_IFN",
                                          "Mig23_LPS" ), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 15",
                            subtitle = "WGS migration")



plist[[31]]


pdf(paste0(directory,"alluvial_plots_pools2-11_13-15.pdf"), 
    width = 12, height = 9)

for(i in 1:27){
  plot(plist[[i]])
}
dev.off()


# when calculating scores for differentiation efficiency preMac vs microglia, 
# do the average of the WGS and scRNA-seq sample props?
# better to correct for this as sc has more uncertainty in cell proportions
