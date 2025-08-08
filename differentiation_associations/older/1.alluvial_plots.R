# alluvial plots iPSC - premacs - microglia
# all pools
# Checking donor proportions and make alluvial plots
# with WGS and vireo estimates at differentiation stages
# all pools
.libPaths(c('/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/',"/software/teamtrynka/ma23/R4.1/libs",
            "/software/R-4.1.0/lib/R/library"))
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
directory = "../../data/results/1.alluvial_plots/"
dir.create(directory, recursive = T)

#### functions ####

load_sc_donors = function( pools = paste0("pool",2:13) ,
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
                                                              pool,"/",file,"/vireoOutput/donor_ids.tsv"),
                                                       header = TRUE)
      
      # if doublets, put doublet in donor_id
      # if unassigned, check which probability is higher
      sc_ncells[[paste0(pool,".",file)]] = sc_ncells[[paste0(pool,".",file)]] %>%
        dplyr::mutate(new_donor_id = case_when(donor_id == "doublet" ~ best_doublet,
                                               donor_id == "unassigned" & prob_max < prob_doublet ~ best_doublet,
                                               donor_id == "unassigned" & prob_max > prob_doublet ~ best_singlet,
                                               .default = donor_id)) %>%
        tidyr::separate_rows(new_donor_id, sep = ",") %>%
        distinct()  %>%
        dplyr::count(new_donor_id) %>%
        dplyr::rename(Line = new_donor_id, Frequency = n)
      
      
      
      
      lines = unlist(strsplit(lines_in_pools[lines_in_pools$Pool == pool,"Lines"],split = ";"))
      sc_ncells[[paste0(pool, ".", file)]] =  sc_ncells[[paste0(pool, ".", file)]] %>%
        dplyr::rows_insert(tibble(Line = lines[!lines %in% sc_ncells[[paste0(pool, ".", file)]]$Line])) %>%
        tidyr::replace_na(list(Frequency = 0))
      
      sc_ncells[[paste0(pool,".",file)]]  = sc_ncells[[paste0(pool,".",file)]] %>%
        dplyr::reframe(Line = Line,
                       prop = Frequency / sum(Frequency),
                       pool = pool, 
                       sample = file, 
                       stage = "microglia")
      
      
    }
  }
  
  sc_ncells = do.call(rbind,Map(cbind,set = names(sc_ncells),sc_ncells))
  return(sc_ncells)
}

alluvial_plot = function(long_df){
  p = ggplot(long_df, aes(y = prop,x = sample,
                          alluvium=Line,stratum = Line,
                          fill = Line, colour = Line,
                          label = Line)) +
    scale_x_discrete(expand = c(.1, .1)) +
    scale_fill_viridis_d(option = "D")+
    scale_color_viridis_d(option = "D")+
    geom_flow(alpha = .6) +
    geom_stratum(alpha = .0) +
    # geom_text(colour = "black", position = position_stack(vjust = 0.5),size = 3) +
    ggfittext::geom_fit_text(stat = "stratum",colour = "black", width = 1/4, min.size = 3) +
    theme_bw() + 
    theme(legend.position = "none",axis.text= element_text(face="bold"),
          axis.title = element_text(face="bold"), plot.title =element_text(face="bold") ) +
    labs(title = paste0("Changing donor proportions for ",unique(long_df$pool))) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  plot(p)
}

#########
# read in info on donors in pools
# from vireo
microglia_prop = load_sc_donors(pools = paste0("pool",2:13)) %>%
  dplyr::select(!set) %>%
  dplyr::mutate(type = "scRNA-seq")
# WGS: iPSC and precursors prop, also migration and phagocytosis assays
wgs_prop=list()

for(pool in  paste0("pool",2:13)){
  files = list.files(paste0("../../data/w/",pool))
  for(fil in files){
    name = paste(stringr::str_split(fil,"_")[[1]][1],
                 stringr::str_split(fil,"_")[[1]][2],
                 stringr::str_split(fil,"_")[[1]][3],sep = "_")
    wgs_prop[[name]] = read.table(paste0("../../data/w/",pool,"/",fil))
    colnames(wgs_prop[[name]]) = c("Line","prop")
    wgs_prop[[name]]$pool = pool
    wgs_prop[[name]]$sample = paste(stringr::str_split(fil,"_")[[1]][2],
                                    stringr::str_split(fil,"_")[[1]][3],sep = "_")
  }
}

wgs_prop = do.call("rbind",wgs_prop)
wgs_prop = wgs_prop %>%
  mutate(sample = if_else(sample %in% c("iPSC_pool","iPSC_WGS","iPSC_D0","D0_iPSC"),"day0_iPSC",sample)) %>%
  mutate(prop = if_else(prop < 0,0,prop)) %>%
  mutate(stage = case_when(str_detect(sample,"phago") ~ "microglia",
                           str_detect(sample,"iPSC")~ "iPSC",
                           str_detect(sample,"migration") ~ "microglia",
                           str_detect(sample,"PreMAC") ~ "preMac",
                           str_detect(sample,"preMAC") ~ "preMac",
                           str_detect(sample,"PreMac") ~ "preMac",
                           str_detect(sample,"UNTR") ~ "microglia",
                           str_detect(sample,"LPS")  ~ "microglia",
                           str_detect(sample,"INF")  ~ "microglia"),
         type = "WGS")
merged = rbind(wgs_prop,microglia_prop)
rm(wgs_prop,microglia_prop)

plist = list()
plist[[1]] = merged %>%
  dplyr::filter(pool == "pool2" & sample %in% c("day0_iPSC","D7-post-thaw_iPSC",
                                                "day29_preMAC","day57_preMAC","day109_preMAC")) %>%
  dplyr::mutate(sample = factor(sample,
                                levels = c("day0_iPSC","D7-post-thaw_iPSC",
                                           "day29_preMAC","day57_preMAC","day109_preMAC"), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 2",
                            subtitle = "During iPSC freeze-thawing and preMAC aging")



plist[[1]]

plist[[2]] = merged %>%
  dplyr::filter(pool == "pool2" & sample %in% c("D7-post-thaw_iPSC",
                                                "day29_preMAC",
                                                "migration_s81", # untreated
                                                "migration_s101", # IFN
                                                "chrmX_WGS",
                                                "AG_10",
                                                "ITMG")) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("D7-post-thaw_iPSC",
                                          "day29_preMAC",
                                          "chrmX_WGS",
                                          "AG_10",
                                          "ITMG",
                                          "migration_s81", # untreated
                                          "migration_s101" # IFN
                                ), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 2",
                            subtitle = "From best match of iPSC and preMACs")



plist[[2]]

plist[[3]]= merged %>%
  dplyr::filter(pool == "pool3" & sample %in% c("day0_iPSC","d35_preMAC",
                                                "h3_diff2_untreated", # untreated
                                                "h3_diff2_IFNg",  # IFN
                                                "h3_diff2_LPS", # LPS
                                                "migration_s135", # untreated
                                                "migration_s147") # LPS
  ) %>%
  dplyr::mutate(sample = factor(sample,
                                levels = c("day0_iPSC","d35_preMAC",
                                           "h3_diff2_untreated", # untreated
                                           "h3_diff2_IFNg",  # IFN
                                           "h3_diff2_LPS", # LPS
                                           "migration_s135", # untreated
                                           "migration_s147"), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 3",
                            subtitle = "From best match of preMAC for migration (day 39)")



plist[[3]]


plist[[4]]= merged %>%
  dplyr::filter(pool == "pool3" & sample %in% c("day0_iPSC","d35_preMAC","d57_preMAC",
                                                "phago_s3D-1", # untreated
                                                "phago_s3E-1",  # IFN
                                                "phago_s3F-1") # LPS
  ) %>%
  dplyr::mutate(sample = factor(sample,
                                levels = c("day0_iPSC","d35_preMAC","d57_preMAC",
                                           "phago_s3D-1","phago_s3E-1","phago_s3F-1",
                                           "migration_s135",
                                           "migration_s147"), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 3",
                            subtitle = "Age of preMAC for phagocytosis samples: day 57")



plist[[4]]


plist[[5]] = merged %>%
  dplyr::filter(pool == "pool4" & sample %in% c("day0_iPSC","d35_preMAC",
                                                "h4_diff2_untreated", # untreated
                                                "h4_diff2_IFNg",  # IFN
                                                "h4_diff2_LPS", # LPS
                                                "migration_s159", # untreated
                                                "migration_s171") # LPS
  ) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("day0_iPSC","d35_preMAC",
                                          "h4_diff2_untreated", # untreated
                                          "h4_diff2_IFNg",  # IFN
                                          "h4_diff2_LPS", # LPS
                                          "migration_s159", # untreated
                                          "migration_s171") # LPS
                                , 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 4",
                            subtitle = "From best match of preMAC for migration (day 42)")




plist[[5]]

plist[[6]] = merged %>%
  dplyr::filter(pool == "pool4" & sample %in% c("day0_iPSC","d35_preMAC","d57_preMAC",
                                                "phago_s4D-1","phago_s4E-1","phago_s4F-1")) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("day0_iPSC","d35_preMAC","d57_preMAC",
                                          "phago_s4D-1","phago_s4E-1","phago_s4F-1"), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 4",
                            subtitle = "Age of preMAC for phagocytosis samples: day 57")



plist[[6]]

plist[[7]] = merged %>%
  dplyr::filter(pool == "pool5" & sample %in% c("day0_iPSC","d35_preMAC",
                                                "P5_diff2_untreated",
                                                "P5_diff2_IFN",
                                                "P5_diff2_LPS",
                                                "phago_s5A-1","phago_s5B-1","phago_s5C-1")) %>%
  dplyr::mutate(sample = factor(sample,
                                levels = c("day0_iPSC","d35_preMAC",
                                           "P5_diff2_untreated",
                                           "P5_diff2_IFN",
                                           "P5_diff2_LPS",
                                           "phago_s5A-1","phago_s5B-1","phago_s5C-1"), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 5",
                            subtitle = "From best match of preMAC for phagocytosis (day 50)")



plist[[7]]

plist[[8]] = merged %>%
  dplyr::filter(pool == "pool6" & sample %in% c("day0_iPSC","d35_preMAC",
                                                "P6_diff2_untreated",
                                                "P6_diff2_IFN",
                                                "P6_diff2_LPS",
                                                "phago_s6A-1","phago_s6B-1","phago_s6C-1")) %>%
  dplyr::mutate(sample = factor(sample,
                                levels = c("day0_iPSC","d35_preMAC",
                                           "P6_diff2_untreated",
                                           "P6_diff2_IFN",
                                           "P6_diff2_LPS",
                                           "phago_s6A-1","phago_s6B-1","phago_s6C-1"), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 6",
                            subtitle = "From best match of preMAC for phagocytosis (day 50)")



plist[[8]]

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
                                levels =c("day0_iPSC","D36_preMAC","D47_preMAC","D54_preMAC",
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
  dplyr::filter(pool == "pool9" & sample %in% c("day0_iPSC","D35_PreMac",
                                                "P9_Diff2_Untreated_A", "P9_Diff2_Untreated_B",
                                                "P9_Diff2_Untreated_C", "P9_Diff2_Untreated_D",
                                                "P9_Diff2_IFN_A"  ,  "P9_Diff2_IFN_B",
                                                "P9_Diff2_IFN_C", "P9_Diff2_LPS_A" , "P9_Diff2_LPS_B",  "P9_Diff2_LPS_C")) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("day0_iPSC","D35_PreMac",
                                          "P9_Diff2_Untreated_A", "P9_Diff2_Untreated_B",
                                          "P9_Diff2_Untreated_C", "P9_Diff2_Untreated_D",
                                          "P9_Diff2_IFN_A"  ,  "P9_Diff2_IFN_B",
                                          "P9_Diff2_IFN_C", "P9_Diff2_LPS_A" , "P9_Diff2_LPS_B",  "P9_Diff2_LPS_C"), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 9",
                            subtitle = "scRNA-seq samples ")



plist[[15]]

plist[[16]] = merged %>%
  dplyr::filter(pool == "pool9" & sample %in% c("day0_iPSC","D35_PreMac","D39_PreMac",
                                                
                                                "Untr_cells","LPS_cells")) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("day0_iPSC","D35_PreMac","D39_PreMac",
                                          
                                          "Untr_cells","LPS_cells"), 
                                ordered = TRUE)) %>% 
  alluvial_plot() + ggtitle(label = "Changing donor proportions for pool 9",
                            subtitle = "WGS migration samples")



plist[[16]]

plist[[17]] = merged %>%
  dplyr::filter(pool == "pool9" & sample %in% c("day0_iPSC","D35_PreMac","D39_PreMac","D50_PreMac",
                                                
                                                "phago_s9A-1",
                                                "phago_s9C-1") # IFN
  ) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("day0_iPSC","D35_PreMac","D39_PreMac","D50_PreMac",
                                          
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
  dplyr::filter(pool == "pool10" & sample %in% c("day0_iPSC","D36_PreMac","D43_PreMac",
                                                 "phago_s10A-1","phago_s10B-1","phago_s10C-1")) %>%
  dplyr::mutate(sample = factor(sample,
                                levels =c("day0_iPSC","D36_PreMac","D43_PreMac",
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



pdf(paste0(directory,"alluvial_plots_pools2-13.pdf"), 
    width = 12, height = 9)

for(i in 1:27){
  plot(plist[[i]])
}
dev.off()

merged %>%
  dplyr::filter(pool %in% c(paste0("pool",2:11),"pool13")) %>%
  
  write.table(.,paste0(directory,"pools2-11_13_changing_props_iPSC_preMacs_microglia_WGS_sc.txt"),
              row.names = FALSE,col.names = TRUE,quote = FALSE,sep = "\t")

# when calculating scores for differentiation efficiency preMac vs microglia, 
# do the average of the WGS and scRNA-seq sample props?
# better to correct for this as sc has more uncertainty in cell proportions
