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
directory = "../../../../OTAR2065_WGS_iPSC_premac_micro/data/results/4.check_donor_proportions/"
dir.create(directory, recursive = T)

#### functions ####

load_sc_donors = function(doublets_unassigned_in_proportion = TRUE,
                          pools = c("pool2", "pool3", "pool4", "pool5", "pool6","pool7","pool8","pool9","pool10","pool11","pool13") ,
                          directory =  "../../../../OTAR2065_scRNA_data_processing/data/"){
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
      
      sc_ncells[[paste0(pool,".",file)]]  = read.table(paste0("../../../../OTAR2065_scRNA_data_processing/data/cellranger/",
                                                              pool,"/",file,"/vireoOutput/summary.tsv"),
                                                       header = TRUE)
      colnames( sc_ncells[[paste0(pool,".",file)]]) = c("Line","Frequency")
      
      
      # fill in missing donors
      n_lines = lines_in_pools[lines_in_pools$Pool == pool,"N_lines"]
      lines = unlist(strsplit(lines_in_pools[lines_in_pools$Pool == pool,"Lines"],split = ";"))
      sc_ncells[[paste0(pool, ".", file)]] =  sc_ncells[[paste0(pool, ".", file)]] %>%
        dplyr::rows_insert(tibble(Line = lines[!lines %in% sc_ncells[[paste0(pool, ".", file)]]$Line])) %>%
        tidyr::replace_na(list(Frequency = 0))
      
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
          dplyr::add_row(Line = "unassigned_plus_doublets",
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
    labs(title = paste0("Changing donor proportions for ",unique(long_df$pool)))
  plot(p)
}

#########
# read in info on donors in pools
# from vireo
microglia_prop = load_sc_donors(doublets_unassigned_in_proportion = FALSE,
                                pools = c("pool2","pool3","pool4",
                                          "pool5","pool6","pool7","pool8","pool11","pool13"))

microglia_prop = microglia_prop %>%
  dplyr::filter(pool %in% c("pool2","pool3","pool4",
                            "pool5","pool6","pool7","pool8","pool11","pool13")) %>%
  dplyr::filter(str_detect(set, "plst|untreated|Untreated|AG|YC|ITMG")) %>%
  dplyr::group_by(pool,Line) %>%
  dplyr::summarise(prop = mean(final_prop),sample="microglia") %>%
  dplyr::ungroup() %>%
  dplyr::relocate(Line,prop,pool,sample)
  

# iPSC and precursors prop
files = list.files("../../../../OTAR2065_WGS_iPSC_premac_micro/data/w/")
wgs_prop=list()
for(fil in files){q
  name = paste(str_split(fil,"_")[[1]][1],
               str_split(fil,"_")[[1]][2],
               str_split(fil,"_")[[1]][3],sep = "_")
  wgs_prop[[name]] = read.table(paste0("../../../../OTAR2065_WGS_iPSC_premac_micro/data/w/",fil))
  colnames(wgs_prop[[name]]) = c("Line","prop")
  wgs_prop[[name]]$pool = str_split(fil,"_")[[1]][1]
  wgs_prop[[name]]$sample = paste(str_split(fil,"_")[[1]][2],str_split(fil,"_")[[1]][3],sep = "_")
}

wgs_prop = do.call("rbind",wgs_prop)
wgs_prop = wgs_prop %>%
  mutate(sample = if_else(sample %in% c("iPSC_pool","iPSC_WGS","iPSC_D0","D0_iPSC"),"day0_iPSC",sample)) %>%
  mutate(prop = if_else(prop < 0,0,prop)) %>%
  mutate(pool = case_when(pool== "P3" ~ "pool3",
                          pool== "P4" ~ "pool4",
                          pool== "P5" ~ "pool5",
                          pool== "P6" ~ "pool6",
                          .default = as.character(pool)))
merged = rbind(wgs_prop,microglia_prop)
rm(wgs_prop,microglia_prop)

p1 = merged %>%
  dplyr::filter(pool == "pool2") %>%
  dplyr::mutate(sample = factor(sample,
                                levels = c("day0_iPSC","D7-post-thaw_iPSC",
                                           "day29_preMAC","day57_preMAC","day109_preMAC",
                                           "microglia"), 
                                ordered = TRUE)) %>% 
  alluvial_plot()

png(paste0(directory,"alluvial_plot_pool2_iPSC_preMac_microglia.png"), 
    width = 9, height = 6, units = "in", res = 400)

p1
dev.off()

p1 = merged %>%
  dplyr::filter(pool == "pool3") %>%
  dplyr::mutate(sample = factor(sample,
                                levels = c("day0_iPSC","d35_preMAC","d57_preMAC","microglia"), 
                                ordered = TRUE)) %>% 
  alluvial_plot()

png(paste0(directory,"alluvial_plot_pool3_iPSC_preMac_microglia.png"), 
    width = 8, height = 6, units = "in", res = 400)

p1
dev.off()

p1 = merged %>%
  dplyr::filter(pool == "pool4") %>%
  dplyr::mutate(sample = factor(sample,
                                levels = c("day0_iPSC","d35_preMAC","d57_preMAC","microglia"), 
                                ordered = TRUE)) %>% 
  alluvial_plot()

png(paste0(directory,"alluvial_plot_pool4_iPSC_preMac_microglia.png"), 
    width = 8, height = 6, units = "in", res = 400)

p1
dev.off()
p1 = merged %>%
  dplyr::filter(pool == "pool5") %>%
  dplyr::mutate(sample = factor(sample,
                                levels = c("day0_iPSC","d35_preMAC","microglia"), 
                                ordered = TRUE)) %>% 
  alluvial_plot()

png(paste0(directory,"alluvial_plot_pool5_iPSC_preMac_microglia.png"), 
    width = 8, height = 6, units = "in", res = 400)

p1
dev.off()
p1 = merged %>%
  dplyr::filter(pool == "pool6") %>%
  dplyr::mutate(sample = factor(sample,
                                levels = c("day0_iPSC","d35_preMAC","microglia"), 
                                ordered = TRUE)) %>% 
  alluvial_plot()

png(paste0(directory,"alluvial_plot_pool6_iPSC_preMac_microglia.png"), 
    width = 8, height = 6, units = "in", res = 400)

p1
dev.off()
p1 = merged %>%
  dplyr::filter(pool == "pool7") %>%
  dplyr::mutate(sample = factor(sample,
                                levels = c("day0_iPSC","D36_preMAC","D47_preMAC","D54_preMAC","microglia"), 
                                ordered = TRUE)) %>% 
  alluvial_plot()

png(paste0(directory,"alluvial_plot_pool7_iPSC_preMac_microglia.png"), 
    width = 8, height = 6, units = "in", res = 400)

p1
dev.off()

p1 = merged %>%
  dplyr::filter(pool == "pool8") %>%
  dplyr::mutate(sample = factor(sample,
                                levels = c("day0_iPSC","D36_preMAC","D40_preMAC","D54_preMAC","microglia"), 
                                ordered = TRUE)) %>% 
  alluvial_plot()

png(paste0(directory,"alluvial_plot_pool8_iPSC_preMac_microglia.png"), 
    width = 8, height = 6, units = "in", res = 400)

p1
dev.off()


p1 = merged %>%
  dplyr::filter(pool == "pool11") %>%
  dplyr::mutate(sample = factor(sample,
                                levels = c("day0_iPSC","day35_preMAC","day39_preMAC","day49_preMAC","microglia"), 
                                ordered = TRUE)) %>% 
  alluvial_plot()

png(paste0(directory,"alluvial_plot_pool11_iPSC_preMac_microglia.png"), 
    width = 8, height = 6, units = "in", res = 400)

p1
dev.off()

p2 = merged %>%
  dplyr::filter(pool == "pool13") %>%
  dplyr::mutate(sample = factor(sample,
                                levels = c("day0_iPSC","day35_preMAC","day43_preMAC","day50_preMAC","microglia"), 
                                ordered = TRUE)) %>% 
  alluvial_plot()

png(paste0(directory,"alluvial_plot_pool13_iPSC_preMac_microglia.png"), 
    width = 8, height = 6, units = "in", res = 400)

p2
dev.off()

merged %>%
  dplyr::filter(pool %in% c(paste0("pool",2:8),"pool11","pool13")) %>%
  
write.table(.,paste0(directory,"pools2-8_11_13_changing_props_iPSC_preMacs_microglia.txt"),
            row.names = FALSE,col.names = TRUE,quote = FALSE,sep = "\t")
