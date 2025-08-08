# Inspecting metadata to find how many lines are kept at the tensorQTL analysis
library(tidyverse)

outDir = "../../data/results/4.Inspect_eQTL_results"

metadata = read_rds(paste0(outDir,"/metadata_list.rds"))
metadata = metadata %>%
  imap_dfr(~mutate(.x, group = .y)) %>%
  mutate(line = donor_id) %>%
  dplyr::filter(group %in% c("60_Not_proliferating_IFN" ,                     
                          "60_Not_proliferating_LPS"   ,                    "60_Not_proliferating_untreated"  ,               "60_Proliferating_LPS"  ,                        
                       "60_Proliferating_untreated"   )) %>%
  select(line,ncells,prolif_treatment,pool)
  
metadata %>%
  pivot_wider(names_from = prolif_treatment,values_from = ncells) %>%
  write_csv(paste0(outDir,"/lines_retained_tensorQTL.csv"))

length(unique(metadata$line)) # 170 lines retained for all conditions, pools 2-10


p = metadata %>%
  ggplot(aes(x=reorder(line,-ncells,FUN = median),y=ncells)) + 
  geom_boxplot() +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),breaks = c(10,50,100,1000,5000,10000,25000,50000),
                     limits = c(10,50000)) + 
  ylab("Number of cells (pseudolog10 axis)") +
  xlab("Line") +
  theme_bw()+
  ggtitle("Number of cells in all proliferation x treatment combinations")+
  coord_flip()

png(paste0(outDir,"/boxplot_number_of_cells_retained_tensorQTL.png"),
    width = 6,height = 17,units = "in",res = 400)

plot(p)

dev.off()
