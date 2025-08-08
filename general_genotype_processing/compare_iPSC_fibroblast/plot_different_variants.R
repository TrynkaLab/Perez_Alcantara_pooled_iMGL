# get numbers of variants that are different

library(tidyverse)

directory = "../../../data/results/fibro"
dir.create(directory, recursive = TRUE)

dirs = list.dirs("../../../data/",recursive = FALSE,full.names = FALSE)

diff_variants = list()
for(dir in dirs){
  diff_variants[[dir]] = read.table(paste0("../../../data/",dir,"/",dir,
                                           ".genotype.hg38.vcf.gz_Nvariants_retained.txt"), 
                                    header = TRUE)
  
}


diff_variants = do.call("rbind",Map(cbind, 
                                    Name = names(diff_variants), diff_variants))

p = diff_variants %>%
  mutate(percent_different_vars = (non_identical_genotype_length/original_length) *100,
         dummy_x = "x") %>%
  ggplot( aes(x = dummy_x, y=percent_different_vars)) +
  geom_boxplot(notch=FALSE) +
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  theme_bw() + 
  ylab("% of variants that differ between fibro and iPSC") + 
  ggtitle("% differences in 10 donors. ~40 million SNPs each") +
  theme(axis.title.x = element_blank(),axis.text.x = element_blank())
  
png(paste0(directory,"/percent_variants_different_fibro_vs_iPSC.png"),
    height = 4,width = 5,res = 400,units = "in")
p
dev.off()

diff_variants %>%
  mutate(percent_different_vars = (non_identical_genotype_length/original_length) *100,
         dummy_x = "x") %>%
  summary(percent_different_vars)
# percent_different_vars       
# Min.   :0.1708                
# 1st Qu.:0.1732         
# Median :0.1753        
# Mean   :0.1843                           
# 3rd Qu.:0.1779                           
# Max.   :0.2690 
