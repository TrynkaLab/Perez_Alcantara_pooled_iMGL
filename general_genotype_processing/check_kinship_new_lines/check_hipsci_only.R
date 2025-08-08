# check hipsci only
.libPaths(c('/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/',"/software/teamtrynka/ma23/R4.1/libs",
            "/software/R-4.1.0/lib/R/library"))
library(tidyverse)
library(ggrepel)
outDir = "../../data/check_kinship_new_lines/"
genotype=read.table("../../data/check_kinship_new_lines/genotype/plink_genotypes/tocheck.genotype.MAF05.eigenvec")
colnames(genotype) = c("donor","line",paste0("genotypePC",1:60))

p2 = ggplot(genotype,aes(x=genotypePC1,y=genotypePC2, label=line))+
  geom_point(alpha=0.7) + 
  geom_text_repel(max.overlaps = 100, col="black",fontface="bold")+
  theme_bw() + 
  ggtitle("Genotype PCs")
pdf(file =  paste0(outDir,"/PC1_PC2_genotype_hipsci.pdf"),width = 10, height = 10)
plot(p2)
dev.off()
# the clones look very far away - need to repeat with clones removed