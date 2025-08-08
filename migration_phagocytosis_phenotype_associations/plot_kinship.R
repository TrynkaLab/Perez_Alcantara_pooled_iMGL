# plotting kinship
# and genotype PCs
.libPaths(c("/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/",
            "/software/teamtrynka/conda/otar2065/lib/R/library",
            "/software/teamtrynka/ma23/R4.1/libs"))
library(patchwork)
library(tidyverse)
library(pheatmap)
library(ggrepel)
outDir = "../../data/results/genotype_kinship_inspection/"
dir.create(outDir)

kinship = read.table("../../../OTAR2065_sc_eQTL/data/kinship/all_pools.genotype.MAF01.hg38.king")
id = read.table("../../../OTAR2065_sc_eQTL/data/kinship/all_pools.genotype.MAF01.hg38.king.id")
colnames(kinship) = id$V1
rownames(kinship) = id$V1
annotation=data.frame(geno_or_wgs=c(rep("genotyped",nrow(kinship)-2),rep("WGS",2)))
rownames(annotation) = colnames(kinship)

pdf(file =  paste0(outDir,"/kinship.pdf"),width = 30, height = 30)
p = pheatmap::pheatmap(kinship,
                       annotation_col = annotation,
                       cutree_cols = 2,
                       annotation_legend=FALSE)
p
dev.off()

# plot only last 30 lines
kinship2=kinship[c((nrow(kinship)-29):nrow(kinship)),
                 c((nrow(kinship)-29):nrow(kinship))]
annotation2 = data.frame(geno_or_wgs=annotation[c((nrow(kinship)-29):nrow(kinship)),])
rownames(annotation2)=rownames(kinship2)
pdf(file =  paste0(outDir,"/kinship_30lines.pdf"),width = 12, height = 10)
p = pheatmap::pheatmap(kinship2,
                       annotation_col = annotation2,
                       cutree_cols = 2,
                       annotation_legend=FALSE)
p
dev.off()

kinship_long = read.table("../../../OTAR2065_sc_eQTL/data/kinship/all_pools.genotype.MAF01.hg38.kin0")
kinship_long %>%
  dplyr::filter(V1 == "curn_3" & V2 == "iukl_1")

# seem too highly related, because they are a different genotyping (WGS) batch from the rest

# PCA

eigen =read.table("../../../OTAR2065_sc_eQTL/data/genotype/plink_genotypes/all_pools.genotype.MAF05.eigenvec")
values =read.table("../../../OTAR2065_sc_eQTL/data/genotype/plink_genotypes/all_pools.genotype.MAF05.eigenval")
eigen = eigen[,-1]
colnames(eigen) = c("line",paste0("PC",1:60))

pdf(file =  paste0(outDir,"/genotypePC1_2.pdf"),width = 15, height = 15)

p = ggplot(eigen,aes(x=PC1,y=PC2,label=line)) +
  geom_point() + 
  ggrepel::geom_text_repel()+
  theme_minimal() + 
  xlab(paste0("PC1: ",round((values$V1[1]/sum(values$V1))*100,2),"% variance")) +
  ylab(paste0("PC1: ",round((values$V1[2]/sum(values$V1))*100,2),"% variance")) + 
  ggtitle("PCA plot of genotypes", subtitle = "no clones, LD-pruned")
p
dev.off()

pdf(file =  paste0(outDir,"/genotypePC2_3.pdf"),width = 15, height = 15)

p = ggplot(eigen,aes(x=PC3,y=PC2,label=line)) +
  geom_point() + 
  ggrepel::geom_text_repel()+
  theme_minimal() +
  xlab(paste0("PC3: ",round((values$V1[3]/sum(values$V1))*100,2),"% variance")) +
  ylab(paste0("PC2: ",round((values$V1[2]/sum(values$V1))*100,2),"% variance")) + 
  ggtitle("PCA plot of genotypes", subtitle = "no clones, LD-pruned")
p
dev.off()

pdf(file =  paste0(outDir,"/genotypePC3_4.pdf"),width = 15, height = 15)

p = ggplot(eigen,aes(x=PC3,y=PC4,label=line)) +
  geom_point() + 
  ggrepel::geom_text_repel()+
  theme_minimal() +
  xlab(paste0("PC3: ",round((values$V1[3]/sum(values$V1))*100,2),"% variance")) +
  ylab(paste0("PC4: ",round((values$V1[4]/sum(values$V1))*100,2),"% variance")) + 
  ggtitle("PCA plot of genotypes", subtitle = "no clones, LD-pruned")
p
dev.off()
