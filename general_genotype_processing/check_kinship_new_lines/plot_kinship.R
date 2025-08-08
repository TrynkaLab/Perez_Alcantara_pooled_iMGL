# plotting kinship
# and genotype PCs
library(patchwork)
library(tidyverse)
library(pheatmap)
library(ggrepel)
outDir = "../../../data/check_kinship_ipmar/"
dir.create(outDir)

kinship = read.table("../../../data/check_kinship_ipmar/ipmar.king")
id = read.table("../../../data/check_kinship_ipmar/ipmar.king.id")
colnames(kinship) = id$V1
rownames(kinship) = id$V1
annotation=data.frame(geno_or_wgs=c(rep("genotyped",nrow(kinship)-3),rep("WGS",3)))
rownames(annotation) = colnames(kinship)

pdf(file =  paste0(outDir,"/kinship.pdf"),width = 30, height = 30)
p = pheatmap::pheatmap(kinship,
                       # annotation_col = annotation,
                       cutree_cols = 2,
                       annotation_legend=FALSE)
p
dev.off()

# plot only last 30 lines
kinship2=kinship[c(220,(nrow(kinship)-29):nrow(kinship)), 
                 c(220,(nrow(kinship)-29):nrow(kinship))]
annotation2 = data.frame(geno_or_wgs=annotation[c(220,(nrow(kinship)-29):nrow(kinship)),])
rownames(annotation2)=rownames(kinship2)
pdf(file =  paste0(outDir,"/kinship_31lines.pdf"),width = 12, height = 10)
p = pheatmap::pheatmap(kinship2,
                       annotation_col = annotation2,
                       annotation_legend=FALSE)
p
dev.off()

kinship_long = read.table("../../../data/check_kinship_ipmar/kinship/ipmar.kin0")
kinship_long %>%
  dplyr::filter(V1 == "iukl_1" & V2 == "curn_3")

# they are no longer highly related after correct pipeline, despite different genotyping (WGS) batch from the rest

# PCA

eigen =read.table("../../../data/check_kinship_ipmar/genotype/all_pools.genotype.MAF05.eigenvec")
values =read.table("../../../data/check_kinship_ipmar/genotype/all_pools.genotype.MAF05.eigenval")
eigen = eigen[,-1]
colnames(eigen) = c("line",paste0("PC",1:60))

pdf(file =  paste0(outDir,"/genotypePC1_2.pdf"),width = 8, height = 8)

p = ggplot(eigen,aes(x=PC1,y=PC2,label=line)) +
  geom_point() + 
  ggrepel::geom_text_repel()+
  theme_minimal() + 
  xlab(paste0("PC1: ",round((values$V1[1]/sum(values$V1))*100,2),"% variance")) +
  ylab(paste0("PC1: ",round((values$V1[2]/sum(values$V1))*100,2),"% variance")) + 
  ggtitle("PCA plot of genotypes", subtitle = "no clones, LD-pruned")
p
dev.off()

pdf(file =  paste0(outDir,"/genotypePC2_3.pdf"),width = 8, height = 8)

p = ggplot(eigen,aes(x=PC3,y=PC2,label=line)) +
  geom_point() + 
  ggrepel::geom_text_repel() +
  theme_minimal() +
  xlab(paste0("PC3: ",round((values$V1[3]/sum(values$V1))*100,2),"% variance")) +
  ylab(paste0("PC2: ",round((values$V1[2]/sum(values$V1))*100,2),"% variance")) + 
  ggtitle("PCA plot of genotypes", subtitle = "no clones, LD-pruned")
p
dev.off()

pdf(file =  paste0(outDir,"/genotypePC3_4.pdf"),width = 8, height = 8)

p = ggplot(eigen,aes(x=PC3,y=PC4,label=line)) +
  geom_point() + 
  ggrepel::geom_text_repel() +
  theme_minimal() +
  xlab(paste0("PC3: ",round((values$V1[3]/sum(values$V1))*100,2),"% variance")) +
  ylab(paste0("PC4: ",round((values$V1[4]/sum(values$V1))*100,2),"% variance")) + 
  ggtitle("PCA plot of genotypes", subtitle = "no clones, LD-pruned")
p
dev.off()
