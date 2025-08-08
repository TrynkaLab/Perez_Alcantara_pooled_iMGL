# PCA plots 1000 genomes
.libPaths(c('/software/teamtrynka/common/R/x86_64-pc-linux-gnu/4.1/',"/software/teamtrynka/ma23/R4.1/libs",
            "/software/R-4.1.0/lib/R/library"))
library(tidyverse)
library(ggrepel)
outDir = "../../data/check_kinship_new_lines/"
# explore genotype PCs with 1k genomes
eur_metadata=read_csv("../../../resources/1000Genomes_eur_ancestry_metadata.csv") %>%
  dplyr::rename(donor=id)

genotype_1k_hipsci=read.table("../../data/check_kinship_new_lines/genotype/plink_genotypes/tocheck.1kgenomes.genotype.MAF05.eigenvec")
names = ifelse(genotype_1k_hipsci$V1==genotype_1k_hipsci$V2,
               yes = genotype_1k_hipsci$V1, 
               no = paste(genotype_1k_hipsci$V1,genotype_1k_hipsci$V2,sep = "_"))
names = ifelse(grepl("NA|HG",names),
               yes =unlist(lapply(str_split(names, "_"), function(x) x[2])),
               no = names)
genotype_1k_hipsci = genotype_1k_hipsci[c(-1:-2)]
colnames(genotype_1k_hipsci) = paste0("genotypePC",1:ncol(genotype_1k_hipsci))
genotype_1k_hipsci$donor = names
genotype_1k_hipsci$origin = ifelse(grepl("NA|HG",genotype_1k_hipsci$donor),
                                   yes = "1000 Genomes", no = "hiPSCi")
genotype_1k_hipsci = merge(genotype_1k_hipsci,eur_metadata,all.x=TRUE, by="donor") %>%
  mutate(pop = if_else(is.na(pop),true = "hiPSCi",false=pop))
genotype_1k_hipsci = genotype_1k_hipsci %>%

  mutate(donor = case_when(origin == "hiPSCi" ~ str_split_fixed(donor,pattern = "_",n=2)[,2],
                           .default = donor)) %>%
  mutate(hipsci_labels = if_else(pop=="hiPSCi",true = donor,false=NA),
         pop_labels = if_else(!duplicated(pop),true = pop,false = NA))

p1 = ggplot(genotype_1k_hipsci,aes(x=genotypePC1,y=genotypePC2,label=hipsci_labels,color=origin))+
  geom_point() + 
  geom_text_repel()+
  theme_bw() + 
  ggtitle("Genotype PCs with 1000 genomes")
p2 = ggplot(genotype_1k_hipsci,aes(x=genotypePC1,y=genotypePC2,label=pop_labels,color=pop,fill=pop))+
  geom_point(alpha=0.7) + 
  geom_text_repel(max.overlaps = 100, col="black",fontface="bold")+
  theme_bw() + 
  ggtitle("Genotype PCs with 1000 genomes")
pdf(file =  paste0(outDir,"/PC1_PC2_genotype_1kgenomes_hipsci.pdf"),width = 10, height = 10)
plot(p2)
dev.off()

pdf(file =  paste0(outDir,"/PC1_PC2_genotype_1kgenomes_hipsci_highlight.pdf"),width = 10, height = 10)
plot(p1 )
dev.off()

# subset 1k genomes to zoom
p1 = genotype_1k_hipsci %>%
  dplyr::filter(pop %in% c("hiPSCi","TSI","IBS","MXL","GBR","CEU","FIN","CLM","PUR","PIL")) %>%
  ggplot(aes(x=genotypePC1,y=genotypePC2,label=pop_labels,color=pop,fill=pop, shape = origin))+
  geom_point(alpha=0.7, size = 2) + 
  theme_bw() + 
  ggtitle("Genotype PCs with 1000 genomes")

p2 = genotype_1k_hipsci %>%
  dplyr::filter(pop %in% c("hiPSCi","TSI","IBS","MXL","GBR","CEU","FIN","CLM","PUR","PIL")) %>%
  dplyr::mutate(origin = factor(origin,levels = c("hiPSCi","1000 Genomes"),ordered = TRUE)) %>%
  ggplot(aes(x=genotypePC1,y=genotypePC2,label=hipsci_labels,color=origin))+
  scale_color_manual(values = c("hiPSCi" = "darkblue", "1000 Genomes" = "lightblue")) +
  geom_point() + 
  geom_text_repel()+
  theme_bw() + 
  ggtitle("Genotype PCs with 1000 genomes")

pdf(file =  paste0(outDir,"/PC1_PC2_genotype_1kgenomes_hipsci_highlight_subset.pdf"),width = 12, height = 7)
plot(p1 + p2)
dev.off()