# Inspect genotype PCs to see if HipSci and IPMAR lines map where they should
library(tidyverse)
library(ggrepel)

pcs = read_tsv("../../data/genotype/plink_genotypes/all_pools.1kgenomes.genotype.MAF05.eigenvec")
hipsci_ipmar = read_tsv("../../data/genotype/all_hipsci_ipmar_donors_donorNames.txt",
                        col_names = FALSE)
thousand_info = read_csv("/lustre/scratch123/hgi/teams/trynka/resources/1000g/1000G/metadata/1000G_sample_information.csv") %>%
  dplyr::rename(IID = Sample)
pcs = pcs %>%
  dplyr::filter(PC2 < -0.0005 & PC1 < 0.0005) %>% # filter out far away ancestries
  dplyr::left_join(thousand_info,by = "IID") %>%
  dplyr::mutate(Sample = case_when(IID %in% hipsci_ipmar$X1 ~ "HipSci / IPMAR",
                                .default = Population),
                labs = case_when(IID %in% c(
                                            "boqx_2","garx_2","sojd_3","yoch_6") ~ IID,
                                 .default = NA),
                alph = case_when(IID %in% hipsci_ipmar$X1 ~ 1,
                                 .default = 0.7)) %>%
  dplyr::filter(!Sample %in% c("PUR","MXL","PEL","ITU","STU","PJL","BEB")) %>% # filter out far away ancestries
  dplyr::mutate(Sample = factor(Sample)) %>%
  dplyr::filter(!is.na(Sample))

pcs$Sample = relevel(pcs$Sample,ref = "HipSci / IPMAR")

pdf("../../data/genotype/all_pools.1kgenomes.genotype.MAF05.PCA.pdf",
    width = 5, height = 4)
pcs %>%
ggplot(aes(x = PC1, y = PC2, col = Sample, label = labs, alpha = alph )) + 
geom_point() +
  geom_text_repel(col = "black") +
   scale_color_brewer(type = "qual",direction = -1) +
  theme_minimal() +
  scale_alpha_continuous(range = c(0.4, 1)) + # doesn't seem to work from inside the tibble
  guides(alpha = "none") # remove legend for alpha

dev.off()