# check gtc and idat files from HipSci for Kaur
library(tidyverse)
lines_in_pools = readr::read_table("../../data/sampleMetadata.txt")

lines_sep = lines_in_pools %>%
  dplyr::filter(sample_ID %in% paste("pool",2:17,sep = "_")) %>%
  dplyr::mutate(donor_IDs = stringr::str_split(donor_IDs,pattern = ";")) %>%
  dplyr::select(donor_IDs) %>%
  unlist() %>%
  unique() %>%
  sort()

lines_no_wgs_ipmar = lines_in_pools %>%
  dplyr::filter(sample_ID %in% "all_pools_no_curn_iukl") %>%
  dplyr::mutate(donor_IDs = stringr::str_split(donor_IDs,pattern = ";")) %>%
  dplyr::select(donor_IDs) %>%
  unlist() %>%
  unique() %>%
  sort()

prs = lines_in_pools %>%
  dplyr::filter(sample_ID %in% "PRS") %>%
  dplyr::mutate(donor_IDs = stringr::str_split(donor_IDs,pattern = ";")) %>%
  dplyr::select(donor_IDs) %>%
  unlist() %>%
  unique() %>%
  sort()

length(lines_sep %in% lines_no_wgs_ipmar) # 261
table(lines_sep %in% lines_no_wgs_ipmar) # 16 missing, IPMAR, curn_3 and iukl_1 (expected)
table(lines_no_wgs_ipmar %in% prs)
setdiff(lines_sep, lines_no_wgs_ipmar)
intersect(lines_sep, lines_no_wgs_ipmar)
intersect(prs, lines_no_wgs_ipmar)

lines_in_pools %>%
  dplyr::filter(sample_ID %in% paste("pool",2:17,sep = "_")) %>%
  dplyr::mutate(donor_IDs = stringr::str_split(donor_IDs,pattern = ";")) %>%
  dplyr::select(donor_IDs) %>%
  unlist() %>%
  unique() %>%
  sort()

# check if lines are present in idat folders
# /lustre/scratch123/hgi/projects/hipsci/releases/data/gtarray/primary_data
files = list.files(path = "/lustre/scratch123/hgi/projects/hipsci/releases/data/gtarray/primary_data",
                   recursive = TRUE,full.names = TRUE)

matched_list = list()
for(line in lines_no_wgs_ipmar){
  matched_list[[line]] =  str_subset(string = files, pattern = fixed(line))[3]
}

hipsci_genotype = do.call("rbind",matched_list)
nrow(hipsci_genotype) == length(lines_no_wgs_ipmar) # all there

matched_list = list()
for(line in lines_no_wgs_ipmar){
  matched_list[[line]] =  str_subset(string = files, pattern = fixed(line))
}
hipsci_genotype = do.call(rbind, lapply(matched_list, as.data.frame))

nrow(hipsci_genotype)/3 == 245 # 3 files per donor

write.table(hipsci_genotype,"../../data/hipsci_genotype_files_all_pools_idat_gtc.txt",
            col.names = FALSE,quote = FALSE,row.names = FALSE)
