# calculate how much variance from metadata is explained by the PCs
library(vegan)
library(tidyverse)
expression_pcs_untreated = read.table("../../data/for_tensorQTL/90/full_metadata_scaled_untreated_Not_proliferating.txt")
expression_pcs_untreated = t(expression_pcs_untreated) %>%
  as.data.frame() %>%
  rownames_to_column("line") %>%
  mutate(across(contains("PC"), as.numeric),
         line = str_replace(line, "\\.", "-"))


donor_metadata = read_csv("../../../OTAR2065_differentiation_efficiency/data/donor_metadata_complete_with_imputed_sex.csv")
metadata_info = read_csv("../../../OTAR2065_differentiation_efficiency/data/metadata_info_hipsci_IPMAR.csv") %>%
  mutate(Disease.Status = case_when(Disease.Status == "Control" ~ "Normal",
                                    Disease.Status == "EOAD" ~ "AD",
                                    Disease.Status == "LOAD" ~ "AD",
                                    .default = "Normal")) 
merged_metadata = donor_metadata %>%
  left_join(metadata_info %>%   select(donor,line), by = "donor") %>% 
  distinct() %>%
  mutate(line = case_when(donor == "Arlene" ~ "Arlene-003",
                          donor == "Bertha" ~ "Bertha-004",
                          donor == "Cindy" ~ "Cindy-005",
                          donor == "Dexter" ~ "Dexter-006",
                          donor == "Fiona" ~ "Fiona-010",
                          donor == "Gilma" ~ "Gilma-009",
                          donor == "Hector" ~ "Hector-011",
                          donor == "Imani" ~ "Imani-012",
                          donor == "Javier" ~ "Javier-013",
                          donor == "Keoni" ~ "Keoni-014",
                          donor == "Mindy" ~ "Mindy-016",
                          donor == "Nestor" ~ "Nestor-017",
                          donor == "Olaf" ~ "Olaf-018",
                          donor == "Qiana" ~ "Qiana-022",
                          donor == "eorc" ~ "eorc_2",
                          donor == "hayt" ~ "hayt_3",
                          donor == "jaqg" ~ "jaqg",
                          donor == "jotn" ~ "jotn_2",
                          donor == "nurh" ~ "nurh",
                          donor == "pitg" ~ "pitg_2",
                          donor == "qecv" ~ "qecv_2",
                          donor == "uict" ~ "uict_1",
                          .default = line),
         Disease.Status = case_when(Disease.Status == "Control" ~ "Normal",
                                    Disease.Status == "EOAD" ~ "AD",
                                    Disease.Status == "LOAD" ~ "AD",
                                    .default = "Normal"))




# Keep only the PC rows/lines we want, in their current order
exprpcs = paste0("PC",1:63)
pc_data <- expression_pcs_untreated %>%
  select(line, exprpcs,genotypePC1,genotypePC2) # this contains expression and genotype PCs fitted
pcs = colnames(pc_data[-1])
# Align metadata to PC order 
metadata_aligned <- pc_data %>%
  select(line) %>%
  left_join(merged_metadata, by = "line")

metadata_aligned %>%
  filter(if_any(-line, is.na)) %>%
  select(line)
# Drop silent Nas (we filtered donors by number of cells before building eQTL input matrix)
vars_meta <- c("Disease.Status", "data_source", "Sex")

keep <- complete.cases(metadata_aligned[, vars_meta]) &
  complete.cases(pc_data[, pcs])

metadata_aligned2 <- metadata_aligned[keep, ]
pc_data2 <- pc_data[keep, ]

# We don't check pool because pools are aggregated in the pseudobulk for eQTL

meta_matrix <- model.matrix(
  ~ Disease.Status + data_source + Sex,
  data = metadata_aligned2
)

rda_res <- rda(
  reformulate(pcs, response = "meta_matrix"),
  data = pc_data2[, pcs]
)

# Extract variance explained:
RsquareAdj(rda_res)

# $r.squared
# [1] 0.892
# adj.r.squared
# 0.832
# Significance:
anova(rda_res, permutations = 999)
# 0.001
