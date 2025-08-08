library(openxlsx)
library(stringr)

donors = read.xlsx("../../data/OTAR2065_Hipsci_lines_in_Pools_july2023.xlsx",sheet = 2, colNames = FALSE)


metadata = data.frame("sample_ID" = unlist(donors[1,]),
                      "number_of_donors" = unlist(donors[2,]),
                      "donor_IDs" = apply(donors[-c(1:2),],2,function(x) paste(na.omit(sort(unlist(x))),collapse = ";")))

# number of unique donors
length(unlist(str_split(metadata[13,3],pattern = ";")))

write.table(
  metadata,
  "../../data/sampleMetadata.txt",
  sep = "\t",
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)
