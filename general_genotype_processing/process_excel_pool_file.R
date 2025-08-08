library(openxlsx)

donors = read.xlsx("../../OTAR2065_Hipsci_lines_in_Pools_160822.xlsx",sheet = 2, colNames = FALSE)

metadata = data.frame("sample_ID" = unlist(donors[1,]),
                      "number_of_donors" = unlist(donors[2,]),
                      "donor_IDs" = apply(donors[-c(1:2),],2,function(x) paste(na.omit(unlist(x)),collapse = ";")))
write.table(
  metadata,
  "../../sampleMetadata.txt",
  sep = "\t",
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)
