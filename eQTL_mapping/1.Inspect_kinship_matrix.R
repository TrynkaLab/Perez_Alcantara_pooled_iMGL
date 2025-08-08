## Investigating kinship matrix for this experiment


args <- commandArgs(trailingOnly = TRUE)

kinship = args[1]
output = args[2]
kinship =  "../../data/kinship/all_pools.genotype.MAF05.hg38.kin0"

message(paste0("Reading kinship matrix from ",kinship))

kinship = read.table(kinship)

message("Note that KING kinship coefficients are scaled such that duplicate samples have kinship 0.5, not 1.\n
First-degree relations (parent-child, full siblings) correspond to ~0.25, second-degree relations correspond to ~0.125, etc")

pdf(file = output, width = 8, height = 8)
    hist(kinship[,6], main = "histogram of pairwise kinship coefficients", 
    xlab = "scaled kinship coefficient")
    
 dev.off()
 
 related  = sum(kinship[,8] > 0.125)
 message(paste0("There are ", related, " individuals with KINSHIP > 0.125"))
 
