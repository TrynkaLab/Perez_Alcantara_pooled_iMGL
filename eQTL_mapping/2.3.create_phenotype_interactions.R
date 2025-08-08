# create phenotype interactions

library(readr)
library(tidyverse)
library(stringr)

output_dir = "../../data/for_tensorQTL/"
dir.create(output_dir)

args = commandArgs(trailingOnly=TRUE)
message(length(args)," arguments provided")
if (length(args)<2) {
  stop("You need to detail 3 arguments: treatment, PC number and phenotype.n",
       call. = FALSE)
} else if (length(args) == 3) {
  treat = args[1]
  pcnumber = as.numeric(args[2])
  phenotype = args[3]
  
}

# treat="untreated"
# pcnumber = 15
# phenotype = "phagocytosis"



if(treat!="IFN"){
  
  scaled = list()
  metadata = list()
for(condition in c("Not_proliferating","Proliferating")){
  metadata[[condition]] = read.table( 
    paste0(output_dir,"/",pcnumber,"/tensorQTL_metadata_sum_sizefactorsNorm_log2_",treat,"_",condition,".txt"), 
    header = TRUE)
  
  scaled[[condition]] = readr::read_table(paste0(output_dir,"/expr_sum_sizefactorsNorm_log2_scaled_centered_",
                                                 treat, "_", condition,".bed.gz"))
  colnames(scaled[[condition]]) = c("chr","start","end","gene_id",colnames(metadata[[condition]]))
  
  # same number of donors
  ncol(metadata[[condition]]) == ncol(scaled[[condition]][-1:-4])
  


message("Creating interaction dataframe")
message("Reading in phenotype data and calculating mean and median per line and treatment, across pools")
if(phenotype == "phagocytosis"){
  mean_pheno = readr::read_csv(paste0("../../../OTAR2065_phenotypic_QTLs/data/results/",
                                      phenotype,
                                      "/1.check_donor_proportions/line_prop_changes_per_well.csv")) %>%
    dplyr::filter(!log_fraction_mean %in% c(NA,Inf, -Inf)) %>%
    dplyr::filter(treatment  ==  treat) %>%
    dplyr::group_by(line) %>%
    dplyr::summarise(scaled_fraction_mean = mean(scaled_log_fraction)) 
  
  median_pheno = readr::read_csv(paste0("../../../OTAR2065_phenotypic_QTLs/data/results/",
                                        phenotype,
                                        "/1.check_donor_proportions/line_prop_changes_per_well.csv")) %>%
    dplyr::filter(!log_fraction_mean %in% c(NA,Inf, -Inf)) %>%
    dplyr::filter(treatment  ==  treat) %>%
    dplyr::group_by(line) %>%
    dplyr::summarise(scaled_fraction_median = median(scaled_log_fraction)) 
  
  # histograms
  
  p1 = ggplot(mean_pheno,aes(x=scaled_fraction_mean)) + 
    geom_histogram(bins = 50) +
    theme_minimal() + 
    ggtitle("Mean scaled fraction: ", phenotype)
  
  p2 = ggplot(median_pheno,aes(x=scaled_fraction_median)) + 
    geom_histogram(bins = 50) +
    theme_minimal() + 
    ggtitle("Median scaled fraction: ", phenotype)
  p1 / p2
  
# Inverse normal transformation
  # https://www.biostars.org/p/80597/
  mean_inverse = mean_pheno %>%
    dplyr::mutate(INT_mean_log_fraction = qnorm((rank(log_fraction_mean,
                                                      na.last="keep")-0.5)/sum(!is.na(log_fraction_mean)))) %>%
    dplyr::select(!log_fraction_mean)
  
  median_inverse = median_pheno %>%
    dplyr::mutate(INT_median_log_fraction = qnorm((rank(log_fraction_median,
                                                      na.last="keep")-0.5)/sum(!is.na(log_fraction_median)))) %>%
    dplyr::select(!log_fraction_median)
  
  p3 = ggplot(mean_inverse,aes(x=INT_mean_log_fraction)) + 
    geom_histogram(bins = 50) +
    theme_minimal() + 
    ggtitle("INT mean log fraction: ", phenotype)
  
  p4 = ggplot(median_inverse,aes(x=INT_median_log_fraction)) + 
    geom_histogram(bins = 50) +
    theme_minimal() + 
    ggtitle("INT median log fraction: ", phenotype)
  (p1+ p2) / (p3 + p4)
  

}
else{
  mean_pheno =readr::read_csv(paste0("../../../OTAR2065_phenotypic_QTLs/data/results/",
phenotype,
"/1.check_line_proportions/line_prop_changes_per_well.csv")) %>%
    dplyr::mutate(treatment = case_when(treatment == "Untreated" ~ "untreated",
                                        .default = treatment)) %>%
    dplyr::filter(!log_fraction_mean %in% c(NA,Inf, -Inf)) %>%
    dplyr::filter(treatment == treat) %>%
    dplyr::filter(condition =="3nM_C5a") %>%
    dplyr::group_by(line) %>%
    dplyr::summarise(scaled_fraction_mean = mean(scaled_log_fraction)) 
  
  median_pheno = readr::read_csv(paste0("../../../OTAR2065_phenotypic_QTLs/data/results/1.check_line_proportions/line_prop_changes_per_well.csv")) %>%
    dplyr::mutate(treatment = case_when(treatment == "Untreated" ~ "untreated",
                                        .default = treatment)) %>%
    dplyr::filter(!log_fraction_mean %in% c(NA,Inf, -Inf)) %>%
    dplyr::filter(treatment == treat) %>%
    dplyr::filter(condition =="3nM_C5a") %>%
    dplyr::group_by(line) %>%
    dplyr::summarise(scaled_fraction_median = median(scaled_log_fraction)) 
}


message("Subsetting donors in expression and covariate files to those in the interaction file")

# filtering interaction files first
mean_pheno = mean_pheno %>%
  dplyr::filter(line %in% colnames(metadata[[condition]]))
median_pheno = median_pheno %>%
  dplyr::filter(line %in% colnames(metadata[[condition]]))

metadata[[condition]] = metadata[[condition]] %>%
  dplyr::select(mean_pheno$line)

scaled[[condition]] = scaled[[condition]] %>%
  dplyr::select("chr","start","end","gene_id",mean_pheno$line)

message("Saving all dataframes")

mean_pheno %>%
  readr::write_tsv(paste0(output_dir,"/",pcnumber,"/interaction_df_mean_",
                          phenotype,"_",
                          treat, "_", condition,".txt"), 
                   col_names = TRUE)

median_pheno %>%
  readr::write_tsv(paste0(output_dir,"/",pcnumber,"/interaction_df_median_",
                          phenotype,"_",
                          treat, "_", condition,".txt"), 
                   col_names = TRUE)

scaled[[condition]] %>%
  dplyr::rename('#chr' = "chr") %>%
  write.table(.,paste0(output_dir,pcnumber,"/expr_sum_sizefactorsNorm_log2_scaled_centered_",phenotype,"_",
                       treat, "_", condition,"_for_interaction.bed"),
              sep = "\t", quote = F, col.names = T, row.names = F)
metadata[[condition]] %>%
write.table(  .,
              paste0(output_dir,pcnumber,"/tensorQTL_metadata_sum_sizefactorsNorm_log2_",phenotype,"_",
                     treat,"_",condition,"_for_interaction.txt"), 
              sep = "\t", quote = F, col.names = T, row.names = T)
}
}else{
  scaled = list()
  metadata = list()
  for(condition in c("Not_proliferating")){
    metadata[[condition]] = read.table( 
      paste0(output_dir,"/",pcnumber,"/tensorQTL_metadata_sum_sizefactorsNorm_log2_",treat,"_",condition,".txt"), 
      header = TRUE)
    
    scaled[[condition]] = readr::read_table(paste0(output_dir,"/expr_sum_sizefactorsNorm_log2_scaled_centered_",
                                                   treat, "_", condition,".bed.gz"))
    colnames(scaled[[condition]]) = c("chr","start","end","gene_id",colnames(metadata[[condition]]))
    
    # same number of donors
    ncol(metadata[[condition]]) == ncol(scaled[[condition]][-1:-4])
    
    
    
    message("Creating interaction dataframe")
    message("Reading in phenotype data and calculating mean and median per line and treatment, across pools")
    if(phenotype == "phagocytosis"){
      mean_pheno = readr::read_csv(paste0("../../../OTAR2065_phenotypic_QTLs/data/results/",
                                          phenotype,
                                          "/1.check_donor_proportions/line_prop_changes_per_well.csv")) %>%
        dplyr::filter(!log_fraction_mean %in% c(NA,Inf, -Inf)) %>%
        dplyr::filter(treatment  ==  treat) %>%
        dplyr::group_by(line) %>%
        dplyr::summarise(scaled_fraction_mean = mean(scaled_log_fraction)) 
      
      median_pheno = readr::read_csv(paste0("../../../OTAR2065_phenotypic_QTLs/data/results/",
                                            phenotype,
                                            "/1.check_donor_proportions/line_prop_changes_per_well.csv")) %>%
        dplyr::filter(!log_fraction_mean %in% c(NA,Inf, -Inf)) %>%
        dplyr::filter(treatment  ==  treat) %>%
        dplyr::group_by(line) %>%
        dplyr::summarise(scaled_fraction_median = median(scaled_log_fraction)) 
    }
    else{
      mean_pheno =readr::read_csv(paste0("../../../OTAR2065_phenotypic_QTLs/data/results/",
phenotype,"/1.check_line_proportions/line_prop_changes_per_well.csv")) %>%
        dplyr::mutate(treatment = case_when(treatment == "Untreated" ~ "untreated",
                                            .default = treatment)) %>%
        dplyr::filter(!log_fraction_mean %in% c(NA,Inf, -Inf)) %>%
        dplyr::filter(treatment == treat) %>%
        dplyr::filter(condition =="3nM_C5a") %>%
        dplyr::group_by(line) %>%
        dplyr::summarise(scaled_fraction_mean = mean(scaled_log_fraction)) 
      
      median_pheno = readr::read_csv(paste0("../../../OTAR2065_phenotypic_QTLs/data/results/",phenotype,"/1.check_line_proportions/line_prop_changes_per_well.csv")) %>%
        dplyr::mutate(treatment = case_when(treatment == "Untreated" ~ "untreated",
                                            .default = treatment)) %>%
        dplyr::filter(!log_fraction_mean %in% c(NA,Inf, -Inf)) %>%
        dplyr::filter(treatment == treat) %>%
        dplyr::filter(condition =="3nM_C5a") %>%
        dplyr::group_by(line) %>%
        dplyr::summarise(scaled_fraction_median = median(scaled_log_fraction)) 
    }
    
    
    message("Subsetting donors in expression and covariate files to those in the interaction file")
    
    # filtering interaction files first
    mean_pheno = mean_pheno %>%
      dplyr::filter(line %in% colnames(metadata[[condition]]))
    median_pheno = median_pheno %>%
      dplyr::filter(line %in% colnames(metadata[[condition]]))
    
    metadata[[condition]] = metadata[[condition]] %>%
      dplyr::select(mean_pheno$line)
    
    scaled[[condition]] = scaled[[condition]] %>%
      dplyr::select("chr","start","end","gene_id",mean_pheno$line)
    
    message("Saving all dataframes")
    
    mean_pheno %>%
      readr::write_tsv(paste0(output_dir,"/",pcnumber,"/interaction_df_mean_",
                              phenotype,"_",
                              treat, "_", condition,".txt"), 
                       col_names = TRUE)
    
    median_pheno %>%
      readr::write_tsv(paste0(output_dir,"/",pcnumber,"/interaction_df_median_",
                              phenotype,"_",
                              treat, "_", condition,".txt"), 
                       col_names = TRUE)
    
    scaled[[condition]] %>%
      dplyr::rename('#chr' = "chr") %>%
      write.table(.,paste0(output_dir,pcnumber,"/expr_sum_sizefactorsNorm_log2_scaled_centered_",phenotype,"_",
                           treat, "_", condition,"_for_interaction.bed"),
                  sep = "\t", quote = F, col.names = T, row.names = F)
    metadata[[condition]] %>%
      write.table(  .,
                    paste0(output_dir,pcnumber,"/tensorQTL_metadata_sum_sizefactorsNorm_log2_",phenotype,"_",
                           treat,"_",condition,"_for_interaction.txt"), 
                    sep = "\t", quote = F, col.names = T, row.names = T)
  }
}
