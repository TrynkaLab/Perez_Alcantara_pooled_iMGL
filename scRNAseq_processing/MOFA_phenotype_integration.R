## MOFA analysis

## Integrating gene expression and phenotypic results
## untreated only to see if it clarifies a bit more of the variation
# setwd("G:/My Drive/server/OTAR2065_differentiation_efficiency/code/otar2065_differentiation_efficiency") # for local computer
library(tidyverse)
library(MOFA2)
library(patchwork)
library(ggrepel)
source("./helpers.R")

outputDir = "../../data/results/7.MOFA_phenotype_integration/"
dir.create(outputDir, recursive = TRUE)

### my data ###

MOFA_data = list()
MOFA_data[["pseudobulk"]] = readr::read_delim("../../data/results/5.prepare_differential_expression_pseudobulks/combined_pseudobulk_raw_counts_treatmentxprolifxidxpool.txt",
                                              delim = "\t") %>%
  dplyr::relocate(gene) %>%
  tibble::column_to_rownames("gene") %>%
  as.data.frame()

# metadata
sex_meta = readr::read_csv("../../data/donor_metadata_complete_with_imputed_sex.csv")
prs = read_tsv("../../../hipsci_genotype_processing/data/prs_ad_bellenguez/hipsci_polygenic_score_AD_Bellenguez_withAPOE.tsv") %>%
  dplyr::select(line, APOE_sum_scaled, prs_scaled, full_PRS) %>%
  dplyr::rename(donor=line) %>%
  dplyr::mutate(line = case_when(donor== "Arlene" ~ "Arlene-003",
                                 donor=="Cindy" ~ "Cindy-005",
                                 donor=="Dexter" ~ "Dexter-006",
                                 donor=="Fiona" ~ "Fiona-010",
                                 donor=="Gilma" ~ "Gilma-009",
                                 donor=="Hector" ~ "Hector-011",
                                 donor=="Imani" ~ "Imani-012",
                                 donor=="Javier" ~ "Javier-013",
                                 donor=="Keoni" ~ "Keoni-014",
                                 donor=="Olaf" ~ "Olaf-018",
                                 donor=="Bertha" ~ "Bertha-004",
                                 donor=="Mindy" ~ "Mindy-016",
                                 donor=="Qiana" ~ "Qiana-022",
                                 donor=="Nestor" ~ "Nestor-017",
                                 .default = donor)) %>%
  dplyr::mutate(donor =  str_split_i(donor,"_",i=1))

microglia_specific_prs = list()
for(treat in c("untreated","LPS","IFN")){
  microglia_specific_prs[[treat]] = read_tsv(paste0("../../../hipsci_genotype_processing/data/prs_ad_bellenguez/PRSice/microglia-specific/hipsci_polygenic_score_AD_Bellenguez_withAPOE_microglia_",
                                                    treat,
                                                    ".tsv"))
  microglia_specific_prs[[treat]]=   microglia_specific_prs[[treat]] %>%
    dplyr::select(line,APOE_sum_scaled,prs_scaled, full_PRS) %>%
    dplyr::mutate(treatment = treat)
}

microglia_specific_prs = do.call("rbind",microglia_specific_prs)
microglia_specific_prs = microglia_specific_prs %>%
  tidyr::pivot_wider(id_cols = line,
                     names_from = treatment,
                     values_from = c(APOE_sum_scaled, prs_scaled, full_PRS),
                     names_glue = "{treatment}_{.value}_microglia_PRS"
  )



pseudobulk_meta = readr::read_delim("../../data/results/5.prepare_differential_expression_pseudobulks/metadata_treatmentxprolifxidxpool.txt",
                                    delim = "\t") %>%
  tibble::column_to_rownames(var = "cols_names") %>%
  dplyr::rename(ncells = count,line= donor_id) %>%
  dplyr::mutate(pseudobulk_colnames=paste(treatment,proliferation_status,line,pool,sep = "_")) %>%
  dplyr::mutate(donor =  if_else(stringr::str_detect(line,"_"),
                                 stringr::str_split_i(line,"_",i=1),
                                 stringr::str_split_i(line,"-",i=1))) 

pseudobulk_meta = pseudobulk_meta %>%
  dplyr::left_join(sex_meta, by = "donor")

# filter out pseudobulks under 100 cells
discarded = pseudobulk_meta$ncells < 100
table(discarded)

MOFA_data[["pseudobulk"]]  = MOFA_data[["pseudobulk"]][,!discarded]
metadata = pseudobulk_meta[!discarded,]
length(unique(metadata$line)) # 194/261, 74%

# filter out proliferating cluster
discarded = metadata$proliferation_status == "Proliferating"
MOFA_data[["pseudobulk"]]  = MOFA_data[["pseudobulk"]][,!discarded]
metadata = metadata[!discarded,]
length(unique(metadata$line)) # 194/261, 74%

# subset to tested genes that have been filtered in the various differential expression analyses
tested_genes = readr::read_rds("../../data/results/5.1.Diff_expr_limma/additive_res.rds")
tested_genes = union(tested_genes$IFNvsUntreated$symbol,union(tested_genes$LPSvsUntreated$symbol,tested_genes$IFNvsLPS$symbol))
MOFA_data[["pseudobulk"]] = MOFA_data[["pseudobulk"]][rownames(MOFA_data[["pseudobulk"]]) %in% tested_genes,]


# normalize expression
dge = edgeR::DGEList(counts=MOFA_data[["pseudobulk"]] )
dge = edgeR::calcNormFactors(dge)
### create design matrix
metadata$treatment = factor(metadata$treatment)
metadata$treatment = relevel(metadata$treatment,ref = "untreated")
metadata$Sex = factor(metadata$Sex)
metadata$Sex = relevel(metadata$Sex,ref = "Female")
metadata$log_ncells = log(metadata$ncells)

mat = model.matrix(~ 0 + treatment,metadata)
mat = as.data.frame(mat)
qr(mat)$rank
ncol(mat)
limma::is.fullrank(mat)
colnames(mat)
mat = mat[,colSums(mat)>0]
colnames(mat)
limma::is.fullrank(mat)

colnames(mat) = make.names(colnames(mat))

v = limma::voom(dge, design = mat, plot=TRUE)
# remove batch effects as recommended, and also covariates we're not interested in
v_corrected = limma::removeBatchEffect(v,batch = metadata$pool,batch2 = metadata$Sex,
                                       covariates =metadata$log_ncells ,design = mat)

# save normalised data
MOFA_data[["pseudobulk"]] = v_corrected
rm(dge,v,mat)
pseudobulk_rownames = rownames(MOFA_data[["pseudobulk"]]) # save this for later

# load phenotype matrices
# I don't include prop_unadjusted_min_value because that has a different distribution to the scaled log fraction
MOFA_data[["phagocytosis"]] = readr::read_csv("../../../OTAR2065_phenotypic_QTLs/data/results/phagocytosis/1.check_line_proportions/line_prop_changes_per_well.csv") %>%
  dplyr::select(line,donor,scaled_log_fraction,pool,treatment,replicate) %>%
  dplyr::filter(line!="sh5y5y") %>%
  dplyr::group_by(line,donor,pool,treatment) %>%
  dplyr::reframe(mean_scaled_log_fraction = mean(scaled_log_fraction)) %>% # mean of replicates
  dplyr::ungroup() %>%
  dplyr::mutate(pseudobulk_colnames=paste(treatment,"Not_proliferating",line,pool,sep = "_")) %>%
  dplyr::select(pseudobulk_colnames,mean_scaled_log_fraction) %>%
  tibble::column_to_rownames("pseudobulk_colnames") %>%
  t(.)
phago_rownames = rownames(MOFA_data[["phagocytosis"]])  # save this for later



# adding columns present in pseudobulk but not in phagocytosis to the phagocytosis object, filling in with NAs
MOFA_data[["phagocytosis"]] = fill_columns_df2_with_df1_NA(MOFA_data[["pseudobulk"]],MOFA_data[["phagocytosis"]])

# same in the other direction
# MOFA_data[["pseudobulk"]] = fill_columns_df2_with_df1_NA(MOFA_data[["phagocytosis"]],MOFA_data[["pseudobulk"]])
MOFA_data[["phagocytosis"]] = MOFA_data[["phagocytosis"]] %>%
  dplyr::select(colnames(MOFA_data[["pseudobulk"]]))
# drop phagocytosis columns not present in pseudobulk

# check columns are identical
identical(colnames(MOFA_data[["phagocytosis"]] ),colnames(MOFA_data[["pseudobulk"]] ) )

# same with migration

# adding to pseudobulk columns present in phagocytosis, filling in with NAs
# MOFA_data[["pseudobulk"]] = MOFA_data[["pseudobulk"]]  %>%
#   dplyr::bind_cols(MOFA_data[["phagocytosis"]][1, 
#                                              setdiff(colnames(MOFA_data[["phagocytosis"]] ),
#                                                      colnames(MOFA_data[["pseudobulk"]] ))]) 


MOFA_data[["migration"]] = readr::read_csv("../../../OTAR2065_phenotypic_QTLs/data/results/migration/1.check_line_proportions/line_prop_changes_per_well.csv") %>%
  dplyr::select(line,donor,scaled_log_fraction,pool,treatment,replicate) %>%
  dplyr::group_by(line,donor,pool,treatment) %>%
  dplyr::reframe(mean_scaled_log_fraction = mean(scaled_log_fraction)) %>% # mean of replicates
  dplyr::ungroup() %>%
  dplyr::mutate(pseudobulk_colnames=paste(treatment,"Not_proliferating",line,pool,sep = "_")) %>% 
  dplyr::select(pseudobulk_colnames,mean_scaled_log_fraction) %>%
  tibble::column_to_rownames("pseudobulk_colnames") %>%
  t(.) 
migr_rownames = rownames(MOFA_data[["migration"]])  # save this for later

# triple fill-in now
MOFA_data[["migration"]] = fill_columns_df2_with_df1_NA(MOFA_data[["pseudobulk"]],MOFA_data[["migration"]])
# MOFA_data[["pseudobulk"]] = fill_columns_df2_with_df1_NA(MOFA_data[["migration"]],MOFA_data[["pseudobulk"]])
# MOFA_data[["phagocytosis"]] = fill_columns_df2_with_df1_NA(MOFA_data[["pseudobulk"]],MOFA_data[["phagocytosis"]])

MOFA_data[["migration"]] = MOFA_data[["migration"]] %>%
  dplyr::select(colnames(MOFA_data[["pseudobulk"]]))

identical(colnames(MOFA_data[["migration"]] ),colnames(MOFA_data[["pseudobulk"]] ) )
identical(colnames(MOFA_data[["phagocytosis"]] ),colnames(MOFA_data[["pseudobulk"]] ) )

# Same for proliferation - only microglia vs premac because we only have microglial expression
MOFA_data[["proliferation"]] = readr::read_csv("../../../OTAR2065_WGS_iPSC_premac_micro/data/results/1.2.scale_proliferation/line_prop_changes_microglia_premac.csv") %>%
  dplyr::mutate(treatment = stringr::str_remove(treatment,"g"),
                donor = if_else(stringr::str_detect(line,"-"), 
                                stringr::str_split_i(line,"-",1),
                                stringr::str_split_i(line,"_",1))) %>% 
  dplyr::select(line,donor,scaled_log_fraction,pool,treatment) %>% # there are several replicates because of different precursor ages (pheno vs scRNA-seq)
  dplyr::group_by(line,donor,pool,treatment) %>%
  dplyr::reframe(mean_scaled_log_fraction = mean(scaled_log_fraction)) %>% # mean of replicates
  dplyr::ungroup() %>%
  dplyr::mutate(pseudobulk_colnames=paste(treatment,"Not_proliferating",line,pool,sep = "_")) %>% 
  dplyr::select(pseudobulk_colnames,mean_scaled_log_fraction) %>%
  tibble::column_to_rownames("pseudobulk_colnames") %>%
  t(.) 
prolif_rownames = rownames(MOFA_data[["proliferation"]])  # save this for later


# triple fill-in now
MOFA_data[["proliferation"]] = fill_columns_df2_with_df1_NA(MOFA_data[["pseudobulk"]],MOFA_data[["proliferation"]])
# MOFA_data[["pseudobulk"]] = fill_columns_df2_with_df1_NA(MOFA_data[["proliferation"]],MOFA_data[["pseudobulk"]])
MOFA_data[["phagocytosis"]] = fill_columns_df2_with_df1_NA(MOFA_data[["pseudobulk"]],MOFA_data[["phagocytosis"]])
MOFA_data[["migration"]] = fill_columns_df2_with_df1_NA(MOFA_data[["pseudobulk"]],MOFA_data[["migration"]])

MOFA_data[["proliferation"]] = MOFA_data[["proliferation"]] %>%
  dplyr::select(colnames(MOFA_data[["pseudobulk"]]))

identical(colnames(MOFA_data[["proliferation"]] ),colnames(MOFA_data[["pseudobulk"]] ) )
identical(colnames(MOFA_data[["phagocytosis"]] ),colnames(MOFA_data[["pseudobulk"]] ) )
identical(colnames(MOFA_data[["migration"]] ),colnames(MOFA_data[["pseudobulk"]] ) )


##### Creating MOFA object ######

# saving new metadata for groups
metadata = data.frame(pseudobulk_colnames = colnames(MOFA_data$pseudobulk)) %>%
  dplyr::as_tibble() %>%
  dplyr::left_join(metadata, by = "pseudobulk_colnames") %>%
  dplyr::filter(stringr::str_detect(pseudobulk_colnames,"Not_proliferating")) %>%
  dplyr::mutate(treatment =  stringr::str_split_i(pseudobulk_colnames,"_",1),
                proliferation_status = paste(stringr::str_split_i(pseudobulk_colnames,"_",2),stringr::str_split_i(pseudobulk_colnames,"_",3), sep = "_"),
                donor = stringr::str_split_i(pseudobulk_colnames,"_",4),
                pool = stringr::str_split_i(pseudobulk_colnames,"_",-1),
                line = stringr::str_replace((stringr::str_replace((stringr::str_replace(pseudobulk_colnames,treatment,"")),
                                                                  "_Not_proliferating_","")),paste0("_",pool),"")) %>%
  dplyr::select(!c(Sex,Age,Disease.Status,data_source,is_sex_imputed)) %>%
  dplyr::left_join(sex_meta)


metadata = metadata %>%
  dplyr::arrange(match(pseudobulk_colnames, colnames(MOFA_data$pseudobulk)))  %>%
  dplyr::left_join(prs %>% dplyr::select(-donor),by="line") %>%
  dplyr::left_join(microglia_specific_prs,by="line")

identical(metadata$pseudobulk_colnames, colnames(MOFA_data$pseudobulk))

#  adding phenotype info

metadata = metadata %>%
  dplyr::left_join(data.frame(pseudobulk_colnames = rownames(t(MOFA_data$phagocytosis)),
                              scaled_phagocytosis = t(MOFA_data$phagocytosis)[,1])) %>%
  dplyr::left_join(data.frame(pseudobulk_colnames = rownames(t(MOFA_data$migration)),
                              scaled_migration = t(MOFA_data$migration)[,1]))  %>%
  dplyr::left_join(data.frame(pseudobulk_colnames = rownames(t(MOFA_data$proliferation)),
                              scaled_proliferation = t(MOFA_data$proliferation)[,1]))

# metadata as data frame
metadata = data.frame(metadata)
rownames(metadata) = colnames(MOFA_data$pseudobulk)

# make all tibbles matrices
MOFA_data$pseudobulk = as.matrix(MOFA_data$pseudobulk)
rownames(MOFA_data$pseudobulk) = pseudobulk_rownames
MOFA_data$phagocytosis = as.matrix(MOFA_data$phagocytosis)
rownames(MOFA_data$phagocytosis) = phago_rownames
MOFA_data$migration = as.matrix(MOFA_data$migration)
rownames(MOFA_data$migration) = migr_rownames
MOFA_data$proliferation = as.matrix(MOFA_data$proliferation)
rownames(MOFA_data$proliferation) = prolif_rownames

###### Training data ########
## creating MOFA
MOFA_data2 = list(MOFA_data$pseudobulk)
names(MOFA_data2) = "RNA"
MOFAobject = create_mofa(MOFA_data2)


### changing model defaults

plot_data_overview(MOFAobject) 
data_opts = get_default_data_options(MOFAobject)
head(data_opts)
# we'll scale and center the matrices
data_opts$scale_views = TRUE

model_opts <- get_default_model_options(MOFAobject)
head(model_opts)
# the authors recommend to leave the default model options

train_opts <- get_default_training_options(MOFAobject)
head(train_opts)
train_opts$convergence_mode = "slow" ## default is fast, for exploratory analysis
## building and training ######

### CONSIDER REMOVING GENES THAT DON'T EXPLAIN LOT OF VARIANCE (HIGH VAR GENES FROM SEURAT) - RETAIN EQTL GENES
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)
MOFAobject@training_options$save_interrupted = TRUE
outfile = file.path(outputDir,"model_pseudobulk_not_proliferating.hdf5")
MOFAobject.trained = run_mofa(MOFAobject, outfile,use_basilisk = TRUE)

MOFAobject.trained = load_model(outfile)

# adding metadata
metadata = metadata %>%
  dplyr::rename(sample = pseudobulk_colnames
  ) %>%
  dplyr::mutate(group = "group1") %>%
  dplyr::select(sample, group, pool, ncells, proliferation_status, line,treatment,
                donor,  APOE_sum_scaled , prs_scaled ,     full_PRS,
                setdiff(colnames(microglia_specific_prs ),"line"),
                scaled_phagocytosis, scaled_migration, scaled_proliferation) %>%   # fix sex errors
  dplyr::mutate(donor = case_when(stringr::str_detect(donor,"-") ~ str_split_i(donor,"-",1),
                                  .default = str_split_i(donor,"_",1))) %>%
  dplyr::left_join(sex_meta)

metadata$log_ncells = log10(metadata$ncells)

write_csv(metadata,paste0(outputDir,"MOFA_metadata_not_proliferating.csv"))
samples_metadata(MOFAobject.trained) = metadata

# save full object
write_rds(MOFAobject.trained,paste0(outputDir,"MOFA_model_not_proliferating.rds"))
MOFAobject.trained = read_rds(paste0(outputDir,"MOFA_model_not_proliferating.rds"))

plot_factor_cor(MOFAobject.trained) # largely uncorrelated, as expected

####  variance decomposition 
p1 = plot_variance_explained(MOFAobject.trained, max_r2=15) + 
  theme(axis.text.x = element_text(angle = 90))

plot(p1)

var_explained = calculate_variance_explained(MOFAobject.trained)
var_explained$r2_per_factor
# RNA
# Factor1  31.2894877
# Factor2  26.3614836
# Factor3   3.1809216
# Factor4   1.1270328
# Factor5   0.9465870
# Factor6   0.8925818
# Factor7   0.6215800
# Factor8   0.5173611
# Factor9   0.4997446
# Factor10  0.4421265
# Factor11  0.3688047
# Factor12  0.3373301
# Factor13  0.3241385
# Factor14  0.2407988
# Factor15  0.1542357
# without fitting pools as groups there are barely any factors that explain variance across phenotypes

## total variance explained by the model
plot_variance_explained(MOFAobject.trained, plot_total = T)[[2]] +
  theme(axis.text.x = element_text(angle = 90))
# expression explains most variance, because it has more factors
#  followed by proliferation, phagocytosis and migration

####### characterizing factors
correlate_factors_with_covariates(MOFAobject.trained, 
                                  covariates = c("Sex","log_ncells","treatment","Age","pool",
                                                 "Disease.Status","donor",
                                                 "APOE_sum_scaled" , "prs_scaled" ,     "full_PRS",
                                                 setdiff(colnames(microglia_specific_prs ),"line"),
                                                 "scaled_phagocytosis", "scaled_migration", "scaled_proliferation"), 
                                  plot="r",transpose = TRUE
)
correlations = correlate_factors_with_covariates(MOFAobject.trained, 
                                                 covariates = c("Sex","log_ncells","treatment","Age","pool",
                                                                "Disease.Status","donor",
                                                                "APOE_sum_scaled" , "prs_scaled" ,     "full_PRS",
                                                                setdiff(colnames(microglia_specific_prs ),"line"),
                                                                "scaled_phagocytosis", "scaled_migration", "scaled_proliferation"), 
                                                 plot="r",return_data = TRUE
)


summary(abs(correlations))

#  Aggregate prettier plots


p_cor = correlations %>%
  as.data.frame() %>%
  dplyr::select(c("Sex","log_ncells","treatment","Age","pool",
                  "Disease.Status","donor",
                  "APOE_sum_scaled" , "prs_scaled" ,     "full_PRS",
                  "untreated_full_PRS_microglia_PRS" ,
                  "scaled_phagocytosis", "scaled_migration", "scaled_proliferation")) %>%
  dplyr::rename("Log(Ncells)" = "log_ncells",
                "Disease status" =  "Disease.Status",
                "Treatment" = "treatment",
                "Pool" = "pool",
                "Donor" = "donor",
                "APOE risk" = "APOE_sum_scaled"  ,
                "AD polygenic risk" = "prs_scaled" ,   
                "AD full PRS" =  "full_PRS" ,
                "Microglia-specific AD full PRS"  = "untreated_full_PRS_microglia_PRS",
                "Phagocytosis phenotype"= "scaled_phagocytosis",
                "Migration phenotype"="scaled_migration",
                "Proliferation phenotype" = "scaled_proliferation"
  ) %>%
  dplyr::mutate(Factor = factor(paste("Factor",1:15),
                                levels = paste("Factor",1:15),
                                ordered = TRUE)) %>%
  tidyr::pivot_longer(cols = -Factor,names_to = "covariates",values_to = "r") %>%
  dplyr::mutate(covariates = factor(covariates,
                                    levels= c("Log(Ncells)" , "Sex","Age", "Treatment","Pool",
                                              "Disease status", "Donor",
                                              "APOE risk" ,
                                              "AD polygenic risk" ,   
                                              "AD full PRS",
                                              "Microglia-specific AD full PRS" ,
                                              "Phagocytosis phenotype",
                                              "Migration phenotype",
                                              "Proliferation phenotype"),
                                    ordered = TRUE)) %>%
  ggplot( aes(y=covariates, x=Factor, fill = r)) + 
  geom_raster() +
  geom_text(aes(label = round(r, digits = 2),
                col = if_else(abs(r) > 0.4, "white", if_else(abs(r) > 0.07, "grey20", NA_character_))  # Nested if_else for color
  )) + 
  scale_color_identity(guide = "none") + 
  scale_fill_gradient2(low = "steelblue",mid = "white", high = "#74121D",name = "R") +
  ggtitle("Correlation between MOFA factors and metadata") +
  theme_minimal() + 
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12)) +
  # facet_grid(cols = vars(`gene list` )) +
  theme(legend.position="right") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), # lables vertical
    strip.text.y = element_blank())   #remove facet bar on y

png(paste0(outputDir,"correlations.png"),width = 12,height = 5)
plot(p_cor)
dev.off()
# inspect factors 3 (phago) 7 (correlated to phago, microglial PRS and prolif pheno), and 4 (migr, prolif) and 8 (phago and microglial PRS)

plot_factor(MOFAobject.trained, 
            factors = 3, 
            color_by = "scaled_phagocytosis"
)

plot_factor(MOFAobject.trained, 
            factors = 7, 
            color_by = "scaled_phagocytosis"
)

plot_factor(MOFAobject.trained, 
            factors =7, 
            color_by = "untreated_full_PRS_microglia_PRS"
)


# features with strong association with the factor are expected to have large absolute values. 
# The sign of the weights indicates the direction of the effect: a positive weights indicates 
# that the feature has higher levels in the cells with positive factor values, and vice-versa.

p2 = plot_weights(MOFAobject.trained,
                  view = "RNA",
                  factor = 3,
                  nfeatures = 10,     # Top number of features to highlight, all negative
                  scale = T           # Scale weights from -1 to 1
) + ggtitle("Factor 3")

p3 = plot_weights(MOFAobject.trained,
                  view = "RNA",
                  factor = 4,
                  nfeatures = 10,     # Top number of features to highlight, all positive
                  scale = T           # Scale weights from -1 to 1
)  + ggtitle("Factor 4")
p4 = plot_weights(MOFAobject.trained,
                  view = "RNA",
                  factor = 8,
                  nfeatures = 10,     
                  scale = T           # Scale weights from -1 to 1
)  + ggtitle("Factor 8")

p5 = plot_weights(MOFAobject.trained,
                  view = "RNA",
                  factor =7,
                  nfeatures = 10,     
                  scale = T           # Scale weights from -1 to 1
)  + ggtitle("Factor 7")

pseudobulk_weights = plot_weights(MOFAobject.trained,
                                  view = "RNA",
                                  factor = "all",
                                  scale = T,return_data = TRUE           # Scale weights from -1 to 1
)

(p2 + p3) / (p4 + p5)



check_gene = v_corrected %>% t() %>% as_data_frame() %>% dplyr::summarise(gene_counts = colSums(.))
check_gene = check_gene %>%
  dplyr::mutate(gene = pseudobulk_rownames) %>%
  dplyr::arrange(desc(gene_counts))
check_gene$rank = 1:nrow(check_gene)
# getting quantiles
fn = ecdf(check_gene$gene_counts)
check_gene$quantile = fn(check_gene$gene_counts)
check_gene = check_gene %>%
  dplyr::left_join(pseudobulk_weights,by = join_by("gene" == "feature"))
p6 = check_gene %>%
  dplyr::group_by(factor) %>%
  slice_max(order_by = abs(value),n = 30) %>%
  dplyr::ungroup() %>%
  dplyr::filter(factor %in% paste0("Factor",1:20)) %>%
  ggplot(aes(x = factor,y = quantile, label = gene)) + 
  geom_point() + 
  geom_text_repel() +
  theme_bw() + 
  ggtitle("Top 30 genes per factor")

print(p6)

plot_top_weights(MOFAobject.trained,
                 view = "RNA",
                 factor = 1,  
                 nfeatures = 10,     
                 scale = T )

plot_top_weights(MOFAobject.trained,
                 view = "RNA",
                 factor = 3,  
                 nfeatures = 10,     
                 scale = T )
plot_top_weights(MOFAobject.trained,
                 view = "RNA",
                 factor = 4,  
                 nfeatures = 10,     
                 scale = T )

plot_top_weights(MOFAobject.trained,
                 view = "RNA",
                 factor = 7,  
                 nfeatures = 10,     
                 scale = T )

plot_top_weights(MOFAobject.trained,
                 view = "RNA",
                 factor = 8,  
                 nfeatures = 10,     
                 scale = T )


plot_data_scatter(MOFAobject.trained, 
                  view = "RNA",
                  factor = 3,  
                  features = 4,
                  sign = "negative",
                  add_lm = TRUE, lm_per_group = FALSE,
                  color_by = "scaled_phagocytosis"
                  
) + labs(y="RNA expression")
plot_data_scatter(MOFAobject.trained, 
                  view = "RNA",
                  factor = 3,  
                  features = 4,
                  sign = "positive",
                  add_lm = TRUE, lm_per_group = FALSE,
                  color_by = "scaled_phagocytosis"
                  
) + labs(y="RNA expression")

plot_data_scatter(MOFAobject.trained, 
                  view = "RNA",
                  factor = 3,  
                  features = 4,
                  sign = "negative",
                  add_lm = TRUE, lm_per_group = FALSE,
                  color_by = "log_ncells"
                  
) + labs(y="RNA expression")

plot_data_heatmap(MOFAobject.trained, 
                  view = "RNA",
                  factor = 7,  
                  features = 25,
                  denoise = TRUE,
                  cluster_rows = TRUE, cluster_cols = TRUE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)

# multi factor scatter plot

p <- plot_factors(MOFAobject.trained, 
                  factors = c(3,7), 
                  color_by = "scaled_phagocytosis",
                  shape_by = "treatment",
                  dot_size = 2.5,
                  show_missing = F, scale = TRUE
)

p <- p + 
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed")

print(p)

p <- plot_factors(MOFAobject.trained, 
                  factors = c(7,8), 
                  color_by = "scaled_phagocytosis",
                  shape_by = "treatment",
                  dot_size = 2.5,
                  show_missing = F, scale = TRUE
)

p <- p + 
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed")

print(p)


factors = plot_factors(MOFAobject.trained, 
                       factors = c(3,7), 
                       color_by = "scaled_phagocytosis",
                       dot_size = 2.5,
                       show_missing = F, scale = TRUE, return_data = TRUE
)

p_f3 = factors %>%
  dplyr::rename(scaled_phagocytosis = color_by,
                factor3=x,
                factor7=y) %>%
  ggplot(aes(y=scaled_phagocytosis,
             x = factor3)) +
  theme_bw() + 
  geom_point() + 
  geom_smooth(method = "lm") +
  ggpubr::stat_cor(label.y = -3.5)+ 
  ggpubr::stat_regline_equation(label.y = -4) +
  xlab("Factor 3") + 
  ylab("Scaled phagocytosis")

plot(p_f3)

p_f7 = factors %>%
  dplyr::rename(scaled_phagocytosis = color_by,
                factor3=x,
                factor7=y) %>%
  ggplot(aes(y=scaled_phagocytosis,
             x = factor7)) +
  theme_bw() + 
  geom_point() + 
  geom_smooth(method = "lm") +
  ggpubr::stat_cor(label.y = -3.5)+ 
  ggpubr::stat_regline_equation(label.y = -4) +
  xlab("Factor 7") + 
  ylab("Scaled phagocytosis")

plot(p_f7)

# PRS (microglia-specific) 

factors = plot_factors(MOFAobject.trained, 
                       factors = c(4,8), 
                       color_by = "untreated_full_PRS_microglia_PRS",
                       dot_size = 2.5,
                       show_missing = F, scale = TRUE, return_data = TRUE
)


p_f4 = factors %>%
  dplyr::rename(microglia_PRS = color_by,
                factor4=x,
                factor8=y) %>%
  # dplyr::group_by(microglia_PRS) %>%  
  # dplyr::summarise(factor11 = mean(factor11)) %>% # averaging factor11 for every identical PRS value, to clarify correlation
  # # this might be controversial
  # dplyr::ungroup() %>%
  
  ggplot(aes(y=microglia_PRS,
             x = factor4)) +
  theme_bw() + 
  geom_point() + 
  geom_smooth(method = "lm") +
  ggpubr::stat_cor(label.y = -3.5)+ 
  ggpubr::stat_regline_equation(label.y = -4) +
  xlab("Factor 4") + 
  ylab("microglia-specific PRS")

plot(p_f4)  # matches correlation table

p_f8 = factors %>%
  dplyr::rename(microglia_PRS = color_by,
                factor4=x,
                factor8=y) %>%
  # dplyr::group_by(microglia_PRS) %>%  
  # dplyr::summarise(factor14 = mean(factor14)) %>% # averaging factor11 for every identical PRS value, to clarify correlation
  # # this might be controversial
  # dplyr::ungroup() %>%
  
  ggplot(aes(y=microglia_PRS,
             x = factor8)) +
  theme_bw() + 
  geom_point() + 
  geom_smooth(method = "lm") +
  ggpubr::stat_cor(label.y = -3.5)+ 
  ggpubr::stat_regline_equation(label.y = -4) +
  xlab("Factor 8") + 
  ylab("microglia-specific PRS")

plot(p_f8)  # matches correlation table


# # # # #  gene set enrichment

library(fgsea)
library(msigdbr)


hallmark = msigdbr::msigdbr(species = "Homo sapiens", category = "H")
# immune = msigdbr::msigdbr(species = "Homo sapiens", subcategory  = "IMMUNESIGDB") #mostly generic useless
wikipath = msigdbr::msigdbr(species = "Homo sapiens", subcategory  = "CP:WIKIPATHWAYS")

hallmark_for_gsea = hallmark %>%
  group_by(gs_name) %>%
  summarise(genes = list(gene_symbol)) %>%
  deframe()
# immune_for_gsea = immune %>%
#   group_by(gs_name) %>%
#   summarise(genes = list(gene_symbol)) %>%
#   deframe()
wikipath_for_gsea = wikipath %>%
  group_by(gs_name) %>%
  summarise(genes = list(gene_symbol)) %>%
  deframe()


## jeremy genes
jeremy_ad_candidates = read_table("../../../resources/Jeremy_medrXiv_AD_candidate_genes.txt")
jeremy_set1_ad_candidates = read_csv("../../../resources/AD_PD_gene_sets_Andrew/Set1_Jeremy_candidates.csv")
jeremy_set3_pd_candidates = read_csv("../../../resources/AD_PD_gene_sets_Andrew/Set3_manually_curated.csv") %>%
  dplyr::filter(str_detect(disease,"PD")) %>%
  pull(gene_name)
# lead phagocytosis sam
sam_crispr = read_tsv("../../../CRISPR/OTAR2065_phagocytosis_CRISPR/data/2024_07_sam_screen_results.tsv") %>%
  dplyr::arrange(FDR) %>%
  slice_head(n = 30) %>%
  pull(id)

#### format my genes

pathways = list()
pathways[["jeremy_ad_candidates"]] = jeremy_ad_candidates$symbol
pathways[["set1_jeremy_ad_candidates"]] = jeremy_set1_ad_candidates$gene_sym # larger
pathways[["set3_jeremy_pd_candidates"]] = jeremy_set3_pd_candidates
pathways[["lead_phagocytosis_sam"]] = sam_crispr 


# working on weights

weights = list()
for (f in paste0("Factor",c(3:9))){
  weights[[f]] = pseudobulk_weights %>%
    dplyr::filter(factor == f) %>%
    dplyr::select(feature,value) %>%
    deframe() 
}


fgsea_results = list()
for(f in paste0("Factor",c(3:9))){
  fgsea_results[[paste("hallmark",f,sep = "_")]] = fgsea(hallmark_for_gsea, weights[[f]], maxSize=500,nproc = 1)
  # fgsea_results[[paste("immune",f,sep = "_")]] = fgsea(immune_for_gsea, weights[[f]], maxSize=500,nproc = 1)
  fgsea_results[[paste("wikipath",f,sep = "_")]] = fgsea(wikipath_for_gsea, weights[[f]], maxSize=500,nproc = 1)
  fgsea_results[[paste("custom",f,sep = "_")]] = fgsea(pathways, weights[[f]], maxSize=500,nproc = 1)
  # fgsea_results[[paste("custom",f,sep = "_")]] = fgsea(pathways, weights[[f]][setdiff(names( weights[[f]]),
  #                                                                                     c("APH1B","FCER1G","HS3ST1","INPP5D",
  #                                                                            "MS4A6A","PICALM","SCIMP","TREM1","TSPAN14"))], maxSize=500,nproc = 1)
  
  
}

sign_results = list()
for(f in  paste0("Factor",c(3:9))){
  sign_results[[paste("hallmark",f,sep = "_")]] = fgsea_results[[paste("hallmark",f,sep = "_")]] %>% as_tibble() %>% dplyr::filter(padj < 0.01) %>% dplyr::arrange(padj )
  # sign_results[[paste("immune",f,sep = "_")]] = fgsea_results[[paste("immune",f,sep = "_")]][order(pval)][padj < 0.0001] # not very useful
  sign_results[[paste("wikipath",f,sep = "_")]] = fgsea_results[[paste("wikipath",f,sep = "_")]] %>% as_tibble() %>% dplyr::filter(padj < 0.01) %>% dplyr::arrange(padj )
  sign_results[[paste("custom",f,sep = "_")]] =  fgsea_results[[paste("custom",f,sep = "_")]] %>% as_tibble() %>% dplyr::filter(padj < 0.05) %>% dplyr::arrange(padj )
  
}


# factor 3 - directly correlated with phagocytosis (cor = 0.41) and also with proliferation
f = "Factor3"
p3_fgsea_hallmark =  plotGseaTable(hallmark_for_gsea[sign_results[[paste("hallmark",f,sep = "_")]]$pathway], 
                                    weights[[f]], 
                                    sign_results[[paste("hallmark",f,sep = "_")]], 
                                    gseaParam = 0.5)
p3_fgsea_hallmark
# positive NES inflammatory response - positive correlation with phagocytosis - same as DGE analysis
# positive NES MYC targets - opposite to DGE analysis
p3_fgsea_custom =  plotGseaTable(pathways[sign_results[[paste("custom",f,sep = "_")]]$pathway], 
                                  weights[[f]], 
                                  sign_results[[paste("custom",f,sep = "_")]], 
                                  gseaParam = 0.5)
p3_fgsea_custom # no significant AD enrichments

p3_fgsea_wikipath =  plotGseaTable(wikipath_for_gsea[sign_results[[paste("wikipath",f,sep = "_")]]$pathway], 
                                    weights[[f]], 
                                    sign_results[[paste("wikipath",f,sep = "_")]], 
                                    gseaParam = 0.5)
p3_fgsea_wikipath # membrane receptors negative NES


# factor 7 - directly correlated with phagocytosis

f = "Factor7"
p7_fgsea_hallmark =  plotGseaTable(hallmark_for_gsea[sign_results[[paste("hallmark",f,sep = "_")]]$pathway], 
                                    weights[[f]], 
                                    sign_results[[paste("hallmark",f,sep = "_")]], 
                                    gseaParam = 0.5)
p7_fgsea_hallmark
# negative NES inflammatory response - inverse correlation with phagocytosis - opposite to factor 3
# positive NES MYC targets - opposite to DGE analysis

p7_fgsea_custom =  plotGseaTable(pathways[sign_results[[paste("custom",f,sep = "_")]]$pathway], 
                                 weights[[f]], 
                                 sign_results[[paste("custom",f,sep = "_")]], 
                                 gseaParam = 0.5)
sign_results[[paste("custom",f,sep = "_")]]$leadingEdge
p7_fgsea_custom #  significant AD enrichments (negative NES - inversely correlated with phagocytosis, directly correlated with higher risk microglial PRS)
# jeremy_ad_candidates      
# [1] "TREM1"      "HLA-DRA"    "MS4A6A"     "HS3ST1"     "HLA-DRB1"   "MS4A4A"     "APOE"       "TREM2"      "PILRA"      "AC090559.1" "INPP5D"     "SCIMP"     
# [13] "CLU"        "TSPAN14"    "PICALM"     "FCER1G"     "APH1B"   
# the problem with this is that there are many of these genes that show little or direct correlation with phagocytosis, so the results can be considered misleading
# most sign oxidative phosphorylation - negative NES, so enriched at negative weights, so directly correlated with high PRS
# other relationships also don't make sense - maybe because correlation is very low (-0.06)
# factor 14 has a higher correlation (-0.11) and shows inflammatory and oxidative phosphorilation in right direction
f = "Factor14"
p14_fgsea_hallmark =  plotGseaTable(hallmark_for_gsea[sign_results[[paste("hallmark",f,sep = "_")]]$pathway], 
                                    weights[[f]], 
                                    sign_results[[paste("hallmark",f,sep = "_")]], 
                                    gseaParam = 0.5)
p14_fgsea_hallmark

p14_fgsea_wikipath =  plotGseaTable(wikipath_for_gsea[sign_results[[paste("wikipath",f,sep = "_")]]$pathway], 
                                    weights[[f]], 
                                    sign_results[[paste("wikipath",f,sep = "_")]], 
                                    gseaParam = 0.5)
p14_fgsea_wikipath

p14_fgsea_custom =  plotGseaTable(pathways[sign_results[[paste("custom",f,sep = "_")]]$pathway], 
                                  weights[[f]], 
                                  sign_results[[paste("custom",f,sep = "_")]], 
                                  gseaParam = 0.5)
p14_fgsea_custom
# significant enrichment in AD candidate genes in high PRS 

# factor 5- directly correlated with phagocytosis (cor = 0.12, but also with log(n_cells) 0.27 so that's not great)
f = "Factor5"
p5_fgsea_hallmark =  plotGseaTable(hallmark_for_gsea[sign_results[[paste("hallmark",f,sep = "_")]]$pathway], 
                                   weights[[f]], 
                                   sign_results[[paste("hallmark",f,sep = "_")]], 
                                   gseaParam = 0.5)
p5_fgsea_hallmark
# negative oxidative phosphorylation - negative NES, so enriched at negative weights, so inversely correlated with high phagocytosis
# high at low phagocytosis

p5_fgsea_custom =  plotGseaTable(pathways[sign_results[[paste("custom",f,sep = "_")]]$pathway], 
                                 weights[[f]], 
                                 sign_results[[paste("custom",f,sep = "_")]], 
                                 gseaParam = 0.5)
p5_fgsea_custom # no significant AD enrichments

p5_fgsea_wikipath =  plotGseaTable(wikipath_for_gsea[sign_results[[paste("wikipath",f,sep = "_")]]$pathway], 
                                   weights[[f]], 
                                   sign_results[[paste("wikipath",f,sep = "_")]], 
                                   gseaParam = 0.5)
p5_fgsea_wikipath
# negative oxidative phosphorylation again

# factors 6 and 7 - inversely correlated (-0.24) and directly correlated (0.18) with phagocytosis

f = "Factor6"
p6_fgsea_hallmark =  plotGseaTable(hallmark_for_gsea[sign_results[[paste("hallmark",f,sep = "_")]]$pathway], 
                                   weights[[f]], 
                                   sign_results[[paste("hallmark",f,sep = "_")]], 
                                   gseaParam = 0.5)
p6_fgsea_hallmark
#  MYC targets positive NES - upregulated in low phagocytosis
# strong positive NES E2F targets (proliferation markers) - upregulated in low phagocytosis
# strong negative NES inflammatory response - upregulated in high phagocytosis
# 
p6_fgsea_custom =  plotGseaTable(pathways[sign_results[[paste("custom",f,sep = "_")]]$pathway], 
                                 weights[[f]], 
                                 sign_results[[paste("custom",f,sep = "_")]], 
                                 gseaParam = 0.5)
p6_fgsea_custom # negative NES - upregulated in high phagocytosis
sign_results[[paste("custom",f,sep = "_")]]$leadingEdge

#jeremy_ad_candidates      
# MS4A6A"     "HLA-DRA"    "TMEM163"    "HS3ST1"     "INPP5D"     "HLA-DRB1"   "TREM2"      "BIN1"       "PICALM"     "NCK2"      
#  "AC090559.1" "ARHGAP45"   "CR1"        "IKZF1"      "RIN3"       "SCIMP"      "CASTOR3"    "PLCG2"      "FCER1G"     "TREM1"     
#  "PTK2B"      "CD2AP"      "ABCA7"      "SORL1"      "MADD"       "TSPAN14"   

p6_fgsea_wikipath =  plotGseaTable(wikipath_for_gsea[sign_results[[paste("wikipath",f,sep = "_")]]$pathway], 
                                   weights[[f]], 
                                   sign_results[[paste("wikipath",f,sep = "_")]], 
                                   gseaParam = 0.5)
p6_fgsea_wikipath

# lots of cell cycle control genes - maybe because scaled proliferation also shows strong negative corr (-0.20)
#  checking factor 7 where there is no correlation with proliferation, but direct correlation with phagoctysosis

f = "Factor7"
p7_fgsea_hallmark =  plotGseaTable(hallmark_for_gsea[sign_results[[paste("hallmark",f,sep = "_")]]$pathway], 
                                   weights[[f]], 
                                   sign_results[[paste("hallmark",f,sep = "_")]], 
                                   gseaParam = 0.5)
p7_fgsea_hallmark
#  MYC targets positive NES - upregulated in high phagocytosis??? contradicts previous factor, so maybe don't focus on MYC at all
# strong positive oxidative phosphorylation - upregulated in high phagocytosis (contradicts factor 5)
#  positive NES inflammatory response - upregulated in high phagocytosis - agrees with previous factor
# 
p7_fgsea_custom =  plotGseaTable(pathways[sign_results[[paste("custom",f,sep = "_")]]$pathway], 
                                 weights[[f]], 
                                 sign_results[[paste("custom",f,sep = "_")]], 
                                 gseaParam = 0.5)
p7_fgsea_custom # no enrichments!

p7_fgsea_wikipath =  plotGseaTable(wikipath_for_gsea[sign_results[[paste("wikipath",f,sep = "_")]]$pathway], 
                                   weights[[f]], 
                                   sign_results[[paste("wikipath",f,sep = "_")]], 
                                   gseaParam = 0.5)
p7_fgsea_wikipath

# again not very useful, agrees with oxidative phosphorilation


# factors 7 and 8- inversely correlated (-0.15) and directly correlated (0.25) with migration


# 8 is strongly inversely correlated (-0.59) with sex, so better leave it alone


f = "Factor8"
p8_fgsea_hallmark =  plotGseaTable(hallmark_for_gsea[sign_results[[paste("hallmark",f,sep = "_")]]$pathway], 
                                   weights[[f]], 
                                   sign_results[[paste("hallmark",f,sep = "_")]], 
                                   gseaParam = 0.5)
p8_fgsea_hallmark

p8_fgsea_custom =  plotGseaTable(pathways[sign_results[[paste("custom",f,sep = "_")]]$pathway], 
                                 weights[[f]], 
                                 sign_results[[paste("custom",f,sep = "_")]], 
                                 gseaParam = 0.5)
p8_fgsea_custom # no enrichments!

p8_fgsea_wikipath =  plotGseaTable(wikipath_for_gsea[sign_results[[paste("wikipath",f,sep = "_")]]$pathway], 
                                   weights[[f]], 
                                   sign_results[[paste("wikipath",f,sep = "_")]], 
                                   gseaParam = 0.5)
p8_fgsea_wikipath








# factor 1 (-0.87 cor with ncells) out of curiosity

weights = list()
for (f in paste0("Factor",c(1))){
  weights[[f]] = pseudobulk_weights %>%
    dplyr::filter(factor == f) %>%
    dplyr::select(feature,value) %>%
    deframe() 
}

for(f in paste0("Factor",c(1))){
  fgsea_results[[paste("hallmark",f,sep = "_")]] = fgsea(hallmark_for_gsea, weights[[f]], maxSize=500,nproc = 1)
  # fgsea_results[[paste("immune",f,sep = "_")]] = fgsea(immune_for_gsea, weights[[f]], maxSize=500,nproc = 1)
  fgsea_results[[paste("wikipath",f,sep = "_")]] = fgsea(wikipath_for_gsea, weights[[f]], maxSize=500,nproc = 1)
  fgsea_results[[paste("custom",f,sep = "_")]] = fgsea(pathways, weights[[f]], maxSize=500,nproc = 1)
  
}

sign_results = list()
for(f in  paste0("Factor",c(1))){
  sign_results[[paste("hallmark",f,sep = "_")]] = fgsea_results[[paste("hallmark",f,sep = "_")]] %>% as_tibble() %>% dplyr::filter(padj < 0.01) %>% dplyr::arrange(padj )
  # sign_results[[paste("immune",f,sep = "_")]] = fgsea_results[[paste("immune",f,sep = "_")]][order(pval)][padj < 0.0001] # not very useful
  sign_results[[paste("wikipath",f,sep = "_")]] = fgsea_results[[paste("wikipath",f,sep = "_")]] %>% as_tibble() %>% dplyr::filter(padj < 0.01) %>% dplyr::arrange(padj )
  sign_results[[paste("custom",f,sep = "_")]] =  fgsea_results[[paste("custom",f,sep = "_")]] %>% as_tibble() %>% dplyr::filter(padj < 0.05) %>% dplyr::arrange(padj )
  
}

f = "Factor1"
p1_fgsea_hallmark =  plotGseaTable(hallmark_for_gsea[sign_results[[paste("hallmark",f,sep = "_")]]$pathway], 
                                   weights[[f]], 
                                   sign_results[[paste("hallmark",f,sep = "_")]], 
                                   gseaParam = 0.5)
p1_fgsea_hallmark

p1_fgsea_custom =  plotGseaTable(pathways[sign_results[[paste("custom",f,sep = "_")]]$pathway], 
                                 weights[[f]], 
                                 sign_results[[paste("custom",f,sep = "_")]], 
                                 gseaParam = 0.5)
p1_fgsea_custom # sam's results have negative NES, so are enriched for highly expressed genes / genes in cells that are more abundant


check_gene %>%
  dplyr::group_by(factor) %>%
  dplyr::filter(gene %in%   sign_results[[paste("custom",f,sep = "_")]]$leadingEdge[[1]]) %>%
  dplyr::ungroup() %>%
  dplyr::filter(factor %in% paste0("Factor",1:10)) %>%
  ggplot(aes(x = factor,y = quantile, label = gene)) + 
  geom_point() + 
  geom_text_repel() +
  theme_bw() + 
  ggtitle("Sam's leading edge in factor 1")

full_sam_crispr = read_tsv("../../../CRISPR/OTAR2065_phagocytosis_CRISPR/data/2024_07_sam_screen_results.tsv") %>%
  dplyr::arrange(FDR)

pdf(paste0(outputDir,"CRISPRko_sam_expression_effects.pdf"),width = 5,height = 4)
check_gene %>%
  dplyr::select(gene,quantile) %>%
  dplyr::distinct() %>%
  dplyr::filter(gene %in%   full_sam_crispr$id) %>%
  dplyr::left_join(full_sam_crispr,by = join_by("gene" == "id")) %>%
  dplyr::arrange(desc(`abs(t)`)) %>%
  ggplot(aes(x = quantile,y =`abs(t)`)) + 
  geom_point() + 
  geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +  # quadratic +
  xlab("Expression percentile") + 
  theme_bw()
dev.off()
# there doesn't seem to be a clear relationship between expression quantile and significance
# but most significant genes by t stat tend to be above 75% percentile



p1_fgsea_wikipath =  plotGseaTable(wikipath_for_gsea[sign_results[[paste("wikipath",f,sep = "_")]]$pathway], 
                                   weights[[f]], 
                                   sign_results[[paste("wikipath",f,sep = "_")]], 
                                   gseaParam = 0.5)
p1_fgsea_wikipath


# check these pathways among strongest DEGs
# pick best factors and plot enrichments with heatmap, highlight some good genes eg candidate genes, TREM1 etc

