
library(tidyverse)
library(mashr)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(grid)
library(sessioninfo)

#-------------------------------------------------------------------------------
#                               3.1 Run mash
#-------------------------------------------------------------------------------

## Set working dir
setwd("/lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_Daianna/")

## Define input, output, and plot dirs
inputdir = paste(getwd(), "input_data", "03_Mash_analysis", "01_run_mash", sep = "/")
outdir = paste(getwd(), "output_data", "03_Mash_analysis", "01_run_mash", sep = "/")
plotdir = paste(getwd(), "plots", "03_Mash_analysis", "01_run_mash", sep = "/")
dir.create(inputdir, recursive = T)
dir.create(outdir, recursive = T)
dir.create(plotdir, recursive = T)

## Input dir
input_dir0 = paste(getwd(), "output_data", "03_Mash_analysis", "00_prepare_input_data", sep = "/")

## Input data
load(paste0(input_dir0, "/B_hat_strong_set.Rdata"), verbose = T)
load(paste0(input_dir0, "/S_strong_set.Rdata"), verbose = T)
load(paste0(input_dir0, "/B_random_set.Rdata"), verbose = T)
load(paste0(input_dir0, "/S_random_set.Rdata"), verbose = T)


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
##                      3.1.1 Set up covariance matrices
##- - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

####################  3.1.1.1 Canonical covariance matrices  ###################

## B_hat with 32K random eQTLs
B_hat_random <- B_random_set[, -c(1:3)] %>% as.matrix()
## S_hat with 32K random eQTLs
S_hat_random <- S_random_set[, -c(1:3)] %>% as.matrix()

## Estimate corr structure in the null tests
mash_data_random_temp = mash_set_data(B_hat_random, S_hat_random)
Vhat = estimate_null_correlation_simple(mash_data_random_temp)

## Mash data object
mash_data_random = mash_set_data(B_hat_random, S_hat_random, V = Vhat)

## Generate canonical cov matrices
U.c = cov_canonical(mash_data_random)   
save(U.c, file = paste0(outdir, "/U.c.Rdata"))
names(U.c)
# [1] "identity"         "effect_IFN"       "effect_LPS"       "effect_untreated" "effect_IFNG_6"    "effect_IFNG_24"   "effect_sLPS_6"    "effect_sLPS_24"  
# [9] "effect_Ctrl_6"    "effect_Ctrl_24"   "equal_effects"    "simple_het_1"     "simple_het_2"     "simple_het_3"   

col_fun = colorRamp2(c(0, 0.5, 1), c("white", "darksalmon", "red"))
labels = c(paste("Microglia", c("IFN", "LPS", "untreated")), paste("Macrophage", c("IFN 6hr", "IFN 24hr", 
                                                                                   "sLPS 6hr", "sLPS 24hr", 
                                                                                   "Ctrl 6hr", "Ctrl 24hr")))
## Heatmaps
h <- list()
i = 1
for(c in U.c){
  h[[i]] <-  grid.grabExpr(draw(ComplexHeatmap::Heatmap(c, cluster_rows = F, 
                               cluster_columns = F, 
                               col = col_fun, 
                               border = T, 
                               name = "Cov", 
                               column_title = "Conditions", 
                               row_title = "Conditions", 
                               row_labels = labels, row_names_gp = gpar(fontsize = 9))))
  i = i +1
}

plot_grid(plotlist = h, align = "hv", nrow = 4)
ggsave(filename = paste0(plotdir, "/Canonical_cov_matrices.pdf"), height = 12, width = 20)


###################  3.1.1.2 Data-driven covariance matrices  ##################

## B_hat with 16K strong eQTLs
B_hat_strong <- as.matrix(B_hat_strong_set[, -c(1:3)])
## S_hat with 16K strong eQTLs
S_hat_strong <- as.matrix(S_strong_set[, -c(1:3)])
mash_data_strong = mash_set_data(B_hat_strong, S_hat_strong, V = Vhat)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
##                3.1.1.2.1 PCA for initial data-driven matrices
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
## SVD PCs
U.pca = cov_pca(data = mash_data_strong, npc = ncol(B_hat_strong), subset = NULL)
save(U.pca, file = paste0(outdir, "/U.pca.Rdata"))
names(U.pca)
# [1] "PCA_1" "PCA_2" "PCA_3" "PCA_4" "PCA_5" "PCA_6" "PCA_7" "PCA_8" "PCA_9" "tPCA" 

h <- list()
PCs <- gsub("_", "", names(U.pca))
i = 1
for(c in U.pca){
  h[[i]] <-  grid.grabExpr(draw(ComplexHeatmap::Heatmap(c, cluster_rows = F, 
                                                        cluster_columns = F, 
                                                        col = col_fun, 
                                                        border = T, 
                                                        name = "Cov", 
                                                        column_title = PCs[i], 
                                                        row_title = "Conditions", 
                                                        row_labels = labels, column_labels = labels,
                                                        row_names_gp = gpar(fontsize = 9),
                                                        column_names_gp = gpar(fontsize = 9))))
  i = i +1
}

plot_grid(plotlist = h, align = "hv", nrow = 2)
ggsave(filename = paste0(plotdir, "/PCA_cov_matrices.pdf"), height = 9.5, width = 26)


## PCA
pca <- prcomp(t(B_hat_strong), center = T)
## % of the variance explained by each PC
pca_vars <- signif(pca$sdev^2/sum(pca$sdev^2)*100, digits = 3)

pca_vars_labs<- paste0(
  "PC", seq(along = pca_vars), "\n", "(",
  pca_vars, "%)")
names(pca_vars_labs) <- paste0("PC", seq(along = pca_vars))
  
pca_data <- melt(pca$x)
colnames(pca_data) <- c("Condition", "PC", "Z")

labels = c(paste("Microglia", c("IFN", "LPS", "untreated")), paste("Macrophage", c("IFN 6hr", "IFN 24hr", 
                                                                                   "sLPS 6hr", "sLPS 24hr", 
                                                                                   "Ctrl 6hr", "Ctrl 24hr")))
names(labels) <- colnames(B_hat_strong)
pca_data$Condition <- sapply(pca_data$Condition, function(c){labels[c]})
pca_data$Condition <- factor(pca_data$Condition, levels = labels)
pca_data$Cells <- sapply(as.vector(pca_data$Condition), function(c){strsplit(c, " ")[[1]][1]})
  
pca_data$PC_lab <- sapply(pca_data$PC, function(pc){pca_vars_labs[pc]})
  
colors <- c("Microglia IFN" = "yellowgreen",
            "Microglia LPS"  = "darkorange1",
            "Microglia untreated" = "gray",
            "Macrophage IFN 6hr" = "yellowgreen",
            "Macrophage IFN 24hr" = "yellowgreen", 
            "Macrophage sLPS 6hr" = "darkorange1",
            "Macrophage sLPS 24hr" = "darkorange1",
            "Macrophage Ctrl 6hr" = "gray",
            "Macrophage Ctrl 24hr" = "gray")

alphas <- c(1,1,1, 0.4, 1, 0.4, 1, 0.4, 1)
shapes <- c("Microglia" = 8, "Macrophage" = 15)

plot <- ggplot(data = pca_data,
               aes(x = PC_lab, y = Z,
                   color = Condition,
                   alpha = Condition)) +
  theme_classic() +
  geom_boxplot(aes(x = PC_lab, y = Z), 
               width= 0.5, color='gray40', 
               linewidth = 0.4, alpha = 0.7, outliers = F) +
  geom_jitter(aes(shape = Cells), size = 2, width = 0.1, stroke = 1) +
  scale_color_manual(values = colors) +
  scale_alpha_manual(values = alphas) +
  scale_shape_manual(values = shapes) +
  labs(x = "PC (% Var Expl)")+
  theme(axis.title = element_text(size = (10)),
        axis.text = element_text(size = (8)),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10))

ggsave(filename = paste0(plotdir, "/PC_boxplots.pdf"), height = 5, width = 10)


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
##              3.1.1.2.2 More refined data-driven matrices with ED  
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
U.ed = cov_ed(mash_data_strong, U.pca, subset = NULL)
save(U.ed, file = paste0(outdir, "/U.ed.Rdata"))
names(U.ed)
# [1] "ED_PCA_1" "ED_PCA_2" "ED_PCA_3" "ED_PCA_4" "ED_PCA_5" "ED_PCA_6" "ED_PCA_7" "ED_PCA_8" "ED_PCA_9" "ED_tPCA" 

h <- list()
ED_PCs <- gsub("_", "", gsub("ED_", "ED-", names(U.ed)))
i = 1
for(c in U.ed){
  h[[i]] <-  grid.grabExpr(draw(ComplexHeatmap::Heatmap(c, cluster_rows = F, 
                                                        cluster_columns = F, 
                                                        col = col_fun, 
                                                        border = T, 
                                                        name = "Cov", 
                                                        column_title = ED_PCs[i], 
                                                        row_title = "Conditions", 
                                                        row_labels = labels, column_labels = labels,
                                                        row_names_gp = gpar(fontsize = 9),
                                                        column_names_gp = gpar(fontsize = 9))))
  i = i +1
}

plot_grid(plotlist = h, align = "hv", nrow = 2)
ggsave(filename = paste0(plotdir, "/ED-PCA_cov_matrices.pdf"), height = 9.5, width = 26)



## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
##                     3.1.2 Fit the mixture MVN model
##- - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#  Important:                     !!!                              *
#  The model fit must be carried out with the large random         *
#  subset of tests for mash to learn and correct for data          *
#  sparsity (not the strong set!!).                                *
## * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

## Fit model using canonical and data-driven cov matrices
m = mash(mash_data_random, c(U.c,U.ed))
save(m, file = paste0(outdir, "/mash_results.Rdata"))

## Estimated mixture proportions for:
## - Cov matrices:
mix_props_cov <- get_estimated_pi(m, dimension = "cov") 
#           null         identity       effect_IFN       effect_LPS effect_untreated    effect_IFNG_6   effect_IFNG_24    effect_sLPS_6 
#   1.804109e-01     0.000000e+00     0.000000e+00     4.420193e-05     0.000000e+00     1.582047e-03     3.098595e-03     8.748427e-03 
# effect_sLPS_24    effect_Ctrl_6   effect_Ctrl_24    equal_effects     simple_het_1     simple_het_2     simple_het_3         ED_PCA_1 
#   3.186076e-03     1.263727e-03     1.785857e-03     5.734390e-03     0.000000e+00     0.000000e+00     0.000000e+00     3.464728e-01 
#       ED_PCA_2         ED_PCA_3         ED_PCA_4         ED_PCA_5         ED_PCA_6         ED_PCA_7         ED_PCA_8         ED_PCA_9 
#   9.172721e-03     9.708697e-04     2.845344e-03     0.000000e+00     0.000000e+00     0.000000e+00     0.000000e+00     2.702814e-05 
#        ED_tPCA 
#   4.346570e-01 

# - Scaling factors
mix_props_grid <- get_estimated_pi(m, dimension = "grid") 
#         null            1            2            3            4            5            6            7            8            9           10 
# 1.804109e-01 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 
#           11           12           13           14           15           16           17           18           19           20           21 
# 6.669764e-03 2.346819e-01 9.902395e-02 7.825014e-04 2.119935e-01 9.531034e-02 1.101759e-01 3.916977e-02 1.890491e-02 2.545462e-03 2.919041e-04 
#           22           23           24 
# 3.924226e-05 0.000000e+00 0.000000e+00 

# - Cov matrix - scaling factor combinations
mix_props_all <- get_estimated_pi(m, dimension = "all") 
length(mix_props_all)
# 577
head(mix_props_all)
#      null         identity.1       effect_IFN.1       effect_LPS.1 effect_untreated.1    effect_IFNG_6.1 
# 0.1804109          0.0000000          0.0000000          0.0000000          0.0000000          0.0000000 
tail(mix_props_all)
# ED_PCA_5.24 ED_PCA_6.24 ED_PCA_7.24 ED_PCA_8.24 ED_PCA_9.24  ED_tPCA.24 
#           0           0           0           0           0           0 

pdf(file = paste0(plotdir, "/mixture_props.pdf"), height = 6, width = 17)
par(mfrow=c(1, 3))

barplot(mix_props_cov, las = 2, cex.names = 0.55)
barplot(mix_props_grid, las = 2, cex.names = 0.55)
barplot(mix_props_all[which(!mix_props_all == 0)], las = 2, cex.names = 0.55)

dev.off()



## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
##                     3.1.3 Extract posterior summaries
##- - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## Compute posterior summaries with previous fit
m2 = mash(mash_data_strong, g = get_fitted_g(m), fixg = TRUE)
save(m2, file = paste0(outdir, "/mash_results_strong_set.Rdata"))

## Same mixture proportions from previous fit? Yes!
mix_props_cov2 <- get_estimated_pi(m2, dimension = "cov") 
mix_props_grid2 <- get_estimated_pi(m2, dimension = "grid") 
mix_props_all2 <- get_estimated_pi(m2, dimension = "all") 

identical(mix_props_cov, mix_props_cov2)
# [1] TRUE
identical(mix_props_grid, mix_props_grid2)
# [1] TRUE
identical(mix_props_all, mix_props_all2)
# [1] TRUE

## Posterior means
posterior_means <- get_pm(m2)
head(posterior_means)
#       effect_IFN  effect_LPS effect_untreated effect_IFNG_6 effect_IFNG_24 effect_sLPS_6 effect_sLPS_24 effect_Ctrl_6 effect_Ctrl_24
# 379  -0.01814949 -0.12455185      -0.02194919   -0.01383641    -0.01231029  -0.016644767   -0.007813214    0.01069049     0.04483669
# 651  -0.11447547 -0.14306458      -0.09183864   -0.15147391    -0.22453657  -0.096388291   -0.213507833   -0.16630274    -0.19940584
# 843  -0.12535929 -0.10170052      -0.07971998   -0.08918569    -0.16708467   0.001663174   -0.154137244   -0.07669196    -0.11421549
# 1828  0.17841173  0.15542892       0.16193262    0.55230209     0.42613783   0.422251504    0.405096527    0.49518213     0.44420754
# 1835  0.13016973  0.08970409       0.11606083    0.48349139     0.35335833   0.324751281    0.327439371    0.40097869     0.36127257
# 2312  0.18414905  0.24268027       0.28918363    0.25895859     0.31203749   0.245432275    0.290233194    0.27335875     0.32438192

## Posterior std deviations
posterior_sd <- get_psd(m2)
head(posterior_sd)
#      effect_IFN effect_LPS effect_untreated effect_IFNG_6 effect_IFNG_24 effect_sLPS_6 effect_sLPS_24 effect_Ctrl_6 effect_Ctrl_24
# 379  0.03101757 0.07203895       0.04585859    0.07655034     0.07627868    0.07800219     0.07606942    0.07544203     0.08687575
# 651  0.03111355 0.03718625       0.03447333    0.07240884     0.07823108    0.07736490     0.07505317    0.06915014     0.07404952
# 843  0.03863761 0.04252295       0.03900228    0.08301305     0.08768072    0.09844505     0.08699975    0.08061219     0.08222545
# 1828 0.05151338 0.05938871       0.05598434    0.11122057     0.09708031    0.10290131     0.09627831    0.10161070     0.09818462
# 1835 0.04867984 0.05774119       0.05389831    0.09367194     0.08075995    0.08362825     0.08088358    0.08268894     0.08146761
# 2312 0.05148597 0.06188698       0.06120811    0.08403689     0.08399301    0.08533290     0.08218749    0.08013631     0.08414305

## Significant eQTLs (signif in at least one condition)
signif_mash_eQTLs_all_conditions <- get_significant_results(m2) 
length(signif_mash_eQTLs_all_conditions)
# [1] 15966

## Sharing of signif effects by sign and magnitude
sharing_effects <- get_pairwise_sharing(m2)
sharing_effects
#                  effect_IFN effect_LPS effect_untreated effect_IFNG_6 effect_IFNG_24 effect_sLPS_6 effect_sLPS_24 effect_Ctrl_6 effect_Ctrl_24
# effect_IFN        1.0000000  0.8131838        0.8358323     0.4907185      0.4879008     0.5201409      0.4931270     0.4912257      0.4711962
# effect_LPS        0.8131838  1.0000000        0.8160840     0.5076780      0.4993712     0.5574593      0.5504439     0.5223232      0.5116435
# effect_untreated  0.8358323  0.8160840        1.0000000     0.5523260      0.5555850     0.5675855      0.5682377     0.5827877      0.5699565
# effect_IFNG_6     0.4907185  0.5076780        0.5523260     1.0000000      0.9660832     0.8902103      0.9009222     0.9260917      0.8973607
# effect_IFNG_24    0.4879008  0.4993712        0.5555850     0.9660832      1.0000000     0.8484505      0.9161796     0.9296914      0.9191386
# effect_sLPS_6     0.5201409  0.5574593        0.5675855     0.8902103      0.8484505     1.0000000      0.8747395     0.8532267      0.8117964
# effect_sLPS_24    0.4931270  0.5504439        0.5682377     0.9009222      0.9161796     0.8747395      1.0000000     0.9416859      0.9368459
# effect_Ctrl_6     0.4912257  0.5223232        0.5827877     0.9260917      0.9296914     0.8532267      0.9416859     1.0000000      0.9824883
# effect_Ctrl_24    0.4711962  0.5116435        0.5699565     0.8973607      0.9191386     0.8117964      0.9368459     0.9824883      1.0000000

pdf(file =paste0(plotdir, "/Pairwise_effect_sharing.pdf"), height = 7, width = 8)
col_fun = colorRamp2(seq(from = 0 , to  = 1, length = 5), 
                      c("mintcream", "lightskyblue1", "lightskyblue2", "lightskyblue3", "lightskyblue4"))
h <- ComplexHeatmap::Heatmap(sharing_effects, 
                         cluster_rows = F, 
                         cluster_columns = F, 
                         col = col_fun, 
                         border = T, 
                         name = "Sharing", 
                         column_title = "Conditions", 
                         row_labels = labels, 
                         column_labels = labels, 
                         row_names_gp = gpar(fontsize = 9), 
                         column_names_gp = gpar(fontsize = 9),
                         cell_fun = function(j, i, x, y, width, height, fill) {
                           grid.text(signif(sharing_effects[i, j], digits = 3), x, y, gp = gpar(fontsize = 9)) } )

h
dev.off()

## Save results
mash_results <- list(posterior_means, posterior_sd, signif_mash_eQTLs_all_conditions, sharing_effects)
save(mash_results, file = paste0(outdir, "/mash_results_list.Rdata"))







## Reproducibility info
options(width = 120)
session_info()
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.1 (2023-06-16)
# os       Ubuntu 22.04.4 LTS
# system   x86_64, linux-gnu
# ui       RStudio
# language (EN)
# collate  en_GB.UTF-8
# ctype    en_GB.UTF-8
# tz       Europe/Belfast
# date     2024-11-28
# rstudio  2024.04.0+735 Chocolate Cosmos (server)
# pandoc   3.1.12.3 @ /opt/view/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [1] CRAN (R 4.3.1)
# ashr                 * 2.2-63    2023-08-21 [1] CRAN (R 4.3.1)
# assertthat             0.2.1     2019-03-21 [1] CRAN (R 4.3.1)
# Biobase              * 2.62.0    2023-10-24 [1] Bioconductor
# BiocGenerics         * 0.48.1    2023-11-01 [1] Bioconductor
# bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.3.1)
# Cairo                  1.6-2     2023-11-28 [1] CRAN (R 4.3.1)
# circlize             * 0.4.16    2024-02-20 [1] CRAN (R 4.3.1)
# cli                    3.6.2     2023-12-11 [1] CRAN (R 4.3.1)
# clue                   0.3-65    2023-09-23 [1] CRAN (R 4.3.1)
# cluster                2.1.6     2023-12-01 [1] CRAN (R 4.3.1)
# codetools              0.2-20    2024-03-31 [1] CRAN (R 4.3.1)
# colorspace             2.1-0     2023-01-23 [1] CRAN (R 4.3.1)
# ComplexHeatmap       * 2.18.0    2023-10-24 [1] Bioconductor
# cowplot              * 1.1.3     2024-01-22 [1] CRAN (R 4.3.1)
# crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.1)
# DelayedArray           0.28.0    2023-10-24 [1] Bioconductor
# digest                 0.6.35    2024-03-11 [1] CRAN (R 4.3.1)
# doParallel             1.0.17    2022-02-07 [1] CRAN (R 4.3.1)
# dplyr                * 1.1.4     2023-11-17 [1] CRAN (R 4.3.1)
# edgeR                  4.0.16    2024-02-18 [1] Bioconductor 3.18 (R 4.3.1)
# fansi                  1.0.6     2023-12-08 [1] CRAN (R 4.3.1)
# farver                 2.1.1     2022-07-06 [1] CRAN (R 4.3.1)
# forcats              * 1.0.0     2023-01-29 [1] CRAN (R 4.3.1)
# foreach                1.5.2     2022-02-02 [1] CRAN (R 4.3.1)
# generics               0.1.3     2022-07-05 [1] CRAN (R 4.3.1)
# GenomeInfoDb         * 1.38.8    2024-03-15 [1] Bioconductor 3.18 (R 4.3.1)
# GenomeInfoDbData       1.2.11    2024-11-06 [1] Bioconductor
# GenomicRanges        * 1.54.1    2023-10-29 [1] Bioconductor
# GetoptLong             1.0.5     2020-12-15 [1] CRAN (R 4.3.1)
# ggplot2              * 3.5.1     2024-04-23 [1] CRAN (R 4.3.1)
# ggrepel                0.9.5     2024-01-10 [1] CRAN (R 4.3.1)
# GlobalOptions          0.1.2     2020-06-10 [1] CRAN (R 4.3.1)
# glue                   1.7.0     2024-01-09 [1] CRAN (R 4.3.1)
# gtable                 0.3.4     2023-08-21 [1] CRAN (R 4.3.1)
# hms                    1.1.3     2023-03-21 [1] CRAN (R 4.3.1)
# invgamma               1.1       2017-05-07 [1] CRAN (R 4.3.1)
# IRanges              * 2.36.0    2023-10-24 [1] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [1] CRAN (R 4.3.1)
# iterators              1.0.14    2022-02-05 [1] CRAN (R 4.3.1)
# labeling               0.4.3     2023-08-29 [1] CRAN (R 4.3.1)
# lattice                0.22-6    2024-03-20 [1] CRAN (R 4.3.1)
# lifecycle              1.0.4     2023-11-07 [1] CRAN (R 4.3.1)
# limma                * 3.58.1    2023-10-31 [1] Bioconductor
# locfit                 1.5-9.9   2024-03-01 [1] CRAN (R 4.3.1)
# lubridate            * 1.9.3     2023-09-27 [1] CRAN (R 4.3.1)
# magick                 2.8.3     2024-02-18 [1] CRAN (R 4.3.1)
# magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.3.1)
# mashr                * 0.2.79    2023-10-18 [1] CRAN (R 4.3.1)
# Matrix                 1.6-5     2024-01-11 [1] CRAN (R 4.3.1)
# MatrixGenerics       * 1.14.0    2023-10-24 [1] Bioconductor
# matrixStats          * 1.2.0     2023-12-11 [1] CRAN (R 4.3.1)
# mixsqp                 0.3-54    2023-12-20 [1] CRAN (R 4.3.1)
# munsell                0.5.1     2024-04-01 [1] CRAN (R 4.3.1)
# mvtnorm                1.2-4     2023-11-27 [1] CRAN (R 4.3.1)
# pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.3.1)
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.1)
# plyr                   1.8.9     2023-10-02 [1] CRAN (R 4.3.1)
# png                    0.1-8     2022-11-29 [1] CRAN (R 4.3.1)
# purrr                * 1.0.2     2023-08-10 [1] CRAN (R 4.3.1)
# R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.1)
# ragg                   1.3.0     2024-03-13 [1] CRAN (R 4.3.1)
# RColorBrewer           1.1-3     2022-04-03 [1] CRAN (R 4.3.1)
# Rcpp                   1.0.12    2024-01-09 [1] CRAN (R 4.3.1)
# RCurl                  1.98-1.14 2024-01-09 [1] CRAN (R 4.3.1)
# readr                * 2.1.5     2024-01-10 [1] CRAN (R 4.3.1)
# reshape2             * 1.4.4     2020-04-09 [1] CRAN (R 4.3.1)
# rjson                  0.2.21    2022-01-09 [1] CRAN (R 4.3.1)
# rlang                * 1.1.3     2024-01-10 [1] CRAN (R 4.3.1)
# rmeta                  3.0       2018-03-20 [1] CRAN (R 4.3.1)
# rstudioapi             0.16.0    2024-03-24 [1] CRAN (R 4.3.1)
# S4Arrays               1.2.1     2024-03-04 [1] Bioconductor 3.18 (R 4.3.1)
# S4Vectors            * 0.40.2    2023-11-23 [1] Bioconductor 3.18 (R 4.3.1)
# scales                 1.3.0     2023-11-28 [1] CRAN (R 4.3.1)
# sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.3.1)
# shape                  1.4.6.1   2024-02-23 [1] CRAN (R 4.3.1)
# SparseArray            1.2.4     2024-02-11 [1] Bioconductor 3.18 (R 4.3.1)
# SQUAREM                2021.1    2021-01-13 [1] CRAN (R 4.3.1)
# statmod                1.5.0     2023-01-06 [1] CRAN (R 4.3.1)
# stringi                1.8.3     2023-12-11 [1] CRAN (R 4.3.1)
# stringr              * 1.5.1     2023-11-14 [1] CRAN (R 4.3.1)
# SummarizedExperiment * 1.32.0    2023-10-24 [1] Bioconductor
# systemfonts            1.0.6     2024-03-07 [1] CRAN (R 4.3.1)
# textshaping            0.3.7     2023-10-09 [1] CRAN (R 4.3.1)
# tibble               * 3.2.1     2023-03-20 [1] CRAN (R 4.3.1)
# tidyr                * 1.3.1     2024-01-24 [1] CRAN (R 4.3.1)
# tidyselect             1.2.1     2024-03-11 [1] CRAN (R 4.3.1)
# tidyverse            * 2.0.0     2023-02-22 [1] CRAN (R 4.3.1)
# timechange             0.3.0     2024-01-18 [1] CRAN (R 4.3.1)
# truncnorm              1.0-9     2023-03-20 [1] CRAN (R 4.3.1)
# tzdb                   0.4.0     2023-05-12 [1] CRAN (R 4.3.1)
# utf8                   1.2.4     2023-10-22 [1] CRAN (R 4.3.1)
# vctrs                  0.6.5     2023-12-01 [1] CRAN (R 4.3.1)
# withr                  3.0.0     2024-01-16 [1] CRAN (R 4.3.1)
# XVector                0.42.0    2023-10-24 [1] Bioconductor
# zlibbioc               1.48.2    2024-03-13 [1] Bioconductor 3.18 (R 4.3.1)
# 
# [1] /opt/view/rlib/R/library
# [2] /software/hgi/installs/softpack/rstudio/.spack/linux-ubuntu22.04-x86_64_v3/gcc-11.4.0/r-4.3.1-bfwldrk76z6f52upk47zepliekn7ayqz/rlib/R/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────




