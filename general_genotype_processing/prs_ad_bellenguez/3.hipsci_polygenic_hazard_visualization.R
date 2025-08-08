#!/usr/bin/env Rscript
library(reshape2)
library(ggplot2)
library(plyr)
library(summarytools)
library(ggpubr)
library(ggExtra)

options(stringsAsFactors=F)

trait = "AD"

## modify for snakemake
input_dir = paste0("../../data/",trait,"/",trait,"_output")
output_dir = paste0(input_dir,"/plots")
hipsci_dir = "../../data/resources"


##
cell_donor_list = list()

for(names in c("cell_lines","donors")){

  cell_donor_list[[names]] = read.table(paste0(input_dir,"/",trait,"_polygenic_hazard.",names,".txt"), header = T)
}
## some cell lines from same donors have slightly different scores! Eg. eorc. Mutations in SNPs that are counted for polygenic risk?

lapply(cell_donor_list,nrow)

# gather everything in sampleID after "-"
cell_donor_list$cell_lines$Line = gsub("^.*?-","",cell_donor_list$cell_lines$sampleID)
## Binning per quartile
cell_donor_list$cell_lines$quartile = rep(NA, nrow(cell_donor_list$cell_lines))
cell_donor_list$cell_lines[cell_donor_list$cell_lines$overallHR_quantile <= 0.25,"quartile"] = "Q1"
cell_donor_list$cell_lines[cell_donor_list$cell_lines$overallHR_quantile <= 0.50 & cell_donor_list$cell_lines$overallHR_quantile > 0.25,"quartile"] = "Q2"
cell_donor_list$cell_lines[cell_donor_list$cell_lines$overallHR_quantile <= 0.75 & cell_donor_list$cell_lines$overallHR_quantile > 0.50,"quartile"] = "Q3"
cell_donor_list$cell_lines[cell_donor_list$cell_lines$overallHR_quantile > 0.75,"quartile"] = "Q4"


## Overall risk hazard per cell donor
ggplot(cell_donor_list$donors, aes(x=overallHR)) +
  geom_histogram(binwidth=1,color="black", fill="white") +
geom_vline(aes(xintercept=median(overallHR)),
           color="blue", linetype="dashed", size=1)


melted = melt(cell_donor_list$donors,id.vars = "donor",
              measure.vars = c("overallHR","apoeHR","polygenicHR"), variable.name = "Score")

p<-ggplot(melted, aes(x=value, fill=Score, color=Score)) +
  geom_histogram(binwidth=0.2,position="identity", alpha=0.5)
p

median <- ddply(melted, "Score", summarise, grp.median=median(value))
# Add mean lines
p+geom_vline(data=median, aes(xintercept=grp.median, color=Score),
             linetype="dashed") +
  theme_bw() + theme(legend.position="bottom")

ggsave(
  paste0(
    output_dir,
    "/hipsci.",
    trait,
    "_polygenic_hazard_histogram.donor.png"
  ),
  device = "png",
  height = 5,
  width = 5,
  units = "in",
  dpi = 400
)

## Boxplot

p <- ggplot(melted, aes(x=Score, y=value)) +
  geom_boxplot() +   theme_bw()

p


ggsave(
  paste0(
    output_dir,
    "/hipsci.",
    trait,
    "_polygenic_hazard_boxplot.donor.png"
  ),
  device = "png",
  height = 5,
  width = 5,
  units = "in",
  dpi = 400
)

## Not all cell lines are European
# Some are diseased
# Some have poor pluripotency

hipsci_lines =  read.delim(paste0(hipsci_dir,"/hipsci_org_all_lines.txt"), header = T, sep = "\t")

# Rename for merging
colnames(hipsci_lines)[1] = "sampleID"

### Neuroseq differentiation information
## Lines pooled in Neuroseq
NeuroSeq_Cell_line_list_forR = read.delim(paste0(hipsci_dir,"/NeuroSeq_Cell_line_list_forR.txt"), header = T, sep = "\t")
colnames(NeuroSeq_Cell_line_list_forR)[1] = "Line"
NeuroSeq_Cell_line_list_forR = NeuroSeq_Cell_line_list_forR[c(-2)]
colnames(NeuroSeq_Cell_line_list_forR)[2] = "Neuroseq_Pool.Status"

## Some failed Neuroseq differentiation
failed_lines =  read.delim(paste0(hipsci_dir,"/Failed_Lines_in_NeuroSeq_recode.txt"), header = T, sep = "\t")
failed_lines = failed_lines[c(-2)]
colnames(failed_lines)[2:5] = c("Neuroseq_Attempted.Pool" , "Neuroseq_Pool.if.repeat.successful", "Neuroseq_Reason.for.Fail",
                               "Neuroseq_Comment")
## differentiation efficiency classification by Julie

diff_efficiency =  read.delim(paste0(hipsci_dir,"/neuroseq_differentiation_efficiency.pool1_17_D52.tsv"), header = T, sep = "\t")
colnames(diff_efficiency) = c("sampleID","Neuroseq_Differentiation_efficiency")

###

## add this info to our lines with PRS
cell_donor_list$cell_lines = merge(cell_donor_list$cell_lines, hipsci_lines, by = "sampleID", all.x = T)

cell_donor_list$cell_lines = merge(cell_donor_list$cell_lines, failed_lines, by = "Line", all.x = T)
cell_donor_list$cell_lines = merge(cell_donor_list$cell_lines, NeuroSeq_Cell_line_list_forR, by = "Line", all.x = T)
cell_donor_list$cell_lines = merge(cell_donor_list$cell_lines, diff_efficiency, by = "sampleID", all.x = T)

## simplifying differentiation pools: differentiated? yes / FAIL / no / ?
cell_donor_list$cell_lines$`Neuroseq_Differentiated` = "no"
cell_donor_list$cell_lines[which(cell_donor_list$cell_lines$`Neuroseq_Pool.Status`=="FAIL"),"Neuroseq_Differentiated"] =   "FAIL"
cell_donor_list$cell_lines[which(cell_donor_list$cell_lines$`Neuroseq_Pool.Status`==""),"Neuroseq_Differentiated"] =   "?"
cell_donor_list$cell_lines[which(cell_donor_list$cell_lines$`Neuroseq_Pool.Status`!="" & cell_donor_list$cell_lines$`Neuroseq_Pool.Status`!="FAIL" & cell_donor_list$cell_lines$`Neuroseq_Pool.Status`!="NA"),"Neuroseq_Differentiated"] =   "yes"


freq(cell_donor_list$cell_lines$Cell.Type, cumul = F)
freq(cell_donor_list$cell_lines$Disease.Status, cumul = F)
# 502 that are Normal
# 1116 for which we have no information
# After removing repeated lines
non_disease_European = cell_donor_list$cell_lines[which(cell_donor_list$cell_lines$Disease.Status=="Normal"),]
length(unique(non_disease_European$sampleID))
# 496 that are Normal
length(unique(non_disease_European$donor))
# 330 donors

# Filter for feeder-free
non_disease_European = non_disease_European[non_disease_European$Culture == "Feeder-free",]
length(unique(non_disease_European$sampleID))
# 496 that are Normal
length(unique(non_disease_European$donor))
# 330 donors



table(cell_donor_list$cell_lines[which(cell_donor_list$cell_lines$Disease.Status=="Normal"),"Predicted.Population"])
 non_disease_European = cell_donor_list$cell_lines[which(cell_donor_list$cell_lines$Disease.Status=="Normal" & cell_donor_list$cell_lines$Predicted.Population=="European"),]
length(unique(non_disease_European$sampleID))
# of the normal, 470 are european
length(unique(non_disease_European$donor))
# 312 donors

freq(cell_donor_list$cell_lines$Sex, cumul = F)
freq(cell_donor_list$cell_lines$Ethnicity, cumul = F)
freq(cell_donor_list$cell_lines$Tissue.Provider, cumul = F)
freq(cell_donor_list$cell_lines$Age, cumul = F)
freq(cell_donor_list$cell_lines$Source.Material, cumul = F)
freq(cell_donor_list$cell_lines$Predicted.Population, cumul = F)
freq(cell_donor_list$cell_lines$Culture, cumul = F)
freq(cell_donor_list$cell_lines$Method.of.derivation, cumul = F)
freq(cell_donor_list$cell_lines$CNV.num.different.regions, cumul = F)
freq(cell_donor_list$cell_lines$Open.access.data, cumul = F) # 0 = open, 1 = restricted
freq(cell_donor_list$cell_lines$Assays.data.available, cumul = F) # 0 = open, 1 = restricted

freq(cell_donor_list$cell_lines$`Neuroseq_Pool.Status`, cumul = F)
# 297 lines differentiated in 21 pools (avg = 14 lines/pool), in which 51 failed differentiation (83% success rate)

freq(cell_donor_list$cell_lines$`Neuroseq_Differentiated`, cumul = F)

check = cell_donor_list$cell_lines

### Possible errors in differentiation scores ####################
## lines in pools without efficiency score
cell_donor_list$cell_lines[cell_donor_list$cell_lines$Neuroseq_Differentiated == "yes" & is.na(cell_donor_list$cell_lines$Neuroseq_Differentiation_efficiency),]

## lines not in pools with efficiency scores
cell_donor_list$cell_lines[cell_donor_list$cell_lines$Neuroseq_Differentiated == "no" & !is.na(cell_donor_list$cell_lines$Neuroseq_Differentiation_efficiency),]

non_disease_European = cell_donor_list$cell_lines[which(cell_donor_list$cell_lines$Disease.Status=="Normal" & cell_donor_list$cell_lines$Predicted.Population=="European"),]
length(unique(non_disease_European$sampleID))
# 470 non-diseased european

# PluriTest Assay compares the transcriptional profile of a sample to an extensive reference set of >450 cell and tissue types:
#   including 223 hESCs, 41 induced pluripotent stem cell lines, somatic cells, and tissues.
#  confirms pluripotency via two separate scores: Pluripotency and Novelty.
# The Pluripotency Score is an indication of how strongly a model-based pluripotency signature is expressed
# in the samples analyzed. The Novelty Score indicates the general model fit for a given sample.
# The Pluripotency Score informs the user on how strongly a model-based pluripotency signature is expressed in the samples analyzed.
# The Novelty Score indicates the general model fit for a given sample: if it is low, the analyzed sample is well represented in the
# current PluriTest data model, if the Novelty Score is high and the sample also has a high Pluripotency Score, further functional
# assays might be necessary determine the identity of the sample

# Pluripotency vs novelty scores vs overall HR quantile


scaterplot <- ggscatter(cell_donor_list$cell_lines, x = "Pluritest.novelty.score", y = "Pluritest.pluripotency.score",
                        size = 3, alpha = 0.6,  color = "overallHR_quantile") +
  scale_color_distiller( palette = "Spectral") +
  border()  + xlab("Novelty score") + ylab("Pluripotency score")

# Marginal density plot of x (top panel) and y (right panel)
xplot <- ggdensity(cell_donor_list$cell_lines, "Pluritest.novelty.score")
yplot <- ggdensity(cell_donor_list$cell_lines, "Pluritest.pluripotency.score")+
  rotate()
# Cleaning the plots
yplot <- yplot + clean_theme()
xplot <- xplot + clean_theme()
# Arranging the plot
ggarrange(xplot, NULL, scaterplot, yplot,
          ncol = 2, nrow = 2,  align = "hv",
          widths = c(2, 1), heights = c(1, 2),
          common.legend = TRUE)



ggsave(
  paste0(
    output_dir,
    "/hipsci.",
    trait,
    "_pluritest_vs_overallHRquantile.cell_lines.png"
  ),
  device = "png",
  height = 8,
  width = 8,
  units = "in",
  dpi = 400
)


# Pluripotency vs novelty scores vs overall HR quaRtile

color_manual = c("#0471A6","#06D6A0", "#F4D35E","#ED254E")
scaterplot <- ggscatter(cell_donor_list$cell_lines, x = "Pluritest.novelty.score", y = "Pluritest.pluripotency.score",
                        size = 3, alpha = 0.6,  color = "quartile") +
  scale_color_manual( values = color_manual) +
  border()  + xlab("Novelty score") + ylab("Pluripotency score")

# Marginal density plot of x (top panel) and y (right panel)
xplot <- ggdensity(cell_donor_list$cell_lines, "Pluritest.novelty.score", fill = "quartile", palette = color_manual)
yplot <- ggdensity(cell_donor_list$cell_lines, "Pluritest.pluripotency.score", fill = "quartile", palette = color_manual)+
  rotate()
# Cleaning the plots
yplot <- yplot + clean_theme()
xplot <- xplot + clean_theme()
# Arranging the plot
ggarrange(xplot, NULL, scaterplot, yplot,
          ncol = 2, nrow = 2,  align = "hv",
          widths = c(2, 1), heights = c(1, 2),
          common.legend = TRUE)



ggsave(
  paste0(
    output_dir,
    "/hipsci.",
    trait,
    "_pluritest_vs_overallHRquaRtile.cell_lines.png"
  ),
  device = "png",
  height = 8,
  width = 8,
  units = "in",
  dpi = 400
)
## Some failed in NeuroSeq diff

## are failed lines evenly distributed among pluripotency/novelty score?

fail_or_not_subset = cell_donor_list$cell_lines[!cell_donor_list$cell_lines$Neuroseq_Differentiated == "no",]

# Ask about possible label errors in differentiated cell lines
freq(fail_or_not_subset$Disease.Status, cumul = F)
freq(fail_or_not_subset$Predicted.Population, cumul = F)


scaterplot <- ggscatter(fail_or_not_subset, x = "Pluritest.novelty.score", y = "Pluritest.pluripotency.score",
                        size = 3, alpha = 0.6,  color = "Neuroseq_Differentiated") +
  scale_color_brewer( palette = "Accent") +
  border()  + xlab("Novelty score") + ylab("Pluripotency score")
# Marginal density plot of x (top panel) and y (right panel)
xplot <- ggdensity(fail_or_not_subset, "Pluritest.novelty.score", fill = "Neuroseq_Differentiated", palette = "Accent")
yplot <- ggdensity(fail_or_not_subset, "Pluritest.pluripotency.score", fill = "Neuroseq_Differentiated", palette = "Accent")+
  rotate()
# Cleaning the plots
yplot <- yplot + clean_theme()
xplot <- xplot + clean_theme()
# Arranging the plot
ggarrange(xplot, NULL, scaterplot, yplot,
          ncol = 2, nrow = 2,  align = "hv",
          widths = c(2, 1), heights = c(1, 2),
          common.legend = TRUE)


ggsave(
  paste0(
    output_dir,
    "/hipsci.",
    trait,
    "_pluritest_vs_differentiated_yesOrFAIL.cell_lines.png"
  ),
  device = "png",
  height = 8,
  width = 8,
  units = "in",
  dpi = 400
)

# boxplots of Pluripotency vs novelty scores by differentiation_status
## only those that were differentiated

melted = melt(fail_or_not_subset,id.vars = c("sampleID","Neuroseq_Differentiated"),
              measure.vars = c("Pluritest.novelty.score"), value.name = "Pluritest.novelty.score")

p1 <- ggplot(melted, aes(x=Neuroseq_Differentiated, y=Pluritest.novelty.score)) +
  geom_boxplot() +   theme_bw()

p1



melted = melt(fail_or_not_subset,id.vars = c("sampleID","Neuroseq_Differentiated"),
              measure.vars = c("Pluritest.pluripotency.score"), value.name = "Pluritest.pluripotency.score")

p2 <- ggplot(melted, aes(x=Neuroseq_Differentiated, y=Pluritest.pluripotency.score)) +
  geom_boxplot() +   theme_bw()

p2
# Arranging the plot
ggarrange(p1,  p2,
          ncol = 2, nrow = , align = "hv",
          widths = c(2, 2), heights = c(1),
          common.legend = TRUE)



ggsave(
  paste0(
    output_dir,
    "/hipsci.",
    trait,
    "_pluritest_vs_diffYesorFAIL.cell_lines.png"
  ),
  device = "png",
  height = 4,
  width = 9,
  units = "in",
  dpi = 400
)





# Pluripotency vs novelty scores vs continuous Differentiation status

# only those that have been differentiated (Neuroseq_Differentiated == yes)


scaterplot1 <- ggscatter(fail_or_not_subset, y = "Pluritest.novelty.score", x = "Neuroseq_Differentiation_efficiency",
                        size = 3, alpha = 0.6) +

  border()  + ylab("Novelty score") + xlab("Neuroseq_Differentiation_efficiency")

scaterplot2 <- ggscatter(fail_or_not_subset, y = "Pluritest.pluripotency.score", x = "Neuroseq_Differentiation_efficiency",
                         size = 3, alpha = 0.6) +

  border()  + ylab("Pluritest.pluripotency.score") + xlab("Neuroseq_Differentiation_efficiency")

# Arranging the plot
ggarrange(scaterplot1,scaterplot2,
          ncol = 2, nrow = 1, align = "hv",
          widths = c(2, 2), heights = c(1),
          common.legend = TRUE)



ggsave(
  paste0(
    output_dir,
    "/hipsci.",
    trait,
    "_pluritest_vs_ContdifferentiationStatus.cell_lines.png"
  ),
  device = "png",
  height = 4,
  width = 9,
  units = "in",
  dpi = 400
)

## Higher novelty scores mean worse fit of model (therefore more unreliable pluripotency scores)

# Setting theshold of scores:

low_novelty = fail_or_not_subset[fail_or_not_subset$Pluritest.novelty.score<1.25,]


scaterplot1 <- ggscatter(low_novelty, y = "Pluritest.novelty.score", x = "Neuroseq_Differentiation_efficiency",
                         size = 3, alpha = 0.6) +

  border()  + ylab("Novelty score") + xlab("Neuroseq_Differentiation_efficiency")

scaterplot2 <- ggscatter(low_novelty, y = "Pluritest.pluripotency.score", x = "Neuroseq_Differentiation_efficiency",
                         size = 3, alpha = 0.6) +

  border()  + ylab("Pluritest.pluripotency.score") + xlab("Neuroseq_Differentiation_efficiency")

# Arranging the plot
ggarrange(scaterplot1,scaterplot2,
          ncol = 2, nrow = 1, align = "hv",
          widths = c(2, 2), heights = c(1),
          common.legend = TRUE)



ggsave(
  paste0(
    output_dir,
    "/hipsci.",
    trait,
    "_pluritest_vs_ContdifferentiationStatus_lowNovelty.cell_lines.png"
  ),
  device = "png",
  height = 4,
  width = 9,
  units = "in",
  dpi = 400
)


#Diff efficiency vs quartiles


melted = melt(fail_or_not_subset,id.vars = c("sampleID","quartile"),
              measure.vars = c("Neuroseq_Differentiation_efficiency"), value.name = "Neuroseq_Differentiation_efficiency")

p1 <- ggplot(melted, aes(x=quartile, y=Neuroseq_Differentiation_efficiency)) +
  geom_boxplot() +   theme_bw()

p1
ggsave(
  paste0(
    output_dir,
    "/hipsci.",
    trait,
    "_ContdifferentiationStatus_by_quartiles.cell_lines.boxplot.png"
  ),
  device = "png",
  height = 4,
  width = 5,
  units = "in",
  dpi = 400
)

scaterplot1 <- ggscatter(fail_or_not_subset, y = "overallHR_quantile", x = "Neuroseq_Differentiation_efficiency",
                         size = 3, alpha = 0.6)
ggsave(
  paste0(
    output_dir,
    "/hipsci.",
    trait,
    "_ContdifferentiationStatus_scatterplot_HRquantile.cell_lines.boxplot.png"
  ),
  device = "png",
  height = 4,
  width = 5,
  units = "in",
  dpi = 400
)
## save the merged cell line table

write.table(cell_donor_list$cell_lines,
            file = paste0(input_dir,"/AD_cell_lines_other_info.txt"),sep = "\t",col.names = T,row.names = F, quote = F)


write.table(non_disease_European,
            file = paste0(input_dir,"/AD_cell_lines_other_info_non_diseased_European.txt"),sep = "\t",col.names = T,row.names = F, quote = F)
