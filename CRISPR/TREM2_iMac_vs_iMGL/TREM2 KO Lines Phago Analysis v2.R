library(dplyr)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(tidyr)
library(stringr)

## load in the data

setwd("~/Cambridge/Validation/TREM2 KO Lines/")


TREM2_KO_001_iMGL = read.csv("TREM2 KO 001/TREM2_KO_001_iMGL_Phago.csv",header=TRUE)
TREM2_KO_001_iMac = read.csv("TREM2 KO 001/TREM2_KO_001_iMac_Phago.csv",header=TRUE)

TREM2_KO_002_iMGL= read.csv("TREM2 KO 002/TREM2_KO_002_iMGL_Phago.csv",header=TRUE)
TREM2_KO_002_iMac = read.csv("TREM2 KO 002/TREM2_KO_002_iMac_Phago.csv",header=TRUE)

TREM2_KO_003_iMGL = read.csv("TREM2 KO 003/TREM2_KO_003_iMGL_Phago.csv",header=TRUE)
TREM2_KO_003_iMac = read.csv("TREM2 KO 003/TREM2_KO_003_iMac_Phago.csv",header=TRUE)

TREM2_KO_004_iMGL = read.csv("TREM2 KO 004/TREM2_KO_004_iMGL_Phago.csv",header=TRUE)
TREM2_KO_004_iMac = read.csv("TREM2 KO 004/TREM2_KO_004_iMac_Phago.csv",header=TRUE)

TREM2_KO_005_iMGL= read.csv("TREM2 KO 005/TREM2_KO_005_iMGL_Phago.csv",header=TRUE)
TREM2_KO_005_iMac = read.csv("TREM2 KO 005/TREM2_KO_005_iMac_Phago.csv",header=TRUE)

TREM2_KO_006_iMGL = read.csv("TREM2 KO 006/TREM2_KO_006_iMGL_Phago.csv",header=TRUE)
TREM2_KO_006_iMac = read.csv("TREM2 KO 006/TREM2_KO_006_iMac_Phago.csv",header=TRUE)

TREM2_KO_007_iMGL = read.csv("TREM2 KO 007/TREM2_KO_007_iMGL_Phago.csv",header=TRUE)
TREM2_KO_007_iMac = read.csv("TREM2 KO 007/TREM2_KO_007_iMac_Phago.csv",header=TRUE)

merged_data=rbind(TREM2_KO_001_iMac,TREM2_KO_001_iMGL,TREM2_KO_002_iMac,TREM2_KO_002_iMGL,TREM2_KO_003_iMac,TREM2_KO_003_iMGL)
colnames(merged_data)=c("Overlay","File Name","Number Singlets","Gate","Number mCherry","mCherry","% mCherry Positive Events","Gating","TREM2 Genotype","Condition","Puro","Cell Type","Repeat")
merged_data$Location="Cambridge"

merged_data_2=rbind(TREM2_KO_004_iMac,TREM2_KO_004_iMGL,TREM2_KO_005_iMac,TREM2_KO_005_iMGL,TREM2_KO_006_iMac,TREM2_KO_006_iMGL,TREM2_KO_007_iMac,TREM2_KO_007_iMGL)
head(merged_data_2)
colnames(merged_data_2)=c("Overlay","File Name","Number Singlets","Gate","Number mCherry","mCherry","% mCherry Positive Events","Gating","TREM2 Genotype","Condition","Puro","Cell Type","Repeat")
merged_data_2$Location="Oxford"

merged_data=rbind(merged_data,merged_data_2)

density(merged_data$`Number Singlets`)
plot(density(merged_data$`Number Singlets`))

merged_data$Sample=paste(
  merged_data$`Cell Type`,
  merged_data$`TREM2 Genotype`,
  merged_data$Condition,
  merged_data$Puro,
  merged_data$Repeat,
  merged_data$Gating,
  sep="_"
)

summary_stats= merged_data %>% 
  group_by(`Cell Type`,Condition,Gating,Puro,Repeat) %>%
  dplyr::summarise(
    n=n(),
    mean=mean(mCherry),
    sd=sd(mCherry)
  )

summary_stats$Repeat=as.factor(summary_stats$Repeat)
summary_stats=summary_stats[!summary_stats$Gating==1,]

### Normalising data to WT lines for each repeat

summary_stats_TREM2_KO_001_iMGL = summary_stats[summary_stats$`Cell Type`=="iMGL" & summary_stats$Repeat==1,]
summary_stats_TREM2_KO_001_iMGL_WT = summary_stats_TREM2_KO_001_iMGL[which(summary_stats_TREM2_KO_001_iMGL$Condition=="WT"),7]
summary_stats_TREM2_KO_001_iMGL$Normalised_WT=summary_stats_TREM2_KO_001_iMGL$mean/summary_stats_TREM2_KO_001_iMGL_WT$mean

summary_stats_TREM2_KO_002_iMGL = summary_stats[summary_stats$`Cell Type`=="iMGL" & summary_stats$Repeat==2,]
summary_stats_TREM2_KO_002_iMGL_WT = summary_stats_TREM2_KO_002_iMGL[which(summary_stats_TREM2_KO_002_iMGL$Condition=="WT"),7]
summary_stats_TREM2_KO_002_iMGL$Normalised_WT=summary_stats_TREM2_KO_002_iMGL$mean/summary_stats_TREM2_KO_002_iMGL_WT$mean

summary_stats_TREM2_KO_003_iMGL = summary_stats[summary_stats$`Cell Type`=="iMGL" & summary_stats$Repeat==3,]
summary_stats_TREM2_KO_003_iMGL_WT = summary_stats_TREM2_KO_003_iMGL[which(summary_stats_TREM2_KO_003_iMGL$Condition=="WT"),7]
summary_stats_TREM2_KO_003_iMGL$Normalised_WT=summary_stats_TREM2_KO_003_iMGL$mean/summary_stats_TREM2_KO_003_iMGL_WT$mean

summary_stats_TREM2_KO_004_iMGL = summary_stats[summary_stats$`Cell Type`=="iMGL" & summary_stats$Repeat==4,]
summary_stats_TREM2_KO_004_iMGL_WT = summary_stats_TREM2_KO_004_iMGL[which(summary_stats_TREM2_KO_004_iMGL$Condition=="WT"),7]
summary_stats_TREM2_KO_004_iMGL$Normalised_WT=summary_stats_TREM2_KO_004_iMGL$mean/summary_stats_TREM2_KO_004_iMGL_WT$mean

summary_stats_TREM2_KO_005_iMGL = summary_stats[summary_stats$`Cell Type`=="iMGL" & summary_stats$Repeat==5,]
summary_stats_TREM2_KO_005_iMGL_WT = summary_stats_TREM2_KO_005_iMGL[which(summary_stats_TREM2_KO_005_iMGL$Condition=="WT"),7]
summary_stats_TREM2_KO_005_iMGL$Normalised_WT=summary_stats_TREM2_KO_005_iMGL$mean/summary_stats_TREM2_KO_005_iMGL_WT$mean

summary_stats_TREM2_KO_006_iMGL = summary_stats[summary_stats$`Cell Type`=="iMGL" & summary_stats$Repeat==6,]
summary_stats_TREM2_KO_006_iMGL_WT = summary_stats_TREM2_KO_006_iMGL[which(summary_stats_TREM2_KO_006_iMGL$Condition=="WT"),7]
summary_stats_TREM2_KO_006_iMGL$Normalised_WT=summary_stats_TREM2_KO_006_iMGL$mean/summary_stats_TREM2_KO_006_iMGL_WT$mean

summary_stats_TREM2_KO_007_iMGL = summary_stats[summary_stats$`Cell Type`=="iMGL" & summary_stats$Repeat==7,]
summary_stats_TREM2_KO_007_iMGL_WT = summary_stats_TREM2_KO_007_iMGL[which(summary_stats_TREM2_KO_007_iMGL$Condition=="WT"),7]
summary_stats_TREM2_KO_007_iMGL$Normalised_WT=summary_stats_TREM2_KO_007_iMGL$mean/summary_stats_TREM2_KO_007_iMGL_WT$mean

summary_stats_TREM2_KO_004_iMGL_ITM = summary_stats[summary_stats$`Cell Type`=="iMGL_ITM" & summary_stats$Repeat==4,]
summary_stats_TREM2_KO_004_iMGL_ITM_WT = summary_stats_TREM2_KO_004_iMGL_ITM[which(summary_stats_TREM2_KO_004_iMGL_ITM$Condition=="WT"),7]
summary_stats_TREM2_KO_004_iMGL_ITM$Normalised_WT=summary_stats_TREM2_KO_004_iMGL_ITM$mean/summary_stats_TREM2_KO_004_iMGL_ITM_WT$mean

summary_stats_TREM2_KO_005_iMGL_ITM = summary_stats[summary_stats$`Cell Type`=="iMGL_ITM" & summary_stats$Repeat==5,]
summary_stats_TREM2_KO_005_iMGL_ITM_WT = summary_stats_TREM2_KO_005_iMGL_ITM[which(summary_stats_TREM2_KO_005_iMGL_ITM$Condition=="WT"),7]
summary_stats_TREM2_KO_005_iMGL_ITM$Normalised_WT=summary_stats_TREM2_KO_005_iMGL_ITM$mean/summary_stats_TREM2_KO_005_iMGL_ITM_WT$mean

summary_stats_TREM2_KO_006_iMGL_ITM = summary_stats[summary_stats$`Cell Type`=="iMGL_ITM" & summary_stats$Repeat==6,]
summary_stats_TREM2_KO_006_iMGL_ITM_WT = summary_stats_TREM2_KO_006_iMGL_ITM[which(summary_stats_TREM2_KO_006_iMGL_ITM$Condition=="WT"),7]
summary_stats_TREM2_KO_006_iMGL_ITM$Normalised_WT=summary_stats_TREM2_KO_006_iMGL_ITM$mean/summary_stats_TREM2_KO_006_iMGL_ITM_WT$mean

summary_stats_TREM2_KO_007_iMGL_ITM = summary_stats[summary_stats$`Cell Type`=="iMGL_ITM" & summary_stats$Repeat==7,]
summary_stats_TREM2_KO_007_iMGL_ITM_WT = summary_stats_TREM2_KO_007_iMGL_ITM[which(summary_stats_TREM2_KO_007_iMGL_ITM$Condition=="WT"),7]
summary_stats_TREM2_KO_007_iMGL_ITM$Normalised_WT=summary_stats_TREM2_KO_007_iMGL_ITM$mean/summary_stats_TREM2_KO_007_iMGL_ITM_WT$mean

summary_stats_TREM2_KO_001_iMac = summary_stats[summary_stats$`Cell Type`=="iMac" & summary_stats$Repeat==1,]
summary_stats_TREM2_KO_001_iMac_WT = summary_stats_TREM2_KO_001_iMac[which(summary_stats_TREM2_KO_001_iMac$Condition=="WT"),7]
summary_stats_TREM2_KO_001_iMac$Normalised_WT=summary_stats_TREM2_KO_001_iMac$mean/summary_stats_TREM2_KO_001_iMac_WT$mean

summary_stats_TREM2_KO_002_iMac = summary_stats[summary_stats$`Cell Type`=="iMac" & summary_stats$Repeat==2,]
summary_stats_TREM2_KO_002_iMac_WT = summary_stats_TREM2_KO_002_iMac[which(summary_stats_TREM2_KO_002_iMac$Condition=="WT"),7]
summary_stats_TREM2_KO_002_iMac$Normalised_WT=summary_stats_TREM2_KO_002_iMac$mean/summary_stats_TREM2_KO_002_iMac_WT$mean

summary_stats_TREM2_KO_003_iMac = summary_stats[summary_stats$`Cell Type`=="iMac" & summary_stats$Repeat==3,]
summary_stats_TREM2_KO_003_iMac_WT = summary_stats_TREM2_KO_003_iMac[which(summary_stats_TREM2_KO_003_iMac$Condition=="WT"),7]
summary_stats_TREM2_KO_003_iMac$Normalised_WT=summary_stats_TREM2_KO_003_iMac$mean/summary_stats_TREM2_KO_003_iMac_WT$mean

summary_stats_TREM2_KO_004_iMac = summary_stats[summary_stats$`Cell Type`=="iMac" & summary_stats$Repeat==4,]
summary_stats_TREM2_KO_004_iMac_WT = summary_stats_TREM2_KO_004_iMac[which(summary_stats_TREM2_KO_004_iMac$Condition=="WT"),7]
summary_stats_TREM2_KO_004_iMac$Normalised_WT=summary_stats_TREM2_KO_004_iMac$mean/summary_stats_TREM2_KO_004_iMac_WT$mean

summary_stats_TREM2_KO_005_iMac = summary_stats[summary_stats$`Cell Type`=="iMac" & summary_stats$Repeat==5,]
summary_stats_TREM2_KO_005_iMac_WT = summary_stats_TREM2_KO_005_iMac[which(summary_stats_TREM2_KO_005_iMac$Condition=="WT"),7]
summary_stats_TREM2_KO_005_iMac$Normalised_WT=summary_stats_TREM2_KO_005_iMac$mean/summary_stats_TREM2_KO_005_iMac_WT$mean

summary_stats_TREM2_KO_006_iMac = summary_stats[summary_stats$`Cell Type`=="iMac" & summary_stats$Repeat==6,]
summary_stats_TREM2_KO_006_iMac_WT = summary_stats_TREM2_KO_006_iMac[which(summary_stats_TREM2_KO_006_iMac$Condition=="WT"),7]
summary_stats_TREM2_KO_006_iMac$Normalised_WT=summary_stats_TREM2_KO_006_iMac$mean/summary_stats_TREM2_KO_006_iMac_WT$mean

summary_stats_TREM2_KO_007_iMac = summary_stats[summary_stats$`Cell Type`=="iMac" & summary_stats$Repeat==7,]
summary_stats_TREM2_KO_007_iMac_WT = summary_stats_TREM2_KO_007_iMac[which(summary_stats_TREM2_KO_007_iMac$Condition=="WT"),7]
summary_stats_TREM2_KO_007_iMac$Normalised_WT=summary_stats_TREM2_KO_007_iMac$mean/summary_stats_TREM2_KO_007_iMac_WT$mean



summary_stats=rbind(
  summary_stats_TREM2_KO_001_iMGL,
  summary_stats_TREM2_KO_001_iMac,
  summary_stats_TREM2_KO_002_iMGL,
  summary_stats_TREM2_KO_002_iMac,
  summary_stats_TREM2_KO_003_iMGL,
  summary_stats_TREM2_KO_003_iMac,
  summary_stats_TREM2_KO_004_iMGL,
  summary_stats_TREM2_KO_004_iMGL_ITM,
  summary_stats_TREM2_KO_004_iMac,
  summary_stats_TREM2_KO_005_iMGL,
  summary_stats_TREM2_KO_005_iMGL_ITM,
  summary_stats_TREM2_KO_005_iMac,
  summary_stats_TREM2_KO_006_iMGL,
  summary_stats_TREM2_KO_006_iMGL_ITM,
  summary_stats_TREM2_KO_006_iMac,
  summary_stats_TREM2_KO_007_iMGL,
  summary_stats_TREM2_KO_007_iMGL_ITM,
  summary_stats_TREM2_KO_007_iMac,
  
)

summary_stats$Condition=paste(summary_stats$Condition,summary_stats$Puro,sep=" ")

summary_stats$Condition
summary_stats$Condition[summary_stats$Condition=="INT 0"]="Intergenic Control"
summary_stats$Condition[summary_stats$Condition=="INT 1"]="Intergenic Control Puro"
summary_stats$Condition[summary_stats$Condition=="KO 0"]="TREM2 KO"
summary_stats$Condition[summary_stats$Condition=="TREM2 0"]="TREM2 Lenti KO"
summary_stats$Condition[summary_stats$Condition=="TREM2 1"]="TREM2 Lenti KO Puro"
summary_stats$Condition[summary_stats$Condition=="VPX 0"]="VPX"
summary_stats$Condition[summary_stats$Condition=="WT 0"]="WT"

summary_stats$Condition=factor(summary_stats$Condition,levels=
                                 c("WT","TREM2 KO","VPX","Intergenic Control","Intergenic Control Puro","TREM2 Lenti KO","TREM2 Lenti KO Puro"))


### Plotted as % phago

cbPalette=c("#999999","#E69F00","#56B4E9","#009E73","#f0e442","#0072b2","#d55e00","#cc79a7","#542788")

### Plotted normalised to WT of each condition

summary_stats_KO = summary_stats[summary_stats$Condition=="TREM2 KO",]
summary_stats_KO$Repeat=factor(summary_stats_KO$Repeat)
summary_stats_KO$grouping = paste(summary_stats_KO$`Cell Type`,summary_stats_KO$Condition,sep=" ")
my_comparisons=list(c("iMac","iMGL"),c("iMac","iMGL_ITM"),c("iMGL","iMGL_ITM"))

TREM2_KO_Lines_Normalised_Phago=ggplot(summary_stats_KO,aes(x=`Cell Type`, y=Normalised_WT))+
  geom_boxplot(width=0.5,
              outlier.shape = NA,
              fill=cbPalette[2:4])+ 
  geom_point(aes(shape=Repeat),
              size=1,
             color="black",
              position=position_dodge(0.2),
             alpha=.8)+
  theme_classic(base_size=7)+
    geom_hline(yintercept = 1,
             color="grey",
             linetype="dashed")+
  scale_y_continuous(trans="log2")+
  scale_shape_manual(values=1:nlevels(summary_stats_KO$Repeat))+
  stat_compare_means(comparisons = my_comparisons, label="p.signif",size=3)+
  theme(plot.title = element_text(face="bold"),axis.title = element_text(face="bold"))+
  labs(y="Log2 Normalised Phagocytosis (to BIONi010-C)",x="Cell Type",title="BIONi010-C-17")

print(TREM2_KO_Lines_Normalised_Phago)

ggsave(plot=last_plot(),"BIONi010C17_Phagocytosis2.pdf",device="pdf",dpi=300,unit="mm",width=90,height=90)

### Only imacs and iMGL

summary_stats_KO = summary_stats[summary_stats$Condition=="TREM2 KO",]
summary_stats_KO = summary_stats_KO[!summary_stats_KO$`Cell Type`=="iMGL_ITM",]
summary_stats_KO$Repeat=factor(summary_stats_KO$Repeat)
summary_stats_KO$grouping = paste(summary_stats_KO$`Cell Type`,summary_stats_KO$Condition,sep=" ")
my_comparisons=list(c("iMac","iMGL"))

TREM2_KO_Lines_Normalised_Phago=ggplot(summary_stats_KO,aes(x=`Cell Type`, y=Normalised_WT))+
  geom_boxplot(width=0.5,
               outlier.shape = NA,
               fill=cbPalette[2:3])+ 
  geom_point(aes(shape=Repeat),
             size=1,
             color="black",
             position=position_dodge(0.2),
             alpha=.8)+
  theme_classic(base_size=7)+
  geom_hline(yintercept = 1,
             color="grey",
             linetype="dashed")+
  scale_y_continuous(trans="log2")+
  scale_shape_manual(values=1:nlevels(summary_stats_KO$Repeat))+
  stat_compare_means(comparisons = my_comparisons, label="p.signif",size=3,method="wilcox.test")+
  theme(plot.title = element_text(face="bold"),axis.title = element_text(face="bold"))+
  labs(y="Log2 Normalised Phagocytosis (to BIONi010-C)",x="Cell Type",title="BIONi010-C-17")

print(TREM2_KO_Lines_Normalised_Phago)

ggsave(plot=last_plot(),"BIONi010C17_Phagocytosis3.pdf",device="pdf",dpi=300,unit="mm",width=90,height=90)
ggsave("TREM2_KO_lines_imacs.jpeg",plot=last_plot(),device="jpeg",dpi=600,unit=("mm"),height=70,width=50)





summary_stats_Lenti = summary_stats[!summary_stats$Condition=="TREM2 KO",]
summary_stats_Lenti = summary_stats_Lenti[!summary_stats_Lenti$Condition=="WT",]

my_comparisons=list(c("VPX","Intergenic Control"),c("VPX","Intergenic Control Puro"),c("VPX","TREM2 KO Puro"),c("Intergenic Control Puro","TREM2 KO Puro"))

TREM2_KO_Lenti_Normalised_Phago=ggplot(summary_stats_Lenti,aes(x=Condition, y=Normalised_WT))+
  geom_boxplot(width=0.4,
               outlier.shape = NA,
               fill=cbPalette[2:6])+ 
  geom_point(aes(shape=Repeat),
             size=1,
             position=position_dodge(0.2),
             alpha=.8)+
  theme_classic(base_size=7)+
    geom_hline(yintercept = 1,
             color="grey",
             linetype="dashed")+
  theme(plot.title = element_text(face="bold"),axis.title = element_text(face="bold"))+
  labs(y="Normalised Phagocytosis (to BIONi010-C iMGL)",x="Condition",title="BIONi010-C TREM2 iMGL Lentiviral KD")

print(TREM2_KO_Lenti_Normalised_Phago)
ggsave("TREM2_KO_lines_Validation_BIONi010-C_TREM2_KD.pdf",plot=last_plot(),device="pdf",dpi=300,,unit=("mm"),height=90,width=180)

#### Only Puro Selected

summary_stats_Lenti_Puro = summary_stats[!summary_stats$Puro==0,]
summary_stats_Lenti_Puro$Condition = str_split_i(summary_stats_Lenti_Puro$Condition," ",1)

summary_stats_Lenti_Puro$Condition[summary_stats_Lenti_Puro$Condition=="Intergenic"]="Intergenic"
summary_stats_Lenti_Puro$Condition[summary_stats_Lenti_Puro$Condition=="TREM2"]="TREM2"

summary_stats_Lenti_Puro$Condition=factor(summary_stats_Lenti_Puro$Condition,levels=
                                 c("TREM2","Intergenic"))



my_comparisons=list(c("TREM2","Intergenic"))

TREM2_KO_Lenti_Puro_Normalised_Phago=ggplot(summary_stats_Lenti_Puro,aes(x=Condition, y=Normalised_WT))+
  geom_boxplot(width=0.4,
               outlier.shape = NA,
               fill=cbPalette[5:6])+ 
  geom_point(aes(shape=Repeat),
             size=1,
             position=position_dodge(0.2),
             alpha=.8)+
  theme_classic(base_size=7)+
  geom_hline(yintercept = 1,
             color="grey",
             linetype="dashed")+
  stat_compare_means(label="p.signif", comparisons=my_comparisons,method="wilcox.test",paired = TRUE,size=3)+
  theme(plot.title = element_text(face="bold"),axis.title = element_text(face="bold"))+
  labs(y="Normalised Phagocytosis (to BIONi010-C iMGL)",x="Perturbation",title="BIONi010-C TREM2 iMGL Lentiviral KO")



plot_grid(TREM2_KO_Lines_Normalised_Phago,TREM2_KO_Lenti_Puro_Normalised_Phago,rel_widths = c(1,1))

ggsave("TREM2_KO_lines_Validation.pdf",plot=last_plot(),device="pdf",dpi=600,unit=("mm"),height=90,width=180)
ggsave("TREM2_KO_lines_Validation.jpeg",plot=last_plot(),device="jpeg",dpi=600,unit=("mm"),height=90,width=180)


summary_stats_Lenti_Puro_lm = lm(Normalised_WT ~ Condition + Repeat, data=summary_stats_Lenti_Puro)
summary(summary_stats_Lenti_Puro_lm)

summary_stats_KO_lm = lm(Normalised_WT ~ `Cell Type` + Repeat, data=summary_stats_KO)
summary(summary_stats_KO_lm)



plot_grid(NULL,NULL,Arrayed_lenti_puro,IC90_ASO_171,NULL,NULL,ncol=2,rel_heights = c(1,1,0.5),labels=c("A","B","C","D","E","F"))



ggsave("../FIGURES/TREM2 ASO IC90 Phagocytosis.jpeg",plot=last_plot(),device="tiff",dpi=300,units=c("mm"),height=90,width=180)


IC90_Figure=plot_grid(NULL,IC90_ASO_WB,NULL,ASO_phago_IC90,ncol=1,rel_heights = c(0.7,1,.5,1),labels="AUTO")
ggsave("~/Paper/Cambridge/Figure X ASO IC90.jpeg",plot=IC90_Figure,device="jpeg",dpi=600,width=10,height=15)
print(IC90_Figure)

