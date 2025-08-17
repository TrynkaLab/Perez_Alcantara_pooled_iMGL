library(dplyr)
library(ggplot2)
library(cowplot)
library(ggpubr)

## load in the data

setwd("~/Cambridge/Validation/TREM2 LV/LV TREM001_TREM002_TREM003 FACS Analysis 09.08.24/")

AOWH_1=read.csv("LVTREM2001_AOWH2_Analysis_1_09.08.24.csv",header=TRUE)
AOWH_2=read.csv("LVTREM2002_AOWH2_Analysis_1_09.08.24.csv",header=TRUE)
AOWH_3=read.csv("LVTREM2003_AOWH2_Analysis_1_09.08.24.csv",header=TRUE)

HEGP_1=read.csv("LVTREM2001_HEGP3_Analysis_1_09.08.24.csv",header=TRUE)
HEGP_2=read.csv("LVTREM2002_HEGP3_Analysis_1_09.08.24.csv",header=TRUE)
HEGP_3=read.csv("LVTREM2003_HEGP3_Analysis_1_09.08.24.csv",header=TRUE)

KOLF_1=read.csv("LVTREM2001_KOLF21S_Analysis_1_09.08.24.csv",header=TRUE)
KOLF_2=read.csv("LVTREM2002_KOLF21S_Analysis_1_09.08.24.csv",header=TRUE)
KOLF_3=read.csv("LVTREM2003_KOLF21S_Analysis_1_09.08.24.csv",header=TRUE)

### add bioner line data

TREM2_KO_001_iMGL = read.csv("~/Cambridge/Validation/TREM2 KO Lines/TREM2 KO 001/TREM2_KO_001_iMGL_Phago.csv",header=TRUE)
TREM2_KO_002_iMGL= read.csv("~/Cambridge/Validation/TREM2 KO Lines/TREM2 KO 002/TREM2_KO_002_iMGL_Phago.csv",header=TRUE)
TREM2_KO_003_iMGL = read.csv("~/Cambridge/Validation/TREM2 KO Lines/TREM2 KO 003/TREM2_KO_003_iMGL_Phago.csv",header=TRUE)

TREM2_KO_iMGL_WT = rbind(TREM2_KO_001_iMGL,TREM2_KO_002_iMGL,TREM2_KO_003_iMGL)
TREM2_KO_iMGL_WT= TREM2_KO_iMGL_WT[TREM2_KO_iMGL_WT$TREM2_Genotype=="WT",]
TREM2_KO_iMGL_WT$Line="BIONi010"
TREM2_KO_iMGL_WT = TREM2_KO_iMGL_WT[,c(1:7,14,10,8,11,13)]


merged_data=rbind(AOWH_1,AOWH_2,AOWH_3,HEGP_1,HEGP_2,HEGP_3,KOLF_1,KOLF_2,KOLF_3,TREM2_KO_iMGL_WT)
head(merged_data)

colnames(merged_data)=c("Overlay","File Name","Number Singlets","Gate","Number mCherry","mCherry","% mCherry Positive Events","Line","Condition","Gating","Puro","Repeat")

density(merged_data$`Number Singlets`)
plot(density(merged_data$`Number Singlets`))

merged_data$Sample=paste(
  merged_data$Line,
  merged_data$Condition,
  merged_data$Gating,
  merged_data$Puro,
  merged_data$Repeat,
  sep="_"
)

summary_stats= merged_data %>% 
  group_by(Line,Condition,Gating,Puro,Repeat) %>%
  dplyr::summarise(
    n=n(),
    mean=mean(mCherry),
    sd=sd(mCherry)
  )

summary_stats$Repeat=as.factor(summary_stats$Repeat)
summary_stats=summary_stats[!summary_stats$Gating==1,]


summary_stats_LV01_AOWH2=summary_stats[summary_stats$Line=="AOWH_2" & summary_stats$Repeat==1,]
summary_stats_LV01_AOWH2_WT=summary_stats_LV01_AOWH2[which(summary_stats_LV01_AOWH2$Condition=="WT"),7]
summary_stats_LV01_AOWH2$Normalised_WT=summary_stats_LV01_AOWH2$mean/summary_stats_LV01_AOWH2_WT$mean

summary_stats_LV02_AOWH2=summary_stats[summary_stats$Line=="AOWH_2" & summary_stats$Repeat==2,]
summary_stats_LV02_AOWH2_WT=summary_stats_LV02_AOWH2[which(summary_stats_LV02_AOWH2$Condition=="WT"),7]
summary_stats_LV02_AOWH2$Normalised_WT=summary_stats_LV02_AOWH2$mean/summary_stats_LV02_AOWH2_WT$mean

summary_stats_LV03_AOWH2=summary_stats[summary_stats$Line=="AOWH_2" & summary_stats$Repeat==3,]
summary_stats_LV03_AOWH2_WT=summary_stats_LV03_AOWH2[which(summary_stats_LV03_AOWH2$Condition=="WT"),7]
summary_stats_LV03_AOWH2$Normalised_WT=summary_stats_LV03_AOWH2$mean/summary_stats_LV03_AOWH2_WT$mean

summary_stats_LV01_HEGP3=summary_stats[summary_stats$Line=="HEGP_3" & summary_stats$Repeat==1,]
summary_stats_LV01_HEGP3_WT=summary_stats_LV01_HEGP3[which(summary_stats_LV01_HEGP3$Condition=="WT"),7]
summary_stats_LV01_HEGP3$Normalised_WT=summary_stats_LV01_HEGP3$mean/summary_stats_LV01_HEGP3_WT$mean

summary_stats_LV02_HEGP3=summary_stats[summary_stats$Line=="HEGP_3" & summary_stats$Repeat==2,]
summary_stats_LV02_HEGP3_WT=summary_stats_LV02_HEGP3[which(summary_stats_LV02_HEGP3$Condition=="WT"),7]
summary_stats_LV02_HEGP3$Normalised_WT=summary_stats_LV02_HEGP3$mean/summary_stats_LV02_HEGP3_WT$mean

summary_stats_LV03_HEGP3=summary_stats[summary_stats$Line=="HEGP_3" & summary_stats$Repeat==3,]
summary_stats_LV03_HEGP3_WT=summary_stats_LV03_HEGP3[which(summary_stats_LV03_HEGP3$Condition=="WT"),7]
summary_stats_LV03_HEGP3$Normalised_WT=summary_stats_LV03_HEGP3$mean/summary_stats_LV03_HEGP3_WT$mean

summary_stats_LV01_KOLF21S=summary_stats[summary_stats$Line=="KOLF2.1S" & summary_stats$Repeat==1,]
summary_stats_LV01_KOLF21S_WT=summary_stats_LV01_KOLF21S[which(summary_stats_LV01_KOLF21S$Condition=="WT"),7]
summary_stats_LV01_KOLF21S$Normalised_WT=summary_stats_LV01_KOLF21S$mean/summary_stats_LV01_KOLF21S_WT$mean

summary_stats_LV02_KOLF21S=summary_stats[summary_stats$Line=="KOLF2.1S" & summary_stats$Repeat==2,]
summary_stats_LV02_KOLF21S_WT=summary_stats_LV02_KOLF21S[which(summary_stats_LV02_KOLF21S$Condition=="WT"),7]
summary_stats_LV02_KOLF21S$Normalised_WT=summary_stats_LV02_KOLF21S$mean/summary_stats_LV02_KOLF21S_WT$mean

summary_stats_LV03_KOLF21S=summary_stats[summary_stats$Line=="KOLF2.1S" & summary_stats$Repeat==3,]
summary_stats_LV03_KOLF21S_WT=summary_stats_LV03_KOLF21S[which(summary_stats_LV03_KOLF21S$Condition=="WT"),7]
summary_stats_LV03_KOLF21S$Normalised_WT=summary_stats_LV03_KOLF21S$mean/summary_stats_LV03_KOLF21S_WT$mean

summary_stats_LV01_BIONi010=summary_stats[summary_stats$Line=="BIONi010" & summary_stats$Repeat==1,]
summary_stats_LV01_BIONi010_WT=summary_stats_LV01_BIONi010[which(summary_stats_LV01_BIONi010$Condition=="WT"),7]
summary_stats_LV01_BIONi010$Normalised_WT=summary_stats_LV01_BIONi010$mean/summary_stats_LV01_BIONi010_WT$mean

summary_stats_LV02_BIONi010=summary_stats[summary_stats$Line=="BIONi010" & summary_stats$Repeat==2,]
summary_stats_LV02_BIONi010_WT=summary_stats_LV02_BIONi010[which(summary_stats_LV02_BIONi010$Condition=="WT"),7]
summary_stats_LV02_BIONi010$Normalised_WT=summary_stats_LV02_BIONi010$mean/summary_stats_LV02_BIONi010_WT$mean

summary_stats_LV03_BIONi010=summary_stats[summary_stats$Line=="BIONi010" & summary_stats$Repeat==3,]
summary_stats_LV03_BIONi010_WT=summary_stats_LV03_BIONi010[which(summary_stats_LV03_BIONi010$Condition=="WT"),7]
summary_stats_LV03_BIONi010$Normalised_WT=summary_stats_LV03_BIONi010$mean/summary_stats_LV03_BIONi010_WT$mean

summary_stats=rbind(
  summary_stats_LV01_AOWH2,
  summary_stats_LV02_AOWH2,
  summary_stats_LV03_AOWH2,
  summary_stats_LV01_HEGP3,
  summary_stats_LV02_HEGP3,
  summary_stats_LV03_HEGP3,
  summary_stats_LV01_KOLF21S,
  summary_stats_LV02_KOLF21S,
  summary_stats_LV03_KOLF21S,
  summary_stats_LV01_BIONi010,
  summary_stats_LV02_BIONi010,
  summary_stats_LV03_BIONi010
)

summary_stats$Condition[summary_stats$Condition=="NT"]="Intergenic"
summary_stats$Condition[summary_stats$Condition=="INT"]="Intergenic"
summary_stats$Condition[summary_stats$Condition=="TREM2"]="TREM2"
summary_stats$Line[summary_stats$Line=="BIONi010"]="BIONi010-C"
summary_stats$Line=factor(summary_stats$Line,levels=c("AOWH_2","HEGP_3","KOLF2.1S","BIONi010-C"))

summary_puro=summary_stats[summary_stats$Puro==1,]
summary_puro$Condition=factor(summary_puro$Condition,levels=c("TREM2","Intergenic","VPX","WT"))
summary_no_puro=summary_stats[summary_stats$Puro==0,]
summary_no_puro$Condition=factor(summary_no_puro$Condition,levels=c("TREM2","Intergenic","VPX","WT"))

### Plotted as % phago

cbPalette=c("#999999","#E69F00","#56B4E9","#009E73","#f0e442","#0072b2","#d55e00","#cc79a7")

plot1=ggplot(summary_no_puro,aes(x=Condition, y=mean))+
 geom_boxplot(width=0.5,size=0.4,fill=cbPalette[2:5],outlier.shape = 0)+ 
 geom_point(
    aes(color=Line,shape=Repeat),
    size=2,
    position=position_dodge(0.8))+
 theme_classic()+
  labs(y="% Phagocytosis",x="sgRNA",title="Without Puro")+
  ylim(0,55)
  
      
plot2=ggplot(summary_puro,aes(x=Condition, y=mean))+
  geom_boxplot(width=0.4,size=0.4,fill=cbPalette[2:3],outlier.shape = FALSE)+ 
  geom_point(
    aes(color=Line,shape=Repeat),
    size=2,
    position=position_dodge(0.5))+
  theme_classic()+
  labs(y="% Phagocytosis",x="sgRNA",title="With Puro")+
  ylim(0,55)

plot_grid(plot1,plot2)

### Plotted normalised to WT of each condition

summary_no_puro_noWT=summary_no_puro[!summary_no_puro$Condition=="WT",]
summary_puro_noWT=summary_puro[!summary_puro$Condition=="WT",]

my_comparisons1=list(c("TREM2","Intergenic"))
my_comparisons2=list(c("TREM2","Intergenic"),c("VPX","Intergenic"),c("TREM2","VPX"))


PHALV001=ggplot(summary_no_puro_noWT,aes(x=Condition, y=Normalised_WT))+
  geom_hline(yintercept = 1,color="grey",linetype="dashed")+
  geom_boxplot(width=0.5,size=0.4,fill=cbPalette[c(5,6,4)],outlier.shape = NA)+ 
  geom_point(
    aes(color=Line,shape=Repeat),
    size=1,
    position=position_dodge(0.2))+
  stat_compare_means(method="anova",size=2)+
  stat_compare_means(comparisons = my_comparisons2,label="p.signif",size=3)+
  theme_classic(base_size=7)+
  ylim(0,4)+
  theme(plot.title = element_text(face="bold"),axis.title = element_text(face="bold"))+
  labs(y="Normalised Phagocytosis (to WT)",x="Perturbation",title="Dual Reporter Assay")



PHALV002=ggplot(summary_puro_noWT,aes(x=Condition, y=Normalised_WT))+
  geom_hline(yintercept = 1,color="grey",linetype="dashed")+
  geom_boxplot(width=0.5,size=0.4,fill=cbPalette[5:6],outlier.shape = NA)+ 
  geom_point(
    aes(color=Line,shape=Repeat),
    size=1,
    position=position_dodge(0.2))+
  geom_line(aes(x=Condition,y=Normalised_WT,group=interaction(Repeat,Line),color=Line),position = position_dodge(0.2))+
  stat_compare_means(comparisons = my_comparisons1,label="p.signif",label.y=4.5,size=3)+
  theme_classic(base_size=7)+
  ylim(0,5)+
  theme(plot.title = element_text(face="bold"),axis.title = element_text(face="bold"))+
  labs(y="Normalised Phagocytosis (to WT)",x="Perturbation",title="Dual Reporter Assay")

Arrayed_lenti_puro = plot_grid(WBLV002,PHALV002)
Arrayed_lenti_no_puro = plot_grid(WBLV001,PHALV001)
print(Arrayed_lenti_puro)
print(Arrayed_lenti_no_puro)


ggsave("Arrayed_TREM2_Lenti_No_Puro_WB_Phago_with_bionio.pdf",plot=Arrayed_lenti_no_puro,device="pdf",dpi=600,unit=("mm"),height=90,width=180)
ggsave("Arrayed_TREM2_Lenti_No_Puro_WB_Phago_with_bionio.jpeg",plot=Arrayed_lenti_no_puro,device="jpeg",dpi=600,unit=("mm"),height=90,width=180)


ggsave("Arrayed_TREM2_Lenti_Puro_WB_Phago_with_bionio.pdf",plot=Arrayed_lenti_puro,device="pdf",dpi=600,unit=("mm"),height=90,width=180)
ggsave("Arrayed_TREM2_Lenti_Puro_WB_Phago_with_bionio.jpeg",plot=Arrayed_lenti_puro,device="jpeg",dpi=600,unit=("mm"),height=120,width=220)



