library(dplyr)
library(ggplot2)
library(cowplot)

## load in the data

setwd("~/Cambridge/Validation/TREM2 ASO/ASO TREM2 004 005 006 FACS Analysis/")

AOWH_1=read.csv("ASOTREM2004_AOWH2_Analysis_1_14.08.24.csv",header=TRUE)
AOWH_2=read.csv("ASOTREM2005_AOWH2_Analysis_1_14.08.24.csv",header=TRUE)
AOWH_3=read.csv("ASOTREM2006_AOWH2_Analysis_1_14.08.24.csv",header=TRUE)

HEGP_1=read.csv("ASOTREM2004_HEGP3_Analysis_1_14.08.24.csv",header=TRUE)
HEGP_2=read.csv("ASOTREM2005_HEGP3_Analysis_1_14.08.24.csv",header=TRUE)
HEGP_3=read.csv("ASOTREM2006_HEGP3_Analysis_1_14.08.24.csv",header=TRUE)

KOLF_1=read.csv("ASOTREM2004_KOLF21S_Analysis_1_14.08.24.csv",header=TRUE)
KOLF_2=read.csv("ASOTREM2005_KOLF21S_Analysis_1_14.08.24.csv",header=TRUE)
KOLF_3=read.csv("ASOTREM2006_KOLF21S_Analysis_1_14.08.24.csv",header=TRUE)

merged_data=rbind(AOWH_1,AOWH_2,AOWH_3,HEGP_1,HEGP_2,HEGP_3,KOLF_1,KOLF_2,KOLF_3)
head(merged_data)

colnames(merged_data)=c("Overlay","File Name","Number Singlets","Gate","Number mCherry","mCherry","% mCherry Positive Events","Condition","Gating","ASO","Repeat","Line")

density(merged_data$`Number Singlets`)
plot(density(merged_data$`Number Singlets`))

merged_data$Sample=paste(
  merged_data$Line,
  merged_data$Condition,
  merged_data$Gating,
  merged_data$ASO,
  merged_data$Repeat,
  sep="_"
)

summary_stats= merged_data %>% 
  group_by(Line,Condition,Gating,ASO,Repeat) %>%
  dplyr::summarise(
    n=n(),
    mean=mean(mCherry),
    sd=sd(mCherry)
  )

summary_stats$Repeat=as.factor(summary_stats$Repeat)
summary_stats=summary_stats[!summary_stats$Gating==1,]


summary_stats_ASO004_AOWH2=summary_stats[summary_stats$Line=="AOWH_2" & summary_stats$Repeat==1,]
summary_stats_ASO004_AOWH2_WT=summary_stats_ASO004_AOWH2[which(summary_stats_ASO004_AOWH2$Condition=="WT"),7]
summary_stats_ASO004_AOWH2$Normalised_WT=summary_stats_ASO004_AOWH2$mean/summary_stats_ASO004_AOWH2_WT$mean

summary_stats_ASO005_AOWH2=summary_stats[summary_stats$Line=="AOWH_2" & summary_stats$Repeat==2,]
summary_stats_ASO005_AOWH2_WT=summary_stats_ASO005_AOWH2[which(summary_stats_ASO005_AOWH2$Condition=="WT"),7]
summary_stats_ASO005_AOWH2$Normalised_WT=summary_stats_ASO005_AOWH2$mean/summary_stats_ASO005_AOWH2_WT$mean

summary_stats_ASO006_AOWH2=summary_stats[summary_stats$Line=="AOWH_2" & summary_stats$Repeat==3,]
summary_stats_ASO006_AOWH2_WT=summary_stats_ASO006_AOWH2[which(summary_stats_ASO006_AOWH2$Condition=="WT"),7]
summary_stats_ASO006_AOWH2$Normalised_WT=summary_stats_ASO006_AOWH2$mean/summary_stats_ASO006_AOWH2_WT$mean

summary_stats_ASO004_HEGP3=summary_stats[summary_stats$Line=="HEGP_3" & summary_stats$Repeat==1,]
summary_stats_ASO004_HEGP3_WT=summary_stats_ASO004_HEGP3[which(summary_stats_ASO004_HEGP3$Condition=="WT"),7]
summary_stats_ASO004_HEGP3$Normalised_WT=summary_stats_ASO004_HEGP3$mean/summary_stats_ASO004_HEGP3_WT$mean

summary_stats_ASO005_HEGP3=summary_stats[summary_stats$Line=="HEGP_3" & summary_stats$Repeat==2,]
summary_stats_ASO005_HEGP3_WT=summary_stats_ASO005_HEGP3[which(summary_stats_ASO005_HEGP3$Condition=="WT"),7]
summary_stats_ASO005_HEGP3$Normalised_WT=summary_stats_ASO005_HEGP3$mean/summary_stats_ASO005_HEGP3_WT$mean

summary_stats_ASO006_HEGP3=summary_stats[summary_stats$Line=="HEGP_3" & summary_stats$Repeat==3,]
summary_stats_ASO006_HEGP3_WT=summary_stats_ASO006_HEGP3[which(summary_stats_ASO006_HEGP3$Condition=="WT"),7]
summary_stats_ASO006_HEGP3$Normalised_WT=summary_stats_ASO006_HEGP3$mean/summary_stats_ASO006_HEGP3_WT$mean

summary_stats_ASO004_KOLF21S=summary_stats[summary_stats$Line=="KOLF2.1S" & summary_stats$Repeat==1,]
summary_stats_ASO004_KOLF21S_WT=summary_stats_ASO004_KOLF21S[which(summary_stats_ASO004_KOLF21S$Condition=="WT"),7]
summary_stats_ASO004_KOLF21S$Normalised_WT=summary_stats_ASO004_KOLF21S$mean/summary_stats_ASO004_KOLF21S_WT$mean

summary_stats_ASO005_KOLF21S=summary_stats[summary_stats$Line=="KOLF2.1S" & summary_stats$Repeat==2,]
summary_stats_ASO005_KOLF21S_WT=summary_stats_ASO005_KOLF21S[which(summary_stats_ASO005_KOLF21S$Condition=="WT"),7]
summary_stats_ASO005_KOLF21S$Normalised_WT=summary_stats_ASO005_KOLF21S$mean/summary_stats_ASO005_KOLF21S_WT$mean

summary_stats_ASO006_KOLF21S=summary_stats[summary_stats$Line=="KOLF2.1S" & summary_stats$Repeat==3,]
summary_stats_ASO006_KOLF21S_WT=summary_stats_ASO006_KOLF21S[which(summary_stats_ASO006_KOLF21S$Condition=="WT"),7]
summary_stats_ASO006_KOLF21S$Normalised_WT=summary_stats_ASO006_KOLF21S$mean/summary_stats_ASO006_KOLF21S_WT$mean

summary_stats=rbind(
  summary_stats_ASO004_AOWH2,
  summary_stats_ASO005_AOWH2,
  summary_stats_ASO006_AOWH2,
  summary_stats_ASO004_HEGP3,
  summary_stats_ASO005_HEGP3,
  summary_stats_ASO006_HEGP3,
  summary_stats_ASO004_KOLF21S,
  summary_stats_ASO005_KOLF21S,
  summary_stats_ASO006_KOLF21S
)

summary_stats$Condition
summary_stats$Condition[summary_stats$Condition==171]="TREM2"
summary_stats$Condition[summary_stats$Condition==192]="TREM2"
summary_stats$Condition[summary_stats$Condition=="171 SC1"]="Scramble 1"
summary_stats$Condition[summary_stats$Condition=="171 SC2"]="Scramble 2"
summary_stats$Condition[summary_stats$Condition=="192 SC2"]="Scramble 2"
summary_stats$Condition[summary_stats$Condition=="192 SC1"]="Scramble 1"

summary_171=summary_stats[summary_stats$ASO==171,]
summary_171$Condition=factor(summary_171$Condition,levels = c("TREM2","Scramble 1","Scramble 2"))

summary_192=summary_stats[summary_stats$ASO==192,]
summary_192$Condition=factor(summary_192$Condition,levels = c("TREM2","Scramble 1","Scramble 2"))


### Plotted as % phago

cbPalette=c("#999999","#E69F00","#56B4E9","#009E73","#f0e442","#0072b2","#d55e00","#cc79a7")
my_comparisons=list(c("Scramble 1","Scramble 2"),c("Scramble 1","TREM2"),c("Scramble 2","TREM2"))


ASO171_Phago=ggplot(summary_171,aes(x=Condition, y=mean))+
  geom_boxplot(width=0.5,size=0.4,fill=cbPalette[2:4],outlier.shape = NA)+ 
  geom_point(
    aes(color=Line,shape=Repeat),
    size=2,
    position=position_dodge(0.8))+
  theme_classic()+
  labs(y="% Phagocytosis",x="ASO",title="TREM2 ASO 171")+
  ylim(0,55)


ASO192_Phago=ggplot(summary_192,aes(x=Condition, y=mean))+
  geom_boxplot(width=0.4,size=0.4,fill=cbPalette[2:4],outlier.shape = NA)+ 
  geom_point(
    aes(color=Line,shape=Repeat),
    size=2,
    position=position_dodge(0.5))+
  theme_classic()+
  labs(y="% Phagocytosis",x="ASO",title="TREM2 ASO 192")+
  ylim(0,55)

plot_grid(ASO171_Phago,ASO192_Phago)

### Plotted normalised to WT of each condition

ASO171_IC50_Phago_Norm=ggplot(summary_171,aes(x=Condition, y=Normalised_WT))+
  geom_boxplot(width=0.5,size=0.4,fill=cbPalette[2:4],outlier.shape = NA)+ 
  geom_point(
    aes(color=Line,shape=Repeat),
    size=1,
    alpha=0.8,
    position=position_dodge(0.5))+
  theme_classic(base_size=7)+
  stat_compare_means(method="anova",size=2)+
  stat_compare_means(comparisons = my_comparisons,label="p.signif",size=3)+
  ylim(0,2)+
  geom_hline(yintercept = 1,color="grey",linetype="dashed")+
  theme(plot.title = element_text(face="bold"),axis.title = element_text(face="bold"))+
  labs(y="Normalised Phagocytosis (to Vehicle)",x="ASO",title="TREM2-171 IC50 Dual Reporter Assay")


ASO192_IC50_Phago_Norm=ggplot(summary_192,aes(x=Condition, y=Normalised_WT))+
  geom_boxplot(width=0.4,size=0.4,fill=cbPalette[2:4],outlier.shape = NA)+ 
  geom_point(
    aes(color=Line,shape=Repeat),
    size=1,
    alpha=0.8,
    position=position_dodge(0.5))+
  theme_classic(base_size = 7)+
  stat_compare_means(method="anova",size=2)+
  stat_compare_means(comparisons = my_comparisons,label="p.signif",size=3)+
  ylim(0,2)+
  geom_hline(yintercept = 1,color="grey",linetype="dashed")+
  theme(plot.title = element_text(face="bold"),axis.title = element_text(face="bold"))+
  labs(y="Normalised Phagocytosis (to Vehicle)",x="ASO",title="TREM2-192 IC50 Dual Reporter Assay")


ASO_phago_IC50=plot_grid(ASO171_IC50_Phago_Norm,ASO192_IC50_Phago_Norm)
print(ASO_phago_IC50)

ggsave("../FIGURES/TREM2_ASO_IC50_Phago.pdf",plot=last_plot(),device="pdf",dpi=600,units=c("mm"),height=90,width=180)
ggsave("../FIGURES/TREM2_ASO_IC50_Phago.jpeg",plot=last_plot(),device="jpeg",dpi=600,units=c("mm"),height=90,width=180)


plot_grid(ASO_IC50_WB,ASO_phago_IC50,ncol=1)
ggsave("../FIGURES/TREM2_ASO_IC50_WB_Phago.pdf",plot=last_plot(),device="pdf",dpi=600,units=c("mm"),height=180,width=180)
ggsave("../FIGURES/TREM2_ASO_IC50_WB_Phago.jpeg",plot=last_plot(),device="jpeg",dpi=600,units=c("mm"),height=180,width=180)


IC50_ASO_171=plot_grid(WB171,ASO171_IC50_Phago_Norm)
plot(IC50_ASO_171)
ggsave("../FIGURES/TREM2_ASO_171_IC50_WB_Phago.pdf",plot=last_plot(),device="pdf",dpi=600,units=c("mm"),height=90,width=180)
ggsave("../FIGURES/TREM2_ASO_171_IC50_WB_Phago.jpeg",plot=last_plot(),device="jpeg",dpi=600,units=c("mm"),height=90,width=180)

IC50_ASO_192=plot_grid(WB192,ASO192_IC50_Phago_Norm)
plot(IC50_ASO_192)
ggsave("../FIGURES/TREM2_ASO_192_IC50_WB_Phago.pdf",plot=last_plot(),device="pdf",dpi=600,units=c("mm"),height=90,width=180)
ggsave("../FIGURES/TREM2_ASO_192_IC50_WB_Phago.jpeg",plot=last_plot(),device="jpeg",dpi=600,units=c("mm"),height=90,width=180)

