library(dplyr)
library(ggplot2)
library(cowplot)
library(ggpubr)

## load in the data

setwd("~/Cambridge/Validation/TREM2 ASO/ASO TREM2 007 008 009 FACS Analysis/")

AOWH_1=read.csv("ASOTREM2007_AOWH2_Analysis_1_05.11.24.csv",header=TRUE)
AOWH_2=read.csv("ASOTREM2008_AOWH2_Analysis_1_07.11.24.csv",header=TRUE)
AOWH_3=read.csv("ASOTREM2009_AOWH2_Analysis_1_12.11.24.csv",header=TRUE)

HEGP_1=read.csv("ASOTREM2007_HEGP3_Analysis_1_05.11.24.csv",header=TRUE)
HEGP_2=read.csv("ASOTREM2008_HEGP3_Analysis_1_07.11.24.csv",header=TRUE)
HEGP_3=read.csv("ASOTREM2009_HEGP3_Analysis_1_12.11.24.csv",header=TRUE)

KOLF_1=read.csv("ASOTREM2007_KOLF2.1S_Analysis_1_05.11.24.csv",header=TRUE)
KOLF_2=read.csv("ASOTREM2008_KOLF2.1S_Analysis_1_07.11.24.csv",header=TRUE)
KOLF_3=read.csv("ASOTREM2009_KOLF2.1S_Analysis_1_12.11.24.csv",header=TRUE)

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


summary_stats_ASO007_AOWH2=summary_stats[summary_stats$Line=="AOWH_2" & summary_stats$Repeat==1,]
summary_stats_ASO007_AOWH2_WT=summary_stats_ASO007_AOWH2[which(summary_stats_ASO007_AOWH2$Condition=="WT"),7]
summary_stats_ASO007_AOWH2$Normalised_WT=summary_stats_ASO007_AOWH2$mean/summary_stats_ASO007_AOWH2_WT$mean

summary_stats_ASO008_AOWH2=summary_stats[summary_stats$Line=="AOWH_2" & summary_stats$Repeat==2,]
summary_stats_ASO008_AOWH2_WT=summary_stats_ASO008_AOWH2[which(summary_stats_ASO008_AOWH2$Condition=="WT"),7]
summary_stats_ASO008_AOWH2$Normalised_WT=summary_stats_ASO008_AOWH2$mean/summary_stats_ASO008_AOWH2_WT$mean

summary_stats_ASO009_AOWH2=summary_stats[summary_stats$Line=="AOWH_2" & summary_stats$Repeat==3,]
summary_stats_ASO009_AOWH2_WT=summary_stats_ASO009_AOWH2[which(summary_stats_ASO009_AOWH2$Condition=="WT"),7]
summary_stats_ASO009_AOWH2$Normalised_WT=summary_stats_ASO009_AOWH2$mean/summary_stats_ASO009_AOWH2_WT$mean

summary_stats_ASO007_HEGP3=summary_stats[summary_stats$Line=="HEGP_3" & summary_stats$Repeat==1,]
summary_stats_ASO007_HEGP3_WT=summary_stats_ASO007_HEGP3[which(summary_stats_ASO007_HEGP3$Condition=="WT"),7]
summary_stats_ASO007_HEGP3$Normalised_WT=summary_stats_ASO007_HEGP3$mean/summary_stats_ASO007_HEGP3_WT$mean

summary_stats_ASO008_HEGP3=summary_stats[summary_stats$Line=="HEGP_3" & summary_stats$Repeat==2,]
summary_stats_ASO008_HEGP3_WT=summary_stats_ASO008_HEGP3[which(summary_stats_ASO008_HEGP3$Condition=="WT"),7]
summary_stats_ASO008_HEGP3$Normalised_WT=summary_stats_ASO008_HEGP3$mean/summary_stats_ASO008_HEGP3_WT$mean

summary_stats_ASO009_HEGP3=summary_stats[summary_stats$Line=="HEGP_3" & summary_stats$Repeat==3,]
summary_stats_ASO009_HEGP3_WT=summary_stats_ASO009_HEGP3[which(summary_stats_ASO009_HEGP3$Condition=="WT"),7]
summary_stats_ASO009_HEGP3$Normalised_WT=summary_stats_ASO009_HEGP3$mean/summary_stats_ASO009_HEGP3_WT$mean

summary_stats_ASO007_KOLF21S=summary_stats[summary_stats$Line=="KOLF2.1S" & summary_stats$Repeat==1,]
summary_stats_ASO007_KOLF21S_WT=summary_stats_ASO007_KOLF21S[which(summary_stats_ASO007_KOLF21S$Condition=="WT"),7]
summary_stats_ASO007_KOLF21S$Normalised_WT=summary_stats_ASO007_KOLF21S$mean/summary_stats_ASO007_KOLF21S_WT$mean

summary_stats_ASO008_KOLF21S=summary_stats[summary_stats$Line=="KOLF2.1S" & summary_stats$Repeat==2,]
summary_stats_ASO008_KOLF21S_WT=summary_stats_ASO008_KOLF21S[which(summary_stats_ASO008_KOLF21S$Condition=="WT"),7]
summary_stats_ASO008_KOLF21S$Normalised_WT=summary_stats_ASO008_KOLF21S$mean/summary_stats_ASO008_KOLF21S_WT$mean

summary_stats_ASO009_KOLF21S=summary_stats[summary_stats$Line=="KOLF2.1S" & summary_stats$Repeat==3,]
summary_stats_ASO009_KOLF21S_WT=summary_stats_ASO009_KOLF21S[which(summary_stats_ASO009_KOLF21S$Condition=="WT"),7]
summary_stats_ASO009_KOLF21S$Normalised_WT=summary_stats_ASO009_KOLF21S$mean/summary_stats_ASO009_KOLF21S_WT$mean

summary_stats=rbind(
  summary_stats_ASO007_AOWH2,
  summary_stats_ASO008_AOWH2,
  summary_stats_ASO009_AOWH2,
  summary_stats_ASO007_HEGP3,
  summary_stats_ASO008_HEGP3,
  summary_stats_ASO009_HEGP3,
  summary_stats_ASO007_KOLF21S,
  summary_stats_ASO008_KOLF21S,
  summary_stats_ASO009_KOLF21S
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
my_comparisons3=list(c("Scramble 1","Scramble 2"),c("Scramble 1","TREM2"),c("Scramble 2","TREM2"))
my_comparisons4=list(c("Scramble 1","Scramble 2"),c("Scramble 1","TREM2"),c("Scramble 2","TREM2"))

ASO171_Phago=ggplot(summary_171,aes(x=Condition, y=mean,fill=Line))+
  geom_boxplot(width=0.7,
               outlier.shape = NA)+ 
  geom_point(aes(group=Line,shape=Repeat),
    size=1,
    position=position_dodge(width=0.7))+
  theme_classic(base_size=7)+
  theme(legend.position = "none")+
  scale_fill_manual(values=cbPalette[2:4],name="Cell Line")+
  labs(y="% Phagocytosis",x="ASO",title="TREM2 ASO 171")+
  ylim(0,55)


ASO192_Phago=ggplot(summary_192,aes(x=Condition, y=mean,fill=Line))+
  geom_boxplot(width=0.7,
               outlier.shape = NA)+ 
  geom_point(aes(group=Line,shape=Repeat),
    size=1,
    position=position_dodge(0.7))+
  theme_classic(base_size=7)+
  scale_fill_manual(values=cbPalette[2:4],name="Cell Line")+
  labs(y="% Phagocytosis",x="ASO",title="TREM2 ASO 192")+
  ylim(0,55)

plot_grid(ASO171_Phago,ASO192_Phago,rel_widths=c(0.8,1))

### Plotted normalised to WT of each condition

ASO171_IC90_Phago_Norm=ggplot(summary_171,aes(x=Condition, y=Normalised_WT))+
  theme_classic(base_size=7)+
  geom_boxplot(width=.5,
              outlier.shape = NA,
              fill=cbPalette[2:4])+ 
  geom_point(aes(color=Line,shape=Repeat),
              size=2,
              position=position_dodge(0.2),
             alpha=.8)+
  stat_compare_means(method="anova",size=4)+
  stat_compare_means(comparisons = my_comparisons3,label="p.signif",size=4)+
  
  theme(plot.title = element_text(face="bold"),axis.title = element_text(face="bold"))+
  geom_hline(yintercept = 1,
             color="grey",
             linetype="dashed")+
  ylim(0,6)+
  labs(y="Normalised Phagocytosis (to Vehicle)",x="ASO",title="TREM2-171 IC90 Dual Reporter Assay")

print(ASO171_IC90_Phago_Norm)

ASO192_IC90_Phago_Norm=ggplot(summary_192,aes(x=Condition, y=Normalised_WT))+
  geom_boxplot(width=0.4,
               outlier.shape = NA,
               fill=cbPalette[2:4])+ 
  geom_point(aes(color=Line,shape=Repeat),
                size=1,
                position=position_dodge(0.2),
             alpha=.8)+
  geom_hline(yintercept = 1,
             color="grey",
             linetype="dashed")+
  stat_compare_means(method="anova",size=2)+
  theme_classic(base_size=7)+
  stat_compare_means(comparisons = my_comparisons4,label="p.signif",size=3)+
  theme(plot.title = element_text(size=7),
        axis.title = element_text(size=7))+
  ylim(0,6)+
  theme(plot.title = element_text(face="bold"),axis.title = element_text(face="bold"))+
  labs(y="Normalised Phagocytosis (to Vehicle)",x="ASO",title="TREM2-192 IC90 Dual Reporter Assay")


ASO_phago_IC90=plot_grid(ASO171_IC90_Phago_Norm,ASO192_IC90_Phago_Norm,rel_widths = c(1,1))

plot(ASO_phago_IC90)

ggsave("../FIGURES/TREM2_ASO_IC90_Phago.pdf",plot=last_plot(),device="pdf",dpi=600,units=c("mm"),height=90,width=180)
ggsave("../FIGURES/TREM2_ASO_IC90_Phago.jpeg",plot=last_plot(),device="jpeg",dpi=600,units=c("mm"),height=90,width=180)

plot_grid(IC90_ASO_WB,ASO_phago_IC90,ncol=1)
ggsave("../FIGURES/TREM2_ASO_IC90_WB_Phago.pdf",plot=last_plot(),device="pdf",dpi=600,units=c("mm"),height=180,width=180)
ggsave("../FIGURES/TREM2_ASO_IC90_WB_Phago.jpeg",plot=last_plot(),device="jpeg",dpi=600,units=c("mm"),height=180,width=180)

IC90_ASO_171=plot_grid(IC90_WB171,ASO171_IC90_Phago_Norm)
plot(IC90_ASO_171)
ggsave("../FIGURES/TREM2_ASO_171_IC90_WB_Phago.pdf",plot=last_plot(),device="pdf",dpi=600,units=c("mm"),height=90,width=180)
ggsave("../FIGURES/TREM2_ASO_171_IC90_WB_Phago.jpeg",plot=last_plot(),device="jpeg",dpi=600,units=c("mm"),height=90,width=180)

IC90_ASO_192=plot_grid(IC90_WB192,ASO192_IC90_Phago_Norm)
plot(IC90_ASO_192)
ggsave("../FIGURES/TREM2_ASO_192_IC90_WB_Phago.pdf",plot=last_plot(),device="pdf",dpi=600,units=c("mm"),height=90,width=180)
ggsave("../FIGURES/TREM2_ASO_192_IC90_WB_Phago.jpeg",plot=last_plot(),device="jpeg",dpi=600,units=c("mm"),height=90,width=180)




