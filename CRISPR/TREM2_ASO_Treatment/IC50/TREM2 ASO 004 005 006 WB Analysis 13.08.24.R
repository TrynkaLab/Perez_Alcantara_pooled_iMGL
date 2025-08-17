library(dplyr)
library(ggplot2)
library(cowplot)

## load in data

setwd("~/Cambridge/Validation/TREM2 ASO/")

results= read.csv("TREM2 ASO 004 005 006 WB Analysis 13.08.24.csv",header=TRUE)
results$Repeat=factor(results$Repeat)
results$ASO=factor(results$ASO)
results_noWT=results[!results$Name=="WT",]

results_noWT$Name[results_noWT$Name=="TREM2 ASO 171"]="TREM2"
results_noWT$Name[results_noWT$Name=="TREM2 ASO 192"]="TREM2"
results_noWT$Name[results_noWT$Name=="Scramble Ctrl 1"]="Scramble 1"
results_noWT$Name[results_noWT$Name=="Scramble Ctrl 2"]="Scramble 2"

results171=results_noWT[results_noWT$ASO==171,]
results171$Name=factor(results171$Name,levels=c("TREM2","Scramble 1","Scramble 2","WT"))

results192=results_noWT[results_noWT$ASO==192,]
results192$Name=factor(results192$Name,levels=c("TREM2","Scramble 1","Scramble 2","WT"))

### Plot Data

cbPalette=c("#999999","#E69F00","#56B4E9","#009E73","#f0e442","#0072b2","#d55e00","#cc79a7")
my_comparisons=list(c("Scramble 1","Scramble 2"),c("Scramble 1","TREM2"),c("Scramble 2","TREM2"))


WB171=ggplot(results171,aes(x=Name, y=Normalised.to.WT))+
  geom_boxplot(width=0.5,size=0.4,fill=cbPalette[2:4],outlier.shape = NA)+ 
  geom_point(
    aes(color=Line,shape=Repeat),
    size=1,
    alpha=.8,
    position=position_dodge(0.5))+
  theme_classic(base_size = 7)+
  stat_compare_means(method="anova",size=2)+
  stat_compare_means(comparisons = my_comparisons,label="p.signif",size=3)+
  geom_hline(yintercept = 1,color="grey",linetype="dashed")+
  labs(y="Relative TREM2 Expression (to Vehicle)",x="ASO",title="TREM2-171 IC50 Protein Levels")+
  theme(plot.title = element_text(face="bold"),axis.title = element_text(face="bold"))+
  ylim(0,2)

print(WB171)

WB192=ggplot(results192,aes(x=Name, y=Normalised.to.WT))+
  geom_boxplot(width=0.5,size=0.4,fill=cbPalette[2:4],outlier.shape = NA)+ 
  geom_point(
    aes(color=Line,shape=Repeat),
    size=1,
    alpha=.8,
    position=position_dodge(0.5))+
  theme_classic(base_size = 7)+
  stat_compare_means(method="anova",size=2)+
  stat_compare_means(comparisons = my_comparisons,label="p.signif",size=3)+
  geom_hline(yintercept = 1,color="grey",linetype="dashed")+
  labs(y="Relative TREM2 Expression (to Vehicle)",x="ASO",title="TREM2-192 IC50 Protein Levels")+
  theme(plot.title = element_text(face="bold"),axis.title = element_text(face="bold"))+
  ylim(0,2)

print(WB192)

ASO_IC50_WB=plot_grid(WB171,WB192)

print(ASO_IC50_WB)

ggsave("FIGURES/TREM2_ASO_IC50_WB.pdf",plot=last_plot(),device="pdf",dpi=600,units=c("mm"),height=90,width=180)


