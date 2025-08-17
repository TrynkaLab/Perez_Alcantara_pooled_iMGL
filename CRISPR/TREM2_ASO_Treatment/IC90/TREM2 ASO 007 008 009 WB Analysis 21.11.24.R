library(dplyr)
library(ggplot2)
library(cowplot)
library(ggbeeswarm)
library(ggpubr)

## load in data
setwd(dir = "~/Cambridge/Validation/TREM2 ASO/")

results= read.csv("TREM2 ASO WB 007 008 009 WB Analysis 21.11.24.csv",header=TRUE)

results$ASO=factor(results$ASO)
results$Condition
results$Condition[results$Condition==171]="TREM2"
results$Condition[results$Condition==192]="TREM2"
results$Condition[results$Condition=="171 SC1"]="Scramble 1"
results$Condition[results$Condition=="171 SC2"]="Scramble 2"
results$Condition[results$Condition=="192 SC2"]="Scramble 2"
results$Condition[results$Condition=="192 SC1"]="Scramble 1"
results$Repeat[results$Repeat=="ASO007"]=1
results$Repeat[results$Repeat=="ASO008"]=2
results$Repeat[results$Repeat=="ASO009"]=3

results_noWT=results[!results$Condition=="WT",]

results171=results_noWT[results_noWT$ASO==171,]
results171$Condition=factor(results171$Condition,levels=c("TREM2","Scramble 1","Scramble 2"))

results192=results_noWT[results_noWT$ASO==192,]
results192$Condition=factor(results192$Condition,levels=c("TREM2","Scramble 1","Scramble 2"))

### Plot Data

cbPalette=c("#999999","#E69F00","#56B4E9","#009E73","#f0e442","#0072b2","#d55e00","#cc79a7")
my_comparisons3=list(c("Scramble 1","Scramble 2"),c("Scramble 1","TREM2"),c("Scramble 2","TREM2"))
my_comparisons4=list(c("Scramble 1","Scramble 2"),c("Scramble 1","TREM2"),c("Scramble 2","TREM2"))


IC90_WB171=ggplot(results171,aes(x=Condition, y=Normalised.to.WT))+
  geom_boxplot(width=0.5,
               outlier.shape = NA,
               fill=cbPalette[2:4])+ 
  geom_point(aes(shape=Repeat,color=Line),
             size=2,
             position = position_dodge(width=0.2),
             alpha=0.8
  )+
  stat_compare_means(method="anova",size=4)+
  stat_compare_means(comparisons = my_comparisons3,label="p.signif",size=4)+
  
  geom_hline(yintercept = 1,
             color="grey",
             linetype="dashed")+
  theme_classic(base_size=7)+
  theme(plot.title = element_text(face="bold"),axis.title = element_text(face="bold"))+
  scale_fill_manual(values=cbPalette[2:4],name="Cell Line")+
  labs(y="Relative TREM2 Expression (to Vehicle)",
       x="ASO",
       title="TREM2-171 IC90 Protein Levels")+
  
  ylim(0,5)

print(IC90_WB171)

IC90_WB192=ggplot(results192,aes(x=Condition, y=Normalised.to.WT))+
  theme_classic(base_size = 7)+
  theme(plot.title = element_text(size=7),
        axis.title = element_text(size=7))+
  geom_boxplot(width=0.5,
               outlier.shape = NA,
               fill=cbPalette[2:4])+ 
  geom_point(aes(shape=Repeat,color=Line),
             size=1,
             position = position_dodge(width=0.2),
             alpha=0.8
             )+
  stat_compare_means(method="anova",size=2)+
  stat_compare_means(comparisons = my_comparisons4,label="p.signif",size=3)+
  geom_hline(yintercept = 1,
             color="grey",
             linetype="dashed")+
  scale_fill_manual(values=cbPalette[2:4],name="Cell Line")+
  labs(y="Relative TREM2 Expression (to Vehicle)",
       x="ASO",
       title="TREM2-192 IC90 Protein Levels")+
  theme(plot.title = element_text(face="bold"),axis.title = element_text(face="bold"))+
  ylim(0,3)

print(IC90_WB192)


IC90_ASO_WB=plot_grid(IC90_WB171,IC90_WB192,rel_widths = c(1,1))
print(IC90_ASO_WB)
ggsave("FIGURES/TREM2_ASO_IC90_WB.pdf",plot=last_plot(),device="pdf",dpi=600,units=c("mm"),height=90,width=180)

