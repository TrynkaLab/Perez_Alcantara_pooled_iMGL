library(dplyr)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(stringr)

## load in data

setwd("~/Cambridge/Validation/TREM2 LV/")

results= read.csv("LV TREM001_TREM002_TREM003_WB_Analysis_16.08.24.csv",header=TRUE)
results=results[!results$Samples=="NO TREM",]
results$Repeat=factor(results$Repeat)

results$Virus[results$Virus=="NT"]="Intergenic"
results$Virus[results$Virus=="T"]="TREM2"

results$Line[results$Line=="KOLF21S"]="KOLF2.1S"
results$Line[results$Line=="AOWH"]="AOWH_2"
results$Line[results$Line=="HEGP"]="HEGP_3"

## remove HEGP_3 rep 1
results=results[-(7:12),]

results$Virus[results$Virus=="PSIV"]="VPX"

results$Virus=factor(results$Virus,levels=c("TREM2","Intergenic","VPX","WT"))

results_noWT=results[!results$Virus=="WT",]
results_noWT_puro=results_noWT[results_noWT$Puro=="1",]
results_noWT_nopuro=results_noWT[results_noWT$Puro=="0",]

results_noWTpSIV=results_noWT[!results_noWT$Virus=="VPX",]
results_noWTpSIV_puro=results_noWTpSIV[results_noWTpSIV$Puro=="1",]
results_noWTpSIV_nopuro=results_noWTpSIV[results_noWTpSIV$Puro=="0",]

### Plot Data 

cbPalette=c("#999999","#E69F00","#56B4E9","#009E73","#f0e442","#0072b2","#d55e00","#cc79a7")
my_comparisons1=list(c("TREM2","Intergenic"))
my_comparisons2=list(c("TREM2","Intergenic"),c("VPX","Intergenic"),c("TREM2","VPX"))

WBLV001=ggplot(results_noWT_nopuro,aes(x=Virus, y=Normalised_to_WT))+
  geom_hline(yintercept = 1,color="grey",linetype="dashed")+
  geom_boxplot(width=0.5,size=0.4,fill=cbPalette[c(5,6,4)],outlier.shape = NA)+ 
  geom_point(
    aes(color=Line,shape=Repeat),
    size=1,
    position=position_dodge(0.2))+
  stat_compare_means(method="anova",size=2)+
  stat_compare_means(comparisons = my_comparisons2,label="p.signif",size=3)+
  theme_classic(base_size=7)+
  labs(y="Relative TREM2 Expression (to WT)",x="Perturbation",title="TREM2 Protein")+
  theme(plot.title = element_text(face="bold"),axis.title = element_text(face="bold"))+
  ylim(0,2)

print(WBLV001)

results_noWT_puro$Virus=factor(results_noWT_puro$Virus,levels=c("TREM2","Intergenic"))

WBLV002=ggplot(results_noWT_puro,aes(x=Virus, y=Normalised_to_WT))+
  geom_hline(yintercept = 1,color="grey",linetype="dashed")+
  geom_boxplot(width=0.5,size=0.4,fill=cbPalette[c(5,6)],outlier.shape = NA)+ 
  geom_point(
    aes(color=Line,shape=Repeat),
    size=1,
    position=position_dodge(0.2))+
  geom_line(aes(x=Virus,y=Normalised_to_WT,group=interaction(Repeat,Line),color=Line),position = position_dodge(0.2))+
  stat_compare_means(comparisons = my_comparisons1,label.y=1.6,label="p.signif",size=3)+
  theme_classic(base_size=7)+
  theme(plot.title = element_text(face="bold"),axis.title = element_text(face="bold"))+
  labs(y="Relative TREM2 Expression (to WT)",x="Perturbation",title="TREM2 Protein")+
  ylim(0,2)

print(WBLV002)

western=plot_grid(WBLV001,WBLV002)
plot(western)

ggsave("LV_WB.jpeg",plot=last_plot(),device="jpeg",width=10,height=5)
plot_grid(NULL,western,NULL,phagocytosis,ncol=1,rel_heights = c(1,1,0.5,1),labels="AUTO")
ggsave("~/Paper/Cambridge/TREM2 Lenti Validation Arrayed.jpeg",plot=last_plot(),device="jpeg",height=12,width=8,dpi=600)


