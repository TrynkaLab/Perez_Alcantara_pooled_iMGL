library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(viridis)
library(cowplot)
library(drc)


##load in data

setwd("~/Cambridge/Validation/TREM2 ASO/")

results=read.csv("WB Analysis Test.csv",header=TRUE)
results$meta=paste(results$ASO,results$Concentration,sep="_")
results$Normalised=results$TREM2/results$Viniculin

TREM2_171_control_average=results[results$Concentration==0.000076300
 & results$ASO=="TREM2-171",]
TREM2_171_control_average=mean(TREM2_171_control_average$Normalised)
TREM2_171=results[results$ASO=="TREM2-171",]
TREM2_171$Normalised_Control=TREM2_171$Normalised/TREM2_171_control_average



TREM2_192_control_average=results[results$Concentration==0.000076300
 & results$ASO=="TREM2-192",]
TREM2_192_control_average=mean(TREM2_192_control_average$Normalised)
TREM2_192=results[results$ASO=="TREM2-192",]
TREM2_192$Normalised_Control=TREM2_192$Normalised/TREM2_192_control_average

results=rbind(TREM2_171,TREM2_192)

summary_stats= results %>% 
  group_by(ASO,Concentration,Experiment) %>%
  dplyr::summarise(
    n=n(),
    mean=mean(Normalised_Control),
    sd=sd(Normalised_Control)
  )

## QC Data

summary_stats

QCplot_171=ggplot(TREM2_171,aes(x=Concentration,y=Normalised_Control,group=ASO))+
  geom_point(aes(shape=Blot,color=Experiment))+
  scale_x_log10(breaks=c(0.00001,0.0001,0.001,0.01,0.1,1,10,20),labels=scales::label_number(drop0trailing=TRUE))+
  theme_classic()+
  labs(title="TREM2 ASO-171 QC",y="Relative TREM2 Protein Levels (of WT)",x=expression(paste("ASO Concentration (",mu,"M)",sep="")),fill="Selection")


print(QCplot_171)


QCplot_192=ggplot(TREM2_192,aes(x=Concentration,y=Normalised_Control,group=ASO))+
  geom_point(aes(shape=Blot,color=Experiment))+
  scale_x_log10(breaks=c(0.00001,0.0001,0.001,0.01,0.1,1,10,20),labels=scales::label_number(drop0trailing=TRUE))+
  theme_classic()+
  labs(title="TREM2 ASO-192 QC",y="Relative TREM2 Protein Levels (of WT)",x=expression(paste("ASO Concentration (",mu,"M)",sep="")),fill="Selection")


print(QCplot_192)


### drop ASO_TREM2002 frin ASI TREN2 WB 001

TREM2_171_clean=TREM2_171[!(TREM2_171$Experiment=="ASO_TREM2_002" & TREM2_171$Blot=="ASO TREM2 WB 001"),]
TREM2_192_clean=TREM2_192[!(TREM2_192$Experiment=="ASO_TREM2_002" & TREM2_192$Blot=="ASO TREM2 WB 001"),]



QCplot_171_clean=ggplot(TREM2_171_clean,aes(x=Concentration,y=Normalised_Control,group=ASO))+
  geom_point(aes(shape=Blot,color=Experiment))+
  scale_x_log10(breaks=c(0.00001,0.0001,0.001,0.01,0.1,1,10,20),labels=scales::label_number(drop0trailing=TRUE))+
  theme_classic()+
  labs(title="TREM2 ASO-171 Clean QC",y="Relative TREM2 Protein Levels (of WT)",x=expression(paste("ASO Concentration (",mu,"M)",sep="")),fill="Selection")


print(QCplot_171_clean)


QCplot_192_clean=ggplot(TREM2_192_clean,aes(x=Concentration,y=Normalised_Control,group=ASO))+
  geom_point(aes(shape=Blot,color=Experiment))+
  scale_x_log10(breaks=c(0.00001,0.0001,0.001,0.01,0.1,1,10,20),labels=scales::label_number(drop0trailing=TRUE))+
  theme_classic()+
  labs(title="TREM2 ASO-192 Clean QC",y="Relative TREM2 Protein Levels (of WT)",x=expression(paste("ASO Concentration (",mu,"M)",sep="")),fill="Selection")

plot_grid(QCplot_171,QCplot_171_clean)
plot_grid(QCplot_192,QCplot_192_clean)

results_cleaned=rbind(TREM2_171_clean,TREM2_192_clean)

summary_stats= results_cleaned %>% 
  group_by(ASO,Concentration,Experiment) %>%
  dplyr::summarise(
    n=n(),
    mean=mean(Normalised_Control),
    sd=sd(Normalised_Control)
  )


colnames(summary_stats)=c("ASO","Concentration","Experiment","n","Normalised_Control","SD")

summary_stats= summary_stats %>% 
  group_by(ASO,Concentration) %>%
  dplyr::summarise(
    n=n(),
    mean=mean(Normalised_Control),
    sd=sd(Normalised_Control)
  )


summary_stats_171=summary_stats[summary_stats$ASO=="TREM2-171",]
summary_stats_192=summary_stats[summary_stats$ASO=="TREM2-192",]

### caculate models

fitted_curve_171= drm(formula= mean ~ Concentration ,
                  data= summary_stats_171,
                  fct = LL.4())

E0_171= fitted_curve_171$coefficients[2]
Einf_171= fitted_curve_171$coefficients[3]
EC50_171= fitted_curve_171$coefficients[4]
H_171= -fitted_curve_171$coefficients[1]

IC50_171= (50/(100-50))^(1/-H_171)*EC50_171
IC90_171= (90/(100-90))^(1/-H_171)*EC50_171

Efficacy_metrics= data.frame(c(E0_171,Einf_171,EC50_171,H_171,IC50_171,IC90_171))
rownames(Efficacy_metrics)=c("E0","Einf","EC50","H","IC50","IC90")

fitted_curve_192= drm(formula= mean ~ Concentration ,
                      data= summary_stats_192,
                      fct = LL.4())

E0_192= fitted_curve_192$coefficients[2]
Einf_192= fitted_curve_192$coefficients[3]
EC50_192= fitted_curve_192$coefficients[4]
H_192= -fitted_curve_192$coefficients[1]

IC50_192= (50/(100-50))^(1/-H_192)*EC50_192
IC90_192= (90/(100-90))^(1/-H_192)*EC50_192

Efficacy_metrics$value=c(E0_192,Einf_192,EC50_192,H_192,IC50_192,IC90_192)

colnames(Efficacy_metrics)=c("TREM2-171","TREM2-192")

Efficacy_metrics

## plot data

TREM2_171_DRC=ggplot(summary_stats_171,aes(x=Concentration,y=mean))+
  geom_point(size=1.2)+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2)+
  theme_classic(base_size = 7)+
  geom_smooth(method=drm,method.args=list(fct=L.4()),se=FALSE,color="#f8766d")+
  scale_x_log10(breaks=c(0.00001,0.0001,0.001,0.01,0.1,1,10,20),labels=scales::label_number(drop0trailing=TRUE))+
  ylim(c(-0.05,1.5))+
  theme(plot.title = element_text(face="bold"),axis.title = element_text(face="bold"))+
  labs(title="TREM2-171 ASO Dose Response",y="Relative TREM2 Protein Levels (to Vehicle)",x=expression(paste("Log ASO Concentration (",mu,"M)",sep="")),fill="Selection")

TREM2_192_DRC=ggplot(summary_stats_192,aes(x=Concentration,y=mean))+
  geom_point(size=1.2)+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=.2)+
  theme_classic(base_size=7)+
  geom_smooth(method=drm,method.args=list(fct=L.4()),se=FALSE,color="#00bfc4")+
  scale_x_log10(breaks=c(0.00001,0.0001,0.001,0.01,0.1,1,10,20),labels=scales::label_number(drop0trailing=TRUE))+
  ylim(c(-0.05,1.5))+
  theme(plot.title = element_text(face="bold"),axis.title = element_text(face="bold"))+
  labs(title="TREM2-192 ASO Dose Response",y="Relative TREM2 Protein Levels (to Vehicle)",x=expression(paste("Log ASO Concentration (",mu,"M)",sep="")),fill="Selection")

TREM2_DRC=plot_grid(TREM2_171_DRC,TREM2_192_DRC,nrow=1)
plot(TREM2_DRC)
ggsave("FIGURES/TREM2_ASO_DRC.pdf",plot=TREM2_DRC,device="pdf",dpi=300,units=c("mm"),height=90,width=180)
ggsave("../FIGURES/TREM2_ASO_DRC.jpeg",plot=TREM2_DRC,device="jpeg",dpi=300,units=c("mm"),height=90,width=180)
