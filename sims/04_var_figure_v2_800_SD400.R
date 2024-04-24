#Program: Create Variance Ratio Figure
#Developed by: BES

library(tidyverse)
library(ggplot2)
theme_set(theme_bw(16))
library(RColorBrewer)
library(ggimage)
library(gridExtra)
library(png)
library(grid)
library(lattice)
library(boot)
library(scales)
library(cowplot)
library(forestplot)

#read in truth
true<- as.data.frame(read.csv("ACE_Truth_SD400.csv"))
ACE.true<-true$x

#bring in results and parse by estimator and model specification
res<-as.data.frame(read.csv("results_N800_SD400.csv"))
res.PI.CS<-res[res$type=="PI,CS",]
res.WTD.CS<-res[res$type=="WTD,CS",]
res.TMLE.CS<-res[res$type=="TMLE,CS",]
res.PI.MO<-res[res$type=="PI,MO",]
res.WTD.MO<-res[res$type=="WTD,MO",]
res.TMLE.MO<-res[res$type=="TMLE,MO",]
res.PI.MW<-res[res$type=="PI,MW",]
res.WTD.MW<-res[res$type=="WTD,MW",]
res.TMLE.MW<-res[res$type=="TMLE,MW",]
res.PI.MB<-res[res$type=="PI,MB",]
res.WTD.MB<-res[res$type=="WTD,MB",]
res.TMLE.MB<-res[res$type=="TMLE,MB",]

#calculate metrics for each scenario (ESE,ASE,SER)
calc_table <- function(ACE,ESse,IFse){
  ESE<-sd(na.omit(ACE))
  #sandwich
  ASE.ES<-mean(na.omit(ESse))
  SER.ES<-ASE.ES/ESE
  #IF
  ASE.IF<-mean(na.omit(IFse))
  SER.IF<-ASE.IF/ESE

  results<-c(ESE, ASE.ES, SER.ES, ASE.IF, SER.IF)
}

#calculate metrics for each estimator and specification
PI.bc<-c("PI,CS", "PI", calc_table(res.PI.CS$ATE,res.PI.CS$ESseATE,res.PI.CS$IFseATE))
WTD.bc<-c("WTD,CS", "Wtd Reg",calc_table(res.WTD.CS$ATE,res.WTD.CS$ESseATE,res.WTD.CS$IFseATE))
TMLE.bc<-c("TMLE,CS", "TMLE",calc_table(res.TMLE.CS$ATE,res.TMLE.CS$ESseATE,res.TMLE.CS$IFseATE))
all.cor<-rbind(PI.bc,WTD.bc,TMLE.bc)

PI.mo<-c("PI,MO", "PI",calc_table(res.PI.MO$ATE,res.PI.MO$ESseATE,res.PI.MO$IFseATE))
WTD.mo<-c("WTD,MO", "Wtd Reg",calc_table(res.WTD.MO$ATE,res.WTD.MO$ESseATE,res.WTD.MO$IFseATE))
TMLE.mo<-c("TMLE,MO", "TMLE",calc_table(res.TMLE.MO$ATE,res.TMLE.MO$ESseATE,res.TMLE.MO$IFseATE))
all.mo<-rbind(PI.mo,WTD.mo,TMLE.mo)

PI.MW<-c("PI,MW", "PI",calc_table(res.PI.MW$ATE,res.PI.MW$ESseATE,res.PI.MW$IFseATE))
WTD.MW<-c("WTD,MW", "Wtd Reg",calc_table(res.WTD.MW$ATE,res.WTD.MW$ESseATE,res.WTD.MW$IFseATE))
TMLE.MW<-c("TMLE,MW", "TMLE",calc_table(res.TMLE.MW$ATE,res.TMLE.MW$ESseATE,res.TMLE.MW$IFseATE))
all.mw<-rbind(PI.MW,WTD.MW,TMLE.MW)

PI.MB<-c("PI,MB", "PI",calc_table(res.PI.MB$ATE,res.PI.MB$ESseATE,res.PI.MB$IFseATE))
WTD.MB<-c("WTD,MB", "Wtd Reg",calc_table(res.WTD.MB$ATE,res.WTD.MB$ESseATE,res.WTD.MB$IFseATE))
TMLE.MB<-c("TMLE,MB", "TMLE",calc_table(res.TMLE.MB$ATE,res.TMLE.MB$ESseATE,res.TMLE.MB$IFseATE))
all.mb<-rbind(PI.MB,WTD.MB,TMLE.MB)

#combine and format
all<-as.data.frame(rbind(all.cor,all.mo,all.mw,all.mb))
colnames(all)<-c('type','estname','ESE', 'ASE.ES', 'SER.ES', 'ASE.IF', 'SER.IF')

data<-merge(res, all, by="type", all.x=TRUE)
data$sd_ratio_ES<-data$ESseATE/as.numeric(data$ESE)
data$sd_ratio_IF<-data$IFseATE/as.numeric(data$ESE)

sand<-data[,c('type','estname','ESseATE','ESE','ASE.ES','SER.ES')]
sand$EST<-"Empirical Sandwich"
colnames(sand)<-c('type','estname','seATE','ESE','ASE','SER','EST')

IF<-data[,c('type','estname','IFseATE','ESE','ASE.IF','SER.IF')]
IF$EST<-"Influence Function"
colnames(IF)<-c('type','estname','seATE','ESE','ASE','SER','EST')

stack<-as.data.frame(rbind(sand,IF))
stack$std_ratio<-stack$seATE/as.numeric(stack$ESE)
stack$Scenario<-factor(data$estname,levels=c('PI','Wtd Reg','TMLE'))
stack$EST2<-as.factor(stack$EST)

CS<-stack[stack$type=='PI,CS'|stack$type=='WTD,CS'|stack$type=='TMLE,CS',]
MO<-stack[stack$type=='PI,MO'|stack$type=='WTD,MO'|stack$type=='TMLE,MO',]
MW<-stack[stack$type=='PI,MW'|stack$type=='WTD,MW'|stack$type=='TMLE,MW',]
MB<-stack[stack$type=='PI,MB'|stack$type=='WTD,MB'|stack$type=='TMLE,MB',]

#create plots for each scenario
var.plot.CS<-ggplot(CS,aes(x=Scenario, y=std_ratio, fill=EST2))+ scale_y_continuous(limits=c(0.5,11),breaks=c(0.5,0.75,0.8,0.9,1,1.1,1.2,1.25,1.5,1.75,2,2.5,3))+geom_violin() + stat_summary(fun = "mean", geom = "crossbar", width = 0.1, color = "black",position = position_dodge(width = .9))+ scale_fill_manual(values = c("#FF3D3D","#57B0FF"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="bottom",legend.title = element_blank())+xlab("Both Correct")+ylab("")+geom_vline(xintercept=c(1.5,2.5), colour='light grey',lwd=1)+ geom_hline(yintercept=1.0, linetype="dashed", color = "black")
#+annotate("text", x=.9, y=-20, label= paste(c("SER:     ES=", round(as.numeric(all[all$type=='PI,CS','SER.ES']),2),"         IF=", round(as.numeric(all[all$type=='PI,CS','SER.IF']),2)), collapse = " "))+annotate("text", x=2.9, y=-20, label= paste(c("         ES=", round(as.numeric(all[all$type=='TMLE,CS','SER.ES']),2),"         IF=", round(as.numeric(all[all$type=='TMLE,CS','SER.IF']),2)), collapse = " "))+annotate("text", x=1.9, y=-20, label= paste(c("         ES=", round(as.numeric(all[all$type=='WTD,CS','SER.ES']),2),"         IF=", round(as.numeric(all[all$type=='WTD,CS','SER.IF']),2)), collapse = " "))
var.plot.CS2<-var.plot.CS+coord_cartesian(ylim = c(0.8, 1.2)) 

var.plot.MO<-ggplot(MO,aes(x=Scenario, y=std_ratio, fill=EST2))+ scale_y_continuous(limits=c(0.5,11),breaks=c(0.5,0.75,0.8,0.9,1,1.1,1.2,1.25,1.5,1.75,2,2.5,3))+geom_violin() + stat_summary(fun = "mean", geom = "crossbar", width = 0.1, color = "black",position = position_dodge(width = .9))+ scale_fill_manual(values = c("#FF3D3D","#57B0FF"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="bottom",legend.title = element_blank())+xlab("Misspecified Outcome")+ylab("")+geom_vline(xintercept=c(1.5,2.5), colour='light grey',lwd=1)+ geom_hline(yintercept=1.0, linetype="dashed", color = "black")
#+annotate("text", x=.9, y=-20, label= paste(c("SER:     ES=", round(as.numeric(all[all$type=='PI,MO','SER.ES']),2),"         IF=", round(as.numeric(all[all$type=='PI,MO','SER.IF']),2)), collapse = " "))+annotate("text", x=2.9, y=-20, label= paste(c("         ES=", round(as.numeric(all[all$type=='TMLE,MO','SER.ES']),2),"         IF=", round(as.numeric(all[all$type=='TMLE,MO','SER.IF']),2)), collapse = " "))+annotate("text", x=1.9, y=-20, label= paste(c("         ES=", round(as.numeric(all[all$type=='WTD,MO','SER.ES']),2),"         IF=", round(as.numeric(all[all$type=='WTD,MO','SER.IF']),2)), collapse = " "))
var.plot.MO2<-var.plot.MO+coord_cartesian(ylim = c(0.8, 1.2)) 

var.plot.MW<-ggplot(MW,aes(x=Scenario, y=std_ratio, fill=EST2))+ scale_y_continuous(limits=c(0.5,11),breaks=c(0.5,0.75,0.8,0.9,1,1.1,1.2,1.25,1.5,1.75,2,2.5,3))+geom_violin() + stat_summary(fun = "mean", geom = "crossbar", width = 0.1, color = "black",position = position_dodge(width = .9))+ scale_fill_manual(values = c("#FF3D3D","#57B0FF"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="bottom",legend.title = element_blank())+xlab("Misspecified Propensity")+ylab("")+geom_vline(xintercept=c(1.5,2.5), colour='light grey',lwd=1)+ geom_hline(yintercept=1.0, linetype="dashed", color = "black")
#+annotate("text", x=.9, y=-20, label= paste(c("SER:     ES=", round(as.numeric(all[all$type=='PI,MW','SER.ES']),2),"         IF=", round(as.numeric(all[all$type=='PI,MW','SER.IF']),2)), collapse = " "))+annotate("text", x=2.9, y=-20, label= paste(c("         ES=", round(as.numeric(all[all$type=='TMLE,MW','SER.ES']),2),"         IF=", round(as.numeric(all[all$type=='TMLE,MW','SER.IF']),2)), collapse = " "))+annotate("text", x=1.9, y=-20, label= paste(c("         ES=", round(as.numeric(all[all$type=='WTD,MW','SER.ES']),2),"         IF=", round(as.numeric(all[all$type=='WTD,MW','SER.IF']),2)), collapse = " "))
var.plot.MW2<-var.plot.MW+coord_cartesian(ylim = c(0.8, 1.2)) 

var.plot.MB<-ggplot(MB,aes(x=Scenario, y=std_ratio, fill=EST2))+ scale_y_continuous(limits=c(0.5,11),breaks=c(0.5,0.75,0.8,0.9,1,1.1,1.2,1.25,1.5,1.75,2,2.5,3))+geom_violin() + stat_summary(fun = "mean", geom = "crossbar", width = 0.1, color = "black",position = position_dodge(width = .9))+ scale_fill_manual(values = c("#FF3D3D","#57B0FF"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="bottom",legend.title = element_blank())+xlab("Both Misspecified")+ylab("")+geom_vline(xintercept=c(1.5,2.5), colour='light grey',lwd=1)+ geom_hline(yintercept=1.0, linetype="dashed", color = "black")
#+annotate("text", x=.9, y=-20, label= paste(c("SER:     ES=", round(as.numeric(all[all$type=='PI,MB','SER.ES']),2),"         IF=", round(as.numeric(all[all$type=='PI,MB','SER.IF']),2)), collapse = " "))+annotate("text", x=2.9, y=-20, label= paste(c("         ES=", round(as.numeric(all[all$type=='TMLE,MB','SER.ES']),2),"         IF=", round(as.numeric(all[all$type=='TMLE,MB','SER.IF']),2)), collapse = " "))+annotate("text", x=1.9, y=-20, label= paste(c("         ES=", round(as.numeric(all[all$type=='WTD,MB','SER.ES']),2),"         IF=", round(as.numeric(all[all$type=='WTD,MB','SER.IF']),2)), collapse = " "))
var.plot.MB2<-var.plot.MB+coord_cartesian(ylim = c(0.8, 1.2)) 

#panel and output
plots<-plot_grid(var.plot.CS2+theme(legend.position = "none",text = element_text(size=12)),
                 var.plot.MO2+theme(legend.position = "none",text = element_text(size=12)),
                 var.plot.MW2+theme(legend.position = "none",text = element_text(size=12)),
                 var.plot.MB2+theme(legend.position = "none",text = element_text(size=12)),
                 nrow=2,ncol=2)

var.plot.CS<-ggplot(CS,aes(x=Scenario, y=std_ratio, fill=EST2))+ scale_y_continuous(limits=c(0.5,11),breaks=c(0.5,0.75,1,1.25,1.5,1.75,2,2.5,3))+geom_violin() + scale_fill_manual(values = c("#FF3D3D","#57B0FF"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="bottom",legend.title = element_blank())+xlab("Both Correct")+ylab("")+geom_vline(xintercept=c(1.5,2.5), colour='light grey',lwd=1)+ geom_hline(yintercept=1.0, linetype="dashed", color = "black")
legend <- get_legend(var.plot.CS)


plot_grid(plots, legend, ncol = 1, rel_heights = c(1, .1))
                 

