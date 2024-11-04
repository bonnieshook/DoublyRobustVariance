#Program: Create summary table based on simulation results
#Developed by: BES

#read in true ACE for bias calculations
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

#calculate metrics for summary table (bias, ASE, ESE, SER, coverage)
calc_table <- function(ACE,ESse,IFse){
  bias<-mean(ACE, na.rm = TRUE)-ACE.true
  relbias<-100*abs(bias)/abs(ACE.true)
  ESE<-sd(na.omit(ACE))
  #sandwich
  ll.ES <- na.omit(ACE)-1.96*na.omit(ESse)
  ul.ES <- na.omit(ACE)+1.96*na.omit(ESse)
  cov.ES <-100*mean((ll.ES <= ACE.true & ul.ES >= ACE.true))
  ASE.ES<-mean(na.omit(ESse))
  SER.ES<-ASE.ES/ESE
  #IF
  ll.IF <- na.omit(ACE)-1.96*na.omit(IFse)
  ul.IF <- na.omit(ACE)+1.96*na.omit(IFse)
  cov.IF <-100*mean((ll.IF <= ACE.true & ul.IF >= ACE.true))
  ASE.IF<-mean(na.omit(IFse))
  SER.IF<-ASE.IF/ESE
  #bootstrap
  ll.boot <- na.omit(ACE)-1.96*na.omit(bootse)
  ul.boot <- na.omit(ACE)+1.96*na.omit(bootse)
  cov.boot <-100*mean((ll.boot <= ACE.true & ul.boot >= ACE.true))
  ASE.boot<-mean(na.omit(bootse))
  SER.boot<-ASE.boot/ESE
  
  results<-c(round(bias,digits=1),"&", round(ESE, digits=1),"&", round(SER.ES, digits=2),"&", round(cov.ES, digits=0),"&", round(SER.boot, digits=2),"&", round(cov.boot, digits=0),"&", round(SER.IF, digits=2),"&", round(cov.IF, digits=0),'\\')
}

#format results for LaTeX table
PI.bc<-c("CS", "&", "PI", "&", calc_table(res.PI.CS$ATE,res.PI.CS$ESseATE,res.PI.CS$IFseATE,res.PI.CS$bootseATE))
WTD.bc<-c("", "&","WTD", "&", calc_table(res.WTD.CS$ATE,res.WTD.CS$ESseATE,res.WTD.CS$IFseATE,res.WTD.CS$bootseATE))
TMLE.bc<-c("", "&","TMLE", "&", calc_table(res.TMLE.CS$ATE,res.TMLE.CS$ESseATE,res.TMLE.CS$IFseATE,res.TMLE.CS$bootseATE))
all.cor<-rbind(PI.bc,WTD.bc,TMLE.bc)

PI.mo<-c("MO", "&","PI", "&", calc_table(res.PI.MO$ATE,res.PI.MO$ESseATE,res.PI.MO$IFseATE,res.PI.MO$bootseATE))
WTD.mo<-c("", "&", "WTD", "&", calc_table(res.WTD.MO$ATE,res.WTD.MO$ESseATE,res.WTD.MO$IFseATE,res.WTD.MO$bootseATE))
TMLE.mo<-c("","&", "TMLE", "&", calc_table(res.TMLE.MO$ATE,res.TMLE.MO$ESseATE,res.TMLE.MO$IFseATE,res.TMLE.MO$bootseATE))
all.mo<-rbind(PI.mo,WTD.mo,TMLE.mo)

PI.MW<-c("MW", "&","PI", "&", calc_table(res.PI.MW$ATE,res.PI.MW$ESseATE,res.PI.MW$IFseATE,res.PI.MW$bootseATE))
WTD.MW<-c("", "&", "WTD", "&", calc_table(res.WTD.MW$ATE,res.WTD.MW$ESseATE,res.WTD.MW$IFseATE,res.WTD.MW$bootseATE))
TMLE.MW<-c("", "&", "TMLE", "&", calc_table(res.TMLE.MW$ATE,res.TMLE.MW$ESseATE,res.TMLE.MW$IFseATE,res.TMLE.MW$bootseATE))
all.mw<-rbind(PI.MW,WTD.MW,TMLE.MW)

PI.MB<-c("MB", "&","PI", "&", calc_table(res.PI.MB$ATE,res.PI.MB$ESseATE,res.PI.MB$IFseATE,res.PI.MB$bootseATE))
WTD.MB<-c("", "&", "WTD", "&", calc_table(res.WTD.MB$ATE,res.WTD.MB$ESseATE,res.WTD.MB$IFseATE,res.WTD.MB$bootseATE))
TMLE.MB<-c("", "&", "TMLE", "&", calc_table(res.TMLE.MB$ATE,res.TMLE.MB$ESseATE,res.TMLE.MB$IFseATE,res.TMLE.MB$bootseATE))
all.mb<-rbind(PI.MB,WTD.MB,TMLE.MB)


#output
all<-rbind(all.cor,all.mo,all.mw,all.mb)
write.csv(all, "SimSummary_n800_SD400.csv")








