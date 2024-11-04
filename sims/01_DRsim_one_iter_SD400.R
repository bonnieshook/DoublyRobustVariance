#Program: 01_DRsim_one_iter_SD400.R
#Purpose: This program generates data under the DGM and analyses it using the 3 point and 2 variance estimators (for a single iteration of the sims)

library(geex)
library(boot)
library(stats)
library(resample)

## collect arguments passed in from SLURM job script
args <- commandArgs(trailingOnly=TRUE)
subdir_name <- as.character(args[1])
num.part <- as.numeric(args[2])
num.boots<-as.numeric(args[3])
sigma<-as.numeric(args[4])


user_home_dir <- "<path>" # top-level directory for the project
dir_path <- paste(user_home_dir,subdir_name,"/",sep="")

sim <- Sys.getenv("SLURM_ARRAY_TASK_ID") 
set.seed(sim)

setwd(user_home_dir)
source('00_Functions.R')
  part.id<-seq(1:num.part)
  
  #covariates
  Z1<-as.vector(rnorm(num.part, mean =155, sd = 7.6)) 
  simdat1<-as.data.frame(cbind(part.id,Z1))
  simdat1$Z1_sq<-(simdat1$Z1-155)^2
  simdat1$Z2<-rbinom(num.part, 1, .25 ) 
  simdat1$Z3<-rbinom(num.part, 1, .75 ) 

  #generate exposure and potential outcomes
  simdat1$e = inv.logit(15-0.1*simdat1$Z1+2.5*simdat1$Z2 - 0.02*simdat1$Z1*simdat1$Z2 - 1*simdat1$Z3 + 0.005*simdat1$Z1*simdat1$Z3)
  simdat1$X<-rbinom(num.part,1,simdat1$e) #anemia
  simdat1$EY0 = 1000+11.5*simdat1$Z1 + 100*simdat1$Z2  - 15*simdat1$Z1*simdat1$Z2   #birth weight under no anemia
  simdat1$EY1 = 1025+11.5*simdat1$Z1 + 100*simdat1$Z2 - 5.5*simdat1$Z1 - 30*simdat1$Z2 + 5*simdat1$Z1*simdat1$Z2  #birth weight under anemia
  simdat1$Y0<- rnorm(num.part, mean =simdat1$EY0, sd=sigma)
  simdat1$Y1<- rnorm(num.part, mean =simdat1$EY1, sd=sigma)
  
  #assign observed outcome
  simdat1$Y<-ifelse(simdat1$X==0,simdat1$Y0,simdat1$Y1)
  
  #create dataset with just observed data
  simdat1.obs<-simdat1[,c('part.id','Z1','Z1_sq','Z2','Z3','X','Y')]
  
  
  #############################################################################################################################################
  ##### Compute each of the three DR estimators on observed data, with 3 variance estimates (influence function and, sandwich, bootstrap)#####
  #############################################################################################################################################
  
  #########################################################
  ##### fit weight models (needed for all estimators) #####
  #########################################################
  
  #correct specification
  prop.model <- glm(X ~ Z1 + Z2 + Z1:Z2 + Z3 + Z1:Z3, family = binomial(), data = simdat1)
  simdat1$PS<-predict(prop.model,simdat1, type="response")
  simdat1$IPTW<-ifelse(simdat1$X==1,simdat1$PS^(-1),(1-simdat1$PS)^(-1))
  
  #misspecified weight model 
  prop.model.mis <- glm(X ~ Z1_sq, family = binomial(), data = simdat1)
  simdat1$PS.mis<-predict(prop.model.mis,simdat1, type="response")
  simdat1$IPTW.mis<-ifelse(simdat1$X==1,simdat1$PS.mis^(-1),(1-simdat1$PS.mis)^(-1))
  
  #create copies of data sets for estimated potential outcomes
  alltrt<-simdat1
  alluntrt<-simdat1
  alltrt$X<-1
  alluntrt$X<-0
  
  ########################################################
  ################Correctly Specified Models##############
  ########################################################
  
  ############################
  #classic (plug-in) AIPW#
  ############################
  
  #fit outcome model (unweighted)
  out.model.unwt <- glm(Y ~ Z1*Z2*X, data = simdat1)
  simdat1$Y0.unwt<-(predict(out.model.unwt, alluntrt, type = "response"))
  simdat1$Y1.unwt<-(predict(out.model.unwt, alltrt, type = "response"))
  
  #get pseudo-potential outcomes
  CM1.PI<-mean(((simdat1$X*simdat1$Y-(simdat1$X-simdat1$PS)*simdat1$Y1.unwt))/simdat1$PS)
  CM0.PI<-mean((((1-simdat1$X)*simdat1$Y+(simdat1$X-simdat1$PS)*simdat1$Y0.unwt))/(1-simdat1$PS))
  ATE.PI<-CM1.PI - CM0.PI
  
  #get 3 standard errors
  DR.PI<-geex_PI(simdat1, X ~ Z1 + Z2 + Z1:Z2 + Z3 + Z1:Z3, Y ~ X*Z1*Z2, CM0.PI, CM1.PI, ATE.PI)
  IFseATE.PI<-IF_Var(simdat1,'PI',ATE.PI)
  bootseATE.PI<-bootstrap_PIAIPW(num.boots, simdat1.obs, X ~ Z1 + Z2 + Z1:Z2 + Z3 + Z1:Z3, Y ~ X*Z1*Z2)
  
  #compile results all 3 variance estimators
  DR.PI2<-as.data.frame(cbind(DR.PI,IFseATE.PI,bootseATE.PI))
  colnames(DR.PI2)<-c('ATE', 'ESseATE','IFseATE','bootseATE')
  DR.PI2$type<-'PI,CS'
  
  
  ##########################
  #weighted regression AIPW#
  ##########################
  
  #fit outcome model (weighted)
  out.model.wt <- glm(Y ~ Z1*Z2*X, weights=IPTW, data = simdat1)
  simdat1$Y0.wt<-(predict(out.model.wt, alluntrt, type = "response"))
  simdat1$Y1.wt<-(predict(out.model.wt, alltrt, type = "response"))
  
  #estimate causal means
  CM0.WTD<-mean(simdat1$Y0.wt)
  CM1.WTD<-mean(simdat1$Y1.wt)
  ATE.WTD<-CM1.WTD - CM0.WTD
  
  #get 3 standard errors
  DR.WTD<-geex_WTD(simdat1, X ~ Z1 + Z2 + Z1:Z2 + Z3 + Z1:Z3, Y ~ Z1*Z2*X, coef(out.model.wt), CM0.WTD, CM1.WTD, ATE.WTD)
  IFseATE.WTD<-IF_Var(simdat1,'WTD',ATE.WTD)
  bootseATE.WTD<-bootstrap_WRAIPW(num.boots, simdat1.obs, X ~ Z1 + Z2 + Z1:Z2 + Z3 + Z1:Z3, Y ~ Z1*Z2*X)
  
  #compile results all 3 variance estimators
  DR.WTD2<-as.data.frame(cbind(DR.WTD,IFseATE.WTD,bootseATE.WTD))
  colnames(DR.WTD2) <-c('ATE', 'ESseATE','IFseATE','bootseATE')
  DR.WTD2$type<-'WTD,CS'
  
  
  ######
  #TMLE#
  ######
  
  #scale the outcome for TMLE
  a<-min(simdat1$Y)
  b<-max(simdat1$Y)
  simdat1$a<-a
  simdat1$b<-b
  simdat1$Y_scaled<-(simdat1$Y-a)/(b-a)

  #fit scaled outcome model (unweighted)
  out.model.sc.unwt <- glm(Y_scaled ~ Z1*Z2*X, data = simdat1)
  simdat1$Y0.sc.unwt<-(predict(out.model.sc.unwt, alluntrt, type = "response"))
  simdat1$Y1.sc.unwt<-(predict(out.model.sc.unwt, alltrt, type = "response"))
  
  #truncate predictions to within the range of the data
  simdat1$Y0.sc.unwt<-ifelse(simdat1$Y0.sc.unwt<min(simdat1$Y_scaled),min(simdat1$Y_scaled),simdat1$Y0.sc.unwt)
  simdat1$Y0.sc.unwt<-ifelse(simdat1$Y0.sc.unwt>max(simdat1$Y_scaled),max(simdat1$Y_scaled),simdat1$Y0.sc.unwt)
  
  simdat1$Y1.sc.unwt<-ifelse(simdat1$Y1.sc.unwt<min(simdat1$Y_scaled),min(simdat1$Y_scaled),simdat1$Y1.sc.unwt)
  simdat1$Y1.sc.unwt<-ifelse(simdat1$Y1.sc.unwt>max(simdat1$Y_scaled),max(simdat1$Y_scaled),simdat1$Y1.sc.unwt)
  
  #fit TMLE targeting model
  target.model.Y0 <- glm(Y_scaled ~ 1, offset=logit(Y0.sc.unwt), weights = (1-X)/(1-PS), 
                         family = binomial(), data = simdat1)
  target.model.Y1 <- glm(Y_scaled ~ 1, offset=logit(Y1.sc.unwt), weights = X/PS, 
                         family = binomial(), data = simdat1)
  
  simdat1$Y0.TMLE<-inv.logit(logit(simdat1$Y0.sc.unwt)+coef(target.model.Y0)[1])*(b-a)+a
  simdat1$Y1.TMLE<-inv.logit(logit(simdat1$Y1.sc.unwt)+coef(target.model.Y1)[1])*(b-a)+a
  CM0.TMLE<-mean(simdat1$Y0.TMLE)
  CM1.TMLE<-mean(simdat1$Y1.TMLE)
  ATE.TMLE<-CM1.TMLE-CM0.TMLE
  
  #get 3 standard errors
  DR.TMLE<-geex_TMLE(simdat1, X ~ Z1 + Z2 + Z1:Z2 + Z3 + Z1:Z3, Y_scaled ~ Z1*Z2*X, coef(target.model.Y0), coef(target.model.Y1), CM0.TMLE, CM1.TMLE, ATE.TMLE)
  IFseATE.TMLE<-IF_Var(simdat1,'TMLE',ATE.TMLE)
  bootseATE.TMLE<-bootstrap_TMLE(num.boots, simdat1.obs, X ~ Z1 + Z2 + Z1:Z2 + Z3 + Z1:Z3, Y_scaled ~ Z1*Z2*X)
  
  #compile results all 3 variance estimators
  DR.TMLE2<-as.data.frame(cbind(DR.TMLE,IFseATE.TMLE,bootseATE.TMLE))
  colnames(DR.TMLE2)<-c('ATE', 'ESseATE','IFseATE','bootseATE')
  DR.TMLE2$type<-'TMLE,CS'
  
  
  ########################################################
  ###########Estimation-Incorrect outcome model###########
  ########################################################
  
  ############################
  #classic (plug-in) AIPW#
  ############################
  
  #fit outcome model (unweighted) 
  out.model.unwt.mis <- glm(Y ~ X + Z1_sq, data = simdat1)
  simdat1$Y0.unwt.mis<-(predict(out.model.unwt.mis, alluntrt, type = "response"))
  simdat1$Y1.unwt.mis<-(predict(out.model.unwt.mis, alltrt, type = "response"))
  
  #get pseudo-potential outcomes
  CM1.PI.mis<-mean(((simdat1$X*simdat1$Y-(simdat1$X-simdat1$PS)*simdat1$Y1.unwt.mis))/simdat1$PS)
  CM0.PI.mis<-mean((((1-simdat1$X)*simdat1$Y+(simdat1$X-simdat1$PS)*simdat1$Y0.unwt.mis))/(1-simdat1$PS))
  ATE.PI.mis<-CM1.PI.mis - CM0.PI.mis
  
  #get 3 standard errors
  DR.PI.mis<-geex_PI(simdat1, X ~ Z1 + Z2 + Z1:Z2 + Z3 + Z1:Z3, Y ~ X+Z1_sq, CM0.PI.mis, CM1.PI.mis, ATE.PI.mis)
  IFseATE.PI.mis<-IF_Var_mis(simdat1,'PI',ATE.PI.mis)
  bootseATE.PI.mis<-bootstrap_PIAIPW(num.boots, simdat1.obs, X ~ Z1 + Z2 + Z1:Z2 + Z3 + Z1:Z3, Y ~ X+Z1_sq)
  
  #compile results all 3 variance estimators
  DR.PI.mis2<-as.data.frame(cbind(DR.PI.mis,IFseATE.PI.mis,bootseATE.PI.mis))
  colnames(DR.PI.mis2)<-c('ATE', 'ESseATE','IFseATE','bootseATE')
  DR.PI.mis2$type<-'PI,MO'
  
  
  ##########################
  #weighted regression AIPW#
  ##########################
  
  #fit outcome model (weighted)
  out.model.wt.mis <- glm(Y ~ X+Z1_sq, weights=IPTW, data = simdat1)
  simdat1$Y0.wt.mis<-(predict(out.model.wt.mis, alluntrt, type = "response"))
  simdat1$Y1.wt.mis<-(predict(out.model.wt.mis, alltrt, type = "response"))
  
  #estimate causal means
  CM0.WTD.mis<-mean(simdat1$Y0.wt.mis)
  CM1.WTD.mis<-mean(simdat1$Y1.wt.mis)
  ATE.WTD.mis<-CM1.WTD.mis - CM0.WTD.mis
  
  #get 3 standard errors
  DR.WTD.mis<-geex_WTD_mis(simdat1, X ~ Z1 + Z2 + Z1:Z2 + Z3 + Z1:Z3, Y ~ X+Z1_sq, coef(out.model.wt.mis), CM0.WTD.mis, CM1.WTD.mis, ATE.WTD.mis)
  IFseATE.WTD.mis<-IF_Var_mis(simdat1,'WTD',ATE.WTD.mis)
  bootseATE.WTD.mis<-bootstrap_WRAIPW(num.boots, simdat1.obs, X ~ Z1 + Z2 + Z1:Z2 + Z3 + Z1:Z3, Y ~ X+Z1_sq)
  
  #compile results all 3 variance estimators
  DR.WTD.mis2<-as.data.frame(cbind(DR.WTD.mis,IFseATE.WTD.mis,bootseATE.WTD.mis))
  colnames(DR.WTD.mis2)<-c('ATE', 'ESseATE','IFseATE','bootseATE')
  DR.WTD.mis2$type<-'WTD,MO'

  
  ######
  #TMLE#
  ######
  
  #fit scaled outcome model (unweighted)
  out.model.sc.unwt.mis <- glm(Y_scaled ~ X+Z1_sq, data = simdat1)
  simdat1$Y0.sc.unwt.mis<-(predict(out.model.sc.unwt.mis, alluntrt, type = "response"))
  simdat1$Y1.sc.unwt.mis<-(predict(out.model.sc.unwt.mis, alltrt, type = "response"))
  
  #truncate predictions to within the range of the data
  simdat1$Y0.sc.unwt.mis<-ifelse(simdat1$Y0.sc.unwt.mis<min(simdat1$Y_scaled),min(simdat1$Y_scaled),simdat1$Y0.sc.unwt.mis)
  simdat1$Y0.sc.unwt.mis<-ifelse(simdat1$Y0.sc.unwt.mis>max(simdat1$Y_scaled),max(simdat1$Y_scaled),simdat1$Y0.sc.unwt.mis)
  
  simdat1$Y1.sc.unwt.mis<-ifelse(simdat1$Y1.sc.unwt.mis<min(simdat1$Y_scaled),min(simdat1$Y_scaled),simdat1$Y1.sc.unwt.mis)
  simdat1$Y1.sc.unwt.mis<-ifelse(simdat1$Y1.sc.unwt.mis>max(simdat1$Y_scaled),max(simdat1$Y_scaled),simdat1$Y1.sc.unwt.mis)
  
  
  #fit targeting models (incorrect outcome model)
  target.model.Y0.mis <- glm(Y_scaled ~ 1, offset=logit(Y0.sc.unwt.mis), weights = (1-X)/(1-PS), 
                             family = binomial(), data = simdat1)
  target.model.Y1.mis <- glm(Y_scaled ~ 1, offset=logit(Y1.sc.unwt.mis), weights = X/PS, 
                             family = binomial(), data = simdat1)
  
  #transform back to range of the data and estimate causal means and ATE
  simdat1$Y0.TMLE.mis<-inv.logit(logit(simdat1$Y0.sc.unwt.mis)+coef(target.model.Y0.mis)[1])*(b-a)+a
  simdat1$Y1.TMLE.mis<-inv.logit(logit(simdat1$Y1.sc.unwt.mis)+coef(target.model.Y1.mis)[1])*(b-a)+a
  CM0.TMLE.mis<-mean(simdat1$Y0.TMLE.mis)
  CM1.TMLE.mis<-mean(simdat1$Y1.TMLE.mis)
  ATE.TMLE.mis<-CM1.TMLE.mis-CM0.TMLE.mis
  
  #get 3 standard errors
  DR.TMLE.mis<-geex_TMLE(simdat1, X ~ Z1 + Z2 + Z1:Z2 + Z3 + Z1:Z3, Y_scaled ~ X+Z1_sq, coef(target.model.Y0.mis), coef(target.model.Y1.mis), CM0.TMLE.mis, CM1.TMLE.mis, ATE.TMLE.mis)
  IFseATE.TMLE.mis<-IF_Var_mis(simdat1,'TMLE',ATE.TMLE.mis)
  bootseATE.TMLE.mis<-bootstrap_TMLE(num.boots, simdat1.obs, X ~ Z1 + Z2 + Z1:Z2 + Z3 + Z1:Z3, Y_scaled ~ X+Z1_sq)
  
  #compile results all 3 variance estimators
  DR.TMLE.mis2<-as.data.frame(cbind(DR.TMLE.mis,IFseATE.TMLE.mis,bootseATE.TMLE.mis))
  colnames(DR.TMLE.mis2)<-c('ATE', 'ESseATE','IFseATE','bootseATE')
  DR.TMLE.mis2$type<-'TMLE,MO'

  
  
  ########################################################
  ###########Estimation-Incorrect weight model###########
  ########################################################
  
  ############################
  #classic (plug-in) AIPW#
  ############################

  #estimate causal means
  CM1.PI.mw<-mean(((simdat1$X*simdat1$Y-(simdat1$X-simdat1$PS.mis)*simdat1$Y1.unwt))/simdat1$PS.mis)
  CM0.PI.mw<-mean((((1-simdat1$X)*simdat1$Y+(simdat1$X-simdat1$PS.mis)*simdat1$Y0.unwt))/(1-simdat1$PS.mis))
  ATE.PI.mw<-CM1.PI.mw - CM0.PI.mw
  
  #get 3 standard errors
  DR.PI.mw<-geex_PI(simdat1, X ~ Z1_sq, Y ~ X*Z1*Z2, CM0.PI.mw, CM1.PI.mw, ATE.PI.mw)
  IFseATE.PI.mw<-IF_Var_mw(simdat1,'PI',ATE.PI.mw)
  bootseATE.PI.mw<-bootstrap_PIAIPW(num.boots, simdat1.obs, X ~ Z1_sq, Y ~ X*Z1*Z2)
  
  #compile results all 3 variance estimators
  DR.PI.mw2<-as.data.frame(cbind(DR.PI.mw,IFseATE.PI.mw,bootseATE.PI.mw))
  colnames(DR.PI.mw2)<-c('ATE', 'ESseATE','IFseATE','bootseATE')
  DR.PI.mw2$type<-'PI,MW'
  
  

  ##########################
  #weighted regression AIPW#
  ##########################
  
  #fit outcome model (correct model but wrong weights) and get pseudo POs
  out.model.wt.mw <- glm(Y ~ X*Z1*Z2, weights=IPTW.mis, data = simdat1)
  simdat1$Y0.wt.mw<-(predict(out.model.wt.mw, alluntrt, type = "response"))
  simdat1$Y1.wt.mw<-(predict(out.model.wt.mw, alltrt, type = "response"))
  
  #estimate causal means
  CM0.WTD.mw<-mean(simdat1$Y0.wt.mw)
  CM1.WTD.mw<-mean(simdat1$Y1.wt.mw)
  ATE.WTD.mw<-CM1.WTD.mw - CM0.WTD.mw
  
  #get 3 standard errors
  DR.WTD.mw<-geex_WTD(simdat1, X ~ Z1_sq, Y ~ X*Z1*Z2, coef(out.model.wt.mw), CM0.WTD.mw, CM1.WTD.mw, ATE.WTD.mw)
  IFseATE.WTD.mw<-IF_Var_mw(simdat1,'WTD',ATE.WTD.mw)
  bootseATE.WTD.mw<-bootstrap_WRAIPW(num.boots, simdat1.obs, X ~ Z1_sq, Y ~ X*Z1*Z2)
  
  #compile results all 3 variance estimators
  DR.WTD.mw2<-as.data.frame(cbind(DR.WTD.mw,IFseATE.WTD.mw,bootseATE.WTD.mw))
  colnames(DR.WTD.mw2)<-c('ATE', 'ESseATE','IFseATE','bootseATE')
  DR.WTD.mw2$type<-'WTD,MW'
  
  
  ######
  #TMLE#
  ######
  
  #fit targeting models (incorrect weight model)
  target.model.Y0.mw <- glm(Y_scaled ~ 1, offset=logit(Y0.sc.unwt), weights = (1-X)/(1-PS.mis), 
                            family = binomial(), data = simdat1)
  target.model.Y1.mw <- glm(Y_scaled ~ 1, offset=logit(Y1.sc.unwt), weights = X/PS.mis, 
                            family = binomial(), data = simdat1)
  
  #estimate causal means and ATE
  simdat1$Y0.TMLE.mw<-inv.logit(logit(simdat1$Y0.sc.unwt)+coef(target.model.Y0.mw)[1])*(b-a)+a
  simdat1$Y1.TMLE.mw<-inv.logit(logit(simdat1$Y1.sc.unwt)+coef(target.model.Y1.mw)[1])*(b-a)+a
  CM0.TMLE.mw<-mean(simdat1$Y0.TMLE.mw)
  CM1.TMLE.mw<-mean(simdat1$Y1.TMLE.mw)
  ATE.TMLE.mw<-CM1.TMLE.mw-CM0.TMLE.mw
  
  #get 3 standard errors
  DR.TMLE.mw<-geex_TMLE(simdat1, X ~ Z1_sq, Y_scaled ~ X*Z1*Z2, coef(target.model.Y0.mw), coef(target.model.Y1.mw), CM0.TMLE.mw, CM1.TMLE.mw, ATE.TMLE.mw)
  IFseATE.TMLE.mw<-IF_Var_mw(simdat1,'TMLE',ATE.TMLE.mw)
  bootseATE.TMLE.mw<-bootstrap_TMLE(num.boots, simdat1.obs, X ~ Z1_sq, Y_scaled ~ X*Z1*Z2)
  
  #compile results all 3 variance estimators
  DR.TMLE.mw2<-as.data.frame(cbind(DR.TMLE.mw,IFseATE.TMLE.mw,bootseATE.TMLE.mw))
  colnames(DR.TMLE.mw2)<-c('ATE', 'ESseATE','IFseATE','bootseATE')
  DR.TMLE.mw2$type<-'TMLE,MW'

  
  ########################################################
  ###########Estimation-Both models Incorrect  ###########
  ########################################################
  
  ############################
  #classic (plug-in) AIPW#
  ############################
  
  #estimate causal means
  CM1.PI.mb<-mean(((simdat1$X*simdat1$Y-(simdat1$X-simdat1$PS.mis)*simdat1$Y1.unwt.mis))/simdat1$PS.mis)
  CM0.PI.mb<-mean((((1-simdat1$X)*simdat1$Y+(simdat1$X-simdat1$PS.mis)*simdat1$Y0.unwt.mis))/(1-simdat1$PS.mis))
  ATE.PI.mb<-CM1.PI.mb - CM0.PI.mb
  
  #get 3 standard errors
  DR.PI.mb<-geex_PI(simdat1, X ~ Z1_sq, Y ~ X+Z1_sq, CM0.PI.mb, CM1.PI.mb, ATE.PI.mb)
  IFseATE.PI.mb<-IF_Var_mb(simdat1,'PI',ATE.PI.mb)
  bootseATE.PI.mb<-bootstrap_PIAIPW(num.boots, simdat1.obs, X ~ Z1_sq, Y ~ X+Z1_sq)
  
  #compile results all 3 variance estimators
  DR.PI.mb2<-as.data.frame(cbind(DR.PI.mb,IFseATE.PI.mb,bootseATE.PI.mb))
  colnames(DR.PI.mb2)<-c('ATE', 'ESseATE','IFseATE','bootseATE')
  DR.PI.mb2$type<-'PI,MB'
  
  
  
  ##########################
  #weighted regression AIPW#
  ##########################
  
  #fit outcome model (weighted)
  out.model.wt.mb <- glm(Y ~ X+Z1_sq, weights=IPTW.mis, data = simdat1)
  simdat1$Y0.wt.mb<-(predict(out.model.wt.mb, alluntrt, type = "response"))
  simdat1$Y1.wt.mb<-(predict(out.model.wt.mb, alltrt, type = "response"))
  
  #estimate causal means
  CM0.WTD.mb<-mean(simdat1$Y0.wt.mb)
  CM1.WTD.mb<-mean(simdat1$Y1.wt.mb)
  ATE.WTD.mb<-CM1.WTD.mb - CM0.WTD.mb
  
  #get 3 standard errors
  DR.WTD.mb<-geex_WTD_mis(simdat1, X ~ Z1_sq, Y ~ X+Z1_sq, coef(out.model.wt.mb), CM0.WTD.mb, CM1.WTD.mb, ATE.WTD.mb)
  IFseATE.WTD.mb<-IF_Var_mb(simdat1,'WTD',ATE.WTD.mb)
  bootseATE.WTD.mb<-bootstrap_WRAIPW(num.boots, simdat1.obs, X ~ Z1_sq, Y ~ X+Z1_sq)
  
  #compile results all 3 variance estimators
  DR.WTD.mb2<-as.data.frame(cbind(DR.WTD.mb,IFseATE.WTD.mb,bootseATE.WTD.mb))
  colnames(DR.WTD.mb2)<-c('ATE', 'ESseATE','IFseATE','bootseATE')
  DR.WTD.mb2$type<-'WTD,MB'
  
  
  ######
  #TMLE#
  ######
  
  #fit targeting models (both models incorrect)
  target.model.Y0.mb <- glm(Y_scaled ~ 1, offset=logit(Y0.sc.unwt.mis), weights = (1-X)/(1-PS.mis), 
                            family = binomial(), data = simdat1)
  target.model.Y1.mb <- glm(Y_scaled ~ 1, offset=logit(Y1.sc.unwt.mis), weights = X/PS.mis, 
                            family = binomial(), data = simdat1)
  
  #estimate causal means and ATE
  simdat1$Y0.TMLE.mb<-inv.logit(logit(simdat1$Y0.sc.unwt.mis)+coef(target.model.Y0.mb)[1])*(b-a)+a
  simdat1$Y1.TMLE.mb<-inv.logit(logit(simdat1$Y1.sc.unwt.mis)+coef(target.model.Y1.mb)[1])*(b-a)+a
  CM0.TMLE.mb<-mean(simdat1$Y0.TMLE.mb)
  CM1.TMLE.mb<-mean(simdat1$Y1.TMLE.mb)
  ATE.TMLE.mb<-CM1.TMLE.mb-CM0.TMLE.mb
  
  #get 3 standard errors 
  DR.TMLE.mb<-geex_TMLE(simdat1, X ~ Z1_sq, Y_scaled ~ X+Z1_sq, coef(target.model.Y0.mb), coef(target.model.Y1.mb), CM0.TMLE.mb, CM1.TMLE.mb, ATE.TMLE.mb)
  IFseATE.TMLE.mb<-IF_Var_mb(simdat1,'TMLE',ATE.TMLE.mb)
  bootseATE.TMLE.mb<-bootstrap_TMLE(num.boots, simdat1.obs, X ~ Z1_sq, Y_scaled ~ X+Z1_sq)
  
  DR.TMLE.mb2<-as.data.frame(cbind(DR.TMLE.mb,IFseATE.TMLE.mb,bootseATE.TMLE.mb))
  colnames(DR.TMLE.mb2)<-c('ATE', 'ESseATE','IFseATE','bootseATE')
  DR.TMLE.mb2$type<-'TMLE,MB'


  
  
  
  ######### combine all results and output #############
  DR.all<-rbind(DR.PI2,DR.WTD2,DR.TMLE2,DR.PI.mis2,DR.WTD.mis2,DR.TMLE.mis2,DR.PI.mw2,DR.WTD.mw2,DR.TMLE.mw2,DR.PI.mb2,DR.WTD.mb2,DR.TMLE.mb2)
  DR.all2<-cbind(DR.all,rep(sim,nrow(DR.all)))
  colnames(DR.all2)<-c('ATE', 'ESseATE','IFseATE','bootseATE','type','sim')
  
  output_filename <- paste(dir_path,"/results_", sim, ".csv", sep="")
  write.csv(DR.all2, output_filename)
  