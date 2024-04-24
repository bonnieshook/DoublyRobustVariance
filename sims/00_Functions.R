#Program: 00_Functions.R
#Purpose: This program contains functions for DR estimators 

library(numDeriv)
library(boot)
library(geex)

##########################################################################################################
#### IF Variance Estimator ####
##########################################################################################################

IF_Var <- function(data,method,est.DR){ 
  
data$method<-method
X<-data$X
Y<-data$Y
ehat<-data$PS
Yhat0<-ifelse(data$method=="PI",data$Y0.unwt,ifelse(data$method=="WTD",data$Y0.wt,ifelse(data$method=="TMLE",data$Y0.TMLE,rep(NA,nrow(data)))))
Yhat1<-ifelse(data$method=="PI",data$Y1.unwt,ifelse(data$method=="WTD",data$Y1.wt,ifelse(data$method=="TMLE",data$Y1.TMLE,rep(NA,nrow(data)))))
IF = (X*Y/ehat - (1-X)*Y/(1-ehat)) - ((X-ehat)/ehat/(1-ehat))*(((1-ehat)*Yhat1)+(ehat*Yhat0))-est.DR
sdATE=sqrt(sum(IF^2)/(length(IF)^2))

return(sdATE)
}

#misspecified outcome
IF_Var_mis <- function(data,method,est.DR){ 
  
  data$method<-method
  X<-data$X
  Y<-data$Y
  ehat<-data$PS
  Yhat0<-ifelse(data$method=="PI",data$Y0.unwt.mis,ifelse(data$method=="WTD",data$Y0.wt.mis,ifelse(data$method=="TMLE",data$Y0.TMLE.mis,rep(NA,nrow(data)))))
  Yhat1<-ifelse(data$method=="PI",data$Y1.unwt.mis,ifelse(data$method=="WTD",data$Y1.wt.mis,ifelse(data$method=="TMLE",data$Y1.TMLE.mis,rep(NA,nrow(data)))))
  IF = (X*Y/ehat - (1-X)*Y/(1-ehat)) - ((X-ehat)/ehat/(1-ehat))*(((1-ehat)*Yhat1)+(ehat*Yhat0))-est.DR
  sdATE=sqrt(sum(IF^2)/(length(IF)^2))
  
  return(sdATE)
}

#misspecified weight
IF_Var_mw <- function(data,method,est.DR){ 
  
  data$method<-method
  X<-data$X
  Y<-data$Y
  ehat<-data$PS.mis
  Yhat0<-ifelse(data$method=="PI",data$Y0.unwt,ifelse(data$method=="WTD",data$Y0.wt.mw,ifelse(data$method=="TMLE",data$Y0.TMLE.mw,rep(NA,nrow(data)))))
  Yhat1<-ifelse(data$method=="PI",data$Y1.unwt,ifelse(data$method=="WTD",data$Y1.wt.mw,ifelse(data$method=="TMLE",data$Y1.TMLE.mw,rep(NA,nrow(data)))))
  IF = (X*Y/ehat - (1-X)*Y/(1-ehat)) - ((X-ehat)/ehat/(1-ehat))*(((1-ehat)*Yhat1)+(ehat*Yhat0))-est.DR
  sdATE=sqrt(sum(IF^2)/(length(IF)^2))
  
  return(sdATE)
}

#both models misspecified 
IF_Var_mb <- function(data,method,est.DR){ 
  
  data$method<-method
  X<-data$X
  Y<-data$Y
  ehat<-data$PS.mis
  Yhat0<-ifelse(data$method=="PI",data$Y0.unwt.mis,ifelse(data$method=="WTD",data$Y0.wt.mb,ifelse(data$method=="TMLE",data$Y0.TMLE.mb,rep(NA,nrow(data)))))
  Yhat1<-ifelse(data$method=="PI",data$Y1.unwt.mis,ifelse(data$method=="WTD",data$Y1.wt.mb,ifelse(data$method=="TMLE",data$Y1.TMLE.mb,rep(NA,nrow(data)))))
  IF = (X*Y/ehat - (1-X)*Y/(1-ehat)) - ((X-ehat)/ehat/(1-ehat))*(((1-ehat)*Yhat1)+(ehat*Yhat0))-est.DR
  sdATE=sqrt(sum(IF^2)/(length(IF)^2))
  
  return(sdATE)
}


##########################################################################################################
#### Classic (Plug-in) AIPW Estimation ####
##########################################################################################################

estfun_PI <- function(data,models){   
  
  X<-data$X
  Y<-data$Y
  
  Xe <- grab_design_matrix(data = data, rhs_formula = grab_fixed_formula(models$e))
  Xm <- grab_design_matrix(data = data, rhs_formula = grab_fixed_formula(models$m))
  
  data0 <- data
  data0$X <- rep(0,nrow(data))
  Xm0 <- grab_design_matrix(data = data0,rhs_formula = grab_fixed_formula(models$m))
  
  data1 <- data
  data1$X <- rep(1,nrow(data))
  Xm1 <- grab_design_matrix(data = data1,rhs_formula = grab_fixed_formula(models$m))
  
  e_pos <- 1:ncol(Xe)
  m_pos <- (max(e_pos) + 1):(max(e_pos) + ncol(Xm))
  
  e_scores <- grab_psiFUN(models$e, data)
  m_scores <- grab_psiFUN(models$m, data)
  
  function(theta){
    p<-length(theta)
    e <- plogis(Xe %*% theta[e_pos]) 
    m0 <- Xm0 %*% theta[m_pos]
    m1 <- Xm1 %*% theta[m_pos]
    
    CM1<-((X*Y - (X-e)*m1)/e)
    CM0<-(((1-X)*Y + (X-e)*m0)/(1-e)) 
    c(e_scores(theta[e_pos]),
      m_scores(theta[m_pos]),
      CM1-theta[p-2], 
      CM0-theta[p-1],
      theta[p-2]- theta[p-1]-theta[p])
  }
}

geex_PI <- function(data, propensity_formula, outcome_formula, CM0.IV, CM1.IV, ATE.IV){
  
  e_model  <- glm(propensity_formula, data = data, family =binomial)
  m_model  <- glm(outcome_formula, data = data)
  models <- list(e = e_model, m=m_model)
  
  geex_resultsdrPI <-m_estimate(
    estFUN = estfun_PI, 
    data   = data, 
    roots = c(coef(e_model),coef(m_model), CM1.IV, CM0.IV, ATE.IV),
    compute_roots = FALSE,
    outer_args = list(models = models))

    DR.PI<-geex_resultsdrPI@estimates[length(geex_resultsdrPI@estimates)]
    seDR.PI <- sqrt(geex_resultsdrPI@vcov[length(geex_resultsdrPI@estimates),length(geex_resultsdrPI@estimates)])
    DR.PI.all<-cbind(DR.PI,seDR.PI)
  
  return(DR.PI.all)
}


##########################################################################################################
#### Weighted Regression AIPW Estimation ####
##########################################################################################################

ScoreWTD<- function(X,Z1,Z2,Y,mu,W){
  c(W*(Y-mu),
    W*Z1*(Y-mu),
    W*Z2*(Y-mu),
    W*X*(Y-mu),
    W*Z1*Z2*(Y-mu),
    W*Z1*X*(Y-mu), 
    W*Z2*X*(Y-mu),
    W*Z1*Z2*X*(Y-mu))}

estfun_WTD <- function(data,models){   
  
  X<-data$X
  Z1<-data$Z1
  Z2<-data$Z2
  Y<-data$Y
  
  Xe <- grab_design_matrix(data = data, rhs_formula = grab_fixed_formula(models$e))
  Xm <- grab_design_matrix(data = data, rhs_formula = grab_fixed_formula(models$m))
  data0 <- data
  data0$X <- rep(0,nrow(data))
  Xm0 <- grab_design_matrix(data = data0,rhs_formula = grab_fixed_formula(models$m))
  
  data1 <- data
  data1$X <- rep(1,nrow(data))
  Xm1 <- grab_design_matrix(data = data1,rhs_formula = grab_fixed_formula(models$m))
  
  e_pos <- 1:ncol(Xe)
  m_pos <- (max(e_pos) + 1):(max(e_pos) + ncol(Xm))
  
  e_scores <- grab_psiFUN(models$e, data)
  
  function(theta){
    p<-length(theta)
    e <- plogis(Xe %*% theta[e_pos]) 
    mu <- Xm %*% theta[m_pos]
    m0 <- Xm0 %*% theta[m_pos]
    m1 <- Xm1 %*% theta[m_pos]
    W<-X*e^(-1)+(1-X)*(1-e)^(-1)    
    
    c(e_scores(theta[e_pos]),
      ScoreWTD(X,Z1,Z2,Y,mu,W),
      m1-theta[p-2], 
      m0-theta[p-1],
      theta[p-2]- theta[p-1]-theta[p])
  }
}

geex_WTD <- function(data, propensity_formula, outcome_formula,coef.WTD, CM0.IV, CM1.IV, ATE.IV){
  
  e_model  <- glm(propensity_formula, data = data, family =binomial)
  m_model  <- glm(outcome_formula, data = data)
  models <- list(e = e_model, m=m_model)
  
  geex_resultsdrWTD <-m_estimate(
    estFUN = estfun_WTD, 
    data   = data, 
    roots = c(coef(e_model),coef.WTD, CM1.IV, CM0.IV, ATE.IV),
    compute_roots = FALSE,
    outer_args = list(models = models))
  
  DR.WTD<-geex_resultsdrWTD@estimates[length(geex_resultsdrWTD@estimates)]
  seDR.WTD <- sqrt(geex_resultsdrWTD@vcov[length(geex_resultsdrWTD@estimates),length(geex_resultsdrWTD@estimates)])
  DR.WTD.all<-cbind(DR.WTD,seDR.WTD)
  
  return(DR.WTD.all)
}



##########################################################################################################
#### Weighted Regression AIPW Estimation (Misspec) ####
##########################################################################################################

ScoreWTD_mis<- function(X,Z1_sq,Y,mu,W){
  c(W*(Y-mu),
    W*X*(Y-mu),
    W*Z1_sq*(Y-mu))
}

estfun_WTD_mis <- function(data,models){   
  
  X<-data$X
  Y<-data$Y
  Z1_sq<-data$Z1_sq
  
  Xe <- grab_design_matrix(data = data, rhs_formula = grab_fixed_formula(models$e))
  Xm <- grab_design_matrix(data = data, rhs_formula = grab_fixed_formula(models$m))
  data0 <- data
  data0$X <- rep(0,nrow(data))
  Xm0 <- grab_design_matrix(data = data0,rhs_formula = grab_fixed_formula(models$m))
  
  data1 <- data
  data1$X <- rep(1,nrow(data))
  Xm1 <- grab_design_matrix(data = data1,rhs_formula = grab_fixed_formula(models$m))
  
  e_pos <- 1:ncol(Xe)
  m_pos <- (max(e_pos) + 1):(max(e_pos) + ncol(Xm))
  
  e_scores <- grab_psiFUN(models$e, data)
  
  function(theta){
    p<-length(theta)
    e <- plogis(Xe %*% theta[e_pos]) 
    mu <- Xm %*% theta[m_pos]
    m0 <- Xm0 %*% theta[m_pos]
    m1 <- Xm1 %*% theta[m_pos]
    W<-X*e^(-1)+(1-X)*(1-e)^(-1)    
    
    c(e_scores(theta[e_pos]),
      ScoreWTD_mis(X,Z1_sq,Y,mu,W),
      m1-theta[p-2], 
      m0-theta[p-1],
      theta[p-2]- theta[p-1]-theta[p])
  }
}

geex_WTD_mis <- function(data, propensity_formula, outcome_formula,coef.WTD, CM0.IV, CM1.IV, ATE.IV){
  
  e_model  <- glm(propensity_formula, data = data, family =binomial)
  m_model  <- glm(outcome_formula, data = data)
  models <- list(e = e_model, m=m_model)
  
  geex_resultsdrWTD <-m_estimate(
    estFUN = estfun_WTD_mis, 
    data   = data, 
    roots = c(coef(e_model),coef.WTD, CM1.IV, CM0.IV, ATE.IV),
    compute_roots = FALSE,
    outer_args = list(models = models))
  
  DR.WTD<-geex_resultsdrWTD@estimates[length(geex_resultsdrWTD@estimates)]
  seDR.WTD <- sqrt(geex_resultsdrWTD@vcov[length(geex_resultsdrWTD@estimates),length(geex_resultsdrWTD@estimates)])
  DR.WTD.all<-cbind(DR.WTD,seDR.WTD)
  
  return(DR.WTD.all)
}



##########################################################################################################
#### TMLE ####
##########################################################################################################

ScoreWTD_tar<- function(Y_sc,mu,W){
  c(W*(Y_sc-mu))
}

estfun_TMLE <- function(data,models){   
  
  X<-data$X
  Y<-data$Y
  Y_sc<-data$Y_scaled
  a<-data$a
  b<-data$b
  
  Xe <- grab_design_matrix(data = data, rhs_formula = grab_fixed_formula(models$e))
  Xm <- grab_design_matrix(data = data, rhs_formula = grab_fixed_formula(models$m))
  data0 <- data
  data0$X <- rep(0,nrow(data))
  Xm0 <- grab_design_matrix(data = data0,rhs_formula = grab_fixed_formula(models$m))
  
  data1 <- data
  data1$X <- rep(1,nrow(data))
  Xm1 <- grab_design_matrix(data = data1,rhs_formula = grab_fixed_formula(models$m))
  
  e_pos <- 1:ncol(Xe)
  m_pos <- (max(e_pos) + 1):(max(e_pos) + ncol(Xm))
  
  e_scores <- grab_psiFUN(models$e, data)
  m_scores <- grab_psiFUN(models$m, data)
  
  function(theta){
    p<-length(theta)
    e <- plogis(Xe %*% theta[e_pos]) 
    m0 <- Xm0 %*% theta[m_pos]
    m1 <- Xm1 %*% theta[m_pos]
    
    mu0 <- inv.logit(theta[max(m_pos)+1]+logit(m0))
    mu1 <- inv.logit(theta[max(m_pos)+2]+logit(m1))
    w0<-(1-X)/(1-e)
    w1<-X/e
    
    CM0<-inv.logit(logit(m0)+theta[max(m_pos)+1])*(b-a)+a
    CM1<-inv.logit(logit(m1)+theta[max(m_pos)+2])*(b-a)+a

    c(e_scores(theta[e_pos]),
      m_scores(theta[m_pos]),
      ScoreWTD_tar(Y_sc,mu0,w0),
      ScoreWTD_tar(Y_sc,mu1,w1),      
      CM1-theta[p-2], 
      CM0-theta[p-1],
      theta[p-2]- theta[p-1]-theta[p])
  }
}

geex_TMLE <- function(data, propensity_formula, outcome_formula, tar0.IV, tar1.IV, CM0.IV, CM1.IV, ATE.IV){
  
  e_model  <- glm(propensity_formula, data = data, family =binomial)
  m_model  <- glm(outcome_formula, data = data)
  models <- list(e = e_model, m=m_model)
  
  geex_resultsdrTMLE <-m_estimate(
    estFUN = estfun_TMLE, 
    data   = data, 
    roots = c(coef(e_model),coef(m_model), tar0.IV, tar1.IV, CM1.IV, CM0.IV, ATE.IV),
    compute_roots = FALSE,
    outer_args = list(models = models))
  
  DR.TMLE<-geex_resultsdrTMLE@estimates[length(geex_resultsdrTMLE@estimates)]
  seDR.TMLE <- sqrt(geex_resultsdrTMLE@vcov[length(geex_resultsdrTMLE@estimates),length(geex_resultsdrTMLE@estimates)])
  DR.TMLE.all<-cbind(DR.TMLE,seDR.TMLE)
  
  return(DR.TMLE.all)
}











