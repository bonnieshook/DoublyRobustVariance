#Program: 00_Functions.R
#Purpose: This program contains functions for DR estimators 

library(numDeriv)
library(boot)
library(geex)
library(resample)

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

### linear outcome model
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

### logit outcome model
estfun_PI_bin <- function(data,models){   
  
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
    m0 <- plogis(Xm0 %*% theta[m_pos])
    m1 <- plogis(Xm1 %*% theta[m_pos])
    
    CM1<-((X*Y - (X-e)*m1)/e)
    CM0<-(((1-X)*Y + (X-e)*m0)/(1-e)) 
    c(e_scores(theta[e_pos]),
      m_scores(theta[m_pos]),
      CM1-theta[p-2], 
      CM0-theta[p-1],
      theta[p-2]- theta[p-1]-theta[p])
  }
}

geex_PI_bin <- function(data, propensity_formula, outcome_formula, CM0.IV, CM1.IV, ATE.IV){
  
  e_model  <- glm(propensity_formula, data = data, family =binomial)
  m_model  <- glm(outcome_formula, data = data, family =binomial)
  models <- list(e = e_model, m=m_model)
  
  geex_resultsdrPI <-m_estimate(
    estFUN = estfun_PI_bin, 
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

#### Linear outcome model
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


#### Logistic outcome model
estfun_WTD_bin <- function(data,models){   
  
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
    mu <- plogis(Xm %*% theta[m_pos])
    m0 <- plogis(Xm0 %*% theta[m_pos])
    m1 <- plogis(Xm1 %*% theta[m_pos])
    W<-X*e^(-1)+(1-X)*(1-e)^(-1)    
    
    c(e_scores(theta[e_pos]),
      ScoreWTD(X,Z1,Z2,Y,mu,W),
      m1-theta[p-2], 
      m0-theta[p-1],
      theta[p-2]- theta[p-1]-theta[p])
  }
}

geex_WTD_bin <- function(data, propensity_formula, outcome_formula,coef.WTD, CM0.IV, CM1.IV, ATE.IV){
  
  e_model  <- glm(propensity_formula, data = data, family =binomial)
  m_model  <- glm(outcome_formula, data = data, family =binomial)
  models <- list(e = e_model, m=m_model)
  
  geex_resultsdrWTD <-m_estimate(
    estFUN = estfun_WTD_bin, 
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

### Linear outcome model
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


### Logit outcome model
estfun_WTD_mis_bin <- function(data,models){   
  
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
    mu <- plogis(Xm %*% theta[m_pos])
    m0 <- plogis(Xm0 %*% theta[m_pos])
    m1 <- plogis(Xm1 %*% theta[m_pos])
    W<-X*e^(-1)+(1-X)*(1-e)^(-1)    
    
    c(e_scores(theta[e_pos]),
      ScoreWTD_mis(X,Z1_sq,Y,mu,W),
      m1-theta[p-2], 
      m0-theta[p-1],
      theta[p-2]- theta[p-1]-theta[p])
  }
}

geex_WTD_mis_bin <- function(data, propensity_formula, outcome_formula,coef.WTD, CM0.IV, CM1.IV, ATE.IV){
  
  e_model  <- glm(propensity_formula, data = data, family =binomial)
  m_model  <- glm(outcome_formula, data = data, family =binomial)
  models <- list(e = e_model, m=m_model)
  
  geex_resultsdrWTD <-m_estimate(
    estFUN = estfun_WTD_mis_bin, 
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

### Linear outcome model
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

### Logit outcome model
estfun_TMLE_bin <- function(data,models){   
  
  X<-data$X
  Y<-data$Y
  Y_sc<-data$Y_scaled
  
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
    m0 <- plogis(Xm0 %*% theta[m_pos])
    m1 <- plogis(Xm1 %*% theta[m_pos])
    
    mu0 <- inv.logit(theta[max(m_pos)+1]+logit(m0))
    mu1 <- inv.logit(theta[max(m_pos)+2]+logit(m1))
    w0<-(1-X)/(1-e)
    w1<-X/e
    
    CM0<-inv.logit(logit(m0)+theta[max(m_pos)+1])
    CM1<-inv.logit(logit(m1)+theta[max(m_pos)+2])
    
    c(e_scores(theta[e_pos]),
      m_scores(theta[m_pos]),
      ScoreWTD_tar(Y_sc,mu0,w0),
      ScoreWTD_tar(Y_sc,mu1,w1),      
      CM1-theta[p-2], 
      CM0-theta[p-1],
      theta[p-2]- theta[p-1]-theta[p])
  }
}

geex_TMLE_bin <- function(data, propensity_formula, outcome_formula, tar0.IV, tar1.IV, CM0.IV, CM1.IV, ATE.IV){
  
  e_model  <- glm(propensity_formula, data = data, family =binomial)
  m_model  <- glm(outcome_formula, data = data, family =binomial)
  models <- list(e = e_model, m=m_model)
  
  geex_resultsdrTMLE <-m_estimate(
    estFUN = estfun_TMLE_bin, 
    data   = data, 
    roots = c(coef(e_model),coef(m_model), tar0.IV, tar1.IV, CM1.IV, CM0.IV, ATE.IV),
    compute_roots = FALSE,
    outer_args = list(models = models))
  
  DR.TMLE<-geex_resultsdrTMLE@estimates[length(geex_resultsdrTMLE@estimates)]
  seDR.TMLE <- sqrt(geex_resultsdrTMLE@vcov[length(geex_resultsdrTMLE@estimates),length(geex_resultsdrTMLE@estimates)])
  DR.TMLE.all<-cbind(DR.TMLE,seDR.TMLE)
  
  return(DR.TMLE.all)
}



##########################################################################################################
#### Bootstrap Estimators ####
##########################################################################################################

### functions to get bootstrap SEs for each estimator of interest

############ Linear outcome model

### Plug-in AIPW
bootstrap_PIAIPW <- function(B, bootdata, propensity_formula, outcome_formula){
  if(B == 0) return(0)
  if(B>0){
    boot.est <- matrix(NaN, nrow = B, ncol = 1)
    datbi<-samp.bootstrap(nrow(bootdata), B)
    for(i in 1:B){
      dati <- bootdata[datbi[,i],]
      skip_to_next=FALSE
      
      ### Fit weight model, get estimated propensity scores, and assign IPTWs
      tryCatch(prop.model  <- glm(propensity_formula, data = dati, family = binomial), error = function(e) { skip_to_next <<- TRUE})
      if(skip_to_next==FALSE) {
        dati$PS<-predict(prop.model,dati, type="response")
        dati$IPTW<-ifelse(dati$X==1,dati$PS^(-1),(1-dati$PS)^(-1))}
      
      #create copies of datasets for estimated potential outcomes
      alltrt.boot<-alluntrt.boot<-dati
      alltrt.boot$X<-1
      alluntrt.boot$X<-0
      
      #fit outcome model (unweighted)
      tryCatch(out.model.unwt  <- glm(outcome_formula, data = dati), error = function(e) { skip_to_next <<- TRUE})
      if(skip_to_next==FALSE) {
        dati$Y0.unwt<-(predict(out.model.unwt, alluntrt.boot, type = "response"))
        dati$Y1.unwt<-(predict(out.model.unwt, alltrt.boot, type = "response"))
        
        #estimate causal means
        CM1.PI<-mean(((dati$X*dati$Y-(dati$X-dati$PS)*dati$Y1.unwt))/dati$PS)
        CM0.PI<-mean((((1-dati$X)*dati$Y+(dati$X-dati$PS)*dati$Y0.unwt))/(1-dati$PS))
        boot.est [i]<-CM1.PI - CM0.PI}
      if(skip_to_next==TRUE) {
        boot.est [i]<-NA}
    }
    boot.se<-sd(boot.est,na.rm=TRUE)
    return(boot.se)
  }}


### Weighted Regression AIPW
bootstrap_WRAIPW <- function(B, bootdata, propensity_formula, outcome_formula){
  if(B == 0) return(0)
  if(B>0){
    boot.est <- matrix(NaN, nrow = B, ncol = 1)
    datbi<-samp.bootstrap(nrow(bootdata), B)
    for(i in 1:B){
      dati <- bootdata[datbi[,i],]
      skip_to_next=FALSE
      
      # Fit weight model, get estimated propensity scores, and assign IPTWs
      tryCatch(prop.model  <- glm(propensity_formula, data = dati, family =binomial), error = function(e) { skip_to_next <<- TRUE})
      if(skip_to_next==FALSE) {
        dati$PS<-predict(prop.model,dati, type="response")
        dati$IPTW<-ifelse(dati$X==1,dati$PS^(-1),(1-dati$PS)^(-1))}
      
      #create copies of datasets for estimated potential outcomes
      alltrt.boot<-alluntrt.boot<-dati
      alltrt.boot$X<-1
      alluntrt.boot$X<-0
      
      #fit outcome model (weighted), get pseudo POs under treatment and no treatment
      tryCatch(out.model.wt  <- glm(outcome_formula, data = dati,weights=IPTW), error = function(e) { skip_to_next <<- TRUE})
      if(skip_to_next==FALSE) {
        dati$Y0.wt<-(predict(out.model.wt, alluntrt.boot, type = "response"))
        dati$Y1.wt<-(predict(out.model.wt, alltrt.boot, type = "response"))
        
        #estimate causal means and ATE
        CM0.WTD<-mean(dati$Y0.wt)
        CM1.WTD<-mean(dati$Y1.wt)
        boot.est [i]<-CM1.WTD - CM0.WTD}
      
      if(skip_to_next==TRUE) {
        boot.est [i]<-NA}
    }
    boot.se<-sd(boot.est,na.rm=TRUE)
    return(boot.se)
  }}


### TMLE
bootstrap_TMLE <- function(B, bootdata, propensity_formula, outcome_formula){
  if(B == 0) return(0)
  if(B>0){
    boot.est <- matrix(NaN, nrow = B, ncol = 1)
    datbi<-samp.bootstrap(nrow(bootdata), B)
    for(i in 1:B){
      dati <- bootdata[datbi[,i],]
      skip_to_next=FALSE
      
      # Fit weight model and get estimated propensity scores
      tryCatch(prop.model  <- glm(propensity_formula, data = dati, family =binomial), error = function(e) { skip_to_next <<- TRUE})
      if(skip_to_next==FALSE) {
        dati$PS<-predict(prop.model,dati, type="response")}
      
      #create copies of datasets for estimated potential outcomes
      alltrt.boot<-alluntrt.boot<-dati
      alltrt.boot$X<-1
      alluntrt.boot$X<-0
      
      #scale the outcome for TMLE
      a.boot<-min(dati$Y)
      b.boot<-max(dati$Y)
      dati$a.boot<-a.boot
      dati$b.boot<-b.boot
      dati$Y_scaled<-(dati$Y-a.boot)/(b.boot-a.boot)
      
      #fit scaled outcome model (unweighted)
      tryCatch(out.model.sc.unwt <- glm(outcome_formula, data = dati), error = function(e) { skip_to_next <<- TRUE})
      if(skip_to_next==FALSE) {
        dati$Y0.sc.unwt<-(predict(out.model.sc.unwt, alluntrt.boot, type = "response"))
        dati$Y1.sc.unwt<-(predict(out.model.sc.unwt, alltrt.boot, type = "response"))
        
        #truncate predictions to within the range of the data
        dati$Y0.sc.unwt<-ifelse(dati$Y0.sc.unwt<min(dati$Y_scaled),min(dati$Y_scaled),dati$Y0.sc.unwt)
        dati$Y0.sc.unwt<-ifelse(dati$Y0.sc.unwt>max(dati$Y_scaled),max(dati$Y_scaled),dati$Y0.sc.unwt)
        
        dati$Y1.sc.unwt<-ifelse(dati$Y1.sc.unwt<min(dati$Y_scaled),min(dati$Y_scaled),dati$Y1.sc.unwt)
        dati$Y1.sc.unwt<-ifelse(dati$Y1.sc.unwt>max(dati$Y_scaled),max(dati$Y_scaled),dati$Y1.sc.unwt)}
      
      #fit TMLE targeting model
      tryCatch(target.model.Y0.boot <- glm(Y_scaled ~ 1, offset=logit(Y0.sc.unwt), weights = (1-X)/(1-PS), 
                                           family = binomial(), data = dati), error = function(e) { skip_to_next <<- TRUE})
      tryCatch(target.model.Y1.boot <- glm(Y_scaled ~ 1, offset=logit(Y1.sc.unwt), weights = X/PS, 
                                           family = binomial(), data = dati), error = function(e) { skip_to_next <<- TRUE})
      if(skip_to_next==FALSE) {
        dati$Y0.TMLE<-inv.logit(logit(dati$Y0.sc.unwt)+coef(target.model.Y0.boot)[1])*(b.boot-a.boot)+a.boot
        dati$Y1.TMLE<-inv.logit(logit(dati$Y1.sc.unwt)+coef(target.model.Y1.boot)[1])*(b.boot-a.boot)+a.boot
        CM0.TMLE<-mean(dati$Y0.TMLE)
        CM1.TMLE<-mean(dati$Y1.TMLE)
        boot.est [i]<-CM1.TMLE-CM0.TMLE}
      
      if(skip_to_next==TRUE) {
        boot.est [i]<-NA}
      
    }
    boot.se<-sd(boot.est,na.rm=TRUE)
    return(boot.se)
  }}


############ Logistic outcome model

### Plug-in AIPW
bootstrap_PIAIPW_bin <- function(B, bootdata, propensity_formula, outcome_formula){
  if(B == 0) return(0)
  if(B>0){
    boot.est <- matrix(NaN, nrow = B, ncol = 1)
    datbi<-samp.bootstrap(nrow(bootdata), B)
    for(i in 1:B){
      dati <- bootdata[datbi[,i],]
      skip_to_next=FALSE
      
      ### Fit weight model, get estimated propensity scores, and assign IPTWs
      tryCatch(prop.model  <- glm(propensity_formula, data = dati, family = binomial), error = function(e) { skip_to_next <<- TRUE})
      if(skip_to_next==FALSE) {
        dati$PS<-predict(prop.model,dati, type="response")
        dati$IPTW<-ifelse(dati$X==1,dati$PS^(-1),(1-dati$PS)^(-1))}
      
      #create copies of datasets for estimated potential outcomes
      alltrt.boot<-alluntrt.boot<-dati
      alltrt.boot$X<-1
      alluntrt.boot$X<-0
      
      #fit outcome model (unweighted)
      tryCatch(out.model.unwt  <- glm(outcome_formula, data = dati, family = binomial), error = function(e) { skip_to_next <<- TRUE})
      if(skip_to_next==FALSE) {
        dati$Y0.unwt<-(predict(out.model.unwt, alluntrt.boot, type = "response"))
        dati$Y1.unwt<-(predict(out.model.unwt, alltrt.boot, type = "response"))
        
        #estimate causal means
        CM1.PI<-mean(((dati$X*dati$Y-(dati$X-dati$PS)*dati$Y1.unwt))/dati$PS)
        CM0.PI<-mean((((1-dati$X)*dati$Y+(dati$X-dati$PS)*dati$Y0.unwt))/(1-dati$PS))
        boot.est [i]<-CM1.PI - CM0.PI}
      if(skip_to_next==TRUE) {
        boot.est [i]<-NA}
    }
    boot.se<-sd(boot.est,na.rm=TRUE)
    return(boot.se)
  }}


### Weighted Regression AIPW
bootstrap_WRAIPW_bin <- function(B, bootdata, propensity_formula, outcome_formula){
  if(B == 0) return(0)
  if(B>0){
    boot.est <- matrix(NaN, nrow = B, ncol = 1)
    datbi<-samp.bootstrap(nrow(bootdata), B)
    for(i in 1:B){
      dati <- bootdata[datbi[,i],]
      skip_to_next=FALSE
      
      # Fit weight model, get estimated propensity scores, and assign IPTWs
      tryCatch(prop.model  <- glm(propensity_formula, data = dati, family =binomial), error = function(e) { skip_to_next <<- TRUE})
      if(skip_to_next==FALSE) {
        dati$PS<-predict(prop.model,dati, type="response")
        dati$IPTW<-ifelse(dati$X==1,dati$PS^(-1),(1-dati$PS)^(-1))}
      
      #create copies of datasets for estimated potential outcomes
      alltrt.boot<-alluntrt.boot<-dati
      alltrt.boot$X<-1
      alluntrt.boot$X<-0
      
      #fit outcome model (weighted), get pseudo POs under treatement and no treatment
      tryCatch(out.model.wt  <- glm(outcome_formula, data = dati,weights=IPTW, family=binomial), error = function(e) { skip_to_next <<- TRUE})
      if(skip_to_next==FALSE) {
        dati$Y0.wt<-(predict(out.model.wt, alluntrt.boot, type = "response"))
        dati$Y1.wt<-(predict(out.model.wt, alltrt.boot, type = "response"))
        
        #estimate causal means and ATE
        CM0.WTD<-mean(dati$Y0.wt)
        CM1.WTD<-mean(dati$Y1.wt)
        boot.est [i]<-CM1.WTD - CM0.WTD}
      
      if(skip_to_next==TRUE) {
        boot.est [i]<-NA}
    }
    boot.se<-sd(boot.est,na.rm=TRUE)
    return(boot.se)
  }}


### TMLE
bootstrap_TMLE_bin <- function(B, bootdata, propensity_formula, outcome_formula){
  if(B == 0) return(0)
  if(B>0){
    boot.est <- matrix(NaN, nrow = B, ncol = 1)
    datbi<-samp.bootstrap(nrow(bootdata), B)
    for(i in 1:B){
      dati <- bootdata[datbi[,i],]
      skip_to_next=FALSE
      
      # Fit weight model and get estimated propensity scores
      tryCatch(prop.model  <- glm(propensity_formula, data = dati, family =binomial), error = function(e) { skip_to_next <<- TRUE})
      if(skip_to_next==FALSE) {
        dati$PS<-predict(prop.model,dati, type="response")}
      
      #create copies of datasets for estimated potential outcomes
      alltrt.boot<-alluntrt.boot<-dati
      alltrt.boot$X<-1
      alluntrt.boot$X<-0
      
      #fit scaled outcome model (unweighted)
      tryCatch(out.model.sc.unwt <- glm(outcome_formula, data = dati, family =binomial), error = function(e) { skip_to_next <<- TRUE})
      if(skip_to_next==FALSE) {
        dati$Y0.sc.unwt<-(predict(out.model.sc.unwt, alluntrt.boot, type = "response"))
        dati$Y1.sc.unwt<-(predict(out.model.sc.unwt, alltrt.boot, type = "response"))}
      
      #fit TMLE targeting model
      tryCatch(target.model.Y0.boot <- glm(Y ~ 1, offset=logit(Y0.sc.unwt), weights = (1-X)/(1-PS), 
                                           family = binomial(), data = dati), error = function(e) { skip_to_next <<- TRUE})
      tryCatch(target.model.Y1.boot <- glm(Y ~ 1, offset=logit(Y1.sc.unwt), weights = X/PS, 
                                           family = binomial(), data = dati), error = function(e) { skip_to_next <<- TRUE})
      if(skip_to_next==FALSE) {
        dati$Y0.TMLE<-inv.logit(logit(dati$Y0.sc.unwt)+coef(target.model.Y0.boot)[1])
        dati$Y1.TMLE<-inv.logit(logit(dati$Y1.sc.unwt)+coef(target.model.Y1.boot)[1])
        CM0.TMLE<-mean(dati$Y0.TMLE)
        CM1.TMLE<-mean(dati$Y1.TMLE)
        boot.est [i]<-CM1.TMLE-CM0.TMLE}
      
      if(skip_to_next==TRUE) {
        boot.est [i]<-NA}
      
    }
    boot.se<-sd(boot.est,na.rm=TRUE)
    return(boot.se)
  }}









