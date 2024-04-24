#Program: 00_Functions.R
#Purpose: This program contains functions for DR estimators 

library(numDeriv)
library(boot)
library(geex)

##########################################################################################################
#### IF Variance Estimator ####
##########################################################################################################

#correct specification
IF_Var <- function(exposure,outcome,prop.score,Yhat0,Yhat1,est.DR){ 
  
X<-exposure
Y<-outcome
ehat<-prop.score
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

# Helper function
gdm <- \(data, model) {
  grab_design_matrix(data = data, rhs_formula = grab_fixed_formula(model))
}

# Estimating function
estfun_WTD <- \(data, models) {
  # Gather necessary data
  X <- data$X
  Xe <- gdm(data, models$e)
  Xm <- gdm(data, models$m)
  data0 <- data1 <- data
  data0$X <- 0
  data1$X <- 1
  Xm0 <- gdm(data0, models$m)
  Xm1 <- gdm(data1, models$m)
  e_scores <- grab_psiFUN(models$e, data)
  m_scores <- grab_psiFUN(models$m, data)
  
  # Compute column indices
  e_pos <- seq_len(ncol(Xe))
  m_pos <- seq_len(ncol(Xm)) + max(e_pos)
  
  
  \(theta){
    
    p <- length(theta)
    e <- plogis(Xe %*% theta[e_pos])
    m0 <- Xm0 %*% theta[m_pos]
    m1 <- Xm1 %*% theta[m_pos]
    W  <- 1/dbinom(X, size = 1, prob = e)
    
    c(e_scores(theta[e_pos]),
      W * m_scores(theta[m_pos]),
      m1 - theta[p - 2],
      m0 - theta[p - 1],
      theta[p - 2] - theta[p - 1] - theta[p]
    )
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











