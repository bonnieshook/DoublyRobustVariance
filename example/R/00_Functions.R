#Program: 00_Functions.R
#Purpose: This program contains functions for DR estimators 

library(numDeriv)
library(boot)
library(geex)
library(resample)

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
#### M-Estimators ####
##########################################################################################################

#### Classic (Plug-in) AIPW Estimation ####

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


#### Weighted Regression AIPW Estimation ####

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


#### TMLE ####

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


##########################################################################################################
#### Bootstrap Estimators ####
##########################################################################################################

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
      
      #fit outcome model (weighted), get pseudo POs under treatement and no treatment
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
















