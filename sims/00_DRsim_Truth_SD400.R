#Program: Determine the ACE empirically based on a sample of 50M obs from the DGM
#Developed by: BES
#Created: 04.19.24

library(boot)
library(stats)

num.part <- 50000000
set.seed(112023)

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
  simdat1$EY0 = 1000+11.5*simdat1$Z1 + 100*simdat1$Z2  - 15*simdat1$Z1*simdat1$Z2   
  simdat1$EY1 = 1025+11.5*simdat1$Z1 + 100*simdat1$Z2 - 5.5*simdat1$Z1 - 30*simdat1$Z2 + 5*simdat1$Z1*simdat1$Z2  
  simdat1$Y0<- rnorm(num.part, mean =simdat1$EY0, sd=400)
  simdat1$Y1<- rnorm(num.part, mean =simdat1$EY1, sd=400)
  
  #define the ACE and output
  ACE<-mean(simdat1$Y1)-mean(simdat1$Y0)
  write.csv(ACE, "ACE_Truth_SD400.csv")
  