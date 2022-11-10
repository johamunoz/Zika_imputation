
library(dplyr)
library("xlsx")
library(mitml)
library(mice)
library(micemd)
library(metafor)
library(metamisc)
library(EpiStats)
library("lme4")
library("logistf")

#######Functions#######
logit <- function(x) {
  log(x/(1-x))
}
inv.logit <- function(x) {
  1/(1+exp(-x))
}


#Functions

#Function to do one stage meta-analysis with log link and random intercept per study, per imputed dataset
f.1ma.r.int<-function(data, outcome_name) {
  for (i in 1:length(unique(data$.imp))) {
    d<-data[data$.imp==i,]
    d$outcome <- d[, outcome_name]
    fit1 <- glmer(outcome ~ zikv_preg + (1 | studyname_fac), 
                  data=d, family = binomial(link = "log"))
    #Store results
    if(i==1) {fit1.coef<-summary(fit1)$coefficients[2,(1:2)]}
    if(i>1) {fit1.coef<-rbind(fit1.coef,summary(fit1)$coefficients[2,(1:2)])}
  }
  colnames(fit1.coef)<-c("log.rr","log.se")
  rownames(fit1.coef)<-NULL
  return(fit1.coef)
}
#Function to pool the results of the one stage meta-analysis using Rubins rules
f.1ma.poolrubin <-function(data,fit1.coef) {
  #Pool with Rubins rules
  pool.rubin<- data.frame(matrix(NA, nrow = 1, ncol = 7))
  colnames(pool.rubin)<-c("log.rr.1ma","within","between","log.var",
                          "log.se.1ma","log.ci.lb.1ma","log.ci.ub.1ma")
  pool.rubin$log.rr.1ma <- mean(fit1.coef$log.rr)
  pool.rubin$within <- mean(fit1.coef$log.se^2)
  pool.rubin$between <- (1 + (1/m)) * var(fit1.coef$log.rr)
  pool.rubin$log.var <- pool.rubin$within + pool.rubin$between
  pool.rubin$log.se.1ma <- sqrt(pool.rubin$log.var)
  pool.rubin$log.ci.lb.1ma <- pool.rubin$log.rr.1ma + qnorm(0.05/2)     * pool.rubin$log.se.1ma
  pool.rubin$log.ci.ub.1ma <- pool.rubin$log.rr.1ma + qnorm(1 - 0.05/2) * pool.rubin$log.se.1ma
  
  rr.outcome<-subset(pool.rubin,select=c(log.rr.1ma,log.se.1ma,log.ci.lb.1ma,log.ci.ub.1ma))
  rr.outcome$rr.1ma<-exp(rr.outcome$log.rr.1ma)
  rr.outcome$rr.ci.lb.1ma<-exp(rr.outcome$log.ci.lb.1ma)
  rr.outcome$rr.ci.ub.1ma<-exp(rr.outcome$log.ci.ub.1ma)
  
  return(rr.outcome)
}

######################################################################################
#################################Load and prepare data################################
######################################################################################

setwd("/Users/jdamen/Documents/Julius/ZIKV analyses/2. Data")
load("imputation2.RData")
setwd("/Users/jdamen/Documents/GitHub/Zika_imputation/2_Code")
source("ZIKV prep.R")

########################Analyses#############


#one stage meta-analysis with log link and random intercept per study, per imputed dataset
fit1.coef<-as.data.frame(f.1ma.r.int(data,"microcephaly_bin"))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(data,fit1.coef)

  
  

c(exp(fit1.coef[1]),exp(fit1.coef[1]-1.96*fit1.coef[2]),exp(fit1.coef[1]+1.96*fit1.coef[2]))






##Objective 1: one-stage meta-analysis per outcome with log link and random intercept per study -> fit 1
#Microcephaly
fit1 <- glmer(microcephaly_bin ~ zikv_preg + (1 | studyname_fac), 
              data=data, family = binomial(link = "log"))
fit1.coef<-summary(fit1)$coefficients[2,]
c(exp(fit1.coef[1]),exp(fit1.coef[1]-1.96*fit1.coef[2]),exp(fit1.coef[1]+1.96*fit1.coef[2])) #RR
#Miscarriage
fit1 <- glmer(miscarriage ~ zikv_preg + (1 | studyname_fac), 
              data=data, family = binomial(link = "log"))
fit1.coef<-summary(fit1)$coefficients[2,]
c(exp(fit1.coef[1]),exp(fit1.coef[1]-1.96*fit1.coef[2]),exp(fit1.coef[1]+1.96*fit1.coef[2])) 
#Fetal loss
fit1 <- glmer(loss ~ zikv_preg + (1 | studyname_fac), 
              data=data, family = binomial(link = "log"))
fit1.coef<-summary(fit1)$coefficients[2,]
c(exp(fit1.coef[1]),exp(fit1.coef[1]-1.96*fit1.coef[2]),exp(fit1.coef[1]+1.96*fit1.coef[2])) 
#Congenital zika syndrome
fit1 <- glmer(czsn ~ zikv_preg + (1 | studyname_fac), 
              data=data, family = binomial(link = "log"))
fit1.coef<-summary(fit1)$coefficients[2,]
c(exp(fit1.coef[1]),exp(fit1.coef[1]-1.96*fit1.coef[2]),exp(fit1.coef[1]+1.96*fit1.coef[2])) 

##Objective 2: one-stage meta-analysis per outcome with log link and random effects on intercept and exposure -> fit 2
#Microcephaly
fit2<-glmer(microcephaly_bin ~ zikv_preg + (1+zikv_preg | studyname_fac), 
            data=data, family = binomial(link = "log"))
fit2.coef<-summary(fit2)$coefficients[2,]
c(exp(fit2.coef[1]),exp(fit2.coef[1]-1.96*fit2.coef[2]),exp(fit2.coef[1]+1.96*fit2.coef[2]))


##Fit 2 + all confounders
#Microcephaly
#Confounders: age, maritalstat
fit2c<-glmer(microcephaly_bin_birth ~ zikv_preg_cent + age + maritalstat + (1+zikv_preg_cent | studyname_fac), 
            data=data, family = binomial(link = "log"))
fit2c.coef<-as.data.frame(summary(fit2c)$coefficients[c(2:nrow(summary(fit2c)$coefficients)),c(1:2)])
colnames(fit2c.coef)<-c("log.rr", "log.rr.se")
fit2c.coef$rr<-exp(fit2c.coef$log.rr)
fit2c.coef$rr.lb<-exp(fit2c.coef$log.rr-1.96*fit2c.coef$log.rr.se)
fit2c.coef$rr.ub<-exp(fit2c.coef$log.rr+1.96*fit2c.coef$log.rr.se)
fit2c.coef


##Fit 2 + effect modifiers (one by one)

##Calculate absolute risks for fit 2: check protocol!

