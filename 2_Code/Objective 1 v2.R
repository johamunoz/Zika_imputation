
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
inv.logit.SE<- function(logit.se, c) {
  logit.se * (c*(1-c))
}
#Function for calculating absolute risk
f.incidence <- function(variable,n.data) {
  a<-summary(variable)
  b<-prop.test(x = a[[2]], n = n.data, correct = TRUE)
  prop<-b$estimate[[1]]
  ci<-b$conf
    return(c(prop,ci))
}
#Function for calculating relative risk
f.rr <- function(rr.table) {
  rr<-(rr.table[2,2]/sum(rr.table[2,]))/(rr.table[1,2]/sum(rr.table[1,]))
  log.rr<-log(rr)
  log.rr.lb<-log.rr-1.96*sqrt((rr.table[2,1]/rr.table[2,2])/sum(rr.table[2,]) +
                                (rr.table[1,1]/rr.table[1,2])/sum(rr.table[1,]))
  log.rr.ub<-log.rr+1.96*sqrt((rr.table[2,1]/rr.table[2,2])/sum(rr.table[2,]) +
                                (rr.table[1,1]/rr.table[1,2])/sum(rr.table[1,]))
  log.rr.se<-(log.rr.ub-log.rr.lb)/(2*1.96)
  rr.lb<-exp(log.rr.lb)
  rr.ub<-exp(log.rr.ub)
  return(c(rr,rr.lb,rr.ub,log.rr,log.rr.lb,log.rr.ub,log.rr.se))
}
#Function to calculate prediction interval using logit transformation
PredInt <- function(fit.rma)
{
  pi <- inv.logit(fit.rma$b + qt(c(0.025, 0.975), df=(fit.rma$k-2))*sqrt(fit.rma$tau2 + fit.rma$se**2))
  return(pi)
} #Function predict.rma uses k-1 degrees of freedom instead of k-2. Therefore written this function
#Function to calculate prediction interval using log transformation
PredIntLog <- function(fit.rma) {
  pi <- exp(fit.rma$b + qt(c(0.025, 0.975), df=(fit.rma$k-2))*sqrt(fit.rma$tau2 + fit.rma$se**2))
  return(pi)
}
#Function to calculate absolute risk of an outcome per study and per imputed dataset
f.abs.perstudy<-function(data, outcome_name) {
  
  dataset.n<-length(data[data$.imp==1,]$studyname)
  #Calculate absolute risk of outcome per study
  for (i in 1:length(unique(data$.imp))) {
    d<-data[data$.imp==i,]
    d$outcome <- d[, outcome_name]
    inc<-as.data.frame(d %>% 
                         group_by(studyname) %>% 
                         summarise(p=f.incidence(outcome,dataset.n)[1],
                                   ci.l=f.incidence(outcome,dataset.n)[2],
                                   ci.u=f.incidence(outcome,dataset.n)[3]))
    if(i==1) {inc.outcome<-inc}
    if(i>1) {inc.outcome<-rbind(inc.outcome,inc)}
  }
  #Replace 0's by 0.000001
  colnames(inc.outcome)<-c("studyname","incidence","ci.l","ci.u")
  inc.outcome$incidence[inc.outcome$incidence==0]<-0.000001
  inc.outcome$ci.l[inc.outcome$ci.l==0]<-0.000001
  
  #Logit transformation
  inc.outcome$logit.incidence<-logit(inc.outcome$incidence)
  inc.outcome$logit.ci.l<-logit(inc.outcome$ci.l)
  inc.outcome$logit.ci.u<-logit(inc.outcome$ci.u)
  inc.outcome$logit.se<-(inc.outcome$logit.ci.u-inc.outcome$logit.ci.l)/(2*1.96)
  inc.outcome
}

#Function to pool the absolute risks per study per imputed dataset using rubins rules
#Results in absolute risk per study
f.abs.poolrubin <-function(data,inc.outcome) {
  pool.rubin<- data.frame(matrix(NA, nrow = max(data$studyname), ncol = 8))
  colnames(pool.rubin)<-c("studyname","logit.abs","within","between","logit.var",
                          "logit.se","logit.ci.lb","logit.ci.ub")
  for (i in 1:max(inc.outcome$studyname,na.rm=T)) {
    a<-inc.outcome[inc.outcome$studyname==i,]
    pool.rubin$studyname[i]<-i
    pool.rubin$logit.abs[i] <- mean(a$logit.incidence)
    pool.rubin$within[i] <- mean(a$logit.se^2)
    pool.rubin$between[i] <- (1 + (1/m)) * var(a$logit.incidence)
    pool.rubin$logit.var[i] <- pool.rubin$within[i] + pool.rubin$between[i]
    pool.rubin$logit.se[i] <- sqrt(pool.rubin$logit.var[i])
    pool.rubin$logit.ci.lb[i] <- pool.rubin$logit.abs[i] + qnorm(0.05/2)     * pool.rubin$logit.se[i]
    pool.rubin$logit.ci.ub[i] <- pool.rubin$logit.abs[i] + qnorm(1 - 0.05/2) * pool.rubin$logit.se[i]
  }
  #Recode studyname
  pool.rubin$studyname<-as.factor(pool.rubin$studyname)
  levels(pool.rubin$studyname)<-studynames
  
  abs.outcome<-pool.rubin
  abs.outcome$incidence<-inv.logit(abs.outcome$logit.abs)*100
  abs.outcome$ci.lb<-inv.logit(abs.outcome$logit.ci.lb)*100
  abs.outcome$ci.ub<-inv.logit(abs.outcome$logit.ci.ub)*100
  
  return(abs.outcome)
}
#Function to perform a meta-analysis of the absolute risks per study
#Results in one overall pooled estimate of absolute risk
f.abs.2s.ma<-function(abs.outcome) {
  #Pool results: two-stage meta-analysis
  fit.rma<-rma(yi = logit.abs, sei = logit.se, method = "REML", test = "knha", 
               data = abs.outcome)
  pool.outcome<-data.frame(matrix(NA, nrow = 1, ncol = 9))
  colnames(pool.outcome)<-c("logit.abs","logit.se","logit.ci.lb","logit.ci.ub",
                            "abs.risk","ci.lb","ci.ub","pi.lb","pi.ub")
  pool.outcome$logit.abs<-fit.rma$beta
  pool.outcome$logit.se<-fit.rma$se
  pool.outcome$logit.ci.lb<-fit.rma$ci.lb
  pool.outcome$logit.ci.ub<-fit.rma$ci.ub
  pool.outcome$abs.risk<-inv.logit(pool.outcome$logit.abs[[1]])*100
  pool.outcome$ci.lb<-inv.logit(pool.outcome$logit.ci.lb)*100
  pool.outcome$ci.ub<-inv.logit(pool.outcome$logit.ci.ub)*100
  #Prediction interval
  PI<-PredInt(fit.rma)
  pool.outcome$pi.lb<-PI[1]*100
  pool.outcome$pi.ub<-PI[2]*100
  return(pool.outcome)
}

#Function to calculate relative risk of an outcome per study and per imputed dataset
f.rel.perstudy<-function(data, outcome_name) {
  for (i in 1:length(unique(data$.imp))) {
    d<-data[data$.imp==i,]
    for (j in 1:length(unique(d$studyname))) {
      d2<-d[d$studyname==j,]
      d2.rrtable<-table(d2$zikv_preg,d2[[outcome_name]])
      #Add continuity correction?
      inc<-(c(unique(d2$studyname),f.rr(d2.rrtable)))
      if(i==1 & j==1) {rr.outcome.all<-inc} else {rr.outcome.all<-rbind(rr.outcome.all,inc)}
    }
  }
  rr.outcome.all
}

#Function to pool the relative risks per study per imputed dataset using rubins rules
#Results in relative risk per study
f.rel.poolrubin <-function(data,rr.outcome.all) {
  #Remove all rows for which RR could not be calculated
  rr.outcome<-rr.outcome.all[!is.na(rr.outcome.all[,2]) & rr.outcome.all[,2]!=0 & !is.infinite(rr.outcome.all[,2]),]
  colnames(rr.outcome)<-c("studyname","rr","ci.l","ci.u","log.rr","log.ci.l","log.ci.u","log.se")
  rr.outcome<-as.data.frame(rr.outcome)
  #Pool with Rubins rules
  pool.rubin<- data.frame(matrix(NA, nrow = length(unique(data$studyname)), ncol = 8))
  colnames(pool.rubin)<-c("studyname","log.rr","within","between","log.var",
                          "log.se","log.ci.lb","log.ci.ub")
  for (i in 1:max(rr.outcome$studyname,na.rm=T)) {
    a<-rr.outcome[rr.outcome$studyname==i,]
    pool.rubin$studyname[i]<-i
    pool.rubin$log.rr[i] <- mean(a$log.rr)
    pool.rubin$within[i] <- mean(a$log.se^2)
    pool.rubin$between[i] <- (1 + (1/m)) * var(a$log.rr)
    pool.rubin$log.var[i] <- pool.rubin$within[i] + pool.rubin$between[i]
    pool.rubin$log.se[i] <- sqrt(pool.rubin$log.var[i])
    pool.rubin$log.ci.lb[i] <- pool.rubin$log.rr[i] + qnorm(0.05/2)     * pool.rubin$log.se[i]
    pool.rubin$log.ci.ub[i] <- pool.rubin$log.rr[i] + qnorm(1 - 0.05/2) * pool.rubin$log.se[i]
  }
  
  #Recode studyname
  pool.rubin$studyname<-as.factor(pool.rubin$studyname)
  levels(pool.rubin$studyname)<-studynames
  
  rr.outcome<-pool.rubin
  rr.outcome$rr<-exp(rr.outcome$log.rr)
  rr.outcome$ci.lb<-exp(rr.outcome$log.ci.lb)
  rr.outcome$ci.ub<-exp(rr.outcome$log.ci.ub)
  
  #Remove studies for which rr could not be calculated
  rr.outcome<-rr.outcome[!is.na(rr.outcome$log.rr),]
  
  return(rr.outcome)
}

#Function to perform a meta-analysis of the relative risks per study
#Results in one overall pooled estimate of relative risk
f.rel.2s.ma<-function(rr.outcome) {
  #Pool results: two-stage meta-analysis
  fit.rma<-rma(yi = log.rr, sei = log.se, method = "REML", test = "knha", 
               data = rr.outcome)
  pool.outcome.rr<-data.frame(matrix(NA, nrow = 1, ncol = 8))
  colnames(pool.outcome.rr)<-c("log.rr","log.rr.se","log.rr.ci.lb","log.rr.ci.ub",
                               "rr","rr.ci.lb","rr.ci.ub","n.studies")
  pool.outcome.rr$log.rr<-fit.rma$beta
  pool.outcome.rr$log.rr.se<-fit.rma$se
  pool.outcome.rr$log.rr.ci.lb<-fit.rma$ci.lb
  pool.outcome.rr$log.rr.ci.ub<-fit.rma$ci.ub
  pool.outcome.rr$rr<-exp(pool.outcome.rr$log.rr[[1]])
  pool.outcome.rr$rr.ci.lb<-exp(pool.outcome.rr$log.rr.ci.lb)
  pool.outcome.rr$rr.ci.ub<-exp(pool.outcome.rr$log.rr.ci.ub)
  pool.outcome.rr$n.studies<-fit.rma$k
  return(pool.outcome.rr)
}

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
#Change directory to save plots
setwd("/Users/jdamen/Documents/Julius/ZIKV analyses/4. Resultaten")

#Data frame to store results
all.results<-data.frame(Outcome=character(),
                        File=numeric(), 
                        User=numeric(), 
                        stringsAsFactors=FALSE)

##################################################################################
###################################Microcephaly###################################
##################################################################################

#Zika-positive women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.zika,"microcephaly_bin")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.zika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20220328 Microcephaly zika positive.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Microcephaly, zika-positive women",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-pool.outcome

#Zika-negative women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.nozika,"microcephaly_bin")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.nozika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20220328 Microcephaly zika negative.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Microcephaly, zika-negative women",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
colnames(pool.outcome)<-c("logit.abs.neg","logit.se.neg","logit.ci.lb.neg","logit.ci.ub.neg",
                          "abs.risk.neg","ci.lb.neg","ci.ub.neg","pi.lb.neg","pi.ub.neg")
results<-cbind(results,pool.outcome)

########Relative risk########

####Two-stage meta-analysis
#calculate relative risk of outcome per study and per imputed dataset
rr.outcome.all<-f.rel.perstudy(data,"microcephaly_bin")
#pool the relative risks per study per imputed dataset using rubins rules, resulting in RR per study
rr.outcome<-f.rel.poolrubin(data,rr.outcome.all)
#Pool the relative risks over the studies, resulting in one pooled summary relative risk over all studies
pool.outcome.rr<-f.rel.2s.ma(rr.outcome)

#Forest plot
#png(file="20220328 Microcephaly RR.png",width=750,height=500,res=100)
metafor::forest(rr.outcome$rr, ci.lb=rr.outcome$ci.lb, ci.ub=rr.outcome$ci.ub, 
                refline = 0, slab = rr.outcome$studyname,
                xlab = "Relative risk", pch = 19, psize=1,
                ylim=(c(-1,nrow(rr.outcome)+3)), cex=1, steps=7, main="Microcephaly",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome.rr$rr, ci.lb=pool.outcome.rr$rr.ci.lb, ci.ub=pool.outcome.rr$rr.ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome.rr)

######One-stage meta-analysis relative risk
#One stage meta-analysis with log link and random intercept per study, per imputed dataset
fit1.coef<-as.data.frame(f.1ma.r.int(data,"microcephaly_bin"))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(data,fit1.coef)

#Store results
results<-cbind(results,pool.outcome.rr.1ma)
results$outcome<-"Microcephaly"
all.results<-results

##################################################################################
###################################Miscarriage###################################
##################################################################################

#Zika-positive women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.zika,"miscarriage")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.zika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20220117 Miscarriage zika positive.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Miscarriage, zika-positive women",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-pool.outcome

#Zika-negative women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.nozika,"miscarriage")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.nozika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20220117 Miscarriage zika negative.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Miscarriage, zika-negative women",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome)

########Relative risk########

#calculate relative risk of outcome per study and per imputed dataset
rr.outcome.all<-f.rel.perstudy(data,"miscarriage")
#pool the relative risks per study per imputed dataset using rubins rules, resulting in RR per study
rr.outcome<-f.rel.poolrubin(data,rr.outcome.all)
#Pool the relative risks over the studies, resulting in one pooled summary relative risk over all studies
pool.outcome.rr<-f.rel.2s.ma(rr.outcome)

#Forest plot
#png(file="20220117 Miscarriage RR.png",width=750,height=500,res=100)
metafor::forest(rr.outcome$rr, ci.lb=rr.outcome$ci.lb, ci.ub=rr.outcome$ci.ub, 
                refline = 0, slab = rr.outcome$studyname,
                xlab = "Relative risk", pch = 19, psize=1,
                ylim=(c(-1,nrow(rr.outcome)+3)), cex=1, steps=7, main="Miscarriage",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome.rr$rr, ci.lb=pool.outcome.rr$rr.ci.lb, ci.ub=pool.outcome.rr$rr.ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome.rr)

######One-stage meta-analysis relative risk
#One stage meta-analysis with log link and random intercept per study, per imputed dataset
fit1.coef<-as.data.frame(f.1ma.r.int(data,"miscarriage"))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(data,fit1.coef)

#Store results
results<-cbind(results,pool.outcome.rr.1ma)
results$outcome<-"Miscarriage"
all.results<-rbind(all.results,results)

##################################################################################
###################################Fetal loss#####################################
##################################################################################

#Zika-positive women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.zika,"loss")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.zika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20220117 Fetal loss zika positive.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Fetal loss, zika-positive women",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-pool.outcome

#Zika-negative women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.nozika,"loss")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.nozika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20220117 Fetal loss zika negative.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Fetal loss, zika-negative women",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome)

########Relative risk########

#calculate relative risk of outcome per study and per imputed dataset
rr.outcome.all<-f.rel.perstudy(data,"loss")
#pool the relative risks per study per imputed dataset using rubins rules, resulting in RR per study
rr.outcome<-f.rel.poolrubin(data,rr.outcome.all)
#Pool the relative risks over the studies, resulting in one pooled summary relative risk over all studies
pool.outcome.rr<-f.rel.2s.ma(rr.outcome)

#Forest plot
#png(file="20220117 Fetal loss RR.png",width=750,height=500,res=100)
metafor::forest(rr.outcome$rr, ci.lb=rr.outcome$ci.lb, ci.ub=rr.outcome$ci.ub, 
                refline = 0, slab = rr.outcome$studyname,
                xlab = "Relative risk", pch = 19, psize=1,
                ylim=(c(-1,nrow(rr.outcome)+3)), cex=1, steps=7, main="Fetal loss",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome.rr$rr, ci.lb=pool.outcome.rr$rr.ci.lb, ci.ub=pool.outcome.rr$rr.ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome.rr)

######One-stage meta-analysis relative risk
#One stage meta-analysis with log link and random intercept per study, per imputed dataset
fit1.coef<-as.data.frame(f.1ma.r.int(data,"loss"))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(data,fit1.coef)

#Store results
results<-cbind(results,pool.outcome.rr.1ma)
results$outcome<-"Fetal loss"
all.results<-rbind(all.results,results)

##################################################################################
########################Congenital zika syndrome##################################
##################################################################################

#Zika-positive women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.zika,"czsn")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.zika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20220117 czs zika positive.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Congenital zika syndrome, zika-positive women",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-pool.outcome

#Forest plot - resized
#png(file="20210930 Congenital Zika Syndrome - resized.png",width=750,height=500,res=100)
#metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
#                refline = 0, slab = abs.outcome$studyname,
#                xlab = "Absolute risk (%)", pch = 19, psize=1,
#                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=8, main="Congenital Zika Syndrome",
#                xlim=(c(-45,55)), alim=(c(0,35)))
#addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
#        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Zika-negative women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.nozika,"czsn")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.nozika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20220117 czs zika negative.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Congenital zika syndrome, zika-negative women",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome)

########Relative risk########

#calculate relative risk of outcome per study and per imputed dataset
rr.outcome.all<-f.rel.perstudy(data,"czsn")
#pool the relative risks per study per imputed dataset using rubins rules, resulting in RR per study
rr.outcome<-f.rel.poolrubin(data,rr.outcome.all)
#Pool the relative risks over the studies, resulting in one pooled summary relative risk over all studies
pool.outcome.rr<-f.rel.2s.ma(rr.outcome)

#Forest plot
#png(file="20220117 czs RR.png",width=750,height=500,res=100)
metafor::forest(rr.outcome$rr, ci.lb=rr.outcome$ci.lb, ci.ub=rr.outcome$ci.ub, 
                refline = 0, slab = rr.outcome$studyname,
                xlab = "Relative risk", pch = 19, psize=1,
                ylim=(c(-1,nrow(rr.outcome)+3)), cex=1, steps=7, main="Congenital zika syndrome",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome.rr$rr, ci.lb=pool.outcome.rr$rr.ci.lb, ci.ub=pool.outcome.rr$rr.ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome.rr)

######One-stage meta-analysis relative risk
#One stage meta-analysis with log link and random intercept per study, per imputed dataset
fit1.coef<-as.data.frame(f.1ma.r.int(data,"czsn"))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(data,fit1.coef)

#Store results
results<-cbind(results,pool.outcome.rr.1ma)
results$outcome<-"CZS"
all.results<-rbind(all.results,results)

##################################################################################
################################Induced abortion##################################
##################################################################################

#Zika-positive women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.zika,"inducedabort")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.zika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20220117 Induced abortion zika positive.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Induced abortion, zika-positive women",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-pool.outcome

#Zika-negative women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.nozika,"inducedabort")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.nozika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20220117 Induced abortion zika negative.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Induced abortion, zika-negative women",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome)

########Relative risk########

#calculate relative risk of outcome per study and per imputed dataset
rr.outcome.all<-f.rel.perstudy(data,"inducedabort")
#pool the relative risks per study per imputed dataset using rubins rules, resulting in RR per study
rr.outcome<-f.rel.poolrubin(data,rr.outcome.all)
#Pool the relative risks over the studies, resulting in one pooled summary relative risk over all studies
pool.outcome.rr<-f.rel.2s.ma(rr.outcome)

#Forest plot
#png(file="20220117 Induced abortion RR.png",width=750,height=500,res=100)
metafor::forest(rr.outcome$rr, ci.lb=rr.outcome$ci.lb, ci.ub=rr.outcome$ci.ub, 
                refline = 0, slab = rr.outcome$studyname,
                xlab = "Relative risk", pch = 19, psize=1,
                ylim=(c(-1,nrow(rr.outcome)+3)), cex=1, steps=7, main="Induced abortion",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome.rr$rr, ci.lb=pool.outcome.rr$rr.ci.lb, ci.ub=pool.outcome.rr$rr.ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome.rr)

######One-stage meta-analysis relative risk
#One stage meta-analysis with log link and random intercept per study, per imputed dataset
fit1.coef<-as.data.frame(f.1ma.r.int(data,"inducedabort"))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(data,fit1.coef)

#Store results
results<-cbind(results,pool.outcome.rr.1ma)
results$outcome<-"Induced abortion"
all.results<-rbind(all.results,results)

##################################################################################
################################Early fetal death#################################
##################################################################################

#Zika-positive women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.zika,"efdeath")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.zika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20220117 Early fetal death zika positive.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Early fetal death, zika-positive women",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-pool.outcome

#Zika-negative women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.nozika,"efdeath")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.nozika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20220117 Early fetal death zika negative.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Early fetal death, zika-negative women",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome)

########Relative risk########

#calculate relative risk of outcome per study and per imputed dataset
rr.outcome.all<-f.rel.perstudy(data,"efdeath")
#pool the relative risks per study per imputed dataset using rubins rules, resulting in RR per study
rr.outcome<-f.rel.poolrubin(data,rr.outcome.all)
#Pool the relative risks over the studies, resulting in one pooled summary relative risk over all studies
pool.outcome.rr<-f.rel.2s.ma(rr.outcome)

#Forest plot
#png(file="20220117 Early fetal death RR.png",width=750,height=500,res=100)
metafor::forest(rr.outcome$rr, ci.lb=rr.outcome$ci.lb, ci.ub=rr.outcome$ci.ub, 
                refline = 0, slab = rr.outcome$studyname,
                xlab = "Relative risk", pch = 19, psize=1,
                ylim=(c(-1,nrow(rr.outcome)+3)), cex=1, steps=7, main="Early fetal death",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome.rr$rr, ci.lb=pool.outcome.rr$rr.ci.lb, ci.ub=pool.outcome.rr$rr.ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome.rr)

######One-stage meta-analysis relative risk
#One stage meta-analysis with log link and random intercept per study, per imputed dataset
fit1.coef<-as.data.frame(f.1ma.r.int(data,"efdeath"))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(data,fit1.coef)

#Store results
results<-cbind(results,pool.outcome.rr.1ma)
results$outcome<-"Early fetal death"
all.results<-rbind(all.results,results)

##################################################################################
#################################Late fetal death#################################
##################################################################################

#Zika-positive women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.zika,"lfdeath")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.zika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20220117 Late fetal death zika positive.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Late fetal death, zika-positive women",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-pool.outcome

#Zika-negative women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.nozika,"lfdeath")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.nozika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20220117 Late fetal death zika negative.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Late fetal death, zika-negative women",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome)

########Relative risk########

#calculate relative risk of outcome per study and per imputed dataset
rr.outcome.all<-f.rel.perstudy(data,"lfdeath")
#pool the relative risks per study per imputed dataset using rubins rules, resulting in RR per study
rr.outcome<-f.rel.poolrubin(data,rr.outcome.all)
#Pool the relative risks over the studies, resulting in one pooled summary relative risk over all studies
pool.outcome.rr<-f.rel.2s.ma(rr.outcome)

#Forest plot
#png(file="20220117 Late fetal death RR.png",width=750,height=500,res=100)
metafor::forest(rr.outcome$rr, ci.lb=rr.outcome$ci.lb, ci.ub=rr.outcome$ci.ub, 
                refline = 0, slab = rr.outcome$studyname,
                xlab = "Relative risk", pch = 19, psize=1,
                ylim=(c(-1,nrow(rr.outcome)+3)), cex=1, steps=7, main="Late fetal death",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome.rr$rr, ci.lb=pool.outcome.rr$rr.ci.lb, ci.ub=pool.outcome.rr$rr.ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome.rr)

######One-stage meta-analysis relative risk
#One stage meta-analysis with log link and random intercept per study, per imputed dataset
fit1.coef<-as.data.frame(f.1ma.r.int(data,"lfdeath"))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(data,fit1.coef)

#Store results
results<-cbind(results,pool.outcome.rr.1ma)
results$outcome<-"Late fetal death"
all.results<-rbind(all.results,results)

##################################################################################
###########################Intrauterine growth restriction########################
##################################################################################

#Zika-positive women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.zika,"igr_curr_preg")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.zika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20220117 IGR zika positive.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Intrauterine growth restriction, zika-positive women",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-pool.outcome

#Zika-negative women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.nozika,"igr_curr_preg")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.nozika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20220117 IGR zika negative.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Intrauterine growth restriction, zika-negative women",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome)

########Relative risk########

#calculate relative risk of outcome per study and per imputed dataset
rr.outcome.all<-f.rel.perstudy(data,"igr_curr_preg")
#pool the relative risks per study per imputed dataset using rubins rules, resulting in RR per study
rr.outcome<-f.rel.poolrubin(data,rr.outcome.all)
#Pool the relative risks over the studies, resulting in one pooled summary relative risk over all studies
pool.outcome.rr<-f.rel.2s.ma(rr.outcome)

#Forest plot
#png(file="20220117 IGR RR.png",width=750,height=500,res=100)
metafor::forest(rr.outcome$rr, ci.lb=rr.outcome$ci.lb, ci.ub=rr.outcome$ci.ub, 
                refline = 0, slab = rr.outcome$studyname,
                xlab = "Relative risk", pch = 19, psize=1,
                ylim=(c(-1,nrow(rr.outcome)+3)), cex=1, steps=7, main="Intrauterine growth restriction",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome.rr$rr, ci.lb=pool.outcome.rr$rr.ci.lb, ci.ub=pool.outcome.rr$rr.ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome.rr)

######One-stage meta-analysis relative risk
#One stage meta-analysis with log link and random intercept per study, per imputed dataset
fit1.coef<-as.data.frame(f.1ma.r.int(data,"igr_curr_preg"))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(data,fit1.coef)

#Store results
results<-cbind(results,pool.outcome.rr.1ma)
results$outcome<-"Intrauterine growth restriction"
all.results<-rbind(all.results,results)

##################################################################################
################################Craniofacial disproportion########################
##################################################################################

#Zika-positive women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.zika,"inf_craniofac_abn_bin")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.zika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20220117 Craniofacial zika positive.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Craniofacial disproportion, zika-positive women",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-pool.outcome

#Zika-negative women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.nozika,"inf_craniofac_abn_bin")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.nozika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20220117 Craniofacial zika negative.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Craniofacial disproportion, zika-negative women",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome)

########Relative risk########

#calculate relative risk of outcome per study and per imputed dataset
rr.outcome.all<-f.rel.perstudy(data,"inf_craniofac_abn_bin")
#pool the relative risks per study per imputed dataset using rubins rules, resulting in RR per study
rr.outcome<-f.rel.poolrubin(data,rr.outcome.all)
#Pool the relative risks over the studies, resulting in one pooled summary relative risk over all studies
pool.outcome.rr<-f.rel.2s.ma(rr.outcome)

#Forest plot
#png(file="20220117 Craniofacial RR.png",width=750,height=500,res=100)
metafor::forest(rr.outcome$rr, ci.lb=rr.outcome$ci.lb, ci.ub=rr.outcome$ci.ub, 
                refline = 0, slab = rr.outcome$studyname,
                xlab = "Relative risk", pch = 19, psize=1,
                ylim=(c(-1,nrow(rr.outcome)+3)), cex=1, steps=7, main="Craniofacial disproportion",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome.rr$rr, ci.lb=pool.outcome.rr$rr.ci.lb, ci.ub=pool.outcome.rr$rr.ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome.rr)

######One-stage meta-analysis relative risk
#One stage meta-analysis with log link and random intercept per study, per imputed dataset
fit1.coef<-as.data.frame(f.1ma.r.int(data,"inf_craniofac_abn_bin"))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(data,fit1.coef)

#Store results
results<-cbind(results,pool.outcome.rr.1ma)
results$outcome<-"Craniofacial disproportion"
all.results<-rbind(all.results,results)

##################################################################################
################################Neuroimaging abnormalities########################
##################################################################################

#Zika-positive women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.zika,"neuroabnormality")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.zika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20220117 Neuroabnormality zika positive.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Neuroimaging abnormalities, zika-positive women",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-pool.outcome

#Zika-negative women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.nozika,"neuroabnormality")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.nozika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20220117 Neuroabnormality zika negative.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Neuroimaging abnormalities, zika-negative women",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome)

########Relative risk########

#calculate relative risk of outcome per study and per imputed dataset
rr.outcome.all<-f.rel.perstudy(data,"neuroabnormality")
#pool the relative risks per study per imputed dataset using rubins rules, resulting in RR per study
rr.outcome<-f.rel.poolrubin(data,rr.outcome.all)
#Pool the relative risks over the studies, resulting in one pooled summary relative risk over all studies
pool.outcome.rr<-f.rel.2s.ma(rr.outcome)

#Forest plot
#png(file="20220117 Neuroabnormality RR.png",width=750,height=500,res=100)
metafor::forest(rr.outcome$rr, ci.lb=rr.outcome$ci.lb, ci.ub=rr.outcome$ci.ub, 
                refline = 0, slab = rr.outcome$studyname,
                xlab = "Relative risk", pch = 19, psize=1,
                ylim=(c(-1,nrow(rr.outcome)+3)), cex=1, steps=7, main="Neuroimaging abnormalities",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome.rr$rr, ci.lb=pool.outcome.rr$rr.ci.lb, ci.ub=pool.outcome.rr$rr.ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome.rr)

######One-stage meta-analysis relative risk
#One stage meta-analysis with log link and random intercept per study, per imputed dataset
fit1.coef<-as.data.frame(f.1ma.r.int(data,"neuroabnormality"))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(data,fit1.coef)

#Store results
results<-cbind(results,pool.outcome.rr.1ma)
results$outcome<-"Neuroimaging abnormalities"
all.results<-rbind(all.results,results)

######################
#write.xlsx2(all.results, "20220117 Results objective 1.xls", col.names = T, row.names = F)
