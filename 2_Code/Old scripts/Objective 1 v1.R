
library(dplyr)
library("xlsx")
library(mitml)
library(mice)
library(micemd)
library(metafor)
library(metamisc)
library(EpiStats)

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
PredInt <- function(fit.rma)
{
  pi <- inv.logit(fit.rma$b + qt(c(0.025, 0.975), df=(fit.rma$k-2))*sqrt(fit.rma$tau2 + fit.rma$se**2))
  return(pi)
} #Function predict.rma uses k-1 degrees of freedom instead of k-2. Therefore written this function
PredIntOE <- function(fit.rma)
{
  pi <- exp(fit.rma$b + qt(c(0.025, 0.975), df=(fit.rma$k-2))*sqrt(fit.rma$tau2 + fit.rma$se**2))
  return(pi)
}

setwd("/Users/jdamen/Documents/Julius/ZIKV analyses/2. Data")
load("imputation2.RData")
data <- complete(merged_imp, "long")
dataset.n<-length(data[data$.imp==1,]$studyname)
rm(merged_imp)
m<-max(data$.imp)
studynames<-c("014-BRA","001-BRA","002-BRA","010-BRA","007-COL",
              "003-GUF","005-ESP","004-ESP","012-TTO","008-USA")

#Create dichotomous outcome variables to calculate incidence
#Microcephaly_bin
data$microcephaly_bin<-data$microcephaly
data$microcephaly_bin[data$microcephaly_bin==2]<-1
data$microcephaly_bin[data$microcephaly_bin==3]<-1
data$microcephaly_bin <- droplevels(data$microcephaly_bin)
#miscarriage (<20 weeks gestation)
data$miscarriage<-0
data$miscarriage[data$bdeath==1 & data$bdeath_ga<20]<-1
data$miscarriage<-as.factor(data$miscarriage)
#Fetal loss (>=20 weeks gestation)
data$loss<-0
data$loss[data$bdeath==1 & data$bdeath_ga>=20]<-1
data$loss<-as.factor(data$loss)
#Early fetal death (20-27 weeks gestation)
data$efdeath<-0
data$efdeath[data$bdeath==1 & data$bdeath_ga>=20 & data$bdeath_ga<28]<-1
data$efdeath<-as.factor(data$efdeath)
#Late fetal death (after 28 weeks gestation)
data$lfdeath<-0
data$lfdeath[data$bdeath==1 & data$bdeath_ga>=28]<-1
data$lfdeath<-as.factor(data$lfdeath)
#Late fetal death (after 28 weeks gestation) with microcephaly
data$lfdeath_micro<-0
data$lfdeath_micro[data$lfdeath==1 & data$bdeath_ga>=28]<-1
data$lfdeath<-as.factor(data$lfdeath)

data.zika<-data[data$zikv_preg==1,]
data.nozika<-data[data$zikv_preg==0,]

########################Analyses#############
#Change directory to save plots
setwd("/Users/jdamen/Documents/Julius/ZIKV analyses/4. Resultaten")

##################################################################################
###################################Microcephaly###################################
##################################################################################

#Calculate incidence of outcome per study
for (i in 1:length(unique(data.zika$.imp))) {
  d<-data.zika[data.zika$.imp==i,]
  inc<-as.data.frame(d %>% 
                       group_by(studyname) %>% 
                       summarise(p=f.incidence(microcephaly,dataset.n)[1],ci.l=f.incidence(microcephaly,dataset.n)[2],ci.u=f.incidence(microcephaly,dataset.n)[3]))
  if(i==1) {inc.outcome<-inc}
  if(i>1) {inc.outcome<-rbind(inc.outcome,inc)}
}
#Replace 0's by 0.000001
colnames(inc.outcome)<-c("studyname","incidence","ci.l","ci.u")
inc.outcome$incidence[inc.outcome$incidence==0]<-0.000001
inc.outcome$ci.l[inc.outcome$ci.l==0]<-0.000001

#Pool with Rubins rules
inc.outcome$logit.incidence<-logit(inc.outcome$incidence)
inc.outcome$logit.ci.l<-logit(inc.outcome$ci.l)
inc.outcome$logit.ci.u<-logit(inc.outcome$ci.u)
inc.outcome$logit.se<-(inc.outcome$logit.ci.u-inc.outcome$logit.ci.l)/(2*1.96)

pool.rubin<- data.frame(matrix(NA, nrow = length(unique(data.zika$studyname)), ncol = 8))
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

#Pool results: two-stage meta-analysis
fit.rma<-rma(yi = logit.abs, sei = logit.se, method = "REML", test = "knha", 
                       data = abs.outcome)
pool.outcome<-data.frame(matrix(NA, nrow = 1, ncol = 7))
colnames(pool.outcome)<-c("logit.abs","logit.se","logit.ci.lb","logit.ci.ub",
                               "abs.risk","ci.lb","ci.ub")
pool.outcome$logit.abs<-fit.rma$beta
pool.outcome$logit.se<-fit.rma$se
pool.outcome$logit.ci.lb<-fit.rma$ci.lb
pool.outcome$logit.ci.ub<-fit.rma$ci.ub
pool.outcome$abs.risk<-inv.logit(pool.outcome$logit.abs[[1]])*100
pool.outcome$ci.lb<-inv.logit(pool.outcome$logit.ci.lb)*100
pool.outcome$ci.ub<-inv.logit(pool.outcome$logit.ci.ub)*100
#Prediction interval
PI<-PredInt(fit.rma)
PI<-cbind(inv.logit(fit.rma$b)[1,1], PI[1], PI[2])*100

#Forest plot
#png(file="20210930 Microcephaly.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,13)), cex=1, steps=7, main="Microcephaly",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

########Relative risk########

for (i in 1:length(unique(data$.imp))) {
  d<-data[data$.imp==i,]
  for (j in 1:length(unique(d$studyname))) {
    d2<-d[d$studyname==j,]
    d2.rrtable<-table(d2$zikv_preg,d2$microcephaly_bin)
    #Add continuity correction?
    inc<-(c(unique(d2$studyname),f.rr(d2.rrtable)))
    if(i==1 & j==1) {rr.outcome<-inc} else {rr.outcome<-rbind(rr.outcome,inc)}
  }
}
rr.outcome
#Remove all rows for which RR could not be calculated
rr.outcome<-rr.outcome[!is.na(rr.outcome[,2]) & rr.outcome[,2]!=0 & !is.infinite(rr.outcome[,2]),]
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

#Pool results: two-stage meta-analysis
fit.rma<-rma(yi = log.rr, sei = log.se, method = "REML", test = "knha", 
             data = rr.outcome)
pool.outcome<-data.frame(matrix(NA, nrow = 1, ncol = 8))
colnames(pool.outcome)<-c("log.rr","log.se","log.ci.lb","log.ci.ub",
                          "rr","ci.lb","ci.ub","n.studies")
pool.outcome$log.rr<-fit.rma$beta
pool.outcome$log.se<-fit.rma$se
pool.outcome$log.ci.lb<-fit.rma$ci.lb
pool.outcome$log.ci.ub<-fit.rma$ci.ub
pool.outcome$rr<-exp(pool.outcome$log.rr[[1]])
pool.outcome$ci.lb<-exp(pool.outcome$log.ci.lb)
pool.outcome$ci.ub<-exp(pool.outcome$log.ci.ub)
pool.outcome$n.studies<-fit.rma$k

#Forest plot
#png(file="20210930 Microcephaly RR.png",width=750,height=500,res=100)
metafor::forest(rr.outcome$rr, ci.lb=rr.outcome$ci.lb, ci.ub=rr.outcome$ci.ub, 
                refline = 0, slab = rr.outcome$studyname,
                xlab = "Relative risk", pch = 19, psize=1,
                ylim=(c(-1,nrow(rr.outcome)+3)), cex=1, steps=7, main="Microcephaly",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome$rr, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

##################################################################################
###################################Miscarriage###################################
##################################################################################

#Calculate incidence of outcome per study
for (i in 1:length(unique(data.zika$.imp))) {
  d<-data.zika[data.zika$.imp==i,]
  inc<-as.data.frame(d %>% group_by(studyname) %>% summarise(p=f.incidence(miscarriage,dataset.n)[1],ci.l<-f.incidence(miscarriage,dataset.n)[2],ci.u<-f.incidence(miscarriage,dataset.n)[3]))
  colnames(inc)<-c("studyname","incidence","ci.l","ci.u")
  if(i==1) {inc.outcome<-inc}
  if(i>1) {inc.outcome<-rbind(inc.outcome,inc)}
}
#Replace 0's by 0.000001
inc.outcome$incidence[inc.outcome$incidence==0]<-0.000001
inc.outcome$ci.l[inc.outcome$ci.l==0]<-0.000001

#Pool with Rubins rules
inc.outcome$logit.incidence<-logit(inc.outcome$incidence)
inc.outcome$logit.ci.l<-logit(inc.outcome$ci.l)
inc.outcome$logit.ci.u<-logit(inc.outcome$ci.u)
inc.outcome$logit.se<-(inc.outcome$logit.ci.u-inc.outcome$logit.ci.l)/(2*1.96)

pool.rubin<- data.frame(matrix(NA, nrow = length(unique(data.zika$studyname)), ncol = 8))
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

#Pool results: two-stage meta-analysis
fit.rma<-rma(yi = logit.abs, sei = logit.se, method = "REML", test = "knha", 
             data = abs.outcome)
pool.outcome<-data.frame(matrix(NA, nrow = 1, ncol = 7))
colnames(pool.outcome)<-c("logit.abs","logit.se","logit.ci.lb","logit.ci.ub",
                          "abs.risk","ci.lb","ci.ub")
pool.outcome$logit.abs<-fit.rma$beta
pool.outcome$logit.se<-fit.rma$se
pool.outcome$logit.ci.lb<-fit.rma$ci.lb
pool.outcome$logit.ci.ub<-fit.rma$ci.ub
pool.outcome$abs.risk<-inv.logit(pool.outcome$logit.abs[[1]])*100
pool.outcome$ci.lb<-inv.logit(pool.outcome$logit.ci.lb)*100
pool.outcome$ci.ub<-inv.logit(pool.outcome$logit.ci.ub)*100
pool.outcome
#Prediction interval
PI<-PredInt(fit.rma)
PI<-cbind(inv.logit(fit.rma$b)[1,1], PI[1], PI[2])*100

#Forest plot
#png(file="20211025 Miscarriage.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,13)), cex=1, steps=7, main="Miscarriage",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

##################################################################################
###################################Fetal loss#####################################
##################################################################################

#Calculate incidence of outcome per study
for (i in 1:length(unique(data.zika$.imp))) {
  d<-data.zika[data.zika$.imp==i,]
  inc<-as.data.frame(d %>% group_by(studyname) %>% summarise(p=f.incidence(loss,dataset.n)[1],ci.l<-f.incidence(loss,dataset.n)[2],ci.u<-f.incidence(loss,dataset.n)[3]))
  colnames(inc)<-c("studyname","incidence","ci.l","ci.u")
  if(i==1) {inc.outcome<-inc}
  if(i>1) {inc.outcome<-rbind(inc.outcome,inc)}
}
#Replace 0's by 0.000001
inc.outcome$incidence[inc.outcome$incidence==0]<-0.000001
inc.outcome$ci.l[inc.outcome$ci.l==0]<-0.000001

#Pool with Rubins rules
inc.outcome$logit.incidence<-logit(inc.outcome$incidence)
inc.outcome$logit.ci.l<-logit(inc.outcome$ci.l)
inc.outcome$logit.ci.u<-logit(inc.outcome$ci.u)
inc.outcome$logit.se<-(inc.outcome$logit.ci.u-inc.outcome$logit.ci.l)/(2*1.96)

pool.rubin<- data.frame(matrix(NA, nrow = length(unique(data.zika$studyname)), ncol = 8))
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

#Pool results: two-stage meta-analysis
fit.rma<-rma(yi = logit.abs, sei = logit.se, method = "REML", test = "knha", 
             data = abs.outcome)
pool.outcome<-data.frame(matrix(NA, nrow = 1, ncol = 7))
colnames(pool.outcome)<-c("logit.abs","logit.se","logit.ci.lb","logit.ci.ub",
                          "abs.risk","ci.lb","ci.ub")
pool.outcome$logit.abs<-fit.rma$beta
pool.outcome$logit.se<-fit.rma$se
pool.outcome$logit.ci.lb<-fit.rma$ci.lb
pool.outcome$logit.ci.ub<-fit.rma$ci.ub
pool.outcome$abs.risk<-inv.logit(pool.outcome$logit.abs[[1]])*100
pool.outcome$ci.lb<-inv.logit(pool.outcome$logit.ci.lb)*100
pool.outcome$ci.ub<-inv.logit(pool.outcome$logit.ci.ub)*100
#Prediction interval
PI<-PredInt(fit.rma)
PI<-cbind(inv.logit(fit.rma$b)[1,1], PI[1], PI[2])*100

#Forest plot
#png(file="20210930 Fetal loss.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,13)), cex=1, steps=7, main="Fetal loss",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

##################################################################################
###################################Fetal loss#####################################
##################################################################################

#Calculate incidence of outcome per study
for (i in 1:length(unique(data.zika$.imp))) {
  d<-data.zika[data.zika$.imp==i,]
  inc<-as.data.frame(d %>% group_by(studyname) %>% summarise(p=f.incidence(czsn,dataset.n)[1],ci.l<-f.incidence(czsn,dataset.n)[2],ci.u<-f.incidence(czsn,dataset.n)[3]))
  colnames(inc)<-c("studyname","incidence","ci.l","ci.u")
  if(i==1) {inc.outcome<-inc}
  if(i>1) {inc.outcome<-rbind(inc.outcome,inc)}
}
#Replace 0's by 0.000001
inc.outcome$incidence[inc.outcome$incidence==0]<-0.000001
inc.outcome$ci.l[inc.outcome$ci.l==0]<-0.000001

#Pool with Rubins rules
inc.outcome$logit.incidence<-logit(inc.outcome$incidence)
inc.outcome$logit.ci.l<-logit(inc.outcome$ci.l)
inc.outcome$logit.ci.u<-logit(inc.outcome$ci.u)
inc.outcome$logit.se<-(inc.outcome$logit.ci.u-inc.outcome$logit.ci.l)/(2*1.96)

pool.rubin<- data.frame(matrix(NA, nrow = length(unique(data.zika$studyname)), ncol = 8))
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

#Pool results: two-stage meta-analysis
fit.rma<-rma(yi = logit.abs, sei = logit.se, method = "REML", test = "knha", 
             data = abs.outcome,control=list(maxiter=200))
pool.outcome<-data.frame(matrix(NA, nrow = 1, ncol = 7))
colnames(pool.outcome)<-c("logit.abs","logit.se","logit.ci.lb","logit.ci.ub",
                          "abs.risk","ci.lb","ci.ub")
pool.outcome$logit.abs<-fit.rma$beta
pool.outcome$logit.se<-fit.rma$se
pool.outcome$logit.ci.lb<-fit.rma$ci.lb
pool.outcome$logit.ci.ub<-fit.rma$ci.ub
pool.outcome$abs.risk<-inv.logit(pool.outcome$logit.abs[[1]])*100
pool.outcome$ci.lb<-inv.logit(pool.outcome$logit.ci.lb)*100
pool.outcome$ci.ub<-inv.logit(pool.outcome$logit.ci.ub)*100
#Prediction interval
PI<-PredInt(fit.rma)
PI<-cbind(inv.logit(fit.rma$b)[1,1], PI[1], PI[2])*100

#Forest plot
#png(file="20210930 Congenital Zika Syndrome.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,13)), cex=1, steps=7, main="Congenital Zika Syndrome",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Forest plot - resized
#png(file="20210930 Congenital Zika Syndrome - resized.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,13)), cex=1, steps=8, main="Congenital Zika Syndrome",
                xlim=(c(-45,55)), alim=(c(0,35)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()