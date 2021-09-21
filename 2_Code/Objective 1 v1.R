
library(dplyr)
library("xlsx")
library(mitml)
library(mice)
library(micemd)
library(metafor)
library(metamisc)

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
f.incidence <- function(variable,n.data) {
  a<-summary(variable)
  b<-prop.test(x = a[[2]], n = n.data, correct = TRUE)
  prop<-b$estimate[[1]]
  ci<-b$conf
    return(c(prop,ci))
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

setwd("C:/Users/jdamen/Documents/Julius centrum/ZIKV analyses/2. Data")
load("merged_imp4.RData")
data <- complete(merged_imp4, "long")
dataset.n<-length(data[data$.imp==1,]$studyname)
rm(merged_imp4)
m<-max(data$.imp)

#Create dichotomous outcome variables to calculate incidence
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
#Late fetal death (???28 weeks gestation)
data$lfdeath<-0
data$lfdeath[data$bdeath==1 & data$bdeath_ga>=28]<-1
data$lfdeath<-as.factor(data$lfdeath)
#Late fetal death (???28 weeks gestation) with microcephaly
data$lfdeath_micro<-0
data$lfdeath_micro[data$lfdeath==1 & data$bdeath_ga>=28]<-1
data$lfdeath<-as.factor(data$lfdeath)

#Calculate incidence of outcome per study
for (i in 1:length(unique(data$.imp))) {
  d<-data[data$.imp==i,]
  inc<-as.data.frame(d %>% group_by(studyname) %>% summarise(p=f.incidence(microcephaly,dataset.n)[1],ci.l<-f.incidence(microcephaly,dataset.n)[2],ci.u<-f.incidence(microcephaly,dataset.n)[3]))
  colnames(inc)<-c("studyname","incidence","ci.l","ci.u")
  if(i==1) {inc.micro<-inc}
  if(i>1) {inc.micro<-rbind(inc.micro,inc)}
}
inc.micro
#Replace 0's by 0.000001
inc.micro$incidence[inc.micro$incidence==0]<-0.000001
inc.micro$ci.l[inc.micro$ci.l==0]<-0.000001

#Pool with Rubins rules
inc.micro$logit.incidence<-logit(inc.micro$incidence)
inc.micro$logit.ci.l<-logit(inc.micro$ci.l)
inc.micro$logit.ci.u<-logit(inc.micro$ci.u)
inc.micro$logit.se<-(inc.micro$logit.ci.u-inc.micro$logit.ci.l)/(2*1.96)

pool.rubin<- data.frame(matrix(NA, nrow = length(unique(data$studyname)), ncol = 8))
colnames(pool.rubin)<-c("studyname","logit.abs","within","between","logit.var",
                                    "logit.se","logit.ci.lb","logit.ci.ub")
for (i in 1:max(inc.micro$studyname,na.rm=T)) {
  a<-inc.micro[inc.micro$studyname==i,]
  pool.rubin$studyname[i]<-i
  pool.rubin$logit.abs[i] <- mean(a$logit.incidence)
  pool.rubin$within[i] <- mean(a$logit.se^2)
  pool.rubin$between[i] <- (1 + (1/m)) * var(a$logit.incidence)
  pool.rubin$logit.var[i] <- pool.rubin$within[i] + pool.rubin$between[i]
  pool.rubin$logit.se[i] <- sqrt(pool.rubin$logit.var[i])
  pool.rubin$logit.ci.lb[i] <- pool.rubin$logit.abs[i] + qnorm(0.05/2)     * pool.rubin$logit.se[i]
  pool.rubin$logit.ci.ub[i] <- pool.rubin$logit.abs[i] + qnorm(1 - 0.05/2) * pool.rubin$logit.se[i]
}
pool.rubin

abs.microcephaly<-pool.rubin
abs.microcephaly$incidence<-inv.logit(abs.microcephaly$logit.abs)*100
abs.microcephaly$ci.lb<-inv.logit(abs.microcephaly$logit.ci.lb)*100
abs.microcephaly$ci.ub<-inv.logit(abs.microcephaly$logit.ci.ub)*100
abs.microcephaly

#Pool results: two-stage meta-analysis
fit.rma<-rma(yi = logit.abs, sei = logit.se, method = "REML", test = "knha", 
                       data = abs.microcephaly)
pool.microcephaly<-data.frame(matrix(NA, nrow = 1, ncol = 7))
colnames(pool.microcephaly)<-c("logit.abs","logit.se","logit.ci.lb","logit.ci.ub",
                               "abs.risk","ci.lb","ci.ub")
pool.microcephaly$logit.abs<-fit.rma$beta
pool.microcephaly$logit.se<-fit.rma$se
pool.microcephaly$logit.ci.lb<-fit.rma$ci.lb
pool.microcephaly$logit.ci.ub<-fit.rma$ci.ub
pool.microcephaly$abs.risk<-inv.logit(pool.microcephaly$logit.abs[[1]])*100
pool.microcephaly$ci.lb<-inv.logit(pool.microcephaly$logit.ci.lb)*100
pool.microcephaly$ci.ub<-inv.logit(pool.microcephaly$logit.ci.ub)*100
pool.microcephaly
#Prediction interval
PI<-PredInt(fit.rma)
PI<-cbind(inv.logit(fit.rma$b)[1,1], PI[1], PI[2])*100

#Forest plot
metafor::forest(abs.microcephaly$incidence, ci.lb=abs.microcephaly$ci.lb, ci.ub=abs.microcephaly$ci.ub, 
                refline = 0, slab = abs.microcephaly$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,13)), cex=1, steps=7, main="Microcephaly",
                xlim=(c(-4,12)), alim=(c(0,6)))
addpoly(x = pool.microcephaly$abs.risk, 
        sei = inv.logit.SE(pool.microcephaly$logit.se,pool.microcephaly$abs.risk), 
        rows=-0, cex=1)#Add pooled
addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
