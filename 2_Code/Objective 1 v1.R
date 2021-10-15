
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

setwd("/Users/jdamen/Documents/Julius/ZIKV analyses/2. Data")
load("imputation2.RData")
data <- complete(merged_imp, "long")
dataset.n<-length(data[data$.imp==1,]$studyname)
rm(merged_imp)
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

########################Analyses#############
#Change directory to save plots
setwd("/Users/jdamen/Documents/Julius/ZIKV analyses/4. Resultaten")

##################################################################################
###################################Microcephaly###################################
##################################################################################

#Calculate incidence of outcome per study
for (i in 1:length(unique(data$.imp))) {
  d<-data[data$.imp==i,]
  inc<-as.data.frame(d %>% group_by(studyname) %>% summarise(p=f.incidence(microcephaly,dataset.n)[1],ci.l<-f.incidence(microcephaly,dataset.n)[2],ci.u<-f.incidence(microcephaly,dataset.n)[3]))
  colnames(inc)<-c("studyname","incidence","ci.l","ci.u")
  if(i==1) {inc.outcome<-inc}
  if(i>1) {inc.outcome<-rbind(inc.outcome,inc)}
}
inc.outcome
#Replace 0's by 0.000001
inc.outcome$incidence[inc.outcome$incidence==0]<-0.000001
inc.outcome$ci.l[inc.outcome$ci.l==0]<-0.000001

#Pool with Rubins rules
inc.outcome$logit.incidence<-logit(inc.outcome$incidence)
inc.outcome$logit.ci.l<-logit(inc.outcome$ci.l)
inc.outcome$logit.ci.u<-logit(inc.outcome$ci.u)
inc.outcome$logit.se<-(inc.outcome$logit.ci.u-inc.outcome$logit.ci.l)/(2*1.96)

pool.rubin<- data.frame(matrix(NA, nrow = length(unique(data$studyname)), ncol = 8))
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
pool.rubin$studyname<-as.character(pool.rubin$studyname)
pool.rubin$studyname[pool.rubin$studyname=="1"]<-"Brazil_BahiaPaudaLima_Costa"
pool.rubin$studyname[pool.rubin$studyname=="2"]<-"Brazil_RiodeJaneiro_CunhaPrata"
pool.rubin$studyname[pool.rubin$studyname=="3"]<-"Brazil_RiodeJaneiro_Joao"
pool.rubin$studyname[pool.rubin$studyname=="4"]<-"Brazil_SP_RibeiraoPreto_Duarte"
pool.rubin$studyname[pool.rubin$studyname=="5"]<-"Colombia_Mulkey"
pool.rubin$studyname[pool.rubin$studyname=="6"]<-"FrenchGuiana_Pomar"
pool.rubin$studyname[pool.rubin$studyname=="7"]<-"Spain_Bardaji"
pool.rubin$studyname[pool.rubin$studyname=="8"]<-"Spain_Soriano"
pool.rubin$studyname[pool.rubin$studyname=="9"]<-"TrinidadTobago_Sohan"
pool.rubin$studyname[pool.rubin$studyname=="10"]<-"USA_Mulkey"

pool.rubin

abs.outcome<-pool.rubin
abs.outcome$incidence<-inv.logit(abs.outcome$logit.abs)*100
abs.outcome$ci.lb<-inv.logit(abs.outcome$logit.ci.lb)*100
abs.outcome$ci.ub<-inv.logit(abs.outcome$logit.ci.ub)*100
abs.outcome

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
#png(file="20210930 Microcephaly.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,13)), cex=1, steps=7, main="Microcephaly",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, 
        sei = inv.logit.SE(pool.outcome$logit.se,pool.outcome$abs.risk), 
        rows=-0, cex=1)#Add pooled
addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()


##################################################################################
###################################Miscarriage###################################
##################################################################################

#Calculate incidence of outcome per study
for (i in 1:length(unique(data$.imp))) {
  d<-data[data$.imp==i,]
  inc<-as.data.frame(d %>% group_by(studyname) %>% summarise(p=f.incidence(miscarriage,dataset.n)[1],ci.l<-f.incidence(miscarriage,dataset.n)[2],ci.u<-f.incidence(miscarriage,dataset.n)[3]))
  colnames(inc)<-c("studyname","incidence","ci.l","ci.u")
  if(i==1) {inc.outcome<-inc}
  if(i>1) {inc.outcome<-rbind(inc.outcome,inc)}
}
inc.outcome
#Replace 0's by 0.000001
inc.outcome$incidence[inc.outcome$incidence==0]<-0.000001
inc.outcome$ci.l[inc.outcome$ci.l==0]<-0.000001

#Pool with Rubins rules
inc.outcome$logit.incidence<-logit(inc.outcome$incidence)
inc.outcome$logit.ci.l<-logit(inc.outcome$ci.l)
inc.outcome$logit.ci.u<-logit(inc.outcome$ci.u)
inc.outcome$logit.se<-(inc.outcome$logit.ci.u-inc.outcome$logit.ci.l)/(2*1.96)

pool.rubin<- data.frame(matrix(NA, nrow = length(unique(data$studyname)), ncol = 8))
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
pool.rubin$studyname<-as.character(pool.rubin$studyname)
pool.rubin$studyname[pool.rubin$studyname=="1"]<-"Brazil_BahiaPaudaLima_Costa"
pool.rubin$studyname[pool.rubin$studyname=="2"]<-"Brazil_RiodeJaneiro_CunhaPrata"
pool.rubin$studyname[pool.rubin$studyname=="3"]<-"Brazil_RiodeJaneiro_Joao"
pool.rubin$studyname[pool.rubin$studyname=="4"]<-"Brazil_SP_RibeiraoPreto_Duarte"
pool.rubin$studyname[pool.rubin$studyname=="5"]<-"Colombia_Mulkey"
pool.rubin$studyname[pool.rubin$studyname=="6"]<-"FrenchGuiana_Pomar"
pool.rubin$studyname[pool.rubin$studyname=="7"]<-"Spain_Bardaji"
pool.rubin$studyname[pool.rubin$studyname=="8"]<-"Spain_Soriano"
pool.rubin$studyname[pool.rubin$studyname=="9"]<-"TrinidadTobago_Sohan"
pool.rubin$studyname[pool.rubin$studyname=="10"]<-"USA_Mulkey"

pool.rubin

abs.outcome<-pool.rubin
abs.outcome$incidence<-inv.logit(abs.outcome$logit.abs)*100
abs.outcome$ci.lb<-inv.logit(abs.outcome$logit.ci.lb)*100
abs.outcome$ci.ub<-inv.logit(abs.outcome$logit.ci.ub)*100
abs.outcome

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
#png(file="20210930 Miscarriage.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,13)), cex=1, steps=7, main="Miscarriage",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, 
        sei = inv.logit.SE(pool.outcome$logit.se,pool.outcome$abs.risk), 
        rows=-0, cex=1)#Add pooled
addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

##################################################################################
###################################Fetal loss#####################################
##################################################################################

#Calculate incidence of outcome per study
for (i in 1:length(unique(data$.imp))) {
  d<-data[data$.imp==i,]
  inc<-as.data.frame(d %>% group_by(studyname) %>% summarise(p=f.incidence(loss,dataset.n)[1],ci.l<-f.incidence(loss,dataset.n)[2],ci.u<-f.incidence(loss,dataset.n)[3]))
  colnames(inc)<-c("studyname","incidence","ci.l","ci.u")
  if(i==1) {inc.outcome<-inc}
  if(i>1) {inc.outcome<-rbind(inc.outcome,inc)}
}
inc.outcome
#Replace 0's by 0.000001
inc.outcome$incidence[inc.outcome$incidence==0]<-0.000001
inc.outcome$ci.l[inc.outcome$ci.l==0]<-0.000001

#Pool with Rubins rules
inc.outcome$logit.incidence<-logit(inc.outcome$incidence)
inc.outcome$logit.ci.l<-logit(inc.outcome$ci.l)
inc.outcome$logit.ci.u<-logit(inc.outcome$ci.u)
inc.outcome$logit.se<-(inc.outcome$logit.ci.u-inc.outcome$logit.ci.l)/(2*1.96)

pool.rubin<- data.frame(matrix(NA, nrow = length(unique(data$studyname)), ncol = 8))
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
pool.rubin$studyname<-as.character(pool.rubin$studyname)
pool.rubin$studyname[pool.rubin$studyname=="1"]<-"Brazil_BahiaPaudaLima_Costa"
pool.rubin$studyname[pool.rubin$studyname=="2"]<-"Brazil_RiodeJaneiro_CunhaPrata"
pool.rubin$studyname[pool.rubin$studyname=="3"]<-"Brazil_RiodeJaneiro_Joao"
pool.rubin$studyname[pool.rubin$studyname=="4"]<-"Brazil_SP_RibeiraoPreto_Duarte"
pool.rubin$studyname[pool.rubin$studyname=="5"]<-"Colombia_Mulkey"
pool.rubin$studyname[pool.rubin$studyname=="6"]<-"FrenchGuiana_Pomar"
pool.rubin$studyname[pool.rubin$studyname=="7"]<-"Spain_Bardaji"
pool.rubin$studyname[pool.rubin$studyname=="8"]<-"Spain_Soriano"
pool.rubin$studyname[pool.rubin$studyname=="9"]<-"TrinidadTobago_Sohan"
pool.rubin$studyname[pool.rubin$studyname=="10"]<-"USA_Mulkey"

pool.rubin

abs.outcome<-pool.rubin
abs.outcome$incidence<-inv.logit(abs.outcome$logit.abs)*100
abs.outcome$ci.lb<-inv.logit(abs.outcome$logit.ci.lb)*100
abs.outcome$ci.ub<-inv.logit(abs.outcome$logit.ci.ub)*100
abs.outcome

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
#png(file="20210930 Fetal loss.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,13)), cex=1, steps=7, main="Fetal loss",
                xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, 
        sei = inv.logit.SE(pool.outcome$logit.se,pool.outcome$abs.risk), 
        rows=-0, cex=1)#Add pooled
addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

##################################################################################
###################################Fetal loss#####################################
##################################################################################

#Calculate incidence of outcome per study
for (i in 1:length(unique(data$.imp))) {
  d<-data[data$.imp==i,]
  inc<-as.data.frame(d %>% group_by(studyname) %>% summarise(p=f.incidence(czsn,dataset.n)[1],ci.l<-f.incidence(czsn,dataset.n)[2],ci.u<-f.incidence(czsn,dataset.n)[3]))
  colnames(inc)<-c("studyname","incidence","ci.l","ci.u")
  if(i==1) {inc.outcome<-inc}
  if(i>1) {inc.outcome<-rbind(inc.outcome,inc)}
}
inc.outcome
#Replace 0's by 0.000001
inc.outcome$incidence[inc.outcome$incidence==0]<-0.000001
inc.outcome$ci.l[inc.outcome$ci.l==0]<-0.000001

#Pool with Rubins rules
inc.outcome$logit.incidence<-logit(inc.outcome$incidence)
inc.outcome$logit.ci.l<-logit(inc.outcome$ci.l)
inc.outcome$logit.ci.u<-logit(inc.outcome$ci.u)
inc.outcome$logit.se<-(inc.outcome$logit.ci.u-inc.outcome$logit.ci.l)/(2*1.96)

pool.rubin<- data.frame(matrix(NA, nrow = length(unique(data$studyname)), ncol = 8))
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
pool.rubin$studyname<-as.character(pool.rubin$studyname)
pool.rubin$studyname[pool.rubin$studyname=="1"]<-"Brazil_BahiaPaudaLima_Costa"
pool.rubin$studyname[pool.rubin$studyname=="2"]<-"Brazil_RiodeJaneiro_CunhaPrata"
pool.rubin$studyname[pool.rubin$studyname=="3"]<-"Brazil_RiodeJaneiro_Joao"
pool.rubin$studyname[pool.rubin$studyname=="4"]<-"Brazil_SP_RibeiraoPreto_Duarte"
pool.rubin$studyname[pool.rubin$studyname=="5"]<-"Colombia_Mulkey"
pool.rubin$studyname[pool.rubin$studyname=="6"]<-"FrenchGuiana_Pomar"
pool.rubin$studyname[pool.rubin$studyname=="7"]<-"Spain_Bardaji"
pool.rubin$studyname[pool.rubin$studyname=="8"]<-"Spain_Soriano"
pool.rubin$studyname[pool.rubin$studyname=="9"]<-"TrinidadTobago_Sohan"
pool.rubin$studyname[pool.rubin$studyname=="10"]<-"USA_Mulkey"

pool.rubin

abs.outcome<-pool.rubin
abs.outcome$incidence<-inv.logit(abs.outcome$logit.abs)*100
abs.outcome$ci.lb<-inv.logit(abs.outcome$logit.ci.lb)*100
abs.outcome$ci.ub<-inv.logit(abs.outcome$logit.ci.ub)*100
abs.outcome

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
pool.outcome
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
addpoly(x = pool.outcome$abs.risk, 
        sei = inv.logit.SE(pool.outcome$logit.se,pool.outcome$abs.risk), 
        rows=-0, cex=1)#Add pooled
addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Forest plot - resized
#png(file="20210930 Congenital Zika Syndrome - resized.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,13)), cex=1, steps=8, main="Congenital Zika Syndrome",
                xlim=(c(-45,55)), alim=(c(0,35)))
addpoly(x = pool.outcome$abs.risk, 
        sei = inv.logit.SE(pool.outcome$logit.se,pool.outcome$abs.risk), 
        rows=-0, cex=1)#Add pooled
addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()