
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


#

######################################################################################
#################################Load and prepare data################################
######################################################################################

setwd("/Users/jdamen/Documents/Julius/ZIKV analyses/2. Data")
load("imputation2.RData")
setwd("/Users/jdamen/Documents/GitHub/Zika_imputation/2_Code")
source("ZIKV prep.R")

########################Analyses#############
#Absolute risk with log link:
fit1 <- glm(microcephaly_bin ~ zikv_preg, 
            data = data, family = binomial(link="log")) #log RR
fit1.coef<-summary(fit1)$coefficients[2,]
c(exp(fit1.coef[1]),exp(fit1.coef[1]-1.96*fit1.coef[2]),exp(fit1.coef[1]+1.96*fit1.coef[2])) #1.060795 -> is dit de RR?

#Logit link
fit1b <- glm(microcephaly_bin ~ zikv_preg, 
            data = data, family = binomial(link="logit"))
fit1b.coef<-summary(fit1b)$coefficients[2,]
c(exp(fit1b.coef[1]),exp(fit1b.coef[1]-1.96*fit1b.coef[2]),exp(fit1b.coef[1]+1.96*fit1b.coef[2])) #0.51584 -> is dit de OR?
#confint(fit1b)

#Random intercept per study
#Log
fit2 <- glmer(microcephaly_bin ~ zikv_preg + (1 | studyname_fac), 
              data=data, family = binomial(link = "log"))
fit2.coef<-summary(fit2)$coefficients[2,]
c(exp(fit2.coef[1]),exp(fit2.coef[1]-1.96*fit2.coef[2]),exp(fit2.coef[1]+1.96*fit2.coef[2])) #1.283628 -> RR met random intercept per studie?
#Logit
fit2b <- glmer(microcephaly_bin ~ zikv_preg + (1 | studyname_fac), 
              data=data, family = binomial(link = "logit"))
fit2b.coef<-summary(fit2b)$coefficients[2,]
c(exp(fit2b.coef[1]),exp(fit2b.coef[1]-1.96*fit2b.coef[2]),exp(fit2b.coef[1]+1.96*fit2b.coef[2])) #0.57106 -> OR met random intercept per studie?

#Random intercept and random slope
#Log
fit3<-glmer(microcephaly_bin ~ zikv_preg + (1+zikv_preg | studyname_fac), 
      data=data, family = binomial(link = "log")) #Geen warnings!?
fit3.coef<-summary(fit3)$coefficients[2,]
c(exp(fit3.coef[1]),exp(fit3.coef[1]-1.96*fit3.coef[2]),exp(fit3.coef[1]+1.96*fit3.coef[2])) #0.8506 -> RR met random intercept en random slope
#Logit
fit3b<-glmer(microcephaly_bin ~ zikv_preg + (1+zikv_preg | studyname_fac), 
            data=data, family = binomial(link = "logit")) 
fit3b.coef<-summary(fit3b)$coefficients[2,]
c(exp(fit3b.coef[1]),exp(fit3b.coef[1]-1.96*fit3b.coef[2]),exp(fit3b.coef[1]+1.96*fit3b.coef[2])) #0.8529 -> OR met random intercept en random slope


#Absolute risks
newdata<-subset(data,select=c(studyname_fac,zikv_preg))
predicted<-newdata
predicted$fit1.pred<-predict(fit1, newdata = newdata, type = "response")
unique(predicted$fit1.pred) #Geen transformatie nodig?
predicted$fit1b.pred<-predict(fit1b, newdata = newdata, type = "response")
unique(predicted$fit1b.pred) #Exact hetzelfde
predicted$fit2.pred<-predict(fit2, newdata = newdata, type = "response")
unique(predicted$fit2.pred)
predicted$fit2b.pred<-predict(fit2b, newdata = newdata, type = "response")
unique(predicted$fit2b.pred) #Vanaf 3e decimaal achter de komma nét iets anders dan fit2
predicted$fit3.pred<-predict(fit3, newdata = newdata, type = "response")
unique(predicted$fit3.pred) #Vrij weinig verschil met fit2
predicted$fit3b.pred<-predict(fit3b, newdata = newdata, type = "response")
unique(predicted$fit3b.pred) #Ook weinig verschil



##Objective 1: one-stage meta-analysis per outcome with log link and random intercept per study -> fit 1
#Microcephaly
fit1 <- glmer(microcephaly_bin ~ zikv_preg + (1 | studyname_fac), 
              data=data, family = binomial(link = "log"))
fit1.coef<-summary(fit1)$coefficients[2,]
c(exp(fit1.coef[1]),exp(fit1.coef[1]-1.96*fit1.coef[2]),exp(fit1.coef[1]+1.96*fit1.coef[2])) 
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
fit2c<-glmer(microcephaly_bin ~ zikv_preg + age + maritalstat + (1+zikv_preg | studyname_fac), 
            data=data, family = binomial(link = "log"))
fit2c.coef<-as.data.frame(summary(fit2c)$coefficients[c(2:nrow(summary(fit2c)$coefficients)),c(1:2)])
colnames(fit2c.coef)<-c("log.rr", "log.rr.se")
fit2c.coef$rr<-exp(fit2c.coef$log.rr)
fit2c.coef$rr.lb<-exp(fit2c.coef$log.rr-1.96*fit2c.coef$log.rr.se)
fit2c.coef$rr.ub<-exp(fit2c.coef$log.rr+1.96*fit2c.coef$log.rr.se)
fit2c.coef


##Fit 2 + effect modifiers (one by one)

##Calculate absolute risks for fit 2: check protocol!

