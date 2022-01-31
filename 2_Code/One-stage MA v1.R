
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


######################################################################################
#################################Load and prepare data################################
######################################################################################

setwd("/Users/jdamen/Documents/Julius/ZIKV analyses/2. Data")
load("imputation2.RData")
data <- complete(merged_imp, "long")
rm(merged_imp)
m<-max(data$.imp)
studynames<-c("014-BRA","001-BRA","002-BRA","010-BRA","007-COL",
              "003-GUF","005-ESP","004-ESP","012-TTO","008-USA")
data$studyname_fac<-as.factor(data$studyname,levels)

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
#data$lfdeath_micro<-0
#data$lfdeath_micro[data$lfdeath==1 & data$bdeath_ga>=28]<-1
#data$lfdeath<-as.factor(data$lfdeath)

data.zika<-data[data$zikv_preg==1,]
data.nozika<-data[data$zikv_preg==0,]

########################Analyses#############

fit1 <- glm(microcephaly_bin ~ zikv_preg, 
            data = data, 
            family = binomial(link="log"))
summary(fit1)
exp(fit1$coefficients[2])

fit2 <- glm(microcephaly_bin ~ 0 + SID + zikv_preg + AGE:SID + GENDER:SID + BILAT_0:SID, 
            data = ds.final, 
            family = binomial(link="log"))
summary(fit2)