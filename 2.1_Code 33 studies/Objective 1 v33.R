
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



######################################################################################
#################################Load and prepare data################################
######################################################################################

setwd("/Users/jdamen/Documents/Julius/ZIKV analyses/2. Data")
data<-read.csv("20221027 zikv_imputed.csv",header=T)
setwd("/Users/jdamen/Documents/GitHub/Zika_imputation/2.1_Code 33 studies")
source("Functions Objective 1.R")
source("ZIKV prep v33.R")

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
inc.outcome<-f.abs.perstudy(data.zika,"microcephaly_bin_birth")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.zika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20221031 Microcephaly zika positive.png",width=750,height=750,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Microcephaly, zika-positive women",
                digits=2L,xlim=(c(-15,85)), alim=(c(0,60)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-pool.outcome

#Zika-negative women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.nozika,"microcephaly_bin_birth")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.nozika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20221031 Microcephaly zika negative.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Microcephaly, zika-negative women",
                digits=2L,xlim=(c(-15,85)), alim=(c(0,60)))
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
rr.outcome.all<-f.rel.perstudy(data,"microcephaly_bin_birth")
#pool the relative risks per study per imputed dataset using rubins rules, resulting in RR per study
rr.outcome<-f.rel.poolrubin(data,rr.outcome.all)
#Pool the relative risks over the studies, resulting in one pooled summary relative risk over all studies
pool.outcome.rr<-f.rel.2s.ma(rr.outcome)

#Store results
results<-cbind(results,pool.outcome.rr)

######One-stage meta-analysis relative risk
#One stage meta-analysis with log link and random intercept per study, per imputed dataset
fit1.coef<-as.data.frame(f.1ma.r.int(data,"microcephaly_bin_birth"))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(data,fit1.coef)

#Forest plot
#png(file="20221031 Microcephaly RR.png",width=750,height=500,res=100)
metafor::forest(rr.outcome$rr, ci.lb=rr.outcome$ci.lb, ci.ub=rr.outcome$ci.ub, 
                refline = 1, slab = rr.outcome$studyname,
                xlab = "Relative risk", pch = 19, psize=1,
                ylim=(c(-1,nrow(rr.outcome)+3)), cex=1, steps=7, main="Microcephaly",
                digits=2L,xlim=(c(-9,11)), alim=(c(0,6)))

addpoly(x = pool.outcome.rr$rr, ci.lb=pool.outcome.rr$rr.ci.lb, ci.ub=pool.outcome.rr$rr.ci.ub,
        rows=0, cex=1)#Add pooled
addpoly(x = pool.outcome.rr.1ma$rr.1ma, ci.lb=pool.outcome.rr.1ma$rr.ci.lb.1ma, ci.ub=pool.outcome.rr.1ma$rr.ci.ub.1ma,
        rows=-1, cex=1)

#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

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
#png(file="20221031 Miscarriage zika positive.png",width=750,height=750,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Miscarriage, zika-positive women",
                digits=2L,xlim=(c(-3,9)), alim=(c(0,6)))
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
#png(file="20221031 Miscarriage zika negative.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Miscarriage, zika-negative women",
                digits=2L,xlim=(c(-3,9)), alim=(c(0,6)))
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
#png(file="20221031 Miscarriage RR.png",width=750,height=500,res=100)
metafor::forest(rr.outcome$rr, ci.lb=rr.outcome$ci.lb, ci.ub=rr.outcome$ci.ub, 
                refline = 1, slab = rr.outcome$studyname,
                xlab = "Relative risk", pch = 19, psize=1,
                ylim=(c(-1,nrow(rr.outcome)+3)), cex=1, steps=7, main="Miscarriage",
                digits=2L,xlim=(c(-9,11)), alim=(c(0,6)))
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
#png(file="20221031 Fetal loss zika positive.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=6, main="Fetal loss, zika-positive women",
                digits=2L,xlim=(c(-15,40)), alim=(c(0,25)))
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
#png(file="20221031 Fetal loss zika negative.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=6, main="Fetal loss, zika-negative women",
                digits=2L,xlim=(c(-15,40)), alim=(c(0,25)))
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
#png(file="20221031 Fetal loss RR.png",width=750,height=500,res=100)
metafor::forest(rr.outcome$rr, ci.lb=rr.outcome$ci.lb, ci.ub=rr.outcome$ci.ub, 
                refline = 1, slab = rr.outcome$studyname,
                xlab = "Relative risk", pch = 19, psize=1,
                ylim=(c(-1,nrow(rr.outcome)+3)), cex=1, steps=7, main="Fetal loss",
                digits=2L,xlim=(c(-9,11)), alim=(c(0,6)))
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
inc.outcome<-f.abs.perstudy(data.zika,"czs")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.zika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20221031 czs zika positive.png",width=750,height=750,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=10, main="Congenital zika syndrome, zika-positive women",
                digits=2L,xlim=(c(-15,70)), alim=(c(0,45)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-pool.outcome

#Zika-negative women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.nozika,"czs")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.nozika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20221031 czs zika negative.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=10, main="Congenital zika syndrome, zika-negative women",
                digits=2L,xlim=(c(-15,70)), alim=(c(0,45)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome)

########Relative risk########

#calculate relative risk of outcome per study and per imputed dataset
rr.outcome.all<-f.rel.perstudy(data,"czs")
#pool the relative risks per study per imputed dataset using rubins rules, resulting in RR per study
rr.outcome<-f.rel.poolrubin(data,rr.outcome.all)
#Pool the relative risks over the studies, resulting in one pooled summary relative risk over all studies
pool.outcome.rr<-f.rel.2s.ma(rr.outcome)

#Forest plot
#png(file="20221031 czs RR.png",width=750,height=500,res=100)
metafor::forest(rr.outcome$rr, ci.lb=rr.outcome$ci.lb, ci.ub=rr.outcome$ci.ub, 
                refline = 0, slab = rr.outcome$studyname,
                xlab = "Relative risk", pch = 19, psize=1,
                ylim=(c(-1,nrow(rr.outcome)+3)), cex=1, steps=7, main="Congenital zika syndrome",
                digits=2L,xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome.rr$rr, ci.lb=pool.outcome.rr$rr.ci.lb, ci.ub=pool.outcome.rr$rr.ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome.rr)

######One-stage meta-analysis relative risk
#First exclude 014-BRA and 023-BRA for this outcome as there is separation
#data.czs<-data[data$studyname_fac!="023-BRA",]
#data.czs<-data.czs[data.czs$studyname_fac!="014-BRA",]
#One stage meta-analysis with log link and random intercept per study, per imputed dataset
fit1.coef<-as.data.frame(f.1ma.r.int(data.czs,"czs"))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(data.czs,fit1.coef)

#Store results
results<-cbind(results,pool.outcome.rr.1ma)
results$outcome<-"CZS"
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
#png(file="20221031 Early fetal death zika positive.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Early fetal death, zika-positive women",
                digits=2L,xlim=(c(-9,11)), alim=(c(0,6)))
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
#png(file="20221031 Early fetal death zika negative.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Early fetal death, zika-negative women",
                digits=2L,xlim=(c(-9,11)), alim=(c(0,6)))
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
#png(file="20221031 Early fetal death RR.png",width=750,height=500,res=100)
metafor::forest(rr.outcome$rr, ci.lb=rr.outcome$ci.lb, ci.ub=rr.outcome$ci.ub, 
                refline = 0, slab = rr.outcome$studyname,
                xlab = "Relative risk", pch = 19, psize=1,
                ylim=(c(-1,nrow(rr.outcome)+3)), cex=1, steps=7, main="Early fetal death",
                digits=2L,xlim=(c(-9,11)), alim=(c(0,6)))
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
#png(file="20221031 Late fetal death zika positive.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Late fetal death, zika-positive women",
                digits=2L,xlim=(c(-9,11)), alim=(c(0,6)))
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
#png(file="20221031 Late fetal death zika negative.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Late fetal death, zika-negative women",
                digits=2L,xlim=(c(-9,11)), alim=(c(0,6)))
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
#png(file="20221031 Late fetal death RR.png",width=750,height=500,res=100)
metafor::forest(rr.outcome$rr, ci.lb=rr.outcome$ci.lb, ci.ub=rr.outcome$ci.ub, 
                refline = 0, slab = rr.outcome$studyname,
                xlab = "Relative risk", pch = 19, psize=1,
                ylim=(c(-1,nrow(rr.outcome)+3)), cex=1, steps=7, main="Late fetal death",
                digits=2L,xlim=(c(-9,11)), alim=(c(0,6)))
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
##############################Postnatal microcephaly##############################
##################################################################################

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
#png(file="20221031 Neuroabnormality zika positive.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Neuroimaging abnormalities, zika-positive women",
                digits=2L,xlim=(c(-9,11)), alim=(c(0,6)))
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
#png(file="20221031 Neuroabnormality zika negative.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Neuroimaging abnormalities, zika-negative women",
                digits=2L,xlim=(c(-9,11)), alim=(c(0,6)))
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
#png(file="20221031 Neuroabnormality RR.png",width=750,height=500,res=100)
metafor::forest(rr.outcome$rr, ci.lb=rr.outcome$ci.lb, ci.ub=rr.outcome$ci.ub, 
                refline = 1, slab = rr.outcome$studyname,
                xlab = "Relative risk", pch = 19, psize=1,
                ylim=(c(-1,nrow(rr.outcome)+3)), cex=1, steps=7, main="Neuroimaging abnormalities",
                digits=2L,xlim=(c(-9,11)), alim=(c(0,6)))
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

##################################################################################
#################################Non-neurologic###################################
##################################################################################


######################
#write.xlsx2(all.results, "20221031 Results objective 1.xls", col.names = T, row.names = F)
