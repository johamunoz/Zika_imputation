
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
library(data.table) 


######################################################################################
#################################Load and prepare data################################
######################################################################################

setwd("/Users/jdamen/Documents/Julius/ZIKV analyses/2. Data")
data<-read.csv("20221027 zikv_imputed.csv",header=T)
setwd("/Users/jdamen/Documents/GitHub/Zika_imputation/2.1_Code 33 studies")
source("Functions Objective 1.R")
source("ZIKV prep v33.R")

#Systematic missings
data_sys <- read.xlsx2("/Users/jdamen/Documents/GitHub/Zika_imputation/1_Input_data/Table 1 outcomes.xlsx",sheetName="Systematic missings")
colnames(data_sys)<-c("variable","001-BRA","002-BRA","003-GUF","004-ESP","005-ESP","006-COL","007-COL",
                      "008-USA","009-GRD","010-BRA","011-BRA","012-TTO","013-BRA","014-BRA",
                      "015-BRA","016-HND","017-USA","018-COL","019-BRA","020-BRA","021-PRI",
                      "022-BRA","023-BRA","024-GTM","025-BRA","026-KEN","027-BRA")
data_sys<-melt(setDT(data_sys), id.vars = c("variable"), variable.name = "studyname",value.name="systematic")
data_sys<-data_sys[systematic==1]
data<-as.data.table(data)
for (i in 1:nrow(data_sys)){
  data[studyname_fac==data_sys[i,]$studyname,(data_sys[i,]$variable):=NA]
}
data<-as.data.frame(data)

data.zika.all<-data[data$zikv_preg==1,]
data.nozika.all<-data[data$zikv_preg==0,]


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

data.zika<-data.zika.all[!is.na(data.zika.all$microcephaly_bin_birth),]
data.nozika<-data.nozika.all[!is.na(data.nozika.all$microcephaly_bin_birth),]
data.1sma<-data[!is.na(data$microcephaly_bin_birth),]

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
fit1.coef<-as.data.frame(f.1ma.r.int(data.1sma,"microcephaly_bin_birth"))
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

data.zika<-data.zika.all[!is.na(data.zika.all$miscarriage),]
data.nozika<-data.nozika.all[!is.na(data.nozika.all$miscarriage),]
data.1sma<-data[!is.na(data$miscarriage),]

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
fit1.coef<-as.data.frame(f.1ma.r.int(data.1sma,"miscarriage"))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(data,fit1.coef)

#Store results
results<-cbind(results,pool.outcome.rr.1ma)
results$outcome<-"Miscarriage"
all.results<-rbind(all.results,results)

##################################################################################
###################################Fetal loss#####################################
##################################################################################

data.zika<-data.zika.all[!is.na(data.zika.all$loss),]
data.nozika<-data.nozika.all[!is.na(data.nozika.all$loss),]
data.1sma<-data[!is.na(data$loss),]

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
fit1.coef<-as.data.frame(f.1ma.r.int(data.1sma,"loss"))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(data,fit1.coef)

#Store results
results<-cbind(results,pool.outcome.rr.1ma)
results$outcome<-"Fetal loss"
all.results<-rbind(all.results,results)

##################################################################################
########################Congenital zika syndrome##################################
##################################################################################

data.zika<-data.zika.all[!is.na(data.zika.all$czs),]
data.nozika<-data.nozika.all[!is.na(data.nozika.all$czs),]
data.1sma<-data[!is.na(data$czs),]

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
#One stage meta-analysis with log link and random intercept per study, per imputed dataset
fit1.coef<-as.data.frame(f.1ma.r.int(data.1sma,"czs"))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(data.1sma,fit1.coef)

#Store results
results<-cbind(results,pool.outcome.rr.1ma)
results$outcome<-"CZS"
all.results<-rbind(all.results,results)


##################################################################################
################################Early fetal death#################################
##################################################################################

data.zika<-data.zika.all[!is.na(data.zika.all$efdeath),]
data.nozika<-data.nozika.all[!is.na(data.nozika.all$efdeath),]
data.1sma<-data[!is.na(data$efdeath),]

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
fit1.coef<-as.data.frame(f.1ma.r.int(data.1sma,"efdeath"))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(data,fit1.coef)

#Store results
results<-cbind(results,pool.outcome.rr.1ma)
results$outcome<-"Early fetal death"
all.results<-rbind(all.results,results)

##################################################################################
#################################Late fetal death#################################
##################################################################################

data.zika<-data.zika.all[!is.na(data.zika.all$lfdeath),]
data.nozika<-data.nozika.all[!is.na(data.nozika.all$lfdeath),]
data.1sma<-data[!is.na(data$lfdeath),]

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
fit1.coef<-as.data.frame(f.1ma.r.int(data.1sma,"lfdeath"))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(data,fit1.coef)

#Store results
results<-cbind(results,pool.outcome.rr.1ma)
results$outcome<-"Late fetal death"
all.results<-rbind(all.results,results)

##################################################################################
##############################Postnatal microcephaly##############################
##################################################################################

data.zika<-data.zika.all[!is.na(data.zika.all$microcephaly_bin_postnatal),]
data.nozika<-data.nozika.all[!is.na(data.nozika.all$microcephaly_bin_postnatal),]
data.1sma<-data[!is.na(data$microcephaly_bin_postnatal),]

#Zika-positive women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.zika,"microcephaly_bin_postnatal")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.zika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20221031 Postnatal microcephaly zika positive.png",width=750,height=750,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Postnatal microcephaly, zika-positive women",
                digits=2L,xlim=(c(-3,9)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-pool.outcome

#Zika-negative women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.nozika,"microcephaly_bin_postnatal")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.nozika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20221031 Postnatal microcephaly zika negative.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Postnatal microcephaly, zika-negative women",
                digits=2L,xlim=(c(-3,9)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome)

########Relative risk########

#calculate relative risk of outcome per study and per imputed dataset
rr.outcome.all<-f.rel.perstudy(data,"microcephaly_bin_postnatal")
#pool the relative risks per study per imputed dataset using rubins rules, resulting in RR per study
rr.outcome<-f.rel.poolrubin(data,rr.outcome.all)
#Pool the relative risks over the studies, resulting in one pooled summary relative risk over all studies
pool.outcome.rr<-f.rel.2s.ma(rr.outcome)

#Forest plot
#png(file="20221031 Postnatal microcephaly RR.png",width=750,height=500,res=100)
metafor::forest(rr.outcome$rr, ci.lb=rr.outcome$ci.lb, ci.ub=rr.outcome$ci.ub, 
                refline = 1, slab = rr.outcome$studyname,
                xlab = "Relative risk", pch = 19, psize=1,
                ylim=(c(-1,nrow(rr.outcome)+3)), cex=1, steps=7, main="Postnatal microcephaly",
                digits=2L,xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome.rr$rr, ci.lb=pool.outcome.rr$rr.ci.lb, ci.ub=pool.outcome.rr$rr.ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome.rr)

######One-stage meta-analysis relative risk
#One stage meta-analysis with log link and random intercept per study, per imputed dataset
fit1.coef<-as.data.frame(f.1ma.r.int(data.1sma,"microcephaly_bin_postnatal"))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(data,fit1.coef)

#Store results
results<-cbind(results,pool.outcome.rr.1ma)
results$outcome<-"Postnatal microcephaly"
all.results<-rbind(all.results,results)

##################################################################################
################################Neuroimaging abnormalities########################
##################################################################################

data.zika<-data.zika.all[!is.na(data.zika.all$neuroabnormality),]
data.nozika<-data.nozika.all[!is.na(data.nozika.all$neuroabnormality),]
data.1sma<-data[!is.na(data$neuroabnormality),]

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
fit1.coef<-as.data.frame(f.1ma.r.int(data.1sma,"neuroabnormality"))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(data,fit1.coef)

#Store results
results<-cbind(results,pool.outcome.rr.1ma)
results$outcome<-"Neuroimaging abnormalities"
all.results<-rbind(all.results,results)

##################################################################################
#################################Non-neurologic###################################
##################################################################################

data.zika<-data.zika.all[!is.na(data.zika.all$nonneurologic),]
data.nozika<-data.nozika.all[!is.na(data.nozika.all$nonneurologic),]
data.1sma<-data[!is.na(data$nonneurologic),]

#Zika-positive women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.zika,"nonneurologic")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.zika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20221031 Non-neurologic abnormalities zika positive.png",width=750,height=750,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Non-neurologic abnormalities, zika-positive women",
                digits=2L,xlim=(c(-3,9)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-pool.outcome

#Zika-negative women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.nozika,"nonneurologic")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.nozika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20221031 Non-neurologic abnormalities zika negative.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Non-neurologic abnormalities, zika-negative women",
                digits=2L,xlim=(c(-3,9)), alim=(c(0,6)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome)

########Relative risk########

#calculate relative risk of outcome per study and per imputed dataset
rr.outcome.all<-f.rel.perstudy(data,"nonneurologic")
#pool the relative risks per study per imputed dataset using rubins rules, resulting in RR per study
rr.outcome<-f.rel.poolrubin(data,rr.outcome.all)
#Pool the relative risks over the studies, resulting in one pooled summary relative risk over all studies
pool.outcome.rr<-f.rel.2s.ma(rr.outcome)

#Forest plot
#png(file="20221031 Non-neurologic abnormalities RR.png",width=750,height=500,res=100)
metafor::forest(rr.outcome$rr, ci.lb=rr.outcome$ci.lb, ci.ub=rr.outcome$ci.ub, 
                refline = 1, slab = rr.outcome$studyname,
                xlab = "Relative risk", pch = 19, psize=1,
                ylim=(c(-1,nrow(rr.outcome)+3)), cex=1, steps=7, main="Non-neurologic abnormalities",
                digits=2L,xlim=(c(-9,11)), alim=(c(0,6)))
addpoly(x = pool.outcome.rr$rr, ci.lb=pool.outcome.rr$rr.ci.lb, ci.ub=pool.outcome.rr$rr.ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome.rr)

######One-stage meta-analysis relative risk
#One stage meta-analysis with log link and random intercept per study, per imputed dataset
fit1.coef<-as.data.frame(f.1ma.r.int(data.1sma,"nonneurologic"))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(data,fit1.coef)

#Store results
results<-cbind(results,pool.outcome.rr.1ma)
results$outcome<-"Non-neurologic abnormalities"
all.results<-rbind(all.results,results)

######################
#write.xlsx2(all.results, "20221031 Results objective 1.xls", col.names = T, row.names = F)

