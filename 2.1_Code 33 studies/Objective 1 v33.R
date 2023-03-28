
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

setwd("/Users/jdamen/Library/CloudStorage/OneDrive-UMCUtrecht/Research/WHO ZIKA/2. Data")
data<-read.csv("20230227 zikv_imputed.csv",header=T)
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

#Exclude studies that are not selected for main analyses
data<-data[data$studyname_fac!="002-BRA" & data$studyname_fac!="008-USA" &  
              data$studyname_fac!="011-BRA" & data$studyname_fac!="013-BRA" & data$studyname_fac!="018-COL",]
data$studyname_fac<-droplevels(data$studyname_fac)

data.zika.all<-data[data$zikv_preg==1,]
data.nozika.all<-data[data$zikv_preg==0,]


########################Analyses#############
#Change directory to save plots
setwd("/Users/jdamen/Library/CloudStorage/OneDrive-UMCUtrecht/Research/WHO ZIKA/4. Resultaten")

#Data frame to store results
all.results<-data.frame(Outcome=character(),
                        File=numeric(), 
                        User=numeric(), 
                        stringsAsFactors=FALSE)

##################################################################################
###################################Microcephaly###################################
##################################################################################

data.zika<-data.zika.all[!is.na(data.zika.all$microcephaly_bin_birth) & (data.zika.all$end_ga_compl>=24 & !is.na(data.zika.all$end_ga_compl)),]
data.nozika<-data.nozika.all[!is.na(data.nozika.all$microcephaly_bin_birth) & (data.nozika.all$end_ga_compl>=24 & !is.na(data.nozika.all$end_ga_compl)),]

#Zika-positive women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.zika,"microcephaly_bin_birth")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.zika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20230328 Microcephaly birth zika positive.png",width=800,height=750,res=100)
#png(file="20230328 Microcephaly birth zika positive SENSITIVITY.png",width=750,height=750,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=6, main="Microcephaly at birth, zika-positive women",
                digits=2L,xlim=(c(-25,80)), alim=(c(0,60)))
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
#png(file="20230328 Microcephaly birth zika negative.png",width=750,height=550,res=100)
#png(file="20230328 Microcephaly birth zika negative SENSITIVITY.png",width=750,height=550,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=6, main="Microcephaly at birth, zika-negative women",
                digits=2L,xlim=(c(-25,80)), alim=(c(0,60)))
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
data.1sma<-data[!is.na(data$microcephaly_bin_birth) & (data$end_ga_compl>=24 & !is.na(data$end_ga_compl)),]
#One stage meta-analysis with log link and random intercept per study, per imputed dataset
fit1.coef<-as.data.frame(f.1ma.r.int(data.1sma,"microcephaly_bin_birth",poisson))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(fit1.coef)

#Forest plot
#png(file="20230328 Microcephaly birth RR.png",width=750,height=600,res=100)
#png(file="20230328 Microcephaly birth RR SENSITIVITY.png",width=750,height=600,res=100)
metafor::forest(rr.outcome$rr, ci.lb=rr.outcome$ci.lb, ci.ub=rr.outcome$ci.ub, 
                refline = 1, slab = rr.outcome$studyname,
                xlab = "Relative risk", pch = 19, psize=1,
                ylim=(c(-1,nrow(rr.outcome)+3)), cex=1, steps=6, main="Microcephaly at birth",
                digits=2L,xlim=(c(-5,10)), alim=(c(0,6)))
addpoly(x = pool.outcome.rr$rr, ci.lb=pool.outcome.rr$rr.ci.lb, ci.ub=pool.outcome.rr$rr.ci.ub,
        rows=0, cex=1)#Add pooled
addpoly(x = pool.outcome.rr.1ma$rr.1ma, ci.lb=pool.outcome.rr.1ma$rr.ci.lb.1ma, ci.ub=pool.outcome.rr.1ma$rr.ci.ub.1ma,
        rows=-1, cex=1)
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome.rr.1ma)
results$outcome<-"Microcephaly at birth"
all.results<-results

##################################################################################
###################################Miscarriage###################################
##################################################################################

data.zika<-data.zika.all[!is.na(data.zika.all$miscarriage),]
data.nozika<-data.nozika.all[!is.na(data.nozika.all$miscarriage),]

#Zika-positive women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.zika,"miscarriage")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.zika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20230328 Miscarriage zika positive.png",width=750,height=750,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=8, main="Miscarriage, zika-positive women",
                digits=2L,xlim=(c(-3,10)), alim=(c(0,7)))
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
#png(file="20230328 Miscarriage zika negative.png",width=750,height=550,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=8, main="Miscarriage, zika-negative women",
                digits=2L,xlim=(c(-3,10)), alim=(c(0,7)))
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

#Store results
results<-cbind(results,pool.outcome.rr)

######One-stage meta-analysis relative risk
data.1sma<-data[!is.na(data$miscarriage),]
#One stage meta-analysis with log link and random intercept per study, per imputed dataset
fit1.coef<-as.data.frame(f.1ma.r.int(data.1sma,"miscarriage",poisson))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(fit1.coef)

#Forest plot
#png(file="20230328 Miscarriage RR.png",width=750,height=400,res=100)
metafor::forest(rr.outcome$rr, ci.lb=rr.outcome$ci.lb, ci.ub=rr.outcome$ci.ub, 
                refline = 1, slab = rr.outcome$studyname,
                xlab = "Relative risk", pch = 19, psize=1,
                ylim=(c(-1,nrow(rr.outcome)+3)), cex=1, steps=6, main="Miscarriage",
                digits=2L,xlim=(c(-5,10)), alim=(c(0,5)))
addpoly(x = pool.outcome.rr$rr, ci.lb=pool.outcome.rr$rr.ci.lb, ci.ub=pool.outcome.rr$rr.ci.ub,
        rows=0, cex=1)#Add pooled
addpoly(x = pool.outcome.rr.1ma$rr.1ma, ci.lb=pool.outcome.rr.1ma$rr.ci.lb.1ma, ci.ub=pool.outcome.rr.1ma$rr.ci.ub.1ma,
        rows=-1, cex=1)
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome.rr.1ma)
results$outcome<-"Miscarriage"
all.results<-rbind(all.results,results)

##################################################################################
###################################Fetal loss#####################################
##################################################################################

data.zika<-data.zika.all[!is.na(data.zika.all$loss) & (data.zika.all$end_ga_compl>=20 & !is.na(data.zika.all$end_ga_compl)),]
data.nozika<-data.nozika.all[!is.na(data.nozika.all$loss) & (data.nozika.all$end_ga_compl>=20 & !is.na(data.nozika.all$end_ga_compl)),]

#Zika-positive women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.zika,"loss")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.zika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20230328 Fetal loss zika positive.png",width=750,height=650,res=100)
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
#png(file="20230328 Fetal loss zika negative.png",width=750,height=550,res=100)
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

#Store results
results<-cbind(results,pool.outcome.rr)

######One-stage meta-analysis relative risk
data.1sma<-data[!is.na(data$loss) & (data$end_ga_compl>=20 & !is.na(data$end_ga_compl)),]
#One stage meta-analysis with log link and random intercept per study, per imputed dataset
fit1.coef<-as.data.frame(f.1ma.r.int(data.1sma,"loss",poisson))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(fit1.coef)

#Forest plot
#png(file="20230328 Fetal loss RR.png",width=750,height=600,res=100)
metafor::forest(rr.outcome$rr, ci.lb=rr.outcome$ci.lb, ci.ub=rr.outcome$ci.ub, 
                refline = 1, slab = rr.outcome$studyname,
                xlab = "Relative risk", pch = 19, psize=1,
                ylim=(c(-1,nrow(rr.outcome)+3)), cex=1, steps=6, main="Fetal loss",
                digits=2L,xlim=(c(-5,10)), alim=(c(0,5)))
addpoly(x = pool.outcome.rr$rr, ci.lb=pool.outcome.rr$rr.ci.lb, ci.ub=pool.outcome.rr$rr.ci.ub,
        rows=0, cex=1)#Add pooled
addpoly(x = pool.outcome.rr.1ma$rr.1ma, ci.lb=pool.outcome.rr.1ma$rr.ci.lb.1ma, ci.ub=pool.outcome.rr.1ma$rr.ci.ub.1ma,
        rows=-1, cex=1)
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome.rr.1ma)
results$outcome<-"Fetal loss"
all.results<-rbind(all.results,results)

##################################################################################
########################Congenital zika syndrome##################################
##################################################################################

data.zika<-data.zika.all[!is.na(data.zika.all$ch_czs) & (data.zika.all$end_ga_compl>=24 & !is.na(data.zika.all$end_ga_compl)),]
data.nozika<-data.nozika.all[!is.na(data.nozika.all$ch_czs) & (data.nozika.all$end_ga_compl>=24 & !is.na(data.nozika.all$end_ga_compl)),]

#Zika-positive women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.zika,"ch_czs")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.zika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20230328 czs study zika positive.png",width=750,height=650,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=10, main="Congenital zika syndrome - study definition, zika-positive women",
                digits=2L,xlim=(c(-15,70)), alim=(c(0,45)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-pool.outcome

#Zika-negative women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.nozika,"ch_czs")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.nozika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20230328 czs study zika negative.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=10, main="Congenital zika syndrome - study definition, zika-negative women",
                digits=2L,xlim=(c(-15,70)), alim=(c(0,45)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome)

########Relative risk########

#calculate relative risk of outcome per study and per imputed dataset
rr.outcome.all<-f.rel.perstudy(data,"ch_czs")
#pool the relative risks per study per imputed dataset using rubins rules, resulting in RR per study
rr.outcome<-f.rel.poolrubin(data,rr.outcome.all)
#Pool the relative risks over the studies, resulting in one pooled summary relative risk over all studies
pool.outcome.rr<-f.rel.2s.ma(rr.outcome)

#Store results
results<-cbind(results,pool.outcome.rr)

######One-stage meta-analysis relative risk
data.1sma<-data[!is.na(data$ch_czs) & (data$end_ga_compl>=24 & !is.na(data$end_ga_compl)),]
#One stage meta-analysis with log link and random intercept per study, per imputed dataset
fit1.coef<-as.data.frame(f.1ma.r.int(data.1sma,"ch_czs",poisson))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(fit1.coef)

#Forest plot
#png(file="20230328 czs study RR.png",width=750,height=500,res=100)
metafor::forest(rr.outcome$rr, ci.lb=rr.outcome$ci.lb, ci.ub=rr.outcome$ci.ub, 
                refline = 1, slab = rr.outcome$studyname,
                xlab = "Relative risk", pch = 19, psize=1,
                ylim=(c(-1,nrow(rr.outcome)+3)), cex=1, steps=6, main="Congenital zika syndrome - study definition",
                digits=2L,xlim=(c(-5,10)), alim=(c(0,5)))
addpoly(x = pool.outcome.rr$rr, ci.lb=pool.outcome.rr$rr.ci.lb, ci.ub=pool.outcome.rr$rr.ci.ub,
        rows=0, cex=1)#Add pooled
addpoly(x = pool.outcome.rr.1ma$rr.1ma, ci.lb=pool.outcome.rr.1ma$rr.ci.lb.1ma, ci.ub=pool.outcome.rr.1ma$rr.ci.ub.1ma,
        rows=-1, cex=1)
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome.rr.1ma)
results$outcome<-"CZS study definition"
all.results<-rbind(all.results,results)

#Export results
#write.xlsx2(all.results, "20230328 Results main outcomes objective 1.xlsx", col.names = T, row.names = F)

#Prepare table for article
table.obj1<-all.results
table.obj1$abs.pos.95<-paste0(sprintf("%.2f",all.results$abs.risk),"% (95% CI ",sprintf("%.2f",all.results$ci.lb),"-",sprintf("%.2f",all.results$ci.ub),")")
table.obj1$abs.neg.95<-paste0(sprintf("%.2f",all.results$abs.risk.neg),"% (95% CI ",sprintf("%.2f",all.results$ci.lb.neg),"-",sprintf("%.2f",all.results$ci.ub.neg),")")
table.obj1$rr.2ma.95<-paste0(sprintf("%.2f",all.results$rr)," (95% CI ",sprintf("%.2f",all.results$rr.ci.lb),"-",sprintf("%.2f",all.results$rr.ci.ub),")")
table.obj1$rr.2ma.95.n<-paste0(sprintf("%.2f",all.results$rr)," (95% CI ",sprintf("%.2f",all.results$rr.ci.lb),"-",sprintf("%.2f",all.results$rr.ci.ub),"); n=",all.results$n.studies)
table.obj1$rr.1ma.95<-paste0(sprintf("%.2f",all.results$rr.1ma)," (95% CI ",sprintf("%.2f",all.results$rr.ci.lb.1ma),"-",sprintf("%.2f",all.results$rr.ci.ub.1ma),")")
table.obj1<-subset(table.obj1,select=c(outcome,abs.pos.95,abs.neg.95,rr.2ma.95.n,rr.1ma.95))

#write.xlsx2(table.obj1, "20230328 Table results main outcomes objective 1.xlsx", col.names = T, row.names = F)

##################################################################################
################################Early fetal death#################################
##################################################################################

data.zika<-data.zika.all[!is.na(data.zika.all$efdeath),]
data.nozika<-data.nozika.all[!is.na(data.nozika.all$efdeath),]

#Zika-positive women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.zika,"efdeath")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.zika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20230328 Early fetal death zika positive.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=6, main="Early fetal death, zika-positive women",
                digits=2L,xlim=(c(-5,10)), alim=(c(0,5)))
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
#png(file="20230328 Early fetal death zika negative.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=6, main="Early fetal death, zika-negative women",
                digits=2L,xlim=(c(-5,10)), alim=(c(0,5)))
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

#Store results
results<-cbind(results,pool.outcome.rr)

######One-stage meta-analysis relative risk
data.1sma<-data[!is.na(data$efdeath) & (data$end_ga_compl>=20 & !is.na(data$end_ga_compl)),]
#One stage meta-analysis with log link and random intercept per study, per imputed dataset
fit1.coef<-as.data.frame(f.1ma.r.int(data.1sma,"efdeath",poisson))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(fit1.coef)

#Forest plot
#png(file="20230328 Early fetal death RR.png",width=750,height=500,res=100)
metafor::forest(rr.outcome$rr, ci.lb=rr.outcome$ci.lb, ci.ub=rr.outcome$ci.ub, 
                refline = 1, slab = rr.outcome$studyname,
                xlab = "Relative risk", pch = 19, psize=1,
                ylim=(c(-1,nrow(rr.outcome)+3)), cex=1, steps=6, main="Early fetal death",
                digits=2L,xlim=(c(-5,10)), alim=(c(0,5)))
#addpoly(x = pool.outcome.rr$rr, ci.lb=pool.outcome.rr$rr.ci.lb, ci.ub=pool.outcome.rr$rr.ci.ub,
#        rows=0, cex=1)#Add pooled
addpoly(x = pool.outcome.rr.1ma$rr.1ma, ci.lb=pool.outcome.rr.1ma$rr.ci.lb.1ma, ci.ub=pool.outcome.rr.1ma$rr.ci.ub.1ma,
        rows=-1, cex=1)
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome.rr.1ma)
results$outcome<-"Early fetal death"
all.results<-rbind(all.results,results)

##################################################################################
#################################Late fetal death#################################
##################################################################################

data.zika<-data.zika.all[!is.na(data.zika.all$lfdeath),]
data.nozika<-data.nozika.all[!is.na(data.nozika.all$lfdeath),]

#Zika-positive women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.zika,"lfdeath")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.zika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20230328 Late fetal death zika positive.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=6, main="Late fetal death, zika-positive women",
                digits=2L,xlim=(c(-5,10)), alim=(c(0,5)))
#addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
#        rows=0, cex=1)#Add pooled
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
#png(file="20230328 Late fetal death zika negative.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=6, main="Late fetal death, zika-negative women",
                digits=2L,xlim=(c(-5,10)), alim=(c(0,5)))
#addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
#        rows=0, cex=1)#Add pooled
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

#Store results
results<-cbind(results,pool.outcome.rr)

######One-stage meta-analysis relative risk
data.1sma<-data[!is.na(data$lfdeath) & (data$end_ga_compl>=28 & !is.na(data$end_ga_compl)),]
#One stage meta-analysis with log link and random intercept per study, per imputed dataset
fit1.coef<-as.data.frame(f.1ma.r.int(data.1sma,"lfdeath",poisson))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(fit1.coef)

#Forest plot
#png(file="20230328 Late fetal death RR.png",width=750,height=500,res=100)
metafor::forest(rr.outcome$rr, ci.lb=rr.outcome$ci.lb, ci.ub=rr.outcome$ci.ub, 
                refline = 1, slab = rr.outcome$studyname,
                xlab = "Relative risk", pch = 19, psize=1,
                ylim=(c(-1,nrow(rr.outcome)+3)), cex=1, steps=6, main="Late fetal death",
                digits=2L,xlim=(c(-5,10)), alim=(c(0,5)))
#addpoly(x = pool.outcome.rr$rr, ci.lb=pool.outcome.rr$rr.ci.lb, ci.ub=pool.outcome.rr$rr.ci.ub,
#        rows=0, cex=1)#Add pooled
addpoly(x = pool.outcome.rr.1ma$rr.1ma, ci.lb=pool.outcome.rr.1ma$rr.ci.lb.1ma, ci.ub=pool.outcome.rr.1ma$rr.ci.ub.1ma,
        rows=-1, cex=1)
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome.rr.1ma)
results$outcome<-"Late fetal death"
all.results<-rbind(all.results,results)

##################################################################################
##############################Postnatal microcephaly##############################
##################################################################################

data.zika<-data.zika.all[!is.na(data.zika.all$microcephaly_bin_postnatal),]
data.nozika<-data.nozika.all[!is.na(data.nozika.all$microcephaly_bin_postnatal),]

#Zika-positive women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.zika,"microcephaly_bin_postnatal")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.zika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20230328 Postnatal microcephaly zika positive.png",width=750,height=750,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Postnatal microcephaly, zika-positive women",
                digits=2L,xlim=(c(-3,9)), alim=(c(0,6)))
#addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
#        rows=0, cex=1)#Add pooled
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
#png(file="20230328 Postnatal microcephaly zika negative.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Postnatal microcephaly, zika-negative women",
                digits=2L,xlim=(c(-3,9)), alim=(c(0,6)))
#addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
#        rows=0, cex=1)#Add pooled
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

#Store results
results<-cbind(results,pool.outcome.rr)

######One-stage meta-analysis relative risk
data.1sma<-data[!is.na(data$microcephaly_bin_postnatal),]
#One stage meta-analysis with log link and random intercept per study, per imputed dataset
fit1.coef<-as.data.frame(f.1ma.r.int(data.1sma,"microcephaly_bin_postnatal",poisson))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(fit1.coef)

#Forest plot
#png(file="20230328 Postnatal microcephaly RR.png",width=750,height=500,res=100)
metafor::forest(rr.outcome$rr, ci.lb=rr.outcome$ci.lb, ci.ub=rr.outcome$ci.ub, 
                refline = 1, slab = rr.outcome$studyname,
                xlab = "Relative risk", pch = 19, psize=1,
                ylim=(c(-1,nrow(rr.outcome)+3)), cex=1, steps=6, main="Postnatal microcephaly",
                digits=2L,xlim=(c(-5,10)), alim=(c(0,5)))
#addpoly(x = pool.outcome.rr$rr, ci.lb=pool.outcome.rr$rr.ci.lb, ci.ub=pool.outcome.rr$rr.ci.ub,
#        rows=0, cex=1)#Add pooled
addpoly(x = pool.outcome.rr.1ma$rr.1ma, ci.lb=pool.outcome.rr.1ma$rr.ci.lb.1ma, ci.ub=pool.outcome.rr.1ma$rr.ci.ub.1ma,
        rows=-1, cex=1)
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome.rr.1ma)
results$outcome<-"Postnatal microcephaly"
all.results<-rbind(all.results,results)

##################################################################################
#############################Fetal microcephaly###################################
##################################################################################

data.zika<-data.zika.all[!is.na(data.zika.all$microcephaly_bin_fet) & (data.zika.all$bdeath==0 & !is.na(data.zika.all$bdeath)),]
data.nozika<-data.nozika.all[!is.na(data.nozika.all$microcephaly_bin_fet) & (data.nozika.all$bdeath==0 & !is.na(data.nozika.all$bdeath)),]

#Zika-positive women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.zika,"microcephaly_bin_fet")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.zika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20230328 Fetal microcephaly zika positive.png",width=750,height=750,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=6, main="Fetal microcephaly, zika-positive women",
                digits=2L,xlim=(c(-15,130)), alim=(c(0,100)))
addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
        rows=0, cex=1)#Add pooled
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-pool.outcome

#Zika-negative women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.nozika,"microcephaly_bin_fet")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.nozika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20230328 Fetal microcephaly zika negative.png",width=750,height=550,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=6, main="Fetal microcephaly, zika-negative women",
                digits=2L,xlim=(c(-15,130)), alim=(c(0,100)))
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
rr.outcome.all<-f.rel.perstudy(data,"microcephaly_bin_fet")
#pool the relative risks per study per imputed dataset using rubins rules, resulting in RR per study
rr.outcome<-f.rel.poolrubin(data,rr.outcome.all)
#Pool the relative risks over the studies, resulting in one pooled summary relative risk over all studies
pool.outcome.rr<-f.rel.2s.ma(rr.outcome)

#Store results
results<-cbind(results,pool.outcome.rr)

######One-stage meta-analysis relative risk
data.1sma<-data[!is.na(data$microcephaly_bin_fet) & (data$end_ga_compl>=24 & !is.na(data$end_ga_compl)),]
#One stage meta-analysis with log link and random intercept per study, per imputed dataset
fit1.coef<-as.data.frame(f.1ma.r.int(data.1sma,"microcephaly_bin_fet",poisson))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(fit1.coef)

#Forest plot
#png(file="20230328 Fetal microcephaly RR.png",width=750,height=600,res=100)
metafor::forest(rr.outcome$rr, ci.lb=rr.outcome$ci.lb, ci.ub=rr.outcome$ci.ub, 
                refline = 1, slab = rr.outcome$studyname,
                xlab = "Relative risk", pch = 19, psize=1,
                ylim=(c(-1,nrow(rr.outcome)+3)), cex=1, steps=6, main="Fetal microcephaly",
                digits=2L,xlim=(c(-5,10)), alim=(c(0,5)))
addpoly(x = pool.outcome.rr$rr, ci.lb=pool.outcome.rr$rr.ci.lb, ci.ub=pool.outcome.rr$rr.ci.ub,
        rows=0, cex=1)#Add pooled
addpoly(x = pool.outcome.rr.1ma$rr.1ma, ci.lb=pool.outcome.rr.1ma$rr.ci.lb.1ma, ci.ub=pool.outcome.rr.1ma$rr.ci.ub.1ma,
        rows=-1, cex=1)
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome.rr.1ma)
results$outcome<-"Fetal microcephaly"
all.results<-results


##################################################################################
################################Neuroimaging abnormalities########################
##################################################################################

data.zika<-data.zika.all[!is.na(data.zika.all$neuroabnormality),]
data.nozika<-data.nozika.all[!is.na(data.nozika.all$neuroabnormality),]

#Zika-positive women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.zika,"neuroabnormality")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.zika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20230328 Neuroabnormality zika positive.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=6, main="Neuroimaging abnormalities, zika-positive women",
                digits=2L,xlim=(c(-5,10)), alim=(c(0,5)))
#addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
#        rows=0, cex=1)#Add pooled
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
#png(file="20230328 Neuroabnormality zika negative.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=6, main="Neuroimaging abnormalities, zika-negative women",
                digits=2L,xlim=(c(-5,10)), alim=(c(0,5)))
#addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
#        rows=0, cex=1)#Add pooled
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

#Store results
results<-cbind(results,pool.outcome.rr)

######One-stage meta-analysis relative risk
data.1sma<-data[!is.na(data$neuroabnormality) & (data$end_ga_compl>=24 & !is.na(data$end_ga_compl)),]
#One stage meta-analysis with log link and random intercept per study, per imputed dataset
fit1.coef<-as.data.frame(f.1ma.r.int(data.1sma,"neuroabnormality",poisson))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(fit1.coef)

#Forest plot
#png(file="20230328 Neuroabnormality RR.png",width=750,height=500,res=100)
metafor::forest(rr.outcome$rr, ci.lb=rr.outcome$ci.lb, ci.ub=rr.outcome$ci.ub, 
                refline = 1, slab = rr.outcome$studyname,
                xlab = "Relative risk", pch = 19, psize=1,
                ylim=(c(-1,nrow(rr.outcome)+3)), cex=1, steps=6, main="Neuroimaging abnormalities",
                digits=2L,xlim=(c(-5,10)), alim=(c(0,5)))
addpoly(x = pool.outcome.rr$rr, ci.lb=pool.outcome.rr$rr.ci.lb, ci.ub=pool.outcome.rr$rr.ci.ub,
        rows=0, cex=1)#Add pooled
addpoly(x = pool.outcome.rr.1ma$rr.1ma, ci.lb=pool.outcome.rr.1ma$rr.ci.lb.1ma, ci.ub=pool.outcome.rr.1ma$rr.ci.ub.1ma,
        rows=-1, cex=1)
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome.rr.1ma)
results$outcome<-"Neuroimaging abnormalities"
all.results<-rbind(all.results,results)

##################################################################################
#################################Non-neurologic###################################
##################################################################################

data.zika<-data.zika.all[!is.na(data.zika.all$nonneurologic),]
data.nozika<-data.nozika.all[!is.na(data.nozika.all$nonneurologic),]

#Zika-positive women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.zika,"nonneurologic")
#Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
abs.outcome<-f.abs.poolrubin(data.zika,inc.outcome)
#Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
pool.outcome<-f.abs.2s.ma(abs.outcome)

#Forest plot
#png(file="20230328 Non-neurologic abnormalities zika positive.png",width=750,height=750,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Non-neurologic abnormalities, zika-positive women",
                digits=2L,xlim=(c(-3,9)), alim=(c(0,6)))
#addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
#        rows=0, cex=1)#Add pooled
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
#png(file="20230328 Non-neurologic abnormalities zika negative.png",width=750,height=500,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Non-neurologic abnormalities, zika-negative women",
                digits=2L,xlim=(c(-3,9)), alim=(c(0,6)))
#addpoly(x = pool.outcome$abs.risk, ci.lb=pool.outcome$ci.lb, ci.ub=pool.outcome$ci.ub,
#        rows=0, cex=1)#Add pooled
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

#Store results
results<-cbind(results,pool.outcome.rr)

######One-stage meta-analysis relative risk
data.1sma<-data[!is.na(data$nonneurologic) & (data$end_ga_compl>=24 & !is.na(data$end_ga_compl)),]
#One stage meta-analysis with log link and random intercept per study, per imputed dataset
fit1.coef<-as.data.frame(f.1ma.r.int(data.1sma,"nonneurologic",poisson))
#Pool with rubins rules
pool.outcome.rr.1ma<-f.1ma.poolrubin(fit1.coef)

#Forest plot
#png(file="20230328 Non-neurologic abnormalities RR.png",width=750,height=500,res=100)
metafor::forest(rr.outcome$rr, ci.lb=rr.outcome$ci.lb, ci.ub=rr.outcome$ci.ub, 
                refline = 1, slab = rr.outcome$studyname,
                xlab = "Relative risk", pch = 19, psize=1,
                ylim=(c(-1,nrow(rr.outcome)+3)), cex=1, steps=6, main="Non-neurologic abnormalities",
                digits=2L,xlim=(c(-5,10)), alim=(c(0,5)))
addpoly(x = pool.outcome.rr$rr, ci.lb=pool.outcome.rr$rr.ci.lb, ci.ub=pool.outcome.rr$rr.ci.ub,
        rows=0, cex=1)#Add pooled
addpoly(x = pool.outcome.rr.1ma$rr.1ma, ci.lb=pool.outcome.rr.1ma$rr.ci.lb.1ma, ci.ub=pool.outcome.rr.1ma$rr.ci.ub.1ma,
        rows=-1, cex=1)
#addpoly(x=PI[,1], ci.lb=PI[,2], ci.ub=PI[,3], rows=-1, cex=1) #Add prediction interval
#dev.off()

#Store results
results<-cbind(results,pool.outcome.rr.1ma)
results$outcome<-"Non-neurologic abnormalities"
all.results<-rbind(all.results,results)

######################
#write.xlsx2(all.results, "20230328 Results objective 1.xls", col.names = T, row.names = F)

