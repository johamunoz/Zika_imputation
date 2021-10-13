#Check imputation parameters before sending jobs to HPC

rm(list=ls())

library(rstudioapi)
library(data.table)
library(dplyr)
library(mitml)
library(mice)
library(micemd)
library(miceadds)
library(Hmisc)
library(VIM)
library(lqmm)

base.dir <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(base.dir)
load("3_Output_data/finaldata.RData")
data<-as.data.table(fdata)

infoselec<-as.data.table(read.csv2("5_Internal_support/Infoselection.csv",sep=";", stringsAsFactors=FALSE, fileEncoding="latin1"))
infoselec<-infoselec[Imputation2>=3,] #only variable selected
infoselec[Vtype=="C"]$Variable #Continuous
infoselec[Vtype=="CAT"]$Variable #Categorical
infoselec[Additionals==3]$Variable # Included in imputation but not as predictors

ncon<-c("age","bdeath_ga","zikv_ga","inf_weight","inf_head_circ_birth","inf_length") # name of continuous variables
nncon<-colnames(data)[!colnames(data)%in%ncon] # name of no continuous variables
ncat<-c("educ","ethnicity","gravidity","maritalstat","microcephaly","parity")#name of categorical variables
nbin<-nncon[!nncon%in%ncat] #name of binary variables
npred<- c("eclampsia","educ","ethnicity","fet_zikv","flavi_alpha_virus","gestdiab","inf_craniofac_abn_bin","maritalstat","tobacco") #name of no predictors


#2. Set methods, prediction and post arguments----
ini <- mice(data, maxit = 0)
meth <- ini$meth
pred <- ini$pred
post <- ini$post

#head(ini$loggedEvents)
#2.1. Assign imputation method to each variable ----
meth[ncon] <- "2l.norm" # normal distribution at study level
meth[c("inf_weight","inf_head_circ_birth")]<-"norm"    #  # Due to convergence problems imputation based on marginal distribution

meth[nbin] <- "2l.bin"#  binomial distribution at study level
meth[c("corticalatrophy","ventriculomegaly","multiplegest","othabnorm","calcifications","storch_patho","hydrocephaly","inducedabort")] <- "logreg" # Due to frequent singular fit warning imputation based on marginal distribution

meth[ncat] <- "2l.pmm" # predictive mean matching at study level
meth[c("inf_length","zikv_ga")]<- "2l.ppm"

#2.2. Post-process the values to constrain normal distributed variables in the range of observable values ----
post["age"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(14, 50))"
post["bdeath_ga"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 50))"
post["zikv_gan"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 65))"
post["inf_weight"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(100, 7000))"
post["inf_head_circ_birth"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(15, 50))"
post["inf_length"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(18, 60))"


#3.3. Set prediction matrix -----
pred <- quickpred(data, minpuc = 0.3) # assignation based on pairwise correlaition
pred[,"studyname"] <- -2 # define the cluster for imputation models at study level.
pred["studyname", "studyname"] <- 0
pred[,npred]<-0 #not assign variables without prediction power into the imputation model
