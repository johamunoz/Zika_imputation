rm(list=ls())

library(data.table)
library(mice)
library(micemd)
library(miceadds)

#0. Load data ----
load("finaldata.RData")

#1. Define type of variables -----
ncon<-c("age","bdeath_ga","zikv_gan","inf_weight","inf_head_circ_birth","inf_length") # name of continuous variables
nncon<-colnames(data)[!colnames(data)%in%ncon] # name of no continuous variables
ncat<-c("gravidity","parity","microcephaly") #name of categorical variables
nbin<-nncon[!nncon%in%ncat] #name of binary variables

#2. Set methods, prediction and post arguments----
ini <- mice(data, maxit = 0)
meth <- ini$meth
pred <- ini$pred
post <- ini$post

#2.1. Assign imputation method to each variable ----
meth[ncon] <- "2l.norm" # normal distribution at study level
meth[c("inf_weight","inf_head_circ_birth","inf_length","zikv_gan")]<-"norm" # Due to convergence problems imputation based on marginal distribution
meth[ncat] <- "2l.pmm" # predictive mean matching at study level
meth[c("na_b")] <- "2l.pmm" # set with ppm at study level 
meth[nbin] <- "2l.bin"#  binomial distribution at study level
meth[c("corticalatrophy","ventriculomegaly","multiplegest","othabnorm","calcifications")] <- "logreg" # Due to frequent singular fit warning imputation based on marginal distribution

#2.2. Post-process the values to constrain normal distributed variables in the range of observable values ----
post["age"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(14, 50))"
post["bdeath_ga"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 50))"
post["zikv_gan"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 65))"
post["inf_weight"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(100, 7000))"
post["inf_head_circ_birth"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(15, 50))"
post["inf_length"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(18, 60))"

#3.3. Set prediction matrix -----
pred <- quickpred(data) # assignation based on pairwise correlaition
pred[, c("bdeath_ga")] <- 0 # as information of prediction is given by nelson allen variable na_b
pred[,"studyname"] <- -2 # define the cluster for imputation models at study level.
pred["studyname", "studyname"] <- 0

#4. Run imputation ----
fun_run<-function(imp,data,pred,meth,post,maxit){
          seed=imp*363+253  
          start<-proc.time()
          micesurv <- mice(data, predictorMatrix = pred,method=meth,post=post, maxit =maxit, m = 1, printFlag = TRUE, seed = seed)
          end<-proc.time()
          return(list(micesurv,end-start))}
