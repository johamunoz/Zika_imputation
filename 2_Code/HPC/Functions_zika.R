#Aim: To set up the imputation methods

rm(list=ls()) # clean environment

# Load packages ---
# Data manipulation package
library(data.table) 

#Imputation packages
library(mice)
library(micemd)
library(miceadds)

# Load dataset and dependencies ----

load("/home/julius_te/jmunoz/Run_zika/finaldata.RData")
data<-as.data.table(fdata)
add_info <- as.data.table(readxl::read_xlsx("/home/julius_te/jmunoz/Run_zika/MasterCodebook_October.xlsx",sheet="237 key")) #CSV file with the
imp_info<-add_info[Imp_obj1=="yes"]


# Assure format of variables
nncon<-imp_info[!Type_var%in%c("Continuous","ID")]$who_name
data[,(nncon):= lapply(.SD, as.factor), .SDcols = nncon] #convert to factors binomial and categorical variables

con_var1<-imp_info[Type_var=="Continuous"&Type_imp==1]$who_name
con_var2<-imp_info[Type_var=="Continuous"&Type_imp==2]$who_name
bin_var1<-imp_info[Type_var=="Binary"&Type_imp==1]$who_name
bin_var2<-imp_info[Type_var=="Binary"&Type_imp==2]$who_name
cat_var1<-imp_info[Type_var=="Categorical"&Type_imp==1]$who_name
cat_var2<-imp_info[Type_var=="Categorical"&Type_imp==2]$who_name


# After checking the histogram we decide to impute this variables as pmm
cat_var2<-c(cat_var2,"end_ga")


#2. Set methods, prediction and post arguments----
ini <- mice(data, maxit = 0)
ini$loggedEvents
meth <- ini$meth
pred <- ini$pred
post <- ini$post


#2.1. Assign imputation method to each variable ----
meth[con_var1] <- "norm" # normal distribution at study level
meth[con_var2] <- "2l.norm" # normal distribution at study level
meth[bin_var1] <- "logreg" # normal distribution at study level
meth[bin_var2] <- "2l.bin" # normal distribution at study level
meth[cat_var1] <- "pmm"#  binomial distribution at study level
meth[cat_var2] <- "2l.pmm" # predictive mean matching at study level


#2.2. Post-process the values to constrain normal distributed variables in the range of observable values ----
post["age"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(14, 55))"
post["zikv_ga"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 46))"
post["ch_weight"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(100, 6000))"
post["end_ga"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 46))"
post["ch_head_circ_birth"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(15, 55))"


#3.3. Set prediction matrix -----
pred <- quickpred(data, minpuc = 0.1) # assignation based on pairwise correlation
pred["childid",] <- 0
pred[,"childid"] <- 0
pred[,"studyimp"] <- -2 # define the cluster for imputation models at study level.


fun_run<-function(imp,data,pred,meth,post,maxit){
  seed=imp*363+253  
  start_time<-proc.time()
  micesurv <- mice(data, predictorMatrix = pred,method=meth,post=post, maxit =maxit, m = 1, printFlag = TRUE, seed = seed)
  end_time<-proc.time()
  time = end_time-start_time
  return(list(micesurv,time))}