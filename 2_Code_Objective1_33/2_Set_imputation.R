#Aim: To set up the imputation methods

rm(list=ls()) # clean environment

#setwd("/Users/jdamen/Documents/GitHub/Zika_imputation") ###Anneke only

# Load packages ---
# Data manipulation package
library(data.table) 
library(here)  # define folder paths

#Imputation packages
library(mice)
library(micemd)
library(miceadds)

# Load dataset and dependencies ----

#load("/home/julius_te/jmunoz/Run_zika/finaldata.RData")
load(here('3_Output_data','finaldata33.RData'))
#add_info <- as.data.table(readxl::read_xlsx("/home/julius_te/jmunoz/Run_zika/MasterCodebook_October.xlsx",sheet="237 key")) #CSV file with the
add_info <- as.data.table(readxl::read_xlsx(here('1_Input_data','MasterCodebook_October.xlsx'),sheet="237 key")) #CSV file with the
add_info <- add_info[Imp_obj1=="yes"]
imp_info<-add_info[order(as.numeric(add_info$Orderimp1)),]

data <- fdata[, imp_info$who_name, with=FALSE] #set order imputation

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
cat_var1<-c(cat_var1,"zikv_ga","gravidity","parity")

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
meth[bin_var2] <- "2l.glm.bin" # normal distribution at study level
meth[cat_var1] <- "pmm"#  binomial distribution at study level
meth[cat_var2] <- "2l.pmm" # predictive mean matching at study level
meth["ch_microcephaly_bin"]<-"~ifelse(ch_microcephaly%in%c(0,3),0,1)"
meth["miscarriage"]<-"~ifelse(birth == '0'& end_ga < 20, 1, 0)"

#2.2. Post-process the values to constrain normal distributed variables in the range of observable values ----
post["age"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(14, 55))"
post["zikv_ga"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 46))"
post["ch_weight"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(100, 6000))"
post["end_ga"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 46))"
post["ch_head_circ_birth"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(15, 55))"
post["microcephaly_bin_postnatal"] <- paste(sep = ";", 'idx <- row.names(imp[["birth"]][imp[["birth"]][,i]==0,])','imp[[j]][idx, i] <- 2')
post["any_abnormality_czs"] <- paste(sep = ";", 'idx <- row.names(imp[["neuroabnormality"]][imp[["neuroabnormality"]][,i]==1,]|
                                                                  imp[["nonneurologic"]][imp[["nonneurologic"]][,i]==1,])','imp[[j]][idx, i] <- 1')
post["parity"] <- "imp[[j]][, i] <- pmin(imp[[j]][, i], data[(!r[, j]) & where[, j], 'gravidity'])"


#2.3. Set prediction matrix -----
pred <- quickpred(data, minpuc = 0.2) # assignation based on pairwise correlation
pred[,"studyimp"] <- -2 # define the cluster for imputation models at study level.
pred["childid",] <- 0
pred[,"childid"] <- 0
pred["ch_microcephaly_bin","ch_microcephaly"] <- 1
pred["ch_microcephaly","ch_microcephaly_bin"] <- 0
pred["zikv_preg","zikv_test_ev"]<-0
pred["zikv_test_ev","zikv_preg"]<-0
pred["end_ga","zikv_ga"]<-1
pred["miscarriage","birth"] <- 1
pred["miscarriage","end_ga"] <- 1
pred[c("end_ga","birth"),"miscarriage"] <- 0


#write.csv(pred,file=here('Documents','GitHub','Zika_imputation','3_Output_data','predmatrix33.csv'))
write.csv(pred,file=here('3_Output_data','predmatrix33.csv'))


fun_run<-function(imp,data,pred,meth,post,maxit){
  seed=imp*363+253  
  start_time<-proc.time()
  micesurv <- mice(data, predictorMatrix = pred,method=meth,post=post, maxit =maxit, m = 1, printFlag = TRUE, seed = seed)
  end_time<-proc.time()
  time = end_time-start_time
  return(list(micesurv,time))}

###END CODE##




