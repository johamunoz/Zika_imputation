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

load(here('3_Output_data','finaldata33.RData'))
add_info <- as.data.table(readxl::read_xlsx(here('1_Input_data','MasterCodebook_October.xlsx'),sheet="237 key")) #CSV file with the
imp_info<-add_info[Final_imputation=="yes"]
data <- fdata

# Assure format of variables
nncon<-imp_info[!Type_var%in%c("Continuous","ID")]$who_name
data[,(nncon):= lapply(.SD, as.factor), .SDcols = nncon] #convert to factors binomial and categorical variables
summary(data)


con_var1<-imp_info[Type_var=="Continuous"&Type_imp==1]$who_name
con_var2<-imp_info[Type_var=="Continuous"&Type_imp==2]$who_name
bin_var1<-imp_info[Type_var=="Binary"&Type_imp==1]$who_name
bin_var2<-imp_info[Type_var=="Binary"&Type_imp==2]$who_name
cat_var1<-imp_info[Type_var=="Categorical"&Type_imp==1]$who_name
cat_var2<-imp_info[Type_var=="Categorical"&Type_imp==2]$who_name
npred<-imp_info[npred==1]$who_name


# After checking the histogram we decide to impute this variables as pmm
cat_var2<-c(cat_var2,"gravidity","end_ga")


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

meth["zikv_preg.end_ga"] <- "~ifelse(zikv_preg == '1', end_ga-38, 0)"
meth["zikv_preg.zikv_ga"] <- "~ifelse(zikv_preg == '1', zikv_ga-20, 0)"
meth["zikv_preg.denv_preg_ever"] <- "~ifelse(zikv_preg == '1'& denv_preg_ever =='1','1','0')"
meth["zikv_preg.chikv_preg_ever"] <- "~ifelse(zikv_preg == '1'& chikv_preg_ever =='1','1','0')"
meth["zikv_preg.comorbid_bin"] <- "~ifelse(zikv_preg == '1'& comorbid_bin =='1','1','0')" 
meth["zikv_preg.comorbid_preg"] <- "~ifelse(zikv_preg == '1'& denv_preg_ever =='1','1','0')"            
meth["zikv_preg.storch_patho"] <- "~ifelse(zikv_preg == '1'& storch_patho =='1','1','0')"     
meth["zikv_preg.arb_symp"] <- "~ifelse(zikv_preg == '1'& arb_symp =='1','1','0')"         
meth["zikv_preg.fever"] <- "~ifelse(zikv_preg == '1'& fever =='1','1','0')"       
meth["zikv_preg.rash"] <- "~ifelse(zikv_preg == '1'& rash =='1','1','0')"
meth["zikv_preg.arthralgia"] <- "~ifelse(zikv_preg == '1'& arthralgia =='1','1','0')"       
meth["zikv_preg.headache"] <- "~ifelse(zikv_preg == '1'& headache =='1','1','0')"    
meth["zikv_preg.arthritis"] <- "~ifelse(zikv_preg == '1'& arthritis =='1','1','0')"

#2.2. Post-process the values to constrain normal distributed variables in the range of observable values ----
post["zikv_assay_ga_1"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 50))"
post["zikv_ga"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 46))"
post["zikv_preg.end_ga"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 46))"
post["zikv_preg.zikv_ga"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 46))"
post["end_ga"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 46))"
post["bmi"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 50))"
post["age"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(14, 55))"
post["ch_weight"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(100, 6000))"
post["ch_head_circ_birth"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(15, 55))"


#3.3. Set prediction matrix -----
pred <- quickpred(data, minpuc = 0.3) # assignation based on pairwise correlation
pred["childid",] <- 0
pred[,"childid"] <- 0
pred[,"studyimp"] <- -2 # define the cluster for imputation models at study level.

varinter=c("end_ga","zikv_ga","denv_preg_ever","chikv_preg_ever","comorbid_bin","comorbid_preg","storch_patho","arb_symp","fever","rash","arthralgia","headache","arthritis")
colpred<-colnames(pred)
predata<-as.data.table(pred)
for(i in varinter) {
  cn<-paste0("zikv_preg.",i)
  exp<-paste0("ifelse(zikv_preg==1&",i,"==1,1,0)")
  expr <- parse(text = paste0(cn, ":=",exp))
  predata[,eval(expr)]
}

pred <- as.matrix(predata)
rownames(pred) <- colpred
pred[c("zikv_preg","end_ga"), c("zikv_preg.end_ga")] <- 0
pred[c("zikv_preg","zikv_ga"), c("zikv_preg.zikv_ga")] <- 0
pred[c("zikv_preg","denv_preg_ever"), c("zikv_preg.denv_preg_ever")] <- 0
pred[c("zikv_preg","chikv_preg_ever"), c("zikv_preg.chikv_preg_ever")] <- 0
pred[c("zikv_preg","comorbid_bin"), c("zikv_preg.comorbid_bin")] <- 0
pred[c("zikv_preg","comorbid_preg"), c("zikv_preg.comorbid_preg")] <- 0
pred[c("zikv_preg","storch_patho"), c("zikv_preg.storch_patho")] <- 0
pred[c("zikv_preg","arb_symp"), c("zikv_preg.arb_symp")] <- 0
pred[c("zikv_preg","fever"), c("zikv_preg.fever")] <- 0
pred[c("zikv_preg","rash"), c("zikv_preg.rash")] <- 0
pred[c("zikv_preg","arthralgia"), c("zikv_preg.arthralgia")] <- 0
pred[c("zikv_preg","headache"), c("zikv_preg.headache")] <- 0
pred[c("zikv_preg","arthritis"), c("zikv_preg.arthritis")] <- 0


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

mice1<-fun_run(imp=1,data=data,pred=pred,meth=meth,post=post,maxit=10)
save(mice1,file=here('3_Output_data','mice1.Rdata'))
# 
mice2<-fun_run(imp=2,data=data,pred=pred,meth=meth,post=post,maxit=10)
save(mice2,file=here('3_Output_data','mice2.Rdata'))
# 
mice3<-fun_run(imp=3,data=data,pred=pred,meth=meth,post=post,maxit=10)
save(mice3,file=here('3_Output_data','mice3.Rdata'))
# 
mice4<-fun_run(imp=4,data=data,pred=pred,meth=meth,post=post,maxit=10)
save(mice4,file=here('3_Output_data','mice4.Rdata'))

mice5<-fun_run(imp=5,data=data,pred=pred,meth=meth,post=post,maxit=10)
save(mice5,file=here('3_Output_data','mice5.Rdata'))
# # 
# mice6<-fun_run(imp=6,data=data,pred=pred,meth=meth,post=post,maxit=10)
# save(mice6,file=here('3_Output_data','mice6.Rdata'))
# 
# mice7<-fun_run(imp=7,data=data,pred=pred,meth=meth,post=post,maxit=10)
# save(mice7,file=here('3_Output_data','mice7.Rdata'))
# 
# mice8<-fun_run(imp=8,data=data,pred=pred,meth=meth,post=post,maxit=10)
# save(mice8,file=here('3_Output_data','mice8.Rdata'))
# 
# mice9<-fun_run(imp=9,data=data,pred=pred,meth=meth,post=post,maxit=10)
# save(mice9,file=here('3_Output_data','mice9.Rdata'))
# 
# mice10<-fun_run(imp=10,data=data,pred=pred,meth=meth,post=post,maxit=10)
# save(mice10,file=here('3_Output_data','mice10.Rdata'))



