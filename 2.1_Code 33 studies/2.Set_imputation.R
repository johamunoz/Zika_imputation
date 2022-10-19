#Aim: To set up the imputation methods

rm(list=ls()) # clean environment

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
add_infoi<-add_info[order(Orderimp)]


base.dir <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(base.dir)
load("3_Output_data/finaldata.RData")
infoselec<-as.data.table(read.csv2("5_Internal_support/Infoselection.csv",sep=";", stringsAsFactors=FALSE, fileEncoding="latin1"))

infoselec<-infoselec[Imputation2>=3,] #only variable selected
infoselec[Vtype=="C"]$Variable #Continuous
infoselec[Vtype=="CAT"]$Variable #Categorical
infoselec[Additionals==3]$Variable # Included in imputation but not as predictors


data<-as.data.table(fdata)
colnames(data)


ncon<-c("age","bdeath_ga","zikv_ga","inf_weight","inf_head_circ_birth","inf_length") # name of continuous variables
nncon<-colnames(data)[!colnames(data)%in%ncon] # name of no continuous variables
data[,(nncon):= lapply(.SD, as.factor), .SDcols = nncon] #convert to factors binomial and categorical variables
data[,studyname:=as.integer(studyname)] # Convert cluster variable studyname to integer 
summary(data)

ncat<-c("educ","ethnicity","gravidity","maritalstat","microcephaly","parity","tobacco")#name of categorical variables
nbin<-nncon[!nncon%in%c(ncat,"studyname","Comb_id")] #name of binary variables
npred<- c("eclampsia","educ","ethnicity","fet_zikv","flavi_alpha_virus","gestdiab","inf_craniofac_abn_bin","maritalstat","tobacco") #name of no predictors


#2. Set methods, prediction and post arguments----
ini <- mice(data, maxit = 0)
meth <- ini$meth
pred <- ini$pred
post <- ini$post

head(ini$loggedEvents)

#2.1. Assign imputation method to each variable ----
meth[ncon] <- "2l.norm" # normal distribution at study level
meth[nbin] <- "2l.bin"#  binomial distribution at study level
meth[ncat] <- "2l.pmm" # predictive mean matching at study level

meth[c("inf_weight","inf_head_circ_birth")]<-"norm"    #  # Due to convergence problems imputation based on marginal distribution
meth[c("inf_length","zikv_ga")]<- "2l.pmm" #zikv_gan does not work as normal
meth[c("corticalatrophy","ventriculomegaly","multiplegest","othabnorm","calcifications","storch_patho","hydrocephaly","inducedabort")] <- "logreg" # Due to frequent singular fit warning imputation based on marginal distribution
meth[c("educ")]<- "pmm"
meth

#2.2. Post-process the values to constrain normal distributed variables in the range of observable values ----
post["age"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(14, 50))"
post["bdeath_ga"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 50))"
post["zikv_ga"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 65))"
post["inf_weight"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(100, 7000))"
post["inf_head_circ_birth"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(15, 50))"
post["inf_length"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(18, 60))"


#3.3. Set prediction matrix -----
pred <- quickpred(data, minpuc = 0.3) # assignation based on pairwise correlaition
pred["Comb_id",] <- 0
pred[,"Comb_id"] <- 0
pred[,"studyname"] <- -2 # define the cluster for imputation models at study level.
pred["studyname", "studyname"] <- 0
pred[,npred]<-0 #not assign variables without prediction power into the imputation model
write.csv(pred,file="5_Internal_support/predmatrix_imp2.csv ")

fun_run<-function(imp,data,pred,meth,post,maxit){
  seed=imp*363+253  
  start<-proc.time()
  micesurv <- mice(data, predictorMatrix = pred,method=meth,post=post, maxit =maxit, m = 1, printFlag = TRUE, seed = seed)
  end<-proc.time()
  return(list(micesurv,end-start))}

fun_run(imp=1,data=data,pred=pred,meth=meth,post=post,maxit=2)
