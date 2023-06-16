
#Passive imputation data provided 09-09-21

rm(list=ls())

# Data manipulation packages
library(rstudioapi) 
library(data.table)
library(dplyr)

# Imputation packages
library(mitml)
library(mice)
library(micemd)
library(miceadds)

# Graphic packages
library(corrplot)
library(ggplot2)
library(ggtext)
library(plotly)

# Field specific package
library(growthstandards)

#0. Load information ----

base.dir <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path)) # set main base working directory
setwd(base.dir)
data<-as.data.table(read.csv('1_Input_data/pilot10_08SEP21_withmetadata.csv', stringsAsFactors=FALSE, fileEncoding="latin1"))


infoexp<-as.data.table(readxl::read_xlsx("1_Input_data/Infoexp.xlsx",sheet="Table")) #Table were are specified the included variables according Expert opinion, also the includes the order in which variables are imputed 
var_inc<-infoexp[Inclusion==1,Variable] #Variables to work with
data<-as.data.table(data[,..var_inc]) #Filter dataset

#1. Initial clean up the dataset -----
#1.0. Check for duplicates ----
data[,comb:=paste0(studyname,mid_original,sep="_")] #create a variable that combine study and patient id
n_occur <- data.table(table(data$comb))
n_occur[n_occur$N>1] #check visually which observations are duplicated
var1<-n_occur[n_occur$N > 1&V1!='_Spain_Soriano'& V1!='_TrinidadTobago_Sohan']$V1 
var1<-n_occur[n_occur$N > 1]$V1
Duplicates<-data[comb%in%var1&multiplegest!=1] # table to report to harmonization team.
write.csv(Duplicates,file='5_Internal_support/Possible_duplicates.csv')

# 1.1. Check format of categorical variables----
char <-   which(sapply(data, is.character)) #which are char
datac<-data[,..char]
summary(datac)
unique(datac$inf_head_circ_age_fu1)
data[,inf_head_circ_age_fu1:=as.numeric(sub(" .*", "",data$inf_head_circ_age_fu1))] # remove years
data[,zikv_elisa_ga_1:=as.numeric(sub(" .*", "",data$zikv_elisa_ga_1))] #remove weeks
summary(data$zikv_prnt_1)


#unique(data$othabnorm_spec) # 39 categories ask if really want to include this
#data[,othabnorm_spec:=NULL] # We temporal remove it from dataset as it is messy
data[zikv_pcr_date_1%in%c("22027","888",""), zikv_pcr_date_1:=NA]
data[zikv_elisa_date_1%in%c(""),zikv_elisa_date_1:=NA]
data[symp_date%in%c(""),symp_date:=NA]
data[,arb_clindiag:=ifelse(arb_clindiag==777,6,arb_clindiag)] # We set the NA value ("777") to a level "6" because we later use it for define CZS variable

#1.2. Set to NA : 666,777,888,999 values----
data[data==""] <-NA
data[data==666] <-  NA
data[data==777] <-  NA
data[data==888] <-  NA
data[data==999] <-  NA
data[data==9999] <-  NA

#1.3 
summary(data$zikv_prnt_1)
summary(data$zikv_prnt_everpos_studydef)
summary(data$zikv_prnt_studydef_1)
summary(data$zikv_prnt_titer_1)
data$zikv_pcr_ga_1
data_test<-data[,c("zikv_elisa_res_1","zikv_elisa_ga_1","zikv_elisa_date_1","zikv_prnt_studydef_1","zikv_pcr_res_1","zikv_pcr_ga_1","studyname")]
summary(data_test)

data_test[,zikv_test:=fcase(zikv_pcr_res_1==1|(zikv_elisa_res_1==1&zikv_prnt_studydef_1==1),"robust",
                            zikv_elisa_res_1==1|zikv_prnt_studydef_1==1, "moderate_limit",
                            zikv_pcr_res_1==0&zikv_elisa_res_1==0,"negative",
                            zikv_pcr_res_1==0|zikv_elisa_res_1==0,"any_negative")]

data_test_sum<-as.data.table(table(pcr=data_test$zikv_pcr_res_1,elisa=data_test$zikv_elisa_res_1,prnt=data_test$zikv_prnt_studydef_1,test=data_test$zikv_test,useNA = "always"))
results_test<-data_test_sum[N>0]

table(data$studyname,data_test$zikv_test)
