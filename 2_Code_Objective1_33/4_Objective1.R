###Aim: Graph for comparing estimations between raw and imputed datasets.

rm(list=ls()) # clean environment

# Load packages ---
# Data manipulation package
library(data.table)
library(dplyr)
library(magrittr)
#setwd("/Users/jdamen/Documents/GitHub/Zika_imputation") #Only for Anneke
library(here)  # define folder paths
library(mice)
# Graphic packages
library(ggplot2)

# Field specific package
library(growthstandards)

source(here('2_Code_Objective1_33','4.1_Functions_Objective1.R'))

# Get original data ----
load(file=here('3_Output_data','ori_data33.RData'))
data_ori$source<-"Original"
data_ori$.imp <- -1

# Get raw data ----
load(file= here('3_Output_data','det_data33.RData')) # fdata=Data_preimputation
data_det$source <- "Deterministic"
data_det$.imp <- 0
data_abort <- data_det[,c("childid","repabort")]

# Get imputation data ----
load("/Users/jmunozav/Desktop/Zika_Jun22/merged_imp.RData")
#load(file= here('3_Output_data','merged_imp.RData')) #Only for Anneke
data_imp <- as.data.table(complete(merged_imp,"long"))
data_imp$source <- "Imputation"
data_imp$.id <- NULL

# Convert factors of 2 levels to numeric if need it
nam_imp<-colnames(data_imp)
nam_nfact<-c(".imp","childid","studyimp","end_ga","ch_weight","ch_head_circ_birth","age","gravidity","parity","ch_microcephaly"," parity","zikv_test_ev","zikv_ga","source")
nam_num<-nam_imp[!nam_imp%in%nam_nfact]
data_imp[, (nam_num) := lapply(.SD, function(x){as.numeric(levels(x))[x]}), .SDcols = nam_num]


# Select only the imputed variables on the original dataset and transform to numeric if it is required
nam_new<-c("arb_preg","birth","end_ga","arb_ever","storch_patho","comorbid_preg","drugs_prescr","zikv_test_ev","neuroabnormality","flavi_alpha_virus", "microcephaly_bin_postnatal","microcephaly_bin_fet","nonneurologic","any_abnormality_czs","zikv_test_ev")
data_ori<-as.data.frame(data_ori[,nam_imp[!nam_imp%in%nam_new], with = FALSE])
namcol <- colnames(data_ori)
nam_num <- namcol[!namcol%in%c("childid","ch_microcephaly","source")]
setDT(data_ori)[, (nam_num) := lapply(.SD, as.numeric), .SDcols = nam_num]
data_ori[,ch_microcephaly:=as.factor(ch_microcephaly)]

# Select only the imputed variables on the deterministic dataset and transform to numeric if it is required
data_det<-as.data.frame(data_det[,nam_imp, with = FALSE])
namcol <- colnames(data_det)
nam_num <- namcol[!namcol%in%c("childid","ch_microcephaly","zikv_test_ev","source")]
setDT(data_det)[, (nam_num) := lapply(.SD, as.numeric), .SDcols = nam_num]


# Merge both datasets ----

study_info <- (readxl::read_xlsx(here('1_Input_data','MasterCodebook_October.xlsx'),sheet="StudyID"))
data_N<-data_ori[, .(N = .N), by = studyimp]
data_all <- dplyr::bind_rows(data_ori,data_det,data_imp)
data_all <- merge(data_all,data_abort,by="childid")
data_all <- merge(data_all,data_N,by="studyimp")
data_all <- merge(data_all,study_info[,c("studyimp","studyname","Included")],by="studyimp")
 


# Create additional variables 
data_all[,bdeath := 1-birth]
data_all[,loss := ifelse(bdeath==1&end_ga>=20,1,0)]
data_all[,efdeath := ifelse(bdeath==1&end_ga>=20&end_ga<28,1,0)] #early fetal death (20-27 weeks gestation)
data_all[,lfdeath := ifelse(bdeath==1&end_ga>=28,1,0)] #late fetal death
data_all[,cend_ga:=ifelse(end_ga>42,42,end_ga)] # modified end_ga used to calculate microcephaly and microchephaly_bin_birth from head circunference
data_all[,ch_sex := as.factor(ifelse(ch_sex== 0, "Male","Female"))]   
data_all[!is.na(ch_sex)&!is.na(cend_ga), hcircm2zscore:=as.numeric(igb_hcircm2zscore(gagebrth = cend_ga*7, hcircm=ch_head_circ_birth,sex=ifelse(ch_sex== 0, "Male","Female")))]
data_all[, microcephaly  := as.factor(ifelse(hcircm2zscore<=-3,2,ifelse(hcircm2zscore<=-2,1,ifelse(hcircm2zscore<=2,0,ifelse(!is.na(hcircm2zscore),3,NA)))))] # given by formula
data_all[, microcephaly_bin_birth:= ifelse(microcephaly%in%c(0,3),0,ifelse(microcephaly%in%c(1,2),1,NA))]
data_all[, who_czs := as.numeric((zikv_test_ev %in% c("Robust","Moderate")| microcephaly_bin_fet==1) & (microcephaly==2 | any_abnormality_czs==1))]


# Separate dataset according to zikv_preg

data_zika <- data_all%>%filter(zikv_preg ==1)
data_nozika <- data_all%>%filter(zikv_preg ==0)


# Microcephaly only population between 24 and 42 gestational age without cases that were reported as abortion

mic_data_all <- data_all%>%filter(cend_ga>=24&cend_ga<=42&repabort!=1)
mic_data_zika <- mic_data_all%>%filter(zikv_preg ==1)
mic_data_nozika <- mic_data_all%>%filter(zikv_preg ==0)


# Miscarriage
mis_data_all <- merge(data_all,study_info[,c("studyimp","Included_mis")],by="studyimp")
mis_data_all[,Included:=Included_mis]
mis_data_zika <- mis_data_all%>%filter(zikv_preg ==1)
mis_data_nozika <- mis_data_all%>%filter(zikv_preg ==0)




# Get the estimates (plots and tables) ----
print_obj1(outcome_name="microcephaly_bin_birth",gentitle = "Microcephaly at birth calculated with head circunference",data_all=mic_data_all,data_zika=mic_data_zika,data_nozika =mic_data_nozika)
print_obj1(outcome_name="ch_microcephaly_bin",gentitle = "Microcephaly at birth from study information",data_all=mic_data_all,data_zika=mic_data_zika,data_nozika =mic_data_nozika)
print_obj1(outcome_name="miscarriage",gentitle = "Miscarriage",data_all=mis_data_all,data_zika=mis_data_zika,data_nozika =mis_data_nozika)
print_obj1(outcome_name="loss",gentitle = "Loss",data_all=data_all,data_zika=data_zika,data_nozika =data_nozika)
print_obj1(outcome_name="ch_czs",gentitle = "Child congenital zika",data_all=data_all,data_zika=data_zika,data_nozika =data_nozika)
print_obj1(outcome_name="who_czs",gentitle = "Child congenital zika (WHO)",data_all=data_all,data_zika=data_zika,data_nozika =data_nozika)





