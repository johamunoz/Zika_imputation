###Aim: Graph for comparing estimations between raw and imputed datasets.

rm(list=ls()) # clean environment

# Load packages ---
# Data manipulation package
library(dplyr)
library(magrittr)
library(here)  # define folder paths
library(mice)
# Graphic packages
library(ggplot2)

source(here('2.1_Code 33 studies','F_graph.R'))

# Get imputation data ----
load("/Users/jmunozav/Desktop/merged_imp.RData")
data_imp <- complete(merged_imp,"long")
data_imp$source <- "Imputation"
data_imp$.id <- NULL

# Convert factors of 2 levels to numeric if need it
data_imp %<>% 
  mutate(across(where(is.factor)&-c(childid,microcephaly,zikv_test_ev,source),
                ~as.numeric(as.character(.x))))

# Get raw data ----
load(file= here('3_Output_data','rawfinaldata33.RData')) # fdata=Data_preimputation
data_raw$source <- "Raw"
data_raw$.imp <- 0

# Select only imputed columns and transform to numeric if it is required
col_imp<-colnames(data_imp)
col_imp<-col_imp[col_imp!="who_czs"]
data_raw %<>% 
          select(all_of(col_imp))%>%
          mutate(across(where(is.factor)&-c(childid,microcephaly,zikv_test_ev,source),
                ~as.numeric(as.character(.x))))%>%
          mutate(across(where(is.character)&-c(childid,microcephaly,zikv_test_ev,source),
                ~as.numeric(.x)))
data_raw%<>%mutate(microcephaly=as.factor(microcephaly))


data_N <- data_raw%>%
  nest(data=-studyimp)%>%
  mutate(meta=map(data,~.x%>%dplyr::summarize(N=n())))%>%
  select(c(studyimp,meta))%>%
  unnest(col=c(meta))


# Get original data ----
load(file = here('3_Output_data','originaldata.RData')) # data 
data_ori$source <- "Ori"
data_ori$.imp <- -1

# Select only imputed columns and transform to numeric if it is required
col_imp <- c(colnames(data_imp))
col_ori <- col_imp[!col_imp %in% c('arb_preg','arb_ever','comorbid_preg','zikv_test_ev','who_czs','neuroabnormality','microcephaly_bin_postnatal','nonneurologic')]
data_ori %<>% 
  select(all_of(col_ori))%>%
  mutate(across(where(is.factor)&-c(childid,microcephaly,source),
                ~as.numeric(as.character(.x))))%>%
  mutate(across(where(is.character)&-c(childid,microcephaly,source),
                ~as.numeric(.x)))



# Merge both datasets ----
data_all <- dplyr::bind_rows(data_ori,data_raw,data_imp)

# Create additional variables 
data_all%<>%
  mutate(bdeath = 1-birth,
         loss = ifelse(bdeath==1&end_ga>=20,1,0),
         efdeath = ifelse(bdeath==1&end_ga>=20&end_ga<28,1,0), #early fetal death (20-27 weeks gestation)
         lfdeath = ifelse(bdeath==1&end_ga>=28,1,0), #late fetal death
         who_czs = as.numeric((zikv_test_ev %in% c("Robust","Moderate")| microcephaly_bin_fet==1) & (microcephaly==2 | any_abnormality_czs==1))) #WHO definition for CZS: Presence of confirmed maternal or fetal ZIKV infection AND (presence of severe microcephaly at birth OR presence of other malformations (including limb contractures, high muscle tone, eye abnormalities, and hearing loss, nose etc.))
         

# Get studynames
study_info <- (readxl::read_xlsx(here('1_Input_data','MasterCodebook_October.xlsx'),sheet="StudyID"))
data_all <- merge(data_all,study_info%>%select(c("studyimp","studyname","Included")),by="studyimp",all.x=TRUE)
data_all <- merge(data_all,data_N,by="studyimp",all.x=TRUE)    


# Separate dataset according to zikv_preg
data_zika <- data_all%>%filter(zikv_preg ==1)
data_nozika <- data_all%>%filter(zikv_preg ==0)



# Microcephaly only population between 24 and 42 gestational age
mdata_all <- data_all%>%filter(end_ga>=24&end_ga<=48)
mdata_zika <- data_zika%>%filter(end_ga>=24&end_ga<=48)
mdata_nozika <- data_nozika%>%filter(end_ga>=24&end_ga<=48)


summary(data_all)


forest_plot_study(data=data_all,outcome_name="microcephaly_bin_birth",plottitle = "Microcephaly at birth,all")
forest_plot_study(data=data_zika,outcome_name="microcephaly_bin_birth",plottitle = "Microcephaly at birth,zika+ mom")
forest_plot_study(data=data_nozika,outcome_name="microcephaly_bin_birth",plottitle = "Microcephaly at birth,zika- mom")

forest_plot_study(data=data_all,outcome_name="zikv_preg",plottitle = "Zika diagnosis in mothers")






