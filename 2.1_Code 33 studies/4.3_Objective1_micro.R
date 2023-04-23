# Aim: Graph for comparing estimations between raw and imputed datasets.

rm(list=ls()) # clean environment

# Load packages ----
# Data manipulation package
library(dplyr)
library(magrittr)
library(here)  # define folder paths
library(mice)


# Load the additional functions ----
  source(here('2.1_Code 33 studies','F_graph.R'))

# Get imputation data ----
  load(file= here('3_Output_data','merged_imp.RData')) # imputed data 
  data_imp <- complete(merged_imp,"long")
  data_imp$source <- "Imputation"
  data_imp$.id <- NULL
 
# Convert factors of 2 levels to numeric if need it
  data_imp %<>% 
    mutate(across(where(is.factor)&-c(childid,microcephaly,zikv_test_ev,source),
                  ~as.numeric(as.character(.x))))
  
  data_imp%<>%
    mutate(miscarriage = ifelse(is.na(birth)|is.na(end_ga),NA,ifelse(birth==0&end_ga<20,1,0)),
           loss = ifelse(is.na(birth)|is.na(end_ga),NA,ifelse(birth==0&end_ga>=20,1,NA)),
           microcephaly_bin_birth:= ifelse(microcephaly%in%c(0,3),0,ifelse(microcephaly%in%c(1,2),1,NA)))
  data_imp<-setDT(data_imp)

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

 
  
  
# Get raw data ----
  load(file= here('3_Output_data','rawfinaldata33.RData')) # fdata=Data_preimputation
  data_raw$source <- "Raw"
  data_raw$.imp <- 0
  col_raw <- c(col_imp)
  
  data_raw %<>% 
    select(all_of(col_imp))%>%
    mutate(across(where(is.factor)&-c(childid,microcephaly,zikv_test_ev,source),
                  ~as.numeric(as.character(.x))))%>%
    mutate(across(where(is.character)&-c(childid,microcephaly,zikv_test_ev,source),
                  ~as.numeric(.x)))

    data_N <- data_raw%>%
              nest(data=-studyimp)%>%
              mutate(meta=map(data,~.x%>%dplyr::summarize(N=n())))%>%
              select(c(studyimp,meta))%>%
              unnest(col=c(meta))
 
 
# Merge source datasets ----
  data_all <- dplyr::bind_rows(data_ori,data_raw,data_imp)

# Create additional variables 
  data_all%<>%
    mutate(efdeath = ifelse(birth==0&end_ga>=20&end_ga<28,1,ifelse(is.na(birth)|is.na(end_ga),NA,0)), #early fetal death (20-27 weeks gestation)
           lfdeath = ifelse(birth==0&end_ga>=28,1,ifelse(is.na(birth)|is.na(end_ga),NA,0))) #late fetal death
  
# Get studynames
  study_info <- (readxl::read_xlsx(here('1_Input_data','MasterCodebook_October.xlsx'),sheet="StudyID"))
  data_all <- merge(data_all,study_info%>%select(c("studyimp","studyname","Included")),by="studyimp",all.x=TRUE)
  data_all <- merge(data_all,data_N,by="studyimp",all.x=TRUE)            


# Separate dataset according to zikv_preg
data_zika <- data_all%>%filter(zikv_preg == 1)
data_nozika <- data_all%>%filter(zikv_preg == 0)


# Get the estimates (plots and tables) ----
print_obj1(outcome_name="microcephaly_bin_birth",gentitle = "Microcephaly at birth")
print_obj1(outcome_name="miscarriage",gentitle = "Miscarriage",dupperRR=700)
print_obj1(outcome_name="loss",gentitle = "Loss")
print_obj1(outcome_name="ch_czs",gentitle = "Child congenital zika")
print_obj1(outcome_name="who_czs",gentitle = "Child congenital zika (WHO)")



# other available outcome_name: 

# "microcephaly_bin_postnatal" (Microcephaly at birth)
# "microcephaly_bin_fet" ("Microcephaly at fetus state")
# "ch_czs" ("Child congenital zika")
# "who_czs" ("Child congenital zika according to WHO definition")
# "neuroabnormality" ("Neuro abnormality")
# "nonneurologic" ("Non neurologic abnormalities")
# "miscarriage" ( "Miscarriage")
# "loss" ("Loss")
#  efdeath" ("Early fetus death")
# "lfdeath" ("Late fetus death")


# Check differences microcephaly among sources---

table(data_imp[studyimp==14]$microcephaly_bin_birth,data_imp[studyimp==14]$.imp,useNA = "always")
table(data_ori[studyimp==14]$microcephaly_bin_birth,useNA = "always")
table(data_raw[studyimp==14]$microcephaly_bin_birth,useNA = "always")
table(data_ori[studyimp==14]$microcephaly_bin_birth,
      data_ori[studyimp==14]$microcephaly,
      useNA = "always")

table(data_raw[studyimp==14]$microcephaly_bin_birth,
      data_raw[studyimp==14]$microcephaly,
      useNA = "always")
