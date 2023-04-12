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

# Get original data ----
  load(file= here('3_Output_data','originaldata.RData')) # data 
  data_ori$source <- "Ori"
  data_ori$.imp <- -1
  
  # Select only imputed columns and transform to numeric if it is required
  col_imp <- colnames(data_imp)
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

  
  data_raw %<>% 
    select(all_of(col_imp))%>%
    mutate(across(where(is.factor)&-c(childid,microcephaly,zikv_test_ev,source),
                  ~as.numeric(as.character(.x))))%>%
    mutate(across(where(is.character)&-c(childid,microcephaly,zikv_test_ev,source),
                  ~as.numeric(.x)))


# Merge both datasets ----
  data_all <- dplyr::bind_rows(data_ori,data_raw,data_imp)

# Create additional variables 
  data_all%<>%
    mutate(bdeath = 1-birth,
           miscarriage = ifelse(bdeath==1&end_ga<20,1,NA),
           loss = ifelse(bdeath==1&end_ga>=20,1,NA),
           efdeath = ifelse(bdeath==1&end_ga>=20&end_ga<28,1,ifelse(is.na(bdeath)|is.na(end_ga),NA,0)), #early fetal death (20-27 weeks gestation)
           lfdeath = ifelse(bdeath==1&end_ga>=28,1,ifelse(is.na(bdeath)|is.na(end_ga),NA,0)), #late fetal death
           microcephaly_bin_birth:= ifelse(microcephaly%in%c(0,3),0,ifelse(microcephaly%in%c(1,2),1,NA))) #microcephaly at birth
  
# Get studynames
  study_info <- (readxl::read_xlsx(here('1_Input_data','MasterCodebook_October.xlsx'),sheet="StudyID"))
  data_all <- merge(data_all,study_info%>%select(c("studyimp","studyname","Included")),by="studyimp",all.x=TRUE)


# Separate dataset according to zikv_preg
data_zika <- data_all%>%filter(zikv_preg ==1)
data_nozika <- data_all%>%filter(zikv_preg ==0)



# Microcephaly at birth (Absolute risk )
MicroAR <- Rpool_studies(data =data_all, # for estimates only on zika+ mother use here: data_zika or zika- mother: data_nozika
                         outcome_name="microcephaly_bin_birth", # it can be used also for: "microcephaly_bin_postnatal","microcephaly_bin_fet","ch_czs","who_czs","neuroabnormality","nonneurologic","miscarriage","loss","efdeath","lfdeath"
                         exposure_name = NA,
                         estimand= "AR",
                         plottitle= "Microcephaly at birth, all mom, logit",
                         t_type= "logit",
                         dupper = NA)

MicroAR$plot # plot 
MicroAR$sum_tdata # plot 

# Microcephaly at birth (Relative risk )

MicroRR<- Rpool_studies(data = data_all,
                        outcome_name ="microcephaly_bin_birth",
                        exposure_name = "zikv_preg",
                        estimand = "RR",
                        plottitle= "Microcephaly at birth, all mom",
                        t_type= "log",
                        mod_type ="binomial",
                        correction = "Sweeting",
                        dupper = 600)
MicroRR$plot # plot 
MicroRR$sum_tdata # plot 



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


