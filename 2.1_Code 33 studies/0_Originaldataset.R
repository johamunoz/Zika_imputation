###Aim: To clean up and to filter out the variables that will be included in the imputation model

rm(list=ls()) # clean environment

# Load packages ---
# Data manipulation package
library(data.table) 
library(here)  # define folder paths
library(rio)   # import STATA file
library(stringr)

# Field specific package
library(growthstandards)

# Graphic packages
library(ggplot2)
library(plotly)
library(mice)


# Load dataset and dependencies ----
data_origin <- as.data.table(import(here('1_Input_data','zikv_033_datasets.dta'))) #  This file you can find it on the dropbox
data_origin[,studyname:=NULL] # as studyname was not assigned to all studies
add_info <- as.data.table(readxl::read_xlsx(here('1_Input_data','MasterCodebook_October.xlsx'),sheet="237 key")) #CSV file with the
study_info <- as.data.table(readxl::read_xlsx(here('1_Input_data','MasterCodebook_October.xlsx'),sheet="StudyID")) #CSV file with the
data <- merge(data_origin, study_info, by="file")
source(here('2.1_Code 33 studies','1.1_Pre_imputation_functions.R'))


# 0. Initial checks ----
# 0.1 Missing data observations----
data[data==""] <-   NA
data[data==555] <-  NA
data[data==666] <-  NA
data[data==888] <-  NA
data[data==999] <-  NA
data[data==9999] <- NA
data[data==999] <-  NA
data[data==777] <-  NA 
data[data=="99/99/9999"]<-NA
data[data=="NA"]<-NA

# 0.2 Check variables class
# Number of observations= # one childID per observation
nrow(data) == length(unique(data$childid))


# Check columns classes

# Check continuous variables
cont_var <- add_info[Type_var == "Continuous" & key_237_variable!="new",]$who_name
type_ga <- sapply(data[,..cont_var], class)
var_ga <- names(type_ga[type_ga=="character"]) # Continuous variables misclassified as character
data[, (var_ga) := lapply(.SD, as.numeric), .SDcols = var_ga]
cont_bound <- cont_bound(add_info,data)
biz_var <- setDT(cont_bound)[Consistent==FALSE,]$who_name  # variables outside the boundaries we set bizarre values as NA
for (var in biz_var){corrvar = correctbiz(var); data[,(var):=corrvar]} # correct bizarre values

# Create additional variables ----
data[,birth := ifelse(!is.na(birth_ga)|!is.na(ch_term),1,NA)] #birth indicator
data[,fet_death := ifelse(inducedabort==1,1,ifelse(birth==1,0,NA))]
data[,fet_death_ga := ifelse(!is.na(endga),endga,ifelse(!is.na(loss_ga),loss_ga,ifelse(!is.na(miscarriage_ga),miscarriage_ga,ifelse(!is.na(inducedabort_ga),inducedabort_ga,NA))))]
data[,birth := ifelse(is.na(birth)&fet_death==1,0,birth)]
data[,end_ga := ifelse(!is.na(endga),endga,ifelse(!is.na(birth_ga),birth_ga,fet_death_ga))] #clean end_ga 

# Microcephaly
data[, microcephaly_bin_fet := ifelse(!is.na(fet_micro),fet_micro,
                                   ifelse(!is.na(fet_micro_diag_tri)|fet_us_micro_tri1==1|fet_us_micro_tri2==1|fet_us_micro_tri3==1,1,NA))]
data[, microcephaly := ch_microcephaly]
data[, microcephaly_bin_birth := ifelse(microcephaly%in%c(0,3),0,ifelse(microcephaly%in%c(1,2),1,ch_microcephaly_bin))]

table(data$microcephaly_bin_birth,data$studyname,useNA = "always")


# Raw total data pre imputation
data_ori <- data
save(data_ori, file =here('3_Output_data','originaldata.RData')) 




