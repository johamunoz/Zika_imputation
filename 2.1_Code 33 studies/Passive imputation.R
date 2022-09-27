
# Data manipulation packages
library(here) 
library(data.table)
library(dplyr)

#Functions

#Function returns 1 if there is any anomaly detected and 0 if there were no anomaly detected across the recorded variables , NA no recorded variables.
checkcon<-function(data,col1){ #
  data<-as.numeric(data)
  data$anyT <- rowSums(data[, .SD, .SDcols = col1], na.rm=T)
  data$allNA <- rowSums(!is.na(data[, .SD, .SDcols = col1]))
  data$Final<- ifelse(data$allNA==0,NA,ifelse(data$anyT>0,1,0))
  return(data$Final)
}


#data <- as.data.table(import(here('Documents','Julius','ZIKV analyses','2. Data','zikv_033_datasets.dta'))) #For Anneke - please leave this line in :-)
data <- as.data.table(import(here('1_Input_data','zikv_033_datasets.dta')))
data[,studycode:=fcase(file=="001","001-BRA",
                       file=="002","002-BRA",
                       file=="003","003-GUF",
                       file=="004","004-ESP",
                       file=="005","005-ESP",
                       file=="006","006-COL",
                       file=="007","007-COL",
                       file=="008","008-USA",
                       file=="009","009-GRD",
                       file=="010","010-BRA",
                       file=="011","011-BRA",
                       file=="012","012-TTO",
                       file=="013","013-BRA",
                       file=="014","014-BRA",
                       file=="015","015-BRA",
                       file=="016","016-HND",
                       file=="017","017-USA",
                       file=="018","018-COL",
                       file=="019","019-BRA",
                       file=="020","020-BRA",
                       file=="021","021-PRI",
                       file=="022","022-BRA",
                       file=="023","023-BRA",
                       file=="024","024-BRA",
                       file=="025","025-BRA",
                       file=="026","026-BRA",
                       file=="027","027-BRA"
)]

#1.1. Set to NA : 666,777,888,999 values----
data[data==""] <-NA
data[data==555] <-  NA
data[data==666] <-  NA
data[data==888] <-  NA
data[data==999] <-  NA
data[data==9999] <-  NA
data[data==999] <-  NA
data[data==777] <-  NA #@Johanna, not sure if we need to keep this line. Maybe better to specify in which variables we would like to recode 777 to missing
data[data=="NA"] <-NA 


#zikv status according to Ricardo
#bdeath and bdeath_ga
#fetal microcephaly
#czs new
#child microcephaly

#Neuroimaging abnormalities other than microcephaly
col1<-c("fet_us_cns_tri2","fet_us_cns_tri3","ch_hydrocephaly","ch_corticalatrophy","ch_calcifications","ch_ventriculomegaly")
data[,neuroabnormality:=checkcon(data=data,col1=col1)]





#Variables that we finally need:
#Exposures: zikv_preg, fet_zikv, zikv_test_ev
#Primary outcomes: bdeath, bdeath_ga, fet_micro, ch_czs, czsn
#Secondary fetal outcomes: igr_curr_preg
#Secondary infant outcomes: ch_microcephaly_bin, ch_weight, ch_craniofac_abn_bin