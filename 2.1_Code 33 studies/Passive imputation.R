
# Data manipulation packages
library(here) 
library(data.table)
library(dplyr)
library(rio)
#Functions

#Function returns 1 if there is any anomaly detected and 0 if there were no anomaly detected across the recorded variables , NA no recorded variables.
checkcon<-function(data,col1){ #
  data[,(col1):= lapply(.SD, as.numeric), .SDcols = col1]
  data$anyT <- rowSums(data[, .SD, .SDcols = col1], na.rm=T)
  data$allNA <- rowSums(!is.na(data[, .SD, .SDcols = col1]))
  data$Final<- ifelse(data$allNA==0,NA,ifelse(data$anyT>0,1,0))
  return(data$Final)
}


#data <- as.data.table(import(here('Documents','Julius','ZIKV analyses','2. Data','zikv_033_datasets.dta'))) #For Anneke - please leave this line in :-)
#data <- as.data.table(import(here('1_Input_data','zikv_033_datasets.dta')))
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
library(growthstandards)

intragrowthzscores <- function(data,birth_ga){
  dataset.zscore = subset(data,is.na(ch_sex) == F&ch_sex != 999)
  lencm2zscore=  igb_lencm2zscore(gagebrth = dataset.zscore$birth_ga*7, lencm=dataset.zscore$ch_length ,
                                  sex = ifelse( dataset.zscore$ch_sex== 0, "Male",ifelse(dataset.zscore$ch_sex== 1, "Female",NA)))  
  wtkg2zscore = igb_wtkg2zscore(gagebrth = dataset.zscore$birth_ga*7, wtkg=dataset.zscore$ch_weight/1000, 
                                sex = ifelse( dataset.zscore$ch_sex== 0, "Male",ifelse(dataset.zscore$ch_sex== 1, "Female",NA)))
  hcircm2zscore = igb_hcircm2zscore(gagebrth = dataset.zscore$birth_ga*7, hcircm=dataset.zscore$ch_head_circ_birth, 
                                    sex = ifelse( dataset.zscore$ch_sex== 0, "Male",ifelse(dataset.zscore$ch_sex== 1, "Female",NA))) 
  list(lencm2zscore,wtkg2zscore,hcircm2zscore)
}
lencm2zscore = intragrowthzscores(data)[[1]]
wtkg2zscore = intragrowthzscores(data)[[2]]
hcircm2zscore = intragrowthzscores(data)[[3]]

data$microcephaly_bin2[is.na(data$ch_sex) ==F] <- ifelse(hcircm2zscore < (-2),1,
                                                          ifelse(hcircm2zscore >= -2 & hcircm2zscore != 999, 0,NA))

data$fet_micro_new<-NA
data$fet_micro_new<-ifelse(data$fet_us_micro_tri1==1 | data$fet_us_micro_tri2==1 | data$fet_us_micro_tri3==1 |
                                data$fet_micro==1 | data$microcephaly_bin2==1,1,NA)
data$temp<-NA
data$temp<-ifelse(data$fet_us_micro_tri1==0 | data$fet_us_micro_tri2==0 | data$fet_us_micro_tri3==0 |
                    data$fet_micro==0 | data$microcephaly_bin2==0,0,NA)
data$fet_micro_new<-ifelse(is.na(data$fet_micro_new) & data$temp==0,0,data$fet_micro_new)


#czs new
#child microcephaly


#Postnatal microcephaly

#Neuroimaging abnormalities other than microcephaly
data$neuroabnormality<-NA
data$neuroabnormality<-ifelse(data$fet_us_cns_tri2==1 | data$fet_us_cns_tri3==1 | data$ch_hydrocephaly==1 |
                                data$ch_corticalatrophy==1 | data$ch_calcifications==1 | data$ch_ventriculomegaly==1,1,NA)
data$temp<-NA
data$temp<-ifelse(data$fet_us_cns_tri2==0 | data$fet_us_cns_tri3==0 | data$ch_hydrocephaly==0 |
                    data$ch_corticalatrophy==0 | data$ch_calcifications==0 | data$ch_ventriculomegaly==0,0,NA)
data$neuroabnormality<-ifelse(is.na(data$neuroabnormality) & data$temp==0,0,data$neuroabnormality)

#Ocular abnormalities
data$ocularabnormality<-NA
data$ocularabnormality<-ifelse(data$fet_us_eyeear_tri2==1 | data$fet_us_eyeear_tri3==1,1,NA)
data$temp<-NA
data$temp<-ifelse(data$fet_us_eyeear_tri2==0 | data$fet_us_eyeear_tri3==0,0,NA)
data$ocularabnormality<-ifelse(is.na(data$ocularabnormality) & data$temp==0,0,data$ocularabnormality)

#Congenital contractures
data$contractures<-NA
data$contractures<-ifelse(data$fet_us_msk_tri2==1 | data$fet_us_msk_tri3==1,1,NA)
data$temp<-NA
data$temp<-ifelse(data$fet_us_msk_tri2==0 | data$fet_us_msk_tri3==0,0,NA)
data$contractures<-ifelse(is.na(data$contractures) & data$temp==0,0,data$contractures)

#Other nonneurologic abnormalities
data$nonneurologic<-NA
data$nonneurologic<-ifelse(data$fet_us_cardio_tri2==1 | data$fet_us_gastro_tri2==1 | data$fet_us_orofac_tri2==1 |
                             data$fet_us_genur_tri2==1 | data$fet_us_cardio_tri3==1 | data$fet_us_gastro_tri3==1 |
                             data$fet_us_orofac_tri3==1 | data$fet_us_genur_tri3==1,1,NA)
data$temp<-NA
data$temp<-ifelse(data$fet_us_cardio_tri2==0 | data$fet_us_gastro_tri2==0 | data$fet_us_orofac_tri2==0 |
                    data$fet_us_genur_tri2==0 | data$fet_us_cardio_tri3==0 | data$fet_us_gastro_tri3==0 |
                    data$fet_us_orofac_tri3==0 | data$fet_us_genur_tri3==0,0,NA)
data$nonneurologic<-ifelse(is.na(data$nonneurologic) & data$temp==0,0,data$nonneurologic)

#BMI
data$bmi<-data$pre_pregweight/((data$height/100)**2)

#Maternal teratogenic drug use

#Maternal vaccination
data$vaccination<-NA
data$vaccination<-ifelse(data$vac_rub==1 | data$vac_vari==1 | data$vac_yf==1,1,NA)
data$temp<-NA
data$temp<-ifelse(data$vac_rub==0 | data$vac_vari==0 | data$vac_yf==0,0,NA)
data$vaccination<-ifelse(is.na(data$vaccination) & data$temp==0,0,data$vaccination)

#Genetic anomalies
data$gen_anomalies<-NA
data$gen_anomalies<-ifelse(data$chromoabn_dx==1 | data$ch_chromoabn==1,1,NA)
data$temp<-NA
data$temp<-ifelse(data$chromoabn_dx==0 | data$ch_chromoabn==0,0,NA)
data$gen_anomalies<-ifelse(is.na(data$gen_anomalies) & data$temp==0,0,data$gen_anomalies)

#Viral load
data$zikv_pcr_vl_1<-as.numeric(data$zikv_pcr_vl_1)

#Concurrent or prior dengue virus
data$denv_preg_ever<-NA
data$denv_preg_ever<-ifelse(data$denv_preg==1 | data$denv_ever==1,1,NA)
data$temp<-NA
data$temp<-ifelse(data$denv_preg==0 | data$denv_ever==0,0,NA)
data$denv_preg_ever<-ifelse(is.na(data$denv_preg_ever) & data$temp==0,0,data$denv_preg_ever)

#Concurrent or prior chik virus
data$chikv_preg_ever<-NA
data$chikv_preg_ever<-ifelse(data$chikv_preg==1 | data$chikv_ever==1,1,NA)
data$temp<-NA
data$temp<-ifelse(data$chikv_preg==0 | data$chikv_ever==0,0,NA)
data$chikv_preg_ever<-ifelse(is.na(data$chikv_preg_ever) & data$temp==0,0,data$chikv_preg_ever)

#Pregnancy related comorbidities
data$comorbid_preg<-NA
data$comorbid_preg<-ifelse(data$pregcomp_bin==1 | data$gestdiab==1 | data$eclampsia==1 |
                             data$preeclampsia==1,1,NA)
data$temp<-NA
data$temp<-ifelse(data$pregcomp_bin==0 | data$gestdiab==0 | data$eclampsia==0 |
                    data$preeclampsia==0,0,NA)
data$comorbid_preg<-ifelse(is.na(data$comorbid_preg) & data$temp==0,0,data$comorbid_preg)

#Variables that we finally need:
#Exposures: zikv_preg, fet_zikv, zikv_test_ev, zikv_ga
#Primary outcomes: bdeath, bdeath_ga, fet_micro, ch_czs, czsn
#Secondary fetal outcomes: igr_curr_prg
#Secondary infant outcomes: ch_microcephaly_bin, ch_weight, ch_craniofac_abn_bin,
#neuroabnormality, ocularabnormality, contractures, nonneurologic
#Confounders: age, educ, maritalstat,ethnicity, bmi, ses, tobacco, drugs_bin, alcohol,
#vaccination, gen_anomalies, zikv_pcr_vl_1, denv_preg_ever, chikv_preg_ever,
#comorbid_bin, comorbid_preg, storch_bin, arb_symp, fever, rash, arthralgia, headache,
#muscle_pain, arthritis, vomiting, abd_pain, bleed, fatigue, sorethroat


data2<-subset(data, select=c(studycode,
                             zikv_preg,fet_zikv, zikv_ga, ch_czs,igr_curr_prg, ch_weight, ch_craniofac_abn_bin,
                             neuroabnormality, ocularabnormality, contractures, nonneurologic,
                             age, educ, maritalstat,ethnicity, bmi, ses, tobacco, drugs_bin, alcohol,
                             vaccination, gen_anomalies, zikv_pcr_vl_1, denv_preg_ever, chikv_preg_ever,
                             comorbid_bin, comorbid_preg, storch_bin, arb_symp, fever, rash, arthralgia, headache,
                             muscle_pain, arthritis, vomiting, abd_pain, bleed, fatigue, sorethroat))
#Change to numeric
data2$studycode<-as.factor(data2$studycode)
data2<-data2 %>% mutate_if(is.character,as.numeric)
#sapply(data2, class)

#Not yet included in subset: zikv_test_ev, bdeath, bdeath_ga,fet_micro, czsn, ch_microcephaly_bin

summary(data2)