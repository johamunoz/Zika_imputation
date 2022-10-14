
#Check https://cran.r-project.org/web/packages/table1/vignettes/table1-examples.html for examples

#First run 1_Pre_imputation.R
library(dplyr)
library(table1)
library(xlsx)

data2<-subset(data, select=c(studycode,birth_ga,
                             zikv_preg,fet_zikv, zikv_ga, ch_czs,igr_curr_prg, ch_microcephaly, ch_weight, ch_craniofac_abn_bin,
                             ocularabnormality, nonneurologic,
                             age, educ, maritalstat,ethnicity, bmi, ses, tobacco, drugs_bin, alcohol, drug_tera,
                             zikv_pcr_vl_1, denv_preg_ever, chikv_preg_ever,
                             comorbid_bin, storch_bin, arb_symp, fever, rash, arthralgia, headache,
                             muscle_pain, arthritis, vomiting, abd_pain, bleed, fatigue, sorethroat,
                             maxbirth_ga,birth,fet_death,fet_death_ga,end_ga,hcircm2zscore,microcephaly_hc,
                             microcephaly,microcephaly_bin_birth,microcephaly_ga,microcephaly_bin_postnatal,
                             neuroabnormality,contractures,cardioabnormality,gastroabnormality,oroabnormality,
                             genurabnormality,any_abnormality_czs,
                             gen_anomalies,zikv_test_ev,czs,flavi_alpha_virus,storch_patho,arb_ever,arb_preg,
                             arb_preg_nz,drugs_prescr,vaccination,comorbid_preg))


data2<-data2 %>% mutate_if(is.character,as.factor)
data2$zikv_preg<-as.factor(data2$zikv_preg)
data2$ocularabnormality<-as.factor(data2$ocularabnormality)
data2$nonneurologic<-as.factor(data2$nonneurologic)
data2$tobacco<-as.factor(data2$tobacco)
data2$drug_tera<-as.factor(data2$drug_tera)
data2$denv_preg_ever<-as.factor(data2$denv_preg_ever)
data2$chikv_preg_ever<-as.factor(data2$chikv_preg_ever)
data2$storch_bin<-as.factor(data2$storch_bin)
data2$birth<-as.factor(data2$birth)
data2$fet_death<-as.factor(data2$fet_death)
data2$microcephaly_hc<-as.factor(data2$microcephaly_hc)
data2$microcephaly_bin_postnatal<-as.factor(data2$microcephaly_bin_postnatal)
data2$neuroabnormality<-as.factor(data2$neuroabnormality)
data2$contractures<-as.factor(data2$contractures)
data2$cardioabnormality<-as.factor(data2$cardioabnormality)
data2$gastroabnormality<-as.factor(data2$gastroabnormality)
data2$oroabnormality<-as.factor(data2$oroabnormality)
data2$genurabnormality<-as.factor(data2$genurabnormality)
data2$any_abnormality_czs<-as.factor(data2$any_abnormality_czs)
data2$gen_anomalies<-as.factor(data2$gen_anomalies)
data2$flavi_alpha_virus<-as.factor(data2$flavi_alpha_virus)
data2$storch_patho<-as.factor(data2$storch_patho)
data2$arb_ever<-as.factor(data2$arb_ever)
data2$arb_preg<-as.factor(data2$arb_preg)
data2$arb_preg_nz<-as.factor(data2$arb_preg_nz)
data2$drugs_prescr<-as.factor(data2$drugs_prescr)
data2$vaccination<-as.factor(data2$vaccination)
data2$comorbid_preg<-as.factor(data2$comorbid_preg)


#data2<-data2 %>% mutate_if(is.character,as.numeric)
#sapply(data2, class)

#Not yet included in subset: zikv_test_ev, bdeath, bdeath_ga,fet_micro, czsn, ch_microcephaly_bin

#

#Descriptives categorical variables
#data2<-data2 %>% mutate_if(is.numeric,as.factor)

#Drop continuous variables
#data2<-subset(data2,select=-c(zikv_ga,ch_weight,age,bmi,zikv_pcr_vl_1))

summary(data2)

as.data.frame(cbind(summary(data2[,2]),summary(data2[,2])/length(data2$studycode),colnames(data2[,2])))

#Exposures
label(data2$zikv_preg)<-"Maternal zika - study definition"
label(data2$fet_zikv)<-"Fetal zika"
label(data2$zikv_test_ev)<-"Maternal zika - Ximenes definition"

#Covariates
label(data2$age)<-"Age"
label(data2$educ)<-"Education"
label(data2$maritalstat)<-"Marital status"
label(data2$ethnicity)<-"Ethnicity"
label(data2$bmi)<-"Body mass index"
label(data2$ses)<-"Socioeconomic status"
label(data2$tobacco)<-"Smoking"
label(data2$drugs_bin)<-"Drug use"
label(data2$alcohol)<-"Alcohol use"
label(data2$drug_tera)<-"Teratogenic drug use"
label(data2$vaccination)<-"Maternal vaccination"
label(data2$gen_anomalies)<-"Genetic anomalies"
label(data2$birth_ga)<-"Gestational age at birth"
label(data2$zikv_ga)<-"Gestational age at zika infection"
label(data2$zikv_pcr_vl_1)<-"Viral load for PCR"
label(data2$denv_preg_ever)<-"Concurrent or prior denv"
label(data2$chikv_preg_ever)<-"Concurrent or prior chikv"
label(data2$comorbid_bin)<-"Comorbidities before pregnancy"
label(data2$comorbid_preg)<-"Pregnancy-related comorbidities"
label(data2$storch_bin)<-"STORCH pathogen infection"
label(data2$arb_symp)<-"Arbovirus-related symptoms"
label(data2$fever)<-"Fever"
label(data2$rash)<-"Rash"
label(data2$arthralgia)<-"Arthralgia"
label(data2$headache)<-"Headache"
label(data2$muscle_pain)<-"Muscle pain"
label(data2$arthritis)<-"Arthritis"
label(data2$vomiting)<-"Vomiting"
label(data2$abd_pain)<-"Abdominal pain"
label(data2$bleed)<-"Bleeding"
label(data2$fatigue)<-"Fatigue"
label(data2$sorethroat)<-"Sore throat"
#label(data2$)<-""


mytable.exp<-table1(~ zikv_preg + fet_zikv + zikv_test_ev,data=data2,excel=1)
mytable.cov<-table1(~ age + educ + maritalstat + ethnicity + bmi + ses + tobacco + drugs_bin + alcohol + drug_tera + vaccination + 
                      gen_anomalies + birth_ga + zikv_ga + zikv_pcr_vl_1 + denv_preg_ever + chikv_preg_ever + comorbid_bin + 
                      comorbid_preg + storch_bin + arb_symp + fever + rash + arthralgia + headache + muscle_pain + arthritis + 
                      vomiting + abd_pain + bleed + fatigue + sorethroat | studycode,data=data2,excel=1)
#write.xlsx(mytable.exp,"/Users/jdamen/Documents/Julius/ZIKV analyses/4. Resultaten/Table 1 exposures.xlsx")
#write.xlsx(mytable.cov,"/Users/jdamen/Documents/Julius/ZIKV analyses/4. Resultaten/Table 1 covariates.xlsx")



