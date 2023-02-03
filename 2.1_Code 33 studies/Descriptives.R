
#Check https://cran.r-project.org/web/packages/table1/vignettes/table1-examples.html for examples

library(dplyr)
library(table1)
library(xlsx)


#Functions
rndr <- function(x, ...) {
  if (is.factor(x) || is.character(x)) {
    c(render.default(x, ...), c(`Overall N`=sum(!is.na(x))))
  } else {
    render.default(x, ...)
  }
}

#Data
data.noimp<-read.csv("/Users/jdamen/Library/CloudStorage/OneDrive-UMCUtrecht/Research/WHO ZIKA/2. Data/20221027 zikv_not_imputed.csv",header=T)

data2<-subset(data.noimp, select=c(studycode,birth_ga,
                             zikv_preg,fet_zikv, zikv_ga, ch_czs,igr_curr_prg, ch_microcephaly, ch_weight, ch_craniofac_abn_bin,
                             ocularabnormality, nonneurologic,
                             age, educ, maritalstat,ethnicity, bmi, ses, tobacco, drugs_bin, alcohol, drug_tera,
                             zikv_pcr_vl_1, denv_preg_ever, chikv_preg_ever,
                             comorbid_bin, storch_bin, arb_symp, fever, rash, arthralgia, headache,
                             muscle_pain, arthritis, vomiting, abd_pain, bleed, fatigue, sorethroat,
                             maxbirth_ga,birth,fet_death,fet_death_ga,end_ga,hcircm2zscore,microcephaly_hc,
                             microcephaly,microcephaly_bin_birth,microcephaly_ga,microcephaly_bin_postnatal,
                             neuroabnormality,contractures,cardioabnormality,gastroabnormality,oroabnormality,
                             genurabnormality,any_abnormality_czs,fet_micro,
                             gen_anomalies,zikv_test_ev,czs,flavi_alpha_virus,storch_patho,arb_ever,arb_preg,
                             arb_preg_nz,drugs_prescr,vaccination,comorbid_preg))
#Exclude studies with selection bias
data2<-data2[data2$studycode!="001-BRA" & data2$studycode!="019-BRA" & data2$studycode!="022-BRA" & data2$studycode!="023-BRA",]

rm(data.noimp)
data2$maritalstat[data2$maritalstat==4]<-NA
data2$ses[data2$ses==3]<-NA

#Create outcome variables
data2$bdeath<-1-data2$birth
#miscarriage (<20 weeks gestation)
data2$miscarriage<-NA
data2$miscarriage[!is.na(data2$bdeath)]<-0
data2$miscarriage[data2$bdeath==1 & data2$end_ga<20]<-1
data2$miscarriage<-as.factor(data2$miscarriage)
#Fetal loss (>=20 weeks gestation)
data2$loss<-NA
data2$loss[!is.na(data2$bdeath)]<-0
data2$loss[data2$bdeath==1 & data2$end_ga>=20]<-1
data2$loss<-as.factor(data2$loss)
#Early fetal death (20-27 weeks gestation)
data2$efdeath<-NA
data2$efdeath[!is.na(data2$bdeath)]<-0
data2$efdeath[data2$bdeath==1 & data2$end_ga>=20 & data2$end_ga<28]<-1
data2$efdeath<-as.factor(data2$efdeath)
#Late fetal death (after 28 weeks gestation)
data2$lfdeath<-NA
data2$lfdeath[!is.na(data2$bdeath)]<-0
data2$lfdeath[data2$bdeath==1 & data2$end_ga>=28]<-1
data2$lfdeath<-as.factor(data2$lfdeath)
#Late fetal death (after 28 weeks gestation) with microcephaly
data2$lfdeath_micro<-NA
data2$lfdeath_micro[!is.na(data2$bdeath) & !is.na(data2$microcephaly_bin_birth)]<-0
data2$lfdeath_micro[data2$lfdeath==1 & data2$microcephaly_bin_birth==1]<-1
data2$lfdeath_micro<-as.factor(data2$lfdeath_micro)

data2<-data2 %>% mutate_if(is.character,as.factor)
#Exposures
data2$zikv_preg<-factor(as.factor(data2$zikv_preg), levels=c(1,0),labels=c("Positive","Negative"))
data2$fet_zikv<-factor(as.factor(data2$fet_zikv), levels=c(1,0),labels=c("Positive","Negative"))
data2$zikv_test_ev <- factor(data2$zikv_test_ev, levels = c("Negative", "Limited", "Moderate","Robust"))

#Outcomes
data2$miscarriage <- factor(data2$miscarriage, levels=c(1,0),labels=c("Yes","No"))
data2$loss <- factor(data2$loss, levels=c(1,0),labels=c("Yes","No"))
data2$microcephaly_bin_birth<-factor(as.factor(data2$microcephaly_bin_birth), levels=c(1,0),labels=c("Yes","No"))
data2$czs<-factor(as.factor(data2$czs), levels=c(1,0),labels=c("Yes","No"))
data2$efdeath <- factor(data2$efdeath, levels=c(1,0),labels=c("Yes","No"))
data2$lfdeath <- factor(data2$lfdeath, levels=c(1,0),labels=c("Yes","No"))
data2$lfdeath_micro <- factor(data2$lfdeath_micro, levels=c(1,0),labels=c("Yes","No"))
data2$igr_curr_prg<-factor(as.factor(data2$igr_curr_prg), levels=c(1,0),labels=c("Yes","No"))
data2$microcephaly_bin_postnatal<-factor(as.factor(data2$microcephaly_bin_postnatal), levels=c(1,0),labels=c("Yes","No"))
data2$ch_craniofac_abn_bin<-factor(as.factor(data2$ch_craniofac_abn_bin), levels=c(1,0),labels=c("Yes","No"))
data2$neuroabnormality<-factor(as.factor(data2$neuroabnormality), levels=c(1,0),labels=c("Yes","No"))
data2$ocularabnormality<-factor(as.factor(data2$ocularabnormality), levels=c(1,0),labels=c("Yes","No"))
data2$contractures<-factor(as.factor(data2$contractures), levels=c(1,0),labels=c("Yes","No"))
data2$nonneurologic<-factor(as.factor(data2$nonneurologic), levels=c(1,0),labels=c("Yes","No"))
data2$any_abnormality_czs<-factor(as.factor(data2$any_abnormality_czs), levels=c(1,0),labels=c("Yes","No"))
#Other outcomes not in table
data2$birth<-as.factor(data2$birth)
data2$fet_death<-as.factor(data2$fet_death)
data2$fet_micro<-as.factor(data2$fet_micro)
data2$microcephaly_hc<-as.factor(data2$microcephaly_hc)
data2$cardioabnormality<-as.factor(data2$cardioabnormality)
data2$gastroabnormality<-as.factor(data2$gastroabnormality)
data2$oroabnormality<-as.factor(data2$oroabnormality)
data2$genurabnormality<-as.factor(data2$genurabnormality)

#Covariates
data2$educ<-factor(as.factor(data2$educ),levels=c(0,1,2,3,4,5),labels=c("No education","Primary school","Secondary school","Some college","Bachelor's degree","Graduate or Professional degree"))
data2$maritalstat<-factor(as.factor(data2$maritalstat), levels=c(0,1,2,3),labels=c("Single","Married/Living as married/Cohabitating","Divorced/Separated","Widowed"))
data2$ethnicity<-factor(as.factor(data2$ethnicity), levels=c(0,1,2,3,4,5),labels=c("Caucasian descent","African descent","East Asian descent","South Asian descent","Indigenous descent","Mixed"))
data2$ses<-factor(as.factor(data2$ses), levels=c(0,1,2),labels=c("Low","Medium","High"))
data2$tobacco<-factor(as.factor(data2$tobacco), levels=c(1,0),labels=c("Yes","No"))
data2$drugs_bin<-factor(as.factor(data2$drugs_bin), levels=c(1,0),labels=c("Yes","No"))
data2$alcohol<-factor(as.factor(data2$alcohol), levels=c(1,0),labels=c("Yes","No"))
data2$drug_tera<-factor(as.factor(data2$drug_tera), levels=c(1,2,0),labels=c("Teratogenic","Risk of teratogenic","Not teratogenic"))
data2$vaccination<-factor(as.factor(data2$vaccination), levels=c(1,0),labels=c("Yes","No"))
data2$gen_anomalies<-factor(as.factor(data2$gen_anomalies), levels=c(1,0),labels=c("Yes","No"))
data2$denv_preg_ever<-factor(as.factor(data2$denv_preg_ever), levels=c(1,0),labels=c("Yes","No"))
data2$chikv_preg_ever<-factor(as.factor(data2$chikv_preg_ever), levels=c(1,0),labels=c("Yes","No"))
data2$comorbid_bin<-factor(as.factor(data2$comorbid_bin), levels=c(1,0),labels=c("Yes","No"))
data2$comorbid_preg<-factor(as.factor(data2$comorbid_preg), levels=c(1,0),labels=c("Yes","No"))
data2$storch_patho<-factor(as.factor(data2$storch_patho), levels=c(1,0),labels=c("Yes","No"))
data2$arb_symp<-factor(as.factor(data2$arb_symp), levels=c(1,0),labels=c("Yes","No"))
data2$fever<-factor(as.factor(data2$fever), levels=c(1,0),labels=c("Yes","No"))
data2$rash<-factor(as.factor(data2$rash), levels=c(1,0),labels=c("Yes","No"))
data2$arthralgia<-factor(as.factor(data2$arthralgia), levels=c(1,0),labels=c("Yes","No"))
data2$headache<-factor(as.factor(data2$headache), levels=c(1,0),labels=c("Yes","No"))
data2$muscle_pain<-factor(as.factor(data2$muscle_pain), levels=c(1,0),labels=c("Yes","No"))
data2$arthritis<-factor(as.factor(data2$arthritis), levels=c(1,0),labels=c("Yes","No"))
data2$vomiting<-factor(as.factor(data2$vomiting), levels=c(1,0),labels=c("Yes","No"))
data2$abd_pain<-factor(as.factor(data2$abd_pain), levels=c(1,0),labels=c("Yes","No"))
data2$bleed<-factor(as.factor(data2$bleed), levels=c(1,0),labels=c("Yes","No"))
data2$fatigue<-factor(as.factor(data2$fatigue), levels=c(1,0),labels=c("Yes","No"))
data2$sorethroat<-factor(as.factor(data2$sorethroat), levels=c(1,0),labels=c("Yes","No"))
#Other covariates not in table
data2$storch_bin<-as.factor(data2$storch_bin)
data2$flavi_alpha_virus<-as.factor(data2$flavi_alpha_virus)
data2$arb_ever<-as.factor(data2$arb_ever)
data2$arb_preg<-as.factor(data2$arb_preg)
data2$arb_preg_nz<-as.factor(data2$arb_preg_nz)
data2$drugs_prescr<-as.factor(data2$drugs_prescr)

#Exclude studies that are only part of sensitivity analyses
#data2<-data2[data2$studycode!="002-BRA" & data2$studycode!="008-USA" & data2$studycode!="011-BRA" & data2$studycode!="013-BRA" & data2$studycode!="018-COL",]

#Exposures
label(data2$zikv_preg)<-"Maternal zika - study definition"
label(data2$fet_zikv)<-"Fetal zika"
label(data2$zikv_test_ev)<-"Maternal zika - Ximenes definition"

#Covariates
label(data2$age)<-"Age (years)"
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
label(data2$birth_ga)<-"Gestational age at birth (weeks)"
label(data2$zikv_ga)<-"Gestational age at zika infection (weeks)"
label(data2$zikv_pcr_vl_1)<-"Viral load for PCR (copies/µL)"
label(data2$denv_preg_ever)<-"Concurrent or prior Dengue virus infection"
label(data2$chikv_preg_ever)<-"Concurrent or prior Chikungunya virus infection"
label(data2$comorbid_bin)<-"Comorbidities before pregnancy"
label(data2$comorbid_preg)<-"Pregnancy-related comorbidities"
label(data2$storch_patho)<-"STORCH pathogen infection"
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

#Outcomes
label(data2$miscarriage)<-"Miscarriage (<20 weeks gestation)"
label(data2$loss)<-"Fetal loss (>=20 weeks gestation)"
label(data2$efdeath)<-"Early fetal death (20-27 weeks gestation)"
label(data2$lfdeath)<-"Late fetal death (>=28 weeks gestation)"
label(data2$lfdeath_micro)<-"Late fetal death with microcephaly"
label(data2$fet_death)<-"Fetus died during pregnancy"
label(data2$fet_death_ga)<-"Gestational age the fetus died (weeks)"
label(data2$fet_micro)<-"Fetal microcephaly"
label(data2$ch_czs)<-"Congenital zika syndrome - study definition"
label(data2$czs)<-"Congenital zika syndrome - WHO definition"
label(data2$igr_curr_prg)<-"Intrauterine growth restriction"
label(data2$microcephaly)<-"Microcephaly at birth"
label(data2$microcephaly_bin_postnatal)<-"Microcephaly after birth"
label(data2$birth)<-"Birth"
label(data2$birth_ga)<-"Gestational age at birth (weeks)"
label(data2$ch_weight)<-"Birth weight (gram)"
label(data2$ch_craniofac_abn_bin)<-"Craniofacial disproportion"
label(data2$neuroabnormality)<-"Neuroimaging abnormality"
label(data2$contractures)<-"Congenital contractures"
label(data2$cardioabnormality)<-"Cardiovascular abnormality"
label(data2$gastroabnormality)<-"Gastrointestinal abnormality"
label(data2$oroabnormality)<-"Oro-facial abnormality"
label(data2$ocularabnormality)<-"Ocular abnormality, congenital deafness or hearing loss"
label(data2$genurabnormality)<-"Genitourinary system abnormality"
label(data2$nonneurologic)<-"Any non-neurologic abnormality"
label(data2$any_abnormality_czs)<-"Any congenital abnormality"
label(data2$end_ga)<-"Gestational age at which the baby was born or died (weeks)"
label(data2$microcephaly_bin_birth)<-"Microcephaly"


#gato<-as.data.table(data.noimp)[,.(numchild=length(unique(childid))) ,by=.(studycode,mid_original)]
#gato[,.(max=max(numchild)),by=.(studycode)]

#mytable.exp<-table1(~ zikv_preg + fet_zikv + zikv_test_ev,data=data2,excel=1)
#mytable.cov<-table1(~ age + educ + maritalstat + ethnicity + bmi + ses + tobacco + drugs_bin + alcohol + drug_tera + vaccination + 
#                      gen_anomalies + birth_ga + zikv_ga + zikv_pcr_vl_1 + denv_preg_ever + chikv_preg_ever + comorbid_bin + 
#                      comorbid_preg + storch_patho + arb_symp + fever + rash + arthralgia + headache + muscle_pain + arthritis + 
#                      vomiting + abd_pain + bleed + fatigue + sorethroat | studycode,data=data2,excel=1)
#write.xlsx(mytable.exp,"/Users/jdamen/Documents/Julius/ZIKV analyses/4. Resultaten/Table 1 exposures.xlsx")
#write.xlsx(mytable.cov,"/Users/jdamen/Documents/Julius/ZIKV analyses/4. Resultaten/Table 1 covariates.xlsx")
#write.xlsx(mytable.out,"/Users/jdamen/Documents/Julius/ZIKV analyses/4. Resultaten/Table 1 outcomes.xlsx")


#Tables: cmd+a -> paste in Word or paste in Excel
#Exposures
#Overall
table1(~ zikv_preg + fet_zikv + zikv_test_ev,data=data2,excel=1)
#Per study
table1(~ zikv_preg + fet_zikv + zikv_test_ev | studycode,data=data2,excel=1)

#Outcomes
#Overall
table1(~ miscarriage + loss + microcephaly_bin_birth + czs + efdeath + lfdeath + lfdeath_micro + igr_curr_prg + microcephaly_bin_postnatal + 
         end_ga + ch_weight + ch_craniofac_abn_bin + neuroabnormality + ocularabnormality + contractures + nonneurologic + 
         any_abnormality_czs,data=data2,excel=1)
#Exclude missings
table1(~ miscarriage + loss + microcephaly_bin_birth + czs + efdeath + lfdeath + lfdeath_micro + igr_curr_prg + microcephaly_bin_postnatal + 
         end_ga + ch_weight + ch_craniofac_abn_bin + neuroabnormality + ocularabnormality + contractures + nonneurologic + 
         any_abnormality_czs | studycode,data=data2,excel=1,
       render.missing = NULL,render.categorical = "FREQ (PCTnoNA%)",render = rndr)
#Per study
table1(~ miscarriage + loss + microcephaly_bin_birth + czs + efdeath + lfdeath + lfdeath_micro + igr_curr_prg + microcephaly_bin_postnatal + 
         end_ga + ch_weight + ch_craniofac_abn_bin + neuroabnormality + ocularabnormality + contractures + nonneurologic + 
         any_abnormality_czs | studycode,data=data2,excel=1)
#By ZIKA status study definition
table1(~ miscarriage + loss + microcephaly_bin_birth + czs + efdeath + lfdeath + lfdeath_micro + igr_curr_prg + microcephaly_bin_postnatal + 
         end_ga + ch_weight + ch_craniofac_abn_bin + neuroabnormality + ocularabnormality + contractures + nonneurologic + 
         any_abnormality_czs | zikv_preg,data=data2,excel=1)
#By ZIKA status Ricardo definition
table1(~ miscarriage + loss + microcephaly_bin_birth + czs + efdeath + lfdeath + lfdeath_micro + igr_curr_prg + microcephaly_bin_postnatal + 
         end_ga + ch_weight + ch_craniofac_abn_bin + neuroabnormality + ocularabnormality + contractures + nonneurologic + 
         any_abnormality_czs | zikv_test_ev,data=data2,excel=1)

#Covariates
#Overall
table1(~ age + educ + maritalstat + ethnicity + bmi + ses + tobacco + drugs_bin + alcohol + drug_tera + vaccination + 
         gen_anomalies + end_ga + zikv_ga + zikv_pcr_vl_1 + denv_preg_ever + chikv_preg_ever + comorbid_bin + 
         comorbid_preg + storch_patho + arb_symp + fever + rash + arthralgia + headache + muscle_pain + arthritis + 
         vomiting + abd_pain + bleed + fatigue + sorethroat,data=data2,excel=1)
#Per study
table1(~ age + educ + maritalstat + ethnicity + bmi + ses + tobacco + drugs_bin + alcohol + drug_tera + vaccination + 
         gen_anomalies + end_ga + zikv_ga + zikv_pcr_vl_1 + denv_preg_ever + chikv_preg_ever + comorbid_bin + 
         comorbid_preg + storch_patho + arb_symp + fever + rash + arthralgia + headache + muscle_pain + arthritis + 
         vomiting + abd_pain + bleed + fatigue + sorethroat | studycode,data=data2,excel=1)



#Original dataset for comparison
library(rio)
library(here) 
data_origin <- import("/Users/jdamen/Documents/Julius/ZIKV analyses/2. Data/zikv_033_datasets.dta")
study_info <- readxl::read_xlsx("/Users/jdamen/Documents/GitHub/Zika_imputation/1_Input_data/MasterCodebook_October.xlsx",sheet="StudyID")
data_origin<- merge(data_origin,study_info,by="file")
table(as.factor(data_origin$ch_microcephaly_bin),as.factor(data_origin$studycode))


#Explore imputed data
setwd("/Users/jdamen/Library/CloudStorage/OneDrive-UMCUtrecht/Research/WHO ZIKA/2. Data")
data<-read.csv("20230202 zikv_imputed.csv",header=T)
setwd("/Users/jdamen/Documents/GitHub/Zika_imputation/2.1_Code 33 studies")
source("Functions Objective 1.R")
source("ZIKV prep v33.R")

data$.imp2<-as.factor(data$.imp)
table1(~ zikv_preg + zikv_test_ev| .imp2,data=data ,excel=1)
table1(~ miscarriage + loss + microcephaly_bin_birth + czs + efdeath + lfdeath + microcephaly_bin_postnatal + 
         end_ga + ch_weight + neuroabnormality + nonneurologic + any_abnormality_czs| .imp2,data=data ,excel=1)
table1(~ age + educ + maritalstat + bmi + drugs_bin + alcohol + end_ga + zikv_ga + denv_preg_ever + chikv_preg_ever + comorbid_bin + 
         comorbid_preg + storch_patho + arb_symp + fever + rash + arthralgia + headache + arthritis | .imp2,data=data ,excel=1)
