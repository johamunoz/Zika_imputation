
#Starts directly after loading imputed dataset
m<-max(data$.imp)
data$studyname<-data$studyimp
studynames<-c("001-BRA","002-BRA","003-GUF","004-ESP","005-ESP","006-COL","007-COL",
              "008-USA","009-GRD","010-BRA","011-BRA","012-TTO","013-BRA","014-BRA",
              "015-BRA","016-HND","017-USA","018-COL","019-BRA","020-BRA","021-PRI",
              "022-BRA","023-BRA","024-GTM","025-BRA","026-KEN","027-BRA")
data$studyname_fac<-factor(data$studyimp)
levels(data$studyname_fac)<-studynames

#Drop variables we do not need
data<-subset(data,select=-c(multiplegest,gravidity,ch_head_circ_birth,cc_hiv,pretermlabor,
                            zikv_igm_res_1,syphilis,zikv_elisa_everpos, zikv_assay_ga_1,zikv_elisa_res_1,
                            drugs_prescr, zikv_pcr_everpos,med_bin))
data$bdeath<-1-data$birth

#Create factors
data$arb_symp<-as.factor(data$arb_symp)
data$birth<-as.factor(data$birth)
data$zikv_preg<-as.factor(data$zikv_preg)
data$microcephaly<-as.factor(data$microcephaly)
data$microcephaly_bin_birth<-as.factor(data$microcephaly_bin_birth)
data$microcephaly_bin_postnatal<-as.factor(data$microcephaly_bin_postnatal)
data$storch_patho<-as.factor(data$storch_patho)
data$rash<-as.factor(data$rash)
data$fever<-as.factor(data$fever)
data$arthralgia<-as.factor(data$arthralgia)
data$comorbid_preg<-as.factor(data$comorbid_preg)
data$alcohol<-as.factor(data$alcohol)
data$headache<-as.factor(data$headache)
data$czs<-as.factor(data$czs)
data$any_abnormality_czs<-as.factor(data$any_abnormality_czs)
data$educ<-as.factor(data$educ)
data$zikv_test_ev<-as.factor(data$zikv_test_ev)
data$neuroabnormality<-as.factor(data$neuroabnormality)
data$drugs_bin<-as.factor(data$drugs_bin)
data$arthritis<-as.factor(data$arthritis)
data$chikv_preg_ever<-as.factor(data$chikv_preg_ever)
data$denv_preg_ever<-as.factor(data$denv_preg_ever)
data$maritalstat<-as.factor(data$maritalstat)
data$comorbid_bin<-as.factor(data$comorbid_bin)
data$fet_micro<-as.factor(data$fet_micro)
data$nonneurologic<-as.factor(data$nonneurologic)
data$bdeath<-as.factor(data$bdeath)

#Create dichotomous outcome variables to calculate incidence
#data$microcephaly_bin <- droplevels(data$microcephaly_bin)
#miscarriage (<20 weeks gestation)
data$miscarriage<-0
data$miscarriage[data$bdeath==1 & data$end_ga<20]<-1
data$miscarriage<-as.factor(data$miscarriage)
#Fetal loss (>=20 weeks gestation)
data$loss<-0
data$loss[data$bdeath==1 & data$end_ga>=20]<-1
data$loss<-as.factor(data$loss)
#Early fetal death (20-27 weeks gestation)
data$efdeath<-0
data$efdeath[data$bdeath==1 & data$end_ga>=20 & data$end_ga<28]<-1
data$efdeath<-as.factor(data$efdeath)
#Late fetal death (after 28 weeks gestation)
data$lfdeath<-0
data$lfdeath[data$bdeath==1 & data$end_ga>=28]<-1
data$lfdeath<-as.factor(data$lfdeath)
#Late fetal death (after 28 weeks gestation) with microcephaly
#data$lfdeath_micro<-0
#data$lfdeath_micro[data$lfdeath==1 & data$microcephaly_bin_birth==1]<-1
#data$lfdeath_micro<-as.factor(data$lfdeath)

#zikv_ga should be NA if there is zikv_preg is 0
data$zikv_ga[data$zikv_preg==0]<-NA

#Trimester of zika infection
data$zikv_tri<-NA
data$zikv_tri[data$zikv_ga<=12]<-1
data$zikv_tri[(data$zikv_ga>12 & data$zikv_ga<=27)]<-2
data$zikv_tri[data$zikv_ga>27]<-3
data$zikv_tri<-as.factor(data$zikv_tri)


