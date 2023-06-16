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

#data_origin <- as.data.table(import(here('Documents','GitHub','Zika_imputation','1_Input_data','zikv_033_datasets.dta')))
#add_info <- as.data.table(readxl::read_xlsx(here('Documents','GitHub','Zika_imputation','1_Input_data','MasterCodebook_October.xlsx'),sheet="237 key")) #CSV file with the
#study_info <- as.data.table(readxl::read_xlsx(here('Documents','GitHub','Zika_imputation','1_Input_data','MasterCodebook_October.xlsx'),sheet="StudyID")) #CSV file with the
#data<- merge(data_origin,study_info,by="file")
#source(here('Documents','GitHub','Zika_imputation','2.1_Code 33 studies','1.1_Pre_imputation_functions.R'))


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
data[data=="88/88/8888"]<-NA
data[data=="NA"]<-NA

# 0.2 Check variables class
# Number of observations= # one childID per observation
nrow(data) == length(unique(data$childid))
# 13991

# 0.3 Data as provided by the harmonization group without any modification
data_ori <- data
save(data_ori, file =here('3_Output_data','ori_data33.RData')) 


# Check columns classes

# Check continuous variables
cont_var <- add_info[Type_var == "Continuous" & key_237_variable!="new",]$who_name
type_ga <- sapply(data[,..cont_var], class)
var_ga <- names(type_ga[type_ga=="character"]) # Continuous variables misclassified as character
data[, (var_ga) := lapply(.SD, as.numeric), .SDcols = var_ga]
cont_bound <- cont_bound(add_info,data)
biz_var <- setDT(cont_bound)[Consistent==FALSE,]$who_name  # variables outside the boundaries we set bizarre values as NA
for (var in biz_var){corrvar = correctbiz(var); data[,(var):=corrvar]} # correct bizarre values



# Convert dates into gestational age (weeks) if it is possible
date_var<-add_info[Type_var == "Date",]$who_name
data[, (date_var) := lapply(.SD, as.Date, format="%d %b %Y"), .SDcols = date_var]
data[,conc_symp:=zikv_symp_date-zikv_symp_ga*7]
data[,conc_assay:=zikv_assay_date_1-zikv_assay_ga_1*7]
data[,conc_pcr:=zikv_pcr_date_1-zikv_pcr_ga_1*7]
data[,conc_elisa:=zikv_elisa_date_1-zikv_elisa_ga_1*7]
data[,conc_date:=apply(data[,c("conc_symp","conc_assay","conc_pcr","conc_elisa")], 1, min, na.rm = TRUE)] #get the min as date of conception

date_concep<-function(var_ga,var_date){
  min<-difftime(data[,get(var_date)],data[,conc_date],units="weeks")
  min<-ifelse(min>0&min<42,min,NA)
  value<-ifelse(is.na(data[,get(var_ga)]),min,data[,get(var_ga)])
}
data[, zikv_symp_ga:=ifelse(!is.na(zikv_symp_ga),zikv_symp_ga,ifelse(zikv_symp_tri==0,13,ifelse(zikv_symp_tri==1,27,ifelse(zikv_symp_tri==2,42,NA))))] # we aprox NA to the max boundary of the trimester
data[, zikv_symp_ga:=date_concep(var_ga="zikv_symp_ga",var_date="zikv_symp_date")] # we further approx with date
data[, zikv_assay_ga_1:=date_concep(var_ga="zikv_assay_ga_1",var_date="zikv_assay_date_1")]
data[, zikv_pcr_ga_1:=date_concep(var_ga="zikv_pcr_ga_1",var_date="zikv_pcr_date_1")]
data[, zikv_elisa_ga_1:=date_concep(var_ga="zikv_elisa_ga_1",var_date="zikv_elisa_date_1")]
data$date_t1_ga <- NA
data[, date_t1_ga:=date_concep(var_ga="date_t1_ga",var_date="date_t1")] #variable of first time visit on ga.
data[, zikv_ga:=ifelse(!is.na(zikv_ga),zikv_ga,ifelse(zikv_tri==0,13,ifelse(zikv_tri==1,27,ifelse(zikv_tri==2,42,NA))))] # we aprox NA to the max boundary of the trimester

# Create additional variables ----
# 1. Fet_death (fetus death variable) and fet_death_ga (time of fetus death) -----

# 1.1. Add miscarriage, loss, loss_etiology, birth in one variable fet_death and fet_death_ga

data[,maxbirth_ga:=fcase(ch_term==1,42,ch_term==2,28,ch_term==3,21,ch_term==4,33,ch_term==5,36,ch_term==6,44,default=NA)] #max ga according to ch_term
data[,birth_ga:=ifelse(is.na(birth_ga),maxbirth_ga,birth_ga)]

data[,birth:=ifelse(!is.na(birth_ga)|!is.na(ch_term)|!is.na(ch_vital_status),1,NA)] #birth indicator
data[,fet_death:=ifelse(inducedabort==1,1,ifelse(birth==1,0,NA))]
data[,fet_death_ga:=ifelse(!is.na(endga),endga,ifelse(!is.na(loss_ga),loss_ga,ifelse(!is.na(miscarriage_ga),miscarriage_ga,ifelse(!is.na(inducedabort_ga),inducedabort_ga,NA))))]
data[,end_ga := ifelse(!is.na(endga),endga,ifelse(!is.na(birth_ga),birth_ga,fet_death_ga))] #clean end_ga 
data[,end_ga :=ifelse(is.na(end_ga),maxbirth_ga,end_ga)] 



#1.2. Use etiology to set fet_death (et 0=live 1=miss 2=loss 3=imp fet_death) 
data[loss_etiology==3&loss==0,c("miscarriage","loss","loss_etiology","birth","fet_death"):= list(0,0,0,1,0)] # 156 study 6
data[loss_etiology==0&(is.na(fet_death)|fet_death==0),c("miscarriage","loss","loss_etiology","birth","fet_death"):= list(0,0,0,1,0)] #birth
data[loss_etiology==1&(is.na(fet_death)|fet_death==1),c("miscarriage","loss","loss_etiology","birth","fet_death"):= list(1,0,1,0,1)] #miscarriage
data[loss_etiology==2&(is.na(fet_death)|fet_death==1),c("miscarriage","loss","loss_etiology","birth","fet_death"):= list(0,0,2,0,1)] #abortion 
data[loss_etiology==3&(is.na(fet_death)|fet_death==1),c("miscarriage","loss","loss_etiology","birth","fet_death"):= list(0,1,3,0,1)] #loss
data[loss_etiology==4,c("miscarriage","loss","loss_etiology","birth","fet_death"):= list(0,0,4,0,1)] #stillbirth 
data[inducedabort==1,c("miscarriage","loss","loss_etiology","birth","fet_death"):= list(0,0,2,0,1)]

data[loss_etiology==1&loss==1&fet_death_ga<20,c("miscarriage","loss"):=list(1,0)]
data[loss_etiology==2&loss==1,c("miscarriage","loss","loss_etiology","birth","fet_death","inducedabort"):= list(0,1,3,0,1,0)] # 6 cases loss missassigned to abortion

data[,repabort:=ifelse(is.na(loss_etiology),0,ifelse(loss_etiology==2,1,0))]  # reported abort
data[ repabort==1,c("miscarriage","loss","loss_etiology","birth","fet_death"):= list(0,0,2,0,1)]
data[fet_death==1&fet_death_ga<20&repabort==0,c("miscarriage","loss","loss_etiology","birth","fet_death"):= list(1,0,1,0,1)]
data[fet_death==1&fet_death_ga>=20&repabort==0,c("miscarriage","loss","loss_etiology","birth","fet_death"):= list(0,1,3,0,1)]
data[loss==1&repabort==0&is.na(loss_etiology),c("miscarriage","loss","loss_etiology","birth","fet_death"):= list(0,1,3,0,1)]
data[miscarriage==1&repabort==0&is.na(loss_etiology),c("miscarriage","loss","loss_etiology","birth","fet_death"):= list(1,0,1,0,1)]
data[loss_etiology==1&birth==1,c("miscarriage","loss","loss_etiology","birth","fet_death"):= list(0,0,0,1,0)]
data[birth==1&is.na(loss_etiology),c("miscarriage","loss","loss_etiology","birth","fet_death"):= list(0,0,0,1,0)]
data[miscarriage==0&loss==0&is.na(birth)&repabort==0&ch_vital_status==0,c("miscarriage","loss","loss_etiology","birth","fet_death"):= list(0,0,0,1,0)]
data[!is.na(data$inf_term),c("miscarriage","loss","loss_etiology","birth","fet_death"):= list(0,0,0,1,0)]

# Check consistency across vital variables
table(losset=data$loss_etiology,loss=data$loss,useNA = "always")
table(losset=data$loss_etiology,miss=data$miscarriage,useNA = "always")
table(losset=data$loss_etiology,induceabort=data$inducedabort,useNA = "always")
table(losset=data$loss_etiology,birth=data$birth,useNA = "always")
table(losset=data$loss_etiology,cvital=data$ch_vital_status,useNA = "always")
table(losset=data$loss_etiology,rep_abort=data$repabort,useNA = "always")
table(cvital=data$ch_vital_status,rep_abort=data$repabort,useNA = "always")
table(cvital=data$ch_vital_status,loss=data$loss,useNA = "always")
table(cvital=data$ch_vital_status,birth=data$birth,useNA = "always") 


# 2. Microcephaly ----

# 2.0 Microcephaly in fetus
data[,microcephaly_bin_fet:=ifelse(!is.na(fet_micro),fet_micro,
                                   ifelse(!is.na(fet_micro_diag_tri)|fet_us_micro_tri1==1|fet_us_micro_tri2==1|fet_us_micro_tri3==1,1,NA))]
table(data$microcephaly_bin_fet, is.na(data$fet_micro_diag_ga),useNA = "always")
table(data$microcephaly_bin_fet, data$fet_micro,useNA = "always")
table(data$microcephaly,useNA = "always")

## 2.1. Microcephaly at birth given by study ----
data[, ch_microcephaly_bin:= ifelse(is.na(ch_microcephaly_bin)&ch_microcephaly%in%c(0,3),0,ifelse(is.na(ch_microcephaly_bin)&ch_microcephaly%in%c(1,2),1,ch_microcephaly_bin))]

data[, microcephaly_ga:=ifelse(!is.na(fet_micro_diag_ga),fet_micro_diag_ga,
                               ifelse(fet_micro_diag_tri==0,13,
                                      ifelse(fet_micro_diag_tri==1,27,
                                             ifelse(fet_micro_diag_tri==2,40,NA))))] 

# 2.2. Postnatal microcephaly ----

data <- micro_postnatal(data) # returns microcephaly_bin_postnatal variable
table(data$loss_etiology,data$microcephaly_bin_postnatal,useNA = "always")

data[, microcephaly_bin_postnatal:=ifelse(loss_etiology%in%c(1,2,3,4),2,microcephaly_bin_postnatal)]


# 3. Abnormalities ----

ncol<-c("fet_us_cns_tri2","fet_us_cns_tri3","ch_hydrocephaly","ch_corticalatrophy","ch_calcifications","ch_ventriculomegaly")
data[,neuroabnormality:=checkcon(data=data,setcol=ncol)]  # Neuroimaging abnormalities

ccol<-c("fet_us_msk_tri2","fet_us_msk_tri3")
data[,contractures:=checkcon(data=data,setcol=ccol)] # Congenital contractures

cacol<-c("fet_us_cardio_tri2","fet_us_cardio_tri3")
data[,cardioabnormality:=checkcon(data=data,setcol=cacol)] # Cardio abnormalities

gacol<-c("fet_us_gastro_tri2","fet_us_gastro_tri3")
data[,gastroabnormality:=checkcon(data=data,setcol=gacol)] # Gastrointestinal abnormalities

orcol<-c("fet_us_orofac_tri2","fet_us_orofac_tri3")
data[,oroabnormality:=checkcon(data=data,setcol=orcol)] # Orofacialintestinal abnormalities

eycol<-c("fet_us_eyeear_tri2","fet_us_eyeear_tri3")
data[,ocularabnormality:=checkcon(data=data,setcol=eycol)] # Ocular abnormalities

gecol<-c("fet_us_genur_tri2","fet_us_genur_tri3")
data[,genurabnormality:=checkcon(data=data,setcol=gecol)] # Genitourinaly abnormalities

noncol<-c("cardioabnormality","gastroabnormality","oroabnormality","genurabnormality")
data[,nonneurologic:=checkcon(data=data,setcol=noncol)] # Non neurological abnormalities

anycol<-c("neuroabnormality","contractures","cardioabnormality","gastroabnormality","oroabnormality","ocularabnormality","genurabnormality","ch_othabnorm","fet_us_bin_tri1","fet_us_bin_tri2","fet_us_bin_tri3")
data[,any_abnormality_czs:=checkcon(data=data,setcol=anycol)]  # Any congenital abnormality excluding microcephaly

gencol<-c("chromoabn_dx","ch_chromoabn")
data[,gen_anomalies:=checkcon(data=data,setcol=gencol)]  # Any congenital abnormality excluding microcephaly


# 4. Zika related test and load

#Zika test with evidence (zikv_test_ev) according to Ricardo Ximene's paper---
data_zik_test_ev <- ziktest_ml(data)  
data <- merge(data,data_zik_test_ev,by="childid",all.x = TRUE)
table( data$zikv_test_ev,data$zikv_preg,useNA = "always")
table( data$zikv_test_ev,useNA = "always")
#Limited Moderate Negative   Robust     <NA> 
#     1     2805      804     2859     7522 

# Include maternal PCR test in robust count
data[,zikv_test_ev_Ricardo_MPCR:=as.factor(ifelse(Maternal_PCR=="Yes","Robust",as.character(zikv_test_ev)))]
table(data$zikv_test_ev_Ricardo_MPCR,useNA = "always") # every cases turns to be robust
data[,zikv_test_ev:=as.factor(zikv_test_ev)]

#Viral load
data$zikv_pcr_vl_1<-as.numeric(data$zikv_pcr_vl_1)


# 5. Exposure to virus or pathogeneus----

#5.1.Concurrent or prior flavi- or alpha virus infection ----
data[,modificated1:=ifelse(arb_clindiag_plus==0,0,ifelse(!is.na(arb_clindiag_plus),1,0))]
data[,modificated2:=ifelse(arb_clindiag!=0&arb_clindiag!=1,0,ifelse(!is.na(arb_clindiag),1,0))] #arb_clindiag==777 to 6 on top
flcol<-c("modificated1","modificated2","denv_ever","chikv_ever")
data[,flavi_alpha_virus:=checkcon(data=data,setcol=flcol)]


#5.2. Intrauterine exposure to storch pathogens----
data[,modificated1:=ifelse(storch==0,0,ifelse(!is.na(storch),1,0))]
stcol<-c("modificated1","storch_bin","toxo","toxo_treat","syphilis","syphilis_treat","varicella","parvo","rubella","cmv","herpes","listeria","chlamydia","gonorrhea","genitalwarts")
data[,storch_patho:=checkcon(data=data,setcol=stcol)]


#5.3.Prior arb virus infection ----
data[,modificated1:=ifelse(zikv_pcr_everpos==1,1,ifelse(zikv_pcr_everpos==0,0,NA))] # 2 indeterminated
parcol<-c("modificated1","zikv_elisa_everpos","denv_ever","chikv_ever")
data[,arb_ever:=checkcon(data=data,setcol=parcol)]

#5.4.Current arb virus infection ----
data[,modificated1:=ifelse(zikv_pcr_res_1==1,1,ifelse(zikv_pcr_res_1==0,0,NA))] # 2 indeterminated
data[,modificated2:=ifelse(zikv_elisa_res_1==1,1,ifelse(zikv_elisa_res_1==0,0,NA))]
data[,modificated3:=ifelse(arb_clindiag==0,0,ifelse(!is.na(arb_clindiag),1,NA))] #arb_clindiag==777 to 6 modified at the beginning= other arbovirus
arcol<-c("modificated1","modificated2","modificated3","zikv_preg","denv_preg","chikv_preg")
data[,arb_preg:=checkcon(data=data,setcol=arcol)]

#5.5. Arb current pregnancy without consider zika
data[,modificated1:=ifelse(arb_clindiag==0|arb_clindiag==1,0,ifelse(!is.na(arb_clindiag),1,NA))] #arb_clindiag==777 to 6 modified at the beginning= other arbovirus
arcoln<-c("modificated1","denv_preg","chikv_preg")
data[,arb_preg_nz:=checkcon(data=data,setcol=arcoln)]

#5.6. Concurrent or prior dengue virus
decol<-c("chikv_preg","chikv_ever")
data[,denv_preg_ever:=checkcon(data=data,setcol=decol)]

#5.7. Concurrent or prior chik virus
chcol<-c("denv_preg","denv_ever")
data[,chikv_preg_ever:=checkcon(data=data,setcol=chcol)]

#5.8 Parity and gravidity
ncol<-c("prev_birthdef_bin","prev_eclampsia","prev_gestdiabet","prev_preeclampsia","prev_pregcomp_oth_bin","prev_pregcomp_oth_spec","prev_pregloss","prev_pretermlabor")
data[,prev_comp:=checkcon(data=data,setcol=ncol)] # previous complications
table(grav=data$gravidity,par=data$parity,useNA = "always")
data[,gravidity:=ifelse(gravidity<parity,parity,gravidity)]
data[,gravidity:=ifelse(is.na(gravidity)&prev_comp==1,1,gravidity)]

# 6. Exposure to drugs and vaccines --- 

#6.1. Maternal prescription drug use----
drcol<-c("med_bin","med_anticonvuls_bin","med_preg_bin","med_fertil_bin")
data[,drugs_prescr:=checkcon(data=data,setcol=drcol)]

#6.2. Maternal vaccination ----
vcol<-c("vac_rub","vac_vari","vac_yf")
data[,vaccination:=checkcon(data=data,setcol=vcol)]

#6.3. Maternal  teratogenic drug use
#Category X is teratogenic.
data[,drug_tera:=ifelse(is.na(med_oth),NA,
                        ifelse(med_oth%in%c("ortho-cyclen","norethindrone (Micronor) 0.35 mg tablet")|med_anticonvuls_bin==1,1,
                               ifelse(med_oth%in%c("Propiltiouracil","Neozine","Povidone-iodine 10% topical solution pyxis"),2,0)))]


#6.4 Preganancy comorbidities ----
corcol<-c("pregcomp_bin","gestdiab","eclampsia","preeclampsia")
data[,comorbid_preg:=checkcon(data=data,setcol=corcol)]


#7. Additional variables
#7.1. BMI
data[, bmi:= pre_pregweight/((height/100)^2)]
data[, bmi:= ifelse(bmi<0|bmi>50,NA,bmi)]

# Format factors and number variables
#data[,microcephaly_bin_birth:=as.numeric(microcephaly_bin_birth)]
#data[,microcephaly :=as.factor(microcephaly )]
data[,pretermlabor :=as.numeric(pretermlabor )]
data[,alcohol :=as.numeric(alcohol)]
data[,ch_czs :=as.numeric(ch_czs)]
data[,microcephaly_bin_postnatal:=as.factor(microcephaly_bin_postnatal)]
data[,microcephaly_bin_fet:=as.numeric(microcephaly_bin_fet)]
data[,miscarriage:=as.numeric(miscarriage)]
data[,ch_microcephaly_bin:=as.numeric(ch_microcephaly_bin)]
data[,ch_microcephaly:=as.factor(ch_microcephaly)]
data[,rash:=as.factor(rash)]
data[,fever:=as.factor(fever)]
data[,cc_hiv:=as.factor(cc_hiv)]
data[,arthralgia:=as.factor(arthralgia)]
data[,arb_symp:=as.factor(arb_symp)]
data[,multiplegest:=as.numeric(multiplegest)]
data[,ch_sex:=as.numeric(ch_sex)]
data[,miscarriage:=as.numeric(miscarriage)]
data[,repabort:=as.numeric(repabort)]

#8  Missing data Plot

var_incl <- add_info[Essential=="yes",]$who_name
var_incl<-var_incl[!var_incl%in% c( "childid","childid_original","fetid_original","fetid","mid","mid_original")]
dataf<-data[,..var_incl]
totalval<-as.data.table(table(dataf$studyname))
colnames(totalval)<-c("studyname","N")
dmatrix<-dataf[, lapply(.SD, function(x) sum(is.na(x))/.N), studyname] #matrix of % of missingness
dmatrix2<-as.data.table(melt(dmatrix,id.vars="studyname"))
dmatrix2<-merge(dmatrix2,totalval,by="studyname")

dmatrix2[,name:=paste0(studyname,"\nN=",N)]
dmatrix2[,miss:=round(value*100,1)]
dmatrix2[,text:=paste0("study: ", studyname, "\n", "variable: ", variable, "\n", "miss%: ",miss)]
dmatrix2[,tmiss:=ifelse(miss==100,"Systematical","Sporadical")]

p<-ggplot(dmatrix2, aes(x=name, y=variable,fill=miss,text=text)) +
  geom_tile() +
  scale_fill_gradientn(colours=c("green","yellow","red")) +
  theme(axis.text.x = element_text(angle = 45,vjust=0.5,size=8),
        axis.text.y = element_text(size=8))+xlab("Study name")+ylab("Variable")+
  labs(fill='%missing') 

ggplotly(p, tooltip="text")

# 9. Flux plot ----
fx<-flux(dataf)
fx$var<-rownames(fx)
fx<-as.data.table(fx)
fx[,ofluxinc:=ifelse(outflux>=0.5,"yes","no")]
fx<-fx[order(-outflux)]
fx[,orderimp:=1:nrow(fx)]
fx<-fx[,c("var","outflux","ofluxinc","orderimp")]

# install java and rJava to use xlsx package to make an interfase between R and xlxs file
library(xlsx)
path<- here('1_Input_data','MasterCodebook_October.xlsx')
if("Outflux" %in% names(getSheets(loadWorkbook(path)))){
  wb <- loadWorkbook(path)
  removeSheet(wb, sheetName = "Outflux")
  saveWorkbook(wb, path)
}

write.xlsx(fx, path, sheetName = "Outflux", append = TRUE,col.names = TRUE,row.names =FALSE) # update Mastercode.xlxs with outflux results.

# 11. Final selected variables ----

# Please refer to the MasterCodebook_October.xlsx on the folder 1_Input_data
add_info <- as.data.table(readxl::read_xlsx(here('1_Input_data','MasterCodebook_October.xlsx'),sheet="237 key")) #CSV file with the
add_infoi <- add_info[order(Orderimp1)]
var_imp <- add_infoi[Imp_obj1=="yes"]$who_name
finc <- study_info[Included==1]$file #included studies

# Raw total data pre imputation

data_det <- data
save(data_det, file =here('3_Output_data','det_data33.RData')) 
#save(fdata, file =here('Documents','GitHub','Zika_imputation','3_Output_data','rawfinaldata33.RData'))


# Data for imputation
fdata <- data[file%in%finc,..var_imp]
save(fdata, file =here('3_Output_data','finaldata33.RData')) 
#save(fdata, file =here('Documents','GitHub','Zika_imputation','3_Output_data','finaldata33.RData'))


