rm(list=ls()) # clean environment

# Load packages ---
# Data manipulation package
  library(data.table) 
  library(here)  # define folder paths
  library(rio)   # import STATA file

# Field specific package
  library(growthstandards)


# Load dataset and dependencies ----
  data <- as.data.table(import(here('1_Input_data','zikv_033_datasets.dta')))
  source(here('2.1_Code 33 studies','1.1_Pre_imputation_functions.R'))
  
  data <- as.data.table(import(here('Documents','Julius','ZIKV analyses','2. Data','zikv_033_datasets.dta'))) 
  ###RUN PREIMPUTATION FUNCTIONS SCRIPT


  
# 0. Initial checks ----
  
# 0.1 Missing data observations----
  data[data==""] <- NA
  data[data==555] <-  NA
  data[data==666] <-  NA
  data[data==888] <-  NA
  data[data==999] <-  NA
  data[data==9999] <-  NA
  data[data==999] <-  NA
# TODO  @Anneke could you please check if there is any valuable information on 777 among your selected variables :) 
  #@Johanna: I checked it and fortunately there are none that we currently use, so we can leave the line in!
  data[data==777] <-  NA 

# 0.2 Check variables class
# Number of observations
  nrow(data) #13992
  length(unique(data$childid))# one childID per observation

# Check columns classes
#TODO @Anneke Can we include on the required list the variables classes so we can check all them here?  
  #@Johanna: yes of course! Can you let me know which variable class you prefer for e.g. categorical variables (with 0 and 1)? Do you want them as numeric?
  # Check ga variables 
  var_ga<-grep(pattern="._ga",x=names(data),value=TRUE)
  type_ga<-sapply(data[,..var_ga], class)
  var_ga<-c("miscarriage_ga","inducedabort_ga","loss_ga","fet_zikv_ga","fet_micro_diag_ga")
  data[, (var_ga) := lapply(.SD, as.numeric), .SDcols = var_ga]


# Study code 
  data[, studycode:= fcase(file=="001","001-BRA",file=="002","002-BRA",file=="003","003-GUF",file=="004","004-ESP",
                           file=="005","005-ESP",file=="006","006-COL",file=="007","007-COL",file=="008","008-USA",
                           file=="009","009-GRD",file=="010","010-BRA",file=="011","011-BRA",file=="012","012-TTO",
                           file=="013","013-BRA",file=="014","014-BRA",file=="015","015-BRA",file=="016","016-HND",
                           file=="017","017-USA",file=="018","018-COL",file=="019","019-BRA",file=="020","020-BRA",
                           file=="021","021-PRI",file=="022","022-BRA",file=="023","023-BRA",file=="024","024-BRA",
                           file=="025","025-BRA",file=="026","026-BRA",file=="027","027-BRA")]


  
# 1. Fet_death (fetus death variable) and fet_death_ga (time of fetus death) -----

# 1.1. Add miscarriage, loss, loss_etiology, birth in one variable fet_death and fet_death_ga
  
# TODO @ Anneke given we did not have a lot of end dates I used for NA cases the maximum ga of the specified ch_term
  #Johanna, that is ok, but I think for the first category, better to take 40 weeks instead of 42. I changed that.
  data[,maxbirth_ga:=fcase(ch_term==1,40,ch_term==2,28,ch_term==3,21,ch_term==4,33,ch_term==5,36,ch_term==6,44,default=NA)] #max ga according to ch_term
  data[,birth_ga:=ifelse(is.na(birth_ga),maxbirth_ga,birth_ga)]

# TODO @ Anneke following Mabel's new notation i defined bdeath as fet_death.. @Johanna, fine!!
  data[,birth:=ifelse(!is.na(birth_ga)|!is.na(ch_term),1,NA)] #birth indicator
  data[,fet_death:=ifelse(inducedabort==1,1,ifelse(birth==1,0,NA))]
  data[,fet_death_ga:=ifelse(!is.na(endga),endga,ifelse(!is.na(loss_ga),loss_ga,ifelse(!is.na(miscarriage_ga),miscarriage_ga,ifelse(!is.na(inducedabort_ga),inducedabort_ga,NA))))]
  data[,end_ga := ifelse(!is.na(endga),endga,ifelse(!is.na(birth_ga),birth_ga,fet_death_ga))] #clean end_ga 
  data[,end_ga :=ifelse(is.na(end_ga),maxbirth_ga,end_ga)] 

#1.2. Use etiology to set fet_death (et 0=live 1=miss 2=loss 3=imp fet_death) 
  data[loss_etiology==0&(is.na(fet_death)|fet_death==0),c("miscarriage","loss","loss_etiology","birth","fet_death"):= list(0,0,0,1,0)] #birth
  data[loss_etiology==1&(is.na(fet_death)|fet_death==1),c("miscarriage","loss","loss_etiology","birth","fet_death"):= list(1,0,1,0,1)] #miscarriage
  data[loss_etiology==2&(is.na(fet_death)|fet_death==1),c("miscarriage","loss","loss_etiology","birth","fet_death"):= list(0,1,2,0,1)] #loss
  data[loss_etiology==3&(is.na(fet_death)|fet_death==1),c("miscarriage","loss","loss_etiology","birth","fet_death"):= list(0,1,2,0,1)] #class as loss
  data[fet_death==1&fet_death_ga<20,c("miscarriage","loss","loss_etiology","birth","fet_death"):= list(1,0,1,0,1)]
  data[fet_death==1&fet_death_ga>=20,c("miscarriage","loss","loss_etiology","birth","fet_death"):= list(0,1,2,0,1)]
  data[loss==1&is.na(loss_etiology),c("miscarriage","loss","loss_etiology","birth","fet_death"):= list(0,1,2,0,1)]
  data[birth==1&loss_etiology==1&ch_vital_status==0,c("miscarriage","loss","loss_etiology","birth","fet_death"):= list(0,0,0,1,0)]
  data[birth==1&is.na(loss_etiology),c("miscarriage","loss","loss_etiology","birth","fet_death"):= list(0,0,0,1,0)]
  data[miscarriage==0&loss==0&is.na(birth)&ch_vital_status==0,c("miscarriage","loss","loss_etiology","birth","fet_death"):= list(0,0,0,1,0)]
  data[!is.na(data$inf_term),c("miscarriage","loss","loss_etiology","birth","fet_death"):= list(0,0,0,1,0)]
  data[inducedabort==1&birth==1,c("miscarriage","loss","loss_etiology","birth","fet_death"):= list(0,1,2,0,1)]
  data[birth==0&ch_vital_status==0,c("miscarriage","loss","loss_etiology","birth","fet_death"):= list(0,0,0,1,0)]

#1.3. Check whether is coherent with inf_alive_birth 0=alive, induce abort 0=No
#TODO @Anneke Could you please ask what is loss_etiology=4?
  #@Johanna it is "Stillbirth or Intrapartum death (death during labor)" 
  checktable<-data.table(table(miss=data$miscarriage,loss=data$loss,et=data$loss_etiology,birth=data$birth,fet_death=data$fet_death,vst=data$ch_vital_status,iabo=data$inducedabort,inft=!is.na(data$ch_term), useNA="always")) #V1 is miscarriage, V2 is loss and N the number of observations
  checktable[N!=0,]



# 2. Microcephaly ----
  
# 2.1. Microcephaly just the moment fetus baby is out!! (microcephaly,microcephaly_bin, microcephayly_ga)
  
  #igb_hcircm2zscore : function (gagebrth, hcircm, sex = "Female")  # birth measurements 

  data[!is.na(ch_sex), hcircm2zscore:=as.numeric(igb_hcircm2zscore(gagebrth = end_ga*7, hcircm=ch_head_circ_birth,sex=ifelse(ch_sex== 0, "Male","Female")))]  
  data[, microcephaly_hc  := ifelse(hcircm2zscore<=-3,2,ifelse(hcircm2zscore<=-2,1,ifelse(hcircm2zscore<=2,0,ifelse(!is.na(hcircm2zscore),3,NA))))] # given by formula
  data[, microcephaly := ifelse(!is.na(microcephaly_hc),microcephaly_hc,ch_microcephaly)] # ch_microcephaly given by hospital so we priorized the result given by circunference
  data[, microcephaly_bin:=ifelse(microcephaly%in%c(0,3),0,ifelse(microcephaly%in%c(1,2),1,ch_microcephaly_bin))]
# TODO @Anneke ask about fet_micro_diag_ga definition (diagnosis!=assessment) and inconsistences with fet_micro_diag_tri I mean can we used us_diagnosis term to aprox the micro_diag_ga?
  data[,  microcephaly_ga:=ifelse(!is.na(fet_micro_diag_ga),fet_micro_diag_ga,
                                ifelse(fet_micro_diag_tri==0,13,
                                ifelse(fet_micro_diag_tri==1,27,
                                ifelse(fet_micro_diag_tri==2,40,NA))))] 


# 2.2. Postnatal microcephaly ----
  data <- micro_postnatal(data) # returns microcephaly_bin_postnatal variable
  table(post=data$microcephaly_bin_postnatal,pre=data$microcephaly_bin)  
  
# 3. Abnormalities ----
  
  ncol<-c("fet_us_cns_tri2","fet_us_cns_tri3","ch_hydrocephaly","ch_corticalatrophy","ch_calcifications","ch_ventriculomegaly")
  data[,neuro_abnormality:=checkcon(data=data,setcol=ncol)]  # Neuroimaging abnormalities
  
  ccol<-c("fet_us_msk_tri2","fet_us_msk_tri3")
  data[,contractures:=checkcon(data=data,setcol=ccol)] # Congenital contractures
  
  cacol<-c("fet_us_cardio_tri2","fet_us_cardio_tri3")
  data[,cardio_abnormality:=checkcon(data=data,setcol=cacol)] # Cardio abnormalities
  
  gacol<-c("fet_us_gastro_tri2","fet_us_gastro_tri3")
  data[,gastro_abnormality:=checkcon(data=data,setcol=gacol)] # Gastrointestinal abnormalities
  
  orcol<-c("fet_us_orofac_tri2","fet_us_orofac_tri3")
  data[,oro_abnormality:=checkcon(data=data,setcol=orcol)] # Orofacialintestinal abnormalities
  
  eycol<-c("fet_us_eyeear_tri2","fet_us_eyeear_tri3")
  data[,ocular_abnormality:=checkcon(data=data,setcol=eycol)] # Ocular abnormalities
  
  gecol<-c("fet_us_genur_tri2","fet_us_genur_tri3")
  data[,genur_abnormality:=checkcon(data=data,setcol=gecol)] # Genitourinaly abnormalities
 
  noncol<-c("cardio_abnormality","gastro_abnormality","oro_abnormality","genur_abnormality")
  data[,nonneuro_abnormality:=checkcon(data=data,setcol=noncol)] # Non neurological abnormalities
 
  anycol<-c("neuro_abnormality","contractures","cardio_abnormality","gastro_abnormality","oro_abnormality","ocular_abnormality","genur_abnormality","ch_othabnorm","fet_us_bin_tri1","fet_us_bin_tri2","fet_us_bin_tri3")
  data[,any_abnormality_czs:=checkcon(data=data,setcol=anycol)]  # Any congenital abnormality excluding microcephaly

  gencol<-c("chromoabn_dx","ch_chromoabn")
  data[,gen_abnormality:=checkcon(data=data,setcol=gencol)]  # Any congenital abnormality excluding microcephaly
  

# 4. Zika related test and load
  
  #Zika test with evidence (zikv_test_ev) according to Ricardo paper---
  data_zik_test_ev <- ziktest_ml(data)  
  
  #Viral load
  data$zikv_pcr_vl_1<-as.numeric(data$zikv_pcr_vl_1)
  

# 5. CZS variable according to WHO definition ----
#TODO@ Anneke how we will calcualte the CZS i mean we will use the zikv_test_ev variable somewhere? 
  #@Johanna YES! I have changed the following lines so we are using it.
  
  #WHO definition for CZS: Presence of confirmed maternal or fetal ZIKV infection AND presence of severe microcephaly at birth AND presence of other malformations (including limb contractures, high muscle tone, eye abnormalities, and hearing loss, nose etc.)
  #Johanna, note that the WHO definition is changed:
  #WHO definition for CZS: Presence of confirmed maternal or fetal ZIKV infection AND (presence of severe microcephaly at birth OR presence of other malformations (including limb contractures, high muscle tone, eye abnormalities, and hearing loss, nose etc.))
  data[,czs2:=ifelse((data$zikv_test_ev=="Robust" | data$zikv_test_ev=="Moderate" | data$fet_zikv==1) & ((data$microcephaly==2) | (data$anyabnormality_czs==1)),1,
                     ifelse(data$zikv_test_ev=="Negative"&data$fet_zikv==0 & data$microcephaly!=2&data$anyabnormality_czs==0,0,NA))] 
  data[,czsn:=ifelse(is.na(ch_czs),czs2,ch_czs)]  

  
      

# 6. Exposure to virus or pathogeneus----
  
  #6.1.Concurrent or prior flavi- or alpha virus infection ----
  data[,modificated1:=ifelse(arb_clindiag_plus==0,0,ifelse(!is.na(arb_clindiag_plus),1,0))]
  data[,modificated2:=ifelse(arb_clindiag!=0&arb_clindiag!=1,0,ifelse(!is.na(arb_clindiag),1,0))] #arb_clindiag==777 to 6 on top
  flcol<-c("modificated1","modificated2","denv_ever","chikv_ever")
  data[,flavi_alpha_virus:=checkcon(data=data,setcol=flcol)]
  
  
  #6.2. Intrauterine exposure to storch pathogens----
  data[,modificated1:=ifelse(storch==0,0,ifelse(!is.na(storch),1,0))]
  stcol<-c("modificated1","storch_bin","toxo","toxo_treat","syphilis","syphilis_treat","varicella","parvo","rubella","cmv","herpes","listeria","chlamydia","gonorrhea","genitalwarts")
  data[,storch_patho:=checkcon(data=data,setcol=stcol)]
  
  
  #6.3.Prior arb virus infection ----
  data[,modificated1:=ifelse(zikv_pcr_everpos==1,1,ifelse(zikv_pcr_everpos==0,0,NA))] # 2 indeterminated
  parcol<-c("modificated1","zikv_elisa_everpos","denv_ever","chikv_ever")
  data[,arb_ever:=checkcon(data=data,setcol=parcol)]
  
  #6.4.Current arb virus infection ----
  data[,modificated1:=ifelse(zikv_pcr_res_1==1,1,ifelse(zikv_pcr_res_1==0,0,NA))] # 2 indeterminated
  data[,modificated2:=ifelse(zikv_elisa_res_1==1,1,ifelse(zikv_elisa_res_1==0,0,NA))]
  data[,modificated3:=ifelse(arb_clindiag==0,0,ifelse(!is.na(arb_clindiag),1,NA))] #arb_clindiag==777 to 6 modified at the beginning= other arbovirus
  arcol<-c("modificated1","modificated2","modificated3","zikv_preg","denv_preg","chikv_preg")
  data[,arb_preg:=checkcon(data=data,setcol=arcol)]
  
  #6.5. Arb current pregnancy without consider zika
  data[,modificated1:=ifelse(arb_clindiag==0|arb_clindiag==1,0,ifelse(!is.na(arb_clindiag),1,NA))] #arb_clindiag==777 to 6 modified at the beginning= other arbovirus
  arcoln<-c("modificated1","denv_preg","chikv_preg")
  data[,arb_preg_nz:=checkcon(data=data,setcol=arcoln)]
  
  #6.6. Concurrent or prior dengue virus
  decol<-c("chikv_preg","chikv_ever")
  data[,denv_preg_ever:=checkcon(data=data,setcol=decol)]
  
  #6.7. Concurrent or prior chik virus
  chcol<-c("denv_preg","denv_ever")
  data[,chikv_preg_ever:=checkcon(data=data,setcol=chcol)]
  

  
# 7. Exposure to drugs and vaccines --- 
  
  #7.1. Maternal prescription drug use----
  drcol<-c("med_bin","med_anticonvuls_bin","med_preg_bin","med_fertil_bin")
  data[,drugs_prescr:=checkcon(data=data,setcol=drcol)]
  
  #7.2. Maternal vaccination ----
  vcol<-c("vac_rub","vac_vari","vac_yf")
  data[,vaccination:=checkcon(data=data,setcol=vcol)]

  #7.3. Maternal  teratogenic drug use
  #Category X is teratogenic.
#TODO @ Anneke no idea what is this,, but thanks God i work with you :), please check it.  
  data[,drug_tera:=ifelse(is.na(med_oth),NA,
                    ifelse(med_oth%in%c("ortho-cyclen","norethindrone (Micronor) 0.35 mg tablet"),1,0))]
  
  #@Johanna can you add another category (number 2) for the following values: "Propiltiouracil" "Neozine" "Povidone-iodine 10% topical solution pyxis"

# 8. Preganancy comorbidities ----
  corcol<-c("pregcomp_bin","gestdiab","eclampsia","preeclampsia")
  data[,comorbid_preg:=checkcon(data=data,setcol=corcol)]
  


