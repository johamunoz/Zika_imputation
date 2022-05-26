
#Passive imputation data provided 09-09-21

rm(list=ls())

# Data manipulation packages
library(rstudioapi) 
library(data.table)
library(dplyr)

# Imputation packages
library(mitml)
library(mice)
library(micemd)
library(miceadds)

# Graphic packages
library(corrplot)
library(ggplot2)
library(ggtext)
library(plotly)

# Field specific package
library(growthstandards)

#0. Load information ----

base.dir <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path)) # set main base working directory
setwd(base.dir)

data<-as.data.table(read.csv('1_Input_data/pilot10_08SEP21_withmetadata.csv', stringsAsFactors=FALSE, fileEncoding="latin1"))

infoexp<-as.data.table(readxl::read_xlsx("1_Input_data/Infoexp.xlsx",sheet="Table")) #Table were are specified the included variables according Expert opinion, also the includes the order in which variables are imputed 
var_inc<-infoexp[Inclusion==1,Variable] #Variables to work with
data<-data[,colnames(data)%in%var_inc,with = FALSE] #Filter dataset

#1. Initial clean up the dataset -----
#1.0. Check for duplicates ----
data[,comb:=paste0(studyname,mid_original,sep="_")] #create a variable that combine study and patient id
n_occur <- data.table(table(data$comb))
n_occur[n_occur$N>1] #check visually which observations are duplicated
var1<-n_occur[n_occur$N > 1&V1!='_Spain_Soriano'& V1!='_TrinidadTobago_Sohan']$V1 
var1<-n_occur[n_occur$N > 1]$V1
Duplicates<-data[comb%in%var1&multiplegest!=1] # table to report to harmonization team.
write.csv(Duplicates,file='5_Internal_support/Possible_duplicates.csv')

# 1.1. Check format of categorical variables----
char <-   which(sapply(data, is.character)) #which are char
datac<-data[,..char]
summary(datac)
unique(datac$inf_head_circ_age_fu1)
data[,inf_head_circ_age_fu1:=as.numeric(sub(" .*", "",data$inf_head_circ_age_fu1))] # remove years
data[,zikv_elisa_ga_1:=as.numeric(sub(" .*", "",data$zikv_elisa_ga_1))] #remove weeks

#unique(data$othabnorm_spec) # 39 categories ask if really want to include this
#data[,othabnorm_spec:=NULL] # We temporal remove it from dataset as it is messy
data[zikv_pcr_date_1%in%c("22027","888",""), zikv_pcr_date_1:=NA]
data[zikv_elisa_date_1%in%c(""),zikv_elisa_date_1:=NA]
data[symp_date%in%c(""),symp_date:=NA]
data[,arb_clindiag:=ifelse(arb_clindiag==777,6,arb_clindiag)] # We set the NA value ("777") to a level "6" because we later use it for define CZS variable

#1.2. Set to NA : 666,777,888,999 values----
data[data==""] <-NA
data[data==666] <-  NA
data[data==777] <-  NA
data[data==888] <-  NA
data[data==999] <-  NA
data[data==9999] <-  NA

data[,Comb_id:=paste(studyname,mid_original,childid_original,sep="-")] #combinated ID from study mother and child id

# 1.3. Set to NA continuos and time variables out of the boundaries
Boundaries<-infoexp[Inclusion==1&Vtype%in%c("C","T"),c("Type","Variable","Vtype","Units","Min","Max")]
Boundaries[,Max:=ifelse(is.na(Max),"Inf",Max)]
boundvar=Boundaries$Variable
for(var in boundvar){
  data$var_bound<-data[,..var]
  Min=as.numeric(Boundaries[Variable==var]$Min)
  Max=as.numeric(ifelse(is.na(Boundaries[Variable==var]$Max),1/0,Boundaries[Variable==var]$Max))
  data[var_bound<=Min|var_bound>=Max,(var):=NA]
}  

#1.4. Check categorical variables----

data[,inf_sex:=ifelse(inf_sex==3,NA,inf_sex)] #3= NA
data[,endga:=ifelse(!is.na(endga),endga,birth_ga)] #end ga has same values as birthga.
data[fet_us_micro_tri2==1|fet_us_micro_tri3==1,microcephaly_bin:=1]
#data[loss_etiology==4,loss_etiology:=NA]

head(data)
#1.5. Check exposures----

#Studies including only women with zikv infection: Brazil_BahiaPaudaLima_Costa Brazil_SP_RibeiraoPreto_Duarte 
#Brazil_RiodeJaneiro_CunhaPrata Colombia_Mulkey FrenchGuiana_Pomar Spain_Soriano TrinidadTobago_Sohan
#zikv_preg can be set to 1 for all records from these studies
#zik_inc<-c("Brazil_RiodeJaneiro_CunhaPrata","Brazil_SP_RibeiraoPreto_Duarte","Brazil_BahiaPaudaLima_Costa","Spain_Soriano","FrenchGuiana_Pomar","TrinidadTobago_Sohan","Colombia_Mulkey")
zik_inc<-c("Brazil_RiodeJaneiro_CunhaPrata","Brazil_SP_RibeiraoPreto_Duarte","TrinidadTobago_Sohan")

data[, .(count = .N,
         ziv_pos=sum(ifelse(zikv_preg==1,1,0),na.rm=TRUE),
         ziv_neg=sum(ifelse(zikv_preg==0,1,0),na.rm=TRUE),
         ziv_na=sum(ifelse(is.na(zikv_preg),1,0),na.rm=TRUE),
         inclusion=ifelse(studyname%in%zik_inc,1,0))
     , by = studyname] 

#                         studyname count ziv_pos ziv_neg ziv_na inclusion
#3:    Brazil_BahiaPaudaLima_Costa    46      13      33      0         1
#1: Brazil_RiodeJaneiro_CunhaPrata    54       1       4     49         1-2 # apply inclusion criteria
#7:       Brazil_RiodeJaneiro_Joao   219     219       0      0         0
#2: Brazil_SP_RibeiraoPreto_Duarte   513     513       0      0         1-2
#9:                Colombia_Mulkey    85      77       0      8         1
#5:             FrenchGuiana_Pomar   291     291       0      0         1
#8:                  Spain_Bardaji   198       4     194      0         0   3 corrected multiple entry
#4:                  Spain_Soriano   278      27     206     45         1
#6:           TrinidadTobago_Sohan   104      83       3     18         1-2  4 removed inclusion criteria
#10:                    USA_Mulkey    95      61      17     17         0




#Studies that includes only woman with zika infection: Check with Diana and Mabell, look previous table
data[studyname%in%zik_inc,zik_preg:=1]
unique(data$studyname)
data[,studycode:=fcase (studyname == "Brazil_BahiaPaudaLima_Costa", "014-BRA",
                        studyname == "Brazil_RiodeJaneiro_CunhaPrata", "001-BRA",
                        studyname == "Brazil_RiodeJaneiro_Joao" , "002-BRA",
                        studyname == "Brazil_SP_RibeiraoPreto_Duarte", "010-BRA",
                        studyname == "Colombia_Mulkey", "007-COL",
                        studyname == "FrenchGuiana_Pomar", "003-GUF",
                        studyname == "Spain_Bardaji" , "005-ESP",
                        studyname == "Spain_Soriano", "004-ESP",
                        studyname == "TrinidadTobago_Sohan", "012-TTO",
                        studyname == "USA_Mulkey", "008-USA")]

#Zika test
data[,zikv_ga_min:=apply(data[,c("zikv_elisa_ga_1","zikv_pcr_ga_1")], 1, min, na.rm = TRUE)]
data[,zikv_ga_min:=ifelse(is.infinite(zikv_ga_min),NA,zikv_ga_min)]
data[,zikv_gan:=ifelse(is.na(zikv_ga),zikv_ga_min,zikv_ga_min)]
data[,zikv_ga:=ifelse(is.na(zikv_ga),zikv_ga_min,zikv_ga_min)]
data[,zikv_tri:=ifelse(is.na(zikv_ga),zikv_tri,ifelse(zikv_ga<13,0,ifelse(zikv_ga<=27,1,2)))]
data[,zikv_test:=fcase(zikv_pcr_res_1==1|(zikv_elisa_res_1==1&zikv_prnt_studydef_1==1),"robust",
                       zikv_elisa_res_1==1|zikv_prnt_studydef_1==1, "moderate_limit",
                       zikv_pcr_res_1==0&zikv_elisa_res_1==0,"negative",
                       zikv_pcr_res_1==0|zikv_elisa_res_1==0,"any_negative")] #according to test paper

#checktable<-as.data.table(table(elisa=data$zikv_elisa_res_1,elisa_everpos=data$zikv_elisa_everpos,pcr=data$zikv_pcr_res_1,zik=data$zikv_preg,arb=data$arb_clindiag,useNA = "always"))
#checktable[N>0] #ask for observations where zik value 0 and other there is 1 aroun 273..which has more priority on information?

#1.4.Check Pregnant woman variables---
data[tobacco==3,tobacco:=NA]

# 2. Microcephaly correction----
#2.0. Correction at baseline value----
data[,microcephaly1:=microcephaly] # given by hospital 
data[!is.na(inf_sex), hcircm2zscore:=as.numeric(igb_hcircm2zscore(gagebrth = birth_ga*7, hcircm=inf_head_circ_birth,sex=ifelse(inf_sex== 0, "Male","Female")))]  
data[, microcephaly2:= ifelse(hcircm2zscore<=-3,2,ifelse(hcircm2zscore<=-2,1,ifelse(hcircm2zscore<=2,0,ifelse(!is.na(hcircm2zscore),3,NA))))] # given by formula
data[, microcephaly:=ifelse(is.na(microcephaly2),microcephaly,microcephaly2)] # we give priority to microcephaly based on head circunference function
#data[is.na(microcephaly),microcephaly:=microcephaly2]
checkmic<-as.data.table(table(micbin=data$microcephaly_bin,mic1=data$microcephaly,mic2=data$microcephaly2,micn=data$microcephalyn,useNA = "always"))
checkmic[N>0]  #check this inconsistency.. i prioririze microcephaly variable
data[,microcephaly_bin:=ifelse(microcephaly%in%c(0,3),0,ifelse(microcephaly%in%c(1,2),1,microcephaly_bin))]

#2.1. Postnatal microcephaly---
for(t in 1:3){
  microvar<-paste0("microcephaly_bin_fu",t)
  var_age<-paste0("inf_head_circ_age_fu",t)
  var_circ<-paste0("inf_head_circ_age_fu",t)
  data$age_fu=data[,..var_age]*30
  data$hcirc_fu=data[,..var_circ]
  data[!is.na(inf_sex), hcircm2zscore_fu:=as.numeric(who_hcircm2zscore(agedays =age_fu, hcircm=hcirc_fu,sex=ifelse(inf_sex== 0, "Male","Female")))]  
  data[,micro_fu:= ifelse(hcircm2zscore_fu<=-3,2,ifelse(hcircm2zscore_fu<=-2,1,ifelse(hcircm2zscore_fu<=2,0,ifelse(!is.na(hcircm2zscore_fu),3,NA))))]
  data[,micro_bin_fu:=ifelse(micro_fu%in%c(0,3),0,ifelse(micro_fu%in%c(1,2),1,NA))]
  data[, (microvar):= micro_bin_fu]
}
data[,postsum:=rowSums(data[,c("microcephaly_bin_fu1","microcephaly_bin_fu2","microcephaly_bin_fu3")],na.rm=T)]
data[,microcephaly_posn:=ifelse(is.na(microcephaly_bin_fu1)&is.na(microcephaly_bin_fu2)&is.na(microcephaly_bin_fu3),NA,ifelse(postsum>0&microcephaly_bin==0,1,0))]


# 3. CZN correction----
#Function returns 1 if there is any anomaly detected and 0 if there were no anomaly detected across the recorded variables , NA no recorded variables.

checkcon<-function(data,col1){ #
  data$anyT <- rowSums(data[, .SD, .SDcols = col1], na.rm=T)
  data$allNA <- rowSums(!is.na(data[, .SD, .SDcols = col1]))
  data$Final<- ifelse(data$allNA==0,NA,ifelse(data$anyT>0,1,0))
  return(data$Final)
}

#Create variable related to neuroimaging abnormalities
data[,modificated1:=ifelse(fet_us_abn_spec_tri1==0,1,NA)] 
data[,modificated2:=ifelse(fet_us_abn_spec_tri2==0,1,NA)] 
data[,modificated3:=ifelse(fet_us_abn_spec_tri3==0,1,NA)] 
col1<-c("modificated1","modificated2","modificated3","hydrocephaly","calcifications","ventriculomegaly","fet_us_cns_tri2","fet_us_cns_tri3")
data[,neuroabnormality:=checkcon(data=data,col1=col1)]


#Create variable related to congenital contractures
data[,modificated1:=ifelse(fet_us_abn_spec_tri1==1,1,NA)] 
data[,modificated2:=ifelse(fet_us_abn_spec_tri2==1,1,NA)] 
data[,modificated3:=ifelse(fet_us_abn_spec_tri3==1,1,NA)] 
col1<-c("modificated1","modificated2","modificated3","fet_us_msk_tri2","fet_us_msk_tri3")
data[,contractures:=checkcon(data=data,col1=col1)]


#Create variable related to cardio abnormalities
data[,modificated1:=ifelse(fet_us_abn_spec_tri1==2,1,NA)] 
data[,modificated2:=ifelse(fet_us_abn_spec_tri2==2,1,NA)] 
data[,modificated3:=ifelse(fet_us_abn_spec_tri3==2,1,NA)]
col1<-c("modificated1","modificated2","modificated3","fet_us_cardio_tri2","fet_us_cardio_tri3")
data[,cardioabnormality:=checkcon(data=data,col1=col1)]


#Create variable related to gastrointestinal abnormalities
data[,modificated1:=ifelse(fet_us_abn_spec_tri1==3,1,NA)] 
data[,modificated2:=ifelse(fet_us_abn_spec_tri2==3,1,NA)] 
data[,modificated3:=ifelse(fet_us_abn_spec_tri3==3,1,NA)] 
col1<-c("modificated1","modificated2","modificated3","fet_us_gastro_tri2","fet_us_gastro_tri3")
data[,gastroabnormality:=checkcon(data=data,col1=col1)]


#Create variable related to orofacialintestinal abnormalities
data[,modificated1:=ifelse(fet_us_abn_spec_tri1==4,1,NA)] 
data[,modificated2:=ifelse(fet_us_abn_spec_tri2==4,1,NA)] 
data[,modificated3:=ifelse(fet_us_abn_spec_tri3==4,1,NA)] 
col1<-c("modificated1","modificated2","modificated3","fet_us_orofac_tri2","fet_us_orofac_tri3")
data[,oroabnormality:=checkcon(data=data,col1=col1)]


#Create variable related to ocular abnormalities
data[,modificated1:=ifelse(fet_us_abn_spec_tri1==5,1,NA)] 
data[,modificated2:=ifelse(fet_us_abn_spec_tri2==5,1,NA)] 
data[,modificated3:=ifelse(fet_us_abn_spec_tri3==5,1,NA)] 
col1<-c("modificated1","modificated2","modificated3","fet_us_eyeear_tri2","fet_us_eyeear_tri3")
data[,ocularabnormality:=checkcon(data=data,col1=col1)]


#Create variable related to genitourinaly abnormalities
data[,modificated1:=ifelse(fet_us_abn_spec_tri1==6,1,NA)] 
data[,modificated2:=ifelse(fet_us_abn_spec_tri2==6,1,NA)] 
data[,modificated3:=ifelse(fet_us_abn_spec_tri3==6,1,NA)] 
col1<-c("modificated1","modificated2","modificated3","fet_us_genur_tri2","fet_us_genur_tri3")
data[,genurabnormality:=checkcon(data=data,col1=col1)]


#Create variable any congenital abnormality excluding microcephaly
col1<-c("neuroabnormality","contractures","cardioabnormality","gastroabnormality","oroabnormality","ocularabnormality","genurabnormality","othabnorm","fet_us_bin_tri1","fet_us_bin_tri2","fet_us_bin_tri3")
data[,anyabnormality_czs:=checkcon(data=data,col1=col1)]
col1<-c(col1,"anyabnormality_czs")

# anyabnormality, nonneurologic: I did not calculated as is the combination of anyabnormality_czs and microcephaly so you can include the anyabnormality in the imputation  
# no sense to add imputation separaterly  neuronormality,ocularnormality, contractures and anyabnormalyty as all are combined  in anyabnormality_czs


#Create a czs variable according to WHO definition
#WHO definition for CZS: Presence of confirmed maternal or fetal ZIKV infection AND presence of severe microcephaly at birth AND presence of other malformations (including limb contractures, high muscle tone, eye abnormalities, and hearing loss, nose etc.)
data[,czs2:=ifelse((data$zikv_preg==1 | data$fet_zikv==1) & (data$microcephaly==2 |data$anyabnormality_czs==1),1,
                   ifelse(data$zikv_preg==0&data$fet_zikv==0 & data$microcephaly!=2&data$anyabnormality_czs==0,0,NA))] 
data[,czsn:=ifelse(is.na(czs),czs2,czs)]
data[,czs2:=NULL]

checktable<-data.table(table(czs=data$czs,czsn=data$czsn,useNA = "always"))
checktable[N>0] #Combinations where the czv calculated differs from the given one. czsn gave information of around 100 observations 


#4. Construct Bdeath (baby death variable) -----

#0. Put miscarriage, loss, loss_etiology, birth in one variable bdeath and bdeath_ga---
data[,birth:=ifelse(!is.na(birth_ga),1,NA)] #birth indicator from bithga
data[,birth2:=ifelse(!is.na(inf_term),1,NA)] #birth indicator from inf_term
data[,bdeath:=ifelse(birth==1|birth2==1,0,NA)]
data[,bdeath_ga:=ifelse(!is.na(endga),endga,ifelse(!is.na(loss_ga),loss_ga,ifelse(!is.na(miscarriage_ga),miscarriage_ga,NA)))]


#Use etiology to set death (et 0=live 1=miss 2=loss 3=imp death), 
data[loss_etiology==0&(is.na(bdeath)|bdeath==0),c("miscarriage","loss","loss_etiology","birth","bdeath"):= list(0,0,0,1,0)] #birth
data[loss_etiology==1&(is.na(bdeath)|bdeath==1),c("miscarriage","loss","loss_etiology","birth","bdeath"):= list(1,0,1,0,1)] #miscarriage
data[loss_etiology==2&(is.na(bdeath)|bdeath==1),c("miscarriage","loss","loss_etiology","birth","bdeath"):= list(0,1,2,0,1)] #loss
data[loss_etiology==3&(is.na(bdeath)|bdeath==1),c("miscarriage","loss","loss_etiology","birth","bdeath"):= list(0,1,2,0,1)] #class as loss
data[bdeath==1&bdeath_ga<20,c("miscarriage","loss","loss_etiology","birth","bdeath"):= list(1,0,1,0,1)]
data[bdeath==1&bdeath_ga>=20,c("miscarriage","loss","loss_etiology","birth","bdeath"):= list(0,1,2,0,1)]
data[loss==1&is.na(loss_etiology),c("miscarriage","loss","loss_etiology","birth","bdeath"):= list(0,1,2,0,1)]
data[birth==1&loss_etiology==1&inf_vital_status==0,c("miscarriage","loss","loss_etiology","birth","bdeath"):= list(0,0,0,1,0)]
data[birth==1&is.na(loss_etiology),c("miscarriage","loss","loss_etiology","birth","bdeath"):= list(0,0,0,1,0)]
data[miscarriage==0&loss==0&is.na(birth)&inf_vital_status==0,c("miscarriage","loss","loss_etiology","birth","bdeath"):= list(0,0,0,1,0)]
data[!is.na(data$inf_term),c("miscarriage","loss","loss_etiology","birth","bdeath"):= list(0,0,0,1,0)]
data[inducedabort==1&birth==1,c("miscarriage","loss","loss_etiology","birth","bdeath"):= list(0,1,2,0,1)]
data[birth==0&inf_vital_status==0,c("miscarriage","loss","loss_etiology","birth","bdeath"):= list(0,0,0,1,0)]


#Coherent with inf_alive_birth 0=alive, induce abort 0=No
checktable<-data.table(table(miss=data$miscarriage,loss=data$loss,et=data$loss_etiology,birth=data$birth,bdeath=data$bdeath,vst=data$inf_vital_status,iabo=data$inducedabort,inft=!is.na(data$inf_term), useNA="always")) #V1 is miscarriage, V2 is loss and N the number of observations
checktable[N!=0,]
data[,birth:=NULL]
data[,birth2:=NULL]



#5. Creation of additional covariates----

#5.1. Maternal prescription drug use----
col1<-c("med_bin","med_anticonvuls_bin","med_preg_bin","med_fertil_bin")
data[,drugs_prescr:=checkcon(data=data,col1=col1)]

#5.2. Maternal vaccination ----
col1<-c("vac_rub_enroll","vac_vari_enroll","vac_yf_enroll")
data[,vaccination:=checkcon(data=data,col1=col1)]

#5.3.Concurrent or prior flavi- or alpha virus infection ----
data[,modificated1:=ifelse(arb_clindiag_plus==0,0,ifelse(!is.na(arb_clindiag_plus),1,0))]
data[,modificated2:=ifelse(arb_clindiag!=0&arb_clindiag!=1,0,ifelse(!is.na(arb_clindiag),1,0))] #arb_clindiag==777 to 6 on top
col1<-c("modificated1","modificated2","denv_ever","chikv_ever")
data[,flavi_alpha_virus:=checkcon(data=data,col1=col1)]

#5.4. Intrauterine exposure to storch pathogens----
data[,modificated1:=ifelse(storch==0,0,ifelse(!is.na(storch),1,0))]
col1<-c("modificated1","storch_bin","toxo","toxo_treat","syphilis","syphilis_treat","varicella","parvo","rubella","cmv","herpes","listeria","chlamydia","gonorrhea","genitalwarts")
data[,storch_patho:=checkcon(data=data,col1=col1)]


#5.5.Prior arb virus infection ----
data[,modificated1:=ifelse(zikv_pcr_everpos==1,1,ifelse(zikv_pcr_everpos==0,0,NA))] # 2 indeterminated
col1<-c("modificated1","zikv_elisa_everpos","denv_ever","chikv_ever")
data[,arb_ever:=checkcon(data=data,col1=col1)]

#5.6.Current arb virus infection ----

data[,modificated1:=ifelse(zikv_pcr_res_1==1,1,ifelse(zikv_pcr_res_1==0,0,NA))] # 2 indeterminated
data[,modificated2:=ifelse(zikv_elisa_res_1==1,1,ifelse(zikv_elisa_res_1==0,0,NA))]
data[,modificated3:=ifelse(arb_clindiag==0,0,ifelse(!is.na(arb_clindiag),1,NA))] #arb_clindiag==777 to 6 modified at the beginning= other arbovirus
col1<-c("modificated1","modificated2","modificated3","zikv_preg","denv_preg","chikv_preg")
data[,arb_preg:=checkcon(data=data,col1=col1)]

#arb current pregnancy without consider zika
data[,modificated1:=ifelse(arb_clindiag==0|arb_clindiag==1,0,ifelse(!is.na(arb_clindiag),1,NA))] #arb_clindiag==777 to 6 modified at the beginning= other arbovirus
col1<-c("modificated1","denv_preg","chikv_preg")
data[,arb_preg_nz:=checkcon(data=data,col1=col1)]


#6. Interactive plot ------
var_inc<-infoexp[order(Order),][!is.na(Inclusion)]$Variable
pldata<-data[, .SD, .SDcols = var_inc]

save(pldata, file = "3_Output_data/pldata.RData")

totalval<-as.data.table(table(data$studycode))
colnames(totalval)<-c("studycode","N")
dmatrix<-pldata[, lapply(.SD, function(x) sum(is.na(x))/.N), studycode] #matrix of % of missingness
dmatrix2<-as.data.table(melt(dmatrix,id.vars="studycode"))
dmatrix2<-merge(dmatrix2,totalval,by="studycode")

dmatrix2[,name:=paste0(studycode,"\nN=",N)]
dmatrix2[,miss:=round(value*100,1)]
dmatrix2[,text:=paste0("study: ", studycode, "\n", "variable: ", variable, "\n", "miss%: ",miss)]
dmatrix2[,tmiss:=ifelse(miss==100,"Systematical","Sporadical")]

p<-ggplot(dmatrix2, aes(x=name, y=variable,fill=miss,text=text)) +
  geom_tile() +
  scale_fill_gradientn(colours=c("green","yellow","red")) +
  theme(axis.text.x = element_text(angle = 45,vjust=0.5,size=8),
        axis.text.y = element_text(size=8))+xlab("Study name")+ylab("Variable")+
  labs(fill='%missing') 

ggplotly(p, tooltip="text")

p2<- ggplot(dmatrix2, aes(x=name, y=variable,fill=tmiss,text=text)) +
  geom_tile() +
  scale_fill_manual(values=c("green","red")) +
  theme(axis.text.x = element_text(angle = 45,vjust=0.5,size=8),
        axis.text.y = element_text(size=8))+xlab("Study name")+ylab("Variable")+
  labs(fill='Type of missing') 

ggplotly(p2, tooltip="text")

#7. Flux plot -----
fx<-flux(data)
fluxplot(data)
outlist<-row.names(fx)[fx$outflux>=0.5]
infoexp[,outflux2:=ifelse(Variable%in%outlist,1,0)]

##8. Select variable according number of studies not completely incomplete----
Nstudies=length(unique(data$studyname))
cols<-colnames(dmatrix)[-1]
dmatrixp<-dmatrix[ , (cols) := lapply(.SD, function(x){ifelse(x<0.8,1,0)}), .SDcols = cols] #Calculate which variables have less htan 80% of missingness
sumstudy<-as.numeric(colSums(dmatrixp[,-1]))
stable<-as.data.table(cbind(Nstudycomp=as.numeric(ifelse(sumstudy>0.3*Nstudies,1,0)),Variable=cols)) #number of studies with complete information above the 30% of total number of studies
infoexp<-merge(infoexp,stable,by="Variable",all = TRUE)

##9. Calculate inclusion criteria -----
infoexp[,Nstudycomp:=as.numeric(Nstudycomp)]
infoexp[,Imputation2 := rowSums(.SD, na.rm = TRUE), .SDcols = c("Inclusion", "Additionals", "outflux2","Nstudycomp")]

write.table(infoexp,file="5_Internal_support/Infoselection.csv",sep=";")
selecvar<-infoexp[order(Order)][Imputation2>=3,]$Variable #Vector with ordered selected variables

fdata<-data[,..selecvar]
save(fdata, file = "3_Output_data/finaldata.RData")


#10. Plot final----
nindexf<-infoexp[Imputation2>2,]
fdata<-data[, .SD, .SDcols = nindexf$Variable]
dmatrix<-fdata[, lapply(.SD, function(x) sum(is.na(x))/.N), studyname]
dmatrix2<-melt(dmatrix,id.vars="studyname")
dmatrix2<-merge(dmatrix2,totalval,by="studyname")
dmatrix2[,name:=fcase(
  studyname=="Brazil_RiodeJaneiro_CunhaPrata","Brazil\nCunhaPrata",
  studyname=="Brazil_SP_RibeiraoPreto_Duarte","Brazil\n Duarte",
  studyname=="Brazil_BahiaPaudaLima_Costa","Brazil\nCosta",
  studyname=="Spain_Soriano" ,"Spain\nSoriano",
  studyname=="FrenchGuiana_Pomar" ,"Fr.Guiana\nPomar",
  studyname=="TrinidadTobago_Sohan" ,"Tri.Tobago\nSohan",
  studyname=="Brazil_RiodeJaneiro_Joao" ,"Brazil\nJoao",
  studyname=="Spain_Bardaji" ,"Spain\nBardaji",
  studyname=="Colombia_Mulkey"  ,"Colombia\nMulkey",
  studyname=="USA_Mulkey"  ,"USA\nMulkey"
)]
dmatrix2[,name:=paste0(name,"\nN=",N)]
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

p2<- ggplot(dmatrix2, aes(x=name, y=variable,fill=tmiss,text=text)) +
  geom_tile() +
  scale_fill_manual(values=c("green","red")) +
  theme(axis.text.x = element_text(angle = 45,vjust=0.5,size=8),
        axis.text.y = element_text(size=8))+xlab("Study name")+ylab("Variable")+
  labs(fill='Type of missing') 

ggplotly(p2, tooltip="text")


####END of CODE##########



# MNAR for microcephaly
table(micro=data$microcephaly_bin,bdeath=data$bdeath,useNA = "always")
sdata<-data[studyname%in%c("Brazil_BahiaPaudaLima_Costa","Brazil_RiodeJaneiro_CunhaPrata","Spain_Bardaji","Spain_Soriano","TrinidadTobago_Sohan"," USA_Mulkey")]
colnames(sdata)
colned=c("studyname","zik_preg","drugs_prescr","inducedabort","age", "parity", "gravidity","inducedabort", "alcohol", "drugs","tobacco", "vaccination","eclampsia","storch_patho","arb_symp","zikv_preg","microcephaly_bin")
sdata<-sdata[,..colned]
data[,ry:=ifelse(is.na(microcephaly_bin),0,1)]
data<-sdata[!is.na(parity)]
data[,firstp:=ifelse(parity==0,1,0)]
library(logistf)

fit5<-glm(ry~parity+age+storch_patho+arb_symp+arb_ever+inducedabort+zikv_preg+arb_symp+storch_patho,data=sdata,family=binomial("logit"))
fit<-logistf(ry~firstp+age+arb_preg_nz+arb_symp+fever+rash+zikv_preg+storch_patho+vaccination,data=data, pl=FALSE,control=logistf.control(maxit=10000, maxstep=100))
fit<-logistf(microcephaly_bin~firstp+age+arb_symp+arb_preg_nz+fever+rash+zikv_preg+storch_patho+vaccination,data=data, pl=FALSE,control=logistf.control(maxit=10000, maxstep=100))

summary(fit)




