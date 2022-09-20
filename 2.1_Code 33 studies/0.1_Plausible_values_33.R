
rm(list=ls()) # clean environment

# Data manipulation packages
library(here) 
library(data.table)
library(dplyr)
library(haven)
require(ggplot2)

#data <- as.data.table(read_dta(here('1_Input_data','zikv_033_datasets_labelandmetada.dta')))
#data <- as.data.table(readxl::read_xlsx(here('1_Input_data','zikv_033_datasets.xlsx'), stringsAsFactors=FALSE, fileEncoding="latin1"))
#data<-as.data.table(read.csv('1_Input_data/pilot10_08SEP21_withmetadata.csv', stringsAsFactors=FALSE, fileEncoding="latin1"))
data <- as.data.table(readxl::read_xlsx(here('1_Input_data','zikv_033_datasets.xlsx'),sheet="Sheet1"))
infoexp <- as.data.table(readxl::read_xlsx(here('1_Input_data','MasterCodebook_Final_June2022 ALL (Repaired).xlsx'),sheet="237 key")) #CSV file with the
var_inc<-as.matrix(infoexp[Essential=="Yes",'WHO variable name'])[,1] #Variables to work with
data<-as.data.table(data[,..var_inc]) #Filter dataset

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

table(data$studycode,useNA = "always")


#1. Initial clean up the dataset -----
#1.0. Check for duplicates ----

data[,file:=ifelse(is.na(file),studycode,file)]
data[,comb:=paste0(file,mid,fetid,childid,sep="_")] #create a variable that combines study and patient id
n_occur <- data.table(table(data$comb))
n_occur[n_occur$N>1] #check visually which observations are duplicated

var1<-n_occur[n_occur$N > 1]$V1
Duplicates<-data[comb%in%var1&multiplegest!=1] # table to report to harmonization team.
write.csv(Duplicates,file='5_Internal_support/Possible_duplicates.csv')


# 1.1. Check format of categorical variables----
char <-   which(sapply(data, is.character)) #which are character
datac<-data[,..char]
#data[zikv_pcr_date_1%in%c("22027","888",""), zikv_pcr_date_1:=NA]
#data[zikv_elisa_date_1%in%c(""),zikv_elisa_date_1:=NA]
#data[symp_date%in%c(""),symp_date:=NA]
#data[,arb_clindiag:=ifelse(arb_clindiag==777,6,arb_clindiag)] # We set the NA value ("777") to a level "6" because we later use it to define CZS variable
#unique(data$age)
#mid_original                    mid                date_t1                    age 
#zikv_symp_date      zikv_assay_date_1        zikv_pcr_date_1      zikv_elisa_date_1 
#zikv_igm_titer_1         fetid_original                  fetid       childid_original 
#childid     ch_head_circ_birth     ch_head_circ_age_1                   file 
#unique(data$travel)

#1.1. Set to NA : 666,777,888,999 values----
data[data==""] <-NA
data[data==555] <-  NA
data[data==666] <-  NA
data[data==888] <-  NA
data[data==999] <-  NA
data[data==9999] <-  NA
data[data==999] <-  NA
#data[data==777] <-  NA
data[data=="NA"] <-NA 


# Plausibility bounds

Boundaries<-infoexp[Essential=="Yes"&Type_var%in%c("Continuous"),c("Tab","WHO variable name","Type_var","Units","Min","Max")]
Boundaries[,Max:=ifelse(is.na(Max),"Inf",Max)]
colnames(Boundaries)<-c("Tab","Var_name","Type_var","Units","Min","Max")
boundvar=Boundaries$Var_name
obs_to_check<-NULL
for(var in boundvar){
  data$var_bound<-data[,..var]
  data$Var_to_check<-var
  Min=as.numeric(Boundaries[Var_name==var]$Min)
  Max=as.numeric(ifelse(is.na(Boundaries[Var_name==var]$Max),1/0,Boundaries[Var_name==var]$Max))
  data_out<-data[var_bound<Min|var_bound>Max]
  obs_to_check<-rbind(obs_to_check,data_out)
}    
save(obs_to_check,file=here('3_Output_data','obs_to_check_33.RData'))




unique(data$weight)

val_check<-obs_to_check[,.(min=min(var_bound),max=max(var_bound), number=.N),by=Var_to_check]

ggplot(data=obs_to_check[Var_to_check=="height",],aes(x=var_bound))+
  geom_histogram(color="darkblue", fill="lightblue")+ theme_minimal()+
  labs(title="Height no plausible values",x="Height(cm)", y = "Count")

ggplot(data=obs_to_check[Var_to_check=="age",],aes(x=var_bound))+
  geom_histogram(color="darkblue", fill="lightblue")+ theme_minimal()+
  labs(title="Age no plausible values",x="Age(years)", y = "Count")

ggplot(data=obs_to_check[Var_to_check=="zikv_symp_sampcol_1",],aes(x=var_bound))+
  geom_histogram(color="darkblue", fill="lightblue")+ theme_minimal()+
  labs(title="zikv_symp_sampcol_1 No plausible values,from symptoms onset to Zika sample ",x="Days", y = "Count")

ggplot(data=obs_to_check[Var_to_check=="zikv_prnt_titer_1",],aes(x=var_bound))+
  geom_histogram(color="darkblue", fill="lightblue")+ theme_minimal()+
  labs(title="zikv_prnt_titer_1 No plausible values",x="zikv_prnt_titer", y = "Count")

ggplot(data=obs_to_check[Var_to_check=="zikv_ga",],aes(x=var_bound))+
  geom_histogram(color="darkblue", fill="lightblue")+ theme_minimal()+
  labs(title="zikv_ga No plausible values",x="weeks", y = "Count")


#4. Construct Bdeath (baby death variable) -----


#0. Put miscarriage, loss, loss_etiology, birth in one variable bdeath and bdeath_ga---
data[,birth:=ifelse(!is.na(birth_ga),1,NA)] #birth indicator from bithga
data[,birth2:=ifelse(!is.na(ch_term),1,NA)] #birth indicator from inf_term=ch_term?
data[,bdeath:=ifelse(birth==1|birth2==1,0,NA)]  # endga?
data[,bdeath_ga:=ifelse(!is.na(loss_ga),loss_ga,ifelse(!is.na(miscarriage_ga),miscarriage_ga,NA))]

data[,end_ga:=ifelse(!is.na(endga),endga,ifelse(!is.na(birth_ga),birth_ga,bdeath_ga))]
#Use etiology to set death (et 0=live 1=miss 2=loss 3=imp death), 
data[loss_etiology==0&(is.na(bdeath)|bdeath==0),c("miscarriage","loss","loss_etiology","birth","bdeath"):= list(0,0,0,1,0)] #birth
data[loss_etiology==1&(is.na(bdeath)|bdeath==1),c("miscarriage","loss","loss_etiology","birth","bdeath"):= list(1,0,1,0,1)] #miscarriage
data[loss_etiology==2&(is.na(bdeath)|bdeath==1),c("miscarriage","loss","loss_etiology","birth","bdeath"):= list(0,1,2,0,1)] #loss
data[loss_etiology==3&(is.na(bdeath)|bdeath==1),c("miscarriage","loss","loss_etiology","birth","bdeath"):= list(0,1,2,0,1)] #class as loss
data[bdeath==1&bdeath_ga<20,c("miscarriage","loss","loss_etiology","birth","bdeath"):= list(1,0,1,0,1)]
data[bdeath==1&bdeath_ga>=20,c("miscarriage","loss","loss_etiology","birth","bdeath"):= list(0,1,2,0,1)]
data[loss==1&is.na(loss_etiology),c("miscarriage","loss","loss_etiology","birth","bdeath"):= list(0,1,2,0,1)]
data[birth==1&loss_etiology==1&ch_vital_status==0,c("miscarriage","loss","loss_etiology","birth","bdeath"):= list(0,0,0,1,0)]
data[birth==1&is.na(loss_etiology),c("miscarriage","loss","loss_etiology","birth","bdeath"):= list(0,0,0,1,0)]
data[miscarriage==0&loss==0&is.na(birth)&ch_vital_status==0,c("miscarriage","loss","loss_etiology","birth","bdeath"):= list(0,0,0,1,0)]
data[!is.na(data$ch_term),c("miscarriage","loss","loss_etiology","birth","bdeath"):= list(0,0,0,1,0)]
data[inducedabort==1&birth==1,c("miscarriage","loss","loss_etiology","birth","bdeath"):= list(0,1,2,0,1)]
data[birth==0&ch_vital_status==0,c("miscarriage","loss","loss_etiology","birth","bdeath"):= list(0,0,0,1,0)]
data[birth==1,c("miscarriage","loss","inducedabort"):=list(0,0,0)]

#Coherent with inf_alive_birth 0=alive, induce abort 0=No
checktable<-data.table(table(miss=data$miscarriage,loss=data$loss,et=data$loss_etiology,birth=data$birth,bdeath=data$bdeath,vst=data$ch_vital_status,iabo=data$inducedabort,inft=!is.na(data$ch_term), useNA="always")) #V1 is miscarriage, V2 is loss and N the number of observations
checktable[N!=0,]
data[,birth:=NULL]
data[,birth2:=NULL]




# PCR test

data$chikv_pcr_res_1
data$denv_pcr_res_1
data$zikv_pcr_date_1
data$zikv_pcr_everpos
data$zikv_pcr_ga_1
data$zikv_pcr_res_1
data$zikv_pcr_tri_1
data$zikv_pcr_vl_1


datapcr<-data[, c("chikv_pcr_res_1","denv_pcr_res_1","zikv_pcr_date_1",
         "zikv_pcr_everpos","zikv_pcr_ga_1","zikv_pcr_res_1",
         "zikv_pcr_tri_1","zikv_pcr_vl_1","end_ga","studycode","comb")]

datanend<-datapcr[!is.na(zikv_pcr_ga_1)&is.na(end_ga),]

ggplot(data=datanend,aes(x=zikv_pcr_ga_1))+
  geom_histogram(color="darkblue", fill="lightblue")+ theme_minimal()+
  labs(title="zikv_pcr_ga without birth or death date",x="weeks", y = "Count")


table(is.na(datapcr$zikv_pcr_date_1),is.na(datapcr$zikv_pcr_ga_1),useNA = "always")
datapcr[,pcr_preg:=ifelse(end_ga>zikv_pcr_ga_1,"during","after")]
table(datapcr$pcr_preg,datapcr$studycode,useNA = "always")





datapcr[,ga:=is.na(zikv_pcr_ga_1)]
datapcr[,date:=is.na(zikv_pcr_date_1)]
datapcr[,tri:=is.na(zikv_pcr_tri_1)]

table(ga=datapcr$ga,date=datapcr$tri)
datapcr[ga==TRUE&date==FALSE,]



