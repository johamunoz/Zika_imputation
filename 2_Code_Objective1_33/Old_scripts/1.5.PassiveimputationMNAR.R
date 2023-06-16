gato<-data[!is.na(inf_head_circ_birth)]
data[]

c("Brazil_RiodeJaneiro_CunhaPrata","Brazil_BahiaPaudaLima_Costa","Spain_Bardaji","Spain_Soriano","TrinidadTobago_Sohan","USA_Mulkey")

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
data<-as.data.table(data[,..var_inc]) #Filter dataset

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

#1.1. Set to NA : 666,777,888,999 values----
data[data==""] <-NA
data[data==666] <-  NA
data[data==777] <-  NA
data[data==888] <-  NA
data[data==999] <-  NA
data[data==9999] <-  NA

data[,Comb_id:=paste(studyname,mid_original,childid_original,sep="-")] #combinated ID from study mother and child id
#1.2. Check outcome variables----
data[,inf_weight:=ifelse(inf_weight>6000|inf_weight<100,NA,inf_weight)]
data[,inf_sex:=ifelse(inf_sex==3,NA,inf_sex)] #3= NA
data[inf_length<18,inf_length:=NA] #only babies with length bigger than 18 cm
data[inf_head_circ_birth==0,inf_head_circ_birth:=NA] #babies with 0 circumference
data[,endga:=ifelse(!is.na(endga),endga,birth_ga)] #end ga has same values as birthga.
data[fet_us_micro_tri2==1|fet_us_micro_tri3==1,microcephaly_bin:=1]
data[loss_etiology==4,loss_etiology:=NA]


#1.3. Check exposures----

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
data[ miscarriage_ga<=0,miscarriage_ga:=NA]# codebook assigned to 666
data[ loss_ga<=0,loss_ga:=NA]# codebook assigned to 666
data[ birth_ga<=0,birth_ga:=NA]# codebook assigned to 666       
data[ zikv_pcr_ga_1<=0,zikv_pcr_ga_1:=NA] # codebook assigned to 666
data[ zikv_elisa_ga_1<=0,zikv_elisa_ga_1:=NA] # codebook assigned to 666
data[ zikv_ga<=0,zikv_ga:=NA]# codebook assigned to 666
data[ symp_ga<=0,symp_ga:=NA]# codebook assigned to 666
data[ arb_clindiag_ga<=0,arb_clindiag_ga:=NA]# codebook assigned to 666


#Zika test
#data[,zikv_ga_min:=apply(data[,c("zikv_elisa_ga_1","zikv_pcr_ga_1")], 1, min, na.rm = TRUE)]
#data[,zikv_ga_min:=ifelse(is.infinite(zikv_ga_min),NA,zikv_ga_min)]
#data[,zikv_gan:=ifelse(is.na(zikv_ga),zikv_ga_min,zikv_ga_min)]
#data[,zikv_gan:=ifelse(is.na(zikv_ga),zikv_ga_min,zikv_ga_min)]
#data[,zikv_ga_min:=NULL]
#data[,zikv_trin:=ifelse(is.na(zikv_gan),zikv_tri,ifelse(zikv_gan<13,0,ifelse(zikv_gan<=27,1,2)))]


#checktable<-as.data.table(table(elisa=data$zikv_elisa_res_1,elisa_everpos=data$zikv_elisa_everpos,pcr=data$zikv_pcr_res_1,zik=data$zikv_preg,arb=data$arb_clindiag,useNA = "always"))
#checktable[N>0] #ask for observations where zik value 0 and other there is 1 aroun 273..which has more priority on information?

#Zika symptoms
#data[,symp_gan:=ifelse(is.na(symp_ga),arb_clindiag_ga,symp_ga)]
#data[,symp_trin:=ifelse(is.na(symp_gan),symp_tri,ifelse(symp_gan<13,0,ifelse(symp_gan<=27,1,2)))]

#checktable<-as.data.table(table(arb=data$arb_symp,fev=data$fever,rash=data$rash,con=data$conjunctivitis,useNA="always"))
#checktable[N>0]   
#table(data$arb_clindiag_plus,data$arb_clindiag,useNA = "always")

#1.4.Check Pregnant woman variables---
data[age<14,age:=NA] #  ask "Spain_Soriano" for this age records, for the moment assigned to NA
data[tobacco==3,tobacco:=NA]


# 2. Microcephaly correction----
data[!is.na(inf_sex), hcircm2zscore:=as.numeric(igb_hcircm2zscore(gagebrth = birth_ga*7, hcircm=inf_head_circ_birth,sex=ifelse(inf_sex== 0, "Male","Female")))]  
data[, microcephaly2:= ifelse(hcircm2zscore<=-3,2,ifelse(hcircm2zscore<=-2,1,ifelse(hcircm2zscore<=2,0,ifelse(!is.na(hcircm2zscore),3,NA))))]
data[is.na(microcephaly),microcephaly:=microcephaly2]
checkmic<-as.data.table(table(micbin=data$microcephaly_bin,mic1=data$microcephaly,mic2=data$microcephaly2,useNA = "always"))
checkmic[N>0]  #check this inconsistency.. i prioririze microcephaly variable
data[,microcephaly_bin:=ifelse(microcephaly%in%c(0,3),0,ifelse(microcephaly%in%c(1,2),1,microcephaly_bin))]
