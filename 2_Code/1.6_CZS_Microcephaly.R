
#Check inconsistences between Microcephaly and CZS 25-03-2022

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


#Zika test
data[,zikv_ga_min:=apply(data[,c("zikv_elisa_ga_1","zikv_pcr_ga_1")], 1, min, na.rm = TRUE)]
data[,zikv_ga_min:=ifelse(is.infinite(zikv_ga_min),NA,zikv_ga_min)]
data[,zikv_gan:=ifelse(is.na(zikv_ga),zikv_ga_min,zikv_ga_min)]
data[,zikv_ga:=ifelse(is.na(zikv_ga),zikv_ga_min,zikv_ga_min)]
data[,zikv_tri:=ifelse(is.na(zikv_ga),zikv_tri,ifelse(zikv_ga<13,0,ifelse(zikv_ga<=27,1,2)))]



#Studies that includes only woman with zika infection: Check with Diana and Mabell, look previous table
data[studyname%in%zik_inc,zik_preg:=1]


#Zika test
data[,zikv_ga_min:=apply(data[,c("zikv_elisa_ga_1","zikv_pcr_ga_1")], 1, min, na.rm = TRUE)]
data[,zikv_ga_min:=ifelse(is.infinite(zikv_ga_min),NA,zikv_ga_min)]
data[,zikv_gan:=ifelse(is.na(zikv_ga),zikv_ga_min,zikv_ga_min)]
data[,zikv_ga:=ifelse(is.na(zikv_ga),zikv_ga_min,zikv_ga_min)]
data[,zikv_tri:=ifelse(is.na(zikv_ga),zikv_tri,ifelse(zikv_ga<13,0,ifelse(zikv_ga<=27,1,2)))]


#checktable<-as.data.table(table(elisa=data$zikv_elisa_res_1,elisa_everpos=data$zikv_elisa_everpos,pcr=data$zikv_pcr_res_1,zik=data$zikv_preg,arb=data$arb_clindiag,useNA = "always"))
#checktable[N>0] #ask for observations where zik value 0 and other there is 1 aroun 273..which has more priority on information?


#1.4.Check Pregnant woman variables---
data[tobacco==3,tobacco:=NA]

# 2. Microcephaly correction----
#2.0. Correction at baseline value----
data[!is.na(inf_sex), hcircm2zscore:=as.numeric(igb_hcircm2zscore(gagebrth = birth_ga*7, hcircm=inf_head_circ_birth,sex=ifelse(inf_sex== 0, "Male","Female")))]  
data[, microcephaly2:= ifelse(hcircm2zscore<=-3,2,ifelse(hcircm2zscore<=-2,1,ifelse(hcircm2zscore<=2,0,ifelse(!is.na(hcircm2zscore),3,NA))))]
data[, microcephalyn:=ifelse(is.na(microcephaly),microcephaly2,microcephaly)]



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
data[,czs2:=ifelse((zikv_preg==1 | fet_zikv==1) & microcephaly==2 & anyabnormality_czs==1,1,
                   ifelse(zikv_preg==0 & fet_zikv==0 & microcephaly!=2& anyabnormality_czs==0,0,NA))] 
data[,czs3:=ifelse((zikv_preg==1 | fet_zikv==1) & microcephalyn==2 & anyabnormality_czs==1,1,
                   ifelse(zikv_preg==0 & fet_zikv==0 & microcephalyn!=2& anyabnormality_czs==0,0,NA))] 
data[,czsn:=ifelse(is.na(czs),czs3,czs)] #combined


###
###Original relationship
##Microcephaly(given), CZS (given)
table(data$microcephaly,data$czs,useNA = "always")

#       0   1 <NA>
#0    886  23  364
#1     20   6   18
#2      1  11    5
#3      0   0    7
#<NA>   2   0  533

# Microcephaly(given), Microcephaly(circunference score)
##Rows vs Columns
table(data$microcephaly,data$microcephaly2,useNA = "always")

####    0   1   2   3 <NA>
#0    586   0   0  64  623
#1      3   7   0   0   34
#2      0   0   6   0   11
#3      0   0   0   0    7
#<NA>  77   4   7  12  435

#So we can observe here that the macrocephaly cases comes from the evaluation of circunference score, 
#basically what we observe is that 64 cases reported as normal have macrocephaly (by circunference) and
#3 cases microchepaly reported cases are normal (by circunference).
# Also you can see that not all the babies can be caracterized with macrocephaly because measurements were
# taken if i am not wrong outside of the timeframe valid to evaluate the circunference. i.e. formula does not
#return any value.


# CSZ(given), CZS(WHO score with original microcephaly value)
table(data$czs,data$czs2,useNA="always")
#       0   1 <NA>
#0      0   0  909
#1      0   9   31
#<NA>  92   2  833

#So basically we observe that the WHO rule only allowed us to classify 2 NA cases as positive CZS cases. 
# At the end we dont do a lot of changes with the WHO rule, but who knows if this will be usefull in the 
# future with the inclusion of new studies

# CSZ(given), CZS(WHO score with combined microcephaly value)

table(data$czs,data$czs3,useNA="always")
#       0   1 <NA>
#0      0   0  909
#1      0   9   31
#<NA>  92   5  830

# So now we classify 5 NA cases as positive CZS cases.

# What we want to see is the relationship between microcephaly and czs as are included before the imputation,
#i.e. czsn, microcephaly$
table(data$microcephalyn,data$czsn,useNA = "always")

#       0   1 <NA>
#0    970  23  357
#1     25   6   17
#2      1  16    7
#3      3   0   16
#<NA>   2   0  433






