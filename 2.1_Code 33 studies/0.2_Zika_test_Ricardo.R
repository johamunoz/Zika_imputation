rm(list=ls()) # clean environment

# ZIKV Test according to Ricardo's paper

# Data manipulation packages
library(here) 
library(data.table)
library(dplyr)
library(haven)
require(ggplot2)

#data<-as.data.table(read.csv("/Users/jdamen/Documents/Julius/ZIKV analyses/2. Data/zikv_033_datasets.csv",header=T)) #Line for Anneke
data <- as.data.table(readxl::read_xlsx(here('1_Input_data','zikv_033_datasets.xlsx'),sheet="Sheet1")) #Line for Johanna
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



# Test modification
data[,bdeath_ga:=ifelse(!is.na(loss_ga),loss_ga,ifelse(!is.na(miscarriage_ga),miscarriage_ga,NA))]
data[,endga:=ifelse(!is.na(endga),endga,ifelse(!is.na(birth_ga),birth_ga,bdeath_ga))]
var_test<-c("childid","endga",paste0("zikv_pcr_ga_",1:10),paste0("zikv_pcr_res_",1:10),
            paste0("zikv_pcr_tri_",1:10),paste0("zikv_pcr_date_",1:10), 
            paste0("zikv_igm_res_",1:10),paste0("zikv_elisa_ga_",1:10),
            paste0("zikv_prnt_titer_",1:2),
            paste0("zikv_prnt_",1:10))

data_test<-as.data.table(data[,..var_test])
colA = paste("zikv_pcr_ga_", 1:10, sep = "")
colB = paste("zikv_pcr_res_", 1:10, sep = "")
colC = paste("zikv_pcr_tri_", 1:10, sep = "")
colD = paste("zikv_pcr_date_", 1:10, sep = "")
colE = paste("zikv_igm_res_", 1:10, sep = "")
colF = paste("zikv_prnt_titer_", 1:2, sep = "")
colG = paste("zikv_prnt_", 1:10, sep = "")
colH = paste("zikv_elisa_ga_", 1:10, sep = "")

data_test_m = melt(data_test, measure = list(colA, colB,colC,colD,colE,colF,colG,colH), 
                   value.name = c("pcr_ga", "pcr_res","pcr_tri","pcr_date","igm_res","prnt_titer","prnt","elisa_ga"),
                   variable.name="visit")
data_test_m[,pcr_res:=ifelse(pcr_res=="FALSE",0,ifelse(pcr_res=="TRUE",1,as.numeric(pcr_res)))]
data_test_m[,pcr_tri:=ifelse(pcr_tri=="FALSE",0,ifelse(pcr_tri=="TRUE",1,as.numeric(pcr_tri)))]
data_test_m[,igm_res:=ifelse(igm_res=="FALSE",0,ifelse(igm_res=="TRUE",1,as.numeric(igm_res)))]
data_test_m[,prnt:=as.numeric(prnt)]
data_test_m[,prnt_0:=ifelse(is.na(prnt),NA,ifelse(prnt%in%c(0),1,0))]
data_test_m[,prntz:=ifelse(!prnt%in%c(0,1),NA,prnt)]
data_test_m[,igmz:=ifelse(!igm_res%in%c(0,1),NA,igm_res)]
data_test_m[,pcr_res0:=ifelse(is.na(pcr_res),NA,ifelse(pcr_res%in%c(0),1,0))]
data_test_m[,pcr_res1:=ifelse(is.na(pcr_res),NA,ifelse(pcr_res%in%c(1),1,0))]
data_test_m[,pcr_res2:=ifelse(is.na(pcr_res),NA,ifelse(pcr_res%in%c(2),1,0))]
data_test_m[,pcr_resNA:=ifelse(is.na(pcr_res),1,0)]

data_test_m[,igm_res0:=ifelse(is.na(igm_res),NA,ifelse(igm_res%in%c(0),1,0))]
data_test_m[,igm_res2:=ifelse(is.na(igm_res),NA,ifelse(igm_res%in%c(2),1,0))]

data_test_m<-data_test_m[order(childid,visit)]
data_test_m[, prnt_titer_20:=ifelse(prnt_titer>=20,1,0)]
data_test_m[, prnt_titer_100:=ifelse(prnt_titer>=100,1,0)]
data_test_m[, prnt_titer_1000:=ifelse(prnt_titer>=1000,1,0)]
data_test_m[, endga:=as.numeric(endga)]
data_test_m[, preg0:=0]
data_test_m[, post1:=endga+4]
data_test_m[, post2:=endga+8]
data_test_m[, post3:=endga+13]
data_test_m[, post6:=endga+26]


testfunc<-function(datav,respv,timev, minv, maxv,name,diff,func){
  datav1<-copy(datav)
  datav1[,min:= get(minv)]
  datav1[,max:= get(maxv)]
  datav1[,time:=get(timev)]
  datav1[,resp:=get(respv)]
  dataf<-datav1[time<=max&time>=min&!is.na(resp),]
  dataf[,resp_lag:=c(NA, resp[-.N]),by=childid]
  dataf[,resp_diff:=resp-resp_lag]
  if(diff==TRUE){
    datares<-dataf[, any(resp_diff==1), by=childid]
  }else if(func=="all"){
    datares<-dataf[, all(resp==1), by=childid] 
  }else if(func=="any"&diff==FALSE){
    datares<-dataf[, any(resp==1), by=childid] 
  }else{
    datares<-dataf[, max(resp), by=childid] 
  }
  colnames(datares)<-c("childid",name)
  datares<-merge(datav,datares,by="childid",all.x = TRUE)
  return(datares)
}

# Robust
robs1<-testfunc(datav=data_test_m,respv="pcr_res1",timev="pcr_ga",minv="preg0", maxv="endga",name="pcr_pos_any_preg",func="any",diff=FALSE)
robs2<-testfunc(datav=robs1,respv="igmz",timev="elisa_ga",minv="preg0", maxv="endga",name="igm_zero_any_preg",func="any",diff=TRUE)
robs3<-testfunc(datav=robs2,respv="prntz",timev="elisa_ga",minv="preg0", maxv="endga",name="prnt_zero_any_preg",func="any",diff=TRUE)
robs40<-testfunc(datav=robs3,respv="igmz",timev="elisa_ga",minv="preg0", maxv="endga",name="igm_pos_any_preg",func="any",diff=FALSE)
robs41<-testfunc(datav=robs40,respv="prnt_titer_20",timev="elisa_ga",minv="preg0", maxv="endga",name="prnt20_any_preg",func="any",diff=FALSE)
robs42<-testfunc(datav=robs41,respv="prntz",timev="elisa_ga",minv="endga", maxv="post6",name="prnt_zero_any_post6",func="any",diff=TRUE)
robs43<-testfunc(datav=robs42,respv="prnt_titer_20",timev="elisa_ga",minv="endga", maxv="post6",name="prnt20_any_post6",func="any",diff=FALSE)
robs43[,robs:=pcr_pos_any_preg|igm_zero_any_preg|prnt_zero_any_preg|(igm_pos_any_preg&(prnt20_any_preg|prnt_zero_any_post6|prnt20_any_post6))]

# Moderated
mod20<-testfunc(datav=robs43,respv="prnt_titer_1000",timev="elisa_ga",minv="preg0", maxv="endga",name="prnt1000_any_preg",func="any",diff=FALSE)
mod21a<-testfunc(datav=mod20,respv="prnt_titer",timev="elisa_ga",minv="preg0", maxv="endga",name="prnt_max_preg",func="max",diff=FALSE)
mod21b<-testfunc(datav=mod21a,respv="prnt_titer",timev="elisa_ga",minv="endga", maxv="post2",name="prnt_max_post2",func="max",diff=FALSE)
mod21b[,prnt_rise_post2:=prnt_max_post2>prnt_max_preg]
mod21b[,prnt_4rise_post2:=prnt_max_post2>4*prnt_max_preg]

mod40a<-testfunc(datav=mod21b,respv="pcr_res2",timev="pcr_ga",minv="preg0", maxv="endga",name="pcr_eq_any_preg",func="any",diff=FALSE)
mod40b<-testfunc(datav=mod40a,respv="igm_res2",timev="elisa_ga",minv="preg0", maxv="endga",name="igm_eq_any_preg",func="any",diff=FALSE)
mod41 <-testfunc(datav=mod40b,respv="prnt_titer_100",timev="elisa_ga",minv="endga", maxv="post3",name="prnt100_any_post3",func="any",diff=FALSE)
mod41[,mod:=igm_pos_any_preg|(prnt1000_any_preg&prnt_rise_post2)|prnt_4rise_post2|(prnt100_any_post3&(pcr_eq_any_preg|igm_eq_any_preg))]

# Limited
lim1<-testfunc(datav=mod41,respv="prnt_titer_100",timev="elisa_ga",minv="preg0", maxv="post6",name="prnt100_any",func="any",diff=FALSE)
lim2<-testfunc(datav=lim1,respv="prntz",timev="elisa_ga",minv="endga", maxv="post3",name="prnt_zero_any_post3",func="any",diff=TRUE)
lim2[,lim:=prnt100_any| prnt_zero_any_post3]

# Negative
neg1<-testfunc(datav=lim2,respv="pcr_res0",timev="pcr_ga",minv="preg0", maxv="endga",name="pcr_neg_any_preg",func="any",diff=FALSE)
neg2<-testfunc(datav=neg1,respv="igm_res0",timev="pcr_ga",minv="preg0", maxv="endga",name="igm_neg_any_preg",func="any",diff=FALSE)
finaltest<-testfunc(datav=neg2,respv="pcr_resNA",timev="pcr_ga",minv="preg0", maxv="endga",name="pcr_na_all_preg",func="all",diff=FALSE)
finaltest[,neg:=pcr_neg_any_preg&igm_neg_any_preg|igm_neg_any_preg&pcr_na_all_preg]

#Colapse all conditions in one variable zikv_test_ev
finaltest[,robs:=as.factor(ifelse(is.na(robs),"NA",robs))]
finaltest[,mod:=as.factor(ifelse(is.na(mod),"NA",mod))]
finaltest[,lim:=as.factor(ifelse(is.na(lim),"NA",lim))]
finaltest[,neg:=as.factor(ifelse(is.na(neg),"NA",neg))]
finaltest[,zikv_test_ev:=ifelse(robs=="TRUE","Robust",ifelse(mod=="TRUE","Moderate",ifelse(lim=="TRUE","Limited",ifelse(neg=="TRUE","Negative",NA))))]
childtest<-finaltest[, head(.SD, 1), by = "childid"]
table(childtest$zikv_test_ev, useNA = "always")

nomeasure<-childtest[is.na(zikv_test_ev),]
navalue<-finaltest[childid%in%nomeasure$childid&pcr_res1==1&!is.na(pcr_ga)&!is.na(endga),]

nomeasure[,zikv_test_ev:=ifelse(robs==TRUE,"Robust",ifelse(mod==TRUE,"Moderate",ifelse(lim==TRUE,"Limited",ifelse(neg==TRUE,"Negative",NA))))]
nomeasure[,navalue:=ifelse(robs==TRUE,"T",ifelse(!is.na(robs),"NA","F"))]



#Check why everyone is missing
table(finaltest$zikv_test_ev,useNA = "always")
108380/(108380+31530)
missing<-finaltest[is.na(zikv_test_ev),]
#Are there women in the dataset with moderate evidence?
#Check if there are women who are igm positive, but do not have seroconversion and no other positive tests
length(data[zikv_igm_res_1==1,]) #First igM is positive - n=4045
igmpos<-igmpos[igmpos$zikv_pcr_everpos==0,] #Exclude with positive pcr
igmpos<-igmpos[is.na(igmpos$zikv_prnt_1) | igmpos$zikv_prnt_1==0,] #Exclude positive prnts
summary(as.factor(igmpos$zikv_igm_res_4)) 
