rm(list=ls()) # clean environment

# Data manipulation packages
library(here) 
library(data.table)
library(dplyr)
library(haven)
require(ggplot2)

data <- as.data.table(readxl::read_xlsx(here('1_Input_data','zikv_033_datasets.xlsx'),sheet="Sheet1"))
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
data[data==777] <-  NA
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


robs1<-testfunc(datav=data_test_m,respv="pcr_res1",timev="pcr_ga",minv="preg0", maxv="endga",name="pcr_pos_any_preg",func="any",diff=FALSE)
robs2<-testfunc(datav=robs1,respv="igmz",timev="elisa_ga",minv="preg0", maxv="endga",name="igm_zero_any_preg",func="any",diff=TRUE)
robs3<-testfunc(datav=robs2,respv="prntz",timev="elisa_ga",minv="preg0", maxv="endga",name="prnt_zero_any_preg",func="any",diff=TRUE)
robs40<-testfunc(datav=robs3,respv="igmz",timev="elisa_ga",minv="preg0", maxv="endga",name="igm_pos_any_preg",func="any",diff=FALSE)
robs41<-testfunc(datav=robs40,respv="prnt_titer_20",timev="elisa_ga",minv="preg0", maxv="endga",name="prnt20_any_preg",func="any",diff=FALSE)
robs42<-testfunc(datav=robs41,respv="prntz",timev="elisa_ga",minv="endga", maxv="post6",name="prnt_zero_any_post6",func="any",diff=TRUE)
robs43<-testfunc(datav=robs42,respv="prnt_titer_20",timev="elisa_ga",minv="endga", maxv="post6",name="prnt20_any_post6",func="any",diff=FALSE)
robs43[,robs:=pcr_pos_any_preg|igm_zero_any_preg|prnt_zero_any_preg|(igm_pos_any_preg&(prnt20_any_preg|prnt_zero_any_post6|prnt20_any_post6))]

mod20<-testfunc(datav=robs43,respv="prnt_titer_1000",timev="elisa_ga",minv="preg0", maxv="endga",name="prnt1000_any_preg",func="any",diff=FALSE)
mod21a<-testfunc(datav=mod20,respv="prnt_titer",timev="elisa_ga",minv="preg0", maxv="endga",name="prnt_max_preg",func="max",diff=FALSE)
mod21b<-testfunc(datav=mod21a,respv="prnt_titer",timev="elisa_ga",minv="endga", maxv="post2",name="prnt_max_post2",func="max",diff=FALSE)
mod21b[,prnt_rise_post2:=prnt_max_post2>prnt_max_preg]
mod21b[,prnt_4rise_post2:=prnt_max_post2>4*prnt_max_preg]

mod40a<-testfunc(datav=mod21b,respv="pcr_res2",timev="pcr_ga",minv="preg0", maxv="endga",name="pcr_eq_any_preg",func="any",diff=FALSE)
mod40b<-testfunc(datav=mod40a,respv="igm_res2",timev="elisa_ga",minv="preg0", maxv="endga",name="igm_eq_any_preg",func="any",diff=FALSE)
mod41 <-testfunc(datav=mod40b,respv="prnt_titer_100",timev="elisa_ga",minv="endga", maxv="post3",name="prnt100_any_post3",func="any",diff=FALSE)
mod41[,mod:=igm_pos_any_preg|(prnt1000_any_preg&prnt_rise_post2)|prnt_4rise_post2|(prnt100_any_post3&(pcr_eq_any_preg|igm_eq_any_preg))]

lim1<-testfunc(datav=mod41,respv="prnt_titer_100",timev="elisa_ga",minv="preg0", maxv="post6",name="prnt100_any",func="any",diff=FALSE)
lim2<-testfunc(datav=lim1,respv="prntz",timev="elisa_ga",minv="endga", maxv="post3",name="prnt_zero_any_post3",func="any",diff=TRUE)
lim2[,lim:=prnt100_any| prnt_zero_any_post3]

neg1<-testfunc(datav=lim2,respv="pcr_res0",timev="pcr_ga",minv="preg0", maxv="endga",name="pcr_neg_any_preg",func="any",diff=FALSE)
neg2<-testfunc(datav=neg1,respv="igm_res0",timev="pcr_ga",minv="preg0", maxv="endga",name="igm_neg_any_preg",func="any",diff=FALSE)
finaltest<-testfunc(datav=neg2,respv="pcr_resNA",timev="pcr_ga",minv="preg0", maxv="endga",name="pcr_na_all_preg",func="all",diff=FALSE)
finaltest[,neg:=pcr_neg_any_preg&igm_neg_any_preg|igm_neg_any_preg&pcr_na_all_preg]
finaltest[,zikv_test_ev:=ifelse(robs==TRUE,"Robust",ifelse(mod==TRUE,"Moderate",ifelse(lim==TRUE,"Limited",ifelse(neg==TRUE,"Negative",NA))))]
finaltest[,zikv_test_ev:=as.factor(zikv_test_ev)]


table(finaltest$robs,finaltest$lim,useNA = "always")
robs1<-testfunc(datav=data_test_m,variable=pcr)
library(tidyr)
library(dplyr)

data_visit<-data_test%>%
  gather(key="full_name",value,zikv_pcr_ga_1:zikv_igm_titer_3)%>%
  separate(full_name, c("variable", "visit"), sep = "_(?!.*_)")%>%
  dcast( mid+endga+visit ~ variable, value.var = "value")
dw %>% 
  gather(v, value, f1.avg:f2.sd)
testout<-function(endga,varvec,timevec,period,transf){
  
}






infoexp <- as.data.table(readxl::read_xlsx(here('1_Input_data','MasterCodebook_Final_June2022 ALL (Repaired).xlsx'),sheet="237 key")) #CSV file with the




var_inc<-as.matrix(infoexp[Essential=="Yes",'WHO variable name'])[,1] #Variables to work with




# date of event (when pregnancy ends)



data<-as.data.table(data[,..var_inc]) #Filter dataset



#3 Test data


data_test<-as.data.table(data[,..var_test]) # data_test





data_test[zikv_pcr_ga_2==TRUE]
data_test$zikv_pcr_res_1 
unique(data$zikv_pcr_res_2)

unique(data$zikv_igm_res_1)
unique(data$zikv_igm_res_2)
unique(data$zikv_igm_res_3)
unique(data$zikv_igm_res_4)
summary(data_test)
unique(data_test$zikv_igm_date_1)
summary(data_test$zikv_pcr_date_2)
data_test[zikv_pcr_ga_2,c("zikv_pcr_ga_2"):=list(NA)]
zikv_pcr_ga_9
zikv_pcr_ga_10
zikv_pcr_res_1
zikv_pcr_res_2
summary(data)
#4. Construct Bdeath (baby death variable) -----

#0. Put miscarriage, loss, loss_etiology, birth in one variable bdeath and bdeath_ga---

data[,birth:=ifelse(!is.na(birth_ga),1,NA)] #birth indicator from bithga
data[,birth2:=ifelse(!is.na(ch_term),1,NA)] #birth indicator from ch_term if no number is specified
data[,bdeath:=ifelse(birth==1|birth2==1,0,NA)]
data[,bdeath_ga:=ifelse(!is.na(loss_ga),loss_ga,ifelse(!is.na(miscarriage_ga),miscarriage_ga,NA))]

data[,endga:=ifelse(!is.na(endga),endga,ifelse(!is.na(birth_ga),birth_ga,bdeath_ga))]
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


# 1.1. Check format of categorical variables----
char <-   which(sapply(data, is.character)) #which are character
datac<-data[,..char]

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

#4. Construct Bdeath (baby death variable) -----


#0. Put miscarriage, loss, loss_etiology, birth in one variable bdeath and bdeath_ga---
data[,birth:=ifelse(!is.na(birth_ga),1,NA)] #birth indicator from bithga
data[,birth2:=ifelse(!is.na(ch_term),1,NA)] #birth indicator from inf_term=ch_term?
data[,bdeath:=ifelse(birth==1|birth2==1,0,NA)]  # endga?
data[,bdeath_ga:=ifelse(!is.na(loss_ga),loss_ga,ifelse(!is.na(miscarriage_ga),miscarriage_ga,NA))]
data[,end_ga:=ifelse(!is.na(birth_ga),birth_ga,bdeath_ga)]
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




