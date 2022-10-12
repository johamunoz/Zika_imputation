# Data manipulation packages
library(here) 
library(data.table)
library(rio)


# ZIKV Test according to Ricardo's paper ---

ziktest_sum<-function(datav,respv,timev, minv, maxv,name,diff,func){
 
  #' Function to summarize at individual level test related variables in a specific period of time
  #' Assumes that childID have a unique ID across all the included studies in the IPD
  #'
  #' @param datav original dataset 
  #' @param respv response variable to summarize
  #' @param timev measurment time of response variable (given in continuos time,i.e. gestational age)
  #' @param minv starting time of time frame
  #' @param maxv end time of time frame 
  #' @param name name given to the return variable
  #' @param diff if calculates lag difference or not (TRUE, FALSE)
  #' @param func specify type  of summary ("all", "any")
  #' @return dataframe that includes new summarized variable 
  
  datav[,min:= get(minv)]
  datav[,max:= get(maxv)]
  datav[,time:=get(timev)]
  datav[,resp:=get(respv)]
  dataf<-datav[time<=max&time>=min&!is.na(resp),]
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


ziktest_ml <- function(data){

  #' Function that returns zikv_test_ev,i.e. a zika test variable with multiple evidence levels ("robust","moderate","limited","negative") 
  #' Refer to Ricardo paper, ensure that end_ga have been defined it before running the function
  #'
  #' @param data original dataset 
  #' @return dataframe that includes new ziktest_ml variable and child id that will be merged with full dataset
  

    # Set dataset with only the zika related variables 
    var_test<-c("childid","end_ga","denv_igm_res_1",paste0("zikv_pcr_ga_",1:10),paste0("zikv_pcr_res_",1:10),
                paste0("zikv_pcr_tri_",1:10),paste0("zikv_pcr_date_",1:10), 
                paste0("zikv_igm_res_",1:10),paste0("zikv_elisa_ga_",1:10),
                paste0("zikv_prnt_titer_",1:2),
                paste0("zikv_prnt_",1:10))
    
    
    data_test <- as.data.table(data[,..var_test])
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
    
    # Calculate additional required variables 
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
    data_test_m[, end_ga:=as.numeric(end_ga)]
    
    data_test_m[, preg0:=0]
    data_test_m[, post1:=end_ga+4]
    data_test_m[, post2:=end_ga+8]
    data_test_m[, post3:=end_ga+13]
    data_test_m[, post6:=end_ga+26]
    
    # Robust
    robs1<-ziktest_sum(datav=data_test_m,respv="pcr_res1",timev="pcr_ga",minv="preg0", maxv="end_ga",name="pcr_pos_any_preg",func="any",diff=FALSE)
    robs2<-ziktest_sum(datav=robs1,respv="igmz",timev="elisa_ga",minv="preg0", maxv="end_ga",name="igm_zero_any_preg",func="any",diff=TRUE)
    robs3<-ziktest_sum(datav=robs2,respv="prntz",timev="elisa_ga",minv="preg0", maxv="end_ga",name="prnt_zero_any_preg",func="any",diff=TRUE)
    robs40<-ziktest_sum(datav=robs3,respv="igmz",timev="elisa_ga",minv="preg0", maxv="end_ga",name="igm_pos_any_preg",func="any",diff=FALSE)
    robs41<-ziktest_sum(datav=robs40,respv="prnt_titer_20",timev="elisa_ga",minv="preg0", maxv="end_ga",name="prnt20_any_preg",func="any",diff=FALSE)
    robs42<-ziktest_sum(datav=robs41,respv="prntz",timev="elisa_ga",minv="end_ga", maxv="post6",name="prnt_zero_any_post6",func="any",diff=TRUE)
    robs43<-ziktest_sum(datav=robs42,respv="prnt_titer_20",timev="elisa_ga",minv="end_ga", maxv="post6",name="prnt20_any_post6",func="any",diff=FALSE)
    robs43[,robs:=pcr_pos_any_preg|igm_zero_any_preg|prnt_zero_any_preg|(igm_pos_any_preg&(prnt20_any_preg|prnt_zero_any_post6|prnt20_any_post6))]
    
    # Moderated
    mod20<-ziktest_sum(datav=robs43,respv="prnt_titer_1000",timev="elisa_ga",minv="preg0", maxv="end_ga",name="prnt1000_any_preg",func="any",diff=FALSE)
    mod21a<-ziktest_sum(datav=mod20,respv="prnt_titer",timev="elisa_ga",minv="preg0", maxv="end_ga",name="prnt_max_preg",func="max",diff=FALSE)
    mod21b<-ziktest_sum(datav=mod21a,respv="prnt_titer",timev="elisa_ga",minv="end_ga", maxv="post2",name="prnt_max_post2",func="max",diff=FALSE)
    mod21b[,prnt_rise_post2:=prnt_max_post2>prnt_max_preg]
    mod21b[,prnt_4rise_post2:=prnt_max_post2>4*prnt_max_preg]
    
    mod40a<-ziktest_sum(datav=mod21b,respv="pcr_res2",timev="pcr_ga",minv="preg0", maxv="end_ga",name="pcr_eq_any_preg",func="any",diff=FALSE)
    mod40b<-ziktest_sum(datav=mod40a,respv="igm_res2",timev="elisa_ga",minv="preg0", maxv="end_ga",name="igm_eq_any_preg",func="any",diff=FALSE)
    mod41 <-ziktest_sum(datav=mod40b,respv="prnt_titer_100",timev="elisa_ga",minv="end_ga", maxv="post3",name="prnt100_any_post3",func="any",diff=FALSE)
#TODO @Anneke check if this is correct    
    # Last condition is not specified in the paper and is related to dengue test results interaction. 
    mod41[,mod:=igm_pos_any_preg|(prnt1000_any_preg&prnt_rise_post2)|prnt_4rise_post2|(prnt100_any_post3&(pcr_eq_any_preg|(igm_eq_any_preg&(denv_igm_res_1!=1|is.na(denv_igm_res_1)))))]
   
    # Limited
    lim1<-ziktest_sum(datav=mod41,respv="prnt_titer_100",timev="elisa_ga",minv="preg0", maxv="post6",name="prnt100_any",func="any",diff=FALSE)
    lim2<-ziktest_sum(datav=lim1,respv="prntz",timev="elisa_ga",minv="end_ga", maxv="post3",name="prnt_zero_any_post3",func="any",diff=TRUE)
    # Last condition is not specified in the paper and is related to dengue test results interaction. 
    lim2[,lim:=prnt100_any| prnt_zero_any_post3|(igm_eq_any_preg&denv_igm_res_1==1)]
    
    # Negative
    neg1<-ziktest_sum(datav=lim2,respv="pcr_res0",timev="pcr_ga",minv="preg0", maxv="end_ga",name="pcr_neg_any_preg",func="any",diff=FALSE)
    neg2<-ziktest_sum(datav=neg1,respv="igm_res0",timev="pcr_ga",minv="preg0", maxv="end_ga",name="igm_neg_any_preg",func="any",diff=FALSE)
    finaltest<-ziktest_sum(datav=neg2,respv="pcr_resNA",timev="pcr_ga",minv="preg0", maxv="end_ga",name="pcr_na_all_preg",func="all",diff=FALSE)
    finaltest[,neg:=pcr_neg_any_preg&igm_neg_any_preg|igm_neg_any_preg&pcr_na_all_preg]
    
    #Colapse all conditions in one variable zikv_test_ev
    finaltest[,robs:=as.factor(ifelse(is.na(robs),"NA",robs))]
    finaltest[,mod:=as.factor(ifelse(is.na(mod),"NA",mod))]
    finaltest[,lim:=as.factor(ifelse(is.na(lim),"NA",lim))]
    finaltest[,neg:=as.factor(ifelse(is.na(neg),"NA",neg))]
    finaltest[,zikv_test_ev:=ifelse(robs=="TRUE","Robust",ifelse(mod=="TRUE","Moderate",ifelse(lim=="TRUE","Limited",ifelse(neg=="TRUE","Negative",NA))))]
    
    childtest<-finaltest[, .SD[c(1)]$zikv_test_ev, by ="childid"]
    colnames(childtest)<-c("childid","zikv_test_ev")
    childtest[,zikv_test_ev:=as.factor(zikv_test_ev)]
    
 return(childtest)
}


# Functions related to CZN correction---

checkcon<-function(data,setcol){ 
  #' Function that summarize if there is any abnormality in a set of variables for a given child
  #' Returns 1: abnormality detected on any of the set of variables.
  #'         0: no abnormality detected across the set of variables.
  #'         NA: no information available on any variable.
  #' @param data original data set 
  #' @param col1 set of variables across which the detection is summarized
  #' @return detection of abnormality in set of variables. 
  
  data[,(setcol):= lapply(.SD, as.numeric), .SDcols = setcol]
  data$anyT <- rowSums(data[, .SD, .SDcols = setcol], na.rm=T)
  data$allNA <- rowSums(!is.na(data[, .SD, .SDcols = setcol]))
  data$Final<- ifelse(data$allNA==0,NA,ifelse(data$anyT>0,1,0))
  return(data$Final)
}


# Function postnatal microcephaly

micro_postnatal<-function(data){ 
  #' Function that summarize if there is any postnatal microcephaly diagnosis 
  #' based on the head circumference
  #' We used igprepost_hcircm2zscore when birth was <=36 weeks and who_hcircm2zscore otherwise
  #' when the functions return inf we set the zcore to NA
  #' @param data original dataset 
  #' @return original dataset + micro_postnatal binary variable
  
  micr_test<-c("childid","ch_sex","fet_death","end_ga",paste0("ch_head_circ_",1:20),paste0("ch_head_circ_age_",1:20))
  data_micr<-data[,..micr_test]
  colA = paste("ch_head_circ_", 1:20, sep = "")
  colB = paste("ch_head_circ_age_", 1:20, sep = "")
  data_micr_m = melt(data_micr, measure = list(colA,colB), 
                     value.name = c("circ", "age"),
                     variable.name="visit")
  
  data_micr_m[,sex:=ifelse(is.na(ch_sex),NA,ifelse(ch_sex== 0, "Male","Female"))]
  data_micr_m[,days:=age*30]
  data_micr_m[!is.na(sex),zscorepreterm:=igprepost_hcircm2zscore(pmagedays=days, hcircm=circ, sex =sex)]
  data_micr_m[!is.na(sex),zscorenorm:=who_hcircm2zscore(agedays=days, hcircm=circ, sex =sex)]
  data_micr_m[,zscore:=ifelse(end_ga<=36,zscorepreterm,zscorenorm)]

  data_micr_m[,zscore:=ifelse(is.infinite(zscore),NA,zscore)]
  data_micr_m[,micro:= ifelse(zscore<=-3,2,ifelse(zscore<=-2,1,ifelse(zscore<=2,0,ifelse(!is.na(zscore),3,NA))))]
  data_micr_m[,micro_bin:=ifelse(micro%in%c(0,3),0,ifelse(micro%in%c(1,2),1,NA))]
  
  post_micro<-data_micr_m[ , max(micro_bin,na.rm=TRUE), by = childid]
  colnames(post_micro)<-c("childid","microcephaly_bin_postnatal")
  post_micro[,microcephaly_bin_postnatal:=ifelse(microcephaly_bin_postnatal=="-Inf",NA,microcephaly_bin_postnatal)]
  data<-as.data.table(merge(data,post_micro,by="childid",all.x = TRUE))
  
  return(data)
}


