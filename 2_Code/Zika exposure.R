


#Script to construct a new zika exposure variable based on the definitions described in the paper by Ximenes et al.

#Load data
setwd("/Users/jdamen/Documents/Julius/ZIKV analyses/2. Data")
zikv<-read.csv("zikv_033_datasets.csv",header=T)
#Recode missings
zikv[zikv==""] <-NA
zikv[zikv==555] <-  NA
zikv[zikv==666] <-  NA
zikv[zikv==777] <-  NA
zikv[zikv==888] <-  NA
zikv[zikv==999] <-  NA
zikv[zikv==9999] <-  NA

zikv$zikv_ximenes<-NA

#Robust evidence of zikv
#Positive PCR
zikv$zikv_ximenes[zikv$zikv_pcr_res_1==1]<-1
zikv$zikv_ximenes[zikv$zikv_pcr_res_2==1]<-1
zikv$zikv_ximenes[zikv$zikv_pcr_res_3==1]<-1
zikv$zikv_ximenes[zikv$zikv_pcr_res_4==1]<-1
zikv$zikv_ximenes[zikv$zikv_pcr_res_5==1]<-1
zikv$zikv_ximenes[zikv$zikv_pcr_res_6==1]<-1
zikv$zikv_ximenes[zikv$zikv_pcr_res_7==1]<-1
zikv$zikv_ximenes[zikv$zikv_pcr_res_8==1]<-1
zikv$zikv_ximenes[zikv$zikv_pcr_everpos==1]<-1
#Seroconversion of IgM
zikv$zikv_ximenes[zikv$zikv_igm_res_1==0 & zikv$zikv_igm_res_2==1 & !is.na(zikv$zikv_igm_res_1) & !is.na(zikv$zikv_igm_res_2)]<-1
zikv$zikv_ximenes[zikv$zikv_igm_res_1==0 & zikv$zikv_igm_res_3==1 & !is.na(zikv$zikv_igm_res_1) & !is.na(zikv$zikv_igm_res_3)]<-1
zikv$zikv_ximenes[zikv$zikv_igm_res_2==0 & zikv$zikv_igm_res_3==1 & !is.na(zikv$zikv_igm_res_2) & !is.na(zikv$zikv_igm_res_3)]<-1
#Seroconversion of PRNT
zikv$zikv_ximenes[zikv$zikv_prnt_1==0 & zikv$zikv_prnt_2==1 & !is.na(zikv$zikv_prnt_1) & !is.na(zikv$zikv_prnt_2)]<-1
zikv$zikv_ximenes[zikv$zikv_prnt_1==0 & zikv$zikv_prnt_3==1 & !is.na(zikv$zikv_prnt_1) & !is.na(zikv$zikv_prnt_3)]<-1
zikv$zikv_ximenes[zikv$zikv_prnt_1==0 & zikv$zikv_prnt_4==1 & !is.na(zikv$zikv_prnt_1) & !is.na(zikv$zikv_prnt_4)]<-1
zikv$zikv_ximenes[zikv$zikv_prnt_1==0 & zikv$zikv_prnt_5==1 & !is.na(zikv$zikv_prnt_1) & !is.na(zikv$zikv_prnt_5)]<-1
zikv$zikv_ximenes[zikv$zikv_prnt_1==0 & zikv$zikv_prnt_6==1 & !is.na(zikv$zikv_prnt_1) & !is.na(zikv$zikv_prnt_6)]<-1
zikv$zikv_ximenes[zikv$zikv_prnt_1==0 & zikv$zikv_prnt_7==1 & !is.na(zikv$zikv_prnt_1) & !is.na(zikv$zikv_prnt_7)]<-1
zikv$zikv_ximenes[zikv$zikv_prnt_1==0 & zikv$zikv_prnt_8==1 & !is.na(zikv$zikv_prnt_1) & !is.na(zikv$zikv_prnt_8)]<-1
zikv$zikv_ximenes[zikv$zikv_prnt_2==0 & zikv$zikv_prnt_3==1 & !is.na(zikv$zikv_prnt_2) & !is.na(zikv$zikv_prnt_3)]<-1
zikv$zikv_ximenes[zikv$zikv_prnt_2==0 & zikv$zikv_prnt_4==1 & !is.na(zikv$zikv_prnt_2) & !is.na(zikv$zikv_prnt_4)]<-1
zikv$zikv_ximenes[zikv$zikv_prnt_2==0 & zikv$zikv_prnt_5==1 & !is.na(zikv$zikv_prnt_2) & !is.na(zikv$zikv_prnt_5)]<-1
zikv$zikv_ximenes[zikv$zikv_prnt_2==0 & zikv$zikv_prnt_6==1 & !is.na(zikv$zikv_prnt_2) & !is.na(zikv$zikv_prnt_6)]<-1
zikv$zikv_ximenes[zikv$zikv_prnt_2==0 & zikv$zikv_prnt_7==1 & !is.na(zikv$zikv_prnt_2) & !is.na(zikv$zikv_prnt_7)]<-1
zikv$zikv_ximenes[zikv$zikv_prnt_2==0 & zikv$zikv_prnt_8==1 & !is.na(zikv$zikv_prnt_2) & !is.na(zikv$zikv_prnt_8)]<-1
zikv$zikv_ximenes[zikv$zikv_prnt_3==0 & zikv$zikv_prnt_4==1 & !is.na(zikv$zikv_prnt_3) & !is.na(zikv$zikv_prnt_4)]<-1
zikv$zikv_ximenes[zikv$zikv_prnt_3==0 & zikv$zikv_prnt_5==1 & !is.na(zikv$zikv_prnt_3) & !is.na(zikv$zikv_prnt_5)]<-1
zikv$zikv_ximenes[zikv$zikv_prnt_3==0 & zikv$zikv_prnt_6==1 & !is.na(zikv$zikv_prnt_3) & !is.na(zikv$zikv_prnt_6)]<-1
zikv$zikv_ximenes[zikv$zikv_prnt_3==0 & zikv$zikv_prnt_7==1 & !is.na(zikv$zikv_prnt_3) & !is.na(zikv$zikv_prnt_7)]<-1
zikv$zikv_ximenes[zikv$zikv_prnt_3==0 & zikv$zikv_prnt_8==1 & !is.na(zikv$zikv_prnt_3) & !is.na(zikv$zikv_prnt_8)]<-1
zikv$zikv_ximenes[zikv$zikv_prnt_4==0 & zikv$zikv_prnt_5==1 & !is.na(zikv$zikv_prnt_4) & !is.na(zikv$zikv_prnt_5)]<-1
zikv$zikv_ximenes[zikv$zikv_prnt_4==0 & zikv$zikv_prnt_6==1 & !is.na(zikv$zikv_prnt_4) & !is.na(zikv$zikv_prnt_6)]<-1
zikv$zikv_ximenes[zikv$zikv_prnt_4==0 & zikv$zikv_prnt_7==1 & !is.na(zikv$zikv_prnt_4) & !is.na(zikv$zikv_prnt_7)]<-1
zikv$zikv_ximenes[zikv$zikv_prnt_4==0 & zikv$zikv_prnt_8==1 & !is.na(zikv$zikv_prnt_4) & !is.na(zikv$zikv_prnt_8)]<-1
zikv$zikv_ximenes[zikv$zikv_prnt_5==0 & zikv$zikv_prnt_6==1 & !is.na(zikv$zikv_prnt_5) & !is.na(zikv$zikv_prnt_6)]<-1
zikv$zikv_ximenes[zikv$zikv_prnt_5==0 & zikv$zikv_prnt_7==1 & !is.na(zikv$zikv_prnt_5) & !is.na(zikv$zikv_prnt_7)]<-1
zikv$zikv_ximenes[zikv$zikv_prnt_5==0 & zikv$zikv_prnt_8==1 & !is.na(zikv$zikv_prnt_5) & !is.na(zikv$zikv_prnt_8)]<-1
zikv$zikv_ximenes[zikv$zikv_prnt_6==0 & zikv$zikv_prnt_7==1 & !is.na(zikv$zikv_prnt_6) & !is.na(zikv$zikv_prnt_7)]<-1
zikv$zikv_ximenes[zikv$zikv_prnt_6==0 & zikv$zikv_prnt_8==1 & !is.na(zikv$zikv_prnt_6) & !is.na(zikv$zikv_prnt_8)]<-1
zikv$zikv_ximenes[zikv$zikv_prnt_7==0 & zikv$zikv_prnt_8==1 & !is.na(zikv$zikv_prnt_7) & !is.na(zikv$zikv_prnt_8)]<-1



summary(as.factor(zikv$zikv_ximenes))


names(zikv2[ , grepl( "zikv_assay" , names(zikv2) ) ])
