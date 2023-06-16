
library(ggplot2)
library(naniar)
library(xlsx)

setwd("/Users/jdamen/Documents/Julius/ZIKV analyses/2. Data")
zikv<-read.csv("zikv_033_datasets.csv",header=T)
#zikv<-read.xlsx("zikv_033_datasets.xlsx", header=T,sheetIndex = 1)

names(zikv)
summary(zikv)

zikv[zikv==""] <-NA
zikv[zikv==555] <-  NA
zikv[zikv==666] <-  NA
zikv[zikv==777] <-  NA
zikv[zikv==888] <-  NA
zikv[zikv==999] <-  NA
zikv[zikv==9999] <-  NA

#Table with percentage of missings per variable
gg_miss_var(zikv, show_pct = TRUE)

missings<-colSums(is.na(zikv)/length(zikv$mid_original)*100) 
missings<-as.data.frame(missings)

#Delete columns without any data
zikv2<-zikv[, which(colMeans(is.na(zikv)) < 1)]
#Deleted columns:
#[1] "n_visit_fu_child_study"    "ch_visit_study_n"          "zikv_igm_pnratio_1"        "zikv_igg_titer_1"          "zikv_igg_pnratio_1"       
#[6] "zikvdenv_prnt_titerdiff_1" "fever_meas_1"              "denv_igm_titer_1"          "chikv_igm_titer_1"         "chikv_igg_titer_1"        
#[11] "treat_fertil"              "ch_chromoabn_test"         "ch_chromoabn_test_oth"     "uti_n"                     "genitalwarts"             
#[16] "ch_rash_dur_1"             "zikv_assay_ga_10"          "zikv_assay_tri_2"          "zikv_assay_tri_3"          "zikv_assay_tri_4"         
#[21] "zikv_assay_tri_5"          "zikv_assay_tri_6"          "zikv_assay_tri_7"          "zikv_assay_tri_8"          "zikv_assay_tri_9"         
#[26] "zikv_assay_tri_10"         "zikv_pcr_ga_9"             "zikv_pcr_ga_10"            "zikv_pcr_tri_3"            "zikv_pcr_tri_4"           
#[31] "zikv_pcr_tri_5"            "zikv_pcr_tri_6"            "zikv_pcr_tri_7"            "zikv_pcr_tri_8"            "zikv_pcr_tri_9"           
#[36] "zikv_pcr_tri_10"           "zikv_pcr_res_9"            "zikv_pcr_res_10"           "zikv_elisa_ga_4"           "zikv_elisa_ga_5"          
#[41] "zikv_elisa_ga_6"           "zikv_elisa_ga_7"           "zikv_elisa_ga_8"           "zikv_elisa_ga_9"           "zikv_elisa_ga_10"         
#[46] "zikv_elisa_tri_4"          "zikv_elisa_tri_5"          "zikv_elisa_tri_6"          "zikv_elisa_tri_7"          "zikv_elisa_tri_8"         
#[51] "zikv_elisa_tri_9"          "zikv_elisa_tri_10"         "zikv_elisa_res_10"         "zikv_igm_res_4"            "zikv_igm_res_5"           
#[56] "zikv_igm_res_6"            "zikv_igm_res_7"            "zikv_igm_res_8"            "zikv_igm_res_9"            "zikv_igm_res_10"          
#[61] "zikv_igg_titer_2"          "zikv_prnt_9"               "zikv_prnt_10" 

#write.csv(missings,file="/Users/jdamen/Documents/Julius/ZIKV analyses/4. Resultaten/20221018 Percentage missings.csv")

names(zikv2)
names(zikv2[ , grepl( "zikv_" , names(zikv2) ) ]) #Overview of variables related to exposure
summary(zikv$weight)

summary(as.factor(zikv$mid))
levels(as.factor(zikv$mid))
summary(as.factor(zikv$mid_original))
levels(as.factor(zikv$mid_original))

#Work with dates
as.Date(as.numeric(zikv$zikv_assay_date_1), origin = "1899-12-30")


