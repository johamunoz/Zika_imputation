###Aim: To describe better the gestational age variables

rm(list=ls()) # clean environment

# Load packages ---
# Data manipulation package
library(data.table) 
library(here)  # define folder paths
library(rio)   # import STATA file
library(stringr)

# Field specific package
library(growthstandards)

# Graphic packages
library(ggplot2)
library(plotly)
library(mice)


# Load dataset and dependencies ----
data_origin <- as.data.table(import(here('1_Input_data','zikv_033_datasets.dta'))) #  This file you can find it on the dropbox
data_origin[,studyname:=NULL] # as studyname was not assigned to all studies
add_info <- as.data.table(readxl::read_xlsx(here('1_Input_data','MasterCodebook_October.xlsx'),sheet="237 key")) #CSV file with the
study_info <- as.data.table(readxl::read_xlsx(here('1_Input_data','MasterCodebook_October.xlsx'),sheet="StudyID")) #CSV file with the
data <- merge(data_origin, study_info, by="file")
source(here('2.1_Code 33 studies','1.1_Pre_imputation_functions.R'))


# 0. Initial checks ----
# 0.1 Missing data observations----
data[data==""] <-   NA
data[data==555] <-  NA
data[data==666] <-  NA
data[data==888] <-  NA
data[data==999] <-  NA
data[data==9999] <- NA
data[data==999] <-  NA
data[data==777] <-  NA 
data[data=="99/99/9999"]<-NA
data[data=="NA"]<-NA

# 0.2 Check variables class
# Number of observations= # one childID per observation
nrow(data) == length(unique(data$childid))

# Check columns classes

# Check continuous variables

cont_var <- add_info[Type_var == "Continuous"& Var_ga=="yes",]$who_name
cont_var<- cont_var[!cont_var%in%c("fet_death_ga", "end_ga", "microcephaly_ga", "date_t1_ga","zikv_ga_n")]
type_ga <- sapply(data[,..cont_var], class)
var_ga <- names(type_ga[type_ga=="character"]) # Continuous variables misclassified as character
data[, (var_ga) := lapply(.SD, as.numeric), .SDcols = var_ga]
cont_bound <- cont_bound(add_info,data)

## Correct unplausible values of ga
data[ arb_clindiag_ga<0,c("childid","arb_clindiag_ga")]# 3 observations
data[	zikv_assay_ga_1<0,c("childid","zikv_assay_ga_1")] # 7 observations
data[	zikv_elisa_ga_1<0,c("childid","zikv_elisa_ga_1")]# 7 observations
data[	zikv_ga<0,c("childid","zikv_ga")] # 1 observations
data[	zikv_pcr_ga_1<0,c("childid","zikv_pcr_ga_1")]# 7 observations
data[	zikv_symp_ga<0,c("childid","zikv_symp_ga")]# 14 observations
data[	zikv_symp_test_1<0,c("childid","zikv_symp_test_1")]# 9observations
  # majority of the outlier values are in the study 004 # Ilich check this observations if you think this is convenient, we assigned to NA

data[,arb_clindiag_ga:=ifelse(arb_clindiag_ga<0.001|arb_clindiag_ga>46,NA,arb_clindiag_ga)]
data[,fet_micro_diag_ga:=ifelse(fet_micro_diag_ga<0.001|fet_micro_diag_ga>46,NA,fet_micro_diag_ga)]
data[,loss_ga:=ifelse(loss_ga<20|loss_ga>46,NA,loss_ga)]
data[,miscarriage_ga:=ifelse(miscarriage_ga<0.001|miscarriage_ga>20,NA,miscarriage_ga)]
data[,zikv_assay_ga_1:=ifelse(zikv_assay_ga_1<0.001|zikv_assay_ga_1>46,NA,zikv_assay_ga_1)]
data[,zikv_elisa_ga_1:=ifelse(zikv_elisa_ga_1<0.001|zikv_elisa_ga_1>46,NA,zikv_elisa_ga_1)]
data[,zikv_ga:=ifelse(zikv_ga<0.001|zikv_ga>46,NA,zikv_ga)]
data[,zikv_pcr_ga_1:=ifelse(zikv_pcr_ga_1<0.001|zikv_pcr_ga_1>46,NA,zikv_pcr_ga_1)]
data[,zikv_symp_ga:=ifelse(zikv_symp_ga<0.001|zikv_symp_ga>46,NA,zikv_symp_ga)]
data[,zikv_symp_test_1:=ifelse(zikv_symp_test_1<0.001|zikv_symp_test_1>46,NA,zikv_symp_test_1)]


# Check date variables
date_var <- add_info[Type_var == "Date"& Var_ga=="yes",]$who_name
data[, (date_var) := lapply(.SD, as.Date, format="%d %b %Y"), .SDcols = date_var]

# Convert dates into gestational age (weeks) if it is possible to get an approx of the conception date
data[,conc_symp:=zikv_symp_date-zikv_symp_ga*7]
data[,conc_assay:=zikv_assay_date_1-zikv_assay_ga_1*7]
data[,conc_pcr:=zikv_pcr_date_1-zikv_pcr_ga_1*7]
data[,conc_elisa:=zikv_elisa_date_1-zikv_elisa_ga_1*7]
data[,conc_date:=apply(data[,c("conc_symp","conc_assay","conc_pcr","conc_elisa")], 1, min, na.rm = TRUE)] #get the min as date of conception


# add an proxy ga  to ga values with NA based on date reported and approx date of conception
date_concep<-function(var_ga,var_date){
  min<-difftime(data[,get(var_date)],data[,conc_date],units="weeks")
  min<-ifelse(min>0&min<42,min,NA)
  value<-ifelse(is.na(data[,get(var_ga)]),min,data[,get(var_ga)])
}

data[,date_frep:=apply(data[,c("date_t1","date_t2","date_t3","date_t4","date_t5","date_t6","date_t7")], 1, min, na.rm = TRUE)] # first date reported of mothers visit
data[,ch_date_frep:=apply(data[,c("ch_date_t1","ch_date_t2","ch_date_t3","ch_date_t4","ch_date_t5","ch_date_t6","ch_date_t7")], 1, min, na.rm = TRUE)] # first date reported of child visit
data[, zikv_symp_ga:=ifelse(!is.na(zikv_symp_ga),zikv_symp_ga,ifelse(zikv_symp_tri==0,13,ifelse(zikv_symp_tri==1,27,ifelse(zikv_symp_tri==2,42,NA))))] # we aprox NA to the max boundary of the trimester
data[, zikv_symp_ga:=date_concep(var_ga="zikv_symp_ga",var_date="zikv_symp_date")] # we further approx with date
data[, zikv_assay_ga_1:=date_concep(var_ga="zikv_assay_ga_1",var_date="zikv_assay_date_1")]
data[, zikv_pcr_ga_1:=date_concep(var_ga="zikv_pcr_ga_1",var_date="zikv_pcr_date_1")]
data[, zikv_elisa_ga_1:=date_concep(var_ga="zikv_elisa_ga_1",var_date="zikv_elisa_date_1")]
data$date_t1_ga<-NA
data[, date_t1_ga:=date_concep(var_ga="date_t1_ga",var_date="date_frep")] #variable of first time visit on ga.

# Check binary variables
bin_var <- add_info[Type_var == "Binary"& Var_ga=="yes",]$who_name
data[, (bin_var) := lapply(.SD, as.factor), .SDcols = bin_var]
summary(data[,..bin_var])
data[,fet_micro_diag_trin:=ifelse(fet_us_micro_tri2==1,1,ifelse(fet_us_micro_tri3==1,2,ifelse(fet_us_micro_tri1==1,0,NA)))]
data[,fet_micro_abn_tri:=ifelse(fet_us_bin_tri2==1|fet_us_bin_tri2==1|fet_us_cns_tri2==1|fet_us_msk_tri2==1|fet_us_cardio_tri2==1|
                                fet_us_gastro_tri2==1|fet_us_orofac_tri2==1|fet_us_eyeear_tri2==1|fet_us_genur_tri2==1,2,
                                ifelse(fet_us_bin_tri3==1| fet_us_cns_tri3==1|fet_us_msk_tri3==1|fet_us_cardio_tri3==1|
                                       fet_us_gastro_tri3==1|fet_us_orofac_tri3==1|fet_us_eyeear_tri3==1|fet_us_genur_tri3==1,3,NA))]
data[,fet_micro_diag_tri:=ifelse(is.na(fet_micro_diag_tri),fet_micro_diag_trin,fet_micro_diag_tri)]
data[, microcephaly_ga:=ifelse(!is.na(fet_micro_diag_ga),fet_micro_diag_ga,
                               ifelse(fet_micro_diag_tri==0,13,
                                      ifelse(fet_micro_diag_tri==1,27,
                                             ifelse(fet_micro_diag_tri==2,40,NA))))] 

# Check categorical variables
cat_var <- add_info[Type_var == "Categorical"& Var_ga=="yes",]$who_name

 # correct tri assignation until 13 weeks trim1(0),27 trim2(1),trim3(2),postbirth()
  data[,zikv_assay_tri_1:=ifelse(zikv_assay_ga_1<=13,0,ifelse(zikv_assay_ga_1<=27,1,ifelse(zikv_assay_ga_1<=46,2,zikv_assay_tri_1)))]
  data[,zikv_pcr_tri_1:=ifelse(zikv_pcr_ga_1<=13,0,ifelse(zikv_pcr_ga_1<=27,1,ifelse(zikv_pcr_ga_1<=46,2,zikv_pcr_tri_1)))]
  data[,zikv_elisa_tri_1:=ifelse(zikv_elisa_ga_1<=13,0,ifelse(zikv_elisa_ga_1<=27,1,ifelse(zikv_elisa_ga_1<=46,2, zikv_elisa_tri_1)))]
  
  data[,zikv_tri:=ifelse(zikv_ga<=13,0,ifelse(zikv_ga<=27,1,ifelse(zikv_ga<=46,2, zikv_tri)))]
  data[,fet_zikv_tri:=ifelse(fet_zikv_ga<=13,0,ifelse(fet_zikv_ga<=27,1,ifelse(fet_zikv_ga<=46,2, fet_zikv_tri)))]
  
  table(fet=data$fet_zikv_tri,mom=data$zikv_tri,useNA="always")
  data[,zikv_tri_n:=ifelse(is.na(zikv_tri),fet_zikv_tri,zikv_tri)] # aprox of zika diagnosis from fet zika diagnosis
  
  data[,zikv_test_tri_n:=apply(data[,c("zikv_assay_tri_1","zikv_elisa_tri_1" ,"zikv_elisa_tri_1")], 1, min, na.rm = TRUE)] # variable first date of test
  data[,zikv_test_tri_n:=ifelse(zikv_test_tri_n==Inf,NA,zikv_test_tri_n)]
  
  table(zikv=data$zikv_tri_n,assay=data$zikv_test_tri_n,useNA="always") # can we approx zikv_diagnosis trimestre with time of any zikv test 
  data[,zikv_tri_n2:=ifelse(is.na(zikv_tri_n),zikv_test_tri_n,zikv_tri_n)]
 
   data[,zikv_ga_n:=ifelse(!is.na(zikv_ga),zikv_ga,ifelse(zikv_tri_n2==0,13,ifelse(zikv_tri_n2==1,27,ifelse(zikv_tri_n2==2,42,NA))))] # we aprox NA to the max boundary of the trimester



# Create additional variables ----
# 1. Fet_death (fetus death variable) and fet_death_ga (time of fetus death) -----

# 1.1. Add miscarriage, loss, loss_etiology, birth in one variable fet_death and fet_death_ga

data[,maxbirth_ga:=fcase(ch_term==1,42,ch_term==2,28,ch_term==3,21,ch_term==4,33,ch_term==5,36,ch_term==6,44,default=NA)] #max ga according to ch_term
data[,birth_ga:=ifelse(is.na(birth_ga),maxbirth_ga,birth_ga)]

data[,birth:=ifelse(!is.na(birth_ga)|!is.na(ch_term),1,NA)] #birth indicator
data[,fet_death:=ifelse(inducedabort==1,1,ifelse(birth==1,0,NA))]
data[,fet_death_ga:=ifelse(!is.na(endga),endga,ifelse(!is.na(loss_ga),loss_ga,ifelse(!is.na(miscarriage_ga),miscarriage_ga,ifelse(!is.na(inducedabort_ga),inducedabort_ga,NA))))]
data[,end_ga := ifelse(!is.na(endga),endga,ifelse(!is.na(birth_ga),birth_ga,fet_death_ga))] #clean end_ga 
data[,end_ga :=ifelse(is.na(end_ga),maxbirth_ga,end_ga)] 



#9  Missing data Plot

table_inc <- add_info[Var_ga=="yes",c("who_name","order_ga")]
table_inc<-table_inc[order(order_ga)]
colnames(table_inc)<-c("variable","order_ga")
var_incl<-table_inc$variable
dataf<-data[,..var_incl]
totalval<-as.data.table(table(dataf$studyname))
colnames(totalval)<-c("studyname","N")
dmatrix<-dataf[, lapply(.SD, function(x) sum(is.na(x))/.N), studyname] #matrix of % of missingness
dmatrix2<-as.data.table(melt(dmatrix,id.vars="studyname"))
dmatrix2<-merge(dmatrix2,totalval,by="studyname")

dmatrix2[,name:=paste0(studyname,"\nN=",N)]
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
