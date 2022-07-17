
rm(list=ls()) # clean environment

# Data manipulation packages
library(here) 
library(data.table)
library(dplyr)
require(ggplot2)

data <- as.data.table(read.csv(here('1_Input_data','zikv_033_datasets.csv'), stringsAsFactors=FALSE, fileEncoding="latin1"))
infoexp <- as.data.table(readxl::read_xlsx(here('1_Input_data','MasterCodebook_Final_June2022 ALL (Repaired).xlsx'),sheet="237 key")) #CSV file with the
var_inc<-as.matrix(infoexp[Essential=="Yes",'WHO variable name'])[,1] #Variables to work with
data<-as.data.table(data[,..var_inc]) #Filter dataset


#1. Initial clean up the dataset -----
#1.0. Check for duplicates ----
data[,comb:=paste0(studyname,mid,fetid,childid,sep="_")] #create a variable that combines study and patient id
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

#1.1. Set to NA : 666,777,888,999 values----
data[data==""] <-NA
data[data==555] <-  NA
data[data==666] <-  NA
data[data==777] <-  NA
data[data==888] <-  NA
data[data==999] <-  NA
data[data==9999] <-  NA
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
save(obs_to_check,file=here('3_Output_data','obs_to_check.RData'))

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


