###Aim: Additional functions for comparing estimations between raw and imputed datasets.


# Load packages ---
# Data manipulation package
library(data.table) 
#setwd("/Users/jdamen/Documents/GitHub/Zika_imputation")
library(here)  # define folder paths
library(metafor)

# Graphic packages
library(ggplot2)

# Load dataset and dependencies ----
load(file =here('3_Output_data','rawfinaldata33.RData')) # fdata=Data_preimputation
study_info <- as.data.table(readxl::read_xlsx(here('1_Input_data','MasterCodebook_October.xlsx'),sheet="StudyID")) #CSV file with the
add_info <- as.data.table(readxl::read_xlsx(here('1_Input_data','MasterCodebook_October.xlsx'),sheet="237 key")) #CSV file with the
data_imp <- as.data.table(read.csv(here('1_Input_data','20221027 zikv_imputed.csv'))) #CSV file with the imputed datasets
data_imp <- merge(data_imp,study_info[,c("studyimp","studyname")],by="studyimp",all.x=TRUE)
data_sys <-as.data.table(readxl::read_xlsx(here('1_Input_data','Table 1 outcomes.xlsx'),sheet="Systematic missings")) #CSV file with the

source(here('2.1_Code 33 studies','FunctionsObjective1mod.R'))


# Give format to variables----

format_data<-function(data,source){
  need_col<-add_info[Figures=="yes"]$who_name
  need_col <- need_col[need_col %in% colnames(data)]
  
  # Drop variables we do not need
  datan<-setDT(data)[,..need_col]
  
  # Create additional variables 
  if(!".imp"%in%colnames(datan)){datan[,'.imp':= 0]}
  if(!"bdeath"%in%colnames(datan)){datan[, bdeath := 1-birth]}
  if(!"miscarriage"%in%colnames(datan)){datan[, miscarriage := ifelse(bdeath==1&end_ga<20,1,0)]} # Miscarriage
  if(!"loss"%in%colnames(datan)){datan[, loss := ifelse(bdeath==1&end_ga>=20,1,0)]}  #Fetal loss (>=20 weeks gestation)
  datan[, efdeath := ifelse(bdeath==1&end_ga>=20&end_ga<28,1,0)]   #Early fetal death (20-27 weeks gestation)
  datan[, lfdeath := ifelse(bdeath==1&end_ga>=28,1,0)]    #Late fetal death (after 28 weeks gestation)
  
  # Convert categorical & binary variables to factors
  need_colfac<-add_info[Figures=="yes"&Type_var%in%c("Categorical","Binary")]$who_name
  datan[, (need_colfac) := lapply(.SD, factor), .SDcols = need_colfac]
  datan[, source:=source]
}

data_imp_f <- format_data(data=data_imp,source="imputation")
data_raw_f <- format_data(data=data_raw,source="raw")

# Separate dataset according to zikv_preg
data_all<-as.data.table(rbind(data_imp_f,data_raw_f))
data.zika<-data_all[zikv_preg==1,]
data.nozika<-data_all[zikv_preg==0,]

# Remove values on systematically missing data
data_sys <-as.data.table(melt(setDT(data_sys), id.vars = c("variable"), variable.name = "studyname",value.name="systematic"))
data_sys<-data_sys[systematic==1]

data_alls<-copy(data_all)
for (i in 1:nrow(data_sys)){
data_alls[studyname==data_sys[i,]$studyname,(data_sys[i,]$variable):=NA]
}
data.zika<-data_alls[zikv_preg==1,]
data.nozika<-data_alls[zikv_preg==0,]


forest_plot_study<-function(data,outcome_name,syst,plottitle){
    # Remove systematical missing studies
    if(syst==TRUE){
      systvector<-as.vector(data_sys[variable==outcome_name]$studyname)
      data<-data[!studyname%in%systvector]
    }
  
    #Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
    inc.outcome<-f.abs.perstudy(data=data,outcome_name)
    
    #For imputed datasets
    #Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
    abs.outcome.imp<-f.abs.poolrubin(data=inc.outcome[.imp!=0])
    abs.outcome.imp[,source:="Imputation"]
    #Pool the absolute risks over the studies, resulting in one pooled summary absolute risk over all studies
    pool.outcome.imp<-f.abs.2s.ma(abs.outcome.imp)
    pool.outcome.imp$source<-"Imputation"
    
    #For raw dataset
    abs.outcome.raw<-f.abs.poolrubin(data=inc.outcome[.imp==0])
    abs.outcome.raw[,source:="Raw"]
    pool.outcome.raw<-f.abs.2s.ma(abs.outcome.raw)
    pool.outcome.raw$source<-"Raw"
    
    #Collapse both datasets
    abs.outcome<-as.data.table(rbind(abs.outcome.imp,abs.outcome.raw))
    pool.outcome<-as.data.table(rbind(pool.outcome.imp,pool.outcome.raw))
    setnames(pool.outcome, "abs.risk", "incidence")
    pool.outcome[,studyname:="TOTAL"]
    
    outcome<-rbind(abs.outcome,pool.outcome,fill=TRUE)
    outcome[,mean:=incidence]
    outcome[,lower:=ci.lb]
    outcome[,upper:=ci.ub]
    outcome[,se:=(ci.ub-ci.lb)/1.96]
    outcome[,cint :=ifelse(is.na(se), "",
                               paste0(sprintf("%.2f(%.2f,%.2f)",
                                       incidence, ci.lb, ci.ub)))]
    
    # Create plot
    dataplot<-outcome[,c("studyname","source","cint","mean","lower","upper")]
    dataplot[,source:=factor(source)]
    dotcols = c("#a6d8f0","#f9b282")
    barcols = c("#008fd5","#de6b35")

    dataplot[,allcint:=paste(studyname,source,cint,sep="-")] 
    dataplot<-dataplot[order(studyname,source),]

    plot<-ggplot(dataplot, aes(x=allcint, y=mean, ymin=lower, ymax=upper,col=source,fill=source)) + 
      geom_linerange(size=2,position=position_dodge(width = 0.5)) +
      geom_point(size=1, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
      facet_grid(studyname~ ., switch = "y",scales="free")+
      scale_fill_manual(values=barcols)+
      scale_color_manual(values=dotcols)+
      scale_y_discrete(limits=c("Imputation","Raw"))+
      xlab("Study name")+ ylab("Absolute risk")+
      ggtitle(plottitle)+
      coord_flip() +
      scale_x_discrete(labels = dataplot$cint)+ 
      facet_grid(studyname~ ., switch = "y",scales="free")+
      theme(strip.placement = "outside")+
      theme(strip.text.y.left = element_text(angle = 0),axis.text.y = element_text(size = 6))
return(plot)}


forest_plot_study(data=data.zika,outcome_name="czs",syst=TRUE,plottitle = "CZS remove systematical")
forest_plot_study(data=data.zika,outcome_name="microcephaly_bin_birth",syst=TRUE,plottitle = "Microcephaly")





