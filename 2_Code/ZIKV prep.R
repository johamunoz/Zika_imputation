
#Starts directly after loading imputed dataset

data <- complete(merged_imp, "long")
rm(merged_imp)
m<-max(data$.imp)
studynames<-c("014-BRA","001-BRA","002-BRA","010-BRA","007-COL",
              "003-GUF","005-ESP","004-ESP","012-TTO","008-USA")
data$studyname_fac<-factor(data$studyname)
levels(data$studyname_fac)<-studynames

#Create dichotomous outcome variables to calculate incidence
#Microcephaly_bin
data$microcephaly_bin<-data$microcephaly
data$microcephaly_bin[data$microcephaly_bin==1]<-1
data$microcephaly_bin[data$microcephaly_bin==2]<-1
data$microcephaly_bin[data$microcephaly_bin==3]<-0 #Macrocephaly
data$microcephaly_bin <- droplevels(data$microcephaly_bin)
#miscarriage (<20 weeks gestation)
data$miscarriage<-0
data$miscarriage[data$bdeath==1 & data$bdeath_ga<20]<-1
data$miscarriage<-as.factor(data$miscarriage)
#Fetal loss (>=20 weeks gestation)
data$loss<-0
data$loss[data$bdeath==1 & data$bdeath_ga>=20]<-1
data$loss<-as.factor(data$loss)
#Early fetal death (20-27 weeks gestation)
data$efdeath<-0
data$efdeath[data$bdeath==1 & data$bdeath_ga>=20 & data$bdeath_ga<28]<-1
data$efdeath<-as.factor(data$efdeath)
#Late fetal death (after 28 weeks gestation)
data$lfdeath<-0
data$lfdeath[data$bdeath==1 & data$bdeath_ga>=28]<-1
data$lfdeath<-as.factor(data$lfdeath)
#Late fetal death (after 28 weeks gestation) with microcephaly
#data$lfdeath_micro<-0
#data$lfdeath_micro[data$lfdeath==1 & data$bdeath_ga>=28]<-1
#data$lfdeath<-as.factor(data$lfdeath)

data.zika<-data[data$zikv_preg==1,]
data.nozika<-data[data$zikv_preg==0,]