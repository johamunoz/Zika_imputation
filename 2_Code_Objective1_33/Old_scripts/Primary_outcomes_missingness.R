
#Passive imputation data provided 09-09-21

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


data<-data[,c("microcephaly_bin","czs","miscarriage","loss","studyname")]

data[data==""] <-NA
data[data==666] <-  NA
data[data==777] <-  NA
data[data==888] <-  NA
data[data==999] <-  NA
data[data==9999] <-  NA


tmicro<-as.data.table(table(data$studyname,data$microcephaly_bin,useNA = "always"))
colnames(tmicro)<-c("Studyname","Response","Number")
tmicro$variable<-"microscephaly_bin"
tczs<-as.data.table(table(data$studyname,data$czs,useNA = "always"))
colnames(tczs)<-c("Studyname","Response","Number")
tczs$variable<-"czs"
tloss<-as.data.table(table(data$studyname,data$loss,useNA = "always"))
colnames(tloss)<-c("Studyname","Response","Number")
tloss$variable<-"loss"
tmisc<-as.data.table(table(data$studyname,data$miscarriage,useNA = "always"))
colnames(tmisc)<-c("Studyname","Response","Number")
tmisc$variable<-"miscarriage"

tprimary<-rbind(tmicro,tczs,tloss,tmisc)
tprimary<-tprimary[!is.na(Studyname),]
head(tprimary)

total<-as.data.table(table(data$studyname))
colnames(total)<-c("Studyname","Total")
tprimary<-merge(tprimary,total,by="Studyname")

tprimary[,name:=paste0(Studyname,"\nN=",Total)]

#Microscephaly plot
ggplot(tprimary[variable=="microscephaly_bin"], aes(fill=Response, y=Number, x=name)) + 
  geom_bar(position="fill", stat="identity")+
  theme(axis.text.x = element_text(angle = 45,vjust=0.5,size=8),
        axis.text.y = element_text(size=8))+xlab("Study name")+ylab("Variable")+
  ylab("Microcephaly binary")

unique(tprimary$variable)
#CZS plot
ggplot(tprimary[variable=="czs"], aes(fill=Response, y=Number, x=name)) + 
  geom_bar(position="fill", stat="identity")+
  theme(axis.text.x = element_text(angle = 45,vjust=0.5,size=8),
        axis.text.y = element_text(size=8))+xlab("Study name")+ylab("Variable")+
  ylab("CZS")

#loss plot
ggplot(tprimary[variable=="loss"], aes(fill=Response, y=Number, x=name)) + 
  geom_bar(position="fill", stat="identity")+
  theme(axis.text.x = element_text(angle = 45,vjust=0.5,size=8),
        axis.text.y = element_text(size=8))+xlab("Study name")+ylab("Variable")+
  ylab("Loss")

#miscarriage plot
ggplot(tprimary[variable=="miscarriage"], aes(fill=Response, y=Number, x=name)) + 
  geom_bar(position="fill", stat="identity")+
  theme(axis.text.x = element_text(angle = 45,vjust=0.5,size=8),
        axis.text.y = element_text(size=8))+xlab("Study name")+ylab("Variable")+
  ylab("Miscarriage")

