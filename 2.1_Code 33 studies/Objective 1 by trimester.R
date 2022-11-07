
library(dplyr)
library("xlsx")
library(mitml)
library(mice)
library(micemd)
library(metafor)
library(metamisc)
library(EpiStats)
library("lme4")
library("logistf")
library(data.table) 
library(ggplot2)


######################################################################################
#################################Load and prepare data################################
######################################################################################

setwd("/Users/jdamen/Documents/Julius/ZIKV analyses/2. Data")
data<-read.csv("20221027 zikv_imputed.csv",header=T)
setwd("/Users/jdamen/Documents/GitHub/Zika_imputation/2.1_Code 33 studies")
source("Functions Objective 1.R")
source("ZIKV prep v33.R")

#Systematic missings
data_sys <- read.xlsx2("/Users/jdamen/Documents/GitHub/Zika_imputation/1_Input_data/Table 1 outcomes.xlsx",sheetName="Systematic missings")
colnames(data_sys)<-c("variable","001-BRA","002-BRA","003-GUF","004-ESP","005-ESP","006-COL","007-COL",
                      "008-USA","009-GRD","010-BRA","011-BRA","012-TTO","013-BRA","014-BRA",
                      "015-BRA","016-HND","017-USA","018-COL","019-BRA","020-BRA","021-PRI",
                      "022-BRA","023-BRA","024-GTM","025-BRA","026-KEN","027-BRA")
data_sys<-melt(setDT(data_sys), id.vars = c("variable"), variable.name = "studyname",value.name="systematic")
data_sys<-data_sys[systematic==1]
data<-as.data.table(data)
for (i in 1:nrow(data_sys)){
  data[studyname_fac==data_sys[i,]$studyname,(data_sys[i,]$variable):=NA]
}
data<-as.data.frame(data)

data.zika.all<-data[data$zikv_preg==1,]
data.nozika.all<-data[data$zikv_preg==0,]


########################Analyses#############
#Change directory to save plots
setwd("/Users/jdamen/Documents/Julius/ZIKV analyses/4. Resultaten")

#Data frame to store results
all.results<-data.frame(Outcome=character(),
                        File=numeric(), 
                        User=numeric(), 
                        stringsAsFactors=FALSE,
                        Trimester=numeric())


##################################################################################
###################################Microcephaly###################################
##################################################################################

data.zika1<-data.zika.all[!is.na(data.zika.all$microcephaly_bin_birth) & data.zika.all$zikv_tri==1,]
data.zika2<-data.zika.all[!is.na(data.zika.all$microcephaly_bin_birth) & data.zika.all$zikv_tri==2,]
data.zika3<-data.zika.all[!is.na(data.zika.all$microcephaly_bin_birth) & data.zika.all$zikv_tri==3,]


#Zika-positive women
#Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
inc.outcome<-f.abs.perstudy(data.zika1,"microcephaly_bin_birth")
abs.outcome1<-f.abs.poolrubin(data.zika1,inc.outcome)
abs.outcome1$trimester<-"First trimester"
inc.outcome<-f.abs.perstudy(data.zika2,"microcephaly_bin_birth")
abs.outcome2<-f.abs.poolrubin(data.zika2,inc.outcome)
abs.outcome2$trimester<-"Second trimester"
inc.outcome<-f.abs.perstudy(data.zika3,"microcephaly_bin_birth")
abs.outcome3<-f.abs.poolrubin(data.zika3,inc.outcome)
abs.outcome3$trimester<-"Third trimester"

result<-rbind(abs.outcome1,abs.outcome2,abs.outcome3)
result$cint<-paste0(sprintf("%.2f",result$incidence),"(",sprintf("%.2f",result$ci.lb),",",sprintf("%.2f",result$ci.ub),")")
result[is.nan(result$logit.var),]$cint<-NA
result$trimester<-as.factor(result$trimester)

result<-result[order(result$studyname,result$trimester),]
# Create plot
#dataplot<-outcome[,c("studyname","source","cint","mean","lower","upper")]
dotcols = c("#a6d8f0","#f9b282","green")
barcols = c("#008fd5","#de6b35","lightgreen")

ggplot(result, aes(x=cint, y=incidence, ymin=ci.lb, ymax=ci.ub,col=trimester,fill=trimester)) + 
  geom_linerange(size=2,position=position_dodge(width = 0.5)) +
  geom_point(size=1, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  scale_fill_manual(values=barcols)+
  scale_color_manual(values=dotcols)+
  xlab("Study name")+ ylab("Absolute risk")+
  ggtitle("test")+
  coord_flip() +
  facet_grid(studyname ~ ., switch = "y",scales="free")+
  theme(strip.placement = "outside")+
  theme(strip.text.y.left = element_text(angle = 0),axis.text.y = element_text(size = 6))


results<-pool.outcome

all.results<-rbind(all.results,results)

#Forest plot
#png(file="20221031 Microcephaly zika positive.png",width=750,height=750,res=100)
metafor::forest(abs.outcome$incidence, ci.lb=abs.outcome$ci.lb, ci.ub=abs.outcome$ci.ub, 
                refline = 0, slab = abs.outcome$studyname,
                xlab = "Absolute risk (%)", pch = 19, psize=1,
                ylim=(c(-1,sum(!is.na(abs.outcome$incidence))+3)), cex=1, steps=7, main="Microcephaly, zika-positive women",
                digits=2L,xlim=(c(-15,85)), alim=(c(0,60)))