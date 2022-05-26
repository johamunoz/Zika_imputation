
rm(list=ls())

# Data manipulation packages
library(rstudioapi) 
library(data.table)
library(dplyr)
library(here) 

# Graphic packages
library(purrr)
library(tidyr)
library(ggplot2)
library(kableExtra)

#0. Load information ----

load(here('3_Output_data','pldata.RData'))
infoexp <- as.data.table(readxl::read_xlsx(here('1_Input_data','Infoexp.xlsx'),sheet="Table")) #CSV file #Table were are specified the included variables according Expert opinion, also the includes the order in which variables are imputed 
var_inc<-infoexp[Inclusion==1,Variable] #Variables to work with
data<-as.data.table(pldata)[,colnames(pldata)%in%var_inc,with = FALSE]
data<-data[,-c("studyname","mid_original","childid_original")]

var_con  <- infoexp[Vtype=="C",Variable] #Variables to work with
datacont <- as.data.table(pldata)[,colnames(pldata)%in%c(var_con,"studycode"),with = FALSE]
long <- melt(setDT(datacont), id.vars = c("studycode"), variable.name = "variable")
long[,group:=paste0(variable,"&",studycode)]
long[,total:=paste0(variable,"&Total")]

myFuncon <- function(x) {
  c(n= sum(!is.na(x)),
    na= sum(is.na(x)),
    mean = mean(x,na.rm=T),
    std = sd(x,na.rm=T),
    min = min(x,na.rm=T), 
    q2.5 = as.numeric(quantile(x, .025,na.rm=T)),
    q25  = as.numeric(quantile(x, .25,na.rm=T)),
    q50  = as.numeric(quantile(x, .5,na.rm=T)),
    q75  = as.numeric(quantile(x, .5,na.rm=T)),
    q97.5 = as.numeric(quantile(x, .975,na.rm=T)),
    max = max(x,na.rm=T))
}

list_study<-tapply(long$value, long$group, myFuncon)
list_total<-tapply(long$value, long$total, myFuncon)
table_con<-do.call(rbind,c(list_study,list_total))
names<-rownames(table_con)
table_con<-data.table(table_con)
table_con[,name:=names]
table_con[table_con=="NaN"] <-  NA
table_con[table_con=="Inf"] <-  NA
table_con[table_con=="-Inf"] <-  NA
table_con[, c("variable","studycode") := tstrsplit(name, "&")]
table_con<-table_con[order(variable,studycode)]
table_con[,name:=NULL]
setcolorder(table_con, c("variable", "studycode", "n","na","mean","std","min", "q2.5","q25","q50","q75","q97.5","max"))
write.csv(table_con,here('6_Tables_graphs', 'Descriptive_continuous'))

#kable format
table_con %>%
  kbl(caption = "Continuos variables",
      digits = c(0,0,0,0,3,3,3,3,3,3,3,3),
      col.names = c("Variable","Study code", "Valid","Missing"," mean","std","min","2.5%","25%","50%","75%","97.5","max")) %>%
  kable_paper(full_width = F) %>%
  column_spec(1, bold = T) %>%
  collapse_rows(columns = 1:2, valign = "top")%>%
  kable_styling(fixed_thead = T,bootstrap_options = c("striped", "hover", "condensed","responsive"))




# Categorical & Binary variables

var_cat <- infoexp[Vtype%in%c("CAT","B"),Variable] #Variables to work with
datacat <- as.data.table(pldata)[,colnames(pldata)%in%c(var_cat,"studycode"),with = FALSE]

longcat <- melt(setDT(datacat), id.vars = c("studycode"), variable.name = "variable")
longcat_study<-longcat[, .N, by=.(studycode,variable,value)]
longcat_total<-longcat_study[, .(N=sum(N,na.rm=T)), by=.(variable,value)]
longcat_total[,studycode:="Total"]
longcatf<-rbind(longcat_study,longcat_total)
longcat_N<-longcatf[, .(total=sum(N,na.rm=T)), by=.(studycode,variable)]
longcatf<-as.data.table(merge(longcatf,longcat_N,by=c("studycode","variable")))
longcatf[,Proportion:=N/total*100]
longcatf[,value:=as.numeric(value)]
longcatf<-longcatf[order(variable,studycode,value)]
longcatf[,total:=NULL]
write.csv(longcatf,here('6_Tables_graphs', 'Descriptive_categorical'))

#Kable format
longcatf1<-dcast(longcatf, variable+value ~ studycode, value.var = c("N", "Proportion"),fun.aggregate = sum)
vec<-apply(expand.grid(c("N","Proportion"), c(unique(longcat$studycode),"Total")), 1, paste, collapse="_")
setcolorder(longcatf1, c("variable","value",vec))

hv=c(1,1,rep(2,11))
names(hv)=c(" "," ",unique(longcat$studycode),"Total")

longcatf1%>%
  kbl(caption="Categorical",
      digits=c(0,0,rep(c(0,1),11)),
      col.names = c("Variable","Level", rep(c("N","%"),11)))%>%
  add_header_above(hv)%>%
  collapse_rows(columns = 1, valign = "top")%>%
  kable_styling(fixed_thead = T,full_width = F,bootstrap_options = c("striped", "hover", "condensed","responsive"))


longcatf %>%
  kbl(caption = "Categorial variables",
      digits = c(0,0,0,0,3),
      col.names = c("Variable","Study code", "Level","N","%")) %>%
  kable_paper(full_width = F) %>%
  column_spec(1, bold = T) %>%
  collapse_rows(columns = 1:2, valign = "top")%>%
  kable_styling(fixed_thead = T,bootstrap_options = c("striped", "hover", "condensed","responsive"))



#Plot continuos

datacont %>%
  select(-pre_pregweight,-zikv_pcr_vl_1,-studycode)%>%
  gather() %>% 
  group_by(key) %>% 
  mutate(mean = mean(value,na.rm=T)) %>% # calculate mean for plotting as well
  ungroup()%>%
  ggplot(aes(value)) +
  geom_density(aes(y = ..count.. * bin), # multiply count by bins
               fill = "blue", alpha = .3, col = NA) + 
  geom_histogram(binwidth = bin, alpha = .5) + # use the same bins here
  geom_vline(aes(xintercept = mean), col = "black",linetype="dashed") + 
  theme_minimal() + 
  labs(y = "Count",x="Value") +
  ggtitle("Histogram of continuous variables")+
  facet_wrap(~ key, scales = "free") 
  
