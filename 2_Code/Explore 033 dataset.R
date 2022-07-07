
library(ggplot2)
library(naniar)

setwd("/Users/jdamen/Documents/Julius/ZIKV analyses/2. Data")
zikv<-read.csv("zikv_033_datasets.csv",header=T)

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

#write.csv(missings,file="/Users/jdamen/Documents/Julius/ZIKV analyses/4. Resultaten/20220706 Percentage missings.csv")

summary(zikv$weight)
