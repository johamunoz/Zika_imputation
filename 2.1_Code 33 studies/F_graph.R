
#Functions necessary for objective 1

library(tidyr)
library(purrr)
library(magrittr)

#######Functions#######
logit <- function(x){x<-ifelse(x==0,1e-10,x);x<-ifelse(x==1,1-1e-10,x);log(x/(1-x))} # better remove NA
logit_inv <- function(x){exp(x)/(1+exp(x))}


#Function to calculate absolute risk of an outcome per study and per imputed dataset
f_sumarize_data <- function(data){
    alpha <- 0.05  
    sum_data<-data%>%summarise(obs_n=n(),obs_sum=sum(outcome,na.rm=TRUE),na_obs=sum(is.na(outcome),na.rm=TRUE))
    ptest = prop.test(sum_data$obs_sum, sum_data$obs_n, conf.level = alpha)
    sum_data%<>%
      mutate(  observed = ptest$estimate,
               lower = ptest$conf.int[[1]],
               upper = ptest$conf.int[[2]])%>%
      mutate( observed = ifelse(observed==0,0.001,observed),
              lower = ifelse(lower==0,0.001,lower),
              upper = ifelse(upper==0,0.001,upper))%>%
      # Logit and se for estimate and observed          
      mutate(logit_obs = logit(observed),
             logit_obs_se = sqrt((logit(upper)-logit(lower))/(2*qnorm(1-alpha/2)*obs_n)))
    return(sum_data)
  }
  

f_studyimp<-function(data, outcome_name) {
  sum_data<- data%>%
             dplyr::select(outcome=!!outcome_name,studyname,'.imp')%>%
             nest(data= - c(studyname,'.imp'))%>%
             mutate(obs_sum = map(data,~f_sumarize_data(data=.x)))%>%
             select(-c(data))%>%
             unnest(col= c(obs_sum))
  return(sum_data)
}


#Function to pool the absolute risks per study per imputed dataset using rubins rules

f_pool <- function(data,est,se){
  alpha <- 0.05
  pdata <- data%>% dplyr::select(est = !!est, se = !!se)
  n <- nrow(pdata)  #number observations
  mu <- mean(pdata$est,na.rm=T) # pool estimate
  w_var <- mean(pdata$se^2,na.rm=T) # within variance
  b_var <- var(pdata$est,na.rm=T) # between variance
  t_var <- w_var + b_var+ b_var/n # total variance
  t_se <- sqrt(t_var) # total standard error
  r <- (b_var + (b_var / n))/ w_var # relative increase variance due to missing values
  v <- (n - 1) * (1 + r^-1)^2 # degrees of freedom
  t <- qt(1-alpha/2, v) #t distribution criteria
  lower <- mu - t_se*t
  upper <- mu + t_se*t
  return(data.frame(cbind(mu=mu,lower=lower,upper=upper)))
}

f_data_pool <-function(data,inc.outcome) {
  data_pool<-data%>%
             filter(.imp!=0)%>%  
             nest(data=-studyname)%>%  
             mutate(pool_sum=map(data,function(x){
               imp_n = nrow(x)
               p_obs=f_pool(data=x,est="logit_obs",se="logit_obs_se")
               observed = logit_inv(p_obs$mu)
               lower = logit_inv(p_obs$lower)
               upper = logit_inv(p_obs$upper)
               return(cbind(p_obs=p_obs,observed=observed,lower=lower,upper=upper))}))%>%
             select(-c(data))%>%
            unnest(col= c(pool_sum))
  return(data_pool)}


forest_plot_study<-function(data,outcome_name,plottitle){
  
  #Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
  inc.outcome<-f_studyimp(data=data,outcome_name)
  
  
  #For raw dataset
  abs.outcome.raw<-inc.outcome%>%filter(.imp==0)%>%select(c("studyname","observed","lower","upper","obs_n","na_obs"))
  abs.outcome.raw$source<-"Raw"
  
  #For imputed datasets
  #Pool the absolute risks over the imputed datasets, resulting in absolute risks per study
  abs.outcome.imp<-f_data_pool(data=inc.outcome)%>%select(c("studyname","observed","lower","upper"))
  abs.outcome.imp$source<-"Imputation"
  abs.outcome.imp%<>%left_join(abs.outcome.raw%>%select(c("studyname","obs_n","na_obs")),by="studyname")
  
  
  #Collapse both datasets
  
  outcome <- dplyr::bind_rows(abs.outcome.imp, abs.outcome.raw)%>%
    mutate(mean = observed,
           se = (upper-lower)/1.96,
           cint=ifelse(is.na(se), "",
                       paste0(sprintf("%.2f(%.2f,%.2f)",
                                      observed, lower, upper))),
           studyname=paste0(studyname),
           Nobs=paste0("N=",obs_n,",NA=",na_obs))
  
  
  
  
  # Create plot
  dataplot<-setDT(outcome)[,c("studyname","source","cint","mean","lower","upper","Nobs")]
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
    xlab("Study name")+ ylab("Absolute risk")+
    ggtitle(plottitle)+
    coord_flip() +
    scale_x_discrete(breaks=levels(factor(dataplot$allcint)),labels = dataplot$cint)+ 
    facet_grid(studyname + Nobs~ ., switch = "y",scales="free")+
    theme(strip.placement = "outside")+
    theme(strip.text.y.left = element_text(angle = 0),axis.text.y = element_text(size = 6))
  return(plot)}
               
