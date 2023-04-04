
#Functions necessary for objective 1

library(tidyr)
library(purrr)
library(magrittr)
library(ggplot2)
library(ggh4x)
library(metafor)

#######Functions#######

trans <-function(e,n,type){
  if (type=="logit") {
    if(any(e == 0) | any(e == n)){
      e <- e + 0.5
      n <- n + 1
    }
    p <- e/n
    y <- log(p/(1 - p))
    v <- 1/e + 1/(n - e)
  }else{ # type arcsine
    p <- e/n
    y <- asin(sqrt(p))
    v <- 1/(4 * n)
  }
  se <- sqrt(v/n)
  return(data.frame(y = y,se = se))
} 


back_trans <- function(x,type){
  if (type=="logit") {
    y <- exp(x)/(1 + exp(x))
  }else{
    y <- (sin(x))^2
    }
  return(y)
} 



# Function to calculate the observed absolute risk, its se and CI (observed, observed_se,lower,upper) and the transformed absolute risk and its se (tr_obs, tr_se)
f_sumarize_data <- function(data,type){
    alpha <- 0.05  
    sum_data<-data%>%summarise(obs_n = n(),
                               obs_sum = sum(outcome,na.rm=TRUE),
                               na_obs = sum(is.na(outcome),na.rm=TRUE),
                               included =first(Included))
    
    if (sum_data$obs_n == sum_data$na_obs){
          observed <- NA
          lower <-  NA
          upper <- NA
          tr_obs <- NA
          tr_se <- NA
    }else{
        
        ptest <- prop.test(sum_data$obs_sum, sum_data$obs_n, conf.level = alpha, correct = TRUE) # exact CI with Yates correction
        trans_observed <- trans(e=sum_data$obs_sum, n=sum_data$obs_n, type = type)
        observed = ptest$estimate
        lower =   ptest$conf.int[[1]]
        upper =   ptest$conf.int[[2]]
        tr_obs = trans_observed$y
        tr_se = trans_observed$se
    }
    sum_data%<>% 
      # in the original scale
      mutate( observed = observed,
              lower = lower,
              upper = upper,
              tr_obs = tr_obs,
              tr_se = tr_se)
    return(sum_data)
  }
  


f_study_imp <- function(data, outcome_name, type) {
  sum_data<- data%>%
             dplyr::select(outcome=!!outcome_name,studyname,Included,'.imp')%>%
             nest(data= - c(studyname,'.imp'))%>%
             mutate(obs_sum = map(data,~f_sumarize_data(data=.x,type=type)))%>%
             select(-c(data))%>%
             unnest(col= c(obs_sum))
  return(sum_data)
}


#Function to pool the absolute risks per study per imputed dataset using rubins rules

f_pool <- function(est, tr_est, tr_se, type){
  
  alpha <- 0.05
  n <- length(tr_est)  #number observations
  mu <- mean(tr_est, na.rm = T) # pool estimate
  w_var <- mean(tr_se^2, na.rm = T) # within variance
  b_var <- var(tr_est, na.rm = T) # between variance
  t_var <- w_var + b_var + b_var/n # total variance
  t_se <- sqrt(t_var) # total standard error
  r <- (b_var + (b_var / n))/ w_var # relative increase variance due to missing values
  v <- (n - 1) * (1 + r^-1)^2 # degrees of freedom
  t <- qt(1-alpha/2, v) #t criteria
  observed <-  mean(est, na.rm = T)
  lower <- back_trans(mu - t_se*t, type = type)
  upper <- back_trans(mu + t_se*t, type = type)
  
  return(data.frame(cbind(observed = observed, lower = lower, upper = upper, tr_obs = mu, tr_se = t_se)))
}


f_data_pool <- function(data,type) {
  data_pool<-data%>%
             filter(.imp>0)%>%  
             nest(data=-studyname)%>%  
             mutate(pool_sum=map(data,function(x){
               imp_n = nrow(x)
               pool_obs=f_pool(est=x$observed,tr_est=x$tr_obs,tr_se=x$tr_se,type=type)}))%>%
    select(-c(data))%>%
            unnest(col= c(pool_sum))
  return(data_pool)}

# Function to calculate the global absolute risk
f_data_total <- function(data,type) {
  
  data_pool <- data%>%
               filter(included==1)
  
  fit.rma <- rma(yi=data_pool$tr_obs,
               sei=data_pool$tr_se,
               method="REML",
               test= "knha")
  
  
  total <- data.frame(  studyname="TOTAL",
                         obs_n = sum(data_pool$obs_n,na.rm = T),
                         obs_sum = sum(data_pool$obs_sum,na.rm = T),
                         na_obs = sum(data_pool$na_obs,na.rm = T),
                         included = sum(data_pool$included,na.rm=T),
                         observed = back_trans(as.numeric(fit.rma$beta),type),
                         lower = back_trans(as.numeric(fit.rma$ci.lb),type),
                         upper = back_trans(as.numeric(fit.rma$ci.ub),type),
                         tr_obs = as.numeric(fit.rma$beta),
                         tr_se = as.numeric(fit.rma$se))
  
  datan <- data%>%rows_insert(total)
  return(datan)}
  

forest_plot_study<-function(data,outcome_name,plottitle,type){
  #Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
  inc.outcome.ind <- f_study_imp(data=data, outcome_name=outcome_name, type= type)
  
  #Estimate global for each dataset (raw and imputed using only the included datasets)
  inc.outcome <- inc.outcome.ind%>%
                 nest(data=-.imp)%>% 
                 mutate( tdata = map(data,~f_data_total(data=.x,type=type)))%>%
                 select(-c(data))%>%
                 unnest(col= c(tdata))%>%
                 mutate(pmiss=na_obs/obs_n*100)  
  
  # Initial num observations
  abs.n <- inc.outcome%>%filter(.imp==0)%>%
           mutate(N=obs_n)%>%
            select("studyname","N","included")
  
  # For original dataset
  abs.outcome.ori <- inc.outcome%>%filter(.imp==-1)%>%
    mutate( source = "Original")%>%
    select(c("studyname","observed","lower","upper","tr_obs","tr_se","source","pmiss"))
  
  # For raw dataset
  abs.outcome.det <- inc.outcome%>%filter(.imp==0)%>%
    mutate( source = "Deterministic")%>%
    select(c("studyname","observed","lower","upper","tr_obs","tr_se","source","pmiss"))
  
  
  # For imputed datases, pool information with Rubins rules
   
   abs.outcome.imp <-f_data_pool(data=inc.outcome,type=type)
   abs.outcome.imp$source <-"Imputation"
   abs.outcome.imp$pmiss <-0.00 
  
  # Merge both datasets
  
  abs.outcome <- dplyr::bind_rows(abs.outcome.imp, abs.outcome.det, abs.outcome.ori)%>%
                 left_join(abs.n, by=c('studyname'))%>%
                  mutate(included = ifelse(is.na(included),"*",""),
                         mean = observed*100,
                         clb = lower*100,
                         cub = upper*100,
                         pmiss = sprintf("%.2f",pmiss))%>%
                  mutate(cint = ifelse(is.na(mean),"   ",paste0(sprintf("%.2f(%.2f,%.2f),",mean, clb, cub))))%>%
                  mutate(allcint = paste0(studyname,included," \n N=",N),
                         cint = paste0(cint," %mis=",pmiss),
                         scint = paste0(studyname,source,cint," %mis=",pmiss))
           
  
  
  # Create plot
  
  abs.outcome$source <- factor(abs.outcome$source,levels=c("Original","Deterministic","Imputation"),ordered=TRUE)
  abs.outcome$studyname <- factor(abs.outcome$studyname,ordered=TRUE)

  dotcols = c("#9DF19D","#C79DF1","#a6d8f0")
  barcols = c("#169C16","#59169C","#008fd5")
  

  
  dataplot <- abs.outcome%>%filter(!(is.na(included)&source=="Imputation"))%>%
              arrange(studyname,source)
  dataplot$scint<-factor(dataplot$scint,levels=dataplot$scint)
  
  plot<-ggplot(dataplot, aes(x=scint, y=mean, ymin=clb, ymax=cub, col=source, fill=source)) + 
    geom_linerange(size=2,position=position_dodge(width = 0.5)) +
    geom_point(size=1, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
    scale_fill_manual(values=barcols)+
    scale_color_manual(values=dotcols)+
    xlab("Study name")+ ylab("Absolute risk")+
    ggtitle(plottitle)+
    coord_flip() +
    scale_x_discrete(breaks=levels(factor(dataplot$scint)),labels = dataplot$cint)+ 
    facet_grid(allcint ~ ., switch = "y",scales="free")+
    theme(strip.placement = "outside")+
    guides(colour = guide_legend(reverse = T),fill = guide_legend(reverse = T))+
    theme(strip.text.y.left = element_text(angle = 0),axis.text.y = element_text(size = 6))+
    labs(caption = "*Excluded study in the imputation and in the total risk estimation")+
    labs(color="Data source",fill="Data source") 
  
  
  return(plot)}
               
