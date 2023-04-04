
#Functions necessary for objective 1

library(tidyr)
library(purrr)
library(magrittr)
library(ggplot2)
library(ggh4x)
library(metafor)

#######Functions#######
logit <- function(x){log(x/(1-x))} # better remove NA
logit_inv <- function(x){exp(x)/(1+exp(x))}


# Function to calculate the observed absolute risk, its se and CI (observed, observed_se,lower,upper) and the logit transformed absolute risk and its se (logit_obs, logit_obs_se)
f_sumarize_data <- function(data){
    alpha <- 0.05  
    sum_data<-data%>%summarise(obs_n = n(),
                               obs_sum = sum(outcome,na.rm=TRUE),
                               na_obs = sum(is.na(outcome),na.rm=TRUE),
                               included =first(Included))
    ptest = prop.test(sum_data$obs_sum, sum_data$obs_n, conf.level = alpha, correct = TRUE) # exact CI with Yates correction
    sum_data%<>% 
      # in the original scale
      mutate( observed = ptest$estimate,
              lower =   ptest$conf.int[[1]],
              upper = pmin(ptest$conf.int[[2]]))%>%
      mutate( observed_se = (upper-lower)/(2*qnorm(1-alpha/2)), # approx se on the original scale based on the CI
              observed_cor=pmax(pmin(ptest$estimate,0.999999),0.000001),# observed corrected when  obs=0 or obs=1 to logit transform  
              lower_cor=pmax(pmin(lower,0.999999),0.000001),
              upper_cor=pmax(pmin(upper,0.999999),0.000001))%>%          
      # Logit and se for estimate and observed          
      mutate( logit_obs = logit(observed_cor),
              logit_obs_se= (logit(upper_cor)-logit(lower_cor))/(2*qnorm(1-alpha/2)))
  
    return(sum_data)
  }
  
alpha=0.05

f_study_imp <- function(data, outcome_name) {
  sum_data<- data%>%
             dplyr::select(outcome=!!outcome_name,studyname,Included,'.imp')%>%
             nest(data= - c(studyname,'.imp'))%>%
             mutate(obs_sum = map(data,~f_sumarize_data(data=.x)))%>%
             select(-c(data))%>%
             unnest(col= c(obs_sum))
  return(sum_data)
}


#Function to pool the absolute risks per study per imputed dataset using rubins rules

f_pool <- function(est, logit_est, logit_se){
  
  alpha <- 0.05
  logit_est <- logit_est[!is.infinite(logit_est)&!is.na(logit_est)]
  logit_se <- logit_se[!is.infinite(logit_se)&!is.na(logit_se)]
  n <- length(logit_est)  #number observations
  mu <- mean(logit_est, na.rm = T) # pool estimate
  w_var <- mean(logit_se^2, na.rm = T) # within variance
  b_var <- var(logit_est, na.rm = T) # between variance
  t_var <- w_var + b_var + b_var/n # total variance
  t_se <- sqrt(t_var) # total standard error
  r <- (b_var + (b_var / n))/ w_var # relative increase variance due to missing values
  v <- (n - 1) * (1 + r^-1)^2 # degrees of freedom
  t <- qt(1-alpha/2, v) #t criteria
  observed <-  mean(est, na.rm = T)
  lower <- logit_inv(mu - t_se*t)
  upper <- logit_inv(mu + t_se*t)
  
  return(data.frame(cbind(observed = observed, lower = lower, upper = upper, logit_obs = mu, logit_obs_se = t_se)))
}


f_data_pool <- function(data) {
  data_pool<-data%>%
             filter(.imp>0)%>%  
             nest(data=-studyname)%>%  
             mutate(pool_sum=map(data,function(x){
               imp_n = nrow(x)
               pool_obs=f_pool(est=x$observed,logit_est=x$logit_obs,logit_se=x$logit_obs_se)}))%>%
    select(-c(data))%>%
            unnest(col= c(pool_sum))
  return(data_pool)}

# Function to calculate the global absolute risk
f_data_total <- function(data) {
  
  data_pool <- data%>%
               filter(included==1)
  
  fit.rma <- rma(yi=data_pool$logit_obs,
               sei=data_pool$logit_obs_se,
               method="REML",
               test= "knha")
  
  
  total <- data.frame(  studyname="TOTAL",
                         obs_n = sum(data_pool$obs_n,na.rm = T),
                         obs_sum = sum(data_pool$obs_sum,na.rm = T),
                         na_obs = sum(data_pool$na_obs,na.rm = T),
                         included = sum(data_pool$included,na.rm=T),
                         observed = logit_inv(as.numeric(fit.rma$beta)),
                         lower = logit_inv(as.numeric(fit.rma$ci.lb)),
                         upper = logit_inv(as.numeric(fit.rma$ci.ub)),
                         observed_se = (logit_inv(as.numeric(fit.rma$ci.ub))- logit_inv(as.numeric(fit.rma$ci.lb)))/(2*qnorm(1-alpha/2)),
                         observed_cor = NA,
                         logit_obs = as.numeric(fit.rma$beta),
                         logit_obs_se = as.numeric(fit.rma$se))
  
  datan <- data%>%rows_insert(total)
  return(datan)}
  

forest_plot_study<-function(data,outcome_name,plottitle){
  #Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
  inc.outcome.ind <- f_study_imp(data=data, outcome_name=outcome_name)
  
  #Estimate global for each dataset (raw and imputed using only the included datasets)
  inc.outcome <- inc.outcome.ind%>%
                 nest(data=-.imp)%>% 
                 mutate( tdata = map(data,~f_data_total(data=.x)))%>%
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
    select(c("studyname","observed","lower","upper","logit_obs","logit_obs_se","source","pmiss"))
  
  # For raw dataset
  abs.outcome.det <- inc.outcome%>%filter(.imp==0)%>%
    mutate( source = "Deterministic")%>%
    select(c("studyname","observed","lower","upper","logit_obs","logit_obs_se","source","pmiss"))
  
  
  # For imputed datases, pool information with Rubins rules
   
   abs.outcome.imp <-f_data_pool(data=inc.outcome)
   abs.outcome.imp$source <-"Imputation"
   abs.outcome.imp$pmiss <-0.00 
  
  # Merge both datasets
  
  abs.outcome <- dplyr::bind_rows(abs.outcome.imp, abs.outcome.det, abs.outcome.ori)%>%
                 left_join(abs.n, by=c('studyname'))%>%
    mutate(mean = observed*100,
           clb = lower*100,
           cub = upper*100,
           pmiss = sprintf("%.2f",pmiss))%>%
    mutate(cint = paste0(sprintf("%.2f(%.2f,%.2f)",mean, clb, cub)))%>%
    mutate(allcint = paste0(studyname," \n N=",N),
           cint = paste0(cint,", %mis=",pmiss),
           scint = paste0(source,cint,", %mis=",pmiss))
           
  
  
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
    theme(strip.text.y.left = element_text(angle = 0),axis.text.y = element_text(size = 6))
  
  return(plot)}
               
