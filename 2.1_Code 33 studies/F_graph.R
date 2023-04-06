
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
  }else { # type arcsine
    p <- e/n
    y <- asin(sqrt(p))
    v <- 1/(4 * n)
  }
  se <- sqrt(v/n)
  return(data.frame(y = y,se = se))
} 


back_trans <- function(x,type){
  if (type == "logit") {
    y <- exp(x)/(1 + exp(x))
  }else if (type == "arcsine"){
    y <- (sin(x))^2
  }else{ # log
      y <- exp(x)
    }
  return(y)
} 

val_na <- function(x){
  if( length(x)==0){
    y <- 0
   }else if(is.infinite(x)|is.na(x)){
      y <- 0
   }else{
    y<-x
  }
     
}



# Function to calculate the observed absolute risk, its se and CI (observed, observed_se,lower,upper) and the transformed absolute risk and its se (tr_obs, tr_se)
f_sumarize_data <- function(data,estimand,type,correction,nacount){
    alpha <- 0.05  
    sum_data <- data%>%summarise(obs_n = n(),
                               obs_sum = sum(outcome,na.rm=TRUE),
                               na_obs = sum(is.na(outcome),na.rm=TRUE),
                               included =first(Included))
    
    if (sum_data$obs_n == sum_data$na_obs){
      observed <- lower <- upper <- tr_obs <- tr_se <- NA
     
    }else{
    
      if (estimand == "RR"){  # Relative Risk "RR"
        
        otable <- as.matrix(table(exp=data$exposure,out=data$outcome))
        colnam <- colnames(otable)
        
        # To count NA values 
        natable <- table(exp = data$exposure, out = data$outcome, useNA = "always")
        colnames(natable) <- paste0("out",colnames(natable))
        rownames(natable) <- paste0("exp",rownames(natable))
        
        na_obs <- natable["expNA","outNA"]
        
        if( nacount ==TRUE){ # if we assume that all NA in outcome are 0 
         
          natable <- natable[rownames(natable)!="expNA",]
          if(is.null(dim(natable))){
            newcol <- sum(natable[c("out0","outNA")],na.rm = T)
            } else{
              newcol <- as.numeric(rowSums( natable[,c("out0","outNA")],na.rm=TRUE)) 
            }
         
            if(is.element("0",colnam)){
              otable[,"0"] <- newcol
            }else{
              otable <- cbind(newcol,otable)
              colnam <- c("0",colnam)
              colnames(otable) <- colnam
            }
        }
        
        sumtab <- data.table(addmargins(otable))
        x1 <- val_na(sumtab[exp==1&out==1,N])
        n1 <- val_na(sumtab[exp==1&out=="Sum",N])
        x2 <- val_na(sumtab[exp==0&out==1,N])
        n2 <- val_na(sumtab[exp==0&out=="Sum",N])
        N <- n1+n1
        noev <- length(sumtab[exp=="Sum"&out==1,N]) == 0
        
        if( n1 == 0 | n2 == 0| noev){
          observed <- lower <- upper <- tr_obs <- tr_se <- NA
        }else{
          c <- 0.5
          if (correction == "Haldan"){
            observed <- (x1+c)*(n2+2*c)/((n1+2*c)*(x2+c))
            tr_se <- sqrt(1/(x1+c)-1/(n1+2*c)+1/(x2+c)-1/(n2+2*c))
          }else if (correction == "Hybrid"){
            observed <- (x1+c)*(n2+c)/((n1+2*c)*(x2+c))
            tr_se <- sqrt(1/(x1+c)-1/(n1+2*c)+1/(x2+c)-1/(n2+c))
          }else{ # no correction correction ="none"
            observed <- (x1)*(n2)/((n1)*(x2))
            tr_se <- sqrt(1/(x1)-1/(n1)+1/(x2)-1/(n2))
          }
        }
        
        tr_obs <- ifelse(observed==0,NA,log(observed))
        lower <- exp(tr_obs-qnorm(1-alpha/2)*tr_se)
        upper <- exp(tr_obs+qnorm(1-alpha/2)*tr_se)
        tr_se <- ifelse(is.infinite(tr_se),NA,tr_se)
        obs_sum <- obs_n - na_obs
        sum_data %<>%
          mutate( obs_sum = obs_sum, # number of observations with complete outcome and exposure.
                  na_obs = na_obs, # all observations with at least one information missing
                  )   
     
      }else{ # Absolute risk "AR"
   
          ptest <- prop.test(sum_data$obs_sum, sum_data$obs_n, conf.level = alpha, correct = TRUE) # exact CI with Yates correction
          trans_observed <- trans(e=sum_data$obs_sum, n=sum_data$obs_n, type = type)
          observed = ptest$estimate
          lower =   ptest$conf.int[[1]]
          upper =   ptest$conf.int[[2]]
          tr_obs = trans_observed$y
          tr_se = trans_observed$se
       }
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
  


f_study_imp <- function(data, outcome_name, exposure_name, estimand,type,correction,nacount) {
  
  if( is.na(exposure_name)){
    sum_data<- data%>%
      dplyr::select(outcome=!!outcome_name,studyname,Included,'.imp')
  }else{
    sum_data<- data%>%
      dplyr::select(outcome=!!outcome_name,exposure=!!exposure_name,studyname,Included,'.imp')
  }
  
  sum_data%<>%nest(data= - c(studyname,'.imp'))%>%
             mutate(obs_sum = map(data,~f_sumarize_data(data = .x,estimand = estimand,type = type, correction =correction, nacount = nacount)))%>%
             select(-c(data))%>%
             unnest(col= c(obs_sum))%>%
             mutate_all(~ifelse(is.nan(.)|is.infinite(.), NA, .))

  return(sum_data)
}


#Function to pool the absolute risks per study per imputed dataset using rubins rules

f_pool <- function(est, tr_est, tr_se, type){
  
  alpha <- 0.05
  mean <- mean(est, na.rm = T) 
  n <- length(tr_est)  #number observations
  mu <- mean(tr_est, na.rm = T) # pool estimate
  w_var <- mean(tr_se^2, na.rm = T) # within variance
  b_var <- var(tr_est, na.rm = T) # between variance
  t_var <- w_var + b_var + b_var/n # total variance
  t_se <- sqrt(t_var) # total standard error
  r <- (b_var + (b_var / n))/ w_var # relative increase variance due to missing values
  v <- (n - 1) * (1 + r^-1)^2 # degrees of freedom
  t <- qt(1-alpha/2, v) #t criteria
  observed <- back_trans(mu, type = type) # mean(est, na.rm = T)
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
  
  fit.rma<-tryCatch(
    expr = {rma(yi=data_pool$tr_obs, sei=data_pool$tr_se, method="REML",test= "knha")},
    error = function(e){ 
      rma(yi=data_pool$tr_obs, sei=data_pool$tr_se, method="REML",test= "knha",control=list(stepadj=0.5))})
  
  
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
  

forest_plot_study<-function(data,outcome_name,exposure_name = NA, estimand, plottitle, type, correction = "none", nacount = TRUE){
  #Calculate the absolute risk of the outcome in every study separate and in every imputed dataset
  inc.outcome.ind <- f_study_imp(data=data, outcome_name=outcome_name, exposure_name=exposure_name, estimand = estimand, type = type, correction= correction, nacount= nacount)
  
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
                 mutate(estimand=estimand)%>%
                  mutate(included = ifelse(included==0,"*",""),
                         mean = ifelse(estimand=="RR",observed,observed*100),
                         clb = ifelse(estimand=="RR",lower,lower*100),
                         cub = ifelse(estimand=="RR",upper,upper*100),
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
  dataplot$scint <- factor(dataplot$scint,levels=dataplot$scint)
  
  plot <- ggplot(dataplot, aes(x=scint, y=mean, ymin=clb, ymax=cub, col=source, fill=source)) + 
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
               
