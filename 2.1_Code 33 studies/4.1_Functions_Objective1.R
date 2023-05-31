# Functions needed for objective 1

# Load packages ----
## Data 
library(here)
library(tidyr)
library(purrr)
library(magrittr)
## Plots
library(ggplot2)
library(ggh4x)
## Modelling
library(metafor) # rma function for 2 step meta-analysis
library(lme4) # glmer function for 1 step meta-analysis
library(sandwich)# sandwich variance poisson
library(merDeriv) # sandwich variance poisson
library(misty) # center exposures

#' Function to transform proportions
#'
#' @param e  number of events
#' @param n  number of observations
#' @param t_type type of transformation ("logit", "arcsine")
#'
#' @return data.frame ( y = transformed point estimate, se = transformed standard error  )
#'

trans <- function(e, n, t_type){
  if (t_type=="logit") {
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
  se <- sqrt(v) # sure https://files.eric.ed.gov/fulltext/EJ1111448.pdf
  return(data.frame(y = y, se = se))
} 


#' Backtransform proportions or RR
#'
#' @param x  transformed value
#' @param t_type transformation type ("logit", "arcsine", "log")
#'
#' @return backtransformed value
#'

back_trans <- function(x, t_type){
  if (t_type == "logit") {
    y <- exp(x)/(1 + exp(x))
  }else if (t_type == "arcsine"){
    y <- (sin(x))^2
  }else{ # log
      y <- exp(x)
    }
  return(y)
} 



# Function to calculate the observed risk, its se and CI (observed, observed_se,lower,upper) and the transformed absolute risk and its se (tr_obs, tr_se)
#' Title
#'
#' @param data data.frame with variables outcome and exposure included
#' @param estimand type of estimand: Absolute risk "AR", Relative risk "RR" 
#' @param t_type type of transformation ("logit", "arcsine", "log")
#' @param correction only for relative risk: when number of events are close to zero or total number of observations ("Haldane","Hybrid","Sweeting","none") as described in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8254431/ and in doi:/10.1002/sim.1761
#'
#' @return data.frame with summary for the given dataset: estimate(observed, lower, upper) and transformed estimate (pointwise estimate: tr_obs, standard error: tr_se)

f_sumarize_data <- function(data,estimand,t_type,correction){
    alpha <- 0.05  
    sum_data <- data%>%summarise(n_obs = n(),
                                 n_ev = sum(outcome,na.rm=TRUE),
                                 n_na = sum(is.na(outcome),na.rm=TRUE),
                                 Included = first(Included),
                                 N=first(N))
    if (sum_data$n_obs != sum_data$n_na){  # not all are NA observations

      if (estimand == "RR"){  # Relative risk "RR"
   
      sum_datar <- data%>%
                  group_by(exposure)%>%
                  summarize(e = sum(outcome,na.rm=T),
                            ena =  sum(is.na(outcome)),
                            n = n(),
                            Included = first(Included),
                            N = first(N))%>%
                  mutate(n_obs = sum(n),
                         n_ev = sum(e,na.rm=T),
                         n_na = sum(n[is.na(exposure)|is.na(e)]))%>%
                  filter(!(is.na(exposure)))%>%
                  filter(ena != n)%>%
                  select(-ena)%>%
                  pivot_wider(names_from = exposure, values_from = c(e,n))
      
        if(nrow(sum_datar)>0){ # data table exist
            e1 <- ifelse(is.null(sum_datar$e_1),NA,sum_datar$e_1)
            n1 <- ifelse(is.null(sum_datar$n_1),NA,sum_datar$n_1)
            e0 <- ifelse(is.null(sum_datar$e_0),NA,sum_datar$e_0)
            n0 <- ifelse(is.null(sum_datar$n_0),NA,sum_datar$n_0)
            sum_data <- sum_datar
     
            if( any(is.na(c(e1,e0,n1,n0)))| e1 == 0 & e0 == 0){ # Studies with no events in either arm# cochrane section-10-4-4-2
               observed <- lower <- upper <- tr_obs <- tr_se <- NA
             
              }else{
              
              if(e1 == 0 | e0 == 0){ #only apply correction when it is required
                  c <- 0.5
                  if (correction == "Haldane"){
                    observed <- (e1+c)*(n0+2*c)/((n1+2*c)*(e0+c))
                    tr_se <- sqrt(1/(e1+c)-1/(n1+2*c)+1/(e0+c)-1/(n0+2*c))
                  }else if (correction == "Hybrid"){
                    observed <- (e1+c)*(n0+c)/((n1+2*c)*(e0+c))
                    tr_se <- sqrt(1/(e1+c)-1/(n1+2*c)+1/(e0+c)-1/(n0+c))
                  }else if (correction == "Sweeting"){
                    R <- n0/n1
                    c0 <- 1/(R+1)
                    c1 <- R/(R+1)
                    observed <- (e1+c1)*(n0+2*c0)/((n1+2*c1)*(e0+c0))
                    tr_se <- sqrt(1/(e1+c1)-1/(n1+2*c1)+1/(e0+c0)-1/(n0+2*c0))
                  }else{ # no correction , correction ="none"
                    observed <- (e1)*(n0)/((n1)*(e0))
                    tr_se <- sqrt(1/(e1)-1/(n1)+1/(e0)-1/(n0))
                  }
              }else{ # no correction is applied
                observed <- (e1)*(n0)/((n1)*(e0))
                tr_se <- sqrt(1/(e1)-1/(n1)+1/(e0)-1/(n0))
                }
             }
            }else{ #data.table does not exist
              observed <- lower <- upper <- tr_obs <- tr_se <- NA
            }
  
        sum_data%<>%
          mutate( observed = observed, 
                  tr_obs = ifelse(observed == 0,NA,log(observed)),
                  tr_se = ifelse(is.infinite(tr_se),NA,tr_se))%>%
          mutate( lower = exp(tr_obs - qnorm(1-alpha/2)*tr_se),
                  upper = exp(tr_obs + qnorm(1-alpha/2)*tr_se))
     
      }else{ # Absolute risk "AR"
   
          ptest <- prop.test(sum_data$n_ev, sum_data$n_obs, conf.level = alpha, correct = TRUE) # exact CI with Yates correction
          trans_observed <- trans(e=sum_data$n_ev, n=sum_data$n_obs, t_type = t_type)
          
          sum_data%<>%
            mutate( observed = ptest$estimate, 
                    tr_obs = trans_observed$y,
                    tr_se = trans_observed$se,
                    lower = ptest$conf.int[[1]],
                    upper = ptest$conf.int[[2]])
          
       }
    }
    
    return(sum_data)
  }
  


#Function to pool the absolute risks per study per imputed dataset using rubins rules

#' @param tr_est  transformed point estimate
#' @param tr_se   transformed standard error
#' @param t_type  type of transformation ("logit", "arcsine", "log")
#'
#' @return data.frame with pooled information in original scale(observed, lower, upper) and transformed scale (tr_obs, tr_se)

f_pool <- function(tr_est, tr_se, t_type){
  
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
  observed <- back_trans(mu, t_type = t_type) # mean(est, na.rm = T)
  lower <- back_trans(mu - t_se*t, t_type = t_type)
  upper <- back_trans(mu + t_se*t, t_type = t_type)
  
  return(data.frame(cbind(observed = observed, lower = lower, upper = upper, tr_obs = mu, tr_se = t_se)))
}


#' Function to get the pooled estimates in all datasets
#'
#' @param data data.frame with all estimates to be pooled
#' @param t_type type of transformation ("logit", "arcsine", "log")
#'
#' @return data.frame with pool information


f_data_pool <- function(data,t_type) {
  data_pool<-data%>%
             filter(.imp>0)%>%  
             nest(data = -studyname)%>%  
             mutate(pool_sum = map(data,function(x){
               imp_n = nrow(x)
               pool_obs = f_pool(tr_est=x$tr_obs,tr_se=x$tr_se,t_type=t_type)
               data.frame(n_obs = mean(x$n_obs,na.rm=T),
                          n_ev = mean(x$n_ev,na.rm=T),
                          Included = mean(x$Included,na.rm=T),
                          nobs_m = mean(x$nobs_m,na.rm=T),
                          nev_m = mean(x$nev_m,na.rm=T),
                          Included_m = mean(x$Included_m,na.rm=T),
                          N = first(x$N),
                          pool_obs)}))%>%
    select(-c(data))%>%
            unnest(col = c(pool_sum))
  return(data_pool)}


#' Total relative risk by using 1-step method via generalized linear mixed models
#'
#' @param data original dataset that includes the imputed and studyname information
#' @param mod_type only for Relative risk, type of model to estimate the model: "binomial","poisson", as given here https://academic.oup.com/aje/article/189/6/508/5812650 when binomial model does not converge we use poisson estimates. 
#'
#' @return Estimate of total relative risk in original scale(pointwise ="observed", lower, upper) and in transformed scale (tr_obs,tr_se)

RR_1step <- function(data, mod_type){
  alpha <- 0.05
  summary  <- data%>%summarize(n_obs=n(),
                   n_ev = sum(outcome,na.rm=TRUE),
                   n_na = sum(Included[is.na(outcome)|is.na(exposure)]),
                   Included = length(unique(studyname)))
  
  data_mod <- data%>% # Included in the model 
    filter(Included == 1)%>%
    filter(!is.na(outcome)&!is.na(exposure))
  
  summary%<>%
    mutate(nobs_m = nrow(data_mod),
           nev_m = sum(data_mod$outcome,na.rm=T),
           Included_m =  length(unique(data_mod$studyname)))

  fit <- NA
  if(mod_type == "binomial"){
    
    data$exposure<-tryCatch(exp={misty::center(data$exposure,type="CWC",cluster=data$studyname)},
                            error=function(e){data$exposure})
    
    fit <- tryCatch( expr={lme4::glmer(outcome ~ exposure + (1|studyname),data=data,family = binomial("log"))},
                     error = function(e){NA})
    
    if(!is.na(fit)){
    tr_obs <- tryCatch(exp={summary(fit)$coefficients[2,1]},error = function(e){NA})
    tr_se <- tryCatch(exp={summary(fit)$coefficients[2,2]},error = function(e){NA})}
  }
   
  if(mod_type == "poisson"| is.na(fit)){ # default option in case binomial do not converge
     fit <- tryCatch( expr={lme4::glmer(outcome ~ exposure + (1|studyname),data=data,family = poisson("log"))},
                      error = function(e){NA})
     if(!is.na(fit)){
       tr_obs <- tryCatch(exp={summary(fit)$coefficients[2,1]},error = function(e){NA})
       tr_se <- tryCatch(exp={sqrt(sandwich(fit, bread = bread(fit, full = TRUE),mean = meat(fit, level = 2))[2,2])},# robust variance
                         error = function(e){NA}) }}
    
  if(is.na(fit)){
      tr_obs <- NA
      tr_se <- NA
  }

    summary%>%
      mutate( observed = exp(tr_obs),
              lower = back_trans(tr_obs - qnorm(1 - alpha/2) * tr_se, t_type="log"),
              upper = back_trans(tr_obs + qnorm(1 - alpha/2) * tr_se,t_type="log"),
              tr_obs = tr_obs,
              tr_se = tr_se)
      }

    
#' Total absolute risk by using 1-step method via generalized linear mixed models
#'
#' @param data summarized data
#' @param t_type type of transformation ("logit", "arcsine", "log")
#'
#' @return Estimate of total absolute risk in original scale(pointwise ="observed", lower, upper) and in transformed scale (tr_obs,tr_se)

AR_1step <- function(data,t_type){
  alpha <- 0.05
  summary  <- data%>%summarize(n_obs = sum(n_obs,na.rm=TRUE),
                               n_ev = sum(n_ev,na.rm=TRUE),
                               n_na = sum(n_na,na.rm=TRUE),
                               Included = sum(Included,na.rm=TRUE))
  
  data_mod <- data%>% # Included in the model 
    filter(Included == 1)%>%
    filter(!is.na(n_obs)&!is.na(n_ev))
  
  summary%<>%
    mutate(nobs_m = sum(data_mod$n_obs,na.rm=T),
           nev_m = sum(data_mod$n_ev,na.rm=T),
           Included_m = sum(data_mod$Included,na.rm=T))

  
  
  fit <- tryCatch( expr={lme4::glmer(cbind(n_ev, n_obs - n_ev) ~ 1 + (1 | studyname), data = data,family = binomial(link = "logit"))},
                   error = function(e){
                     tryCatch( expr = {lme4::glmer(cbind(n_ev, n_obs - n_ev) ~ 1 + (1 | studyname), data = data,family = binomial(link = "logit"),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))},
                               error = function(e){NA})
                               })
                  
  
  if(!is.na(fit)){
    tr_obs <- as.numeric(summary(fit)$coefficients[1])
    tr_se <- as.numeric(summary(fit)$coefficients[2])
  }else{
    tr_obs <- NA
    tr_se <- NA
  }
  
  summary%>%
    mutate( observed = back_trans(tr_obs,t_type),
            lower = back_trans(tr_obs - qnorm(1-alpha/2)*tr_se,t_type),
            upper = back_trans(tr_obs + qnorm(1-alpha/2)*tr_se,t_type),
            tr_obs = tr_obs,
            tr_se = tr_se)
}
    


# Total risk (absolute and relative) by using 2-step method via the rma funcion in the package metafor
#'
#' @param data summarized data of each imputed-studyname dataset
#' @param t_type type of transformation ("logit", "arcsine", "log")
#'
#' @return Estimate of total absolute risk in original scale(pointwise ="observed", lower, upper) and in transformed scale (tr_obs,tr_se)

R_2step <- function(data, t_type){
  alpha<-0.05
  data_pool <- data%>%
               filter(Included==1)
  
  data_mod <- data%>% # Included in the model 
    filter(Included==1)%>%
    filter(!is.na(tr_obs)&!is.na(tr_se))
  
   data_2step <- data.frame( n_obs = sum(data_pool$n_obs,na.rm = T),
                             n_ev = sum(data_pool$n_ev,na.rm = T),
                             n_na = sum(data_pool$n_na,na.rm = T),
                             Included = sum(data_pool$Included,na.rm=T),
                             nobs_m = sum(data_mod$n_obs,na.rm=T),
                             nev_m = sum(data_mod$n_ev,na.rm=T),
                             Included_m = sum(data_mod$Included,na.rm=T))
  
   # knha only applied in datasets with >10 studies DOI: 10.1002/jrsm.1316
    fit.rma <- tryCatch(
      expr = { metafor::rma(yi=data_mod$tr_obs, sei=data_mod$tr_se, method="REML",test= "knha")},
      error = function(e){
        tryCatch( expr = {metafor::rma(yi=data_mod$tr_obs, sei=data_mod$tr_se, method="REML",test= "knha",control=list(stepadj=0.5))},
                  error = function(e){NA})
      })
        
  
  if(all(is.na(fit.rma))){
    observed<-lower<-upper<-tr_obs<-tr_se<-NA
  }else{
    observed = back_trans(as.numeric(fit.rma$beta),t_type)
    lower = back_trans(as.numeric(fit.rma$ci.lb),t_type)
    upper = back_trans(as.numeric(fit.rma$ci.ub),t_type)
    tr_obs = as.numeric(fit.rma$beta)
    tr_se = as.numeric(fit.rma$se)
  
  }
  

  data_2step%>%mutate( observed = observed,
                        lower = lower,
                        upper = upper,
                        tr_obs = tr_obs,
                        tr_se = tr_se)
}
  
  
  

#'  Function to get the summarized risk (relative or absolute) estimates (sumamary table "sum_tdata" and plot "plot) from a dataset with imputed and studyname information 
#'
#' @param data dataset with imputed and studyname information where the information source ( original data "original", data after deterministic imputation "Deterministic", complete data after imputation "Imputation") is provided 
#' @param outcome_name name of the outcome variable to estimate
#' @param exposure_name name of the exposure variable, in case of absolute risk = NA
#' @param estimand type of estimand: Absolute risk "AR", Relative risk "RR"
#' @param plottitle title of the plot
#' @param t_type type of transformation ("logit", "arcsine", "log")
#' @param correction only for relative risk: when number of events are close to zero or total number of observations ("Haldane","Hybrid","Sweeting","none") as described in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8254431/ and in doi:/10.1002/sim.1761
#' @param mod_type only for relative risk, type of model to estimate the model: "binomial","poisson", as given here https://academic.oup.com/aje/article/189/6/508/5812650 when binomial model does not converge we use poisson estimates. 
#' @param dupper maximum upper confidence interval value to display in the plot, default NA.
#'
#' @return list with "plot_all", plot that shows the risk estimates at study and data information level.
#'                   "plot_imp", plot that shows the risk estimates at study level only for the imputed dataset.
#'                   "sum_tdata"  summary table at study and information level
#'                   "data_stuimp" summary table at study and imputation level (.imp=-1 "original data", .imp=0 "data after deterministic imputation", .imp>0 "imputed datasets")

Rpool_studies <- function(data, outcome_name, exposure_name = NA, estimand, plottitle, t_type, correction = "Haldane", mod_type ="binomial", dupper = NA){
  
  # Select only variables we need.
  if( is.na(exposure_name)){ # AR
    data%<>%dplyr::select(outcome=!!outcome_name, studyname, Included, N,'.imp')
   }else{ #RR
    data%<>%dplyr::select(outcome=!!outcome_name, exposure=!!exposure_name, studyname, Included,N,'.imp')}
  
  
  #Summarize the risk in every study separate and in every imputed dataset
    data_sum <- data%>%nest(data= - c(studyname,'.imp'))%>%
                   mutate(obs_sum = map(data,~f_sumarize_data(data = .x, estimand = estimand, t_type = t_type, correction = correction)))%>%
                   select(-c(data))%>%
                   unnest(col= c(obs_sum))%>%
                   mutate_all(~ifelse(is.nan(.)|is.infinite(.), NA, .))
    
  
  #Estimate global estimator for each dataset via 2-step estimator
    data_2step <- data_sum %>%
                  nest(data = -'.imp')%>%
                  mutate(total2 = map(data,~R_2step(data =.x, t_type = t_type)))%>%
                  select(-c(data))%>%
                  unnest(col= c(total2))%>%
                  mutate(studyname="TOTAL-2step")
    
    

  #Estimate global estimator for each dataset via 1-step estimator
    if (estimand =="RR"){
      # 1-step
      data_1step <- data%>%
                    filter(Included==1)%>%
                    nest(data = -.imp)%>% 
                    mutate(pool = map( data, ~RR_1step(data=.,mod_type)))%>%
                    select(-c(data))%>%
                    unnest(col= c(pool))%>%
                    mutate(studyname="TOTAL-1step")
   
    }else{ # estimand ="AR"
      
      data_1step <- data_sum%>%
                   filter(Included==1)%>%
                   nest(data = -.imp)%>% 
                   mutate(pool = map( data, ~AR_1step(data=.,t_type)))%>%
                   select(-c(data))%>%
                   unnest(col= c(pool))%>%
                   mutate(studyname="TOTAL-1step")
    
    }
        
    
    # Rbind the summary dataset with the values of global estimates of the 1_step and 2_step methods --    
    tdata_sum <- dplyr::bind_rows(data_sum, data_1step,data_2step)%>%
                    mutate(pmiss = n_na/n_obs*100,
                           outcome = outcome_name,
                           exposure = exposure_name,
                           estimand = estimand) 
        
        
    # For original dataset
    tdata_ori <- tdata_sum%>%filter(.imp == -1)%>%
      mutate( source = "Original")
      

    # For dataset after applying deterministic imputation
    tdata_det <- tdata_sum%>%filter(.imp == 0)%>%
      mutate( source = "Deterministic")
    
    
    # For imputed datases, pool information with Rubins rules
    
    tdata_imp <- f_data_pool(data = tdata_sum, t_type = t_type)
    tdata_imp$source <- "Imputation"
    tdata_imp$pmiss <- 0.00   
    tdata_imp$estimand <- estimand
        
        

  
  # Merge both datasets
  
  sum_tdata <- dplyr::bind_rows( tdata_ori, tdata_det, tdata_imp)%>%
                 mutate(estimand=estimand)%>%
                 mutate_all(~ifelse(is.nan(.), NA, .))%>%
                 mutate( mean = ifelse(estimand=="RR",observed,observed*100),
                         clb = ifelse(estimand=="RR",lower,lower*100),
                         cub = ifelse(estimand=="RR",upper,upper*100),
                         pmiss = sprintf("%.2f",pmiss))%>%
                 select(-.imp)
  
  plot_data <- sum_tdata %>%
               mutate(Sincluded = ifelse(Included==0,"*",""),
                      cint0 = ifelse(is.na(mean),"   ",paste0(sprintf("%.2f(%.2f,%.2f)",mean, clb, cub))),
                      cint = paste0(cint0,","))%>%
                  mutate(allcint = paste0(studyname,Sincluded,ifelse(is.na(N),"",paste0(" \n N=",N))),
                         pcint = paste0(cint," %mis=",pmiss),
                         scint = paste0(studyname,source,cint," %mis=",pmiss))
           
  
  # Create plot
  
  plot_data$source <- factor(plot_data$source,levels=c("Original","Deterministic","Imputation"),ordered=TRUE)
  plot_data$studyname <- factor(plot_data$studyname,ordered=TRUE)

  dotcols = c("#9DF19D","#C79DF1","#a6d8f0")
  barcols = c("#169C16","#59169C","#008fd5")
  

  
  dataplot <- plot_data%>%
              arrange(studyname,source)
  dataplot$scint <- factor(dataplot$scint,levels=dataplot$scint)
  xlabv <- ifelse(estimand == "RR","Relative risk","Absolute risk")
  int<-ifelse(estimand == "RR",1,0)
  
  if(!is.na(dupper)){ # if a upper limit to display is given
    dataplot%<>%mutate(cub = ifelse(cub > dupper, dupper, cub),
                       mean =ifelse(mean > dupper, dupper, mean))
  }
  
  
  
  # Plot all source of information 
  
  plot_all <- ggplot(dataplot, aes(x=scint, y=mean, ymin=clb, ymax=cub, col=source, fill=source)) + 
    geom_linerange(size=2,position=position_dodge(width = 0.5)) +
    geom_point(size=1, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
    scale_fill_manual(values=barcols)+
    scale_color_manual(values=dotcols)+
    xlab("Study name")+ ylab(xlabv)+
    ggtitle(plottitle)+
    coord_flip() +
    geom_hline(yintercept=int,linetype="dashed", color = "black", size=0.2)+
    scale_x_discrete(breaks=levels(factor(dataplot$scint)),labels = dataplot$pcint)+ 
    facet_grid(allcint ~ ., switch = "y",scales="free")+
    theme(strip.placement = "outside")+
    guides(colour = guide_legend(reverse = T),fill = guide_legend(reverse = T))+
    theme(strip.text.y.left = element_text(angle = 0,size=5),axis.text.y = element_text(size = 5))+
    labs(caption = "*Excluded study in the imputation and in the total risk estimation")+
    labs(color="Data source",fill="Data source")+
    theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
  
  
  
# Plot only imputed datasets
    dataplot%<>%
               filter(source=="Imputation")%>%
               filter(!is.na(mean))%>%
               mutate(pcint=cint,
                      type= ifelse(studyname%in%c("TOTAL-1step","TOTAL-2step"),"0","1"),
                      N = ifelse(is.na(N),"",paste0("N=",N)),
                      y=interaction(N,studyname, sep = "&"))

    dataplot$y <- factor(dataplot$y, levels=rev(levels(dataplot$y)))
  

    
  plot_imp <- ggplot(dataplot, aes(x = mean, y = y)) +
    theme_classic()+
    geom_vline(xintercept=int,linetype="dashed", color = "red", size=0.5)+
    scale_color_manual(values=c("blue","black"))+
    scale_shape_manual(values=c(18,15))+
    scale_size_manual(values=c(5,2))+
    geom_point(aes(x=mean,shape=type,color=type,size=type)) +
    geom_linerange(aes(xmin=clb, xmax=cub,color=type))+
    guides(y = guide_axis_nested(delim = "&"),
           y.sec = guide_axis_manual(
             breaks = dataplot$y,
             labels = dataplot$cint0)) +
    theme(axis.text.y.left = element_text(margin = margin(r = 5, l = 5),size=7,hjust=0),
          axis.text.y.right = element_text(margin = margin(r = 5, l = 5),size=7,hjust=1),
          ggh4x.axis.nesttext.y = element_text(margin = margin(r = 5, l = 5),face="bold",hjust=0),
          ggh4x.axis.nestline = element_blank())+
    ylab("Study name")+ xlab(xlabv)+
    ggtitle(plottitle)+
    theme(legend.position="none",plot.title = element_text(hjust = 0.5))  
    


  return(list(plot_all = plot_all,
              plot_imp = plot_imp,
              sum_tdata = sum_tdata,
              data_stuimp = tdata_sum))}




#' Print all the output of objective 1
#'
#' @param outcome_name  Name of the outcome variable
#' @param title general title to display in graphs
#' @param title upper limit to display in RR graph
#'
#' @return all output is saved it on the folder here("6_Tables_graphs","Objective1")

print_obj1 <- function(outcome_name,gentitle,dupperRR=NA,data_all,data_zika,data_nozika){
  
  # Absolute risk ----
  ARall <- Rpool_studies(data = data_all, # for estimates only on zika+ mother use here: data_zika or zika- mother: data_nozika
                         outcome_name = outcome_name, # it can be used also for: "microcephaly_bin_postnatal","microcephaly_bin_fet","ch_czs","who_czs","neuroabnormality","nonneurologic","miscarriage","loss","efdeath","lfdeath"
                         exposure_name = NA,
                         estimand = "AR",
                         plottitle = paste0(gentitle),
                         t_type = "logit",
                         dupper = NA)
  
  ggsave(filename=here("6_Tables_graphs","Objective1",paste0(outcome_name,"_ARall.jpg")), plot=ARall$plot_all, width=12, height=9, units="in")
  ggsave(filename=here("6_Tables_graphs","Objective1",paste0(outcome_name,"_ARimp.jpg")), plot=ARall$plot_imp, width=12, height=9, units="in")
  
  
  ARpos <- Rpool_studies(data = data_zika, # for estimates only on zika+ mother use here: data_zika or zika- mother: data_nozika
                         outcome_name = outcome_name, # it can be used also for: "microcephaly_bin_postnatal","microcephaly_bin_fet","ch_czs","who_czs","neuroabnormality","nonneurologic","miscarriage","loss","efdeath","lfdeath"
                         exposure_name = NA,
                         estimand = "AR",
                         plottitle = paste0(gentitle," in ZIKV+ pregnant woman"),
                         t_type = "logit",
                         dupper = NA)
  
  ggsave(filename=here("6_Tables_graphs","Objective1",paste0(outcome_name,"_pos_ARall.jpg")), plot=ARpos$plot_all, width=12, height=9, units="in")
  ggsave(filename=here("6_Tables_graphs","Objective1",paste0(outcome_name,"_pos_ARimp.jpg")), plot=ARpos$plot_imp, width=12, height=9, units="in")
  
  ARneg <- Rpool_studies(data = data_nozika, # for estimates only on zika+ mother use here: data_zika or zika- mother: data_nozika
                         outcome_name = outcome_name, # it can be used also for: "microcephaly_bin_postnatal","microcephaly_bin_fet","ch_czs","who_czs","neuroabnormality","nonneurologic","miscarriage","loss","efdeath","lfdeath"
                         exposure_name = NA,
                         estimand = "AR",
                         plottitle = paste0(gentitle," in ZIKV- pregnant woman"),
                         t_type = "logit",
                         dupper = NA)
  
  ggsave(filename=here("6_Tables_graphs","Objective1",paste0(outcome_name,"_neg_ARall.jpg")), plot=ARneg$plot_all, width=12, height=9, units="in")
  ggsave(filename=here("6_Tables_graphs","Objective1",paste0(outcome_name,"_neg_ARimp.jpg")), plot=ARneg$plot_imp, width=12, height=9, units="in")
  
  # Relative risk ----
  RRall<- Rpool_studies(data = data_all,
                        outcome_name = outcome_name,
                        exposure_name = "zikv_preg",
                        estimand = "RR",
                        plottitle = paste0(gentitle),
                        t_type = "log",
                        mod_type ="binomial",
                        correction = "Haldane",
                        dupper = dupperRR)
  
  ggsave(filename=here("6_Tables_graphs","Objective1",paste0(outcome_name,"_RRall.jpg")), plot=RRall$plot_all, width=12, height=9, units="in")
  ggsave(filename=here("6_Tables_graphs","Objective1",paste0(outcome_name,"_RRimp.jpg")), plot=RRall$plot_imp, width=12, height=9, units="in")
  
  
  ALLlist<-list(ARall=ARall,ARpos= ARpos,ARneg= ARneg,RRall=RRall)

  save(ALLlist,file=here("6_Tables_graphs","Objective1",paste0(outcome_name,"_list.RData")))
}

               
