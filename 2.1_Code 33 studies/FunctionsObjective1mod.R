
#0. General transformations ----
logit <- function(x) {log(x/(1-x))}
inv.logit <- function(x) {1/(1+exp(-x))}
inv.logit.SE<- function(logit.se, c) {logit.se * (c*(1-c))}

#1. Absolute risk functions ----
#1.1. Outcome at study level---- 
f.incidence <- function(data,outcome_name) {
  a<-summary(data[,get(outcome_name)])
  b<-prop.test(x = a[[2]], n = sum(a), correct = TRUE)
  prop<-b$estimate[[1]]
  ci<-b$conf
  return (list(incidence = prop,
               ci.l=ci[[1]],
               ci.u=ci[[2]]))
}


#1.2.Outcome per study and per imputed dataset ----
f.abs.perstudy<-function(data, outcome_name) {
  
  inc.outcome<-setDT(data)[, f.incidence(data=.SD, outcome_name=outcome_name), by = list(studyname,.imp)]

  #Replace 0's by 0.000001
  
  inc.outcome$incidence[inc.outcome$incidence==0]<-0.000001
  inc.outcome$ci.l[inc.outcome$ci.l==0]<-0.000001
  inc.outcome$ci.u[inc.outcome$ci.u==0]<-0.000001
  
  #Replace 1's by 0.999
  
  inc.outcome$incidence[inc.outcome$incidence==1]<-0.9999
  inc.outcome$ci.l[inc.outcome$ci.l==1]<-0.9999
  inc.outcome$ci.u[inc.outcome$ci.u==1]<-0.9999
  
  #Logit transformation
  inc.outcome$logit.incidence<-logit(inc.outcome$incidence)
  inc.outcome$logit.ci.l<-logit(inc.outcome$ci.l)
  inc.outcome$logit.ci.u<-logit(inc.outcome$ci.u)
  inc.outcome$logit.se<-(inc.outcome$logit.ci.u-inc.outcome$logit.ci.l)/(2*1.96)
  inc.outcome
}

#1.3. Pool ----
f.pool<-function(data,m){
  logit.abs <- mean(data$logit.incidence)
  within <- mean(data$logit.se^2)
  if(m > 1){
  between <- (1 + (1/m)) * var(data$logit.incidence)}else{
  between<-0}
  logit.var <- within + between
  logit.se <- sqrt(logit.var)
  logit.ci.lb <- logit.abs + qnorm(0.05/2)     * logit.se
  logit.ci.ub <- logit.abs + qnorm(1 - 0.05/2) * logit.se
  incidence<-inv.logit(logit.abs)*100
  ci.lb<-inv.logit(logit.ci.lb)*100
  ci.ub<-inv.logit(logit.ci.ub)*100
  return(list(logit.abs=logit.abs,within=within,between=between,
              logit.var=logit.var,logit.se=logit.se,logit.ci.lb=logit.ci.lb,logit.ci.ub=logit.ci.ub,
              incidence=incidence,ci.lb=ci.lb,ci.ub=ci.ub))
}

f.abs.poolrubin <-function(data) {
  data<-setDT(data)
  abs.outcome<-data[, f.pool(data=.SD, m <- length(unique(data$.imp))), by = studyname]
  return(abs.outcome)
}

#1.4. Overall ----
f.abs.2s.ma<-function(abs.outcome) {
  #Pool results: two-stage meta-analysis
  fit.rma<-rma(yi = logit.abs, sei = logit.se, method = "REML", test = "knha", 
               data = abs.outcome)
  pool.outcome<-data.frame(matrix(NA, nrow = 1, ncol = 9))
  colnames(pool.outcome)<-c("logit.abs","logit.se","logit.ci.lb","logit.ci.ub",
                            "abs.risk","ci.lb","ci.ub","pi.lb","pi.ub")
  pool.outcome$logit.abs<-fit.rma$beta
  pool.outcome$logit.se<-fit.rma$se
  pool.outcome$logit.ci.lb<-fit.rma$ci.lb
  pool.outcome$logit.ci.ub<-fit.rma$ci.ub
  pool.outcome$abs.risk<-inv.logit(pool.outcome$logit.abs[[1]])*100
  pool.outcome$ci.lb<-inv.logit(pool.outcome$logit.ci.lb)*100
  pool.outcome$ci.ub<-inv.logit(pool.outcome$logit.ci.ub)*100
  #Prediction interval
  PI<-PredInt(fit.rma)
  pool.outcome$pi.lb<-PI[1]*100
  pool.outcome$pi.ub<-PI[2]*100
  return(pool.outcome)
}
