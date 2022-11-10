

fit1 <- glmer(outcome ~ zikv_preg_cent + age_cent + educ_cent + maritalstat_cent + 
                bmi_cent + drugs_bin_cent + alcohol_cent +
                (1+zikv_preg_cent | studyname_fac),
              data=d, family = binomial)



f.1ma.cov<-function(data, outcome_name,family.link) {
  for (i in 1:length(unique(data$.imp))) {
    d<-data[data$.imp==i,]
    d$outcome <- as.numeric(levels(d[, outcome_name]))[d[, outcome_name]]
    fit1 <- glmer(outcome ~ zikv_preg_cent + age_cent + educ_cent + maritalstat_cent + 
                    bmi_cent + drugs_bin_cent + alcohol_cent + (1+zikv_preg_cent | studyname_fac),
                  data=d, family = family.link) #family.link can be binomial(link = "log"), binomial, poisson or others
    #Store results
    if(i==1) {fit1.coef<-summary(fit1)$coefficients[-1,(1:2)]}
    if(i>1) {fit1.coef<-rbind(fit1.coef,summary(fit1)$coefficients[-1,(1:2)])}
  }
  colnames(fit1.coef)<-c("log.rr","log.se")
  rownames(fit1.coef)<-NULL
  return(fit1.coef)
}
f.1ma.cov(data,"microcephaly_bin_birth",poisson)
