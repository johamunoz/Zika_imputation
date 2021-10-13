source("~/Run_zika/Functions_zika.R")

############# SET SEED ###############
args <- commandArgs(trailingOnly=TRUE)
imp <- as.numeric(args[1])

wd <-  "/home/julius_te/jmunoz/Results_z1/"
fnsuffix <- ".RData"
fname <- paste0(wd, "Imp.", imp, fnsuffix)

if(file.exists(fname)) {
  stop(paste("Result already available for simulation", imp))
} 

# Output HPC
 output<- fun_run(imp=imp,data=data,pred=pred,meth=meth,post=post,maxit=20)
 save(output, file = fname)

