library(mice)
library(here)
library(data.table)
library(dplyr)
library(ggplot2)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])[[1]]
}


upfunc<-function(path){
  load(path)
  output<-output[[1]]
  }

mergeMice <- function (imp) {
  merged_imp <- NULL
  for (n in 1:length(imp)){
    if (is.null(merged_imp)){
      merged_imp <- imp[[n]]
    }else{
      counter <- merged_imp$m
      # Update counter
      merged_imp$m <- merged_imp$m + imp[[n]]$m
      # Rename chains
      dimnames(imp[[n]]$chainMean)[[3]] <-sprintf("Chain %d", (counter + 1):merged_imp$m)
      dimnames(imp[[n]]$chainVar)[[3]] <-sprintf("Chain %d", (counter + 1):merged_imp$m)
      # Merge chains
      merged_imp$chainMean <- abind::abind(merged_imp$chainMean,imp[[n]]$chainMean)
      merged_imp$chainVar <- abind::abind(merged_imp$chainVar, imp[[n]]$chainVar)
      
      for (nn in names(merged_imp$imp)){
        # Non-imputed variables are not in the data.frame format but are null
        
        if (!is.null(imp[[n]]$imp[[nn]])){
          
          colnames(imp[[n]]$imp[[nn]]) <- (counter + 1):merged_imp$m
          merged_imp$imp[[nn]] <- cbind(merged_imp$imp[[nn]],imp[[n]]$imp[[nn]])
          
        }
      }
    }
  }
  return(merged_imp)
}

path<-file.path(here("3_Output_data"))
#path<-file.path(here('Documents','GitHub','Zika_imputation',"3_Output_data"))
list_files <- list.files(path,pattern="^mice*",full.names=TRUE)
imp <- lapply(list_files, loadRData)
merged_imp <- mergeMice(imp)

merged_imp$loggedEvents
plot(merged_imp)
densityplot(merged_imp)

#load("/Users/jdamen/Library/CloudStorage/OneDrive-UMCUtrecht/Research/WHO ZIKA/2. Data/merged_imp.RData")
#data.imp<-complete(merged_imp,action="long")
#write.csv(data.imp,"/Users/jdamen/Library/CloudStorage/OneDrive-UMCUtrecht/Research/WHO ZIKA/2. Data/20230227 zikv_imputed.csv", row.names = FALSE)
