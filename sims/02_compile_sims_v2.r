#Program: 02_compile_sims_v2.R
#Purpose: This program compiles sim results from each iteration (stored as separate files) into a single file

library(tidyverse)

## collect arguments passed in from SLURM job script
args <- commandArgs(trailingOnly=TRUE)
num.sims <- as.numeric(args[1])
subdir_name <- as.character(args[2])

# set working directory
user_home_dir <- "<location>" # top-level directory for the project
dir_path <- paste(user_home_dir,subdir_name,"/",sep="")
setwd(dir_path)


# get all results and put in one dataframe
list_of_all_results <- vector(mode = "list", length = num.sims)

for(i in 1:num.sims){
  skip_to_next<-FALSE
  # read in all results from results directory
  tryCatch(list_of_all_results[[i]] <- read.csv(file = paste("results_", i, ".csv", sep = "")
  )[,-1], error=function(e){skip_to_next <- TRUE})
  if(skip_to_next==TRUE) {
    list_of_all_results[[i]]<-matrix(rep(NA,6),nrow=1)
    skip_to_next<-FALSE}
    }
  
  all_results <- bind_rows(list_of_all_results)
  setwd(user_home_dir)
  write.csv(all_results, file = paste("results_", subdir_name, ".csv", sep = ""))
  




