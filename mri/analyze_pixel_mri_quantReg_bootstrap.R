# Fig 4 - MRI pixelwise, quantile regression 
# Data comes from script: process_data_for_R.m

# ------------ Quantile Regression with Clustering (boot.rq) ------------ 
# RUN: source('/Users/linacolucci/Documents/GitHub/MRI/mri/analyze_pixel_mri_quantReg_bootstrap.R')

# Prep
rm(list=ls())
library(quantreg)
library(dplyr)
library(tictoc)

# ------ USER INPUTS ---------
data = read.csv('/Paper/data_for_R/mri_pixelwise_muscle_RelAmp2_subset100.csv') 
quantileList <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9) #seq(0.05, 0.95, by=0.05)
nReps <- 1000 # number bootstrap replications
savePath <- '/Results3' #this folder and a 'tmp' sub-folder will be created if they do not already exist
# ----------------------------

# Suppress warnings
# oldw <- getOption("warn")
# options(warn = -1)
# options(warn = oldw) #turn on default warning, which is 0

# Check folder existence & create if necessary
if (!dir.exists(savePath)) {
  dir.create(savePath)
  cat('Result folder was created:',savePath,'\n')
}
tmpSavePath <- file.path(savePath,'tmp')
if (!dir.exists(tmpSavePath)) {
  dir.create(tmpSavePath)
  cat('Temp. sub-folder was created:',tmpSavePath,'\n')
}

# Prep data
data$AM[data$PM==1] <- 0 
data$AM[data$PM==0] <- 1 
data$HC[data$HD==1] <- 0 
data$HC[data$HD==0] <- 1 
y <- data$Values
clus <- data$SubjectNum

# Analysis
masterList <- vector('list',length(quantileList))
resultTable <- data.frame()
iList=1
jCounter = 1
tic("total")
for (q in quantileList) {
  tic("loop iteration")
  
  # run bootstrapped quantile reg.
  modelHDPM <- rq(y ~ PM * HD, tau = q, method = "fn", data)
  resultHDPM <- summary(modelHDPM, se = "boot", tau = q, R = nReps, bsmethod = "wild", cluster = clus)
  
  modelHCAM <- rq(y ~ AM * HC, tau = q, method = "fn", data)
  resultHCAM <- summary(modelHCAM, se = "boot", tau = q, R = nReps, bsmethod = "wild", cluster = clus)
  
  # store raw results in list
  tmpHDPM <- list(levels = 'HDPM', quantile = resultHDPM$tau, coefs = resultHDPM$coefficients, formula = resultHDPM$call, R = nReps)
  masterList[[iList]] <- tmpHDPM
  
  tmpHCAM <- list(levels = 'HCAM', quantile = resultHCAM$tau, coefs = resultHCAM$coefficients, formula = resultHCAM$call, R = nReps)
  masterList[[iList+1]] <- tmpHCAM
  
  # store important results in dataframe
  resultTable <- rbind(resultTable, data.frame(quantile = resultHDPM$tau, level = 'diff_PrePost_HC', value = resultHDPM$coefficients[2,1], p = resultHDPM$coefficients[2,4]))
  resultTable <- rbind(resultTable, data.frame(quantile = resultHDPM$tau, level = 'diff_HCHD_AM',    value = resultHDPM$coefficients[3,1], p = resultHDPM$coefficients[3,4]))
  resultTable <- rbind(resultTable, data.frame(quantile = resultHDPM$tau, level = 'interaction',     value = resultHDPM$coefficients[4,1], p = resultHDPM$coefficients[4,4]))
  
  resultTable <- rbind(resultTable, data.frame(quantile = resultHCAM$tau, level = 'diff_PrePost_HD', value = resultHCAM$coefficients[2,1], p = resultHCAM$coefficients[2,4]))
  resultTable <- rbind(resultTable, data.frame(quantile = resultHCAM$tau, level = 'diff_HCHD_PM',    value = resultHCAM$coefficients[3,1], p = resultHCAM$coefficients[3,4]))
  # resultTable <- rbind(resultTable, data.frame(quantile = resultHCAM$tau, level = 'interaction_HCAM',     value = resultHCAM$coefficients[4,1], p = resultHCAM$coefficients[4,4])) #same as interaction HDPM
  
  # Save temporary variables (in case something crashes)
  save(tmpHCAM, tmpHDPM, file=paste(tmpSavePath,'/tau',gsub('[.]','',q),'.RData',sep=''))
  write.csv(resultTable, file=paste(tmpSavePath,'/upToTau',gsub('[.]','',q),'.csv',sep=''))
  
  # Print progress
  cat("Finished loop", jCounter, "out of", length(quantileList), "( quantile", q,")", "\n")
  iList<- iList+2
  jCounter <- jCounter + 1
  toc()
}
print("Finished analysis!")
toc()

# Save results
save(masterList,resultTable, file=paste(savePath,'/results.RData',sep=''))
write.csv(resultTable, file=paste(savePath,'/resultSummary.csv',sep=''))
print("Results have been saved!")
print("----- DONE -----")


