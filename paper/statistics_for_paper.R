## Statistics for Paper
#  Data comes from Matlab script: process_data_for_R.m and analysis_bi_vs_weight.m (for BI)

# Prep
rm(list=ls())
library(dplyr)
library(perm)
library(EnvStats)
library(quantreg)

# ----- Demographics (Table 1) -----
data = read.csv('/Users/linacolucci/Documents/Thesis/Paper/data_for_R/demographics.csv')
data <- data[data$RecordID!='HDMRI-04',] #remove HDMRI04
whichCol <- 'Race'
HC <- data[data$HD==0,][[whichCol]]
HD <- data[data$HD==1,][[whichCol]]

# Welch Two Sample t-test
t.test(HC, HD, var.equal=FALSE, paired=FALSE) # use if 20+ samples each group. makes assumptions about underlying distribution

# *** Permutation test (it's not actually exact), p-vals calcualted from two-sample (or paired) perm test with X number of permutations
permTS(HC, HD, method='exact.mc', control=permControl(nmc=10^5-1)) # will get slightly diff p each time. Good test even if small sample size. Highest power out of non-parametric tests (Mann-Whitney, etc. makes more underlying assumptions, on ranks, etc.) 

# Categorical Var: Binomial test
#binom.test(HC,HD) #??? doesn't work. This test would be used to tell if, for example, gender distribution is different from 50% (random) or not
fisher.test(HC, HD) # Fisher's Exact Test
# Cat Var: Logistic Regression
#model <- glm(Race ~ . , family=binomial(link='logit'),data)
#summary(model)

# *** What about HC AM vs PM? e.g. WeightPre vs WeightPost (one-sample permutation test)
AM <- data$WeightPre[data$HD==1] # HD==0 for HC, HD==1 for HD
PM <- data$WeightPost[data$HD==1]
diff <- data$WeightChange[data$HD==1]
#t.test(AM, PM, var.equal=FALSE, paired=TRUE)
oneSamplePermutationTest(diff, n.permutations=1e5)


# ----- MRI: Pixel-wise Change (Fig. 2) -----
#   --- for a single ROI (experimenting) ---- 
# data = read.csv('/Users/linacolucci/Documents/Thesis/Paper/data_for_R/mri_pixelwise_allROIs_RelAmp2_change.csv') # just Relative Amp 2 data
# whichROI <- 'muscle_anterior'
# HC <- data[data$ROI==whichROI & data$HD==0,]$Change
# HD <- data[data$ROI==whichROI & data$HD==1,]$Change
# permTS(HC, HD, method='exact.mc', control=permControl(nmc=10^5-1)) 
#   --- loop through all ROIs (for paper) ----- 
data = read.csv('/Users/linacolucci/Documents/Thesis/Paper/data_for_R/mri_pixelwise_allROIs_allParameters_change.csv')
results <- data.frame()
paramList <- unique(data$Parameter)
nParams <- length(paramList)
roiList <- unique(data$ROI)
nROIs <- length(roiList)
for (i in 1:nParams) {
  whichParam <- paramList[i]
  for (j in 1:nROIs) {
    whichROI <- roiList[j]
    HC <- data[data$Parameter==whichParam & data$ROI==whichROI & data$HD==0,]$Change
    HD <- data[data$Parameter==whichParam &data$ROI==whichROI & data$HD==1,]$Change
    permResult <- permTS(HC, HD, method='exact.mc', control=permControl(nmc=10^5-1)) 
    results <- rbind(results, data.frame(ROI = whichROI, parameter = whichParam, pVal = permResult$p.value, meanDiffs = permResult$estimate))
  }
}
rownames(results) <- c()
write.csv(results, file = paste('/Users/linacolucci/Documents/Thesis/Paper/figures/MRI_pixelwise', '/mri_pixelwise_lastIntegralDiff_permTest_Rstudio_pVals.csv',sep = ''),row.names = F)
# 1-exp results (for bone only)
data = read.csv('/Users/linacolucci/Documents/Thesis/Paper/data_for_R/mri_pixelwise_1exp_allROIs_allParameters_change.csv')
whichROI <- 'bone'
HC <- data[data$ROI==whichROI & data$HD==0,]$Value
HD <- data[data$ROI==whichROI & data$HD==1,]$Value
permTS(HC, HD, method='exact.mc', control=permControl(nmc=10^5-1)) 

# ----- MRI: Pixel-wise Muscle AM vs PM (Fig 4B) -----
# Raw Data: all am/pm rel. amp 2 values within muscle
data = read.csv('/Users/linacolucci/Documents/Thesis/Paper/data_for_R/mri_pixelwise_muscle_RelAmp2.csv') 
# can't figure out how to load pixel values as a list within a dataframe column, going to just load separate columns instead
hc_AM = read.csv('/Users/linacolucci/Documents/Thesis/Paper/data_for_R/mri_pixelwise_muscle_RelAmp2_HC_AM.csv') 
hc_PM = read.csv('/Users/linacolucci/Documents/Thesis/Paper/data_for_R/mri_pixelwise_muscle_RelAmp2_HC_PM.csv') 
hd_AM = read.csv('/Users/linacolucci/Documents/Thesis/Paper/data_for_R/mri_pixelwise_muscle_RelAmp2_HD_AM.csv') 
hd_PM = read.csv('/Users/linacolucci/Documents/Thesis/Paper/data_for_R/mri_pixelwise_muscle_RelAmp2_HD_PM.csv') 
# ???? also ask about comparing HC AM vs PM, and HD AM vs PM
t.test(hc_AM, hc_PM, var.equal=FALSE, paired=FALSE)
# Quantile Regression with Clustering (boot.rq)
data = read.csv('/Users/linacolucci/Documents/Thesis/Paper/data_for_R/mri_pixelwise_muscle_RelAmp2_subset.csv') 
x <- model.matrix(~ data$PM * data$HD, data)
y <- data$Values
clus <- data$SubjectNum
boot.rq(x,y,tau=0.5, R=200, cluster=clus)

model1 <- rq(y ~ PM * HD, tau = 0.5, method = "fn", data)
summary(model1, se = "boot", tau = 0.5, R = 1000, bsmethod = "wild", cluster = clus)


# ------- MRI ROI Muscle All (Fig 5 B,C and Table 2)                                                                                                                                                                                                                                                # ----- MRI: ROI Average Muscle All (Table 2) -----
data = read.csv('/Users/linacolucci/Documents/Thesis/Paper/data_for_R/mri_roi_muscleAll.csv')
am_HC = data$AM[data$HD==0]
am_HD = data$AM[data$HD==1] 
pm_HC = data$PM[data$HD==0]
pm_HD = data$PM[data$HD==1]
diff_HC = data$Change[data$HD==0]
diff_HD = data$Change[data$HD==1]
# exact permutation
permTS(am_HC, am_HD, method='exact.mc', control=permControl(nmc=10^5-1))
permTS(pm_HC, pm_HD, method='exact.mc', control=permControl(nmc=10^5-1))
permTS(diff_HC, diff_HD, method='exact.mc', control=permControl(nmc=10^5-1))
# one-sample permutation (HC: AM vs PM, for example)
oneSamplePermutationTest(diff_HC, n.permutations=1e5)
oneSamplePermutationTest(diff_HD, n.permutations=1e5)


# ----- MR Sensor: R2 (Table 5) -----
data = read.csv('/Users/linacolucci/Documents/Thesis/Paper/data_for_R/mrSensor_R2.csv')
am_HC = data$AM[data$HD==0]
am_HD = data$AM[data$HD==1] 
pm_HC = data$PM[data$HD==0]
pm_HD = data$PM[data$HD==1]
diff_HC = data$Change[data$HD==0]
diff_HD = data$Change[data$HD==1]
# exact permutation 
permTS(am_HC, am_HD, method='exact.mc', control=permControl(nmc=10^5-1))
permTS(pm_HC, pm_HD, method='exact.mc', control=permControl(nmc=10^5-1))
permTS(diff_HC, diff_HD, method='exact.mc', control=permControl(nmc=10^5-1))
# one-sample permutation (HC: AM vs PM, for example)
oneSamplePermutationTest(diff_HC, n.permutations=1e5)
oneSamplePermutationTest(diff_HD, n.permutations=1e5)


# ------- MRI Small Voxels (Fig 8)
data = read.csv('/Users/linacolucci/Documents/Thesis/Paper/data_for_R/mri_rois_smallVoxels.csv')
# --- anterior 1 ----
whichROI <- 'voxel_lateral2'
temp <- filter(data, ROI == whichROI)
am_HC <- temp$AM[temp$HD==0]
am_HD <- temp$AM[temp$HD==1]
pm_HC <- temp$PM[temp$HD==0]
pm_HD <- temp$PM[temp$HD==1]
diff_HC <- temp$Change[temp$HD==0]
diff_HD <- temp$Change[temp$HD==1]
# exact permutation 
permTS(am_HC, am_HD, method='exact.mc', control=permControl(nmc=10^5-1))
permTS(pm_HC, pm_HD, method='exact.mc', control=permControl(nmc=10^5-1))
permTS(diff_HC, diff_HD, method='exact.mc', control=permControl(nmc=10^5-1))
# one-sample permutation (HC: AM vs PM, for example)
oneSamplePermutationTest(diff_HC, n.permutations=1e5)
oneSamplePermutationTest(diff_HD, n.permutations=1e5)


# ---- BI data ------
rm(list=ls())

data = read.csv('/Users/linacolucci/Documents/Paper/figures/Bioimpedance/BIdata.csv')
listType <- as.character(unique(data$Type)); 
listWholeBody <- c(0,1)
# --- 
resultsHCvsHD <- data.frame()
resultsAMvsPM <- data.frame()
for (i in 1:length(listWholeBody)){
  whichWholeBody <- listWholeBody[i]
  for (j in 1:length(listType)){
    whichType <- listType[j]
    
    # Prep Data for Stas
    subData <- filter(data, Type==whichType & WholeBody == whichWholeBody)
    if (all(is.nan(subData$Pre))) {next} #go on to next loop iteration if all values are NaN (i.e. for whole body TBW_calfLength)
    am_HC <- subData$Pre[subData$HD==0]
    am_HD <- subData$Pre[subData$HD==1]
    pm_HC <- subData$Post[subData$HD==0]
    pm_HD <- subData$Post[subData$HD==1]
    diff_HC <- subData$Change[subData$HD==0]
    diff_HD <- subData$Change[subData$HD==1]
    
    am_HC <- am_HC[!is.na(am_HC)]
    am_HD <- am_HD[!is.na(am_HD)]
    pm_HC <- pm_HC[!is.na(pm_HC)]
    pm_HD <- pm_HD[!is.na(pm_HD)]
    diff_HC <- diff_HC[!is.na(diff_HC)]
    diff_HD <- diff_HD[!is.na(diff_HD)]
    
    ## ---- HC vs HD ----
    # Welch Two Sample t-test
    welchAM <- t.test(am_HC, am_HD, var.equal=FALSE, paired=FALSE,na.action=na.omit) # use if 20+ samples each group. makes assumptions about underlying distribution
    welchPM <- t.test(pm_HC, pm_HD, var.equal=FALSE, paired=FALSE,na.action=na.omit) # use if 20+ samples each group. makes assumptions about underlying distribution
    welchDiff <- t.test(diff_HC, diff_HD, var.equal=FALSE, paired=FALSE,na.action=na.omit) # use if 20+ samples each group. makes assumptions about underlying distribution
    tempWelchRow <- data.frame(Type=whichType, Test = 'Welch_TwoSample', WholeBody=whichWholeBody, Pre_HCvsHD=welchAM$p.value, Post_HCvsHD=welchPM$p.value, Change_HCvsHD=welchDiff$p.value)
    # Permutation test (it's not actually exact), p-vals calcualted from two-sample (or paired) perm test with X number of permutations
    permAM <- permTS(am_HC, am_HD, method='exact.mc', control=permControl(nmc=10^5-1),na.action=na.omit) # will get slightly diff p each time. Good test even if small sample size. Highest power out of non-parametric tests (Mann-Whitney, etc. makes more underlying assumptions, on ranks, etc.) 
    permPM <- permTS(pm_HC, pm_HD, method='exact.mc', control=permControl(nmc=10^5-1),na.action=na.omit)
    permDiff <- permTS(diff_HC, diff_HD, method='exact.mc', control=permControl(nmc=10^5-1),na.action=na.omit)
    tempPermRow <- data.frame(Type=whichType, Test = 'Permutation', WholeBody=whichWholeBody, Pre_HCvsHD=permAM$p.value, Post_HCvsHD=permPM$p.value, Change_HCvsHD=permDiff$p.value)
    # Join
    resultsHCvsHD <- rbind(resultsHCvsHD, tempPermRow)
    resultsHCvsHD <- rbind(resultsHCvsHD, tempWelchRow)
    
    ## ---- AM vs PM ---- 
    # Welch paired t-test
    welchHC <- t.test(am_HC, pm_HC, var.equal=FALSE, paired=TRUE)
    welchHD <- t.test(am_HD, pm_HD, var.equal=FALSE, paired=TRUE)
    tempWelch1SampleRow <- data.frame(Type=whichType, Test = 'Welch_Paired', WholeBody=whichWholeBody, HC_AMvsPM=welchHC$p.value, HD_AMvsPM=welchHD$p.value)
    # One-sample permutation test
    permHC <- oneSamplePermutationTest(diff_HC, n.permutations=1e5)
    permHD <- oneSamplePermutationTest(diff_HD, n.permutations=1e5)
    tempPerm1SampleRow <- data.frame(Type=whichType, Test='Permutation',WholeBody=whichWholeBody, HC_AMvsPM=permHC$p.value, HD_AMvsPM=permHD$p.value)
    # Join results
    resultsAMvsPM <- rbind(resultsAMvsPM, tempPerm1SampleRow)
    resultsAMvsPM <- rbind(resultsAMvsPM, tempWelch1SampleRow)

    # Delte temp vars
    rm(subData, whichType,tempWelchRow,tempPermRow,tempWelch1SampleRow, tempPerm1SampleRow)
    
  }
}
# Save
savePath = '/Users/linacolucci/Documents/Paper/figures/Bioimpedance'
write.csv(resultsHCvsHD, file=paste(savePath,'/BI_stats_HCvsHD.csv',sep=''))
write.csv(resultsAMvsPM, file=paste(savePath,'/BI_stats_AMvsPM.csv',sep=''))






