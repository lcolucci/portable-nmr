%% Aggregate QuantReg Results - Supplemental Table (RCE Cluster Analysis) 

% Load Data
resultsA = importQuantRegResults('/Paper/RCE_results/results/resultSummary.csv'); %0.1-0.7
resultsB = importQuantRegResults('/Paper/RCE_results/results2/tmp/upToTau09.csv'); %0.9
resultsC = importQuantRegResults('/Paper/RCE_results/results3/resultSummary.csv'); %0.6
resultsD = importQuantRegResults('/Paper/RCE_results/results4/resultSummary.csv'); %0.7-0.8

% Aggregate Data
results = vertcat(resultsA, resultsB, resultsD(resultsD.quantile==0.800001,:))


%% Process Results
% Select just quantiles that have interaction < 0.05 (happens to be all of them)
whichQuantiles = results(results.level=='interaction' & results.p < 0.05,:).quantile; 
resultsSubset = results(ismember(results.quantile,whichQuantiles),:); 

% Difference btw HC and HD AM
HCHD_AM = resultsSubset(resultsSubset.level == 'diff_HCHD_AM',:); 

% Difference btw HC and HD PM
HCHD_PM = resultsSubset(resultsSubset.level == 'diff_HCHD_PM',:); 

% Difference btw Pre and Post HC
PrePost_HC = resultsSubset(resultsSubset.level == 'diff_PrePost_HC',:); 

% Difference btw Pre and Post HD
PrePost_HD = resultsSubset(resultsSubset.level == 'diff_PrePost_HD',:); 


