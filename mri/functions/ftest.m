%% F-test
%  This scripts performs an F-test to determine which multi-exponential
%  model is best for the data
%
%  INPUTS
%       time - vector
%       data - vector with T2 decay data (same length as time)
%       whichModels - structure with 0/1 for whether or not to perform a
%                     particular type of multi-exp model 
%       startVals - structure with arrays in each multi-exp model field that
%                   contain starting values for each multi-exp fit
%
%  OUTPUTS
%       modelComparisons - table with results of model comparison. Each row
%                          is a different head-to-head model comparison with 
%                          corresponding f-ratio, and p-value
%       fittingResults - structure with the fitting result of each
%                        multi-exp model. A table is contained within each
%                        field of the structure.
%
%  2018 LColucci

function [modelComparisons,fittingResults] = ftest(time,data,whichModels,startVals)
    
% Calculate sum of squares and degrees of freedom for each model 
[ss, df, fittingResults] = calculate_ss_df(time, data, whichModels, startVals); 

% Calculate f-test comparisons between nested models
modelComparisons = table();
modelComparisons.Model1 = {'exp1', 'exp2', 'exp1',       'exp2',       'exp3',       'exp1offset', 'exp2offset', 'exp3', 'exp4'}'; % model 1 is the simpler of the two being compared
modelComparisons.Model2 = {'exp2', 'exp3', 'exp1offset', 'exp2offset', 'exp3offset', 'exp2',       'exp3',       'exp4', 'exp5'}' ; 
for i = 1:height(modelComparisons)
    try 
        modelComparisons.f_ratio(i) =  ( (ss.(modelComparisons.Model1{i}) - ss.(modelComparisons.Model2{i})) / (df.(modelComparisons.Model1{i}) - df.(modelComparisons.Model2{i})) ) / ( ss.(modelComparisons.Model2{i}) / df.(modelComparisons.Model2{i}) ); 
        modelComparisons.p_val(i) = 1 - fcdf(modelComparisons.f_ratio(i), (df.(modelComparisons.Model1{i}) - df.(modelComparisons.Model2{i})), df.(modelComparisons.Model2{i})); 
    catch % if that model wasn't calculated, skip it
        modelComparisons.f_ratio(i) = NaN; 
        modelComparisons.p_val(i) = NaN; 
    end
end