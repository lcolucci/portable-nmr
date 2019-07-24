%% Determine optimal model to fit to data
%  This script runs 6 types of fits: 1-exp, 1-exp with offset, 2-exp, 2-exp with offset, 3-exp, 3-exp with offset
%  Applies the F-test to determine which is the optimal model for this particular data
%
%   INPUTS
%       time - define time array
%       data - define data array (same size as time)
%       ignore - number of initial points to skip (if any)
%       whichModels - put a 1 next to the type of models you want to run on this dataset
%       startVals - (optional) user-defined starting values for each type of fit
%
%   OUTPUTS
%       modelComparisons - table that runs head-to-head comparison of the models and determines which one is a better fit.
%                          p-value < 0.05 means the more complex model (Model2) should be chosen.
%       fittingResults - if you want to explore the actual fits that were given to each model,  


clearvars -except pixelT2Decays result_1exp data_table data
%% ------- USER INPUTS ----------
%data = data_table;
%data = data(data.AMPM=='am' | data.AMPM=='pm',:); 
%x = nmr.Time{1,1}; %data.Time{1,1}; %8:8:8*32; 
%y = nmr.T2Decay{1,1}; %data.T2Decay{1,1}; %pixelT2Decays.HDMRI01b.PM.slice1.decay{4092,1} ; %pixelT2Decays.HC01.AM.slice1.decay{183,1}; 
ignore = 1; % 0 = doesn't skip any points. 1 = skips 1st point and starts on 2nd point. 

% Which models do you want to compare? 
whichModels.do1exp = 1; % 1 or something else
whichModels.do1expoffset = 0;
whichModels.do2exp = 1;
whichModels.do2expoffset = 0;
whichModels.do3exp = 1;
whichModels.do3expoffset = 0;
whichModels.do4exp = 1; 
whichModels.do5exp = 1; 

% Which starting values do you want to specify for each fit (if any)? 
startVals.exp5 = [2 2 0.5 20 8 40 10 100 2 300];
startVals.exp4 = [2 2 8 35 10 100 2 300]; 
startVals.exp3offset = [4 20 15 70 4 275 5]; %[1000 20 1000 100 500 400 5]; % optional. empty square bracket [] or enter in 7 starting vals. 
startVals.exp3 = [2 2 8 35 10 100]; %[1000 20 1000 100 500 400]; % [] or 6 vals 
startVals.exp2offset = [600 60 1600 250 5]; % [] or 5 vals
startVals.exp2 = [100 40 1000 250]; % [] or 4 vals
startVals.exp1offset = []; % [] or 3 vals
startVals.exp1 = []; % [] or 2 vals
%% ------------------------------
master = table(); 
whichRow = 0; 
subjectList = unique(data.Subject); 
nSubjects = length(subjectList); 
for iSubject=1:nSubjects
    tempData = data(data.Subject == subjectList(iSubject),:); 
    AMPMlist = unique(tempData.AMPM); 
    nAMPM = length(AMPMlist); 
    for jAMPM = 1:nAMPM
        whichRow = whichRow + 1; 
        
        % Prep Data
        x = tempData(tempData.AMPM == AMPMlist(jAMPM),:).Time{:}; 
        y = tempData(tempData.AMPM == AMPMlist(jAMPM),:).T2Decay{:}; 
        
        % Ignore initial points
        x = x(ignore+1:end); 
        y = y(ignore+1:end); 
        
        % Calculate sum of squares and degrees of freedom for each model 
        [ss, df, fittingResults] = calculate_ss_df(x, y, whichModels, startVals); 

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
        
        % Aggregate results in master
        master.Subject(whichRow) = subjectList(iSubject); 
        master.AMPM(whichRow) = AMPMlist(jAMPM); 
        % 1- vs 2-exp
        master.F_ratio_1vs2(whichRow) = modelComparisons.f_ratio(1); 
        master.p_val_1vs2(whichRow) = modelComparisons.p_val(1); 
        % 2- vs 3-exp
        master.F_ratio_2vs3(whichRow) = modelComparisons.f_ratio(2); 
        master.p_val_2vs3(whichRow) = modelComparisons.p_val(2); 
        % 3- vs 4-exp
        master.F_ratio_3vs4(whichRow) = modelComparisons.f_ratio(8); 
        master.p_val_3vs4(whichRow) = modelComparisons.p_val(8); 
        % 4- vs 5-exp
        master.F_ratio_4vs5(whichRow) = modelComparisons.f_ratio(9); 
        master.p_val_4vs5(whichRow) = modelComparisons.p_val(9);
        
        % Delete temp vars
        clear modelComparisons fittingResults ss df x y 
    end
end

% Prep data
x = x(ignore+1:end); 
y = y(ignore+1:end); 

% Calculate sum of squares and degrees of freedom for each model 
[ss, df, fittingResults] = calculate_ss_df(x, y, whichModels, startVals); 

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

