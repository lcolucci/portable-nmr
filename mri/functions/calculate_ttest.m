%% calculate_ttest
%  This function takes a table of 1-exp or 2-exp fit results for a single ROI and 
%  summarizes it into a table of AM, PM, AM-PM Change results with a
%  p-value comparing HC vs HD values. 
%
%  This function is called in the script: "analyze_ROI_multiexp_results.m"
%       
%
%   Change - positive value = increase from AM to PM
%            negative value = decrease from AM to PM
%
% 2018-07-24 Implemented ttest2 with unequal variances for everything (i.e.
%            a Welch test instead of Student t-test)


function master = calculate_ttest(data, whichParamsList, doGOF)

%% Set Checks and Defaults
if length(unique(data.Subject))==1 
    master = []; 
    return
end
if length(unique(data.HD))==1
    master = []; 
    return
end 
if ~iscell(whichParamsList) || isempty(whichParamsList) || isempty(whichParamsList{1})
    error('''whichParamsList'' must be a cell array of strings')
end
if ~exist('whichParamsList') && any(strcmpi(dataSubset1.Properties.VariableNames,'relax2'))
    whichParamsList = {'relax1','relax2','relative_amp2'}; 
end
if ~exist('whichParamsList') && ~any(strcmpi(dataSubset1.Properties.VariableNames,'relax2'))
    whichParamsList = {'relax1'}; 
end
if ~exist('doGOF')
    doGOF = 1; % default is to include the rmse/r2 columns 
end
nParams = length(whichParamsList);

%% 'Change' table: Calculate AM-PM Differences
change = varfun(@(x) x(2)-x(1), data,'InputVariables',whichParamsList, 'GroupingVariables','Subject');
change.GroupCount = []; % Delete 'GroupCount' column 
try
    change.Properties.VariableNames = {'Subject',whichParamsList{:}}; 
catch
    change.Properties.VariableNames = {'subject',whichParamsList{:}}; 
end
change.AMPM(:) = categorical({'Change'}); 
ids_HD = contains(cellstr(change.Subject), 'HD'); 
change.HD(ids_HD)=1; % add 'HD' column 

%% Vertically stack 'data' and 'change'
% Find column IDs of the data columns I want to keep  
try colSubject = find(strcmpi(data.Properties.VariableNames, 'Subject')); catch colSubject = find(strcmpi(data.Properties.VariableNames, 'subject')); end
colAMPM = find(strcmpi(data.Properties.VariableNames, 'AMPM')); 
colHD = find(strcmpi(data.Properties.VariableNames, 'HD')); 
for ii=1:nParams
    colParameters(ii) = find(strcmpi(data.Properties.VariableNames, whichParamsList{ii})); 
end
% Vertically stack tables
dataJoined = vertcat(data(:,[colSubject, colAMPM, colHD, colParameters]), change);

%% Loop through each parameter and calculate HC vs HD p-value
master = table(); % initialize
for iParam=1:nParams
    
    % Which parameter to analyze? 
    whichParam = whichParamsList{iParam}; 
    whichParamColID = find(strcmpi(dataJoined.Properties.VariableNames, whichParam)); 
    
    % Table with just 1 Variable
    dataSubset = dataJoined(:,[1:3,whichParamColID]); % long-format table
    tTestTable = unstack(dataSubset, whichParam,'AMPM'); %convert to wide format
    
    % Add r2 and rmse values
    if doGOF==1
        tTestTable.rmse_AM = data.rmse(data.AMPM=='AM');
        tTestTable.rmse_PM = data.rmse(data.AMPM=='PM');
        tTestTable.r2_AM = data.r2(data.AMPM=='AM');
        tTestTable.r2_PM = data.r2(data.AMPM=='PM');
    end
    
    % Add ROI and whichParam columns 
    tTestTable.ROI(:) = data.ROI(1);
    tTestTable.Parameter(:) = {whichParam}; 
    
    % Calculate p-values for HC vs HD
    % Welch Test: two-sample, two-sided t-test with unequal variances
    [h, p_AM] =     ttest2( tTestTable.AM(tTestTable.HD==0),     tTestTable.AM(tTestTable.HD==1)    , 'VarType','unequal');
    [h, p_PM] =     ttest2( tTestTable.PM(tTestTable.HD==0),     tTestTable.PM(tTestTable.HD==1)    , 'VarType','unequal');
    [h, p_Change] = ttest2( tTestTable.Change(tTestTable.HD==0), tTestTable.Change(tTestTable.HD==1), 'VarType','unequal');
    
    % Add p-value results to end of table
    ttestRow = height(tTestTable)+1; 
    tTestTable.Subject(ttestRow) = {'HC vs HD, t-test p-value'}; 
    tTestTable.HD(ttestRow) = NaN; 
    tTestTable.AM(ttestRow) = p_AM; 
    tTestTable.PM(ttestRow) = p_PM; 
    tTestTable.Change(ttestRow) = p_Change; 
    if doGOF==1
        tTestTable.rmse_AM(ttestRow) = NaN; 
        tTestTable.rmse_PM(ttestRow) = NaN; 
        tTestTable.r2_AM(ttestRow) = NaN; 
        tTestTable.r2_PM(ttestRow) = NaN; 
    end
%     tTestTable.ROI(ttestRow) = NaN; 
    tTestTable.Parameter(ttestRow) = {' - '}; 
    
    % Store t-test results for this parameter in master table
    master.ROI(iParam) = data.ROI(1);
    master.Parameter(iParam) = {whichParam}; 
    master.tTest(iParam) = {tTestTable}; 
    
    % Clear temporary variables
    clear tTestTable dataSubset whichParam whichParamColID p* h
end


