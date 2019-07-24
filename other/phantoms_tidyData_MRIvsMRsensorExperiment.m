%% Tidy Phantom Data from MRI vs MR Sensor Experiment
% This script is specifically made to tidy the MRI vs MR Sensor experiment
% data, where T2 measurements were collected on both phantoms and meat
% sample with both the MRI and MR sensor on the same day. 
%

clear

%% ---- USER INPUTS ----
doSave = 1; 
savePath = '/Volumes/CimaLabLargeFiles/Hydration/MRI_Study/ProcessedData/phantoms'; 
%% ---------------------

%% Initialize
master = table(); % initialize results table
whichRow = 0; 

%% MR Sensor - Phantoms
load('/Phantoms_MRIvsMR/Leg/Results/rawData_AM_phantoms_allIndividualTrials.mat')
load('/Phantoms_MRIvsMR/Leg/Results/rawData_AM_phantoms_allIndividualTrials_Time.mat')
data = dataMaster; 

fnSamples = fieldnames(data); 
nSamples = length(fnSamples); 
for iSample = 1:nSamples
    subStrings = strsplit(fnSamples{iSample},'_'); 
    
    tempData = data.(fnSamples{iSample}); 
    fnTrials = data.(fnSamples{iSample}).Properties.VariableNames; 
    nTrials = length(fnTrials); 
    for jTrial = 1:nTrials
        whichRow = whichRow +1; 
        
        master.Sample(whichRow) = subStrings(1); 
        master.Sensor(whichRow) = {'MR Sensor'}; 
        master.nDummies(whichRow) = subStrings(end); 
        master.Trial(whichRow) = fnTrials(jTrial); 
        master.T2Decay(whichRow) = {data.(fnSamples{iSample}).(fnTrials{jTrial})}; 
        master.Time(whichRow) = {time(:,iSample)};
        %master.Slice(whichRow) = NaN; 
    end
end
clear data

%% MR Sensor - Meat
load('/Volumes/CimaLabLargeFiles/lcolucci/2017-11-12_Phantoms_MRIvsMR/Leg/Results/rawData_PM_meat_allIndividualTrials.mat')
load('/Volumes/CimaLabLargeFiles/lcolucci/2017-11-12_Phantoms_MRIvsMR/Leg/Results/rawData_PM_meat_allIndividualTrials_Time.mat')
data = dataMaster; 

fnSamples = fieldnames(data); 
nSamples = length(fnSamples); 
for iSample = 1:nSamples
    subStrings = strsplit(fnSamples{iSample},'_'); 
    
    tempData = data.(fnSamples{iSample}); 
    fnTrials = data.(fnSamples{iSample}).Properties.VariableNames; 
    nTrials = length(fnTrials); 
    for jTrial = 1:nTrials
        whichRow = whichRow +1; 
        
        master.Sample(whichRow) = subStrings(1); 
        master.Sensor(whichRow) = {'MR Sensor'}; 
        master.nDummies(whichRow) = subStrings(end); 
        master.Trial(whichRow) = fnTrials(jTrial); 
        master.T2Decay(whichRow) = {data.(fnSamples{iSample}).(fnTrials{jTrial})}; 
        master.Time(whichRow) = {time(:,iSample)}; 
        %master.Slice(whichRow) = {NaN}; 
    end
end
clear data

%% MRI - Phantoms
load('/Volumes/CimaLabLargeFiles/lcolucci/2017-11-12_Phantoms_MRIvsMR/MRI/Results/mri_ROI_means_phantoms.mat')
data = mri_ROImeans; 

fnTrials = fieldnames(data); 
nTrials = length(fnTrials); 
for iTrial = 1:nTrials 
    
    fnTE = fieldnames(data.(fnTrials{iTrial})); %nDummies
    nTE = length(fnTE); 
    
    for jTE = 1:nTE
        
        fnSamples = fieldnames(data.(fnTrials{iTrial}).(fnTE{jTE}));
        nSamples = length(fnSamples); 
        for kSample=1:nSamples
            
            fnSlices = fieldnames(data.(fnTrials{iTrial}).(fnTE{jTE}).(fnSamples{kSample}));
            nSlices = length(fnSlices); 
            for lSlice = 1:nSlices
                whichRow = whichRow +1; 
                subStrings = strsplit(fnSamples{kSample},'_'); 
                
                master.Sample(whichRow) = subStrings(1); 
                master.Sensor(whichRow) = {'MRI'}; 
                master.nDummies(whichRow) = fnTE(jTE); 
                master.Trial(whichRow) = fnTrials(iTrial); 
                master.Slice(whichRow) = {fnSlices(lSlice)}; 
                master.T2Decay(whichRow) = {data.(fnTrials{iTrial}).(fnTE{jTE}).(fnSamples{kSample}).(fnSlices{lSlice})};
                if contains(fnTE{jTE},'shorterTE')
                    TE = 8; 
                    master.Time(whichRow) = {TE:TE:32*TE}; 
                elseif contains(fnTE{jTE},'longerTE')
                    TE = 25.5; 
                    master.Time(whichRow) = {TE:TE:32*TE};
                else
                    disp('No time array was set')
                end
            end
        end
    end
end
clear data

%% MRI - Meat
load('/Volumes/CimaLabLargeFiles/lcolucci/2017-11-12_Phantoms_MRIvsMR/MRI/Results/mri_ROI_means_meat.mat')
data = mri_ROImeans; 

fnSamples = fieldnames(data); 
nSamples = length(fnSamples); 
for iSample = 1:nSamples
    whichRow = whichRow +1; 
    subStrings = strsplit(fnSamples{iSample},'_');
    
    master.Sample(whichRow) = subStrings(1); 
    master.Sensor(whichRow) = {'MRI'}; 
    master.nDummies(whichRow) = {'shorterTE'}; 
    master.Trial(whichRow) = {'trial1'}; 
    master.Slice(whichRow) = subStrings(end); 
    master.T2Decay(whichRow) = {data.(fnSamples{iSample})};
    TE = 8; 
    master.Time(whichRow) = {TE:TE:32*TE}; 
end
clear data


%% Clean up the dataset

% Convert to categorical variables 
master.Sample = categorical(master.Sample); 
master.Sensor = categorical(master.Sensor); 
master.nDummies = categorical(master.nDummies); 
master.Trial = categorical(master.Trial); 
%master.Slice = categorical(master.Slice); 

% re-name 'oil' to 'oil Phantom'
master.Sample(master.Sample=='oil') = categorical({'oilPhantom'}); 
master.Sample(master.Sample=='porkbellyfat') = categorical({'fat'}); 


% average all time points together
uniqueRowCombos =unique(master(:,1:3), 'rows'); 
nUniqueRowCombos = height(uniqueRowCombos); 
for iRow = 1:nUniqueRowCombos
    whichRow = whichRow +1;
    
    tempData = master(master.Sample==uniqueRowCombos.Sample(iRow) & master.Sensor==uniqueRowCombos.Sensor(iRow) & master.nDummies==uniqueRowCombos.nDummies(iRow),:);
    if uniqueRowCombos.Sensor(iRow)=='MRI' && height(tempData) > 1
        tempData = tempData(strcmp([tempData.Slice{:}],'avg'),:); 
    end
    nMeas = height(tempData); 
    aggregateT2 = []; aggregateTime = []; 
    for j=1:nMeas
        
        if size(tempData.T2Decay{j},1) > size(tempData.T2Decay{j},2)
            aggregateT2 = horzcat(aggregateT2, tempData.T2Decay{j});
            aggregateTime = horzcat(aggregateTime, tempData.Time{j});
            longDim = 2;
        else
            aggregateT2 = vertcat(aggregateT2, tempData.T2Decay{j});
            aggregateTime = vertcat(aggregateTime, tempData.Time{j});
            longDim = 1; 
        end
    end
    avgT2 = mean(aggregateT2,longDim);
    avgTime = mean(aggregateTime,longDim);
    
    % Put results in master
    master.Sample(whichRow) = uniqueRowCombos.Sample(iRow); 
    master.Sensor(whichRow) = uniqueRowCombos.Sensor(iRow); 
    master.nDummies(whichRow)=uniqueRowCombos.nDummies(iRow); 
    master.T2Decay(whichRow) = {avgT2}; 
    master.Time(whichRow) = {avgTime}; 
    master.Trial(whichRow) = categorical({'Average'}); 
    
    clear aggregateT2 avgT2 aggregateTime avgTime
end
    
    
%% Save
if doSave==1
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end
    
    save(fullfile(savePath, 'phantoms_T2decays.mat'), 'master')
end


