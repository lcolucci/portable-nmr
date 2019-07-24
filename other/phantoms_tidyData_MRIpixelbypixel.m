%% Tidy phantom data for MRI pixel-by-pixel Results
%  This script loads in existing pixel-by-pixel MRI results for both meat
%  and phantom data and organizes it into a master table. The result of
%  this script is loaded into phantoms_visualize_MRIvsMRsensorExperiment.m
% 

doSave = 1; 
savePath = '/Volumes/CimaLabLargeFiles/Hydration/MRI_Study/ProcessedData/phantoms';

pixelbypixel = table(); 
whichRow = 0; 
%% Meat 1exp
meat_1exp = load('/Volumes/CimaLabLargeFiles/lcolucci/2017-11-12_Phantoms_MRIvsMR/MRI/Results/meat_1exp.mat'); 
data = meat_1exp.resultsMeat_1exp; 
fnSample = fieldnames(data); 
nSamples = length(fnSample); 
for iSample = 1:nSamples
    substrings = strsplit(fnSample{iSample},'_');
    fnSlices = fieldnames(data.(fnSample{iSample})); 
    nSlices = length(fnSlices); 
    for jSlice = 1:nSlices
        whichRow = whichRow+1; 
        pixelbypixel.Sample(whichRow) = substrings(1); 
        pixelbypixel.Sensor(whichRow) = {'MRI'}; 
        pixelbypixel.nDummies(whichRow) = {'shorterTE'}; 
        pixelbypixel.Slice(whichRow) = fnSlices(jSlice);
        pixelbypixel.Fit(whichRow) = {'1exp'}; 
        pixelbypixel.Data(whichRow) = {data.(fnSample{iSample}).(fnSlices{jSlice})}; 
    end
    
    clear substrings fnSlices nSlices 
end
clear data

%% Meat 2exp
meat_2exp = load('/Volumes/CimaLabLargeFiles/lcolucci/2017-11-12_Phantoms_MRIvsMR/MRI/Results/meat_2exp.mat');
data = meat_2exp.resultsMeat_2exp; 
fnSample = fieldnames(data); 
nSamples = length(fnSample); 
for iSample = 1:nSamples
    substrings = strsplit(fnSample{iSample},'_');
    fnSlices = fieldnames(data.(fnSample{iSample})); 
    nSlices = length(fnSlices); 
    for jSlice = 1:nSlices
        whichRow = whichRow+1; 
        pixelbypixel.Sample(whichRow) = substrings(1); 
        pixelbypixel.Sensor(whichRow) = {'MRI'}; 
        pixelbypixel.nDummies(whichRow) = {'shorterTE'}; 
        pixelbypixel.Slice(whichRow) = fnSlices(jSlice);
        pixelbypixel.Fit(whichRow) = {'2exp'}; 
        pixelbypixel.Data(whichRow) = {data.(fnSample{iSample}).(fnSlices{jSlice})}; 
    end
    clear substrings fnSlices nSlices 
end
clear data

%% Phantoms AM 1exp 
load('/Volumes/CimaLabLargeFiles/lcolucci/2017-11-12_Phantoms_MRIvsMR/MRI/Results/allResults_1exp_phantoms.mat')
data = allResults_1exp; 
fnSample = fieldnames(data); 
nSamples = length(fnSample); 
for iSample = 1:nSamples
    sampleSubstrings = strsplit(fnSample{iSample},'_');
    
    fnTimePts = fieldnames(data.(fnSample{iSample})); 
    nTimePts = length(fnTimePts); 
    for jTimePt = 1:nTimePts
        timePtSubstrings = strsplit(fnTimePts{jTimePt}, '_'); 
        
        fnSlices = fieldnames(data.(fnSample{iSample}).(fnTimePts{jTimePt}));
        nSlices = length(fnSlices); 
        for kSlice = 1:nSlices
            whichRow = whichRow+1; 
            pixelbypixel.Sample(whichRow) = sampleSubstrings(1); 
            pixelbypixel.TimePt(whichRow) = timePtSubstrings(1); 
            pixelbypixel.Sensor(whichRow) = {'MRI'}; 
            pixelbypixel.nDummies(whichRow) = timePtSubstrings(end); %{'shorterTE'}; 
            pixelbypixel.Slice(whichRow) = fnSlices(kSlice);
            pixelbypixel.Fit(whichRow) = {'1exp'}; 
            pixelbypixel.Data(whichRow) = {data.(fnSample{iSample}).(fnTimePts{jTimePt}).(fnSlices{kSlice})}; 
        end
        
        clear timePtSubstrings fnSlices nSlices
    end
    
    clear sampleSubstrings fnTimePts nTimePts
end
clear data  

%% 2exp: Oil and agarShort 
load('/Volumes/CimaLabLargeFiles/lcolucci/2017-11-12_Phantoms_MRIvsMR/MRI/Results/allResults_2exp_PM_oil_agarShort.mat')
data = allResults_2exp;
fnSample = fieldnames(data); 
nSamples = length(fnSample); 
for iSample = 1:nSamples
    sampleSubstrings = strsplit(fnSample{iSample},'_');
    
    fnTimePts = fieldnames(data.(fnSample{iSample})); 
    nTimePts = length(fnTimePts); 
    for jTimePt = 1:nTimePts
        timePtSubstrings = strsplit(fnTimePts{jTimePt}, '_'); 
        
        fnSlices = fieldnames(data.(fnSample{iSample}).(fnTimePts{jTimePt}));
        nSlices = length(fnSlices); 
        for kSlice = 1:nSlices
            whichRow = whichRow+1; 
            pixelbypixel.Sample(whichRow) = sampleSubstrings(1); 
            pixelbypixel.TimePt(whichRow) = timePtSubstrings(1); 
            pixelbypixel.Sensor(whichRow) = {'MRI'}; 
            pixelbypixel.nDummies(whichRow) = timePtSubstrings(end); %{'shorterTE'}; 
            pixelbypixel.Slice(whichRow) = fnSlices(kSlice);
            pixelbypixel.Fit(whichRow) = {'2exp'}; 
            pixelbypixel.Data(whichRow) = {data.(fnSample{iSample}).(fnTimePts{jTimePt}).(fnSlices{kSlice})}; 
        end
        
        clear timePtSubstrings fnSlices nSlices
    end
    
    clear sampleSubstrings fnTimePts nTimePts
end
clear data 

%% 2exp: agarLong and legSensorPhantom
load('/Volumes/CimaLabLargeFiles/lcolucci/2017-11-12_Phantoms_MRIvsMR/MRI/Results/allResults_2exp_PM_agarLong_legPhantom.mat')
data = allResults_2exp;
fnSample = fieldnames(data); 
nSamples = length(fnSample); 
for iSample = 1:nSamples
    sampleSubstrings = strsplit(fnSample{iSample},'_');
    
    fnTimePts = fieldnames(data.(fnSample{iSample})); 
    nTimePts = length(fnTimePts); 
    for jTimePt = 1:nTimePts
        timePtSubstrings = strsplit(fnTimePts{jTimePt}, '_'); 
        
        fnSlices = fieldnames(data.(fnSample{iSample}).(fnTimePts{jTimePt}));
        nSlices = length(fnSlices); 
        for kSlice = 1:nSlices
            whichRow = whichRow+1; 
            pixelbypixel.Sample(whichRow) = sampleSubstrings(1); 
            pixelbypixel.TimePt(whichRow) = timePtSubstrings(1); 
            pixelbypixel.Sensor(whichRow) = {'MRI'}; 
            pixelbypixel.nDummies(whichRow) = timePtSubstrings(end); %{'shorterTE'}; 
            pixelbypixel.Slice(whichRow) = fnSlices(kSlice);
            pixelbypixel.Fit(whichRow) = {'2exp'}; 
            pixelbypixel.Data(whichRow) = {data.(fnSample{iSample}).(fnTimePts{jTimePt}).(fnSlices{kSlice})}; 
        end
        
        clear timePtSubstrings fnSlices nSlices
    end
    
    clear sampleSubstrings fnTimePts nTimePts
end
clear data 

%% Save
if doSave==1
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end
    
    save(fullfile(savePath,'phantoms_mri_pixelbypixel.mat'), 'pixelbypixel')
end
