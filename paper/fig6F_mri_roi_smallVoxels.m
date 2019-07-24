%% ------ Fig 6F, Tables S9-S12: MRI small voxels -------
%
%  
%  PREVIOUS SCRIPTS
%       data comes from script: analyze_ROI_multiexp_results.m 
%            NOTE: pVals are for a two-sample, two-tailed t-test with unequal vars (Welch test) 
%
%  OUTPUTS
%
%
%  Lina A. Colucci, 2018


clear; close all

%% ----- USER INPUTS ----
whichParameter = 'relative_amp2'; 
pathToMriFolder = '/Volumes/CimaLabLargeFiles/Hydration'; % path to Hydration folder in the drive (folder that contains 'MRI_Study' folder)
% ---
doSave = 1; % 1 if you want to save figure
savePath = '/Users/linacolucci/Documents/Thesis/Paper/figures/MRI_ROI';
%% ---------------------

%% Load Data
load(fullfile(pathToMriFolder,'MRI_Study/ProcessedData/mri/ROIs/ROIs_2expFits_tTestResults.mat'))
data = tTestResults; 
data.Parameter = categorical(data.Parameter); 

listROIs = {'voxel_anterior1','voxel_anterior2','voxel_lateral1','voxel_lateral2'}; 
nROIs = length(listROIs); 
for iROIs = 1:nROIs
    % Subset
    temp = data(data.ROI==listROIs(iROIs) & data.Parameter==whichParameter, :); 
    
    % Data
    tempData = temp.tTest{1,1}; 
    tempData = tempData(:,1:5); 
    
    % Save
    if doSave==1  
        if ~(7==exist(savePath)) % exist = 7 when it's a folder
            mkdir(savePath)
        end 
        writetable(tempData, fullfile(savePath, sprintf('fig8_MRI_%s_%s_%s_TableSummary.csv',listROIs{iROIs},whichParameter,temp.FitType{1})))
    end
end





