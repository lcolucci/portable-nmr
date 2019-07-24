%% process_ROImeans
%  This script loads in 'ROImeans_table.mat' and fits the decay of each ROI 
%  to 1- and 2-exp fits. 
%  
%  PREVIOUS SCRIPT
%       data is coming from: process_pixelData_to_ROImeans.m
%
%  INPUTS
%
%
%  OUTPUTS
%       ROIs_2expFits.mat and ROIs_1expFits.mat - structure containing a fit table for each ROI 
%           Table columns are: 
%           relax1 | relax1std | relax2 | relax2std | etc. | r2 | sse | etc. | startVals | Subject | AMPM | ROI | HD
%           
%           One row per subject/AMPM combination
%
%           NOTE: std = 95% CI of the fit (Matlab's cftool)
%
%  Lina A. Colucci, 2018


clear 
%% ---- USER INPUTS ----
time = [8:8:8*32]'; % Time array: [8:8:8*32]' or [25.5:25.5:25.5*32]';
ignore = 1;         % Skip 1st point or not? 0 = use all points, 1 = skip 1st point, etc. 
StartPoint = [1500 60 1000 170]; % For 2-exp fit, what are starting points? 
ROIlist = {''};     % Which ROI do you want to analyze? If empty, then all ROIs will be analyzed 
pathToMriStudy = '/Volumes/CimaLabLargeFiles/Hydration'; % Path to folder that contains 'MRI_Study' folder (contains the ROImeans data)
doSave = 0;         % 1 to save
savePath = '/Volumes/CimaLabLargeFiles/Hydration/MRI_Study/ProcessedData/mri/ROIs'; % Results will be saved in this folder

%% ---------------------

%% Load Data
load(fullfile(pathToMriStudy, 'MRI_Study/ProcessedData/mri/ROIs/ROImeans_table.mat'))
data = ROImeans_table; % <--- if name changes, need to adjust this

%% Time: ignore initial points 
tempTime = time(1+ignore:end); 

%% Loop through ROIs
if isempty(ROIlist) || isempty(ROIlist{:})
    ROIlist = unique(data.ROI); 
end
nROIs = length(ROIlist); 
for iROIs = 1:nROIs
     
    % Subset Data on 1 ROI 
    if iscategorical(ROIlist)
        whichROI = char(ROIlist(iROIs));
    elseif iscell(ROIlist)
        whichROI = ROIlist{iROIs}; 
    else
        error('''whichROI'' must be a cell or categorical variable')
    end
    dataSubset = data(data.AllSlices==1 & data.ROI == whichROI,:); %default is to process AllSlices==1
    
    % Looping through All Datasets
    nRows = height(dataSubset); 
    tempResult_2exp_All = table(); 
    tempResult_1exp_All = table(); 
    for jRows = 1:nRows

        % Select the right T2 Decay
        tempAmp = dataSubset.T2Decay{jRows}; 

        % Ignore initial points
        tempAmp = tempAmp(1+ignore:end); 

        % Do 1- and 2-exp Fittings
        tempResult_2exp = t2_fitting_2exp(tempTime, tempAmp, StartPoint);
        tempResult_1exp = t2_fitting_1exp(tempTime, tempAmp);

        % Calculate Relative Amp columns
        tempResult_2exp.relative_amp1 = 100.* tempResult_2exp.amp1 ./ nansum([tempResult_2exp.amp1, tempResult_2exp.amp2]); 
        tempResult_2exp.relative_amp2 = 100.* tempResult_2exp.amp2 ./ nansum([tempResult_2exp.amp1, tempResult_2exp.amp2]); 

        % Add Columns to Fitting Results
        tempResult_2exp.Subject = dataSubset.Subject(jRows); 
        tempResult_2exp.AMPM = dataSubset.AMPM(jRows); 
        tempResult_2exp.ROI = dataSubset.ROI(jRows); 
        
        tempResult_1exp.Subject = dataSubset.Subject(jRows); 
        tempResult_1exp.AMPM = dataSubset.AMPM(jRows); 
        tempResult_1exp.ROI = dataSubset.ROI(jRows); 
        
        % HD 
        if contains(cellstr(dataSubset.Subject(jRows)),'hd','IgnoreCase',true)
            tempResult_2exp.HD = 1; 
            tempResult_1exp.HD = 1; 
        else
            tempResult_2exp.HD = 0;
            tempResult_1exp.HD = 0; 
        end

        % Concatenate All Fittings in a table
        tempResult_2exp_All = vertcat(tempResult_2exp_All, tempResult_2exp) ; 
        tempResult_1exp_All = vertcat(tempResult_1exp_All, tempResult_1exp) ; 
        
        % Clear temporary variables
        clear tempResult_1exp tempResult_2exp tempAmp 
    end
    
    % Store results from this ROI in master structure 
    results_2exp.(whichROI) = tempResult_2exp_All;   
    results_1exp.(whichROI) = tempResult_1exp_All;  
    
    % Clear temporary variables
    clear whichROI dataSubset nRows tempResult_2exp_All tempResult_1exp_All
end

%% Save
if doSave==1
    save(fullfile(savePath, 'ROIs_2expFits.mat'), 'results_2exp')
    save(fullfile(savePath, 'ROIs_1expFits.mat'), 'results_1exp')
end

