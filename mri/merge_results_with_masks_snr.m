%%  Merge_Results_with_Masks_SNR
%   This script addts to the existing pixel-by-pixel fit results. 
%   This script can (1) merge the fit results with the ROI masks so we know 
%   which mask each pixel is in, (2) calculate the SNR or each pixel by
%   loading in the 'noiseSummary' that can the noise level for each scan. 
%
%   USER INPUTS
%       whichResults - which pixel-by-pixel results do you want to modify? '1exp' or '2exp'
%       pathToMriStudy - Path to folder that contains 'MRI_Study' folder (this is needed to load 1-exp and 2-exp results, and mask files) 
%       doSNR - do you want to add a column for pixel-by-pixel SNR? 1 if yes
%       pathToNoiseSummary - path to noiseSummary.mat variable 
%       doMaskMerge - do you want to add a column for each ROI showing whether this pixel appears in that ROI? 1 if yes
%       pathToWorkbook - Path to config.xlsx (contains paths to load mask files) 
%       sheetName - Which sheet name in config.xlsx to load?
%       doSave - do you want to save these updates results? 1 if yes
%
%   OUTPUTS
%       if doSave=1, updated 'results_1exp' or 'results_2exp' will overwrite the original files
%       if doSNR=1, same result but with an additional column for 'SNR' 
%                   SNR = maxVal / noise 
%                   noise = the Gaussian-corrected noise calculated in process_noise.m script
%       if doMaskMerge=1 same result but with additinoal columns for each Mask that existed for that patient. 
%                   The columns are named after the Mask name and are populated with 1 or NaN to indicate if that pixel is located in that Mask
%                   If a column with maskName already exists, that column is deleted first and then a the new one is added

clear 
%% --- USER INPUTS ---
% --- 
whichResults = 'pixelData'; % '1exp' or '2exp' or 'pixelData', load pixelDecays, 1-, or 2-exp pixel-by-pixel results
pathToMriStudy = '/Volumes/CimaLabLargeFiles/Hydration'; % Path to folder that contains 'MRI_Study' folder (contains 1-exp and 2-exp results, and mask files) - '/Volumes/CimaLabLargeFiles/Hydration'
% ---
doSNR = 0; % Calculate SNR for each pixel? 1 for yes. (only applies to '1exp' and '2exp', not 'pixelData')
pathToNoiseSummary = '/Volumes/CimaLabLargeFiles/Hydration/MRI_Study/ProcessedData/mri/mri_noise_summary.mat';
% ---
doMaskMerge = 1; % Calculate SNR for each pixel? 1 for yes
pathToWorkbook = '/GitHub/MRI/mri/config.xlsx'; % Path to config.xlsx (contains paths to load mask files) 
sheetName = 'MRI_Study_T2';                             % Which sheet name in config.xlsx to load?
% ---
doSave = 1; %1 for yes. NOTE: will overwrite and replace the original results (mri_pixelByPixel_1exp.mat, or mri_pixelByPixel_2exp.mat)
%% ---- END USER INPUTS ----

%% Load Data 
% Load config.xlsx file
config = load_config_file(pathToWorkbook, sheetName); 

% Load noise summary
if doSNR==1, load(pathToNoiseSummary), end

% Load pixel-by-pixel fit results
if strcmpi(whichResults, '1exp')
    pathToData = 'MRI_Study/ProcessedData/mri/pixel_by_pixel/mri_pixelByPixel_1exp.mat'; %if folder structure/names change, then change this path
    load(fullfile(pathToMriStudy, pathToData)) 
    resultName = who('result*'); result = eval(resultName{1}); 
elseif strcmpi(whichResults, '2exp')
    pathToData = 'MRI_Study/ProcessedData/mri/pixel_by_pixel/mri_pixelByPixel_2exp.mat'; 
    load(fullfile(pathToMriStudy, pathToData))
    resultName = who('result*'); result = eval(resultName{1}); 
elseif strcmpi(whichResults, 'pixelData')
    pathToData = 'MRI_Study/ProcessedData/mri/pixel_by_pixel/mri_t2PixelDecays.mat'; 
    load(fullfile(pathToMriStudy, pathToData))
    result = pixelT2Decays; 
end

%% Analysis
fnSubjects = fieldnames(result); % fn=fieldname
nSubjects = length(fnSubjects); 
whichRowSize = 0; 
for iSubjects=1:nSubjects
    whichSubject = fnSubjects{iSubjects}; 
    fnAMPM = fieldnames(result.(whichSubject));
    
    %% Loop through AM/PM
    for jAmPm=1:length(fnAMPM)
        whichAmPm = fnAMPM{jAmPm}; 
        whichConfigRow = find(config.subject==whichSubject & config.AMPM==whichAmPm); 
        
        % Directory of masks
        pathToMaskFolder = fullfile(pathToMriStudy, config.pathMasks{whichConfigRow}); 
        pathToMaskFolderSearchNiftis = fullfile(pathToMaskFolder,'*.nii');
        maskList = dir(pathToMaskFolderSearchNiftis); 
        
        % Find noise mask and delete it from list
        noiseRow = find(contains({maskList.name},'noise','IgnoreCase',true)); 
        maskList(noiseRow,:) = [];         
        
        %% Loop through all masks 
        nMasks = length(maskList);
        for kMasks = 1:nMasks
            % Load mask
            mask = load_nifti(fullfile(pathToMaskFolder,maskList(kMasks).name)); 
            mask = mask.vol; 
            
            % Define maskName
            maskName = maskList(kMasks).name; 
            maskName = strrep(maskName,'.nii',''); 
            
            %% Loop through slices
            fnSlices = fieldnames(result.(fnSubjects{iSubjects}).(fnAMPM{jAmPm})); 
            nSlices = size(mask,3); 
            for lSlices = 1:nSlices 
                sliceName = strcat('slice',num2str(lSlices)); 
                if ~strcmpi(sliceName, fnSlices{lSlices})
                    error('There is something wrong with the slice names/numbering. Look into it. lSlice number does not match up with fnSlices{lSlices}.')
                end
                tempMask = mask(:,:,lSlices); % select 1 slice
                lookupTable = table(); % create fresh lookup table for each slice
                [lookupTable.row, lookupTable.col] = find(tempMask ==1); % find row/col value of where mask=1
                lookupTable.(maskName)(:) = 1;             
                               
                %% Calculate SNR
                if doSNR==1
                    % Find corresponding noise level for this Subject/AMPM
                    whichRowNoiseSummary = find(noiseSummary.subject == fnSubjects{iSubjects} & noiseSummary.AMPM == fnAMPM{jAmPm}); 
                    noiseLevel = noiseSummary.noise(whichRowNoiseSummary); 
                    % Calculate SNR
                    result.(fnSubjects{iSubjects}).(fnAMPM{jAmPm}).(sliceName).SNR = result.(fnSubjects{iSubjects}).(fnAMPM{jAmPm}).(sliceName).maxVal ./ noiseLevel; 
                end
                                
                %% Merge Mask and Result
                if doMaskMerge==1
                   % If Column with this maskName already exists, delete it first
                   colNames = result.(fnSubjects{iSubjects}).(fnAMPM{jAmPm}).(sliceName).Properties.VariableNames; 
                   if any(strcmpi(colNames, maskName))
                       result.(fnSubjects{iSubjects}).(fnAMPM{jAmPm}).(sliceName).(maskName) = []; % delete existing column
                   end
                   
                   % Join data with lookup table to identify which pixels belong to this mask
                   result.(fnSubjects{iSubjects}).(fnAMPM{jAmPm}).(sliceName) = outerjoin(result.(fnSubjects{iSubjects}).(fnAMPM{jAmPm}).(sliceName), lookupTable,'MergeKeys',1); 
                end
                                
                clear whichRowNoiseSummary noiseLevel
            end
        end
    end     
end

%% Save
if doSave==1
    if strcmpi(whichResults, '1exp')
        result_1exp = result; 
        save(fullfile(pathToMriStudy, pathToData), 'result_1exp') %will overwrite the existing data
    end
    if strcmpi(whichResults, '2exp')
        result_2exp = result; 
        save(fullfile(pathToMriStudy, pathToData), 'result_2exp') %will overwrite the existing data
    end
    if strcmpi(whichResults, 'pixelData')
        pixelT2Decays = result; 
        save(fullfile(pathToMriStudy, pathToData), 'pixelT2Decays') %will overwrite the existing data
    end
end
    
