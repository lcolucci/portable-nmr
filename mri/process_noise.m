%% Process Noise
%  This script calculates the noise level for each MRI scan and puts the
%  results in a summary table. 
%
%   INPUTS
%       whichIdToAnalyzeLoop - vector of patients to analyze. Take values from 'id' column of config.xlsx (sheet=MRI_Study_T2)
%       whichImageToAnalyzeLoop - cell array of strings
%       pathToMriStudy - string path to folder that contains the 'MRI_Study' folder
%       pathToWorkBook - string path to config.xlsx file
%       sheetName - which sheet of the config.xlsx file to look at?
%       doSave - do you want to save results? 1 if so. 
%       savePath - string. Results will be saved in this folder. 
%
%   OUTPUTS
%       If doSave==1, 'noiseMaster' and 'noiseSummary' will be saved tosavePath
%
%       noiseMaster 
%               Columns: subject | AMPM | slice | timept | noiseVector | mean | median | std
%               noiseVector = vector with all noise pixel values for that slice/timept
%               mean/median/std = the calculated statistics on noiseVector
%
%       noiseSummary
%               Columns: subject | AMPM | noiseRayleigh | noise
%               noiseRayleigh = the stdev of the noiseVector across all slices/timePts
%               noise         = the noise corrected for the Rayleigh distribution in the background of any magnitude MRI image
%                               (aka Gaussian noise)
%                               noise = noiseRayleigh / sqrt(2 - pi/2) 
%                               REF: https://web.stanford.edu/class/rad229/Notes/B4-SNR.pdf
%                                    Nishimura book... and many others


clear; close all

%% --- USER INPUTS ---
whichIdToAnalyzeLoop = [29:30];                          % Look at config.xlsx to find ID of the subjects you want to analyze
whichImageToAnalyzeLoop = {'legShorterTE'};              % (cell array) Which column of config.xlsx to analyze? e.g. {'legShorterTE'}, {'fingerLongerTE','additionalScan01'},etc.
% ---
pathToMriStudy = '/Users/linacolucci/Documents/'; % '/Volumes/CimaLabLargeFiles/Hydration'; % Path to folder that contains 'MRI_Study' folder
pathToWorkbook = '/Users/linacolucci/Documents/GitHub/MRI/mri/config.xlsx'; % Path to config.xlsx 
sheetName = 'MRI_Study_T2';                             % Which sheet name in config.xlsx to load?
% ---
doSave = 0; % 1 to save
savePath = '/Volumes/CimaLabLargeFiles/Hydration/MRI_Study/ProcessedData/mri'; 
%% --- END OF USER INPUTS ---

%% Calculate detailed table of noise values
config = load_config_file(pathToWorkbook, sheetName); % Load config excel file
whichMaskToProcess = 'noise'; % define noise mask
whichRow = 0; noiseMaster = table(); 
nIDsInLoop = length(whichIdToAnalyzeLoop); % an ID points to a Subject+AMPM combination
for iIDs = 1:nIDsInLoop
    try 
    % Pick the right variables from loop 
    whichIdToAnalyze = whichIdToAnalyzeLoop(iIDs); 
    whichConfigRow = find(config.id==whichIdToAnalyze);    
    [whichImageToAnalyze] = checkinputs_static_or_looping(iIDs, whichImageToAnalyzeLoop); %  Go from 'inputLoop' to 'input'  
%   fprintf('\nCurrently analyzing patientID #%d (%s %s), which is #%d out of %d patients to analyze. \n',config.id(whichConfigRow), config.subject{whichConfigRow}, config.AMPM{whichConfigRow}, iIDs, nIDsInLoop)

    % Load image (nifti)
    pathToImageNifti = fullfile(pathToMriStudy, config.pathNiftis{whichConfigRow}, config.(whichImageToAnalyze){whichConfigRow});
    image = load_nifti(pathToImageNifti); image = image.vol; 

    % Load mask (nifti)
    pathToMask = fullfile(pathToMriStudy, config.pathMasks{whichConfigRow},strcat(whichMaskToProcess,'.nii'));
    mask = load_nifti(pathToMask); mask = mask.vol; 
    mask = logical(mask); 

    %% Loop through slices
    nSlices = size(image, 3); 
    for jSlice = 1:nSlices
%         fprintf('Processing Slice %d of %d (%s) \n',jSlice,nSlices,datestr(now))
        
        %% Loop through TimePts
        whichRowTimePts = 0; nTimePts = size(image,4); 
        for kTimePt=1:nTimePts 
            % Go to next row with each loop iteration
            whichRow = whichRow+1;

            % Select Image/Noise for 1 slice, 1 timept
            imageSlice = image(:,:,jSlice,kTimePt); % 1 slice, 1 timePts
            noiseVector = imageSlice(squeeze(mask(:,:,jSlice))); % apply mask to image
           
            % Store results in master table
            noiseMaster.subject(whichRow) = (config.subject(whichConfigRow)); 
            noiseMaster.AMPM(whichRow) = (config.AMPM(whichConfigRow)); 
            noiseMaster.slice(whichRow) = jSlice; 
            noiseMaster.timept(whichRow) = kTimePt; 
            noiseMaster.noiseVector(whichRow) = {noiseVector}; 
            noiseMaster.mean(whichRow) = nanmean(noiseVector); 
            noiseMaster.median(whichRow) = nanmedian(noiseVector); 
            noiseMaster.std(whichRow) = nanstd(noiseVector); 

            clear imageSlice noiseVector 
        end
        
    end  
    catch
       fprintf('\n***NOTE***\nPatientID #%d (%s %s) was skipped. \n',config.id(whichConfigRow), config.subject{whichConfigRow}, config.AMPM{whichConfigRow})
    end
end
% Sort rows
noiseMaster = sortrows(noiseMaster, [1,2]); 
% <<<<<<<<<<<< not sure if there are pros/cons to now converting subject/AMPM columns to categorical type??

% Save
if doSave==1
    save(fullfile(savePath, 'mri_noise_all.mat'),'noiseMaster')
end


%% Summary Noise Table (average across all slices and timePts)

% Find # unique subject/AMPM combos
tempSubjectAMPM = noiseMaster(:,1:2); %select just subject/AMPM columns
[uniqueCombinations,ia,ic] = unique(tempSubjectAMPM, 'rows');
nCombinations = height(uniqueCombinations);

% Aggregate all noise vectors for a patient/AMPM combination and calculate stdev
whichRow = 0; noiseSummary = table(); 
for iCombos = 1:nCombinations
    whichRow = whichRow + 1; 
    
    % Subset data
    rowCriteria = noiseMaster.subject == uniqueCombinations.subject(iCombos) & noiseMaster.AMPM == uniqueCombinations.AMPM(iCombos); 
    tempData = noiseMaster(rowCriteria,:); 
    
    % Aggregate all noise vectors
    noiseVectorAll = []; 
    nSlices = height(tempData); 
    for jRows = 1:nSlices
        noiseVectorAll = vertcat(noiseVectorAll, tempData.noiseVector{jRows}); % I should pre-allocate but not huge performance hit so I'll leave it for now
    end
    
    % Put results in summary table
    noiseSummary.subject(whichRow) = (uniqueCombinations.subject(iCombos)); 
    noiseSummary.AMPM(whichRow) = (uniqueCombinations.AMPM(iCombos));
    noiseSummary.noiseRayleigh(whichRow) = nanstd(noiseVectorAll); 
end
noiseSummary.noise(:) = noiseSummary.noiseRayleigh ./ sqrt((2 - pi/2)); %correction for rayleigh noise in background
noiseSummary = sortrows(noiseSummary, [1,2]); 
% <<<<<<<<<<<< not sure if there are pros/cons to now converting subject/AMPM columns to categorical type??

% Save
if doSave==1
    save(fullfile(savePath, 'mri_noise_summary.mat'),'noiseSummary')
end

