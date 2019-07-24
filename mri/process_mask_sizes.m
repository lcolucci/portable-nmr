%%  Process Mask Sizes
%   This script calculates the # of pixels of each mask for each patient 
%
%   USER INPUTS
%       whichIdToAnlayzeLoop = vector of patients to analyze. Take values from 'id' column of config.xlsx (sheet=MRI_Study_T2)
%       pathToWorkBook - string path to config.xlsx file
%       sheetName - which sheet of the config.xlsx file to look at?
%       pathToMriStudy - string path to folder that contains the 'MRI_Study' folder
%       doSave - do you want to save results? 1 if so. 
%       savePath - path to folder where 'maskSizes' will be saved
%
%   OUTPUTS
%       maskSizes - a table with # of pixels in each ROI for each subject
%           Columns: subject | AMPM | slice | {names of ROIs} 
%
clear 

%% --- USER INPUTS ---
whichIdToAnalyzeLoop = [1:28];  % Look at config.xlsx to find IDs of the subjects you want to analyze
% ---
pathToWorkbook = '/Users/linacolucci/Documents/GitHub/MRI/mri/config.xlsx'; % Path to config.xlsx (contains paths to load mask files) 
sheetName = 'MRI_Study_T2';                             % Which sheet name in config.xlsx to load?
pathToMriStudy = '/Volumes/CimaLabLargeFiles/Hydration'; % Path to folder that contains 'MRI_Study' folder (contains mask files) 
% ---
doSave = 1; %1 for yes
savePath = '/Volumes/CimaLabLargeFiles/Hydration/MRI_Study/ProcessedData/mri'; % This is where 'maskSizes' will get saved (doMaskSize). 
                                                                               % results from doMaskMerge & doSNR will overwrite and replace the original results (mri_pixelByPixel_1exp.mat, or mri_pixelByPixel_2exp.mat)
%% ---- END USER INPUTS ----

%% Load config.xlsx file
config = load_config_file(pathToWorkbook, sheetName);

%% Analysis
nIDsInLoop = length(whichIdToAnalyzeLoop); 
whichRowSize = 0; maskSizes = table(); 
for iIDs=1:nIDsInLoop
    whichConfigRow = find(config.id==whichIdToAnalyzeLoop(iIDs)); 
        
    % Directory of masks
    pathToMaskFolder = fullfile(pathToMriStudy, config.pathMasks{whichConfigRow}); 
    pathToMaskFolderSearchNiftis = fullfile(pathToMaskFolder,'*.nii');
    maskList = dir(pathToMaskFolderSearchNiftis); 

    % Find noise mask and delete it from list
    noiseRow = contains({maskList.name},'noise','IgnoreCase',true); 
    maskList(noiseRow,:) = [];         

    %% Loop through all masks found in folder
    nMasks = length(maskList);
    for kMasks = 1:nMasks
        % Load mask
        mask = load_nifti(fullfile(pathToMaskFolder,maskList(kMasks).name)); 
        mask = mask.vol; 

        % Define maskName
        maskName = maskList(kMasks).name; 
        maskName = strrep(maskName,'.nii',''); 

        %% Loop through slices
        nSlices = size(mask,3); 
        if nSlices ~= 4 % warning sign
            fprintf('\n*** NOTE *** \n%s %s %s\nThis mask does not have 4 slices as expected.\nMake sure that this is okay and that the slice number that''s getting recorded in maskSizes table is accurate.\n********\n',config.subject(whichConfigRow), config.AMPM(whichConfigRow),maskName)
        end
        for lSlices = 1:nSlices 

            % Select 1 slice of mask    
            tempMask = mask(:,:,lSlices); 

            % Calculate Mask Size
            whichRowSize = whichRowSize +1; 
            maskSizes.subject(whichRowSize) = config.subject(whichConfigRow); 
            maskSizes.AMPM(whichRowSize) = config.AMPM(whichConfigRow);
            maskSizes.slice(whichRowSize) = lSlices; % <-- if fprintf NOTE gets activateda above, make sure slice # is what you want it to be 
                                                     % i.e. there will be issues if there is only slice # 1,3,4. It'll get recorded as 1,2,3. 
            maskSizes.mask(whichRowSize) = {maskName}; 
            maskSizes.nPixels(whichRowSize) = sum(sum(tempMask));

            clear tempMask
        end
        clear maskName mask
    end
    clear  whichConfigRow pathToMaskFolder pathToMaskFolderSearchNiftis maskList noiseRow
end
%% Convert Table Wide Format
maskSizes = unstack(maskSizes,'nPixels','mask');

%% Save
if doSave==1
    save(fullfile(savePath, 'mask_sizes.mat'), 'maskSizes')
end
    