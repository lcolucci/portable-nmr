%% Process MRI Data
%  + This is the heavy work-horse function that analyzes raw mri T2 data in a pixelwise manner. Once this script runs, don't ever need 
%    to go back to raw mri data ever again. Everything can be obtained from the resulting files that it generates. 
%  + It generates mri_pixelByPixel_[1 or 2]exp.mat, mri_t2PixelDecays.mat
%  + This script is a wrapper around the function 'pixel_by_pixel_analysis.m', which does pixel-by-pixel fittings (e.g. 1-exp, 2-exp, ILT fits, etc.)
%  + This script loops through multiple MRI datasets, gives each dataset to the function 'pixel_by_pixel_analysis.m', and then organizes the
%    results in a big structure which contains results for all patients.
%  
%   INPUTS
%       whichIdToAnalyzeLoop - vector of patients to analyze. Take values from 'id' column of config.xlsx (sheet=MRI_Study_T2)
%       pathSaveBase - string. Results will be saved in an additional sub-folder within this path. 
%       whichImageToAnalyzeLoop - cell array of strings
%       whichMaskToProcessLoop - cell array of strings. **This code should really only ever be run on the largest possible mask. 
%                                                         No need to run it again after that. 
%                                                         Exception to this rule would be if you want to do pixel-by-pixel ILT on a small ROI (this analysis takes forever to run).
%       timeLoop - cell array of vectors (time arrays)
%       ignore - # of initial points to ignore when doing T2 fittings (should be = 1, unless explicity testing something different)
%       startVals1ExpLoop - cell array of 1x2 vectors 
%       startVals2ExpLoop - cell array of 1x4 vectors
%       pathToMriStudy - string path to folder that contains the 'MRI_Study' folder
%       pathToWorkBook - string path to config.xlsx file
%       sheetName - which sheet of the config.xlsx file to look at?
%
%       doPixelT2Decays - save the individual T2 decay for each pixel
%       do1ExpFit - pixel-by-pixel 1exp fit
%       do2ExpFit - pixel-by-pixel 2exp fit
%       doSave - do you want to save results? 1 if so. 
%
%   OUTPUTS
%       If doSave==1, results from whichever 'do*' fit options will be
%       saved in the following folder structure: 
%
%       PATHSAVEBASE
%           MASK_NAME
%               mri_pixelByPixel_2exp.mat
%               mri_t2PixelDecays.mat
%               ...(results for each type of analysis for all patient IDs)
%
%               ID1
%               ID2
%                   mri_pixelByPixel_1exp_HC01_PM.mat
%                   mri_pixelByPixel_2exp_HC01_PM.mat
%                   mri_t2PixelDecays_HC01_PM.mat
%                   ...(results for each type of analysis for a single patient ID)
%               ... (folders for each ID that was analyzed. Contains the results from just that one patient ID. This is saved along the way so results are not lost) 
% 
%
%  Lina A. Colucci, 2018     
%
%  Potential Future To-Dos 
%   + Change how results are saved. From a deep structure to a table format: Subject | PM | Slice | Results_2exp (table)
%
%

clear; close all

%% --- USER INPUTS ---
whichIdToAnalyzeLoop = [30];                          % Look at config.xlsx to find ID of the subjects you want to analyze
pathSaveBase = '/Data/Processed_Data';   % Base path to save results to. Results will be saved to an additional sub-folder named after the mask
whichImageToAnalyzeLoop = {'legShorterTE'};              % (cell array) Which column of config.xlsx to analyze? e.g. {'legShorterTE'}, {'fingerLongerTE','additionalScan01'},etc.
whichMaskToProcessLoop = {'all'};                  % (cell array) Script will loop through each pixel of this ROI. Generally, do the largest possible ROI, e.g. {'leg_whole'} or {'all'}
timeLoop = {[8:8:8*32]};                                 % (cell array) time array, cell array of length 1 or same # of vectors as IDs in 'whichIdToAnalyzeLoop' e.g. {[8:8:8*32]}, {[25.5:25.5:25.5*32]}
ignore = 1;                                              % # of initial points to ignore when doing the fittings (generally 1)
% ----
doPixelT2Decays = 1;   % 1 or something else
doSaturatedPixels = 1; % 1 or something else
do1ExpFit = 1;         % 1 or something else
do2ExpFit = 0;         % 1 or something else
do1ExpFitOffset=0;     % 1 or something else
do2ExpFitOffset=0;     % 1 or something else
doILT=0;               % 1 or something else
doSave = 1;   % Save results? 1 or something else
% ---
startVals1ExpLoop = {[2000 80]}; 
startVals2ExpLoop = {[1500 50 1500 210]}; 
startVals1ExpOffsetLoop = {[2000, 50, 50]}; 
startVals2ExpOffsetLoop = {[1500 55 1500 200 50]};
% ---
pathToMriStudy = '/Documents'; % Path to folder that contains 'MRI_Study' folder
pathToWorkbook = '/GitHub/MRI/mri/config.xlsx'; % Path to config.xlsx 
sheetName = 'MRI_Study_T2';                             % Which sheet name in config.xlsx to load?

%% --- END OF USER INPUTS ---

%%  Checks
%      Check that these variables are cell arrays and have a length of 1 or length(whichIdToAnalyzeLoop). If not, throw an error.
%                   >>  variables: timeLoop, whichMaskToProcessLoop, whichImageToAnalyzeLoop, startVals**Loop
nIDsInLoop = length(whichIdToAnalyzeLoop); 
checkinputs_cells(nIDsInLoop,whichImageToAnalyzeLoop,whichMaskToProcessLoop,timeLoop,startVals1ExpLoop,startVals2ExpLoop,startVals1ExpOffsetLoop,startVals2ExpOffsetLoop)

% Load config excel file
config = load_config_file(pathToWorkbook, sheetName);

%% Loop 
for iMask = 1:length(whichMaskToProcessLoop)
    whichMaskToProcess = whichMaskToProcessLoop{iMask}; 
    for jPatient = 1:nIDsInLoop
        %try 
        whichIdToAnalyze = whichIdToAnalyzeLoop(jPatient); 
        whichConfigRow = find(config.id==whichIdToAnalyze);
        
        fprintf('\nCurrently analyzing patientID #%d (%s %s), which is #%d out of %d patients to analyze. \n',config.id(whichConfigRow), config.subject{whichConfigRow}, config.AMPM{whichConfigRow}, jPatient, nIDsInLoop)
        
        %% Check: use the same input for all patients or looping input for each patient? 
        %  Go from 'inputLoop' to 'input'
        [whichImageToAnalyze,TE,startVals1Exp,startVals2Exp,startVals1ExpOffset,startVals2ExpOffset] = checkinputs_static_or_looping(jPatient, whichImageToAnalyzeLoop,timeLoop,startVals1ExpLoop,startVals2ExpLoop,startVals1ExpOffsetLoop,startVals2ExpOffsetLoop);
        
        %% Load image (nifti)
        pathToImageNifti = fullfile(pathToMriStudy, config.pathNiftis{whichConfigRow}, config.(whichImageToAnalyze){whichConfigRow});
        try 
            image = load_nifti(pathToImageNifti); image = image.vol; 
        catch
            image = niftiread(pathToImageNifti); % if not processed with freesurfer
        end 

        %% Load mask (nifti)
        pathToMask = fullfile(pathToMriStudy, config.pathMasks{whichConfigRow},strcat(whichMaskToProcess,'.nii'));
        mask = load_nifti(pathToMask); mask = mask.vol; 

        
        %% Do Analyses: Loop through all slices and pixels in this patient/mask combo
        clear pixelT2Decays_singlePatient result1exp_singlePatient result2exp_singlePatient saturatedPixels_singlePatient
        [pixelT2Decays_singlePatient, saturatedPixels_singlePatient, result1exp_singlePatient, result2exp_singlePatient,result1expOffset_singlePatient, result2expOffset_singlePatient,resultILT_singlePatient] = ...
                                                pixel_by_pixel_analysis(image, mask, TE, ignore, ...
                                                doPixelT2Decays,doSaturatedPixels, do1ExpFit, do2ExpFit,do1ExpFitOffset,do2ExpFitOffset,doILT,...
                                                startVals1Exp,startVals2Exp, startVals1ExpOffset, startVals2ExpOffset); 
        
        %% Put results from a single patient into a larger structure for all patients
        nameSubject = config.subject{whichConfigRow};
        nameAMPM = config.AMPM{whichConfigRow};        
        
        if (doPixelT2Decays == 1)
            pixelT2Decays.(nameSubject).(nameAMPM) = pixelT2Decays_singlePatient; 
        end
        if (doSaturatedPixels == 1)
            saturatedPixels.(nameSubject).(nameAMPM) = saturatedPixels_singlePatient;
        end
        if (do1ExpFit == 1)
            result_1exp.(nameSubject).(nameAMPM) = result1exp_singlePatient; 
        end
        if (do2ExpFit == 1)
            result_2exp.(nameSubject).(nameAMPM) = result2exp_singlePatient;
        end
        if (do1ExpFitOffset == 1)
            result_1expOffset.(nameSubject).(nameAMPM) = result1expOffset_singlePatient; 
        end
        if (do2ExpFitOffset == 1)
            result_2expOffset.(nameSubject).(nameAMPM) = result2expOffset_singlePatient;
        end
        if (doILT==1)
            result_ILT.(nameSubject).(nameAMPM) = resultILT_singlePatient; 
        end
        
        %% Temporary save of individual patient results (in case there's an error before the script finishes running, all the intermediate data gets saved)
        if (doSave==1) && (length(whichIdToAnalyzeLoop) > 1)
            % Create save folder & subfolder named after the mask
            if 7~=exist(fullfile(pathSaveBase, whichMaskToProcess, sprintf('ID%d',whichIdToAnalyze)), 'dir') %exist==7 if name is a folder 
                mkdir(fullfile(pathSaveBase, whichMaskToProcess, sprintf('ID%d',whichIdToAnalyze)));
            end

            % Save
            if (doPixelT2Decays==1)
                save(fullfile(pathSaveBase,whichMaskToProcess, sprintf('ID%d',whichIdToAnalyze), sprintf('mri_t2PixelDecays_%s_%s.mat',nameSubject,nameAMPM)), 'pixelT2Decays_singlePatient')
            end
            if (doSaturatedPixels==1)
                save(fullfile(pathSaveBase,whichMaskToProcess, sprintf('ID%d',whichIdToAnalyze), sprintf('mri_saturatedPixels_%s_%s.mat',nameSubject,nameAMPM)), 'saturatedPixels_singlePatient')
            end
            if (do1ExpFit==1)
                save(fullfile(pathSaveBase,whichMaskToProcess, sprintf('ID%d',whichIdToAnalyze), sprintf('mri_pixelByPixel_1exp_%s_%s.mat',nameSubject,nameAMPM)), 'result1exp_singlePatient')
            end  
            if (do2ExpFit==1)
                save(fullfile(pathSaveBase,whichMaskToProcess, sprintf('ID%d',whichIdToAnalyze), sprintf('mri_pixelByPixel_2exp_%s_%s.mat',nameSubject,nameAMPM)), 'result2exp_singlePatient')
            end    
            if (do1ExpFitOffset==1)
                save(fullfile(pathSaveBase,whichMaskToProcess, sprintf('ID%d',whichIdToAnalyze), sprintf('mri_pixelByPixel_1expOffset_%s_%s.mat',nameSubject,nameAMPM)), 'result1expOffset_singlePatient')
            end  
            if (do2ExpFitOffset==1)
                save(fullfile(pathSaveBase,whichMaskToProcess, sprintf('ID%d',whichIdToAnalyze), sprintf('mri_pixelByPixel_2expOffset_%s_%s.mat',nameSubject,nameAMPM)), 'result2expOffset_singlePatient')
            end  
            if (doILT==1)
                save(fullfile(pathSaveBase,whichMaskToProcess, sprintf('ID%d',whichIdToAnalyze), sprintf('mri_pixelByPixel_ILT_%s_%s.mat',nameSubject,nameAMPM)), 'resultILT_singlePatient')
            end  
        end
        
        clear nameSubject nameAMPM whichIdToAnalyze whichConfigRow whichImageToAnalyze TE startVals1Exp startVals2Exp pathToImageNifti image pathToMask mask 
%         catch
%             fprintf('********** ---ID#%d could not be processed--- ********** \n',whichIdToAnalyzeLoop(jPatient))
%         end
    end    
    
    %% Order fields alphabetically and save
    if (doSave==1)
        % Create save folder & subfolder named after the mask
        if (7~=exist(fullfile(pathSaveBase, whichMaskToProcess), 'dir')) %exist==7 if name is a folder 
            mkdir(fullfile(pathSaveBase, whichMaskToProcess)); 
        end
        
        % Save
        if (doPixelT2Decays==1)
            pixelT2Decays = orderfields(pixelT2Decays); 
            save(fullfile(pathSaveBase, whichMaskToProcess, 'mri_t2PixelDecays.mat'), 'pixelT2Decays')
        end
        if (do1ExpFit==1)
            result_1exp = orderfields(result_1exp); 
            save(fullfile(pathSaveBase, whichMaskToProcess, 'mri_pixelByPixel_1exp.mat'), 'result_1exp')
        end  
        if (do2ExpFit==1)
            result_2exp = orderfields(result_2exp); 
            save(fullfile(pathSaveBase, whichMaskToProcess, 'mri_pixelByPixel_2exp.mat'), 'result_2exp')
        end    
        if (do1ExpFitOffset==1)
            result_1expOffset = orderfields(result_1expOffset); 
            save(fullfile(pathSaveBase, whichMaskToProcess, 'mri_pixelByPixel_1expOffset.mat'), 'result_1expOffset')
        end  
        if (do2ExpFitOffset==1)
            result_2expOffset = orderfields(result_2expOffset); 
            save(fullfile(pathSaveBase, whichMaskToProcess, 'mri_pixelByPixel_2expOffset.mat'), 'result_2expOffset')
        end 
        if (doILT==1)
            result_ILT = orderfields(result_ILT); 
            save(fullfile(pathSaveBase, whichMaskToProcess, 'mri_pixelByPixel_ILT.mat'), 'result_ILT','-v7.3')
        end 
    end
       
end

