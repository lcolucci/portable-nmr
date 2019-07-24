%% Heatmap Pixel-by-Pixel Results
%   This is a wrapper around the function 'visualize_pixelbypixel_heatmap.m' that applies the function to several
%   patients and saves the figures as well
%
%   INPUTS
%       load results - load in results_1exp or results_2exp 
%       results - pick which fitting type you want to plot (results_1exp or results_2exp)
%       whichParameters - which parameters from the results_2exp/1exp tables to plot? Up to 2 parameters at a time (i.e. relax2, r2, etc.)
%       savePath - path to folder where figures should be saved
%       doSave - do you want to save the figures? saved as both .png and .eps
%       pathGithubFolder - path to where the MRI Github folder (eventually to be renamed PhD) is located. This is so the colormaps can be loaded.
%
%   OUTPUTS
%       figures for each subject/AMPM combination saved in the specified folder
%

%% 
clear; close all 

% Load results

% Load colormaps
pathGithubFolder = '/GitHub'; % --- USER --- path of computer to the 'MRI' github folder
colormapBaseFolder = fullfile(pathGithubFolder, 'mri/functions/colormaps'); 
colormapList = dir(fullfile(colormapBaseFolder, '*.mat')); 
for iCMaps = 1:length(colormapList)
    load(fullfile(colormapBaseFolder, colormapList(iCMaps).name))
end

%% ------- USER INPUTS -------
results = result_2exp; 
whichParameters = {'relative_amp2'}; % which parameters to plot in heatmaps? up to 2 parameters, i.e. {'relax1','relax2'}, {'relax1','r2'}, or {'r2'}
colorLimits = {[0 100]}; % cell array, i.e. {[0 500], [150 1000]}, {[0 1000]}, or {} <-- default limits will set in
colormaps = {cmap_amp2}; %relaxation times: cmap_relax1, cmap_relax1b, cmap_relax2d
% --- Save ---
savePath = 'Thesis/Figures/MRI/Pixel_by_Pixel/Results_2exp/Heatmaps/amp2'; % include '2exp' or '1exp' in the sub-folder name 
doSave = 1; % 1 or something else
%% -------------------------------


%% Loop through subjects to generate figures
fnSubjects = fieldnames(results); % fn = fieldnames 
nSubjects = length(fnSubjects); 
for iSubjects = 1:nSubjects
    fnAMPM = fieldnames(results.(fnSubjects{iSubjects})); 
    nAMPM = length(fnAMPM); 
    
    for jAMPM = 1:nAMPM
        tempResults = results.(fnSubjects{iSubjects}).(fnAMPM{jAMPM}); %structure at 'slice' level
        fnSlices = fieldnames(tempResults); 
        nSlices = length(fnSlices); 
        
        %% Calculate relative amp 
        
        %% Check data and delete unwanted pixels <<<<<< checks happen in 'visualize_pixelbypixel_heatmap'

        %% Visualize Heatmaps 
        optionalTitle = sprintf('%s %s',fnSubjects{iSubjects},fnAMPM{jAMPM}); 
        visualize_pixelbypixel_heatmap(tempResults, whichParameters,colormaps, optionalTitle,colorLimits); %<--- I SHOULD APPLY CHECKS WITHIN THIS FUNCTION. AT LEAST GIVE OPTION OF HAVING A CHECK. 
        
        %% Save
        if doSave==1
            if ~(7==exist(savePath)) % exist = 7 when it's a folder
                mkdir(savePath)
            end
%             saveas(gcf, fullfile(savePath, sprintf('heatmap_%s_%s.png', fnSubjects{iSubjects}, fnAMPM{jAMPM})))
            saveas(gcf, fullfile(savePath, sprintf('heatmap_%s_%s.eps', fnSubjects{iSubjects}, fnAMPM{jAMPM})),'epsc')
        end
    end
end

