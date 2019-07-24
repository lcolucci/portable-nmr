%%  Analyze ROIs 1-exp and 2-exp Fit Results
%   This script takes multi-exp fitting results and processes them by
%   (1) visualizing AM vs PM results in figures (i.e. "line plots" for each individual's AM and PM values), 
%   and/or (2) doing HC vs HD t-tests and summarizing the results in tables
%
%   PREVIOUS SCRIPTS
%       data comes from script: process_ROImeans.m
%
%   USER INPUTS
%       whichResults - '1exp' or '2exp' - will load the appropriate MRI ROI results
%       whichParamsList - which parameters to plot / analyze? cell array of strings
%       doTTest - do t-test analysis and summary tables? 
%       doFigures - generate plots of fit results? 
%       pathToMriStudy - path to folder containing 'MRI Folder', enables fit results to be automatically loaded
%       doSave - do you want to save results? 
%       savePath - where to save results? 
%
%   OUTPUTS
%       if doSave==1, then tTestResults and/or figures will be saved in the savePath
%           tTestResults - table with columns ROI | Parameter | tTest (table) | FitType (1exp/2exp)
%                          tTest = table with the following columns: 
%                                   Subject 
%                                   HD (0/1) 
%                                   AM (values of Param} 
%                                   PM (values of Param) 
%                                   Change (PM-AM values, + is increase from AM to PM, - is decrease from AM to PM) 
%                                   rmse AM (rmse of the AM fit for this subject and ROI) 
%                                   rmse PM ( " " PM fit)
%                                   r2 AM (r2 of the AM fit for this subject & ROI) 
%                                   r2 PM (" " PM fit)
%                                   ROI (which ROI is this table for)
%                                   Parameter (which paramter is this table for, i.e. 'relax1', 'relative_amp2')
%                                   
%           figures - 1 figure per ROI with # subplots = # of parameters in 'whichParamsList)
%                     Each subplot has AM,PM on x-axis with scatter plot of the value on y-axis. 
%                     Each subject's point is labeled 

clear, close all

%% ---- USER INPUTS ----
whichResults = '2exp'; %'1exp' or '2exp'
whichParamsList = {'relax1','relax2','relative_amp2'}; % Which parameters do you want to plot / analyze? i.e. {'relax1','relax2','relative_amp2'}
doTTest = 1; % Calculate tables with t-tests? 
doFigures = 0; % Plot figures for each ROI? 
pathToMriStudy = '/Volumes/CimaLabLargeFiles/Hydration'; % Path to folder that contains 'MRI_Study' folder (contains ROI_2expFits, etc.) 
% ---
doSave=1; % 1 for yes
savePath = '/Volumes/cimalablargefiles/Hydration/MRI_Study/ProcessedData/mri/ROIs'; % Results will be saved in this folder
%% ---------------------

%% Load data
if strcmpi(whichResults,'1exp')
    load(fullfile(pathToMriStudy,'MRI_Study/ProcessedData/mri/ROIs/ROIs_1expFits.mat'))
    data = results_1exp; 
end
if strcmpi(whichResults,'2exp')
    load(fullfile(pathToMriStudy,'MRI_Study/ProcessedData/mri/ROIs/ROIs_2expFits.mat'))
    data = results_2exp; 
end
  
%% Loop through ROIs
listROIs = fieldnames(data); 
nROIs = length(listROIs); 
tTestResults = table(); 
for iROIs =3:nROIs
    
    %% Pick 1 ROI
    whichROI = listROIs{iROIs}; 
    dataSubset1 = data.(whichROI);
    
    %% Calculate t-test summary table for this 1 ROI 
    if doTTest==1
        try
        temp_tTest = calculate_ttest(dataSubset1, whichParamsList); 
        if isempty(temp_tTest)
            % if there was some reason t-test couldn't be done, skip
        else
            tTestResults = vertcat(tTestResults, temp_tTest); % master table, aggregate results across all ROIs
        end
        catch
        end
    end
    
    %% Create plots
    if doFigures == 1
        plot_multiexp_results(dataSubset1, whichParamsList) ;
        mtit(sprintf('%s',whichROI),'fontsize',22,'xoff',0,'yoff',0.03,'Interpreter','none');
        
        % Save
        if doSave==1  
            if ~(7==exist(savePath)) % exist = 7 when it's a folder
                mkdir(savePath)
            end 
%           saveas(fig, fullfile(savePath,sprintf('results_2exp_%s.png',whichROI))) % figure is too big to be saved normally, so I'm just going to save as .eps
                                                                                    % If I want to debug the .png thing later, look at this answer: https://www.mathworks.com/matlabcentral/answers/102382-how-do-i-specify-the-output-sizes-of-jpeg-png-and-tiff-images-when-using-the-print-function-in-mat
            saveas(gcf, fullfile(savePath, sprintf('results_%s_%s.eps',whichResults, whichROI)),'epsc')
        end
    end
   
    % Clear intermediate variables
    clear temp_tTest whichROI dataSubset1
end

% Add column to show which fit type this is 
if doTTest==1
    tTestResults.FitType(:) = {whichResults};
end
%% Save 
if doSave==1 && doTTest==1
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end 
    save(fullfile(savePath, sprintf('ROIs_%sFits_tTestResults.mat',whichResults)), 'tTestResults')
end

