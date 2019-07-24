%% ------- Figure 6A,C, Table S8: MR Sensor Relative Amp 2 ------
% Script used to be called plot_mrSensor_boxplots.m (located in mr_sensors folder)
%
% This script creates all figures/tables for pure MR sensor results: 
%       fig 6 C: 4 boxplots
%       fig 6 A: 2 boxplots
%       table S8: summary of MR sensor results 
%       supplemental fig: "line plot" for MR sensor (AM-PM values for each individual)
%
% PREVIOUS SCRIPTS
%       data comes from script: analyze_forced3expFits.m
%
%
% OUTPUTS
%       Fig 6 Subplot C: 4 boxplots
%       Fig 6 Subplot A: 2 boxplots
%       Table S8: summary of MR sensor results 
%
% STATISTICS
%       the accurate statistics comes from R script (statistics_for_paper.R)
%       make sure to populate p-values from R script into final figures for publication
%       have not found a reliable permutation test function for Matlab
%
% Lina A. Colucci, 2018 


clear

%% ----- USER INPUTS ----
whichParameter = 'relative_amp2'; % Which parameter to plot / create table for? 
pathToMriStudy = '/Volumes/CimaLabLargeFiles/Hydration'; % Path to folder that contains 'MRI_Study' folder 
% ---
doSave = 1; % 1 if you want to save figure
savePath = '/Users/linacolucci/Documents/Thesis/Paper/figures/MRsensor';
%% ---------------------


%% Load and prep data. Calculate uniqueID and AM-PM change. 
load(fullfile(pathToMriStudy, 'MRI_Study/ProcessedData/mr_sensor_leg/legSensor_forced3exp_40ms_250ms.mat'))
data = result; 
data.HD(contains(data.subject,'HD'))=1;    
dataChange = varfun(@diff, data, 'GroupingVariables','subject','InputVariables',{'relax','relative_amp1','relative_amp2','relative_amp3'}); 
dataChange.HD(contains(dataChange.subject,'HD'))=1;
dataChange.SubjectName = categorical(dataChange.subject);
data.subject = categorical(data.subject); 
data.Properties.VariableNames{24} = 'Subject';
data.uniqueID = findgroups(data(:,[25,26])); % form unique grouping from PM and HD columns
data.AMPM(data.PM==1) = {'PM'}; % calculate_ttest function needs a column called 'AMPM'
data.AMPM(data.PM==0) = {'AM'};
data.AMPM = categorical(data.AMPM); % and it needs to be categorical
data.ROI(:) = categorical({'mrSensor'}); % calculate_ttest function needs a column called 'ROI'

%% Define colors
keynoteBlue = [0 162 255]./255;
keynoteOrange = [255 147 0]./255;
grey = [138 138 138]./255;
ishaRed = [142 51 49]./255; % 157 43 44  
ishaBlue = [32 25 77]./255; % 45 38 90
ishaGreen = [47 85 59]./255;
ishaTeal = [98 144 128]./255; % 85 165 139
scienceRed = [197, 22, 29]./255;

%% -------------- Subplot 6C: four boxplots,  AM: HC vs HD.  PM: HC vs HD ------------------
%  this is identical to script from fig5_tabl2_mri_roi.m 

f = customBoxplot4(data, whichParameter);
ylabel('R_2 (%)')

% Save
if doSave==1  
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end 
%   saveas(fig, fullfile(savePath,sprintf('results_2exp_%s.png',whichROI))) % figure is too big to be saved normally, so I'm just going to save as .eps
    saveas(gcf, fullfile(savePath, sprintf('fig7A_MRsensor_%s_AMvsPM_boxplots4.eps',whichParameter)),'epsc')
end

%% ------------ Subplot 6A: two boxplots, Change: HC vs HD --------------
whichParameter2 = strcat('diff_',whichParameter); 
f = customBoxplot2(dataChange, whichParameter2, dataChange.(whichParameter2), dataChange.HD);
ylabel('R_2 (%)')

% Save
if doSave==1  
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end 
%   saveas(fig, fullfile(savePath,sprintf('results_2exp_%s.png',whichROI))) % figure is too big to be saved normally, so I'm just going to save as .eps
    saveas(gcf, fullfile(savePath, sprintf('fig7B_MRsensor_ROI_%s_Change_boxplots2.eps',whichParameter)),'epsc')
end

%% --------- Table S8: Summary of MR Sensor R2 Results ---------
% Code taken from script process_mrSensor_tableForPaper.m 
% Columns of summary: 
%       Subject | HD | AM | PM | Change and corresponding p-values for each indicator

summary = calculate_ttest(data, {whichParameter});
% summary.tTest{1,1}(end,:) = []; % delete pVals (last row)
summary = summary.tTest{1,1}(:,1:5); %delete cols after 'Change'

% Save as csv
if doSave==1  
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end
    writetable(summary,fullfile(savePath, sprintf('table5_MRsensor_%s_TableSummary.csv',whichParameter)))
end
