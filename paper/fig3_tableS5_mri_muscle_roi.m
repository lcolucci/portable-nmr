%% ------- Figure 3 and Table S5: MRI ROI 2-exp Muscle_All ------
%% Fig 3  - boxplots and AM-PM "line" plots
%% Table S5 - summary with p-vals
%
% This script generates the figures/tables corresponding to 
% MRI Muscle_All ROI 2-exp long relative amplitude (R_L). 
%
% PREVIOUS SCRIPT
%       plot data comes from script: process_ROImeans.m
%       tTest data for Table 2 comes from script: analyze_ROI_multiexp_results.m 
%               NOTE: pVals are for a two-sample, two-tailed t-test with unequal vars)  
%       real statistics comes from R script: statistics_for_paper.R
%
% INPUTS
%       whichROI - which ROI do you want to plot data for? i.e. muscle_all, subcu_all, muscle_lateral, leg_whole, etc. 
%       whichParameter - which parameter do you want to plot data for? relative_amp2, relax1, relax2
%       pathToDrive - path to CimaLabLargeFiles folder (different than most scripts that ask for path to Hydration project folder)
%       doSave - 1 to save figures and table
%       savePath - where to save outputs
%
% OUTPUTS
%       Subplot A - "line plot" with AM and PM values for each individual subject
%       Subplot B - figure with 4 total boxplots: 2 groups (AM and PM), and 2 boxplots within each group (HC and HD)
%       Subplot C - figure with 2 total boxplots: AM-PM change in HC subjects, AM-PM change in HD subjects
%       Table 2 - table with one row for each subject and columns: Subject | AM | PM | Change
%
% STATISTICS
%       the accurate statistics comes from R script
%       make sure to populate p-values from R script into final figures for publication
%       there is no reliable permutation test function for Matlab
%
% Lina A. Colucci, 2018


clear; close all

%% ----- USER INPUTS ----
whichROI = 'muscle_all'; 
whichParameter = 'relative_amp2'; 
pathToDrive = '/Volumes/CimaLabLargeFiles'; % path to CimaLabLargeFiles folder in the drive
% ---
doSave = 1; % 1 if you want to save figure
savePath = '/Users/linacolucci/Documents/Thesis/Paper/figures/MRI_ROI';
%% ---------------------


%% Load and prep data
load(fullfile(pathToDrive,'Hydration/MRI_Study/ProcessedData/mri/ROIs/ROIs_2expFits.mat'))
data = results_2exp.(whichROI); 
data.uniqueID = findgroups(data(:,[21,23])); % form unique grouping from AMPM and HD columns


%% --------- Subplot A: "Line" plot of AM vs PM value for each individual subject --------------------
% this is script based on analyze_ROI_multiexp_results.m (which is more flexible and allows generation of lots of different ROIs at once)

% define colors
ishaRed = [142 51 49]./255; % 157 43 44  
scienceRed = [197, 22, 29]./255;

% Create plot
plot_multiexp_results(data, whichParameter, 0, 22, 2, {'w',scienceRed},{'k','k'}) ;
set(gca, 'linewidth',2,'XColor','k','YColor','k'); box off
ylabel('R_L (%)','Interpreter','tex')
title(sprintf('%s',whichROI),'Interpreter','none')
%mtit(sprintf('%s',whichROI),'fontsize',22,'xoff',0,'yoff',0.03,'Interpreter','none');

% Save
if doSave==1  
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end 
%   saveas(fig, fullfile(savePath,sprintf('results_2exp_%s.png',whichROI))) % figure is too big to be saved normally, so I'm just going to save as .eps
    saveas(gcf, fullfile(savePath, sprintf('fig5A_MRI_ROI_%s_%s_IndivLinePlots.eps',whichParameter, whichROI)),'epsc')
end


%% -------------- Subplot B: four boxplots,  AM: HC vs HD.  PM: HC vs HD ------------------
f = customBoxplot4(data, whichParameter);
ylabel('R_L (%)')
ylim([11 45])

% Save
if doSave==1  
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end 
%   saveas(fig, fullfile(savePath,sprintf('results_2exp_%s.png',whichROI))) % figure is too big to be saved normally, so I'm just going to save as .eps
    saveas(gcf, fullfile(savePath, sprintf('fig5B_MRI_ROI_%s_%s_AMvsPMboxplots.eps',whichParameter, whichROI)),'epsc')
end

%% -------------- Subplot C: two boxplots, Change: HC vs HD --------------
% Prep data
changeData = varfun(@diff, data,'InputVariables',{'relax1','relax2','relative_amp2'}, 'GroupingVariables','Subject'); 
changeData.HD = contains(string(changeData.Subject),'HDMRI'); 
x = changeData.(strcat('diff_',whichParameter)); 
groups = changeData.HD; 

% Plot
figure
h=axes; 
boxplot(x, groups)
set(gcf, 'color', 'w');
set(gca, 'fontsize',32,'linewidth',2); box off

% customize colors and lines
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects 
boxesHC = [a(2),a(4),a(6),a(8),a(10),a(12),a(14)]; 
set(boxesHC, 'Color','k','LineWidth',3,'LineStyle','-')
boxesHD = [a(1),a(3),a(5),a(7),a(9),a(11),a(13)]; 
set(boxesHD, 'Color','k','LineWidth',3,'LineStyle','-')

% labels
ylabel('\Delta R_L (%)')
set(gca,'xticklabel',{'HC','HD'},'XColor','k','YColor','k')

% calculate significance
diff_hc = changeData.diff_relative_amp2(changeData.HD==0); 
diff_hd = changeData.diff_relative_amp2(changeData.HD==1); 
pVal = permutationTest(diff_hc, diff_hd, 1e5-1, 'exact',1); %% I don't trust that these perm tests are working

% add significance stars to plot
model_series = [quantile(diff_hc, 0.95), quantile(diff_hd, 0.95)]; 
sigstar2(model_series, get(gca,'XTickLabel'),[pVal],0,'boxplot',h)


% Save
if doSave==1  
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end 
%   saveas(fig, fullfile(savePath,sprintf('results_2exp_%s.png',whichROI))) % figure is too big to be saved normally, so I'm just going to save as .eps
    saveas(gcf, fullfile(savePath, sprintf('fig5C_MRI_ROI_%s_%s_ChangeBoxplots.eps',whichParameter, whichROI)),'epsc')
end


%% ---------- Table 2 - Summary of Individual Values and p-vals -------------------
% NOTE: the pVals are Welch Test values

% Load and prep data
clear data
load(fullfile(pathToDrive,'Hydration/MRI_Study/ProcessedData/mri/ROIs/ROIs_2expFits_tTestResults.mat'))
data = tTestResults; 
data.Parameter = categorical(data.Parameter); 
data = data.tTest(data.ROI==whichROI & data.Parameter==whichParameter); 
data = data{1,1}; 
data = data(:,1:5); % delete everything after 'Change'

% Save as csv
if doSave==1  
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end 
    writetable(data, fullfile(savePath, sprintf('table2_MRI_ROI_%s_%s_TableSummary.csv',whichParameter,whichROI)))
end


