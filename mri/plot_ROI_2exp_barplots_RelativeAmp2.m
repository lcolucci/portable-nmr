%% Bar graphs of MRI ROI 2-exp (indiv. subjects): Relative Amp 2 of Muscular Tissue
%  This script loads in the MRI ROI 2-exp data and plots relative amplitude
%  2 values for each subject as a barplot. AM and PM bars for each subject
%  are grouped together. 
%
%  INPUTS
%       doSave - 1 if you want to save figures
%       savePath - folder in which to save figures
%       pathToDrive - path to CimaLabLargeFiles folder in the drive
%
%  OUTPUTS
%       barplot of relative amp 2 for each subjects. Gets saved to
%       savePath (if doSave=1 as poth .eps and .fig 
%
% L.Colucci 2018

clear; close all 

%% ----- USER INPUTS ---- 
doSave = 0; % 1 if you want to save figure
savePath = '/Users/linacolucci/Documents/Thesis/Thesis/Figures/MRI';
pathToDrive = '/Volumes/CimaLabLargeFiles'; % path to CimaLabLargeFiles folder in the drive
%% ---------------------

%% Load and prep data
load(fullfile(pathToDrive,'Hydration/MRI_Study/ProcessedData/mri/ROIs/ROIs_2expFits.mat'))
data = results_2exp.muscle_all; 
% data = data(data.HD==1,:); %<<<<<< If you want just HD's to be plotted, uncomment this line

%% Re-shape data to proper form for a group bar plot
fnSubjects = unique(data.Subject); %fn = fieldname
nSubjects = length(fnSubjects); 
plotData = []; % initialize empty matrix
for iSubject = 1:nSubjects
    tempData = data(data.Subject == fnSubjects(iSubject), :); 
    plotData = [plotData; ...
                tempData.relative_amp2(tempData.AMPM=='AM') tempData.relative_amp2(tempData.AMPM=='PM')]; 
end

%% Plot
figure
bar_handle = bar(fnSubjects, plotData,'LineStyle','none'); 
set(gcf, 'color', 'w','Position',[203, 195, 1196, 760]); set(gca, 'fontsize',23,'linewidth',2); box off
ylabel('Relative Amplitude 2 (%)')
title('MRI Relative Amplitude 2 of Muscular Tissue')
legend('AM','PM'); legend boxoff
keynoteOrange = [255 147 0]./255;
keynoteBlue = [0 162 255]./255; 
bar_handle(1).FaceColor = keynoteBlue; 
bar_handle(2).FaceColor = keynoteOrange;

%% Save 
if doSave==1
    % make folder if it doesn't exist
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end
    
    saveas(gcf, fullfile(savePath, 'bar_mri_relativeAmp2_AMPM.eps'),'epsc')
    saveas(gcf, fullfile(savePath, 'bar_mri_relativeAmp2_AMPM.fig'))
end

