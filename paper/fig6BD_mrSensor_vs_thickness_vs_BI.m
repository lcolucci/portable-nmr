%% Fig 6 - B and D - MR Sensor vs Thickness and BI 
%% B: Plot MR Change in R2 vs BI Change in Re
%% D: Plot Subcu Thickness vs MR Relative Amplitude 3 
%  This script loads in the thicknesses dataset and MR sensor forced 3-exp
%  results and plots total thickness vs. relative amplitude 3. Total
%  thickness = subcu + skin thickness. 
%  
%  PREVIOUS SCRIPTS
%       thickness data comes from script: process_subcu_thickness.m
%       MR sensor data comes from script: analyze_forced3expFits.m    
%       BI data comes from script: process_bioimpedance.m 
%
%  INPUTS
%       doSave - 1 if you want to save figures
%       savePath - folder in which to save figures
%       pathToDrive - path to CimaLabLargeFiles folder in the drive
%
%  OUTPUTS
%       figure of total thickness vs. relative amp 3. Gets saved in
%       savePath if doSave=1 as poth .eps and .fig 
%
% L.Colucci 2018

clear; close all 

%% ----- USER INPUTS ----
pathToDrive = '/Volumes/CimaLabLargeFiles'; % path to CimaLabLargeFiles folder in the drive
posVector = [560   528   1*560   1*420]; % position vecotr 
markerSize = 16; 
fontSize = 26; %for labels and axes
% ---
doSave = 0; % 1 if you want to save figure
savePath = '/Users/linacolucci/Documents/Thesis/Paper/figures/MRsensor';
%% ---------------------

%% Load Datasets
% MR Sensor
load(fullfile(pathToDrive,'Hydration/MRI_Study/ProcessedData/mr_sensor_leg/legSensor_forced3exp_40ms_250ms.mat'))
result.SubjectName = categorical(result.subject); 
data_MRsensor = result; 
change_MRsensor = varfun(@diff, result, 'GroupingVariables','subject','InputVariables',{'relax','relative_amp1','relative_amp2','relative_amp3'}); 
change_MRsensor.SubjectName = categorical(change_MRsensor.subject); 
clear result
% Thickness
load(fullfile(pathToDrive,'Hydration/MRI_Study/ProcessedData/other/thicknesses_summary_mean_std.mat'))
data_thickness = data_summary; 
clear data_summary
% BI 
load(fullfile(pathToDrive,'Hydration/MRI_Study/ProcessedData/bioimpedance/BI_Averages.mat'))
data_BI = bi(bi.WholeBody==0,:); 
data_BI.SubjectName = categorical(data_BI.SubjectName); 
 
%% ---- Fig 6B: MR Sensor Change in R2 vs. BI Change in Re -----
%  code came from script: plot_BI_vs_mrSensor.m

% Prep BI data
list_subjects = unique(data_BI.SubjectName); 
BI_change_results = table(); 
BI_results = table(); 
% Need to have this loop, instead of just applying varfun because some
% subjects don't have a 'last' value and so varfun returns an error when that's the case
for ii=1:length(list_subjects)
    % subset data, just 1 subject
    tempData = data_BI(data_BI.SubjectName==list_subjects(ii),:); 
    
%     % Find pre and post time points
%     tempData = tempData(tempData.TimePt~=90 & tempData.TimePt~=98,:); 
%     [amTimePt, whichRowAM] = min(tempData.TimePt);
%     [pmTimePt, whichRowPM] = max(tempData.TimePt); %(tempData.TimePt~=90 & tempData.TimePt~=98)) 
%     % subset data, just pre and post timePts
%     tempData = tempData([whichRowAM, whichRowPM],:); 
    
    % subset data
    tempData = tempData(tempData.TimePt_Last == 1 | tempData.TimePt_First == 1, :); 
    
    % Calculate pre-post difference
    tempResult = varfun(@diff, tempData, 'GroupingVariables','SubjectName','InputVariables',{'Rinf','Re','Ri','Rcentre'}); 
    
    % Aggregate with master
    BI_change_results = vertcat(BI_change_results, tempResult); 
    BI_results = vertcat(BI_results, tempData); 
    
    % Delete temp var
     clear tempResult tempData whichRow* 

end

% Plot: Change in Leg BI vs MR Sensor
param1 = 'diff_Re'; 
plotData = join(BI_change_results,change_MRsensor,'Keys',{'SubjectName','GroupCount'});
figure; 
plot(-1.*plotData.(param1), plotData.diff_relative_amp2,'ko','MarkerSize',markerSize,'MarkerFaceColor','k')
l = lsline; 
l.LineWidth = 3; 
xlabel('Calf BI: \DeltaECF (\Omega)')
ylabel('MR Sensor: \DeltaR_2 (%)')
rsquared(-1.*plotData.(param1), plotData.diff_relative_amp2)
text(-26, 0, sprintf('r^2 = %.3f', rsquared(-1.*plotData.(param1), plotData.diff_relative_amp2)),'FontSize',fontSize)
set(gcf, 'color', 'w','Position',posVector);
set(gca, 'fontsize',fontSize,'linewidth',2,'XColor','k','YColor','k'); box off 
%title(sprintf('MR Sensor Change R_2 vs BI Change R_e'))
xlim([-30, -5])
ylim([-10 3])

% Save
if doSave==1
    % make folder if it doesn't exist
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end
    saveas(gcf, fullfile(savePath, 'fig7D_mrSensorChangeR2_vs_BIchangeRe.eps'),'epsc')
end



%% ---- Fig 6D: MR Sensor R3 vs. Thickness -----
% Came from script: plot_thickness_vs_mrSensor.m

% Prep Data
% Select Subset
data_thickness = data_thickness(data_thickness.Tissue == 'total',:); % select the total thickness data
% Join thickness and MR sensor data into single table
data_MRsensor.Subject_Name = categorical(data_MRsensor.subject); % change to categorical to match 'subject_name' of thickness data (need this for outerjoin)
data = outerjoin(data_thickness, data_MRsensor);
% Select just AM data
plotData = data(data.PM==0,:); %Select just AM

% Plot
figure; 
plot(plotData.mean_Thickness_mm, plotData.relative_amp3,'ko','MarkerSize',markerSize,'MarkerFaceColor','k')
set(gcf, 'color', 'w','Position',posVector); %[74, 68, 1550, 850]
set(gca, 'fontsize',fontSize,'linewidth',2,'XColor','k','YColor','k'); box off 
l = lsline;
l.LineWidth = 3; 
xlabel('Thickness (mm)')
ylabel('R_3 (%)')
text(4, 33, sprintf('r^2 = %.3f', rsquared(plotData.mean_Thickness_mm, plotData.relative_amp3)),'FontSize',fontSize)
%title(sprintf('MR Sensor Pre R_3 vs Thickness, r^2=%.4f',rsquared(plotData.mean_Thickness_mm, plotData.relative_amp3)))
%title({'Relative Amplitude_3 of the AM MRI Scan',' vs.', 'Total Thickness of Skin and Subcutaneous Tissue',sprintf('r^2=%.4f',rsquared(plotData.mean_Thickness_mm, plotData.relative_amp3))})

% Save
if doSave==1
    % make folder if it doesn't exist
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end
    
    saveas(gcf, fullfile(savePath, 'fig7C_mrSensorR3_vs_thickness.eps'),'epsc')
    %saveas(gcf, fullfile(savePath, 'fig7C_mrSensorR3_vs_thickness.fig'))
end


