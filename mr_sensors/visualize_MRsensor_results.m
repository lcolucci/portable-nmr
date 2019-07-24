%% 	Visualize Results of MR Sensor
%	A type of scatter plot to visualize MR Sensor Data
%  		x-axis has 'Pre' and 'Post' markings
%		y-axis has RA2 values for each subject
%	Same as fig.3A except for MR sensor data instead of MRI data

%% --- USER INPUTS ---
pathToMriStudy = '/Volumes/CimaLabLargeFiles/Hydration'; % Path to folder that contains 'MRI_Study' folder (contains ROI_2expFits, etc.) 
whichParamsList = {'relative_amp2'}; % Which parameters do you want to plot / analyze? i.e. {'relax1','relax2','relative_amp2'}
doSave = 0; 
savePath = '/Users/linacolucci/Documents/Thesis/Thesis/Figures/MR_Sensor'; 
%% -------------------


%% Load and process MR Sensor Data
load(fullfile(pathToMriStudy,'MRI_Study/ProcessedData/mr_sensor_leg/legSensor_forced3exp_40ms_250ms.mat'))
data = result; 
data.Properties.VariableNames{24} = 'Subject'; %re-name 'subject' to 'Subject'
data.HD(contains(data.Subject,'HC')) = 0; 
data.HD(contains(data.Subject,'HD')) = 1; 

%% Create plots
plot_multiexp_results(data, whichParamsList) ;
ylabel('Relative Amplitude 2 (%)')
title('MR Sensor')
% Save
if doSave==1  
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end 
    saveas(gcf, fullfile(savePath, 'results_MRsensor_Amp2.eps'),'epsc')
end
