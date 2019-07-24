% Make Supplemental Figure S13: Raw T2 decay plots for HC and HD subjects
%
% L.Colucci 2019

%clear; close all 

%% ----- USER INPUTS ----
% Load data_table
load('/Volumes/cimalablargefiles/Hydration/MRI_Study/ProcessedData/mr_sensor_leg/legSensor_T2Decays_table.mat')
whichHC = 30; %row corresponding to HC subject
whichHD = 66; %row corresponding to HD subject
% ---
doSave = 1; % 1 if you want to save figure
savePath = '/Documents/STM/Figures';
%% ---------------------

%% Define Plot Appearance 
posVector = [560   528   1*560   1*420]; % position vecotr 
markerSize = 16; 
fontSize = 26; %for labels and axes

% Define Colors
keynoteBlue = [0 162 255]./255;
keynoteOrange = [255 147 0]./255;
grey = [138 138 138]./255;
ishaRed = [142 51 49]./255; % 157 43 44  
ishaBlue = [32 25 77]./255; % 45 38 90
ishaGreen = [47 85 59]./255;
ishaTeal = [98 144 128]./255; % 85 165 139
scienceRed = [197, 22, 29]./255;

%% HC Plot 
% Define data
data_hc = data_table.T2Decay{whichHC,1}; 
time_hc = data_table.Time{whichHC,1}; 

% Plot
figure; plot(time_hc, data_hc, 'Color','k', 'LineWidth',2)
set(gcf, 'color', 'w','Position',posVector); %[74, 68, 1550, 850]
set(gca, 'fontsize',fontSize,'linewidth',2,'XColor','k','YColor','k'); box off 
xlabel('Time (ms)')
ylabel('Amplitude (a.u.)')
title('Healthy Control (HC)')
xlim([-1 max(time_hc)])
ylim([-2 max(data_hc)])

% Save
if doSave==1
    % make folder if it doesn't exist
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end
    saveas(gcf, fullfile(savePath, 'suppl_NMRsensor_RawT2_HC.eps'),'epsc')
end

%% HD Plot 
% Define data
data_hd = data_table.T2Decay{whichHD,1}; 
time_hd = data_table.Time{whichHD,1}; 

% Plot
figure; plot(time_hd, data_hd, 'Color',ishaRed, 'LineWidth',2)
set(gcf, 'color', 'w','Position',posVector); 
set(gca, 'fontsize',fontSize,'linewidth',2,'XColor','k','YColor','k'); box off 
xlabel('Time (ms)')
ylabel('Amplitude (a.u.)')
title('Dialysis Patient (HD)')
xlim([-1 max(time_hd)])
ylim([-2 max(data_hd)])

% Save
if doSave==1
    % make folder if it doesn't exist
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end
    saveas(gcf, fullfile(savePath, 'suppl_NMRsensor_RawT2_HD.eps'),'epsc')
end

%% Analyze SNR values of NMR data 
load('/Volumes/cimalablargefiles/Hydration/MRI_Study/ProcessedData/mr_sensor_leg/legSensor_forced3exp_40ms_250ms.mat')

snr_data = table(result.subject, result.PM, result.snr);
snr_data.Properties.VariableNames = {'subject','PM','snr'};
snr_data.HD(contains(snr_data.subject, 'HDMRI')) = 1; 

snr_mean = mean(snr_data.snr) %75.7964
snr_std = std(snr_data.snr) %26.5781

hc_mean = mean(snr_data.snr(snr_data.HD == 0)) % 80.4728
hc_std = std(snr_data.snr(snr_data.HD == 0)) %24.5560
hc_min = min(snr_data.snr(snr_data.HD == 0)) %42.7225
hc_max = max(snr_data.snr(snr_data.HD == 0)) %128.1371

hd_mean = mean(snr_data.snr(snr_data.HD == 1)) %71.1200
hd_std = std(snr_data.snr(snr_data.HD == 1)) %28.5838
hd_min = min(snr_data.snr(snr_data.HD == 1)) %18.4088
hd_max = max(snr_data.snr(snr_data.HD == 1)) %126.2594




