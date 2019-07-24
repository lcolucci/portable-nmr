%% Plot MRI pixel-by-pixel histogram of muscular vs subcutaneous tissue
%  This script loads in the MRI pixel-by-pixel bi-exponential fit results
%  and plots a histogram of a single patient's subcu vs muscular tissue
%  values. 
%
%  Code comes from script: plot_pixel_by_pixel_histogram_muscle_vs_subcu.m
%                          (cosmetic differences btw this and original script)
%
%  PREVIOUS SCRIPTS
%       data comes from script: process_mri.m 
%
%  NOTE: takes a while to load the data b/c its large. Caching might be smart. 
%
%  INPUTS
%       whichSubject - which subject's data do you want to plot? 
%       whichAMPM - which time point do you want to plot? 
%       doSave - 1 if you want to save figures
%       savePath - folder in which to save figures
%       pathToDrive - path to CimaLabLargeFiles folder in the drive
%       whichParams - don't need to touch. Which parameters to plot? 
%       histBinEdges - don't need to touch. What binning to use for histogram? (ms)
%
%   OUTPUTS
%       overlaid histogram of muscular and subcutaneous tissue relaxation times 
%       for 1 subject. Histogram saved to savePath (if doSave=1 as poth .eps and .fig 
%       Title of histogram: 
%               MRI Pixel-by-Pixel 2-exp Fit: Histogram of \tau Values in Muscle and Subcu Tissue
%
%
%   L.Colucci 2018

clearvars -except results
close all

%% ----- USER INPUTS ----
whichSubject = 'HC01'; % which subject to plot?
whichAMPM = 'AM'; % which dataset to plot? 
posVector = [560   528   1*600   1*420]; % size of plot (same size as fig. 7 subplots)
fontSize = 26; 
% Save
doSave = 1; % 1 if you want to save figure
savePath = '/Users/linacolucci/Desktop'; % where to save plot? 
% Path
pathToDrive = '/Volumes/CimaLabLargeFiles'; % path to CimaLabLargeFiles folder in the drive
% Histogram Parameters
whichParams ={'relax1','relax2'}; % which parameters to plot in histogram? {'relative_amp2'} or {'relax1','relax2'}
histBinEdges = 1:1:1000; % bin size/length, i.e. x-axis of histograms
%% ---------------------

%% Define colors
keynoteBlue = [0 162 255]./255;
keynoteOrange = [255 147 0]./255;
grey = [138 138 138]./255;
ishaRed = [142 51 49]./255; % 157 43 44  
ishaBlue = [32 25 77]./255; % 45 38 90
ishaGreen = [47 85 59]./255;
ishaTeal = [98 144 128]./255; % 85 165 139
scienceRed = [197, 22, 29]./255;

colorsList = {scienceRed,grey}; %{'b',[.8 .8 .8],[1 .6431 .1098]}; % colors to use in histograms/plots <<<<<<<<<<<<<<<<<<<<<<<<<

%% Load data
if exist('results','var')~=1 %if not yet loaded
    load(fullfile(pathToDrive, 'Hydration/MRI_Study/ProcessedData/mri/pixel_by_pixel/mri_pixelByPixel_2exp.mat'))
    results = result_2exp; 
end

%% Select data (1 subject, 1 AMPM) 
tempResults = results.(whichSubject).(whichAMPM); %structure at 'slice' level
fnSlices = fieldnames(tempResults); 
nSlices = length(fnSlices); 

%% Concatenate data from all slices into 1 long-format table
data = table(); %initialize empty table
for kSlices = 1:nSlices
    tempData = tempResults.(fnSlices{kSlices}); % data from one slice  
    data = vertcat(data, tempData); 
end 

% *** IMPORTANT *** The order of the next 3 sections matters. 
%   1. 'checkfits' must go first because the 99% RMSE cutoff needs to happen with all the dataset for this subject (no subsetting yet. ROI = 'all')
%   2. 'relative amplitudes' goes after 'checkfits' because otherwise you have to add 'relative_amp1=NaN' and 'relative_amp2=NaN' to all of the criteria in checkfits
%   3. 'Subset histogram' could switch places with 'relative amplitudes' but if it goes after then it gives ability to subset data based on relative amplitude values

%% Apply fitting criteria (delete pixel results that don't meet criteria)
data = checkfits(data); %check_1and2exp(data); %   

%% Add 'relative amplitude' components if 2exp fit
if any(strcmp('relax2',data.Properties.VariableNames)) %if it's a 2exp fit, calculate relative amplitudes
    data.relative_amp1 = 100* data.amp1 ./ nansum([data.amp1, data.amp2],2); 
    data.relative_amp2 = 100* data.amp2 ./ nansum([data.amp1, data.amp2],2);
end

%% Apply a subset to select specific tissues
dataMuscle = data(data.muscle_all==1,:); 
dataSubcu = data(data.subcu_all==1,:); 

%% Create vector with the parameter to plot (if more than one param, i.e. {'relax1', 'relax2'}, then vertically concatenate parameter values first)
histogramParam_Muscle = [];  
for kParams = 1:length(whichParams)
    histogramParam_Muscle = vertcat(histogramParam_Muscle, dataMuscle.(whichParams{kParams})); 
end

histogramParam_Subcu = [];  
for kParams = 1:length(whichParams)
    histogramParam_Subcu = vertcat(histogramParam_Subcu, dataSubcu.(whichParams{kParams}));
end

%% Plot
figure
histogram(histogramParam_Muscle, histBinEdges,'FaceColor',colorsList{1}, 'Normalization','probability','EdgeColor','none','FaceAlpha',.9); % Muscular Tissue
hold on; 
histogram(histogramParam_Subcu, histBinEdges,'FaceColor',colorsList{2}, 'Normalization','probability','EdgeColor','none','FaceAlpha',.9); % Subcutaneous Tissue
legend('Muscular Tissue','Subcutaneous Tissue'); legend boxoff 
xlabel('\tau_S and \tau_L (ms)')
ylabel('Relative Probability') 
set(gcf, 'color', 'w','Position',posVector)
set(gca, 'fontsize',fontSize,'linewidth',2); box off % Position',[77, 403, 1400, 550]
%title({'MRI Pixel-by-Pixel 2-exp Fit:','Histogram of \tau Values in Muscle and Subcu Tissue'})
xlim([0 300]) %<<<<<<<<<<<<<<< NOTE: There are some values 500+ms but not many

%% Save 
if doSave==1
    % make folder if it doesn't exist
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end
    
    saveas(gcf, fullfile(savePath, sprintf('fig1C_histogram_muscle_vs_subcu_%s_%s.svg',whichSubject, whichAMPM)))
    %saveas(gcf, fullfile(savePath, sprintf('fig1C_histogram_muscle_vs_subcu_%s_%s.fig',whichSubject, whichAMPM)))
end
