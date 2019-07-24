%% ----------------------------------- Fig 1D-F -------------------------------
%% MRI pixel-wise data  |  AM-PM Change in each ROI (Bar Graphs)  |  Comparison of HC vs HD 
% 
%   This script generates one figure containing many bar graphs for a particular indicator (i.e. R_L, tau_S, etc.) (user picks indicator). 
%   There are up to 10 ROIs (user inputs) and one set of bar graphs for each ROI.
%   One set has 2 bars: one for the AM-PM change for HC subjects, the other for the AM-PM change for HD subjects
%   There is a significance star above the HC and HD bars. The user can pick which statistical test to use. 
%   
%   To generate the figure for the paper, run this script 3 times each time
%   for a different indicator (i.e. relax1, relax2, relative_amp2). 
%
%   PREVIOUS SCRIPTS
%       data comes from script: process_pixelbypixel_averages.m
%       real statistics comes from R script: statistics_for_paper.R
%
%   INPUTS
%       whichDataset - do you want to plot the pixel-wise 2exp fit or pixel-wise 1exp fit dataset? 'pixelByPixel_2exp' or 'pixelByPixel1exp'
%                      data comes from script: process_pixelbypixel_averages.m
%       whichParam - which indicator do you want to plot? relax1, relax2, relative_amp2
%       whichPval - which statistical test to compare HC and HD? 
%                   equalVar = ttest2 two-sided with equal variances.
%                   unequalVar = ttest2 two-sided with unequal variances.
%                   permTest = two-sample, two-sided exact permutation test with 1e5-1 repetitions
%       pathToMriStudy = path to folder that contains 'MRI_Study' folder, which contains all the data
%       whichROIs = specify which ROIs to visualize in bar graph
%       which1exp = specify which ROIs, if any, if you want to use the 1-exp results instead of 2-exp. If any ROIs are specified here,
%                   it'll only appear in the relax1 plots. It will be empty in any of the other plots (i.e relax2 or relative_amp2)
%       doSave = 1 to save plots as .eps
%       savePath = if doSave=1, specify path to which figures will be saved
%
%   OUTPUTS
%       eps figure gets saved in the savePath folder
%       'data' (i.e. the table that is loaded) is saved as a csv file. This table goes into supplemental materials. 
%
%   Lina A. Colucci, 2018

clear

%% --- USER INPUTS ---
whichDataset = 'pixelByPixel_2exp'; %'pixelByPixel_2exp' or 'pixelByPixel1exp' or STILL TO BE implemented: 'roi1exp' or 'roi2exp'
whichParam = 'relax2'; %relax1, relax2, relative_amp2
whichPval = 'unequalVar'; % which p-val do you want to plot? 'permTest', 'unequalVar' or 'equalVar'
                        % equalVar - two-sample, two-tailed t-test with equal variances btw HC and HD subjects. ttest2(Mean_HC, Mean_HD)
                        % unequalVar - two-sample, two-tailed t-test with unequal variances btw HC and HD subjects. ttest2(Mean_HC, Mean_HD,'Vartype','unequal')
                        % permTest - two-sample, two-tailed exact permutation test with 10e5-1 repetitions. permutationTest(Mean_HC, Mean_HD, 1e5-1, 'exact',1)
pathToMriStudy = '/Volumes/CimaLabLargeFiles/Hydration'; % Path to folder that contains 'MRI_Study' folder 
whichROIs = {'leg_whole','bone','marrow','subcu_all','muscle_all','muscle_anterior','muscle_deepposterior','muscle_gastrocnemius','muscle_lateral','muscle_soleus'}; 
which1exp = {'bone'}; %if you want an ROI to have 1exp results, write it here. only relax1 result will appear in bar graph
% ---
doSave = 1; 
savePath = '/Users/linacolucci/Documents/Thesis/Paper/figures/MRI_pixelwise'; 
%% --------------------

%% Load Dataset
if strmatch(whichDataset,'pixelByPixel_2exp')
    load(fullfile(pathToMriStudy, 'MRI_Study/ProcessedData/mri/pixel_by_pixel/mri_pixelByPixel_2exp_summary_AMPMchange.mat'))
    data = master; 
end
if strmatch(whichDataset, 'roi2exp')
    load(fullfile(pathToMriStudy,'MRI_Study/ProcessedData/mri/ROIs/mri_roi_2exp_summary_AMPMchange.mat'))
    data = master; 
end
if strmatch(whichDataset,'pixelByPixel_1exp')
    load(fullfile(pathToMriStudy, 'MRI_Study/ProcessedData/mri/pixel_by_pixel/mri_pixelByPixel_1exp_summary_AMPMchange.mat'))
    data = master; 
end
if ~isempty(which1exp)
    load(fullfile(pathToMriStudy, 'MRI_Study/ProcessedData/mri/pixel_by_pixel/mri_pixelByPixel_1exp_summary_AMPMchange.mat'))
    data1exp = master; 
end

%% Subset parameter and ROIs 
barData = data(data.Parameter==whichParam,:); 
barData = barData(ismember(unique(barData.ROI), categorical(whichROIs)), :); 
barData = sortrows(barData, 'ROI'); 

% Sort table in specified ROI order 
% [a, newOrder] = sort(categorical(whichROIs)); 
% test = barData(newOrder,:); %% sorting tabel in order specified by 'whichROIs' is not working perfectly. muscle_all in wrong place. Don't know why. 

if ~isempty(which1exp) %if params other than relax1, NaN. 
    n1exp = length(which1exp); 
    for i = 1:n1exp
        whichRow2exp = find(barData.ROI == which1exp{i}); 
        whichRow1exp = find(data1exp.ROI == which1exp{i}); 
        if strcmpi(whichParam,'relax1') 
            barData(whichRow2exp,:) = data1exp(whichRow1exp,:); 
        else
            barData(whichRow2exp,3:end) = {NaN}; 
        end
    end
end

% Data to be plotted as a bar graph
model_series = [barData.Mean_HC, barData.Mean_HD];
model_series = -1 .* model_series; %decrease should be negative
%Data to be plotted as the error bars
model_error =  [barData.Std_HC, barData.Std_HD];
% Other inputs to bar graph
names = cellstr(barData.ROI); 
xLab = ''; 
if contains(whichParam, 'relax')
    yLab = 'Change (ms)'; 
end
if contains(whichParam, 'amp')
    yLab = 'Change (%)';
end
legendLabs = {'HC','HD'};

% Plot
[f, h] = bar_with_error(model_series, model_error, names, xLab, yLab, legendLabs); 
title(sprintf('%s',whichParam),'Interpreter','none')

set(f,'Position',[ 73   255   856   612])

% Add significance bar 
if strcmpi(whichPval, 'unequalVar')
    sigstar2(model_series, get(gca,'XTickLabel'), barData.pVals_unequalVar,0,'bar') 
end
if strcmpi(whichPval, 'equalVar')
    sigstar2(model_series, get(gca,'XTickLabel'), barData.pVals,0,'bar') 
end
if strcmpi(whichPval, 'permTest')
    sigstar2(model_series, get(gca,'XTickLabel'), barData.pVals_permTest,0,'bar') 
end

% Save
if doSave==1 
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)   
    end
    % save figure
    saveas(gcf, fullfile(savePath, sprintf('mri_pixelwise_lastIntegralDiff_ROIs_%s_%s_%s.eps',whichParam, whichPval, whichDataset)),'epsc') 
    % save table as a csv (for supplemental materials) 
    writetable(data, fullfile(savePath, 'mri_pixelwise_lastIntegralDiff_HCvsHD.csv'))
end

