%% Visualize which tissue is changing: Bar Graphs
%   This scripts loads in AM-to-PM change data and plots bar graphs with significance lines/stars. 
%   Bar graph has tissue groups and a separate bar for HC vs HD.
%       

clear

%% --- USER INPUTS ---
whichDataset = 'pixelByPixel_2exp'; %'pixelByPixel_2exp' or 'pixelByPixel1exp' or STILL TO BE implemented: 'roi1exp' or 'roi2exp'
whichParam = 'relative_amp2'; %relax1, relax2, relative_amp2
whichPval = 'unequalVar'; %'unequalVar' or 'equalVar' (which ttest2 p-val do you want to plot?) 
pathToMriStudy = '/Volumes/CimaLabLargeFiles/Hydration'; % Path to folder that contains 'MRI_Study' folder 
whichROIs = {'leg_whole','bone','marrow','subcu_all','muscle_all','muscle_anterior','muscle_deepposterior','muscle_gastrocnemius','muscle_lateral','muscle_soleus'}; 
which1exp = {'bone'}; %if you want an ROI to have 1exp results, write it here. only relax1 result will appear in bar graph
doSave = 0; 
savePath = 'Paper/figures/bargraphs_changing_rois'; 
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

set(f,'Position',[ 73   339   851   528])

% Add significance bar 
if strcmpi(whichPval, 'unequalVar')
    sigstar2(model_series, get(gca,'XTickLabel'), barData.pVals_unequalVar,0) 
end
if strcmpi(whichPval, 'equalVar')
    sigstar2(model_series, get(gca,'XTickLabel'), barData.pVals,0) 
end 

% Save
if doSave==1 
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)   
    end
    saveas(gcf, fullfile(savePath, sprintf('changingROI_%s_%s_%s.eps',whichParam, whichPval, whichDataset)),'epsc')
end