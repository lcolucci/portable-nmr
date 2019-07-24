%% Process MRI ROI Average: Which Tissue is Changing? 

clear; close all


%% ------- USER INPUTS -------
% - - - Which Fitting Results?
whichData = 'result_2exp'; % 'result_1exp' or 'result_2exp'
whichROIs = {'leg_whole','bone_tibia','marrow_tibia','subcu_all','muscle_all','muscle_anterior','muscle_deepposterior','muscle_gastrocnemius','muscle_lateral','muscle_soleus'}; 
% - - -
% Paths to data 
pathToMriStudy = '/Volumes/CimaLabLargeFiles/Hydration'; % Path to folder that contains 'MRI_Study' folder (contains raw images and fit results, mri_pixelByPixel_1exp.mat, etc.)
% - - - Save
doSave = 0; % 1 to save
savePath = '/Volumes/CimaLabLargeFiles/Hydration/MRI_Study/ProcessedData/mri/ROIs'; 
%% -------------------------------

%% Load data
if strmatch(whichData, 'result_2exp')
    load(fullfile(pathToMriStudy, 'MRI_Study/ProcessedData/mri/ROIs/ROIs_2expFits.mat'))
    data = results_2exp; 
end

%% Process and re-format
masterDifferenceValue = table(); 
masterDifferenceStd = table(); 
masterTtest = table(); 
roiList = fieldnames(data); 
nROIs = length(roiList); 
for iROI = 1:nROIs
    tempData = data.(roiList{iROI});
    
    % Subtract AM and PM
    difference = varfun(@diff, tempData,'InputVariables',{'relax1','relax2','relative_amp2'}, 'GroupingVariables',{'Subject','HD','ROI'});
    
    % Mean/Stdev of HC and HD
    differenceValue = varfun(@mean, difference, 'InputVariables',{'diff_relax1','diff_relax2','diff_relative_amp2'},'GroupingVariables',{'HD','ROI'});
    differenceStd= varfun(@std, difference, 'InputVariables',{'diff_relax1','diff_relax2','diff_relative_amp2'},'GroupingVariables',{'HD','ROI'});
    %differenceValueStd = join(differenceValue, differenceStd); 
    
    % t-test HC vs HD
    try 
        temp_ttest = table(); 
        [tstat, temp_ttest.pVals_relax1] = ttest2(difference.diff_relax1(difference.HD==0),difference.diff_relax1(difference.HD==1));
        [tstat, temp_ttest.pVals_relax2] = ttest2(difference.diff_relax2(difference.HD==0),difference.diff_relax2(difference.HD==1));
        [tstat, temp_ttest.pVals_relative_amp2] = ttest2(difference.diff_relative_amp2(difference.HD==0),difference.diff_relative_amp2(difference.HD==1));
        temp_ttest.ROI = categorical(roiList(iROI));
    catch
    end
    
    % Aggregate into master
    masterDifferenceValue = vertcat(masterDifferenceValue, differenceValue); 
    masterDifferenceStd = vertcat(masterDifferenceStd, differenceStd); 
    masterTtest = vertcat(masterTtest, temp_ttest); 
    
    % Delete temp vars
    clear tempData difference differenceValueStd temp_ttest differenceValue differenceStd
    
end

% Delete 'Group Count' column
masterDifferenceValue.GroupCount = []; %delete group count column
masterDifferenceStd.GroupCount = []; %delete group count column

% Tidy data
masterDifferenceValue = tidydata(masterDifferenceValue, 'Mean', {'mean_diff_relax1','mean_diff_relax2','mean_diff_relative_amp2'}, {'relax1','relax2','relative_amp2'}); 
masterDifferenceStd = tidydata(masterDifferenceStd, 'Std', {'std_diff_relax1','std_diff_relax2','std_diff_relative_amp2'}, {'relax1','relax2','relative_amp2'}); 

% Merge Mean and Std into Master
masterDifference = join(masterDifferenceValue, masterDifferenceStd); 

% Merge with p-values (t-test results)
masterTtest2 = stackdata(masterTtest, 'pVals', {'pVals_relax1','pVals_relax2','pVals_relative_amp2'}, {'relax1','relax2','relative_amp2'});

% Merge MeanStd with pVals
master = join(masterTtest2,masterDifference,'Keys',{'ROI','Parameter'}); 


%% Save results
if doSave==1 
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end
    save(fullfile(savePath, 'mri_roi_2exp_summary_AMPMchange.mat'),'master')
end


%% Delete unnecessary ROIs for paper figure (phantoms, voxels, etc.) 
master(master.ROI=='all'|master.ROI=='landmark'|master.ROI=='legSensorROI'|master.ROI=='muscle_voxel'|...
    master.ROI=='muscle_voxel_15'|master.ROI=='phantom_long'|master.ROI=='phantom_short'|master.ROI=='phantom_oil'|...
    master.ROI=='subcu_voxel'|master.ROI=='voxel_anterior1'|master.ROI=='voxel_anterior2'|master.ROI=='voxel_lateral1'|...
    master.ROI=='voxel_lateral2'|master.ROI=='marrow_fibia',:)=[]; 


%% Bar Graphs
whichParam = 'relax1'; %relax1, relax2, relative_amp2
barData = master(master.Parameter==whichParam,:); 


% Data to be plotted as a bar graph
model_series = [barData.Mean_HC, barData.Mean_HD];
%Data to be plotted as the error bars
model_error =  [barData.Std_HC, barData.Std_HD];
% Other inputs to bar graph
names = cellstr(barData.ROI); 
xLab = ''; 
yLab = 'Change (%)';
legendLabs = {'HC','HD'};

% Plot
bar_with_error(model_series, model_error, names, xLab, yLab, legendLabs)
title(sprintf('%s',whichParam),'Interpreter','none')



%% FUNCTIONS

% tidydata 
%       performs 'stackdata', then unstacks so that HC and HD have separate columns
function unstackedData = tidydata(data, valuename, oldnames, newnames)    
    stackeddata = stackdata(data, valuename, oldnames, newnames);
    
    unstackedData = unstack(stackeddata, valuename,'HD');
    unstackedData.Properties.VariableNames = {'ROI','Parameter',strcat(valuename,'_HC'),strcat(valuename,'_HD')}; 
end

% stackdata
%       stacks data from horizontal form to long form
function stackeddata = stackdata(data, valuename, oldnames, newnames)
    stackeddata = stack(data,oldnames);
    try 
        stackeddata.Properties.VariableNames = {'HD','ROI','Parameter',valuename}; 
    catch
        stackeddata.Properties.VariableNames = {'ROI','Parameter',valuename}; 
    end
    stackeddata.Parameter = renamecats(stackeddata.Parameter, oldnames,newnames);
end
