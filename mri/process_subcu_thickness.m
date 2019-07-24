%% Process Subcu Thicknesses
%  This script loads the excel sheet with subcu/skin thickness measurements
%  and processes it into (1) a tidy master datasheet (data) that can be used to
%  explore/plot the data, as well as (2) a tidy summary datasheet (data_summary) 
%  that will be loaded and used in future MRI/LegSensor analyses. 
%
%  Previous Scripts
%   No previous scripts. Data comes from excel workbook: subcu_thickness.xlsx
%   Data was acquired by taking length measurements on each MRI scan using OsiriX Lite
%
%  Inputs
%   pathToWorkbook = path to where the excel sheet is located
%   sheetName = the name of the excel sheet to load in 
%   savePath = the folder where results should be saved to
%   doSave = do you want the results to be saved? (0 | 1) 
%
%  Outputs
%   data = tidy, long-form table with ALL thickness measurements
%   data_summary = table with mean/stdev thickness measurements of skin/subcu/total for each subject
%
%
%  Intermediate Analysis Code (at Bottom of Script)
%   Additional code is included at the bottom of this script to generate
%   the types of figures/intermediate analyses that will be useful for the thesis. 
%

% --- USER INPUTS ---
pathToWorkbook = '/Volumes/CimaLabLargeFiles/Hydration/MRI_Study/Data/subcu_thickness.xlsx'; %'/Users/linacolucci/Desktop/subcu_thickness.xlsx'; 
sheetName ='thicknesses_matlab';
savePath = '/Volumes/CimaLabLargeFiles/Hydration/MRI_Study/ProcessedData'; 
doSave = 1; % 1 or something else
% -------------------

%% Load data
data_raw = load_subcuthickness(pathToWorkbook, sheetName); % <--- If excel sheet changes structure, must update this function 

%% Arrange data into tidy long form (master dataset)
% Re-order the raw data into long-form
    thicknessColumnIDs = find(contains(data_raw.Properties.VariableNames, 'Slice')); %IDs of all the columns with thickness data in them
    nThicknessColStart = min(thicknessColumnIDs); 
    nThicknessColEnd = max(thicknessColumnIDs); 
    data = stack(data_raw, nThicknessColStart:nThicknessColEnd, 'NewDataVariableName','Thickness_mm','IndexVariableName','Slice'); %long-form table
    data.Slice = cellstr(data.Slice); %convert to string cell instead of categorical 
% Convert 'Slice' column into 2 separate columns (i.e. 'Slice1_a' into '1' and 'a')
    sliceColumnSplit = regexp(data.Slice, '_','split'); %split 'Slice1' and 'b' into different cells
    temp = vertcat(sliceColumnSplit{:}); % put cells into a matrix with 2 columns
    temp = table(temp(:,1),temp(:,2),'VariableNames',{'Slice','Repetition'}); 
    temp.Slice = str2double(erase(temp.Slice,'Slice')); %remove 'Slice' and convert to numeric values
% Now attach the 'temp' table with the 2 new columns into the 'data' table
    data.Slice = []; %delete column
    data = horzcat(data(:,1:end-1), temp,data(:,end)); % <--- this is the tidy, master dataset
% Clear intermediate variables that are no longer needed
    clear temp sliceColumnSplit nThickness* thicknessColumnIDs

    
%% Summary Table: Average/Stdev of Thicknesses (master table): average all thicknesses across AM and PM for each subject
% Calculate mean and stdev for each subject and tissue type
    data_summary = grpstats(data, {'Subject_Name','Tissue'},{'mean','std'},'DataVars','Thickness_mm'); 
% Temporary table that calculates total thickness as sum of skin + subcu thicknesses
    tempTotalThick = varfun(@sum, data_summary, 'InputVariables',{'mean_Thickness_mm','std_Thickness_mm'},'GroupingVariables','Subject_Name');
    tempTotalThick.Tissue(:) = categorical({'total'}); 
    tempTotalThick.Properties.VariableNames{3} = 'mean_Thickness_mm';
    tempTotalThick.Properties.VariableNames{4} = 'std_Thickness_mm';
    tempTotalThick = [tempTotalThick(:,1), tempTotalThick(:,end), tempTotalThick(:,(2:end-1))];
% Add tempTotalThick table results into the master summary table
    data_summary = vertcat(data_summary, tempTotalThick);
    data_summary = sortrows(data_summary,'Subject_Name');
    data_summary.Properties.RowNames = {};
    data_summary.GroupCount = []; 
% Clear intermediate variables that are no longer needed 
    clear tempTotalThick 

%% Save Results
if doSave==1
    save(fullfile(savePath, 'thicknesses_allmeasurements.mat'), 'data')
    save(fullfile(savePath, 'thicknesses_summary_mean_std.mat'), 'data_summary')
end

%% ===============================================================================
%% ---- Code for Intermediate Analyses that will be useful for thesis figures ----

% Table: Calculate separate AM/PM thickness means and stdevs for each subject
ampm = grpstats(data, {'Subject_Name','AMPM','Tissue'},{'mean','std'},'DataVars','Thickness_mm'); %<-- better than varfun b/c it calculates multiples stats at once
ampm_skin = ampm(ampm.Tissue =='skin',:); % 
ampm_subcu = ampm(ampm.Tissue =='subcu',:); 

% Boxplot of Thickness Measurements
tempData = data(data.Tissue == 'subcu',:); % USER: pick whatever selection criteria for the tempData
figure; boxplot(tempData.Thickness_mm, {tempData.Subject_Name, tempData.AMPM}) %<-- can include more than 1 variable (i.e. more than just 'Subject_Name') in 2nd entry of boxplot
title('Boxplot of Subcu Thicknesses for Each Subject')
ylabel('Thickness (mm)')
set(gcf, 'color', 'w');set(gca, 'FontSize', 16);

% Table: Calculate AM-PM Differences to see if any major discrepancies exist
ampm_differences = varfun(@(x) x(1)-x(2), ampm, 'InputVariables',{'mean_Thickness_mm'},'GroupingVariables',{'Subject_Name','Tissue'});

% Table: Perform ttest to see if differences btw AM-PM thickness values are statistically significant or not
ampm_differences_ttest = varfun(@(x) ttest(x(1:length(x)/2),x(length(x)/2+1:end)), data, 'InputVariables',{'Thickness_mm'},'GroupingVariables',{'Subject_Name','Tissue'});






