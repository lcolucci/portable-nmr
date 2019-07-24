%% Process ROI_means into HC vs HD means
%  This script loads ROImeans_table and averages decays together to produce
%  an average HC AM, HC PM, HD AM, and HD PM decay for each ROI. 
%
%  PREVIOUS SCRIPTS
%       data comes from script: process_mri.m (for MRI) OR process_t2decays_table.m (for MR Sensor)
%
% INPUTS
%   load ROImeans_table.mat
%   doSave - 1 to save 'results' output
%   savePath - file path for where to save 'results' 
%
% OUTPUTS
%   results - table with average T2 decay for HC, HD, AM, PM for each ROI 
%
% Lina Colucci 2018 

clear
%% ---- USER INPUTS -----
load('/Volumes/CimaLabLargeFiles/Hydration/MRI_Study/ProcessedData/mr_sensor_leg/legSensor_T2Decays_table.mat')
data = data_table; % change if loaded data has a different name
% Save
doSave = 1; % 1 to save
savePath = '/Volumes/CimaLabLargeFiles/Hydration/MRI_Study/ProcessedData/mr_sensor_leg/Groups_HCvsHD';% '/Volumes/CimaLabLargeFiles/Hydration/MRI_Study/ProcessedData/mri/ROIs'; 
%% ----------------------

%% Initialize empty results table
results = table(); 

%% HC AM
% Select correct subset of data
try % MRI data
    tempData1 = data((data.AMPM == 'AM'|data.AMPM=='am') & data.AllSlices==1 & contains(cellstr(data.Subject),'HC'), :);
catch % MR Sensor data
    tempData1 = data((data.AMPM == 'AM'|data.AMPM=='am') & contains(cellstr(data.Subject),'HC'), :);
end
% Loop 
ROIlist = unique(tempData1.ROI); 
nROI = length(ROIlist); 
lastRow = height(results); 
for iROI = 1:nROI 
    tempData2 = tempData1(tempData1.ROI == ROIlist(iROI),:); 
    
    % Average T2 decays from each subject together
    allDecays = []; 
    for jSubjects = 1:height(tempData2)
        singleDecay = tempData2.T2Decay{jSubjects,1}(:); 
        if size(singleDecay,1) > size(singleDecay,2)
            singleDecay = singleDecay';
        end
        if ~isempty(allDecays) && size(singleDecay,2) ~= size(allDecays,2) 
            singleDecay(size(allDecays,2)+1:end)=[]; 
        end
        allDecays = vertcat(allDecays, singleDecay); 
    end
    avgAllDecays = nanmean(allDecays,1); 
    
    % Aggregate results
    results.Subject(lastRow + iROI) = categorical({'HC'}); 
    results.AMPM(lastRow + iROI) = categorical({'AM'}); 
    results.ROI(lastRow + iROI) = ROIlist(iROI); 
    results.T2Decay(lastRow + iROI) = {avgAllDecays(:)}; 
    try results.Time(lastRow + iROI) = tempData2.Time(jSubjects-1); catch, end
    
    % Delete temp variables
    clear tempData2 allDecays avgAllDecays
end
clear tempData1

%% HC PM
% Select correct subset of data
try % MRI data
    tempData1 = data((data.AMPM == 'PM'|data.AMPM=='pm') & data.AllSlices==1 & contains(cellstr(data.Subject),'HC'), :);
catch % MR Sensor data
    tempData1 = data((data.AMPM == 'PM'|data.AMPM=='pm') & contains(cellstr(data.Subject),'HC'), :);
end
% Loop 
ROIlist = unique(tempData1.ROI); 
nROI = length(ROIlist); 
lastRow = height(results); 
for iROI = 1:nROI 
    tempData2 = tempData1(tempData1.ROI == ROIlist(iROI),:); 
    
    % Average T2 decays from each subject together
    allDecays = []; 
    for jSubjects = 1:height(tempData2)
        singleDecay = tempData2.T2Decay{jSubjects,1}(:); 
        if size(singleDecay,1) > size(singleDecay,2)
            singleDecay = singleDecay';
        end
        if ~isempty(allDecays) && size(singleDecay,2) ~= size(allDecays,2) 
            singleDecay(size(allDecays,2)+1:end)=[]; 
        end
        allDecays = vertcat(allDecays, singleDecay); 
    end
    avgAllDecays = nanmean(allDecays,1); 
    
    % Aggregate results
    results.Subject(lastRow + iROI) = categorical({'HC'}); 
    results.AMPM(lastRow + iROI) = categorical({'PM'}); 
    results.ROI(lastRow + iROI) = ROIlist(iROI); 
    results.T2Decay(lastRow + iROI) = {avgAllDecays(:)}; 
    try results.Time(lastRow + iROI) = tempData2.Time(jSubjects-1); catch, end
    
    % Delete temp variables
    clear tempData2 allDecays avgAllDecays
end
clear tempData1

%% HD AM
% Select correct subset of data
try % MRI data
    tempData1 = data((data.AMPM == 'AM'|data.AMPM=='am') & data.AllSlices==1 & contains(cellstr(data.Subject),'HD'), :);
catch % MR Sensor data
    tempData1 = data((data.AMPM == 'AM'|data.AMPM=='am') & contains(cellstr(data.Subject),'HD'), :);
end% Loop 
ROIlist = unique(tempData1.ROI); 
nROI = length(ROIlist); 
lastRow = height(results); 
for iROI = 1:nROI 
    tempData2 = tempData1(tempData1.ROI == ROIlist(iROI),:); 
    
    % Average T2 decays from each subject together
    allDecays = []; 
    for jSubjects = 1:height(tempData2)
        singleDecay = tempData2.T2Decay{jSubjects,1}(:); 
        
        %% Special exceptions of MR sensor data collected with wrong TE
        if tempData2.ROI(jSubjects) == 'legSensor' && strcmpi(cellstr(tempData2.Subject(jSubjects)),'HDMRI01') 
            % I'm using row # 7 (HDMRI05) as the one with correct time array. Change if this is no longer true. 
            modifiedDecay = nan(size(tempData2.Time{7,1}));
            [intersection, ia, ib] = intersect(tempData2.Time{jSubjects,1} , tempData2.Time{7,1}  );
            modifiedDecay(ib) = tempData2.T2Decay{jSubjects,1}(ia);
            singleDecay = modifiedDecay; 
        end
        if tempData2.ROI(jSubjects) == 'legSensor' && strcmpi(cellstr(tempData2.Subject(jSubjects)),'HDMRI02') 
            % I'm using row # 7 (HDMRI05) as the one with correct time array. Change if this is no longer true. 
            modifiedDecay = nan(size(tempData2.Time{7,1}));
            [intersection, ia, ib] = intersect(tempData2.Time{jSubjects,1} , tempData2.Time{7,1}  );
            modifiedDecay(ib) = tempData2.T2Decay{jSubjects,1}(ia);
            singleDecay = modifiedDecay; 
        end
        if tempData2.ROI(jSubjects) == 'legSensor' && strcmpi(cellstr(tempData2.Subject(jSubjects)),'HDMRI03') 
            % Need to add NaNs to the end
            nNaNs = length(allDecays) - length(singleDecay);
            singleDecay = [singleDecay;nan(nNaNs,1)]; 
        end
        
        %% Decay is a row vector. If column, flip it. 
        if size(singleDecay,1) > size(singleDecay,2)
            singleDecay = singleDecay';
        end
        %% If singleDecay for this subject is shorter than the other ones, cut off the end of the singleDecay to make it the same length. 
        if ~isempty(allDecays) && size(singleDecay,2) ~= size(allDecays,2) 
            singleDecay(size(allDecays,2)+1:end)=[]; 
        end
        allDecays = vertcat(allDecays, singleDecay); 
    end
    avgAllDecays = nanmean(allDecays,1); 
    
    % Aggregate results
    results.Subject(lastRow + iROI) = categorical({'HD'}); 
    results.AMPM(lastRow + iROI) = categorical({'AM'}); 
    results.ROI(lastRow + iROI) = ROIlist(iROI); 
    results.T2Decay(lastRow + iROI) = {avgAllDecays(:)}; 
    try results.Time(lastRow + iROI) = tempData2.Time(jSubjects-1); catch, end
    
    % Delete temp variables
    clear tempData2 allDecays avgAllDecays
end
clear tempData1

%% HD PM
% Select correct subset of data
try % MRI data
    tempData1 = data((data.AMPM == 'PM'|data.AMPM=='pm') & data.AllSlices==1 & contains(cellstr(data.Subject),'HD'), :);
catch % MR Sensor data
    tempData1 = data((data.AMPM == 'PM'|data.AMPM=='pm') & contains(cellstr(data.Subject),'HD'), :);
end
% Loop
ROIlist = unique(tempData1.ROI); 
nROI = length(ROIlist); 
lastRow = height(results); 
for iROI = 1:nROI 
    tempData2 = tempData1(tempData1.ROI == ROIlist(iROI),:); 
    
    % Average T2 decays from each subject together
    allDecays = []; 
    for jSubjects = 1:height(tempData2)
        singleDecay = tempData2.T2Decay{jSubjects,1}(:); 
        
        %% Special exceptions of MR sensor data collected with wrong TE
        if tempData2.ROI(jSubjects) == 'legSensor' && strcmpi(cellstr(tempData2.Subject(jSubjects)),'HDMRI01') 
            % I'm using row # 7 (HDMRI05) as the one with correct time array. Change if this is no longer true. 
            modifiedDecay = nan(size(tempData2.Time{7,1}));
            [intersection, ia, ib] = intersect(tempData2.Time{jSubjects,1} , tempData2.Time{7,1}  );
            modifiedDecay(ib) = tempData2.T2Decay{jSubjects,1}(ia);
            singleDecay = modifiedDecay; 
        end
        if tempData2.ROI(jSubjects) == 'legSensor' && strcmpi(cellstr(tempData2.Subject(jSubjects)),'HDMRI02') 
            % I'm using row # 7 (HDMRI05) as the one with correct time array. Change if this is no longer true. 
            modifiedDecay = nan(size(tempData2.Time{7,1}));
            [intersection, ia, ib] = intersect(tempData2.Time{jSubjects,1} , tempData2.Time{7,1}  );
            modifiedDecay(ib) = tempData2.T2Decay{jSubjects,1}(ia);
            singleDecay = modifiedDecay; 
        end
        if tempData2.ROI(jSubjects) == 'legSensor' && strcmpi(cellstr(tempData2.Subject(jSubjects)),'HDMRI03') 
            % Need to add NaNs to the end
            nNaNs = length(allDecays) - length(singleDecay);
            singleDecay = [singleDecay;nan(nNaNs,1)]; 
        end
        
        if size(singleDecay,1) > size(singleDecay,2)
            singleDecay = singleDecay';
        end
        if ~isempty(allDecays) && size(singleDecay,2) ~= size(allDecays,2) 
            singleDecay(size(allDecays,2)+1:end)=[]; 
        end
        allDecays = vertcat(allDecays, singleDecay); 
    end
    avgAllDecays = nanmean(allDecays,1); 
    
    % Aggregate results
    results.Subject(lastRow + iROI) = categorical({'HD'}); 
    results.AMPM(lastRow + iROI) = categorical({'PM'}); 
    results.ROI(lastRow + iROI) = ROIlist(iROI); 
    results.T2Decay(lastRow + iROI) = {avgAllDecays(:)}; 
    try results.Time(lastRow + iROI) = tempData2.Time(jSubjects-1); catch, end
    
    % Delete temp variables
    clear tempData2 allDecays avgAllDecays
end
clear tempData1

%% Save
if doSave==1
    save(fullfile(savePath,'ROImeans_HCHD_AMPM.mat'),'results')
end
