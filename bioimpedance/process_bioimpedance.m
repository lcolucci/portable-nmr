%% Bioimpedance - Clean & Tidy Master Dataset 
%  This script takes raw BI data and turns it into a tidy master BI table
%
%  INPUTS:  savePath
%           add paths to new BI data as it becomes available
%           path to load demogrpahics data (which contains leg circumferences, etc.)
%           custom edits to timePt labels: the "Tidy the Data" section, may need new code if the new datasets have timePt labels that aren't described in code 
%  OUTPUTS: 2 figures consisting of 4 subplots (Whole Body vs Leg, TBW vs ECF)
%           BI_AllIndividualTrials.mat - tidy table of BI data for each individual trial
%           BI_Averages.mat - tidy table of BI data where trials at the same time pt are averaged together
%           BI_BadTrials.mat - tidy table showing which BI trials were deemed 'bad' and not included in the two other tidy tables
%
%  Rzero = Re = Resistivity associated with ECF (lower frequencies don't pass through cells)
%  Rinf = Resistivity associated with TBW (higher frequencies pass through everything)
%  Leg = electrodes are on lateral side of right calf
%  Whole Body = electrodes are on right hand to right foot
%
%  Sub-function calculateLegVolumesBI.m 
%           Re-calculates the leg BI volumes with proper equations 
%
%  2018-01-22 Lina Colucci 
%   

clear

%% Set Save Path
savePath = '/Volumes/cimalablargefiles/Hydration/MRI_Study/ProcessedData/bioimpedance'; 

%% Load paths to .xlsx BI data
basepath = '/Volumes/CimaLabLargeFiles/Hydration/'; % path to folder that contains 'MRI_Study' folder
BIpaths.HC01 = fullfile(basepath, 'MRI_Study/Data/HC001/Bioimpedance/batch.xlsx'); 
BIpaths.HC01b = fullfile(basepath, 'MRI_Study/Data/HC01b/Bioimpedance/hc01b_batch.xlsx'); 
BIpaths.HC02 = fullfile(basepath, 'MRI_Study/Data/HC002/Bioimpedance/hc002_batch.xlsx'); 
BIpaths.HC03 = fullfile(basepath, 'MRI_Study/Data/HC003/Bioimpedance/HC003_batch.xlsx'); 
BIpaths.HC04 = fullfile(basepath, 'MRI_Study/Data/HC004/Bioimpedance/hc004_batch.xlsx'); 
BIpaths.HC05 = fullfile(basepath, 'MRI_Study/Data/HC005/Bioimpedance/hc005_batch.xlsx'); 
BIpaths.HC06 = fullfile(basepath, 'MRI_Study/Data/HC06/Bioimpedance/hc06_batch.xlsx'); 
BIpaths.HDMRI01 =  fullfile(basepath, 'MRI_Study/Data/HDMRI01/Bioimpedance/hdmri01_batch.xlsx'); 
BIpaths.HDMRI01b = fullfile(basepath, 'MRI_Study/Data/HDMRI01b/Bioimpedance/batch_hd01b.xlsx'); 
BIpaths.HDMRI02 = fullfile(basepath, 'MRI_Study/Data/HDMRI02/Bioimpedance/hdmri02_batch.xlsx'); 
BIpaths.HDMRI02b = fullfile(basepath, 'MRI_Study/Data/HDMRI02b/Bioimpedance/hdmri02b_batch.xlsx'); 
BIpaths.HDMRI03 = fullfile(basepath, 'MRI_Study/Data/HDMRI03/Bioimpedance/hdmri03_batch.xlsx'); 
BIpaths.HDMRI04b =  fullfile(basepath, 'MRI_Study/Data/HDMRI04b/Bioimpedance/hdmri04b_batch.xlsx'); 
BIpaths.HDMRI05 = fullfile(basepath, 'MRI_Study/Data/HDMRI05/Bioimpedance/hdmri05_batch.xlsx'); 

%% Load Demographics (need for calculating leg segmental volumes) 
load('/Volumes/cimalablargefiles/Hydration/MRI_Study/ProcessedData/redcap/demographics.mat');

%% Time Point Numbers
% If additional ticks are added here, make sure to add them to uniqueSortedTimePts (search for *** symbols)
preTick = -1; connectTick = -0.5; disconnectTick = 98; postTick = 99; continuousTick = 90; bolusTick = 80; % USER 

%% Lookup table
lookup = table({'HC01','HC01b','HC02','HC03','HC04','HC05','HC06','HDMRI01','HDMRI01b','HDMRI02','HDMRI02b','HDMRI03','HDMRI04b','HDMRI05'}',...
                [101, 101, 102, 103, 104, 105,106, 1,1,2,2,3,40,5]',...
                [1,2,1,1,1,1,1,1,2,1,2,1,1,1]',...
                [0,0,0,0,0,0,0,1,1,1,1,1,1,1]'); 
lookup.Properties.VariableNames = {'SubjectName','Subject','VisitNum','HD'}; 


%% ---------------------- AUTOMATED CODE (NO USER INPUTS NEEDED) -----------------------------------
% Loop
fieldnames_subjects = fieldnames(BIpaths); 
masterAverages = table();
badTrials = table(); 
masterIndividualTrials = table(); 
for ii= 1:length(fieldnames_subjects)
    
    %% Load raw data from excel files
    data = importBIdata(BIpaths.(fieldnames_subjects{ii})); 
    
    %% Custom edits to the raw data (doing this here so I don't have to make change to the original excel files)
    if strcmp(fieldnames_subjects(ii),'HDMRI05')
        data.Filename = strrep(data.Filename,'mri05-pst-lg-0008.mfu','mri05-pst-wb-0008.mfu'); %These values are wb measurements, not leg ones
        data.Filename = strrep(data.Filename,'mri05-pst-lg-0009.mfu','mri05-pst-wb-0009.mfu'); % " " 
        data.Filename = strrep(data.Filename,'mri05-pst-lg-0010.mfu','mri05-pst-wb-0010.mfu'); % " "
        data.Filename = strrep(data.Filename,'mri-05-','mri05-'); %remove dash
    end
        
    %% Tidy the data
    data.Gender = []; % delete Gender b/c it's the only 'text' column, it causes issues, and it's unncessary since all subjects are Male 
    data.SubjectName(:) = fieldnames_subjects(ii); % Add subject name (for lookup table reference)
    data.Date = datetime(data.Date); % Convert Time column to time format (instead of text)
    data.WholeBody = contains(data.Filename,'wb'); % whole body = 1. leg measurement = 0. 
    for jj=1:height(data)% Time Point column
        tempFilenameAsSeparatedCells = strsplit(data.Filename{jj},'-'); %tempFilenameAsSeparatedCells is a 1x4 cell: subject, timept, wb_lg, trial
        
        %Pre
        if contains(tempFilenameAsSeparatedCells{2}, 'pre','IgnoreCase',1)
            data.TimePt(jj) = preTick; 
        % Post
        elseif contains(tempFilenameAsSeparatedCells{2}, 'post','IgnoreCase',1) | contains(tempFilenameAsSeparatedCells{2}, 'end','IgnoreCase',1) | contains(tempFilenameAsSeparatedCells{2}, 'pst','IgnoreCase',1)
            data.TimePt(jj) = postTick; 
        else
            data.TimePt(jj) = NaN; 
        end
        
        % Trial number within TimePt
        data.Trial(jj) = str2double(strrep(tempFilenameAsSeparatedCells{4},'.mfu','')); % Remove '.mfu' and convert to number

        clear tempString
    end
    
    % Merge with lookup table (to create subject #, visit #, and HD columns) 
    data = join(data, lookup);
    
    % Re-organize columns so identifying info appears first
    colID_SubjectName = find(strcmp(data.Properties.VariableNames,'SubjectName')); 
    colID_Subject = find(strcmp(data.Properties.VariableNames,'Subject')); 
    colID_Trial = find(strcmp(data.Properties.VariableNames,'Trial')); 
    data = [data(:,colID_SubjectName), data(:,[colID_Subject:end]), data(:,[colID_SubjectName+1:colID_Trial]),data(:,1:colID_SubjectName-1)]; 
    
    % Re-organize so that columns that don't need stdev appear at the end 
    colID_Flow = find(strcmp(data.Properties.VariableNames,'Flow')); %from Flow to Reject
    colID_Reject = find(strcmp(data.Properties.VariableNames,'Reject'));%from Flow to Reject
    colID_Td = find(strcmp(data.Properties.VariableNames,'Td'));%keep Td column
    colID_Total = find(strcmp(data.Properties.VariableNames,'Total'));%from Total to Ignored
    colID_Ignored = find(strcmp(data.Properties.VariableNames,'Ignored'));%from Total to Ignored
    colID_Height = find(strcmp(data.Properties.VariableNames,'Height'));%from Height to Age
    colID_Age = find(strcmp(data.Properties.VariableNames,'Age'));%from Height to Age

    data = [data(:,1:colID_Flow-1), data(:,colID_Td), data(:,colID_Ignored+1:colID_Height-1), data(:,colID_Age+1:end), data(:,colID_Flow:colID_Reject), data(:,colID_Total:colID_Ignored), data(:,colID_Height:colID_Age)];
    colID_BMI = find(strcmp(data.Properties.VariableNames,'BMI'));%move BMI to end to be next to height/weight/etc.
    data = [data(:,1:colID_BMI-1), data(:,colID_BMI+1:end), data(:,colID_BMI)]; 
        
    % Calculate Xmax (the X value that corresponds to Rcentre. It is the max point of the semi-circle) 
    data.Xmax = sqrt(data.Zchar.^2 - data.Rcentre.^2); % Zchar = sqrt( Rcenter^2 + Xmax^2) 
    colID_Xcentre = find(strcmp(data.Properties.VariableNames,'Xcentre'));  
    data = [data(:,1:(colID_Xcentre-1)), data(:,end), data(:,colID_Xcentre:end-1)]; % reorganize so that Xmax is next to Xcentre
    
    % Convert Subject Name to categorical 
    data.SubjectName = categorical(data.SubjectName); 
    
    % Create columns to identify Last and First timePts within each subject 
    uniqueSortedTimePts = sort(unique(data.TimePt)); 
    uniqueSortedTimePts(uniqueSortedTimePts == connectTick | uniqueSortedTimePts == disconnectTick | uniqueSortedTimePts == continuousTick | uniqueSortedTimePts ==  bolusTick) = []; % *** Add tick types here if they are added above 
    lastTimePt = uniqueSortedTimePts(end); 
    secondToLastTimePt = uniqueSortedTimePts(end-1); 
    firstTimePt = uniqueSortedTimePts(1); 
    secondTimePt = uniqueSortedTimePts(2); 
    data.TimePt_Last = (data.TimePt == lastTimePt); 
    data.TimePt_SecondToLast = (data.TimePt == secondToLastTimePt); 
    data.TimePt_First = (data.TimePt == firstTimePt); 
    data.TimePt_Second = (data.TimePt == secondTimePt); 
    
    %% Identify & Delete Bad Measurements
    %   Criteria: 
    %       1. Rinf <= 0 (a negative resistivity doesn't make physical sense)
    %       2. Mcap <= 0 (negative membrane capacitance also doesn't make sense)
    %       3. Td <= 0 

    % Put the rows for bad measurements in a separate table 
    tempBadTrials = data(data.Rinf <= 0 | data.Td == 0 | data.Mcap <=0,:); 
    badTrials = vertcat(badTrials, tempBadTrials); 
    % Delete the bad rows from main data table
    data(data.Rinf <=0,:) = []; 
    data(data.Td == 0,:) =[];
    data(data.Mcap <=0,:) = []; 
    
    %% Plot for checking individual trials 
    %  This section plots all trial results as Resistivity vs. Reactance 
    %  Whole Body
    data_wb = data(data.WholeBody ==1,:); 
    data_wb(isnan(data_wb.TimePt),:) = []; 
    data_wb(data_wb.TimePt==98,:) = []; 
    data_wb(data_wb.TimePt==-0.5,:) = []; 
    figure
    scatter(data_wb.Rcentre, data_wb.Xmax); 
    text(data_wb.Rcentre, data_wb.Xmax, strcat('timePt',num2str(data_wb.TimePt), '_',num2str(data_wb.Trial)))
    xlabel('R_{center}'); ylabel('X_c max')
    %xlim([300 700]); ylim([20 80]); 
    set(gcf, 'color', 'w');set(gca, 'FontSize', 16);
    title(sprintf('%s Whole Body BI', fieldnames_subjects{ii}))
     
    %  Leg
    data_lg = data(data.WholeBody ==0,:); 
    data_lg(isnan(data_lg.TimePt),:) = []; 
    data_lg(data_lg.TimePt==98,:) = []; 
    data_lg(data_lg.TimePt==-0.5,:) = []; 
    figure
    scatter(data_lg.Rcentre, data_lg.Xmax); 
    text(data_lg.Rcentre, data_lg.Xmax, strcat('timePt',num2str(data_lg.TimePt), '_',num2str(data_lg.Trial)))
    xlabel('R_{center}'); ylabel('X_c max')
    set(gcf, 'color', 'w');set(gca, 'FontSize', 16);
    title(sprintf('%s Leg Segmental BI', fieldnames_subjects{ii}))
    
    %% Merge with circumferences data
    % Prep 'Lengths' table
    colsSubjectName = find(strcmpi(demographics.Properties.VariableNames,'SubjectName'));
    colsCalf = find(contains(demographics.Properties.VariableNames,'Calf'));
    lengths = demographics(:,[colsSubjectName,colsCalf]); 
    lengths.C1 = mean([lengths.CalfMajorCircumferencePre, lengths.CalfMajorCircumferencePost],2,'omitnan'); 
    lengths.C2 = mean([lengths.CalfMinorCircumferencePre, lengths.CalfMinorCircumferencePost],2,'omitnan'); 
    % Merge BI and Lengths Datasets
    data = join(data, lengths, 'Keys','SubjectName');

    %% Re-calculate Leg Segmental TBW, ECF, ICF results 
    % Add columns that are added to leg data
    data_wb = data(data.WholeBody ==1,:); 
    data_wb.ECF_calfLength(:) = NaN; 
    data_wb.ICF_calfLength(:) = NaN; 
    data_wb.TBW_calfLength(:) = NaN; 
    data_wb.ECF1_calfLength(:) = NaN; 
    data_wb.ICF1_calfLength(:) = NaN; 
    data_wb.TBW1_calfLength(:) = NaN; 
    % Calculate new leg volumnes
    data_lg = data(data.WholeBody ==0,:);
    data_lg = calculateLegVolumesBI(data_lg); 
    % Merge data_wb and data_lg again into one dataset
    data = vertcat(data_wb, data_lg);     
    
    %% Average & Stdev of Trials  
    % Calculate averages and stdevs 
    colID_Date = find(strcmp(data.Properties.VariableNames,'Date'));  
    data_mean = varfun(@mean, data, 'InputVariables',[colID_Date:width(data)],'GroupingVariables', {'SubjectName','Subject','VisitNum','HD','TimePt','WholeBody'}); % <-- I'm doing an average of all the values, might not be the right thing. Maybe have to delete values that are a ceratin % away from the mean/median or something like that (?) 
    data_std = varfun(@std, data, 'InputVariables',[colID_Date:width(data)],'GroupingVariables', {'SubjectName','Subject','VisitNum','HD','TimePt','WholeBody'}); 
    % Rename columns
    colID_mean_Date = find(strcmp(data_mean.Properties.VariableNames,'mean_Date')); 
    data_mean.Properties.VariableNames(colID_mean_Date:end) = data.Properties.VariableNames(colID_Date:end); 
    data_std.Properties.VariableNames(colID_mean_Date:end) = strcat(data.Properties.VariableNames(colID_Date:end),'_std'); 
    % Join mean and std tables
    colID_mean_FFM1 = find(strcmp(data_mean.Properties.VariableNames,'FFM1')); 
    colID_mean_ECFcalfLength = find(strcmp(data_mean.Properties.VariableNames,'ECF_calfLength'));
    data_mean_std_JustColsWithStdev = zipFastener(data_mean(:,[colID_mean_Date:colID_mean_FFM1,colID_mean_ECFcalfLength:end]), data_std(:,[colID_mean_Date:colID_mean_FFM1,colID_mean_ECFcalfLength:end])); % Just cols I want to have a stdev for
    data_mean_std_AllCols = [data_mean(:,1:colID_mean_Date-1), data_mean_std_JustColsWithStdev, data_mean(:,colID_mean_FFM1+1:(colID_mean_ECFcalfLength-1))]; %vertcat all columns
    
    %% Put results from this subject into master data tables
    masterAverages = vertcat(masterAverages,data_mean_std_AllCols); 
    masterIndividualTrials = vertcat(masterIndividualTrials, data); 
end


%% Save results
bi_individuals = masterIndividualTrials; 
bi = masterAverages; 
save(fullfile(savePath, 'BI_AllIndividualTrials.mat'), 'bi_individuals')
save(fullfile(savePath, 'BI_Averages.mat'), 'bi')
save(fullfile(savePath, 'BI_BadTrials.mat'), 'badTrials')

