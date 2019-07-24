%% Demographics (REDCap data analysis)
clear

%% USER INPUTS
% CSVs to REDCap data: all non-PHI data with labels
pathHD = '/MRI_Study/ProcessedData/redcap/MRIHDStudy_DATA_LABELS_2018-07-09_1451.csv';
pathHC = '/MRI_Study/ProcessedData/redcap/HCMRINMRStudy_DATA_LABELS_2018-09-27_1210.csv'; 
% Save
doSave = 1; 
savePath = '/Volumes/cimalablargefiles/Hydration/MRI_Study/ProcessedData/redcap'; 
%%

%% Manually define subject ID lookuptable
idLookup = table(); 
% idLookup.ID = [1, 1, 2, 2, 3, 4, 4, 5,101,102,103,104,105,101,106]'; 
% idLookup.Visit = [1, 2, 1, 2, 1, 1, 2, 1,1,1,1,1,1,2,1]'; 
% idLookup.HD = [1 1 1 1 1 1 1 1 0 0 0 0 0 0 0]'; 
idLookup.RecordID = categorical({'HDMRI-01', 'HDMRI-01B', 'HDMRI-02', 'HDMRI-02b','HDMRI-03','HDMRI-04','HDMRI-04b','HDMRI-05','HC001','HC002','HC003','HC004','HC005','HC-001b','HC006'}'); 


%% Load data
dataHD = importRedcapHD(pathHD); 
dataHC = importRedcapHC(pathHC); 

%% Select subset
dataHD2 = dataHD(:,[1,61,63:69,86:87,109:111,113:119,121:124,126,128,131,134:139,141,144:146,148,152:166]);    
dataHC2 = dataHC(:,[2,4,5:7,8:15,17,19:20,22,24:28,30:46,62:71]);     
    
    
%% paper table subset
% dataHD3 = dataHD2(:,[1,3,5,6,9,35,39:43,46:48]); 
dataHD3 = dataHD2(:,[1,3,5,6,9,35,39, 31:34,37:38,40:49]); 
dataHD3.Properties.VariableNames = {'RecordID','Age','Race','Ethnicity','BMI','WeightPre','WeightPost','CalfLength','CalfElectrodeSpacing','CalfMajorCircumferencePre','CalfMinorCircumferencePre','CalfMajorCircumferencePost','CalfMinorCircumferencePost','SodiumPre','BUNpre','CreatininePre','WBCpre','HemoglobinPre','HematocritPre','PlateletsPre','OsmolalityPre','BNPpre','AlbuminPre'}; 
dataHC3 = dataHC2(:,[1,23,25:26,30,28,29,34:36,38,37,39:49]);
dataHC3.Properties.VariableNames = {'RecordID','Age','Race','Ethnicity','BMI','WeightPre','WeightPost','CalfLength','CalfElectrodeSpacing','CalfMajorCircumferencePre','CalfMinorCircumferencePre','CalfMajorCircumferencePost','CalfMinorCircumferencePost','SodiumPre','BUNpre','CreatininePre','WBCpre','HemoglobinPre','HematocritPre','PlateletsPre','OsmolalityPre','BNPpre','AlbuminPre'}; 

data = vertcat(dataHD3,dataHC3); 
data = join(data,idLookup); 
data.WeightChange = data.WeightPre - data.WeightPost; 
data.WeightChangePercent = data.WeightChange ./ data.WeightPre .*100; 

%% Re-name and add other identifiers
% Subject look-up table
subject_lookup = table(); 
subject_lookup.RecordID = categorical({'HDMRI-01','HDMRI-01B','HDMRI-02','HDMRI-02b','HDMRI-03','HDMRI-04','HDMRI-04b',...
    'HDMRI-05','HC-001b','HC001','HC002','HC003','HC004','HC005','HC006'})'; 
subject_lookup.SubjectName = categorical({'HDMRI01','HDMRI01b','HDMRI02','HDMRI02b','HDMRI03','HDMRI04','HDMRI04b',...
    'HDMRI05','HC01b','HC01','HC02','HC03','HC04','HC05','HC06'})'; 
subject_lookup.Subject = [1 1 2 2 3 4 4 5 101 101 102 103 104 105 106]'; 
subject_lookup.VisitNum = [1 2 1 2 1 1 2 1 1 2 1 1 1 1 1]'; 
subject_lookup.HD = [1 1 1 1 1 1 1 1 0 0 0 0 0 0 0]'; 

% Join subject_lookup table and demographics data
demographics = join(subject_lookup, data, 'Keys','RecordID'); 

%% Save
if doSave==1
    % make folder if it doesn't exist
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end
    save(fullfile(savePath, 'demographics.mat'),'demographics')
end


%% Calculate two-sample t-test p values comparing HCs to HDs for certain values
whichCols = [2,5:14,18:19]; 
for i=whichCols
    tempData = data(:,i); 
    [a p(i) ] = ttest2(table2array(tempData(data.HD==1,1)), table2array(tempData(data.HD==0,1)));
    data.(i)(16,1) = p(i);
end
data = vertcat(data, table(p)); 


%% Save as csv for Rstudio
data(end,:) = [];
writetable(data,'/Users/linacolucci/Documents/Thesis/Paper/data_for_R/demographics.csv')


