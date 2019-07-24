%% Aggregate procesed MR sensor data into a table or structure
%% Generate mrSensor_T2Decays_table.mat and *_struct.mat
%  This script loads in T2 Decay data from the MR sensor and organizes it
%  into a single matlab variable: both structure and table formats
%
%  PREVIOUS SCRIPTS
%       data comes from script: process_mrSensor_rawData.m
% 
%  OUTPUTS
%       data_table and data_struct - all processed MR sensor data from all 
%       subjects and time points in a single table or structure
%
% Lina A. Colucci, 2017

%% --- USER INPUTS ---
doSave = 1; %1 to save
savePath = '/Volumes/CimaLabLargeFiles/Hydration/MRI_Study/ProcessedData/mr_sensor_leg'; 
%% --------------------

%% Leg Sensor Data - organize it into a structure
%  This code was already written in "explore_forced3exp.m" script
baseFolder = '/Volumes/CimaLabLargeFiles'; 
% HC01
temp = load(fullfile(baseFolder, 'HC001/Leg/ProcessedData/CPMG_0hr.mat'));
data_struct.HC01.am = temp.tempAmp; clear temp
temp = load(fullfile(baseFolder, 'HC001/Leg/ProcessedData/CPMG_4hr.mat'));
data_struct.HC01.pm = temp.tempAmp; clear temp
temp = load(fullfile(baseFolder, 'HC001/Leg/ProcessedData/time.mat'));
data_struct.HC01.time = temp.time; clear temp


% HC02
temp = load(fullfile(baseFolder, 'HC002/Leg/ProcessedData/CPMG_0hr_legB.mat')); 
data_struct.HC02.am = temp.tempAmp; clear temp
temp = load(fullfile(baseFolder, 'HC002/Leg/ProcessedData/CPMG_4hr_leg.mat')); 
data_struct.HC02.pm = temp.tempAmp; clear temp
temp = load(fullfile(baseFolder, 'HC002/Leg/ProcessedData/time.mat')); 
data_struct.HC02.time = temp.time; clear temp


% HC03
temp = load(fullfile(baseFolder, 'HC003/Leg/ProcessedData/CPMG_leg_0hr.mat')); 
data_struct.HC03.am = temp.tempAmp; clear temp
temp = load(fullfile(baseFolder, 'HC003/Leg/ProcessedData/CPMG_leg_3hr.mat')); 
data_struct.HC03.pm = temp.tempAmp; clear temp
temp = load(fullfile(baseFolder, 'HC003/Leg/ProcessedData/time.mat')); 


% HC04
temp = load(fullfile(baseFolder, 'HC004/Leg/ProcessedData/CPMG_0hr_leg_real.mat')); 
data_struct.HC04.am = temp.tempAmp; clear temp
temp = load(fullfile(baseFolder, 'HC004/Leg/ProcessedData/CPMG_4hr_Leg.mat')); 
data_struct.HC04.pm = temp.tempAmp; clear temp
temp = load(fullfile(baseFolder, 'HC004/Leg/ProcessedData/time.mat')); 
data_struct.HC04.time = temp.time; clear temp


% HC05
temp = load(fullfile(baseFolder, 'HC005/Leg/ProcessedData/CPMG_0hr_Leg.mat')); 
data_struct.HC05.am = temp.tempAmp; clear temp
temp = load(fullfile(baseFolder, 'HC005/Leg/ProcessedData/CPMG_4hr_Leg.mat')); 
data_struct.HC05.pm = temp.tempAmp; clear temp
temp = load(fullfile(baseFolder, 'HC005/Leg/ProcessedData/time.mat')); 
data_struct.HC05.time = temp.time; clear temp


% HDMRI01
temp = load(fullfile(baseFolder, 'HDMRI01/Leg/ProcessedData/CPMG_Pre.mat')); 
data_struct.HDMRI01.am = temp.tempAmp; clear temp
temp = load(fullfile(baseFolder, 'HDMRI01/Leg/ProcessedData/CPMG_Post.mat')); 
data_struct.HDMRI01.pm = temp.tempAmp; clear temp
temp = load(fullfile(baseFolder, 'HDMRI01/Leg/ProcessedData/time.mat')); 
data_struct.HDMRI01.time = temp.time; clear temp


% HDMRI02
temp = load(fullfile(baseFolder, 'HDMRI02/Leg/ProcessedData/Pre.mat'));
data_struct.HDMRI02.am = temp.tempAmp; clear temp
temp = load(fullfile(baseFolder, 'HDMRI02/Leg/ProcessedData/End.mat'));
data_struct.HDMRI02.pm = temp.tempAmp; clear temp
temp = load(fullfile(baseFolder, 'HDMRI02/Leg/ProcessedData/Pre_time.mat'));
data_struct.HDMRI02.time = temp.tempTime; clear temp


% HDMRI03
temp = load(fullfile(baseFolder, 'HDMRI03/Leg/ProcessedData/cpmg_pre.mat')); 
data_struct.HDMRI03.am = temp.tempAmp; clear temp
temp = load(fullfile(baseFolder, 'HDMRI03/Leg/ProcessedData/cpmg_End.mat')); 
data_struct.HDMRI03.pm = temp.tempAmp; clear temp
temp = load(fullfile(baseFolder, 'HDMRI03/Leg/ProcessedData/time.mat')); 
data_struct.HDMRI03.time = temp.time; clear temp


%HDMRI04b
temp = load(fullfile(baseFolder, 'HDMRI04b/Leg/ProcessedData/cpmg_pre.mat')); 
data_struct.HDMRI04b.am = temp.tempAmp; clear temp
temp = load(fullfile(baseFolder, 'HDMRI04b/Leg/ProcessedData/cpmg_Post.mat')); 
data_struct.HDMRI04b.pm = temp.tempAmp; clear temp
temp = load(fullfile(baseFolder, 'HDMRI04b/Leg/ProcessedData/time.mat')); 
data_struct.HDMRI04b.time = temp.time; clear temp


%---- Prospective Data ---
baseFolder = '/Volumes/CimaLabLargeFiles'; 
temp = load(fullfile(baseFolder, 'HC01b/Leg/ProcessedData/cpmg_Pre_hd.mat')); 
data_struct.HC01b.am = temp.tempAmp; clear temp
temp = load(fullfile(baseFolder, 'HC01b/Leg/ProcessedData/cpmg_4hr_hd.mat')); 
data_struct.HC01b.pm = temp.tempAmp; clear temp
temp = load(fullfile(baseFolder, 'HC01b/Leg/ProcessedData/cpmg_Pre_hd_time.mat')); 
data_struct.HC01b.time = temp.tempTime; clear temp


temp = load(fullfile(baseFolder, 'HC06/Leg/ProcessedData/cpmg_0hr_hd.mat')); 
data_struct.HC06.am = temp.tempAmp; clear temp; 
temp = load(fullfile(baseFolder, 'HC06/Leg/ProcessedData/cpmg_4hr_hd.mat')); 
data_struct.HC06.pm = temp.tempAmp; clear temp; 
temp = load(fullfile(baseFolder, 'HC06/Leg/ProcessedData/cpmg_0hr_hd_time.mat')); 
data_struct.HC06.time = temp.tempTime; clear temp; 


temp = load(fullfile(baseFolder, 'HDMRI01b/Leg/ProcessedData/cpmg_Pre_hd.mat')); 
data_struct.HDMRI01b.am = temp.tempAmp; clear temp
temp = load(fullfile(baseFolder, 'HDMRI01b/Leg/ProcessedData/cpmg_Post_hd.mat')); 
data_struct.HDMRI01b.pm = temp.tempAmp; clear temp
temp = load(fullfile(baseFolder, 'HDMRI01b/Leg/ProcessedData/cpmg_Pre_hd_time.mat')); 
data_struct.HDMRI01b.time = temp.tempTime; clear temp

temp = load(fullfile(baseFolder, 'HDMRI02b/Leg/ProcessedData/cpmg_Pre_hd.mat')); 
data_struct.HDMRI02b.am = temp.tempAmp; clear temp
temp = load(fullfile(baseFolder, 'HDMRI02b/Leg/ProcessedData/cpmg_Post_hd.mat')); 
data_struct.HDMRI02b.pm = temp.tempAmp; clear temp
temp = load(fullfile(baseFolder, 'HDMRI02b/Leg/ProcessedData/cpmg_Pre_hd_time.mat')); 
data_struct.HDMRI02b.time = temp.tempTime; clear temp


temp = load(fullfile(baseFolder, 'HDMRI05/Leg/ProcessedData/cpmg_Pre_hd.mat'));
data_struct.HDMRI05.am = temp.tempAmp; clear temp; 
temp = load(fullfile(baseFolder, 'HDMRI05/Leg/ProcessedData/cpmg_Post_hd.mat'));
data_struct.HDMRI05.pm = temp.tempAmp; clear temp; 
temp = load(fullfile(baseFolder, 'HDMRI05/Leg/ProcessedData/cpmg_Pre_hd_time.mat'));
data_struct.HDMRI05.time = temp.tempTime; clear temp; 


data_struct = orderfields(data_struct); 

%% Table Form - organize data from structure into table format

data_table = table(); 
fnSubjects = fieldnames(data_struct); 
nSubjects = length(fnSubjects); 
whichRow = 0; 
for iSubject = 1:nSubjects
    fnTimePts = fieldnames(data_struct.(fnSubjects{iSubject})); 
    nTimePts = length(fnTimePts); 
    
    for jTimePts = 1:nTimePts
        if strcmpi(fnTimePts{jTimePts},'time')
            continue
        end
        whichRow = whichRow + 1; 
        data_table.Subject(whichRow) = categorical(fnSubjects(iSubject)); 
        data_table.AMPM(whichRow) = categorical(fnTimePts(jTimePts)); 
        data_table.ROI(whichRow) = categorical({'legSensor'}); 
        data_table.T2Decay(whichRow) = {data_struct.(fnSubjects{iSubject}).(fnTimePts{jTimePts})}; 
        data_table.Time(whichRow) = {data_struct.(fnSubjects{iSubject}).time}; 
    end
end

%% Save
if doSave==1
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end
    
    save(fullfile(savePath, 'legSensor_T2Decays_struct.mat'), 'data_struct')
    save(fullfile(savePath, 'legSensor_T2Decays_table.mat'), 'data_table')
end
    
