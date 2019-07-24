%% Multi-Exponential Fitting of T2 Data from MR Sensor (analysis_MultiExpT2Fit)  
%  This file processes data (1Exp and 3Exp) from a single patient and outputs graphs
%  illustrating the goodness of the fits and a matrix with the fitting
%  results. 
%
%  This script used to be called analysis_processNMR.m (in 'Dialysis' folder)
% 
%  2018 NOTE: I don't use the fits from this script anymore. I do use the
%  processing. Purpose of this script is now to: 
%       Process raw data from each subject and each time point. Data from all trials are averaged together
%       (straight averaging). From now on can use data from these .mat files rather than having
%       to touch raw data again.
%
%  PREVIOUS SCRIPTS
%       none. This script takes raw data. 
%       It requires reading from config_mrSensor.csv for file locations (or other config files)
%
%  Outputs: 
%      + Plots_StdDevofTrials.png -- Plot the standard deviation that you get from 
%                                     averaging together the 3 (or however many) trials
%                                     were collected for that data point.
%                                     High std. dev. (normal=below 3 or 4) shows that trials
%                                     are very different from each other. Take time to look 
%                                     at them closely. Might not want to use all
%                                     trial data for the final analysis.                                     
%      + Plots_T2Decays.png -- For each time point, plot raw data overlayed
%                              with 1-, 2-, and 3-exp fitting
%      + fitResults_Leg_1Exp.csv -- Plots results for 1-exp fitting for all
%                                   time points. SNR, Rsquare, SSE, RMSE,
%                                   StartingVals, Fitting Results, nTrials,
%                                   StdDev Btw Trials, TE, Who Did Analysis, etc.
%      + fitResults_Leg_2exp.csv -- Ditto for 2-exponential fitting
%      + fitResults_Leg_3exp.csv -- Ditto for 3-exponential fitting
%
%  Matlab Variable Names
%      + time -- (nrows=nEchoes, ncols=nTimePts) Time array custom generated based on the TE for each time point  
%      + amp -- (nrows=nEchoes, ncols=nTimePts) Amplitude array for each time point 
%      + save1exp, save2exp, save3exp -- the cell arrays that are saved into a .csv file (explained in 'Outputs' section above)
%      + stdev -- (nrows=nEchoes, ncols=nTimePts) the data that is plotted in Figure 2 (Plots_StdDevofTrials.png explained above)
%      + trialTE -- (nrows=nTrials, ncols=nTimePts) for each time point there are 3 (or more) trials taken. This matrix shows the TE used in each of those trials. 
%      + rawAmps -- structure containing data from the (usually) 3 trials that then get averaged together to make the amp array for a particular time point: 'rawAmps(1).name'
%                   each field in the structure is a time point (in chronological order). Each field contains a 3 column matrix (nrows=nEchoes, ncols=nTrials).
%
%  Lina A. Colucci 2015
%
%  2015-09-09 The code only works for a config.csv file with only 1 patient
%  2015-09-15 The code pulls in the 1 patient from the 'User Inputs' (still only 1 patient but automatically picks from config file)
%  2015-09-28 The code generates 'tidy data' outputs. Added subject, timePt, Fitting, Sensor to the saved results
%  2015-10-31 Improved help section 
%  2015-11-03 Added 'Finger' Sensor analysis capability. Choose whether analyzing 'leg' (up to 3exp fit) or 'finger' (up to 5exp fit) sensor data. 
%  2016-08-09 Made code compatible with MRI sensor data
%  2016-08-10 Added capability of analyzing phantom data. Added more columns to 'importconfigfile' function
%  2016-11-14 Added LegPhantomPilot (only up to 2exp fit) and LegPhantom (up to 3exp fit), added columns in config file
%
%  FUTURE TO DOs
%       + make a single config xlsx file with different tabs for the different studies (rather than a different csv file for each study)
%       + modernize the whole script by using tables (I didn't know about them at the time so read in data all as individual vectors. Way messier than necessary) 
%       + modernize how I save results (one single table or structure instead of separate .mat file for each time point data)
%  --------------------------------------------------------------------------------
clear; clc; close all; workingDirectory = pwd; 

%% -------------------------- USER INPUTS ----------------------------- %%
% Which patient(s) do you want to analyze? (Enter in their ID number)
patient = 7; 

% Is this the normal HD study or the MRI dialysis study?
study = 'MRI'; %options: 'HD' (pilot study) or 'MRI' or 'Other'. Different file paths and config files are used for each. 
% if 'Other', must fill out these next variables: 
importPath = '/Users/linacolucci/Documents/GitHub/Dialysis/config_mri.csv'; %where config.csv file is located
baseFolder = '/Volumes/CimaLabLargeFiles/lcolucci/MRI_Study/Data/HDMRI'; %path to folder that holds subfolders whose names are listed in legTimeA, legTimeB, etc. 
savePath = '/Volumes/CimaLabLargeFiles/lcolucci/MRI_Study/Analysis/HDHMRI01b/Leg'; %where you want results saved

% Analyze leg or finger data? 
whichSensor = 'Leg'; %Options: "Leg" or "Finger" or "FingerPhantom" or "LegPhantom" or "LegPhantomPilot" (NMR MOUSE leg sensor used during first 30 HD patients)

% How many initial points to ignore during fitting? 0 means starting on pt #1. 1 means starting on pt #2, etc. 
ignore = 1; % ignore for finger and fingerphantom = 2
            % ignore for leg = 0 or 1 
            % ignore for legphantompilot = 1

% How many dummy echoes are there and what is the TE? 
nDummies = 3; % Make sure that all TEs on the patient are the same (error message will pop up if they are not). 
%TE = .1; % (ms) i.e. 65us = 0.065ms OR 300us=0.3ms

%Which trials to skip (if any)?
skip = NaN; %This is a single numer or an array [2,5,8] that tells the code which trial numbers to skip and not do a fitting for. 
            
% Starting values for fittings? %If something here use it to overwrite, otherwise just use whatever is in function code. 
startPoint1Exp = [9 35]; %e.g. [AmpVal RelaxVal] No comma in between; %leg: #25[3 80]; #6 [8 50] %LegPhantomNew [9 35]
startPoint2Exp = [1 5 9 100]; %finger: 50 80 100 200 %LegPhantomNew [1 1 9 30]
startPoint3Exp = [5 35 10 100 5 300]; %leg Old [2 10 2 50 2 100] [1 15 3 60 1 200] %newLeg [5 35 10 100 5 300] %finger: 100 50 200 100 100 200 %Matt's [100 10 100 100 50 1000] %LegPhantomNew [1 1 7 30 2 65] [0.5 7 10 30 1.5 0.9] 
startPoint4Exp = [200 30 300 80 80 200 10 500]; %Matt's [20 5 40 20 40 100 5 350]
startPoint5Exp = [20 5 50 30 100 80 40 150 15 300]; %Matt's: [5 5 80 25 200 60 100 150 30 400] 

% Do you want to save the resuts?
doSave = 'Yes'; %Options: Yes, yes, No, no

%Who is running this analysis?
name = 'LAC'; 

%% -------------------------------------------------------------------- %%


%% Connect
% Connect to folders that contain scripts I need
path(path, '/Users/linacolucci/Documents/1Lab/Dialysis Study/Analysis/Dialysis/functions'); 
path(path, fullfile(workingDirectory, 'functions')); 

%  Connect to the folders
% ***!!!!*** Any time there are additional columns added to config.csv, this line must change ****
switch study 
    case 'HD' 
        importPath = fullfile(workingDirectory,'config_pilot.csv');
    case 'MRI'
        importPath = fullfile(workingDirectory, 'config_mrSensor.csv'); 
    case {'Other'}
        importPath = importPath;     
end
[ID,filename,numPtsLeg,legTimeA,legTimeB,legTimeC,legTimeD,legTimeE,legTimeF,legTimeG,legTimeH,legTimeI,legTimeJ,legTimeK, legTimeL, legTimeM, legTimeN, legTimeO, legTimeP, legTimeQ, numPtsFinger,fingerTimeA,fingerTimeB,fingerTimeC,fingerTimeD,fingerTimeE,fingerTimeF,fingerTimeG,fingerTimeH,fingerTimeI,fingerTimeJ, fingerTimeK, fingerTimeL, fingerTimeM, fingerTimeN, fingerTimeO, fingerTimeP, fingerTimeQ] = importconfigfile2(importPath); 

% Analyze the right patient. Generate the right file names. 
rowNumToAnalyze = find(ID==patient); %find row # of subject you analyze
filename = filename{rowNumToAnalyze};
nTimePtsLeg = numPtsLeg(rowNumToAnalyze);
legTimeA = legTimeA{rowNumToAnalyze};
legTimeB = legTimeB{rowNumToAnalyze};
legTimeC = legTimeC{rowNumToAnalyze}; 
legTimeD = legTimeD{rowNumToAnalyze};
legTimeE = legTimeE{rowNumToAnalyze};
legTimeF = legTimeF{rowNumToAnalyze};
legTimeG = legTimeG{rowNumToAnalyze};
legTimeH = legTimeH{rowNumToAnalyze};
legTimeI = legTimeI{rowNumToAnalyze};
legTimeJ = legTimeJ{rowNumToAnalyze};
legTimeK = legTimeK{rowNumToAnalyze};
legTimeL = legTimeL{rowNumToAnalyze};
legTimeM = legTimeM{rowNumToAnalyze};
legTimeN = legTimeN{rowNumToAnalyze};
legTimeO = legTimeO{rowNumToAnalyze};
legTimeP = legTimeP{rowNumToAnalyze};
legTimeQ = legTimeQ{rowNumToAnalyze};
nTimePtsFinger = numPtsFinger(rowNumToAnalyze);
fingerTimeA = fingerTimeA{rowNumToAnalyze};
fingerTimeB = fingerTimeB{rowNumToAnalyze};
fingerTimeC = fingerTimeC{rowNumToAnalyze};
fingerTimeD = fingerTimeD{rowNumToAnalyze};
fingerTimeE = fingerTimeE{rowNumToAnalyze};
fingerTimeF = fingerTimeF{rowNumToAnalyze};
fingerTimeG = fingerTimeG{rowNumToAnalyze};
fingerTimeH = fingerTimeH{rowNumToAnalyze};
fingerTimeI = fingerTimeI{rowNumToAnalyze};
fingerTimeJ = fingerTimeJ{rowNumToAnalyze};
fingerTimeK = fingerTimeK{rowNumToAnalyze};
fingerTimeL = fingerTimeL{rowNumToAnalyze};
fingerTimeM = fingerTimeM{rowNumToAnalyze};
fingerTimeN = fingerTimeN{rowNumToAnalyze};
fingerTimeO = fingerTimeO{rowNumToAnalyze};
fingerTimeP = fingerTimeP{rowNumToAnalyze};
fingerTimeQ = fingerTimeQ{rowNumToAnalyze};

switch study 
    case 'HD' 
        baseFolder = '/Volumes/cimalab/Hydration/MGH/Lenovo_NMR/MGH';      
    case 'MRI'
        baseFolder = '/Volumes/cimalab/Hydration/MRI_Study/Data'; 
    case 'MRI_Lenovo' %When using Lenovo laptop for this analysis
        baseFolder = 'Z:\Hydration\MRI_Study\Data'; 
    case 'HD_Lenovo' %When using Lenovo laptop for this analysis
        baseFolder = 'Z:\Hydration\MGH\Lenovo_NMR\MGH'; 
    case 'Other'
        baseFolder = baseFolder;     
end 
fileFormat = 'data2.csv'; 

% Generate a structure so that folders can be referenced in the loop & define nTimePts
switch whichSensor 
    case {'Leg','LegPhantom','LegPhantomPilot'}
        field = 'name'; 
        value = {legTimeA,legTimeB, legTimeC,legTimeD,legTimeE,legTimeF,legTimeG,legTimeH,legTimeI, legTimeJ,legTimeK, legTimeL, legTimeM, legTimeN, legTimeO, legTimeP, legTimeQ}; 
        sensorTime = struct(field, value); 
        nTimePts = nTimePtsLeg; 
        whichSensorCat = 'Leg'; %sensor category. For the purposes of file paths, Leg phantom data is in same Leg folder
    case {'Finger','FingerPhantom'}
        field = 'name'; 
        value = {fingerTimeA, fingerTimeB, fingerTimeC, fingerTimeD, fingerTimeE, fingerTimeF, fingerTimeG, fingerTimeH,fingerTimeI,fingerTimeJ, fingerTimeK, fingerTimeL, fingerTimeM, fingerTimeN, fingerTimeO, fingerTimeP, fingerTimeQ}; 
        sensorTime = struct(field, value); 
        nTimePts = nTimePtsFinger; 
        whichSensorCat = 'Finger'; 
end

%% Load and Generate 'amp' Array
% This section loads in data from each time point, averages the trials of
% each time point together, and outputs a matrix ('amp') where each column
% is the amplitude data for each time point. It also outputs a matrix
% ('stdev') that is the standard deviation of the 3 amplitudes from the 3 
% trials collected per time point
% 
% See rawAmps (i.e. amplitude for each of the 3 trials before they are
% averaged together) by accessing: e.g. 'rawAmps(2).name' 
[amp, stdev, nEchoes, trials, rawAmps] = generateAmpArray(nTimePts, 1, 2, baseFolder, filename, whichSensorCat,sensorTime, fileFormat, study)

%% Find TE and Generate 'time' Array
[time, trialTE] = generateTimeArray(nTimePts, nEchoes, baseFolder, filename, whichSensorCat, sensorTime, study); 
    
%% Analysis 
[results1Exp, results2Exp, results3Exp, results4Exp, results5Exp] = analysisT2Decay(nTimePts, time, amp, ignore, skip, whichSensor, startPoint1Exp, startPoint2Exp, startPoint3Exp, startPoint4Exp, startPoint5Exp);

%% Create a String that Identifies Each Time Point (i.e Pre, 0hr, 1hr,
% etc.) by removing 'CPMG_Leg_' and '_hd' from the 'legTime' variable
sensorTimeString = generateTimePtStringFromFilename(sensorTime, nTimePts, study); %******* WHAT ABOUT FOR FINGER? OR PHANTOM?

%% Look at Individual Trials (before averaging)
plotIndividualTrials(nTimePts, time, baseFolder, filename, whichSensorCat,sensorTime, fileFormat, study);

%% Plots
% Plot the T2 Decays for Each Time Point 
% BEWARE! If code order changes here, change the 'Export' section as well, i.e. '-f1' and '-f2'
plotT2DecaysAllTimePts(nTimePts, ignore, time, amp, results1Exp, results2Exp, results3Exp, results4Exp, results5Exp, sensorTime, filename, whichSensor)

% Residuals Plot
plotStdDevofTrials(nTimePts, stdev, sensorTime, filename)

%% Export

% Move into working directory 
cd(workingDirectory); % Move into pwd if not already. pwd contains the scripts & functions that are backed into GitHub

switch study
    case 'HD'
        basePath = '/Volumes/cimalab/Hydration/MGH/outputs'; 
        
        % Create 'outputs' folder if it doesn't exist
        if ~exist(basePath,'dir')
            mkdir(basePath);
        end

        % Create 'outputs' folder for this patient if it doesn't exist
        if ~exist(fullfile(basePath,filename,whichSensor),'dir')
            mkdir(fullfile(basePath,filename,whichSensor))
        end

        % Define Save Path within the 'outputs' folder
        savePath = fullfile(basePath, filename,whichSensor); 

    case 'MRI'
        % Create 'Analysis' folder if it doesn't exist
        if ~exist(fullfile('/Volumes/cimalab/Hydration/MRI_Study/Analysis',filename),'dir')
            mkdir(fullfile('/Volumes/cimalab/Hydration/MRI_Study/Analysis',filename));
        end
        
        % Within 'Analysis' folder, create 'Leg' and 'Finger' subfolders if they don't exist. 
        if ~exist(fullfile('/Volumes/cimalab/Hydration/MRI_Study/Analysis',filename,whichSensor),'dir')
            mkdir(fullfile('/Volumes/cimalab/Hydration/MRI_Study/Analysis',filename,whichSensor));
        end
       
        % Define Save Path within the 'Analysis' folder
        savePath = fullfile('/Volumes/cimalab/Hydration/MRI_Study/Analysis',filename,whichSensor); 
    
    case 'Other'
        savePath = savePath; 
        if ~exist(savePath, 'dir')
            mkdir(savePath)
        end
end

% Add Column: Time this file was generated 
date(1:nTimePts, 1)= {datestr(now, 'yyyy-mm-dd HH:MM:SS')}; 
% Add Column: Subject ID 
subject(1:nTimePts,1)=patient; 
% Add a proper time point column
timePt(1:nTimePts, 1)=sensorTimeString';
% Add Column: Fitting Type 
fitting1(1:nTimePts,1)=1; 
fitting2(1:nTimePts,1)=2; 
fitting3(1:nTimePts,1)=3; 
fitting4(1:nTimePts,1)=4;
fitting5(1:nTimePts,1)=5; 
% Add Column: Sensor Type 
switch whichSensor
    case {'Leg','LegPhantom','LegPhantomPilot'}
        sensor(1:nTimePts,1)={'leg'};
    case {'Finger','FingerPhantom'}
        sensor(1:nTimePts,1)={'finger'}; 
end
% Add Column: User
user(1:nTimePts, 1) = {name};
% Add Column: TE
teArray = time(1,:)';
if nDummies==0
    teArray = time(1,:)'; % if 0 dummy echoes, then TE is just the 1st value in time Array
else
    teArray = time(1,:)'/(nDummies+1);
    %teArray(1:nTimePts,1) = TE; % otherwise, I need to input the TE
end
% Add Column: Median of Std Deviation btw the Trials
stdevTrials = median(stdev)'; 
% Add Column: # points that were ignored in fitting
ignoreArray(1:nTimePts, 1) = ignore; 
% Add Column: # of dummy Echoes (nDummies)
nDummiesArray(1:nTimePts,1)= nDummies; 
% Add Column: skipped trials
%skipArray(1:nTimePts,1) = skip; %Can't figure out how to add this if it's an array (not just a single value)

% Create header strings
header1ExpResults = {'subject','timePt','nMeas','snr','rSquare','sse','rmse','startValAmp1','startValRelax1', 'amp','ampStd','relax','relaxStd',...
    'nTrials','stdevBtwTrials','te','nDummies','sensor','fitting','ignore','date', 'user'}; 
header2ExpResults = {'subject','timePt','nMeas','snr','rSquare','sse','rmse','startValAmp1','startValRelax1', 'startValAmp2','startValRelax2', 'amp1','ampStd1','relax1',...
    'relaxStd1', 'amp2','ampStd2','relax2','relaxStd2',...
    'nTrials','stdevBtwTrials','te','nDummies','sensor','fitting','ignore','date', 'user'};
header3ExpResults = {'subject','timePt','nMeas','snr','rSquare','sse','rmse','startValAmp1','startValRelax1', 'startValAmp2','startValRelax2', 'startValAmp3','startValRelax3',...
    'amp1','ampStd1','relax1','relaxStd1', 'amp2','ampStd2','relax2','relaxStd2', 'amp3','ampStd3','relax3','relaxStd3',...
    'nTrials','stdevBtwTrials','te','nDummies','sensor','fitting','ignore','date', 'user'};
header4ExpResults = {'subject','timePt','nMeas','snr','rSquare','sse','rmse','startValAmp1','startValRelax1', 'startValAmp2','startValRelax2', 'startValAmp3','startValRelax3','startValAmp4','startValRelax4',...
    'amp1','ampStd1','relax1','relaxStd1', 'amp2','ampStd2','relax2','relaxStd2', 'amp3','ampStd3','relax3','relaxStd3','amp4','ampStd4','relax4','relaxStd4'...
    'nTrials','stdevBtwTrials','te','nDummies','sensor','fitting','ignore','date', 'user'};
header5ExpResults = {'subject','timePt','nMeas','snr','rSquare','sse','rmse','startValAmp1','startValRelax1', 'startValAmp2','startValRelax2', 'startValAmp3','startValRelax3','startValAmp4','startValRelax4','startValAmp5','startValRelax5',...
    'amp1','ampStd1','relax1','relaxStd1', 'amp2','ampStd2','relax2','relaxStd2', 'amp3','ampStd3','relax3','relaxStd3','amp4','ampStd4','relax4','relaxStd4','amp5','ampStd5','relax5','relaxStd5',...
    'nTrials','stdevBtwTrials','te','nDummies','sensor','fitting','ignore','date', 'user'};

% Value to reference in cell arrays
saveConstant = 11; % for the existing headers, the value is 11. There are 11 additional columns beyond the ones in "results" (which go from nMeas-->relaxStd).
                   % if you add more columns, need to add the number of new columns to this number. 

% Initialize the blank cell array to put results into
save1exp = cell(size(results1Exp,1)+1, size(results1Exp,2)+saveConstant); 
save2exp = cell(size(results2Exp,1)+1, size(results2Exp,2)+saveConstant); 

% add the headers
save1exp(1,1:end) = header1ExpResults; 
save2exp(1,:) = header2ExpResults; 

% add the subject and time pt
save1exp(2:end,1:2) = horzcat(num2cell(subject), timePt); 
save2exp(2:end,1:2) = horzcat(num2cell(subject), timePt); 

% add trials, stdev btw trials, sensor, fitting, date, user
save1exp(2:end, end-(saveConstant-3):end) = horzcat(num2cell(trials), num2cell(stdevTrials), num2cell(teArray), num2cell(nDummiesArray), sensor, num2cell(fitting1), num2cell(ignoreArray), date, user); %Edit "-7" if you change the 'Save' info (# should be 1 less than length of array)
save2exp(2:end, end-(saveConstant-3):end) = horzcat(num2cell(trials), num2cell(stdevTrials), num2cell(teArray), num2cell(nDummiesArray), sensor, num2cell(fitting2), num2cell(ignoreArray), date, user);

% add fitting results
save1exp(2:end, 3:end-(saveConstant-2)) = num2cell(results1Exp); %Edit "-8" if you change the 'Save' info (should be one above the value in '%add trials')
save2exp(2:end, 3:end-(saveConstant-2)) = num2cell(results2Exp);

switch whichSensor
    case {'Finger','Leg','LegPhantom'}
        %initialize blank
        save3exp = cell(size(results3Exp,1)+1, size(results3Exp,2)+saveConstant); 
        %add headers
        save3exp(1,:) = header3ExpResults; 
        %add subject and time pt
        save3exp(2:end,1:2) = horzcat(num2cell(subject), timePt); 
        % add trials, stdev btw trials, sensor, fitting, date, user
        save3exp(2:end, end-(saveConstant-3):end) = horzcat(num2cell(trials), num2cell(stdevTrials), num2cell(teArray), num2cell(nDummiesArray), sensor, num2cell(fitting3), num2cell(ignoreArray), date, user);
        % add fitting results
        save3exp(2:end, 3:end-(saveConstant-2)) = num2cell(results3Exp);
end

switch whichSensor %do all 'Finger' steps at once for the 4- and 5-exp (so that I don't have to do switches every time)
    case 'Finger'
        %initialize blank
        save4exp = cell(size(results4Exp,1)+1, size(results4Exp,2)+saveConstant); 
        save5exp = cell(size(results5Exp,1)+1, size(results5Exp,2)+saveConstant);         
        %add headers
        save4exp(1,:) = header4ExpResults; 
        save5exp(1,:) = header5ExpResults; 
        %add subject and time pt
        save4exp(2:end,1:2) = horzcat(num2cell(subject), timePt);
        save5exp(2:end,1:2) = horzcat(num2cell(subject), timePt);
        % add trials, stdev btw trials, sensor, fitting, date, user
        save4exp(2:end, end-(saveConstant-3):end) = horzcat(num2cell(trials), num2cell(stdevTrials), num2cell(teArray), num2cell(nDummiesArray), sensor, num2cell(fitting4), num2cell(ignoreArray), date, user); %Edit "-7" if you change the 'Save' info (# should be 1 less than length of array)
        save5exp(2:end, end-(saveConstant-3):end) = horzcat(num2cell(trials), num2cell(stdevTrials), num2cell(teArray), num2cell(nDummiesArray), sensor, num2cell(fitting5), num2cell(ignoreArray), date, user);
        % add fitting results
        save4exp(2:end, 3:end-(saveConstant-2)) = num2cell(results4Exp); %Edit "-8" if you change the 'Save' info (should be one above the value in '%add trials')
        save5exp(2:end, 3:end-(saveConstant-2)) = num2cell(results5Exp);
end

% Save the results
if any(strcmp(doSave, {'Yes','YES','yes','y','Y'}))
    switch whichSensor
        case {'Leg','LegPhantom','LegPhantomPilot'}
            cell2csv(fullfile(savePath,  'Leg_fitResults_1Exp.csv'), save1exp, ',')
            cell2csv(fullfile(savePath,  'Leg_fitResults_2Exp.csv'), save2exp, ',')
            
            switch whichSensor
                case {'Leg','LegPhantom'}
                    cell2csv(fullfile(savePath,  'Leg_fitResults_3Exp.csv'), save3exp, ',')
            end
            
            print(fullfile(savePath,  'Leg_Plots_IndividualTrials'), '-f1','-dpng') % BEWARE! If code order changes, these might not refer to right figures any more
            print(fullfile(savePath,  'Leg_Plots_T2Decays'), '-f2','-dpng') 
            print(fullfile(savePath,  'Leg_Plots_StdDevofTrials'),'-f3','-dpng')
            
            savefig(figure(1),fullfile(savePath,  'Leg_Plots_IndividualTrials'))% Save as Matlab figure
            savefig(figure(2),fullfile(savePath,  'Leg_Plots_T2Decays')) 
            savefig(figure(3),fullfile(savePath,  'Leg_Plots_StdDevofTrials'))
            
            disp('Results and plots have been saved')
        case {'Finger','FingerPhantom'}
            cell2csv(fullfile(savePath,  'Finger_fitResults_1Exp.csv'), save1exp, ',')
            cell2csv(fullfile(savePath,  'Finger_fitResults_2Exp.csv'), save2exp, ',')

            switch whichSensor
                case 'Finger'
                    cell2csv(fullfile(savePath,  'Finger_fitResults_3Exp.csv'), save3exp, ',')
                    cell2csv(fullfile(savePath,  'Finger_fitResults_4Exp.csv'), save4exp, ',')
                    cell2csv(fullfile(savePath,  'Finger_fitResults_5Exp.csv'), save5exp, ',')            
            end
            
            print(fullfile(savePath,  'Finger_Plots_IndividualTrials'), '-f1','-dpng') % BEWARE! If code order changes, these might not refer to right figures any more
            print(fullfile(savePath,  'Finger_Plots_T2Decays'), '-f2','-dpng') 
            print(fullfile(savePath,  'Finger_Plots_StdDevofTrials'),'-f3','-dpng')

            savefig(figure(1),fullfile(savePath,  'Finger_Plots_IndividualTrials')) % Save as Matlab figure
            savefig(figure(2),fullfile(savePath,  'Finger_Plots_T2Decays')) 
            savefig(figure(3),fullfile(savePath,  'Finger_Plots_StdDevofTrials'))
            
            disp('Results and plots have been saved')
    end
    % Save amplitude as MAT files
    for i=1:size(amp,2)
        varSavePath = fullfile(savePath, 'ProcessedData');
        if ~exist(varSavePath, 'dir')
            mkdir(varSavePath); 
        end
        tempAmp = amp(:,i);
        save(fullfile(varSavePath, strcat(timePt{i},'.mat')), 'tempAmp')
        tempTime = time(:,i); 
        save(fullfile(varSavePath, strcat(timePt{i},'_time','.mat')), 'tempTime')
    end
    
else 
    disp('Results and plots were not saved. Edit "save" variable if you would like to save.')
end

