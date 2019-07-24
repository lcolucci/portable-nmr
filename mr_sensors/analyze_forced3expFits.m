%% Force Leg Sensor Data into 3 exp
%  Based on script in old "Dialysis" folder called "explore_forced3exp.m" -
%  which has a lot of additional code for exploring the fit results
%
%  This script loads in leg sensor data for all patients then
%  applies a forced 3-exp fit to each dataset and saves result.
%
%  It's possible to switch out different fitting types for experimentation purposes. 
%
%  PREVIOUS SCRIPTS
%       data comes from script: aggregate_t2decays_table.m 
%
%  INPUTS
%       tau1 - enforced value of 1st relaxation time
%       tau3 - enforced value of 3rd relaxation time
%
%  OUTPUTS
%       result - is a table with results from the forced 3exp fit. 
%               'legSensor_forced3exp_[tau1]ms_[tau3]ms.mat'
%               Has the following columns; 
%               relax | relaxStd | amp1 | etc. | r2 | sse | etc. | model |
%               lowerBounds | upperBounds | tau1 | tau3 | subject | PM
%       mrSensor_R2 - abbreviated version of the table, that is fed into R for 
%               statistical tests, and shown as table S8 of the STM paper. 
%  
%  NOTE: if fits are not converging, try changing the starting values (need
%        to go into the fitting function)

clear

%% --- USER INPUTS ---
pathToMriStudy = '/Volumes/CimaLabLargeFiles/Hydration'; % Path to folder that contains 'MRI_Study' folder 
tau1 = 40;  % fixed tau 1 (ms). Default = 40ms
tau3 = 250; % fixed tau 3 (ms). Default = 250ms. 
% ---
doSave = 0; 
savePath = '/Volumes/cimalablargefiles/Hydration/MRI_Study/ProcessedData/mr_sensor_leg'; 
%%

%% Load data table
load(fullfile(pathToMriStudy, 'MRI_Study/ProcessedData/mr_sensor_leg/legSensor_T2Decays_table.mat')); % called data_table 
data = data_table; 

%% Define the Equation (fix tau1 and tau3, float everything else)
equation = sprintf('a*exp(-x/%d)+b*exp(-x/r)+c*exp(-x/%d)', tau1, tau3); 

%% Fitting
result = table(); 
listSubjects = unique(data.Subject); 
nSubjects = length(listSubjects); 
listAMPM = categorical({'am','pm'});
nAMPM = length(listAMPM); 
for iSubject = 1:nSubjects
    for jAMPM = 1:nAMPM
        % Prep data
        temp = data(data.Subject == listSubjects(iSubject) & data.AMPM==listAMPM(jAMPM),:);
        x = temp.Time{1,1}; 
        y = temp.T2Decay{1,1}; 
        
        % Fitting
        tempResult1 = t2_fitting_3exp_floatAllAmpsAnd1Relax(x,y,equation); 
        
        % Format results
        tempResult1.tau1 = tau1; 
        tempResult1.tau3 = tau3; 
        tempResult1.subject = temp.Subject;
        if listAMPM(jAMPM)=='pm'
            tempResult1.PM = 1; 
        else
            tempResult1.PM = 0; 
        end
        
        % Add result to master 
        result = vertcat(result, tempResult1); 
        
        clear tempResult1 x y temp
    end
end


%% Save
if doSave==1 
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)   
    end
    %writetable(result, fullfile(savePath, sprintf('legSensor_forced3exp_%dms_%dms.csv',tau1,tau3))) % .csv
    save(fullfile(savePath, sprintf('legSensor_forced3exp_%dms_%dms.mat',tau1,tau3)), 'result') % .mat
end
        
  
%% Process MR sensor data into table S8
%  Convert from legSensor_forced3exp_40ms_250ms.mat into summary table of R2. 
%  This table is fed into statics_for_paper.R
%
%  Table columns:   AM | PM | Change and corresponding p-values for each indicator
%                   
%  This table is loaded into R for further statistics (permutation and Welch tests) and then put into supplemental materials

%% Load data
load('/Volumes/cimalablargefiles/Hydration/MRI_Study/ProcessedData/mr_sensor_leg/legSensor_forced3exp_40ms_250ms.mat')

%% Edit 
result.Properties.VariableNames{24} = 'Subject';
result.HD = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1]';
result.AMPM(result.PM==1) = {'PM'};
result.AMPM(result.PM==0) = {'AM'};
result.AMPM = categorical(result.AMPM); 
result.ROI(:) = categorical({'mrSensor'});

mrSensor_R2 = calculate_ttest(result, {'relative_amp2'});


%% To save as csv for Rstudio
mrSensor_R2.tTest{1,1}(end,:) = [];
writetable(mrSensor_R2.tTest{1,1},'/Users/linacolucci/Documents/Thesis/Paper/data_for_R/mrSensor_R2.csv')
    

