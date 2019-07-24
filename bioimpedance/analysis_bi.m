%% BI Analysis 
%  This script produces scatterplots and boxplots of BI data. 
%
%  Input: 
%        Tidy BI datatable (either BI_Averages.mat or BI_Individual.mat) from previous script (process_bioimpedance.m) 
%  Outputs: 
%        Multiple plots and BI_change datatable, which calculates difference between pre-and-post BI values
%
%       1. Scatter Plots: creates scatter plots of specified columns 
%                         data - tidy datatable to plot
%                         paramYlist - list of column names to plot on y-axis
%                         paramXlist - list of column names to plot on x-axis
%            As many scatterplots will be saved as (# paramYlist) x (# paramXlist)
%       2. 4 Boxplots: Pre vs Post
%            A figure will be created for each dataset specified by 'loopWholeBody' (and 'loopValue'). Each figure has 4 boxplots representing the pre- and post-BI values for both HC and HD participants. 
%       3. 2 Boxplots: The pre-to-post change in BI values for both HC and HD participants
%            Same as above but the pre-to-post change is plotted so there are only 2 boxplots instead of 4
%       4. Data formatting for p-value calculations in R
%            A table that formats BI data for p-value calculations in R (Welch tests and permutation tests; fed into script 'statistics_for_paper.R')
%
%  Lina A. Colucci 


clear; close all

%% ----- USER INPUTS ----
% Plot specifications
posVector = [560   528   1*560   1*420]; % position vector 
markerSize = 16; 
fontSize = 26; %for labels and axes
% Save
doSave = 1; 
savePath = '/Paper/figures/Bioimpedance'; % save path for raw figures 
savePathTidyData = ''; %save path for tidy data (to be uploaded in R for calculating stats)
% ----- BIOIMPEDANCE ----
load('/Volumes/cimalablargefiles/Hydration/MRI_Study/ProcessedData/bioimpedance/BI_Averages.mat')
replaceSecondtoLast = 1; % HC02, HC06, HD01b do not have last BI timePts, only second-to-last ones
                         % 0 if you want to leave those subjects out of the analysis
                         % 1 if you want to use the second-to-last timePt and include those subjects

% ----- DEMOGRAPHICS ----
load('/Volumes/cimalablargefiles/Hydration/MRI_Study/ProcessedData/redcap/demographics.mat')
%% ------------------


%% Load Data
% BI
dataBI = bi; clear bi
dataBI.SubjectName = categorical(dataBI.SubjectName); 
if replaceSecondtoLast==1
    dataBI.TimePt_Last(30) = 1;
    %bi.TimePt_SecondToLast(30) = 0;
    dataBI.TimePt_Last(90) = 1;
    %bi.TimePt_SecondToLast(90) = 0;
    dataBI.TimePt_Last(64) = 1;
    %bi.TimePt_SecondToLast(64) = 0;
end

% Demographics 
dataDemographics = demographics; clear demographics 


%% Calculate BI pre-post difference

% Prep BI data
list_subjects = unique(dataBI.SubjectName); 
BI_change_results = table(); 
BI_results = table(); 
% Need to have this loop, instead of just applying varfun because some
% subjects don't have a 'last' value and so varfun returns an error when that's the case
for jj=1:length(list_subjects)

    % subset data, just 1 subject
    tempData = dataBI(dataBI.SubjectName==list_subjects(jj),:); 
    list_wholeBody = unique(tempData.WholeBody); 
    
    for kk=1:length(list_wholeBody)

        % subset data further (just pre/post for either wholeBody or leg)
        tempData2 = tempData(tempData.WholeBody == list_wholeBody(kk),:); 
        tempData2 = tempData2(tempData2.TimePt_Last == 1 | tempData2.TimePt_First == 1, :); 

        % Calculate pre-post difference
        % positive number means value increased from pre to post. negative number means value decreased from pre to post. 
        tempResultDiff = varfun(@diff, tempData2, 'GroupingVariables',{'SubjectName','WholeBody'},'InputVariables',{'Rinf','Re','Ri','Rcentre','TBW','TBW1','ECF','ECF1','ICF','ICF1','TBW_calfLength','TBW1_calfLength','ECF_calfLength','ECF1_calfLength','ICF_calfLength','ICF1_calfLength'}); 

        % Calculate % difference: (pre-post)/pre*100
        % positive number means value increased from pre to post. negative number means value decreased from pre to post. 
        try 
            tempResultPercent = varfun(@(x) 100.*(x(2)-x(1))/x(1), tempData2, 'GroupingVariables',{'SubjectName','WholeBody'},'InputVariables',{'Rinf','Re','Ri','Rcentre','TBW','TBW1','ECF','ECF1','ICF','ICF1','TBW_calfLength','TBW1_calfLength','ECF_calfLength','ECF1_calfLength','ICF_calfLength','ICF1_calfLength'});
            tempResultPercent.Properties.VariableNames = {'SubjectName','WholeBody','GroupCount','percentDiff_Rinf','percentDiff_Re','percentDiff_Ri','percentDiff_Rcentre','percentDiff_TBW','percentDiff_TBW1','percentDiff_ECF','percentDiff_ECF1','percentDiff_ICF','percentDiff_ICF1','percentDiff_TBW_calfLength','percentDiff_TBW1_calfLength','percentDiff_ECF_calfLength','percentDiff_ECF1_calfLength','percentDiff_ICF_calfLength','percentDiff_ICF1_calfLength'};

            % Merge difference and percentDifference results
            tempResultAggregate = join(tempResultDiff, tempResultPercent,'Keys',{'SubjectName','WholeBody','GroupCount'}); 

            % Aggregate with master
            BI_change_results = vertcat(BI_change_results, tempResultAggregate); 

        catch
            % if subject doesn't have a 'last timePt' value, error will occur
        end

        % Aggregate clean, raw results (not differences/percent differences)
        BI_results = vertcat(BI_results, tempData2); 

        % Delete temp var
         clear tempResultDiff tempResultPercent tempResultAggregate tempData2 whichRow* 
         
    end
    
    clear tempData
end


%% Aggregate BI and demographics
data = join(BI_change_results, dataDemographics,'Keys','SubjectName');

%% Scatter Plot: BI vs Weight Change
% Specify plotting data
paramYlist = {'diff_ECF','percentDiff_ECF','diff_ECF1','percentDiff_ECF1', 'diff_TBW','percentDiff_TBW','diff_TBW1','percentDiff_TBW1','diff_Re','percentDiff_Re','diff_Rinf','percentDiff_Rinf','diff_ECF_calfLength','percentDiff_ECF_calfLength','diff_ECF1_calfLength','percentDiff_ECF1_calfLength', 'diff_TBW_calfLength','percentDiff_TBW_calfLength','diff_TBW1_calfLength','percentDiff_TBW1_calfLength','diff_Re_calfLength','percentDiff_Re_calfLength','diff_Rinf_calfLength','percentDiff_Rinf_calfLength'}; 
paramXlist = {'WeightChange','WeightChangePercent'};

% Plot
listWholeBody = unique(data.WholeBody); 
for ii = 1:length(listWholeBody)
    plotData = data(data.WholeBody==listWholeBody(ii),:); 
    for jj=1:length(paramXlist)
        paramX = paramXlist{jj}; 
        for kk=1:length(paramYlist)
            paramY = paramYlist{kk}; 
            try
                figure
                plot(-1.*plotData.(paramX), plotData.(paramY),'ko','MarkerSize',markerSize,'MarkerFaceColor','k')
                l = lsline; 
                l.LineWidth = 3; 
                xlabel(sprintf('%s',paramXlist{jj}),'Interpreter','none')
                ylabel(sprintf('%s',paramYlist{kk}),'Interpreter','none')
                xl = xlim; yl = ylim;
                textX = xl(1) + abs(xl(1)-xl(2))/2; 
                textY = yl(1) + abs(yl(1)-yl(2))*.9; 
                text(textX, textY, sprintf('r^2 = %.3f', rsquared(-1.*plotData.(paramX), plotData.(paramY))),'FontSize',fontSize)
                set(gcf, 'color', 'w','Position',posVector);
                set(gca, 'fontsize',fontSize,'linewidth',2,'XColor','k','YColor','k'); box off 
                title(sprintf('Whole Body = %d',listWholeBody(ii)))

                % Save
                if doSave==1
                    % make folder if it doesn't exist
                    if ~(7==exist(fullfile(savePath,'ScatterPlots'))) % exist = 7 when it's a folder
                        mkdir(fullfile(savePath,'ScatterPlots'))
                    end
                    saveas(gcf, fullfile(savePath, 'ScatterPlots',sprintf('wholeBody%d_%s_vs_%s.eps',listWholeBody(ii),paramXlist{jj},paramYlist{kk})),'epsc')
                    %saveas(gcf, fullfile(savePath, sprintf('wholeBody%d_%s_vs_%s.png',listWholeBody(ii),paramXlist{jj},paramYlist{kk})),'png')
                end
            
            catch 
                % Whole body measurements won't have '*_calfLength' variables to plot
            end
        end
    end
end


%% ------------ Boxplots -------

%% 4 BOXPLOTS: PRE VS POST
% Looping vars
loopWholeBody = [1 1 0 0 1 1 0 0 0 0]; 
loopValue = {'Re','Rinf','Re','Rinf','ECF','TBW','ECF','TBW','ECF_calfLength','TBW_calfLength'}; 
loopTitle = {'Whole Body R_e', 'Whole Body R_{inf}', 'Leg R_e', 'Leg R_{inf}', 'Whole Body ECF', 'Whole Body TBW', 'Leg ECF', 'Leg TBW', 'Leg ECF CalfLength', 'Leg TBW CalfLength'}; 
loopYlabel = {'Resistance (\Omega)','Resistance (\Omega)','Resistance (\Omega)','Resistance (\Omega)','Fluid (L)','Fluid (L)','Fluid (L)','Fluid (L)','Fluid (L)','Fluid (L)'}; 

for ii=1:length(loopWholeBody)
    
    dataPlot = dataBI(dataBI.WholeBody == loopWholeBody(ii),:);
    dataPlot = dataPlot(dataPlot.TimePt_Last == 1 | dataPlot.TimePt_First == 1, :); 
    colID_HD = find(strcmpi(dataPlot.Properties.VariableNames, 'HD')); 
    colID_TimePtLast = find(strcmpi(dataPlot.Properties.VariableNames, 'TimePt_Last')); 
    dataPlot.uniqueID = findgroups(dataPlot(:,[colID_TimePtLast,colID_HD])); % group by TimePt_Last and HD to have groups in right numerical order 

    customBoxplot4(dataPlot,loopValue{ii},dataPlot.(loopValue{ii}), dataPlot.uniqueID)
    ylabel(sprintf('%s',loopYlabel{ii}))
    title(sprintf('%s',loopTitle{ii}))

    % Save
    if doSave==1
        % make folder if it doesn't exist
        if ~(7==exist(fullfile(savePath,'PrevsPost'))) % exist = 7 when it's a folder
            mkdir(fullfile(savePath,'PrevsPost'))
        end
        saveas(gcf, fullfile(savePath, 'PrevsPost',sprintf('%s.eps',loopTitle{ii})),'epsc')
    end
end


%% 2 BOXPLOTS: CHANGE FOR HC VS HD

%Looping vars
loopWholeBody = [1 1 1 1 0 0 0 0 0 0]; 
loopValue = {'diff_Re', 'diff_Rinf','diff_ECF','diff_TBW',...
             'diff_Re', 'diff_Rinf','diff_ECF','diff_TBW', 'diff_ECF_calfLength','diff_TBW_calfLength'};
loopTitle = {'Whole Body \DeltaR_e','Whole Body \DeltaR_{inf}','Whole Body \DeltaECF','Whole Body \DeltaTBW',...
             'Leg \DeltaR_e','Leg \DeltaR_{inf}','Leg \DeltaECF','Leg \DeltaTBW','Leg Calf Length \DeltaECF','Leg Calf Length \DeltaTBW'};
loopYlabel = {'\DeltaResistance (\Omega)', '\DeltaResistance (\Omega)', '\DeltaFluid (L)','\DeltaFluid (L)',...
             '\DeltaResistance (\Omega)', '\DeltaResistance (\Omega)', '\DeltaFluid (L)','\DeltaFluid (L)', '\DeltaFluid (L)','\DeltaFluid (L)'}; 
    

dataBoxplot = BI_change_results; 
dataBoxplot.uniqueID = contains(string(dataBoxplot.SubjectName),'HD'); 
dataBoxplot.HD = dataBoxplot.uniqueID; 

for ii = 1:length(loopWholeBody)
    
    dataPlot = dataBoxplot(dataBoxplot.WholeBody == loopWholeBody(ii),:); 
    customBoxplot2(dataPlot, loopValue{ii}, dataPlot.(loopValue{ii}), dataPlot.uniqueID);
    ylabel(sprintf('%s',loopYlabel{ii}))
    title(sprintf('%s',loopTitle{ii}))
    
    % Save
    if doSave==1
        % make folder if it doesn't exist
        if ~(7==exist(fullfile(savePath,'Change'))) % exist = 7 when it's a folder
            mkdir(fullfile(savePath,'Change'))
        end
        saveas(gcf, fullfile(savePath, 'Change',sprintf('%s.eps',strrep(loopTitle{ii}, '\Delta','Change in '))),'epsc')
    end
end


%% Format Data for p-Val calculations in R

% Format Change data
change = BI_change_results;
colIDs_percentDiff = find(contains(change.Properties.VariableNames,'percentDiff')); 
change(:,colIDs_percentDiff) = []; % delete 'percentDiff' columns
change(:,'GroupCount') = []; 
change(:,'diff_Rcentre') = [];
change = stack(change,{'diff_Rinf','diff_Re','diff_Ri','diff_TBW','diff_ECF','diff_TBW1','diff_ECF1','diff_ICF','diff_ICF1','diff_TBW_calfLength','diff_ECF_calfLength','diff_TBW1_calfLength','diff_ECF1_calfLength','diff_ICF_calfLength','diff_ICF1_calfLength'},'NewDataVariableName','Change','IndexVariableName','Type');
change.Type = categorical(replace(string(change.Type),'diff_',''));

% Format Pre data
pre = dataBI(dataBI.TimePt_First==1,:); 
pre = pre(:,[1,4,6,24,26,28,36,38,40,42,44,46, 56,58,60,62,64,66]); 
pre = stack(pre, {'Rinf','Re','Ri','TBW','TBW1','ECF','ECF1','ICF','ICF1','ECF_calfLength','ICF_calfLength','TBW_calfLength','ECF1_calfLength','ICF1_calfLength','TBW1_calfLength'},'NewDataVariableName','Pre','IndexVariableName','Type'); 

% Format Post data
post = dataBI(dataBI.TimePt_Last==1,:); 
post = post(:,[1,4,6,24,26,28,36,38,40,42,44,46, 56,58,60,62,64,66]); 
post = stack(post, {'Rinf','Re','Ri','TBW','TBW1','ECF','ECF1','ICF','ICF1','ECF_calfLength','ICF_calfLength','TBW_calfLength','ECF1_calfLength','ICF1_calfLength','TBW1_calfLength'},'NewDataVariableName','Post','IndexVariableName','Type'); 

% Merge data frames
BIdataTidy = join(post, pre); 
BIdataTidy = join(BIdataTidy, change); 
BIdataTidy = BIdataTidy(:,[1:4 6 5 end]); %Put 'Pre' in front of 'Post'


% Save data
if doSave==1
    % make folder if it doesn't exist
    if ~(7==exist(fullfile(savePath))) % exist = 7 when it's a folder
        mkdir(fullfile(savePath))
    end
    save(fullfile(savePath, 'BIdata.mat'),'BIdataTidy')
    writetable(BIdataTidy, fullfile(savePath, 'BIdata.csv'))
end

