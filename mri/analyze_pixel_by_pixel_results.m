%% Visualize pixel-by-pixel fittings
%   
%   FEATURES
%       + Hightlight pixels of a particular relaxation time or amplitude
%       + Plot histograms of certain values 
%       + Deal with 1exp and 2exp fits 

% Initialization 
clearvars -except result_2exp,%close all
keynoteOrange = [255 147 0]./255;
keynoteBlue = [0 162 255]./255; 
lineThickness = 3; 
allPixels = []; 


%% ------- USER INPUTS -------
% - - - Which Fitting Results?
whichResults = 'result_2exp'; % 'result_1exp' or 'result_2exp'
% - - - Which histograms? 
doHistogram = 1; 
doCDF = 1; 
doDiff = 1; % plot difference on CDF plots? plot integral of CDF diff? 
whichParams ={'relax2'};    % {'relax1','relax2'}; % {'relative_amp2'}; 
histBinEdges = 1:1:1000; %1000;
masterTitle = 'Relax 2 of Whole Leg'; % can be empty '' or a string that will appear in every histogram/CDF/integral title, e.g. 'Relax1 & Relax2' (descriptor of what's being plotted is good idea)
xName = 'Relax 2 (ms)'; % x-label of histograms, 'T2 Relaxation Time (ms)' or 'Relative Amplitude 2 (%)' <<<<<<<<<<<<<<<<<<<<<<<<<
colorsList = {keynoteBlue, keynoteOrange}; %{'b',[.8 .8 .8],[1 .6431 .1098]}; % colors to use in histograms/plots <<<<<<<<<<<<<<<<<<<<<<<<<
histogramSubsetCriteria = 'data.leg_whole==1'; % empty '' or some subsetting criteria, i.e. 'data.muscle_all==1'
% - - - Highlight Pixels 
doHighlight = 0; 
highlightSubsetCriteria = 'data.relax2 >70 & data.relax2<170';%'data.criteria_rmseCutoff==1 | data.criteria_max==1 | data.criteria_min==1 | data.criteria_stdNaN==1'; %i.e. 'data.relax2>170' will turn into: data = data(data.relax2>170, :); %anything in addition to histogramSuset above? 
figureTitle = 'Pixels with Relaxation Times btw 70-170ms'; %'Pixels that were Deleted by Criteria'; 
% - - -
% Config File for Loading Images
pathToMriStudy = '/Volumes/CimaLabLargeFiles/Hydration'; % Path to folder that contains 'MRI_Study' folder (contains raw images and fit results, mri_pixelByPixel_1exp.mat, etc.)
pathToWorkbook = '/mri/config.xlsx'; % Path to config.xlsx 
sheetName = 'MRI_Study_T2';                             % Which sheet name in config.xlsx to load?
whichImageToAnalyze = 'legShorterTE'; % for the highlight pixels on image, which image to use? (which column of the config.xlsx sheet?) 
timePt = 2; % for the highlight pixels on image, which timePt of the image to use? 
% - - - Save
doSave = 0; % 1 to save
savePath = '/Thesis/Figures/MRI/Pixel_by_Pixel/Results_2exp/legSensorROI'; 
%% -------------------------------

fontSize = 22; 

% Load fit results
if strcmpi(whichResults, 'result_1exp')
    load(fullfile(pathToMriStudy, 'MRI_Study/ProcessedData/mri/pixel_by_pixel/mri_pixelByPixel_1exp.mat'))
    results = eval(whichResults); 
elseif strcmpi(whichResults, 'result_2exp')
    load(fullfile(pathToMriStudy, 'MRI_Study/ProcessedData/mri/pixel_by_pixel/mri_pixelByPixel_2exp.mat')); % new results
    results = eval(whichResults); 
%     load('/Volumes/CimaLabLargeFiles/Hydration/MRI_Study/ProcessedData/mri/archive/pixel_by_pixel/mri_pixel_by_pixel_2exp.mat'); % original results
%     results = allResults_2exp; results = results.leg_whole  ; 
else
    error('The ''whichResults'' string does not match either ''result_1exp'' or ''result_2exp''.')
end

% Load config excel file
config = load_config_file(pathToWorkbook, sheetName);

fnSubjects = fieldnames(results); % fn = fieldnames 
nSubjects = length(fnSubjects); 
aggregate = table(); aggregateID = 0; 
for iSubjects = 1:nSubjects 
    fnAMPM = fieldnames(results.(fnSubjects{iSubjects})); 
    nAMPM = length(fnAMPM); 
    
    amPmLoopNumber = 0; 
    for jAMPM = 1:nAMPM
        amPmLoopNumber = amPmLoopNumber + 1; 
        aggregateID = aggregateID + 1; 
        
        %% Select data (1 subject, 1 AMPM) 
        tempResults = results.(fnSubjects{iSubjects}).(fnAMPM{jAMPM}); %structure at 'slice' level
        fnSlices = fieldnames(tempResults); 
        nSlices = length(fnSlices); 
    
        %% Concatenate data from all slices into 1 long-format table
        data = table(); %initialize empty table
        for kSlices = 1:nSlices
            tempData = tempResults.(fnSlices{kSlices}); % data from one slice  
            data = vertcat(data, tempData); 
        end 
        
        % *** IMPORTANT *** The order of the next 3 sections matters. 
        %   1. 'checkfits' must go first because the 99% RMSE cutoff needs to happen with all the dataset for this subject (no subsetting yet. ROI = 'all')
        %   2. 'relative amplitudes' goes after 'checkfits' because otherwise you have to add 'relative_amp1=NaN' and 'relative_amp2=NaN' to all of the criteria in checkfits
        %   3. 'Subset histogram' could switch places with 'relative amplitudes' but if it goes after then it gives ability to subset data based on relative amplitude values
        
        %% Apply fitting criteria (delete pixel results that don't meet criteria)
        data = checkfits(data); %check_1and2exp(data); %   
        
        %% Add 'relative amplitude' components if 2exp fit
        if any(strcmp('relax2',data.Properties.VariableNames)) %if it's a 2exp fit, calculate relative amplitudes
            data.relative_amp1 = 100* data.amp1 ./ nansum([data.amp1, data.amp2],2); 
            data.relative_amp2 = 100* data.amp2 ./ nansum([data.amp1, data.amp2],2);
        end
        
        %% Apply a subset to the data before generating histograms
        if exist('histogramSubsetCriteria') && ~isempty(histogramSubsetCriteria)
            data = data(eval(histogramSubsetCriteria),:);
        end

        %% Straight Normalized Histograms
        if doHistogram ==1
            % HistogramParam: Create vector with the parameter to plot (if more than one param, i.e. {'relax1', 'relax2'}, then vertically concatenate parameter values first)
            if length(whichParams) ==1
                histogramParam = data.(whichParams{1}); 
            else
                histogramParam = []; % <<<<< is there a way to pre-allocate? would it be worth it? have to have another loop. 
                for kParams = 1:length(whichParams)
                    histogramParam = vertcat(histogramParam, data.(whichParams{kParams})); %<<<<<<<<<<< what happens to NaN values? 
                end
            end

            % Plot 
            if amPmLoopNumber==1 
                figHist = figure; 
            end
            figure(figHist)
            hold on;
            histogram(histogramParam, histBinEdges,'FaceColor',colorsList{amPmLoopNumber}, 'Normalization','probability','EdgeColor','none','FaceAlpha',1); %<<<<<<<< MAKE NORMALIZATION TYPE A VARIABLE AND AUTO PUT IT IN TITLE
            xlabel(sprintf('%s', xName))
            ylabel('Sum of Bar Heights = 1') 
            set(gcf, 'color', 'w'); set(gca, 'fontsize',fontSize,'linewidth',4) % Position',[77, 403, 1400, 550]
            title(sprintf('%s %s',fnSubjects{iSubjects}, masterTitle)) %
            legend(fnAMPM); legend boxoff 
%             ylim([0 0.04]); xlim([0 100]); 
            
            if doSave==1  
                if ~(7==exist(savePath)) % exist = 7 when it's a folder
                    mkdir(savePath)
                end
                saveas(gcf, fullfile(savePath, sprintf('histogram_%s.png', fnSubjects{iSubjects})))
                saveas(gcf, fullfile(savePath, sprintf('histogram_%s.svg', fnSubjects{iSubjects})),'svg') %eps format does not support transparency
            end
            
            % Aggregate results
             allPixels = vertcat(allPixels, histogramParam); 
%            prctile(histogramParam,[90, 95, 97.5, 98, 99])
            aggregate.subject(aggregateID) = fnSubjects(iSubjects); 
            aggregate.AMPM(aggregateID) = fnAMPM(jAMPM);
            aggregate.bins(aggregateID) = {histBinEdges}; 
            histogramParam(isnan(histogramParam)) = []; %delete nan values
            aggregate.values(aggregateID) = {histogramParam}; 
        end
        
        %% CDF 
        if doCDF == 1
            % calculate CDF
            startIntegralID = 1; % at which point to start the integration?
            endIntegralID = length(histBinEdges); % where to end the integration? Makes the script more generalizable if I ever want to edit integration bounds. 
            histCDF{amPmLoopNumber} = ecdf2(histogramParam, histBinEdges, startIntegralID, endIntegralID); 
            % plot 
            if amPmLoopNumber ==1
                figCDF = figure; 
            end
            figure(figCDF); hold on; 
            plot(histBinEdges(startIntegralID+1:endIntegralID), histCDF{amPmLoopNumber}, 'Color',colorsList{amPmLoopNumber},'LineWidth',lineThickness)
            if doDiff ==1   % difference btw PM-AM CDFs
                if length(histCDF) == 2
                    diffCDF = histCDF{2} - histCDF{1}; 
                    plot(histBinEdges(startIntegralID+1:endIntegralID), diffCDF,'k','LineWidth',lineThickness)
                    legend({fnAMPM{:},'Difference PM - AM'},'Location','Best'); legend boxoff
                    aggregate.cdfDiff(aggregateID) = {diffCDF}; 
                end
            else
                legend({fnAMPM{:}},'Location','Best'); legend boxoff
            end
            set(gcf, 'color', 'w'); set(gca, 'fontsize',fontSize)
            xlabel(sprintf('%s',xName))
            title(sprintf('%s %s CDF',fnSubjects{iSubjects},masterTitle)) %
            
            aggregate.cdf(aggregateID) = {histCDF{amPmLoopNumber}}; 
            
            %% Integral of CDF Diff
            if doDiff==1
                if exist('diffCDF')
                    integralCdfDiff = cumtrapz(histBinEdges(startIntegralID+1:endIntegralID), diffCDF); 
                    figIntegralDiff = figure; 
                    plot(histBinEdges(startIntegralID+1:endIntegralID), integralCdfDiff,'k','LineWidth',lineThickness)
                    set(gcf, 'color', 'w'); set(gca, 'fontsize',fontSize)
                    xlabel(sprintf('%s',xName))
                    title(sprintf('%s %s\nIntegral of Difference btw AM-PM CDFs',fnSubjects{iSubjects}, masterTitle))

                    if doSave==1
                        saveas(figIntegralDiff, fullfile(savePath, sprintf('integral_cdfDiff_%s.png', fnSubjects{iSubjects})))
                        saveas(figIntegralDiff, fullfile(savePath, sprintf('integral_cdfDiff_%s.eps', fnSubjects{iSubjects})),'epsc')
                    end
                    
                    aggregate.integralDiff(aggregateID) = {integralCdfDiff}; 
                end
            end
            
            if doSave==1  
                if ~(7==exist(savePath)) % exist = 7 when it's a folder
                    mkdir(savePath)
                end
%                 saveas(figCDF, fullfile(savePath, sprintf('cdf_%s.png', fnSubjects{iSubjects})))
                saveas(figCDF, fullfile(savePath, sprintf('cdf_%s.eps', fnSubjects{iSubjects})),'epsc')
            end
            
        end
        
        
        %% Highlight Pixels on Image
        if doHighlight == 1
            % Subset the data (in addition to histogramSubsetCriteria) 
            if exist('highlightSubsetCriteria') && ~isempty(highlightSubsetCriteria)
                data = data(eval(highlightSubsetCriteria),:);
            end

            % Load image (nifti)
            whichConfigRow = find(config.subject==fnSubjects{iSubjects} & config.AMPM == fnAMPM{jAMPM});
            pathToImageNifti = fullfile(pathToMriStudy, config.pathNiftis{whichConfigRow}, config.(whichImageToAnalyze){whichConfigRow});
            image = load_nifti(pathToImageNifti); image = image.vol; 
            % Generate title
            optionalTitle = sprintf('%s %s \n %s',fnSubjects{iSubjects},fnAMPM{jAMPM}, figureTitle); 
            % Highlight pixels on image
            figHighlight = visualize_highlightpixels(data, image, timePt,optionalTitle); 

            % Save
            if doSave==1
                if ~(7==exist(savePath)) % exist = 7 when it's a folder
                    mkdir(savePath)
                end
%                 saveas(gcf, fullfile(savePath, sprintf('highlight_%s_%s.png', fnSubjects{iSubjects}, fnAMPM{jAMPM})))
                saveas(gcf, fullfile(savePath, sprintf('highlight_%s_%s.eps', fnSubjects{iSubjects}, fnAMPM{jAMPM})),'epsc')
            end
        end
        clear whichConfigRow pathToImageNifti image 
    end
    clear histCDF diffCDF figCDF figHist figIntegralDiff integralCdfDiff 
end

%% Overlay of Plots
% convert subject/AMPM to categorical
aggregate.AMPM = categorical(aggregate.AMPM); 
aggregate.subject = categorical(aggregate.subject); 
ids_HD = contains(cellstr(aggregate.subject), 'HDMRI'); 
aggregate.HD(ids_HD) = 1;

%% Overlaid CDF
% overlay of AM CDF
data = aggregate(aggregate.AMPM=='AM',:);
figure1 = figure; 
subplot(1,2,1)
for ii=1:height(data)
    if data.HD(ii)==1, linestyle = '--'; end
    if data.HD(ii)==0, linestyle = '-'; end
    x = data.bins{ii};
    plot(x(2:end), data.cdf{ii},linestyle,'LineWidth',2); hold on; 
end
set(gcf, 'color', 'w'); set(gca, 'fontsize',fontSize)
xlabel(sprintf('%s',xName))
title(sprintf('AM CDF: %s ',masterTitle))
axes1 = axes('Parent',figure1,...
    'Position',[0.530845110771581 0.11 0.334659090909091 0.815]);
hold(axes1,'on');
% overlay of PM CDF
data = aggregate(aggregate.AMPM=='PM',:);
axes2 = subplot(1,2,2);
set(axes2, 'Position', [0.5241    0.1100    0.3347    0.8150])
for ii=1:height(data)
    if data.HD(ii)==1, linestyle = '--'; end
    if data.HD(ii)==0, linestyle = '-'; end
    x = data.bins{ii};
    plot(x(2:end), data.cdf{ii},linestyle,'LineWidth',2); hold on; 
end
legend1 = legend(cellstr(data.subject),'Location','eastoutside'); legend boxoff
set(legend1,'Position',[0.884453781512605 0.263095238095238 0.0945378151260501 0.608333333333333]);
set(gcf, 'color', 'w','Position',[26 534 1190 420]); set(gca, 'fontsize',fontSize)
xlabel(sprintf('%s',xName))
title(sprintf('PM CDF: %s ',masterTitle))
if doSave==1 
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end
    saveas(gcf, fullfile(savePath, 'cdf_All_AMPM.eps'),'epsc')
    saveas(gcf, fullfile(savePath, 'cdf_All_AMPM.png'))
end

%% Overlay of CDF Diff
data = aggregate(aggregate.AMPM=='PM',:);
figure; 
for ii=1:height(data)
    if data.HD(ii)==1, linestyle = '--'; end
    if data.HD(ii)==0, linestyle = '-'; end
    x = data.bins{ii};
    plot(x(2:end), data.cdfDiff{ii},linestyle,'LineWidth',2); hold on; 
end
legend(cellstr(data.subject),'Location','eastoutside'); legend boxoff
set(gcf, 'color', 'w','Position',[35   106   800   420]); set(gca, 'fontsize',fontSize)
xlabel(sprintf('%s',xName))
title(sprintf('AM-PM CDF Difference: %s ',masterTitle))
if doSave==1 %<<<<<<<<<<<<<<<<<<<<< should figure out a way to tell if there are going to be 2 loops, only save on 2nd time 
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end
    saveas(gcf, fullfile(savePath, 'cdfDiff_All.eps'),'epsc')
    saveas(gcf, fullfile(savePath, 'cdfDiff_All.png'))
end

%% overlaid integral of CDF difference
data = aggregate(aggregate.AMPM=='PM',:);
figure; 
for ii=1:height(data)
    x = data.bins{ii}; 
    if data.HD(ii)==1, linestyle = '--'; end
    if data.HD(ii)==0, linestyle = '-'; end
    plot(x(2:end), data.integralDiff{ii},linestyle, 'LineWidth',2); hold on
end
legend(cellstr(data.subject),'Location','eastoutside'); legend boxoff
set(gcf, 'color', 'w','Position',[35   106   800   420]); set(gca, 'fontsize',fontSize)
xlabel(sprintf('%s',xName))
title(sprintf('Integral of CDF Diff.: %s ',masterTitle))
if doSave==1 %<<<<<<<<<<<<<<<<<<<<< should figure out a way to tell if there are going to be 2 loops, only save on 2nd time 
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end
    saveas(gcf, fullfile(savePath, 'integral_cdfDiff_All.eps'),'epsc')
    saveas(gcf, fullfile(savePath, 'integral_cdfDiff_All.png'))
end

% t-test of final integral CDF difference
tempData = aggregate(aggregate.AMPM=='PM' & aggregate.HD==0, :); 
for ii=1:height(tempData)
    tempVector = tempData.integralDiff{ii}; 
    integral_HC(ii) = tempVector(end); 
end
tempData = aggregate(aggregate.AMPM=='PM' & aggregate.HD==1, :); 
for ii=1:height(tempData)
    tempVector = tempData.integralDiff{ii}; 
    integral_HD(ii) = tempVector(end); 
end
[p h] = ttest2(integral_HC, integral_HD)

% Compare to BI 
load('/Volumes/CimaLabLargeFiles/Hydration/MRI_Study/ProcessedData/bioimpedance/BI_Averages.mat') %called 'masterAverages'
% -- Leg BI: Change Pre-Post
data = masterAverages(masterAverages.WholeBody == 0, :); % whole body or leg?
list_subjects = unique(data.SubjectName); 
BIresults = table(); 
for ii=1:length(list_subjects)
    % subset data
    tempDataForOneSubject = data(strcmp(data.SubjectName,list_subjects{ii}),:); 
    
    % Create post_minus_pre variables
    BIresults.Subject(ii) = list_subjects(ii);
    if (sum(tempDataForOneSubject.TimePt_Last==1) + sum(tempDataForOneSubject.TimePt_First==1)) < 2
        BIresults.Rcentre_change(ii) = NaN; 
        BIresults.Rzero_change(ii) = NaN; 
        BIresults.Rinf_change(ii) = NaN;
    else
        BIresults.Rcentre_change(ii) = tempDataForOneSubject.Rcentre(tempDataForOneSubject.TimePt_Last==1) - tempDataForOneSubject.Rcentre(tempDataForOneSubject.TimePt_First==1);  
        BIresults.Rzero_change(ii) = tempDataForOneSubject.Rzero(tempDataForOneSubject.TimePt_Last==1) - tempDataForOneSubject.Rzero(tempDataForOneSubject.TimePt_First==1);  
        BIresults.Rinf_change(ii) = tempDataForOneSubject.Rinf(tempDataForOneSubject.TimePt_Last==1) - tempDataForOneSubject.Rinf(tempDataForOneSubject.TimePt_First==1);
    end
end
BIresults.mri_integral(1:7) = integral_HC'; 
BIresults.mri_integral(8:14) = integral_HD'; 
% Plot ECF
figure; 
plot(BIresults.Rzero_change, BIresults.mri_integral,'o'); lsline
xlabel('Change in Leg BI Rzero (%)'); ylabel(sprintf('%s',xName));
set(gcf, 'color', 'w'); set(gca, 'fontsize',16)
title(sprintf('Leg BI Change in Rzero (ECF): %s ',masterTitle))
text(BIresults.Rzero_change, BIresults.mri_integral, BIresults.Subject)
disp(sprintf('Rzero r2 = %.3f',rsquared(BIresults.Rzero_change, BIresults.mri_integral)))
% Plot TBW
figure; 
plot(BIresults.Rinf_change, BIresults.mri_integral,'o')
xlabel('Change in Leg BI Rinf (%)'); ylabel(sprintf('%s',xName));
set(gcf, 'color', 'w'); set(gca, 'fontsize',16)
title(sprintf('Leg BI Change in Rinf (TBW): %s ',masterTitle))
text(BIresults.Rinf_change, BIresults.mri_integral, BIresults.Subject)
disp(sprintf('Rinf r2 = %.3f',rsquared(BIresults.Rinf_change, BIresults.mri_integral)))

% Compare to Leg Sensor
load('/Volumes/CimaLabLargeFiles/Hydration/MRI_Study/ProcessedData/mr_sensor_leg/legSensor_forced3exp_40ms_250ms.mat'); % called 'result'
legSensorChange = varfun(@diff, result, 'GroupingVariables','subject','InputVariables',{'relax','relative_amp1','relative_amp2','relative_amp3'}); 
legSensorChange.mri_integral(1:7) = integral_HC'; 
legSensorChange.mri_integral(8:14) = integral_HD'; 
figure; 
plot(legSensorChange.diff_relative_amp2, legSensorChange.mri_integral,'o'); lsline
text(legSensorChange.diff_relative_amp2, legSensorChange.mri_integral, legSensorChange.subject)
xlabel('Leg Sensor Change in Relative Amp 2 (%)'), ylabel(sprintf('%s',xName));
set(gcf, 'color', 'w'); set(gca, 'fontsize',16)
title(sprintf('Leg Sensor Change Rel. Amp. 2: %s ',masterTitle))
disp(sprintf('Leg Sensor Change Rel. Amp. 2 r2 = %.3f',rsquared(legSensorChange.diff_relative_amp2, legSensorChange.mri_integral)))


%% --- NOTES ---
%% 3 ways to plot CDF 
%         figure; hold on
%         % 1 - cdfplot - don't know how to extract values 
%         cdfplot(histogramParam); hold on; 
%         % 2 - ecdf - x and f are same length as length of histogramParam, so can't subtract AM and PM from each other
%         [f x] = ecdf(histogramParam); 
%         plot(x, f)
%         % 3 - histcounts and cumtrapz - more work but behaves as desired
%         f2 = histcounts(histogramParam, histBinEdges); 
%         plot(histBinEdges(2:end), cumtrapz(histBinEdges(2:end), f2)./max(cumtrapz(histBinEdges(2:end), f2)))
