%% ---------- Generate Master Data for Bar Graph (pixel-by-pixel results) ---------------
%  This script processes raw pixelwise results into summary tables:
%       aggregate - contains all possible pixelwise information in oneplace
%       master - summarizes difference between AM-PM for HC and HD populations (i.e. the last value of the integralDiff vector)
%  Note that the checkfits(data) function is applied to data before calculating any of the results. 
%
%  Script is adapted from analyze_pixel_by_pixel_results.m (edited to not generate any figures)
%   
%  PREVIOUS SCRIPT
%       data comes from: process_mri.m 
%
%  INPUTS
%       whichResults - loads in either mri_pixelByPixel_2exp.mat or mri_pixelByPixel_1exp.mat
%                   data comes from: process_mri.m 
%       whichParams - which indicators to process? all or a subset of: relax1, relax2, relative_amp2
%       whichHistBinEdges - for the chosen indicators, what are histogram bins that should be used? 
%                   For relaxation times 1-1000ms by 1ms increments is good. 
%                   For relative_amplitude, 1-100% by 1% increments is good. 
%       whichROIs - which ROIs to generate results for? 
%       pathToMriStudy - path to folder that contains 'MRI_Study' folder, which contains all the data
%       doSave = 1 to save plots as .eps
%       savePath = if doSave=1, specify path to which figures will be saved
%   
%   OUTPUTS
%       aggregate - saved as mri_pixelByPixel_[1|2]exp_summary_all.mat
%                   massive master table containig all possible info. Has the following columns: 
%                       Subject - subject name. HC01, HDMRI02b, etc.
%                       HD - 0 or 1
%                       ROI - name of ROI. leg_whole, muscle_lateral, subcu_all, etc. 
%                       Parameter - relax1, relax2, or relative_amp2
%                       AMPM - AM or PM
%                       Bins - user-specified vector of bins (i.e. the x-axis if plotting cdf, cdfDiff, integralDiff)
%                              the vector that would be fed into histogram(data, bins) to generate a histogram with the specified bin edges
%                       Values - array of all pixel values 
%                       cdf - the result of taking the integral of histogram(data, bins)
%                             ecdf2(Values, Bins, Bins(1), Bins(end))
%                       cdfDiff - the difference between AM and PM cdfs for a particular subject/ROI/parameter 
%                                 cdf(AM) - cdf(PM) 
%                       integralDiff - the integral of cdfDiff
%                                      cumtrapz(Bins(2:end), cdfDiff)
%   
%       master - saved as mri_pixelByPixel_[1|2]exp_summary_AMPMchange.mat
%                summary table from 'aggregate' that compares the last integralDiff value of HC vs HD
%                groups. The last integralDiff value is the total AM-PM change.
%                This is the only single value that makes sense to compare for pixel-wise data. 
%
%                Table the following columns:
%                       ROI - name of ROI. leg_whole, muscle_lateral, subcu_all, etc. 
%                       Parameter - relax1, relax2, or relative_amp2
%                       Mean_HC - average of last integralDiff values for HC subjects
%                       Mean_HD -          "        "                 for HD subjects
%                       Std_HC - stdev of last integralDiff values for HC subjects
%                       Std_HD -       "             "             for HD subjects
%                       pVals - two-sample, two-tailed t-test with equal variances btw HC and HD subjects
%                               ttest2(Mean_HC, Mean_HD)
%                       pVals_unequalVar - two-sample, two-tailed t-test with unequal variances btw HC and HD subjects
%                               ttest2(Mean_HC, Mean_HD,'Vartype','unequal')
%                       pVals_permTest - two-sample, two-tailed exact permutation test with 10e5-1 repetitions
%                               permutationTest(Mean_HC, Mean_HD, 1e5-1, 'exact',1)
%                               *** NOTE: I downloaded this function off of Mathworks file exchange ***
%                               https://www.mathworks.com/matlabcentral/fileexchange/63276-permutation-test
%                               I checked its results against R's perm package: permTS(HC, HD, method='exact.mc', control=permControl(nmc=10^5-1)), which I verified with the statistician
%                               I got nearly identical results, though sometimes the thousandth decimal value would be different
%                               I figured this was good enough, though note that it is a file-exchange function and not a built-in, validated function
%
%
%
% Lina A. Colucci 2018 



% Initialization 
clear

%% ------- USER INPUTS -------
% - - - Which Fitting Results?
whichResults = 'result_1exp'; % 'result_1exp' or 'result_2exp'
% - - - Which histograms? 
whichParams ={'relax1'} %,'relax2','relative_amp2'};    % {'relax1','relax2'}; % {'relative_amp2'}; 
whichHistBinEdges = {[1:1:1000],[1:1:1000],[1:1:100]}; 
whichROIs = {'leg_whole','bone','marrow','subcu_all','muscle_all','muscle_anterior','muscle_deepposterior','muscle_gastrocnemius','muscle_lateral','muscle_soleus'}; 
% - - -
% Paths to data 
pathToMriStudy = '/Volumes/CimaLabLargeFiles/Hydration'; % Path to folder that contains 'MRI_Study' folder (contains raw images and fit results, mri_pixelByPixel_1exp.mat, etc.)
% - - - Save
doSave = 1; % 1 to save
savePath = '/Volumes/CimaLabLargeFiles/Hydration/MRI_Study/ProcessedData/mri/pixel_by_pixel'; 
%% -------------------------------

fontSize = 22; 

% Turn off warning that rows are being added 
%       Warning: The assignment added rows to the table, but did not assign values to all of the table's
%       existing variables. Those variables have been extended with rows containing default values. 
warning('off','MATLAB:table:RowsAddedExistingVars')

% Load fit results
if strcmpi(whichResults, 'result_1exp')
    load(fullfile(pathToMriStudy, 'MRI_Study/ProcessedData/mri/pixel_by_pixel/mri_pixelByPixel_1exp.mat'))
    results = eval(whichResults); 
elseif strcmpi(whichResults, 'result_2exp')
    load(fullfile(pathToMriStudy, 'MRI_Study/ProcessedData/mri/pixel_by_pixel/mri_pixelByPixel_2exp.mat')); % new results
    results = eval(whichResults); 
else
    error('The ''whichResults'' string does not match either ''result_1exp'' or ''result_2exp''.')
end


%% Generate 'aggregate' master dataset 
%  If 'aggregate' exists, skip this step 
fnSubjects = fieldnames(results); % fn = fieldnames 
nSubjects = length(fnSubjects); 
aggregate = table(); aggregateID = 0; 
nParams = length(whichParams); 
nROIs = length(whichROIs); 
for iParam = 1:nParams
    param = string(whichParams{iParam}); 
    histBinEdges = whichHistBinEdges{iParam}; 
    
    for iROI = 1:nROIs
        roi = string(whichROIs{iROI}); 

        for iSubjects = 1:nSubjects 
            fnAMPM = fieldnames(results.(fnSubjects{iSubjects})); 
            nAMPM = length(fnAMPM); 

            amPmLoopNumber = 0; 
            for jAMPM = 1:nAMPM
                try 
                amPmLoopNumber = amPmLoopNumber + 1; 
                aggregateID = aggregateID + 1; 
                
                % Build 'aggregate' results table 
                aggregate.Subject(aggregateID) = categorical(fnSubjects(iSubjects));
                if contains(fnSubjects(iSubjects), 'HC')
                    aggregate.HD(aggregateID) = 0; 
                else
                    aggregate.HD(aggregateID) = 1;
                end
                aggregate.ROI(aggregateID) = categorical(roi); 
                aggregate.Parameter(aggregateID) = categorical(param); 
                aggregate.AMPM(aggregateID) = categorical(fnAMPM(jAMPM));
                aggregate.Bins(aggregateID) = {histBinEdges}; 

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

                %% Select ROI subset to the data before generating histograms
                data = data(data.(char(roi))==1,:);

                %% Straight Normalized Histograms
                    % HistogramParam: Create vector with the parameter to plot (if more than one param, i.e. {'relax1', 'relax2'}, then vertically concatenate parameter values first)
                    if length(param) ==1 %************ FIGURE THIS OUT LATER. WILL I EVER WANT TO JOIN RELAX1 AND RELAX2? **********
                        histogramParam = data.(char(param)); 
                    else
                        histogramParam = []; % <<<<< is there a way to pre-allocate? would it be worth it? have to have another loop. 
                        for kParams = 1:length(param)
                            histogramParam = vertcat(histogramParam, data.(char(param{kParams}))); %<<<<<<<<<<< what happens to NaN values? 
                        end
                    end
                    histogramParam(isnan(histogramParam)) = []; %delete nan values
        
                    % Put into master structure 
                    aggregate.Values(aggregateID) = {histogramParam}; 

                %% CDF 
                    % calculate CDF
                    startIntegralID = 1; % at which point to start the integration?
                    endIntegralID = length(histBinEdges); % where to end the integration? Makes the script more generalizable if I ever want to edit integration bounds. 
                    histCDF{amPmLoopNumber} = ecdf2(histogramParam, histBinEdges, startIntegralID, endIntegralID); 
                    if length(histCDF) == 2
                        diffCDF = histCDF{2} - histCDF{1}; 
                        aggregate.cdfDiff(aggregateID) = {diffCDF}; 
                    end   
                    aggregate.cdf(aggregateID) = {histCDF{amPmLoopNumber}}; 

                    %% Integral of CDF Diff
                    if exist('diffCDF')
                        integralCdfDiff = cumtrapz(histBinEdges(startIntegralID+1:endIntegralID), diffCDF); 
                        aggregate.integralDiff(aggregateID) = {integralCdfDiff}; 
                    end 
                    
                catch
                end
                
                clear histogramParam
                
                end
                % Delete AM/PM temp variables 
                clear integralCdfDiff histCDF diffCDF 
                
                clear whichConfigRow pathToImageNifti image 
                
            end
        clear histCDF diffCDF figCDF figHist figIntegralDiff integralCdfDiff 

    end
end

% Delete empty rows in 'aggregate'
aggregate(cellfun(@(x) isempty(x), aggregate.Values),:) = []; 

% Save 'aggregate'
if doSave==1 
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end
    if strcmpi(whichResults, 'result_1exp')
        save(fullfile(savePath, 'mri_pixelByPixel_1exp_summary_all.mat'),'aggregate')
    end
    if strcmpi(whichResults, 'result_2exp')
        save(fullfile(savePath, 'mri_pixelByPixel_2exp_summary_all.mat'),'aggregate')
    end
end

%% Generate 'indiv' dataset 
%  Contains the last integral value for each subject and parameter
%  (pixel-by-pixel quantification of AM-to-PM change) 
clearvars -except aggregate doSave savePath whichResults results

data = aggregate(aggregate.AMPM=='PM',:);
indiv = data(:,1:4); 
indiv.Value(:) = NaN; 
for i=1:height(data) 
    indiv.Value(i) = data.integralDiff{i,1}(end); 
end

% Calculate statistics HC vs HD
ttestResults = table(); 
whichRow = 0; 
roiList = unique(indiv.ROI); 
nROIs = length(roiList); 
for iROI = 1:nROIs
    tempData = indiv(indiv.ROI==roiList(iROI),:); 
    paramList = unique(tempData.Parameter); 
    nParams = length(paramList); 
    for jParam = 1:nParams
        whichRow = whichRow + 1; 
        tempData2 = tempData(tempData.Parameter==paramList(jParam),:); 
   
        try 
            ttestResults.ROI(whichRow) = roiList(iROI); 
            ttestResults.Parameter(whichRow) = paramList(jParam); 
            [tstat, ttestResults.pVals(whichRow)] = ttest2(tempData2.Value(tempData2.HD==0), tempData2.Value(tempData2.HD==1)); 
            [tstat, ttestResults.pVals_unequalVar(whichRow)] = ttest2(tempData2.Value(tempData2.HD==0), tempData2.Value(tempData2.HD==1),'Vartype','unequal'); 
            [ttestResults.pVals_permTest(whichRow), observeddifference, effectsize] = permutationTest(tempData2.Value(tempData2.HD==0), tempData2.Value(tempData2.HD==1), 1e5-1, 'exact',1);
        catch
        end
    end
end

%% master average dataset 
avgMean = table();
avgStd = table(); 
whichRow = 0; 
hdList = unique(indiv.HD); 
for iHD = 1:length(hdList)
    tempData = indiv(indiv.HD==hdList(iHD),:); 
    roiList = unique(tempData.ROI); 
    nROIs = length(roiList); 
    
    for jROI = 1:nROIs
        tempData2 = tempData(tempData.ROI==roiList(jROI),:); 
        paramList = unique(tempData2.Parameter); 
        nParams = length(paramList); 
        
        for kParam = 1:nParams
            whichRow = whichRow + 1; 
            
            tempData3 = tempData2(tempData2.Parameter==paramList(kParam),:); 
            tempValue = mean(tempData3.Value); 
            tempStd = std(tempData3.Value); 
            
            avgMean.HD(whichRow) = hdList(iHD); 
            avgMean.ROI(whichRow) = roiList(jROI); 
            avgMean.Parameter(whichRow) = paramList(kParam); 
            avgMean.Mean(whichRow) = tempValue; 
            
            avgStd.HD(whichRow) = hdList(iHD); 
            avgStd.ROI(whichRow) = roiList(jROI); 
            avgStd.Parameter(whichRow) = paramList(kParam);
            avgStd.Std(whichRow) = tempStd; 
            
            % clear temporary variables
            clear tempValue tempStd tempData3
        end
        clear tempData2
    end
    clear tempData 
end

% Unstack data so that 'HC' and 'HD' have separate columns
avgMean = unstackdata(avgMean, 'Mean');
avgStd = unstackdata(avgStd, 'Std');

% Join mean and stdev values
avg = join(avgMean, avgStd); 

% Join values with p-Values
master = join(avg, ttestResults); 

% Save 'master' dataset from which bargraphs will be created
if doSave==1 
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end
    if strcmpi(whichResults, 'result_1exp')
        save(fullfile(savePath, 'mri_pixelByPixel_1exp_summary_AMPMchange.mat'),'master')
    end
    if strcmpi(whichResults, 'result_2exp')
        save(fullfile(savePath, 'mri_pixelByPixel_2exp_summary_AMPMchange.mat'),'master')
    end
end

%% FUNCTIONS

% unstackdata 
%       unstacks data so that HC and HD values have separate columns
function unstackedData = unstackdata(data, valuename)        
    unstackedData = unstack(data, valuename,'HD');
    unstackedData.Properties.VariableNames = {'ROI','Parameter',strcat(valuename,'_HC'),strcat(valuename,'_HD')}; 
end
