%% PROCESS PIXEL DATA TO ROI MEANS 
%  This code takes raw T2 decays for each pixel (pixelT2Decays) and  
%  averages them together into an average T2 decay for each ROI
%
%  PREVIOUS SCRIPTS
%       data comes from script: process_mri.m
%
%  INPUT: 
%       pixelData = structure with the 32-pt decay for each pixel
%         ** Make sure the pixelData has columns showing which pixels appear in each ROI. 
%            If it doesn't first run pixelData through script called "merge_results_with_masks_snr.m"
%
%  OUTPUT: 
%       ROImeans = table or structure that has the 32-pt average decay curve for each ROI for each subject/scan
%                  ILT/multi-exp fits/etc. can be done with these decay curves           
%        
%       Organization of the structure: 
%     - ROImeans
%           + ROIs (leg_whole, muscle_all, etc.)
%               ++ Subjects (HC01, HDMRI04b, etc.)
%                   +++ AM/PM
%                       ++++ Slices (slice1, slice2, slice3, slice4, avg)
%                           +++++ The 32-pt decay as a vector
%
% After drawing more ROIs, run 'merge_results_with_masks_snr.m" and
% this script again to generate master 'ROImeans_table.mat' table

clear 

%% USER INPUTS
pathToMriStudy = '/Volumes/CimaLabLargeFiles/Hydration'; % Path to folder that contains 'MRI_Study' folder (contains pixelT2Decays data) 
doSave = 1; 
savePath = '/Volumes/CimaLabLargeFiles/Hydration/MRI_Study/ProcessedData/mri/ROIs'; 
%% --------

%% Load data
pathToData = 'MRI_Study/ProcessedData/mri/pixel_by_pixel/mri_t2PixelDecays.mat'; 
load(fullfile(pathToMriStudy, pathToData))
data = pixelT2Decays; 

%% Automated Loop 
% Initialize master table/structure
nRows = 3000; 
ROImeans_table = table(cell(nRows,1), cell(nRows,1),nan(nRows,1), nan(nRows,1), cell(nRows,1), cell(nRows,1)); % initialize 3000
ROImeans_table.Properties.VariableNames = {'Subject','AMPM','Slice','AllSlices','ROI','T2Decay'};
ROImeans_struct = struct(); 

% Identify subjects to loop through 
fnSubjects = fieldnames(data); %fn = fieldname
nSubjects = length(fnSubjects);
whichRow = 0; 
for iSubject = 1:nSubjects
    
    fnAMPM = fieldnames(data.(fnSubjects{iSubject})); 
    nAMPM = length(fnAMPM); 
    
    tic        
    for kAMPM=1:nAMPM

        fnSlices = fieldnames(data.(fnSubjects{iSubject}).(fnAMPM{kAMPM}));
        nSlices = length(fnSlices); 
        pixelDecays_allSlices = struct(); 
        
        for lSlice = 1:nSlices
                
            % Find List of ROIs in this slice to loop through 
            listROIs = data.(fnSubjects{iSubject}).(fnAMPM{kAMPM}).(fnSlices{lSlice}).Properties.VariableNames; %list of all column names for this one subject & slice
            decayColumnNumber = find(contains(listROIs, 'decay')); %which column in 'decay' in?
            listROIs = listROIs(decayColumnNumber+1:end); % ROIs appear after the 'decay' column. Delete columns before ROIs so just ROIs left in this list. 
            nROIs = length(listROIs); 
            
            for jROI = 1:nROIs

                % Subset Data: Select just ROI pixels from a single slice               
                data_subset_1Slice = data.(fnSubjects{iSubject}).(fnAMPM{kAMPM}).(fnSlices{lSlice}); %subset just 1 slice
                try 
                    data_subset_1Slice_1ROI = data_subset_1Slice(data_subset_1Slice.(listROIs{jROI})==1,:); % subset just pixels in 1 ROI
                    whichRow = whichRow +1; % only add to row count when the above dataset actually exists
                catch % if this ROI doesn't exist for this subject/slice, go on to next loop iteration
                    continue
                end   
                
                % Aggregate T2 decays from all pixels into a matrix                 
                pixelDecays = cell2mat(data_subset_1Slice_1ROI.decay); % NOTE: There are occasionally empty decays so pixelDecays has fewer rows than data_subset_1Slice_1ROI. This happens when the hand-drawn ROI contains a pixel that does not have a t2Decay associated with it (i.e. it was not processed by the process_mri script)
                
                % Calculate average pixelDecay for this slice/ROI
                pixelDecays_average = nanmean(pixelDecays,1); 
                
                % Add data from this slice to table that will consolidate data from all slices
                if isempty(fieldnames(pixelDecays_allSlices)) ||  ~any(strcmp(fieldnames(pixelDecays_allSlices),listROIs{jROI}))
                    pixelDecays_allSlices.(listROIs{jROI}) = []; % initialize empty matrix if first time adding data to this ROI
                end
                pixelDecays_allSlices.(listROIs{jROI}) = vertcat(pixelDecays_allSlices.(listROIs{jROI}), pixelDecays); %keep a matrix of all decays from all slices in this ROI
                
                % Store T2 Decays from this one slice into master table
                ROImeans_table.Subject(whichRow) = fnSubjects(iSubject);
                ROImeans_table.AMPM(whichRow) = fnAMPM(kAMPM);
                ROImeans_table.Slice(whichRow) = str2double(regexp(fnSlices{lSlice}, '\d+','match')); %fnSlices(lSlice);
                ROImeans_table.AllSlices(whichRow) = 0; 
                ROImeans_table.ROI(whichRow) = listROIs(jROI);
                ROImeans_table.T2Decay(whichRow) = {pixelDecays_average};
                
                % Store results  " "  into master structure
                ROImeans_struct.(listROIs{jROI}).(fnSubjects{iSubject}).(fnAMPM{kAMPM}).(fnSlices{lSlice}) = pixelDecays_average; 
                
                % Delete temporary variables
                clear pixelDecays pixelDecays_average
            end
        end
        
        % Average all slices together - Loop through pixelDecays_allSlices and average data from all slices together 
        listROIs = fieldnames(pixelDecays_allSlices); nROIs = length(listROIs); % these should be identical to last time 'listROIs' and 'nROIs' were calculated but just in case, let's re-calculate them here
        for jROI = 1:nROIs
            
            % Isolate the desired ROI dataset
            tempData = pixelDecays_allSlices.(listROIs{jROI}); 
            
            % Average all pixel decays for all slices
            pixelDecays_allSlices_average = nanmean(tempData,1); 
            
            % Add result to structure
            ROImeans_struct.(listROIs{jROI}).(fnSubjects{iSubject}).(fnAMPM{kAMPM}).avg = pixelDecays_allSlices_average; 
            
            % Add result to table
            whichRow = whichRow +1; 
            ROImeans_table.Subject(whichRow) = fnSubjects(iSubject);
            ROImeans_table.AMPM(whichRow) = fnAMPM(kAMPM);
            ROImeans_table.Slice(whichRow) = 99; % 99 = All slices
            ROImeans_table.AllSlices(whichRow) = 1;    
            ROImeans_table.ROI(whichRow) = listROIs(jROI); 
            ROImeans_table.T2Decay(whichRow) = {pixelDecays_allSlices_average};
            
            % Delete temporary variables 
            clear tempData  pixelDecays_allSlices_average
        end
    end
    toc
end

%% Format ROImeans_table
ROImeans_table(any(ismissing(ROImeans_table),2),:) = []; % delete empty rows (will be based on Slice or AllSlices being NaN) 
ROImeans_table.Subject = categorical(ROImeans_table.Subject); % convert to categorical variables
ROImeans_table.AMPM = categorical(ROImeans_table.AMPM); 
ROImeans_table.ROI = categorical(ROImeans_table.ROI); 

%% Save
if doSave==1
    save(fullfile(savePath,'ROImeans_table.mat'), 'ROImeans_table') %will overwrite the existing data
    save(fullfile(savePath,'ROImeans_struct.mat'), 'ROImeans_struct') %will overwrite the existing data
end

