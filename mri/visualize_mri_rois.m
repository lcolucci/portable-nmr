%% Visualize ROIs on Images

clear
%% ------ USER INPUTS --------
whichIdToAnalyzeLoop = [1];   % Look at config.xlsx to find ID of the subjects you want to analyze
maskList = {'muscle_anterior','muscle_deepposterior','muscle_gastrocnemius','muscle_lateral','muscle_soleus'}; 
overlayAllMasks = 1;            % 0 for separate figure per mask, 1 for all masks on a single figure
whichSlice = 1;                % 99 for plotting all slices (subplots), 1-4 for plotting only a specific slice (no subplots)
% --- path to config and images
pathToMriStudy = '/Volumes/CimaLabLargeFiles/Hydration'; % Path to folder that contains 'MRI_Study' folder (contains raw images and fit results, mri_pixelByPixel_1exp.mat, etc.)
pathToWorkbook = '/GitHub/MRI/mri/config.xlsx'; % Path to config.xlsx 
sheetName = 'MRI_Study_T2';                             % Which sheet name in config.xlsx to load?
whichImageToAnalyze = 'legShorterTE'; % Which image to visualize? (which column of the config.xlsx sheet?) 
timePt = 2; % Which timePt of the image to visualize (timePt 2 is the brightest)? 
% --- save
doSave=1; 
savePath = 'MuscleROIs'; 
%% ----- END USER INPUTS ------

%% Load config excel file
config = load_config_file(pathToWorkbook, sheetName);

%% Set params
keynoteOrange = [255 147 0]./255;
keynoteBlue = [0 162 255]./255; 
scienceRed = [197, 22, 29]./255;


nIDsInLoop = length(whichIdToAnalyzeLoop); 
for iPatient = 1:nIDsInLoop
    
    whichIdToAnalyze = whichIdToAnalyzeLoop(iPatient); 
    whichConfigRow = find(config.id==whichIdToAnalyze);
        
    %% Load image (nifti)
    pathToImageNifti = fullfile(pathToMriStudy, config.pathNiftis{whichConfigRow}, config.(whichImageToAnalyze){whichConfigRow});
    image = load_nifti(pathToImageNifti); image = image.vol; 
    imageScaled = image(:,:,:,timePt)./4096; %scale the image to be from 0 to 1 instead of 0 to 4096 (2^12 bit encoding in MRI scanner)

    %% Loop Through Masks and Slices
    if overlayAllMasks ==1
        figure % One figure for all masks
    end
    nMasks = length(maskList); 
    for jMask = 1:nMasks
        
         %% Load mask (nifti)
         if ~isempty(maskList{1})
            whichMask = maskList{jMask}; 
            pathToMask = fullfile(pathToMriStudy, config.pathMasks{whichConfigRow},strcat(whichMask,'.nii'));
            mask = load_nifti(pathToMask); mask = mask.vol; 
         else
             whichMask = '';
         end

        %%  Visualize
        figure(1)
        %figure(2)
        if whichSlice==99
            nSlices = size(imageScaled,3); 
            for kSlice = 1:nSlices
                figure(1)
                subplot(nSlices,1,kSlice)
                imshow(imageScaled(:,:,kSlice)); hold on
            end
        else
            figure(1)
            imshow(imageScaled(:,:,whichSlice)); hold on
            %figure(2)
            %imshow(imageScaled(:,:,whichSlice))
        end
%         title(sprintf('Slice %d',whichSlice),'fontweight','bold','fontsize',12)
            
        %% Overlay Mask
        if ~isempty(maskList{1})
            figure(1)
            hold on; contour(mask(:,:,whichSlice),'Fill','on')
        end
        figure(1)
        set(gcf, 'color', 'w', 'Position', [221   47 1056  884])
        mtit(gcf, sprintf('%s %s %s',config.subject{whichConfigRow}, config.AMPM{whichConfigRow}, whichMask), 'fontsize',18,'xoff',0,'yoff',0.02,'Interpreter', 'none');
        %figure(2)
        %set(gcf, 'color', 'w', 'Position', [221   47 1056  884])
        
        %% Save
        if doSave==1
            if ~(7==exist(savePath)) % exist = 7 when it's a folder
                mkdir(savePath)
            end
            saveas(figure(1), fullfile(savePath, sprintf('visualizeROI_%s_%s_%s.eps', config.subject{whichConfigRow}, config.AMPM{whichConfigRow}, whichMask)),'epsc')
            %saveas(figure(2), fullfile(savePath, sprintf('visualizeSubject_%s_%s_%s.png', config.subject{whichConfigRow}, config.AMPM{whichConfigRow}, whichMask)))
            close all 
        end
    end
end

