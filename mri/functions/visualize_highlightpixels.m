%% Visualize: Highlight Pixels
% INPUT
%       feed it data that's already been checked/filtered/subset as desired
%       data as a long table (rather than structure)
%    
%  3/30/18 LAC - Added capability that if no pixels meet the criteria (i.e. 'data' table is empty), 
%                just plot the image with title '0 pixels'

function figHighlight = visualize_highlightpixels(data, image, timePt,optionalTitle) %color? 

%% Checks  
keynoteBlue = [0 162 255]./255; 
keynoteOrange = [255 147 0]./255;
highlightColor = keynoteOrange; %'cyan'; 
if ~exist('timePt')
    timePt = 2; %visualize timePt #2 unless otherwise specified
end

% This filter makes it so that only the pixels that are contributing to the histogram will be highlighted (instead of the NaN ones)
% try
%         data = data(~(isnan(data.relax1) & isnan(data.relax2)),:); % 2exp fit
% catch
%         data = data(~(isnan(data.relax1)),:); % 1exp fit
% end
        
%% Scale the Image to be from range 0-1 (instead of 0-4096)
imageScaled = image(:,:,:,timePt)./4096; %scale the image to be from 0 to 1 instead of 0 to 4096 (2^12 bit encoding in MRI scanner)
        
%% Loop through slices and plot highlighted pixels
if isempty(data)
    nSlices = 4; % default
else
    nSlices = 4; % default
    sliceList=[1,2,3,4];
%     sliceList = unique(data.slice); 
%     sliceList(isnan(sliceList)) = []; %delete nan values
%     nSlices = length(sliceList); 
end
figHighlight = figure; 
for iSlice=1:nSlices
    if isempty(data) %If no pixels meet this criteria, just plot the image
        subplot(nSlices,1,iSlice)   
        imshow(imageScaled(:,:,iSlice))
        ylabel(sprintf('Slice %d',iSlice),'fontweight','bold','fontsize',12)
        title(sprintf('0 pixels')); set(gca, 'FontSize', 14);
    else
        % Select 1 Slice at a Time
        tempData = data(data.slice==sliceList(iSlice),:);  
        nPixels = height(tempData); % # of pixels in this slice

        % Create Array of Coordinates for Which Pixels to Highlight
        highlightCoordinates = []; %reset for each slice
        for jPixels=1:nPixels                                                       
            highlightCoordinates = vertcat(highlightCoordinates, [tempData.col(jPixels) tempData.row(jPixels) 1 1]); % [col row height width] <-- I want a 1x1 pixel highlight
        end

        % 'insertShape' creates the actual image+highlight object
        annotatedPointsSlice = insertShape(imageScaled(:,:,iSlice), 'FilledRectangle', highlightCoordinates,'Color',highlightColor); 

        subplot(nSlices,1,iSlice)   
        imshow(annotatedPointsSlice); % show image
        ylabel(sprintf('Slice %d',sliceList(iSlice)),'fontweight','bold','fontsize',12)
        title(sprintf('%d pixels', nPixels)); set(gca, 'FontSize', 14);

        clear tempData annotatedPointsSlice
    end
end

% Format Figure
set(gcf, 'color', 'w', 'Position', [607    53   390   900])
if exist('optionalTitle')
    mtit(gcf, sprintf('%s',optionalTitle), 'fontsize',18,'xoff',0,'yoff',0.02);
    % suptitle(
end




 