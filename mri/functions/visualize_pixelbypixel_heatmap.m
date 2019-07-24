%% visualize pixel-by-pixel heatmap
%  This function takes a fitting result for a single subject and plots heatmaps of the specified
%  parameters. It creates a separate subplot for each slice. 
%
%  There is a wrapper function (visualize_heatmaps.m) that loops this
%  function across multiple subjects and saves the resulting figures. 
%
%   INPUTS
%       fitData - 1exp or 2exp fit results (structures) at the level where there the only sub-structures are 'slices' 
%                 fitData
%                       slice1
%                       ...
%                       slice4
%                           < table of results including relax1, r2, etc. >
%       whichParamsToPlot - names of params that will be visualized as heatmaps (up to 2) - (cell arrays) i.e. {'relax1','r2'}
%       --- OPTIONAL INPUTS ---
%       colormaps - STRONGLY ENCOURAGED to specify this parameter. Cell array of 1 or 2 colormaps
%       optionalTitle - title of the figure
%       colorlimits - cell array of 1 or 2 vectors specifying lower + upper bounds of the colorbars, i.e. {[0 150], [0 1000]}
%
%   OUTPUTS
%       heatmap - a figure plotting a heatmap of 'whichParamsToPlot' for every slice in the image
%                 min figure subplot dimensions = 1x1 if 1 parameter and 1 slice exist
%                 max figure subplot dimensions = nx2 if 2 parameters and n slices exist
%
%
%   FUTURE TO DO'S
%   + should fitData already have checks applied to it??? or should it be applied here? 
%       

function heatmap = visualize_pixelbypixel_heatmap(fitData, whichParamsToPlot, colormaps, optionalTitle, colorlimits)                                               

%% Checks
if ~iscell(whichParamsToPlot) || ~ischar(whichParamsToPlot{1})
    error('whichParamsToPlot must be a cell or cell array of characters.')
end
if length(whichParamsToPlot)==1
    whichParam1 = whichParamsToPlot{1}; 
    nCols = 1; 
elseif length(whichParamsToPlot)==2
    whichParam1 = whichParamsToPlot{1}; 
    whichParam2 = whichParamsToPlot{2}; 
    nCols = 2; 
elseif length(whichParamsToPlot)==3
    whichParam1 = whichParamsToPlot{1}; 
    whichParam2 = whichParamsToPlot{2}; 
    whichParam3 = whichParamsToPlot{3}; 
    nCols = 3; 
else
    error('whichParamsToPlot must have a length of 1, 2, or 3.')
end
if ~isempty('colormaps')
    colormap1 = colormaps{1};
    if length(colormaps)==2
        colormap2 = colormaps{2}; 
    end
else
    colormap1 = colormap('parula'); 
    colormap2 = colormap('jet'); 
end
if exist('colorlimits') && isempty(colorlimits)
    clear colorlimits
end    
%% Loop through slices and plot heatmaps
fnSlices = fieldnames(fitData); % fn = fieldnames
nSlices = length(fnSlices); 
heatmap = figure; 

for iSlices = 1:nSlices 
    
    % Select data for just one slice
    tempData = fitData.(fnSlices{iSlices}); 
    % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< apply checks???
%     tempData = checkfits(tempData); 
%     tempData = tempData(tempData.leg_whole==1,:); 
%     
    % Calculate relative_amp2
    try
    tempData.relative_amp1 = tempData.amp1 ./ (tempData.amp1 + tempData.amp2) .* 100; 
    tempData.relative_amp2 = tempData.amp2 ./ (tempData.amp1 + tempData.amp2) .* 100; 
    catch
        % if amp2 doesn't exist
    end
    
    % Populate matrix with heatmap data
    heatmap1 = nan(192,144); %nRow, nCol); 
    heatmap2 = nan(192,144); %nRow, nCol);
    for jPixels = 1:height(tempData)
        row = tempData.row(jPixels); 
        col = tempData.col(jPixels); 
        heatmap1(row, col) = tempData.(whichParam1)(jPixels); 
        if nCols==2
            heatmap2(row, col) = tempData.(whichParam2)(jPixels);
        end
    end

    % set colorbar limits (lower and upper bounds) 
    if exist('colorlimits')
        colorlimits1 = colorlimits{1}; 
        if length(colorlimits)==2
            colorlimits2 = colorlimits{2}; 
        end
    end
    if ~exist('colorlimits')
        if any(strcmp('relax2',tempData.Properties.VariableNames)) % default caxis for visualizing 2exp results
            colorlimits1 = [0 150]; 
            colorlimits2 = [0 1000]; 
        else
        	colorlimits1 = [0 500]; % default caxis for visualizing 1exp results. This looks good at either [0 500] or [0 1000] 
        end
    end

    % Plot
    % 1st parameter
    if nCols==1
        subplot(nSlices,1,iSlices)
        imagesc(heatmap1); 
        colorbar; colormap(gca, colormap1); caxis(colorlimits1)
    elseif nCols==2
        subplot(nSlices,2,2*iSlices-1)
        imagesc(heatmap1); 
        colorbar; colormap(gca, colormap1); caxis(colorlimits1)
    elseif nCols==3
        subplot(nSlices,3,3*iSlices-1)
        imagesc(heatmap1); 
        colorbar; colormap(gca, colormap1); caxis(colorlimits1)
    end
    title(sprintf('%s',whichParam1),'Interpreter','none'); 
    ylabel(sprintf('%s \n ',fnSlices{iSlices}),'fontweight','bold','fontsize',12)
    set(gca,'DataAspectRatio',[1 1 1]); set(gca,'ytick',[],'xtick',[]); set(gca, 'FontSize', 12);
    
    %2nd parameter
    if nCols==2
        subplot(nSlices,2,2*iSlices)
        imagesc(heatmap2); 
        title(sprintf('%s',whichParam2),'Interpreter','none'); 
        colormap(gca, colormap2); caxis(colorlimits2)
        colorbar; set(gca,'DataAspectRatio',[1 1 1]); set(gca,'ytick',[],'xtick',[]); set(gca, 'FontSize', 12);
    elseif nCols==3
        subplot(nSlices,3,3*iSlices)
        imagesc(heatmap2); 
        title(sprintf('%s',whichParam2),'Interpreter','none'); 
        colormap(gca, colormap2); caxis(colorlimits2)
        colorbar; set(gca,'DataAspectRatio',[1 1 1]); set(gca,'ytick',[],'xtick',[]); set(gca, 'FontSize', 12);
    end
    
    %3rd parameter
    if nCols==3
        subplot(nSlices,3,3*iSlices)
        imagesc(heatmap3); 
        title(sprintf('%s',whichParam3),'Interpreter','none'); 
        colormap(gca, colormap2); caxis(colorlimits2)
        colorbar; set(gca,'DataAspectRatio',[1 1 1]); set(gca,'ytick',[],'xtick',[]); set(gca, 'FontSize', 12);
    end
    
end

% Format figure
set(gcf, 'color', 'w'); 

if nCols==1
    set(gcf, 'Position', [607    53   394   902])
elseif nCols==2
    set(gcf, 'Position',[97    53   482   902])
end

if exist('optionalTitle')
%     suptitle(sprintf('%s',optionalTitle))
    mtit(gcf, sprintf('%s',optionalTitle), 'fontsize',18,'xoff',0,'yoff',.05);
end

