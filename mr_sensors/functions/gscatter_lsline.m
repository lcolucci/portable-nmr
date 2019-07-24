%% gscatter_lsline
%  This function augments Matlab's built-in 'gscatter' function by adding
%  best-fit linear lines and displaying the equation and r2 value in the legend 
%  for each group of gscatter. 
%
%  There is an optional ability to add a text label to each point in the
%  scatter plot as well
%
%  INPUTS 
%       x = array of x values - REQUIRED
%       y = array of y values (same size as x) - REQUIRED
%       groups = variable that defines how x & y will be grouped - REQUIRED
%                can be a categorical variable, numeric vector, character array, string array, or cell array of character vectors
%                Alternatively, group can be a cell array containing several grouping variables (such as {g1 g2 g3})
%       --- Optional Inputs ---
%       markerColors = character vector (i.e. 'k', 'krbm') or string of 3-column arrays (i.e. [0 0 0.5], {[0 0.5 0.5], [.12, .17, 0]})
%       markerShapes = single or multiple character vector or string of symbols (i.e. 'x','xo*s')
%       markerSizes = single value or array of values (i.e. 10, [10,5,20])
%       xName = x-axis label
%       yName = y-axis label
%       textLabels = character vectors (same size as x & y) that will appear next to each marker
%
%  OUTPUT
%       a scatter plot figure with best fit lines. White background. 16pt font.
%
%  Lina A. Colucci (on top of Matlab's built in gscatter function)

function gscatter_lsline(x, y, groups, markerColors, markerShapes, markerSizes,xName, yName, textLabels)

%% Set defaults if certain inputs are empty
if ~exist('markerColors')
    markerColors = 'k';
end
if ~exist('markerShapes')
    markerShapes = 'o'; 
end
if ~exist('markerSizes')
    markerSizes = 7; 
end
if ~exist('xName')
    xName = ''; 
end
if ~exist('yName')
    yName = ''; 
end
%% Perform Checks

% Add this later

%% Do plotting 
figure; 
gscatter(x, y, groups, markerColors, markerShapes, markerSizes,'on',xName,yName); hold on
legendObject = get(legend); 

% Calculate r2 for each fit
r2 = rsquared_groups(x, y, groups); 

% Plot best fit lines and find equations for them
bestFitLines = lsline; %linear best fit line for each group in gscatter
nGroups = size(bestFitLines,2); %number of lines (one for each group)
for iGroup=1:nGroups
    fits{iGroup} = polyfit(get(bestFitLines(iGroup),'xdata'), get(bestFitLines(iGroup),'ydata'), 1); 
end
% Update the legend
for iGroup = 1:nGroups %<-- need this separate loop b/c of issue with 'fits'. It gets populated in reverse order so whole 'fits' array needs to be generated before updating the legend
    legendObject.String(iGroup) = strcat(legendObject.String(iGroup), sprintf(': %.2fx + %.2f, r2 = %.3f',fits{nGroups+1-iGroup}(1), fits{nGroups+1-iGroup}(2), r2{iGroup})); %update legend entries
                                                                                                         % ^^for some reason, fits get populated in reverse order so need to read them in reverse order to have be in proper location for legend
end
legend(legendObject.String,'Location','Best'); legend boxoff
set(gcf, 'color', 'w');set(gca, 'FontSize', 16);

% Add marker labels (optional) 
if exist('textLabels')
    if ~ischar(char(textLabels))
        error('The ''textLabels'' are not char and cannot be converted to char format.')
    end
    text(x, y,char(textLabels))
else
end

% Set interpreter to 'none' (so that things with '_' don't appear as subscript)
set(0, 'DefaultLegendInterpreter', 'none')
set(0, 'DefaultTextInterpreter','none')
 
% set(groot, 'DefaultAxesTickLabelInterpreter', 'none') <-- globally turn off interpreter
 

