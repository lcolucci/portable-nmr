%% Plot Multiexp Results
%   This script generates "line plot" with AM and PM values for each individual subject
%   It is used for results of a multi-exp fit
%   It creates 1 figure with as many subplots as there are elements in 'whichParamsList'
%
%   Function is called in the following scripts, for example,
%   analyze_ROI_multiexp_results.m, visualize_multiexp_results.m, etc. 
%
%   Lina A. Colucci, 2018 

function figRoiResults = plot_multiexp_results(data, whichParamsList,includeText, markerSize, lineWidth, listFaceColors,listEdgeColors) 

% Checks
if strcmp(class(whichParamsList),'char')
    whichParamsList = {whichParamsList}; % convert to cell array
end

% Define aesthetics for {HC,HD} plots
keynoteOrange = [255 147 0]./255;
keynoteBlue = [0 162 255]./255; 
listMarkers = {'o--', '^-'};
if exist('listFaceColors','var')~=1 % 1 if a variable in workspace
    listFaceColors = {'w',keynoteBlue}; 
end
if exist('listEdgeColors','var')~=1
    listEdgeColors = {'k',keynoteBlue}; 
end
if exist('markerSize','var')~=1 
    markerSize = 22; %16
end
if exist('lineWidth','var')~=1
    lineWidth = 2; %1
end
if exist('includeText','var')~=1 
    includeText= 1; 
end
fontSize=32; 


%% Loop through parameters
nParams = length(whichParamsList);
figRoiResults = figure;
for jParam = 1:nParams

    % Pick 1 Parameter
    whichParam = whichParamsList{jParam}; 

    subplot(1,nParams,jParam)
    fnSubjects = unique(data.Subject);
    nSubjects = length(fnSubjects); 

    for kSubject = 1:nSubjects

        % Select 1 Subject
        data.Subject = categorical(data.Subject); 
        dataSubset = data(data.Subject==fnSubjects(kSubject), :); 

        % Select plot colors
        if dataSubset.HD(1)==1
            marker = listMarkers{2}; 
            faceColor = listFaceColors{2}; 
            edgeColor = listEdgeColors{2}; 
        else
            marker = listMarkers{1}; 
            faceColor = listFaceColors{1}; 
            edgeColor = listEdgeColors{1}; 
        end

        % Plot
        try 
        valueAM = dataSubset.(whichParam)(dataSubset.AMPM=='AM',:); 
        valuePM = dataSubset.(whichParam)(dataSubset.AMPM=='PM',:);
        catch
            valueAM = dataSubset.(whichParam)(dataSubset.PM==0,:); 
            valuePM = dataSubset.(whichParam)(dataSubset.PM==1,:); 
        end
        plot([0 1], [valueAM, valuePM], marker,'MarkerFaceColor',faceColor,'MarkerEdgeColor', edgeColor, 'Color',edgeColor, 'MarkerSize',markerSize,'LineWidth',lineWidth); hold on

        % Add 'Subject' Text Label 
        if includeText==1
            text(-.35, dataSubset.(whichParam)(1), char(dataSubset.Subject(1)))
        end
    end
    ylabel(sprintf('%s',whichParam),'Interpreter','none')
    xlim([-.5 1.5])
    xticks([0 1]); xticklabels({'Pre','Post'})
    set(gca, 'FontSize', fontSize);

end
set(gcf, 'color', 'w', 'Position', [73   37  550.*nParams  900])
%mtit(sprintf('%s',whichROI),'fontsize',22,'xoff',0,'yoff',0.03,'Interpreter','none')


% original 
% set(gcf, 'color', 'w','Position',[560    62   640   886]);set(gca, 'FontSize', 16);
