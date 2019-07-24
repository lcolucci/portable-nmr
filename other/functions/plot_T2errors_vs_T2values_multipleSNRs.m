%%  plot_T2errors_vs_T2values_multipleSNRs
%   This script is called in syntheticDecays_explore_T2limits_for_TE_SNR.m
%   It plots T2_error vs. T2 values for multiple different SNRs

% dataStruct = table that has a column 
% dataStructColumn = which column of dataStruct to loop through? 
% xData = which column of dataStruct.dataStructColumn table to plot on x-axis? 't2_div_te'
% yData = which column of dataStruct.dataStructColumn table to plot on y-axis? 'errorT2' or 'errorT2_percent'
% yName = y-axis label. 'T_2 measurement error (%)'
% legendData = data for legend. 
% doSecondXAxis = 1. for yes (will be 't2' below 't2_div_te') 
% xNames = x-axis label. '{'Ratio T2/TE','T2 Relaxation Time (ms)'}'

function [fig, ax, h] = plot_T2errors_vs_T2values_multipleSNRs(dataStruct, dataStructColumn, xData, yData, yName, legendData, doSecondXAxis, xNames)

%% Checks and Defaults
if ~exist('xData') || isempty(xData)
    xData = 't2_div_te'; 
end
if ~exist('yData') || isempty(yData)
    xData = 'errorT2_percent'; 
end
if ~exist('xNames') || isempty(xNames)
    xNames = '{''Ratio T2/TE'',''T2 Relaxation Time (ms)''}'; 
end
if ~exist('yName') || isempty(yName)
    yName = 'T_2 measurement error (%)'; 
end
if ~exist('legendData') || isempty(legendData)
    legendData = dataStruct.TargetSNR; 
end
if ~exist('doSecondAxis') || isempty(doSecondAxis)
    doSecondXAxis = 1; 
end

dataToLoop = dataStruct.(dataStructColumn); 

markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
nRows = size(dataToLoop,1); 
fig = figure; 
for ii=1:nRows
    hold on
    plot(dataToLoop{ii,1}.(xData), dataToLoop{ii,1}.(yData), strcat(markers{ii},'-'),'LineWidth',1.5,'MarkerSize',12);  
end
set(gcf, 'color', 'w','Position',[64 49 1124 906]); set(gca, 'fontsize',16)
legend(cellstr(num2str(legendData, 'SNR=%-d'))); legend boxoff    
ax = xlabel(eval(xNames));
ylabel(sprintf('%s',yName)); 

if doSecondXAxis == 1 %default is to plot T2 relaxation time 
    % Add an additional x-label below current one
    xTicks = get(gca, 'xtick');
    yTicks = get(gca,'ytick');
    minY = min(yTicks);
    if strcmp(yData,'errorT2')
        VerticalOffset = .6;
        HorizontalOffset = 1;
    end
    if strcmp(yData,'errorT2_percent')
        VerticalOffset = 1;
        HorizontalOffset = 1;
    end
    for xx = 1:length(xTicks)
        h(xx) = text(xTicks(xx)-HorizontalOffset, minY - VerticalOffset, sprintf('%.0fms',dataToLoop{1,1}.TE(1) .* xTicks(xx)), 'Fontsize',16); % Create a text box at every Tick label position
    end
    set(ax, 'Units','Normalized','position', [0.5 -.06 0]);  % shift the y label to the left by 20
%     set(ax, 'position', get(ax,'position')-[0,.65,0]);  % shift the y label down by 0.65
                                                      % If I ever want to work with Normalized units: set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]); %https://stackoverflow.com/questions/14611259/how-to-adjust-the-distance-between-the-y-label-and-the-y-axis-in-matlab
end