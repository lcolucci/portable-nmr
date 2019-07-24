%% Grouped Bar Graph with One-Sided Error Bars
%  Based on accepted answer to question: https://www.mathworks.com/matlabcentral/answers/102220-how-do-i-place-errorbars-on-my-grouped-bar-graph-using-function-errorbar-in-matlab

function [f,h] = bar_with_error(model_series, model_error, names, xLab, yLab, legendLabs)
% define color
scienceRed = [197, 22, 29]./255;

% Creating axes and the bar graph
f = figure; 
ax = axes;
h = bar(model_series,'BarWidth',1);

% Set custom color for each bar face if desired
h(1).FaceColor = 'white';
h(1).LineWidth = 2; 
h(2).FaceColor = scienceRed; %'yellow';
%h(2).LineStyle = 'none'; 
h(2).LineWidth = 2;

% Properties of the bar graph as required
ax.YGrid = 'on';
ax.GridLineStyle = '-';
ax.XTickLabelRotation = 45; 
xticks(ax,[1:size(model_series,1)]);

% Naming each of the bar groups
xticklabels(ax,names);
set(gca,'TickLabelInterpreter','none')

% X and Y labels
xlabel(sprintf('%s',xLab));
ylabel(sprintf('%s',yLab)); %'Change in Relative Amplitude 2 (%)');

% Creating a legend and placing it outside the bar plot
lg = legend(legendLabs,'AutoUpdate','off'); %'HC','HD'
lg.Location = 'Best';
lg.Orientation = 'Horizontal';
legend boxoff
hold on;

% Finding the number of groups and the number of bars in each group
ngroups = size(model_series, 1);
nbars = size(model_series, 2);

% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));

% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    for j=1:ngroups
        if model_series(j,i) < 0
            errorbar(x(j), model_series(j,i), model_error(j,i),0, 'k', 'linestyle', 'none','linewidth',2);
        else
            errorbar(x(j), model_series(j,i), 0, model_error(j,i), 'k', 'linestyle', 'none','linewidth',2);
        end
    end
end

% General figure formatting
set(gcf, 'color', 'w','Position',[80, 150, 1427, 805]); 
set(gca, 'fontsize',32,'linewidth',2); box off

