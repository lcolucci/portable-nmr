%% Custom function to produce figure with 2 boxplots: CHANGE HC vs HD
%  This is used for MRI data (fig 5) and MR Sensor data (fig 7) 
%
%  SYNTAX
%       either: 
%               customBoxplot2(~ANYTHING~, ~ANYTHING~, x, groups)
%
%       or: 
%               customBoxplot2(data, whichParameter)      ** (if there's a column called 'uniqueID')
%
%  INPUT
%       built to load "change" dataset (i.e. AM-PM difference)
%
%
%  Lina A. Colucci, 2018

function f = customBoxplot2(data, whichParameter, x, groups)

% Prep data
if exist('x','var')~=1 % 1 if a variable in workspace
    x = data.(whichParameter); 
end
if exist('groups','var')~=1
    groups = data.uniqueID; 
end

% Plot
f = figure;
h=axes; 
boxplot(x, groups)
set(gcf, 'color', 'w');
set(gca, 'fontsize',32,'linewidth',2); box off

% customize colors and lines
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects 
boxesHC = [a(2),a(4),a(6),a(8),a(10),a(12),a(14)]; 
set(boxesHC, 'Color','k','LineWidth',3,'LineStyle','-')
boxesHD = [a(1),a(3),a(5),a(7),a(9),a(11),a(13)]; 
set(boxesHD, 'Color','k','LineWidth',3,'LineStyle','-')

% labels
ylabel('\Delta R_2 (%)')
set(gca,'xticklabel',{'HC','HD'},'XColor','k','YColor','k')

% calculate significance
diff_hc = data.(whichParameter)(data.HD==0);
diff_hd = data.(whichParameter)(data.HD==1);
% pVal = permutationTest(diff_hc, diff_hd, 1e5-1, 'exact',1); %% I don't trust that these perm tests are working
[~, pVal] = ttest2(diff_hc, diff_hd, 'VarType','unequal');

% add significance stars to plot
model_series = [quantile(diff_hc, 0.95), quantile(diff_hd, 0.95)]; 
sigstar2(model_series, get(gca,'XTickLabel'),[pVal],0,'boxplot',h)


