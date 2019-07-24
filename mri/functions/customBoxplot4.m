%% Custom function to produce figure with 4 boxplots: PRE HC vs HD, POST HC vs HD
%  This is used for MRI data (fig 5) and MR Sensor data (fig 7) 
%
%  SYNTAX
%       either: 
%               customBoxplot4(~ANYTHING~, ~ANYTHING~, x, groups)
%
%       or: 
%               customBoxplot4(data, whichParameter)      ** (if column is called 'uniqueID')
%
% INPUTS
%       data - table that contains data for boxplots. 
%              if x and groups ARE specified then, it's ignored
%                       i.e. customBoxplot4(~, ~, x, groups)
%              if x and groups are NOT specified, then data must have columns called 'whichParameter' and 'uniqueID'
%                       i.e. customBoxplot4(data, whichParameter)
%
%       whichParameter - name of column in data that will become 'x' vector
%
%       x - optional. vector of data points to visualize as boxplot
%
%       groups - optional. vector (same size as x) with which data point belongs in which group
%                the boxplots will appear in order of these IDs
%                it should be: 
%                   am_hc = 1 
%                   am_hd = 2 
%                   pm_hc = 3 
%                   pm_hd = 4 
%                uniqueID = findgroups(data(:,[PM,HD])) - columns should be AM = first (either 0 or AM), HC = first (either 0 or HC) 
%                                                         (when ordered alphanumerically)          
%                
%                
% OUTPUTS
%       f - handle to figure with 4 total boxplots: 2 groups (AM and PM), and 2 boxplots within each group (HC and HD)
%           NOTE: significance stars are based on Matlab's Welch test. Need to do Permutation Test in R for the true pVal 
%                 (usually doesn't change the # of stars)
%
% Lina A. Colucci, 2018
%

function f = customBoxplot4(data, whichParameter, x, groups)

% Prep data
if exist('x','var')~=1 % 1 if a variable in workspace
    x = data.(whichParameter); 
end
if exist('groups','var')~=1
    groups = data.uniqueID; 
end
positions = [1 1.2 1.6 1.8]; % 3 3.5]; 

% Plot
f = figure; 
h = axes; 
boxplot(x, groups, 'positions', positions); 
set(gca, 'xtick', [mean(positions(1:2)) mean(positions(3:4))], 'fontsize',32,'linewidth',2); box off
set(gcf, 'color', 'w');

%customize colors and lines
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects 
boxesHC = [a(12),a(10),a(6),a(8),a(28),a(26),a(24),a(22),a(20),a(18),a(16),a(14)]; % for some reason adding outline colors look completely diff than desired color. 
set(boxesHC, 'Color','k','LineWidth',3,'LineStyle','-')
boxesHD = [a(5),a(7),a(9),a(11),a(13),a(15),a(17),a(19),a(21),a(23),a(25),a(27)]; 
set(boxesHD, 'Color','k','LineWidth',3,'LineStyle','-')
%ylim([11 50]) %<--- TURN ON WHEN GENERATING MR SENSOR PLOTS


% labels
set(gca,'xticklabel',{'Pre','Post'},'XColor','k','YColor','k')
ylabel('R_2 (%)')

% calculate significance 
am_hc = data.(whichParameter)(data.uniqueID==1); 
am_hd = data.(whichParameter)(data.uniqueID==2); 
pm_hc = data.(whichParameter)(data.uniqueID==3); 
pm_hd = data.(whichParameter)(data.uniqueID==4); 
% pVals_am = permutationTest(am_hc, am_hd, 1e5-1,'exact',1); % Don't trust Matlab's permutation test. Need to run test in R
% pVals_pm = permutationTest(pm_hc, pm_hd, 1e5-1,'exact',1); % exact permutation test with 1e5-1 reps
[~, pVals_am] = ttest2(am_hc, am_hd, 'VarType','unequal');
[~, pVals_pm] = ttest2(pm_hc, pm_hd, 'VarType','unequal');

% add significance starts to plot
model_series = [quantile(data.(whichParameter)(data.uniqueID==1), 0.95), quantile(data.(whichParameter)(data.uniqueID==2), 0.95); quantile(data.(whichParameter)(data.uniqueID==3), 0.95), quantile(data.(whichParameter)(data.uniqueID==4), 0.95)]; 
sigstar2(model_series, get(gca,'XTickLabel'), [pVals_am, pVals_pm],0,'boxplot',h)

% legend and plot details
hLegend = legend(findall(gca,'Tag','Box'), {'HC','HD'},'Location','Best'); legend boxoff

