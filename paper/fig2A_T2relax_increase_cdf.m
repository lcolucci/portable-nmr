%% Generate the relaxation times are elevated cdf plot (fig 2A)
%  Code based on script: visualize_relax_increase.m
%
%  PREVIOUS SCRIPTS
%       data comes from script: process_pixelbypixel_averages.m (which has
%                   code for plotting individual subjects and HC vs HD averages)
%                   this script plots 95% CI instead of Stdev
%
%  OUTPUTS
%       title: MRI pixelwise 2-exp: \tau_L of the Whole Leg'



clear 

%% ----- USER INPUTS ----
posVector = [560   528   1*560   1*420]; % size of plot (same size as fig. 7 subplots)
fontSize = 26; 
% Save
doSave = 1; % 1 if you want to save figure
savePath = '/Users/linacolucci/Documents/Thesis/Paper/figures/MRI_pixelwise'; % where to save plot? 
% Path
pathToMriFolder = '/Volumes/CimaLabLargeFiles/Hydration'; % path to CimaLabLargeFiles folder in the drive
%% ----------------------

%% Load Data
load(fullfile(pathToMriFolder, 'MRI_Study/ProcessedData/mri/pixel_by_pixel/mri_pixelByPixel_2exp_summary_all.mat'))

%% Define Colors
keynoteBlue = [0 162 255]./255;
keynoteOrange = [255 147 0]./255;
grey = [138 138 138]./255;
ishaRed = [142 51 49]./255; % 157 43 44  
ishaBlue = [32 25 77]./255; % 45 38 90
ishaGreen = [47 85 59]./255;
ishaTeal = [98 144 128]./255; % 85 165 139
scienceRed = [197, 22, 29]./255;

%% Calculate HC mean/std/95%CI vs HD mean/std/95%CI
z = 1.96; % z-stat for 95% CI

% All
data = aggregate(aggregate.ROI=='leg_whole' & aggregate.Parameter=='relax2' & aggregate.AMPM=='AM',:); 
nRows = height(data); 
cdfAll = []; 
for i=1:nRows
    cdfAll = vertcat(cdfAll, data.cdf{i,1}); 
end

cdfAvg = mean(cdfAll,1); 
cdfStd = std(cdfAll,1); 
cdfCI = z .* cdfStd ./ sqrt(height(data)); % CI =  z * stdev / sqrt(n)
cdfBins = data.Bins{1,1}(2:end); 

%% Plotting
% Prep data for plotting
x = [cdfBins, fliplr(cdfBins)]; 
y = [cdfAvg - cdfCI, fliplr(cdfAvg + cdfCI)]; 

% Plot
% population average (including the individually-plotted subjects)
figure; hold on
fill(x, y, grey, 'FaceAlpha',0.4,'EdgeColor','none')
plot(cdfBins, cdfAvg, 'Color',grey, 'LineWidth',3)
% individual subjects
plot(cdfBins, data(data.Subject=='HDMRI01',:).cdf{1,1},'LineWidth',3, 'Color',ishaRed)
plot(cdfBins, data(data.Subject=='HDMRI01b',:).cdf{1,1},'--','LineWidth',3, 'Color',ishaBlue)
plot(cdfBins, data(data.Subject=='HDMRI02b',:).cdf{1,1},'-.','LineWidth',3, 'Color',ishaGreen)


ylabel('Cumulative Probability')
xlabel('\tau_L (ms)')
set(gcf, 'color', 'w','Position',posVector); 
set(gca, 'fontsize',fontSize,'linewidth',3,'XColor','k','YColor','k'); box off
ylim([0 1])
legend('95% CI', 'Mean', 'HD 1','HD 1b','HD 2b', 'Location', 'Best')
legend boxoff
xlim([0 600])



% Save 
if doSave==1
    % make folder if it doesn't exist
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end
    saveas(gcf, fullfile(savePath, sprintf('fig3A_mri_pixelwise_relaxation.svg')))
    %saveas(gcf, fullfile(savePath, sprintf('fig2B_mri_pixelwise_cdf_post.eps')),'epsc')
end
