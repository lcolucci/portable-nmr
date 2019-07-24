%% ---- Fig 2D-F - MRI Pixelwise CDF and Integral of CDF Difference ----
%
%  Code comes from script: visualize_relativeAmp_elevation_etc.m 
%                         (with cosmetic, organizational changes, and made errors 95% CI instead of stdev)
%  
%  PREVIOUS SCRIPTS
%       data comes from script: process_pixelbypixel_averages.m
%
%  OUTPUTS
%       General description for all plots: MRI Pixel-wise 2-exp: Relative Amplitude 2 (R_L) of Muscular Tissue
%       Subplot A - Pre: CDF of HC and HD (mean +/- 95% CI)
%       Subplot B - Post: CDF of HC and HD (mean +/- 95% CI)
%       Subplot C - Change: Integral of AM-PM CDF of HC and HD (mean +/- 95% CI)
%       Supplemental (2x) - Like subplots A and B but showing each individual
%                           subject rather than population averages
%
%  Lina A. Colucci, 2018

clear; close all

%% ----- USER INPUTS ----
posVector = [560   528   1*560   1*420]; % size of plot (same size as fig. 7 subplots)
fontSize = 26; 
% Save
doSave = 1; % 1 if you want to save figure
savePath = '/Users/linacolucci/Documents/Thesis/Paper/figures/MRI_pixelwise'; % where to save plot? 
% Path
pathToMriFolder = '/Volumes/CimaLabLargeFiles/Hydration'; % path to CimaLabLargeFiles folder in the drive
%% ----------------------

%% Load data
load(fullfile(pathToMriFolder, 'MRI_Study/ProcessedData/mri/pixel_by_pixel/mri_pixelByPixel_2exp_summary_all.mat'))

%% Set colors
keynoteBlue = [0 162 255]./255;
keynoteOrange = [255 147 0]./255;
grey = [138 138 138]./255;
ishaRed = [142 51 49]./255; % 157 43 44  
ishaBlue = [32 25 77]./255; % 45 38 90
ishaGreen = [47 85 59]./255;
ishaTeal = [98 144 128]./255; % 85 165 139
scienceRed = [197, 22, 29]./255;

%% Prep Data
z = 1.96; % z-stat for 95% CI
%n = 7; % sample size

% average HC AM
data = aggregate(aggregate.ROI=='muscle_all' & aggregate.Parameter=='relative_amp2' & aggregate.AMPM=='AM' & aggregate.HD==0,:); 
nRows = height(data); 
cdfAll = []; 
for i=1:nRows
    cdfAll = vertcat(cdfAll, data.cdf{i,1}); 
end
cdfAvg_HC_AM = mean(cdfAll,1); 
cdfStd_HC_AM = std(cdfAll,1); 
cdfCI_HC_AM = z .* cdfStd_HC_AM ./ sqrt(height(data)); % CI =  z * stdev / sqrt(n)
cdfBins = data.Bins{1,1}(2:end); 

clear data cdfAll

% average HC PM
data = aggregate(aggregate.ROI=='muscle_all' & aggregate.Parameter=='relative_amp2' & aggregate.AMPM=='PM' & aggregate.HD==0,:); 
nRows = height(data); 
cdfAll = []; 
for i=1:nRows
    cdfAll = vertcat(cdfAll, data.cdf{i,1}); 
end
cdfAvg_HC_PM = mean(cdfAll,1); 
cdfStd_HC_PM = std(cdfAll,1); 
cdfCI_HC_PM = z .* cdfStd_HC_PM ./ sqrt(height(data)); % CI =  z * stdev / sqrt(n)
clear data cdfAll

% average HD AM
data = aggregate(aggregate.ROI=='muscle_all' & aggregate.Parameter=='relative_amp2' & aggregate.AMPM=='AM' & aggregate.HD==1,:); 
nRows = height(data); 
cdfAll = []; 
for i=1:nRows
    cdfAll = vertcat(cdfAll, data.cdf{i,1}); 
end
cdfAvg_HD_AM = mean(cdfAll,1); 
cdfStd_HD_AM = std(cdfAll,1); 
cdfCI_HD_AM = z .* cdfStd_HD_AM ./ sqrt(height(data)); % CI =  z * stdev / sqrt(n)
clear data cdfAll

% average HD PM
data = aggregate(aggregate.ROI=='muscle_all' & aggregate.Parameter=='relative_amp2' & aggregate.AMPM=='PM' & aggregate.HD==1,:); 
nRows = height(data); 
cdfAll = []; 
for i=1:nRows
    cdfAll = vertcat(cdfAll, data.cdf{i,1}); 
end
cdfAvg_HD_PM = mean(cdfAll,1); 
cdfStd_HD_PM = std(cdfAll,1); 
cdfCI_HD_PM = z .* cdfStd_HD_PM ./ sqrt(height(data)); % CI =  z * stdev / sqrt(n)
clear data cdfAll



%% --- Subplot A (AM) and B (PM): Avg HC vs HD CDFs ---
% Prep data for plotting
x = [cdfBins, fliplr(cdfBins)]; 
y_HC_AM = [cdfAvg_HC_AM - cdfCI_HC_AM, fliplr(cdfAvg_HC_AM + cdfCI_HC_AM)]; 
y_HD_AM = [cdfAvg_HD_AM - cdfCI_HD_AM, fliplr(cdfAvg_HD_AM + cdfCI_HD_AM)]; 
y_HC_PM = [cdfAvg_HC_PM - cdfCI_HC_PM, fliplr(cdfAvg_HC_PM + cdfCI_HC_PM)]; 
y_HD_PM = [cdfAvg_HD_PM - cdfCI_HD_PM, fliplr(cdfAvg_HD_PM + cdfCI_HD_PM)]; 

% ------- AM --------
figure; hold on
fill(x, y_HC_AM, grey, 'FaceAlpha',0.4,'EdgeColor','none')
plot(cdfBins, cdfAvg_HC_AM, 'Color',grey, 'LineWidth',3)

fill(x, y_HD_AM, ishaRed, 'FaceAlpha',0.4,'EdgeColor','none')
plot(cdfBins, cdfAvg_HD_AM, 'Color',ishaRed, 'LineWidth',3)

legend('HC 95% CI','HC Mean','HD 95% CI','HD Mean','Location','Best'); legend boxoff
ylabel('Cumulative Probability')
xlabel('R_L (%)')
title('Pre')
set(gcf, 'color', 'w','Position',posVector); set(gca, 'fontsize',fontSize,'linewidth',2,'XColor','k','YColor','k'); box off
ylim([0 1])
set(gcf,'Position',posVector)

% Save 
if doSave==1
    % make folder if it doesn't exist
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end
    saveas(gcf, fullfile(savePath, sprintf('fig2A_mri_pixelwise_cdf_pre.svg')))
    % saveas(gcf, fullfile(savePath, sprintf('fig2A_mri_pixelwise_cdf_pre.eps')),'epsc')
end

% ------- PM -------
figure; hold on
fill(x, y_HC_PM, grey, 'FaceAlpha',0.4,'EdgeColor','none')
plot(cdfBins, cdfAvg_HC_PM, 'Color',grey, 'LineWidth',3)

fill(x, y_HD_PM, ishaRed, 'FaceAlpha',0.4,'EdgeColor','none')
plot(cdfBins, cdfAvg_HD_PM, 'Color',ishaRed, 'LineWidth',3)

legend('HC 95% CI','HC Mean','HD 95% CI','HD Mean','Location','Best'); legend boxoff
ylabel('Cumulative Probability')
xlabel('R_L (%)')
title('Post')
set(gcf, 'color', 'w','Position',posVector); set(gca, 'fontsize',fontSize,'linewidth',2,'XColor','k','YColor','k'); box off
ylim([0 1])

% Save 
if doSave==1
    % make folder if it doesn't exist
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end
    saveas(gcf, fullfile(savePath, sprintf('fig2B_mri_pixelwise_cdf_post.svg')))
    %saveas(gcf, fullfile(savePath, sprintf('fig2B_mri_pixelwise_cdf_post.eps')),'epsc')
end

%% ---- Subplot C: Cumulative Change (Integral of CDF Diff) ------
% average HC 
data = aggregate(aggregate.ROI=='muscle_all' & aggregate.Parameter=='relative_amp2' & aggregate.AMPM=='PM' & aggregate.HD==0,:); 
data(data.Subject=='HC02',:) = []; 
data(data.Subject=='HC06',:) = []; 
nRows = height(data); 
cdfAll = []; 
for i=1:nRows
    cdfAll = vertcat(cdfAll, data.integralDiff{i,1}); 
end
integralAvg_HC = mean(cdfAll,1); 
integralStd_HC = std(cdfAll,1); 
integralCI_HC = z .* integralStd_HC ./ sqrt(height(data)); % CI =  z * stdev / sqrt(n)
cdfBins = data.Bins{1,1}(2:end); 
clear data cdfAll

% average HD
data = aggregate(aggregate.ROI=='muscle_all' & aggregate.Parameter=='relative_amp2' & aggregate.AMPM=='PM' & aggregate.HD==1,:); 
nRows = height(data); 
cdfAll = []; 
for i=1:nRows
    cdfAll = vertcat(cdfAll, data.integralDiff{i,1}); 
end
integralAvg_HD = mean(cdfAll,1); 
integralStd_HD = std(cdfAll,1); 
integralCI_HD = z .* integralStd_HD ./ sqrt(height(data)); % CI =  z * stdev / sqrt(n)
cdfBins = data.Bins{1,1}(2:end); 
clear data cdfAll

% Prep data
x = [cdfBins, fliplr(cdfBins)]; 
y_HC = [integralAvg_HC - integralCI_HC, fliplr(integralAvg_HC + integralCI_HC)]; 
y_HD = [integralAvg_HD - integralCI_HD, fliplr(integralAvg_HD + integralCI_HD)]; %multiply by -1 so that HDs "decrease"

% Plot
figure; hold on
fill(x, -1.*y_HC, grey, 'FaceAlpha',0.4,'EdgeColor','none')
plot(cdfBins, -1.*integralAvg_HC, 'Color',grey, 'LineWidth',4)

fill(x, -1.*y_HD, ishaRed, 'FaceAlpha',0.4,'EdgeColor','none')
plot(cdfBins, -1.*integralAvg_HD, 'Color',ishaRed, 'LineWidth',4)

% Add individual subjects - HC02 and HC06
data = aggregate(aggregate.ROI=='muscle_all' & aggregate.Parameter=='relative_amp2' & aggregate.AMPM=='PM' & aggregate.Subject=='HC02',:); 
plot(cdfBins, -1.*data.integralDiff{1,1},'--','Color',ishaGreen, 'LineWidth',4)
data = aggregate(aggregate.ROI=='muscle_all' & aggregate.Parameter=='relative_amp2' & aggregate.AMPM=='PM' & aggregate.Subject=='HC06',:); 
plot(cdfBins, -1.*data.integralDiff{1,1},'-.','Color',ishaBlue, 'LineWidth',4)

legend('HC 95% CI','HC Mean','HD 95% CI','HD Mean','HC 2', 'HC 6','Location','Best'); legend boxoff
ylabel('Cumulative Probability')
xlabel('R_L (%)')
title('AM-to-PM Change')
set(gcf, 'color', 'w','Position',posVector); set(gca, 'fontsize',fontSize,'linewidth',2,'XColor','k','YColor','k'); box off

% Save 
if doSave==1
    % make folder if it doesn't exist
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end
    saveas(gcf, fullfile(savePath, sprintf('fig2C_mri_pixelwise_cdf_change.svg')))
    %saveas(gcf, fullfile(savePath, sprintf('fig2C_mri_pixelwise_cdf_change.eps')),'epsc')
end


%% ----- Supplemental: CDF of Individual Subjects ------
% ----- AM -----
data = aggregate(aggregate.ROI=='muscle_all' & aggregate.Parameter=='relative_amp2' & aggregate.AMPM=='AM' & aggregate.HD==0,:); 
nRows = height(data); 
figure
for i=1:nRows
    plot(data.Bins{i,1}(2:end), data.cdf{i,1},'LineWidth',3); hold on
end
data = aggregate(aggregate.ROI=='muscle_all' & aggregate.Parameter=='relative_amp2' & aggregate.AMPM=='AM' & aggregate.HD==1,:); 
nRows = height(data);
for i=1:nRows
    plot(data.Bins{i,1}(2:end), data.cdf{i,1},'--','LineWidth',4); hold on
end
legend(cellstr(unique(aggregate.Subject)),'Location','Best'); legend boxoff
ylabel('Cumulative Probability')
xlabel('R_L (%)')
title('Pre MRI Pixelwise 2-exp: R_L of Muscular Tissue')
set(gcf, 'color', 'w','Position',2.*posVector); set(gca, 'fontsize',fontSize,'linewidth',3,'XColor','k','YColor','k'); box off
ylim([0 1])

% Save 
if doSave==1
    % make folder if it doesn't exist
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end
    saveas(gcf, fullfile(savePath, sprintf('supplemental_mri_pixelwise_cdf_pre_individuals.eps')),'epsc')
end
clear data 

% ----- PM -----
data = aggregate(aggregate.ROI=='muscle_all' & aggregate.Parameter=='relative_amp2' & aggregate.AMPM=='PM' & aggregate.HD==0,:); 
nRows = height(data); 
figure
for i=1:nRows
    plot(data.Bins{i,1}(2:end), data.cdf{i,1},'LineWidth',3); hold on
end
data = aggregate(aggregate.ROI=='muscle_all' & aggregate.Parameter=='relative_amp2' & aggregate.AMPM=='PM' & aggregate.HD==1,:); 
nRows = height(data);
for i=1:nRows
    plot(data.Bins{i,1}(2:end), data.cdf{i,1},'--','LineWidth',4); hold on
end
legend(cellstr(unique(aggregate.Subject)),'Location','Best'); legend boxoff
ylabel('Cumulative Probability')
xlabel('R_L (%)')
title('Post MRI Pixelwise 2-exp: R_L of Muscular Tissue')
set(gcf, 'color', 'w','Position',2.*posVector); set(gca, 'fontsize',fontSize,'linewidth',3,'XColor','k','YColor','k'); box off
ylim([0 1])

% Save 
if doSave==1
    % make folder if it doesn't exist
    if ~(7==exist(savePath)) % exist = 7 when it's a folder
        mkdir(savePath)
    end
    saveas(gcf, fullfile(savePath, sprintf('supplemental_mri_pixelwise_cdf_post_individuals.eps')),'epsc')
end