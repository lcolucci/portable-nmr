%% Generate figures for paper - all of the figures with MRI Relative Amplitude

% Set colors
keynoteBlue = [0 162 255]./255;
keynoteOrange = [255 147 0]./255;
grey = [138 138 138]./255;
ishaRed = [142 51 49]./255; % 157 43 44  
ishaBlue = [32 25 77]./255; % 45 38 90
ishaGreen = [47 85 59]./255;
ishaTeal = [98 144 128]./255; % 85 165 139
scienceRed = [197, 22, 29]./255;

% Load data
load('/Volumes/cimalablargefiles/Hydration/MRI_Study/ProcessedData/mri/pixel_by_pixel/mri_pixelByPixel_2exp_summary_all.mat')

%% average HC AM
data = aggregate(aggregate.ROI=='muscle_all' & aggregate.Parameter=='relative_amp2' & aggregate.AMPM=='AM' & aggregate.HD==0,:); 
nRows = height(data); 
cdfAll = []; 
for i=1:nRows
    cdfAll = vertcat(cdfAll, data.cdf{i,1}); 
end
cdfAvg_HC_AM = mean(cdfAll,1); 
cdfStd_HC_AM = std(cdfAll,1); 
cdfBins = data.Bins{1,1}(2:end); 

clear data cdfAll

%% average HC PM
data = aggregate(aggregate.ROI=='muscle_all' & aggregate.Parameter=='relative_amp2' & aggregate.AMPM=='PM' & aggregate.HD==0,:); 
nRows = height(data); 
cdfAll = []; 
for i=1:nRows
    cdfAll = vertcat(cdfAll, data.cdf{i,1}); 
end
cdfAvg_HC_PM = mean(cdfAll,1); 
cdfStd_HC_PM = std(cdfAll,1); 
clear data cdfAll

%% average HD AM
data = aggregate(aggregate.ROI=='muscle_all' & aggregate.Parameter=='relative_amp2' & aggregate.AMPM=='AM' & aggregate.HD==1,:); 
nRows = height(data); 
cdfAll = []; 
for i=1:nRows
    cdfAll = vertcat(cdfAll, data.cdf{i,1}); 
end
cdfAvg_HD_AM = mean(cdfAll,1); 
cdfStd_HD_AM = std(cdfAll,1); 
clear data cdfAll

%% average HD PM
data = aggregate(aggregate.ROI=='muscle_all' & aggregate.Parameter=='relative_amp2' & aggregate.AMPM=='PM' & aggregate.HD==1,:); 
nRows = height(data); 
cdfAll = []; 
for i=1:nRows
    cdfAll = vertcat(cdfAll, data.cdf{i,1}); 
end
cdfAvg_HD_PM = mean(cdfAll,1); 
cdfStd_HD_PM = std(cdfAll,1); 
clear data cdfAll

%% Plots

%% ============================================ %%
%% Individual plots
% AM
data = aggregate(aggregate.ROI=='muscle_all' & aggregate.Parameter=='relative_amp2' & aggregate.AMPM=='AM' & aggregate.HD==0,:); 
nRows = height(data); 
figure
for i=1:nRows
    plot(data.Bins{i,1}(2:end), data.cdf{i,1},'LineWidth',3); hold on
end
data = aggregate(aggregate.ROI=='muscle_all' & aggregate.Parameter=='relative_amp2' & aggregate.AMPM=='AM' & aggregate.HD==1,:); 
nRows = height(data);
for i=1:nRows
    plot(data.Bins{i,1}(2:end), data.cdf{i,1},'--','LineWidth',3); hold on
end
legend(cellstr(unique(aggregate.Subject)),'Location','Best'); legend boxoff
ylabel('Cumulative Probability')
xlabel('Relative Amplitude_{Long} (%)')
title('Relative Amplitude_{Long} of Muscular Tissue')
set(gcf, 'color', 'w'); set(gca, 'fontsize',25); box off
ylim([0 1])
set(gcf,'Position',[ -1415  230  740   533])

%% ============================================ %%
%% --- Avg AM HC vs HD---
% Prep data for plotting
x = [cdfBins, fliplr(cdfBins)]; 
y_HC_AM = [cdfAvg_HC_AM - cdfStd_HC_AM, fliplr(cdfAvg_HC_AM + cdfStd_HC_AM)]; 
y_HD_AM = [cdfAvg_HD_AM - cdfStd_HD_AM, fliplr(cdfAvg_HD_AM + cdfStd_HD_AM)]; 
y_HC_PM = [cdfAvg_HC_PM - cdfStd_HC_PM, fliplr(cdfAvg_HC_PM + cdfStd_HC_PM)]; 
y_HD_PM = [cdfAvg_HD_PM - cdfStd_HD_PM, fliplr(cdfAvg_HD_PM + cdfStd_HD_PM)]; 

% AM
figure; hold on
fill(x, y_HC_AM, grey, 'FaceAlpha',0.4,'EdgeColor','none')
plot(cdfBins, cdfAvg_HC_AM, 'Color',grey, 'LineWidth',3)

fill(x, y_HD_AM, ishaRed, 'FaceAlpha',0.4,'EdgeColor','none')
plot(cdfBins, cdfAvg_HD_AM, 'Color',ishaRed, 'LineWidth',3)

legend('Healthy Controls Stdev','Healthy Controls Mean','Hemodialysis Subjects Stdev','Hemodialysis Subjects Mean','Location','Best'); legend boxoff
ylabel('Cumulative Probability')
xlabel('Relative Amplitude_{Long} (%)')
title('Pre: Relative Amplitude_{Long} of Muscular Tissue')
set(gcf, 'color', 'w'); set(gca, 'fontsize',25); box off
ylim([0 1])
set(gcf,'Position',[ -1415  230  740   533])

% PM
figure; hold on
fill(x, y_HC_PM, grey, 'FaceAlpha',0.4,'EdgeColor','none')
plot(cdfBins, cdfAvg_HC_PM, 'Color',grey, 'LineWidth',3)

fill(x, y_HD_PM, ishaRed, 'FaceAlpha',0.4,'EdgeColor','none')
plot(cdfBins, cdfAvg_HD_PM, 'Color',ishaRed, 'LineWidth',3)

legend('Healthy Controls Stdev','Healthy Controls Mean','Hemodialysis Subjects Stdev','Hemodialysis Subjects Mean','Location','Best'); legend boxoff
ylabel('Cumulative Probability')
xlabel('Relative Amplitude_{Long} (%)')
title('Post: Relative Amplitude_{Long} of Muscular Tissue')
set(gcf, 'color', 'w'); set(gca, 'fontsize',25); box off
ylim([0 1])
set(gcf,'Position',[ -1415  230  740   533])

%% ============================================ %%
%% Cumulative Change
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
cdfBins = data.Bins{1,1}(2:end); 
clear data cdfAll

% Prep data
x = [cdfBins, fliplr(cdfBins)]; 
y_HC = [integralAvg_HC - integralStd_HC, fliplr(integralAvg_HC + integralStd_HC)]; 
y_HD = [integralAvg_HD - integralStd_HD, fliplr(integralAvg_HD + integralStd_HD)]; 

% Plot
figure; hold on
fill(x, y_HC, grey, 'FaceAlpha',0.4,'EdgeColor','none')
plot(cdfBins, integralAvg_HC, 'Color',grey, 'LineWidth',4)

fill(x, y_HD, ishaRed, 'FaceAlpha',0.4,'EdgeColor','none')
plot(cdfBins, integralAvg_HD, 'Color',ishaRed, 'LineWidth',4)

% Add individual subjects - HC02 and HC06
data = aggregate(aggregate.ROI=='muscle_all' & aggregate.Parameter=='relative_amp2' & aggregate.AMPM=='PM' & aggregate.Subject=='HC02',:); 
plot(cdfBins, data.integralDiff{1,1},'--','Color',ishaGreen, 'LineWidth',4)
data = aggregate(aggregate.ROI=='muscle_all' & aggregate.Parameter=='relative_amp2' & aggregate.AMPM=='PM' & aggregate.Subject=='HC06',:); 
plot(cdfBins, data.integralDiff{1,1},'-.','Color',ishaBlue, 'LineWidth',4)

legend('Healthy Controls Stdev','Healthy Controls Mean','Hemodialysis Subjects Stdev','Hemodialysis Subjects Mean','HC2','HC6','Location','Best'); 
legend boxoff
xlabel('Relative Amplitude_{Long} (%)')
title('Change in Relative Amplitude_{Long} of Muscular Tissue')
set(gcf, 'color', 'w'); set(gca, 'fontsize',25); box off
set(gcf,'Position',[ -1415  230  740   533])