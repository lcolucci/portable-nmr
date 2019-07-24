%% Generate the relaxation times are elevated cdf plot

clear 

%% Load Data
load('/Volumes/cimalablargefiles/Hydration/MRI_Study/ProcessedData/mri/pixel_by_pixel/mri_pixelByPixel_2exp_summary_all.mat')
keynoteBlue = [0 162 255]./255;
keynoteOrange = [255 147 0]./255;
grey = [138 138 138]./255;
ishaRed = [142 51 49]./255; % 157 43 44  
ishaBlue = [32 25 77]./255; % 45 38 90
ishaGreen = [47 85 59]./255;
ishaTeal = [98 144 128]./255; % 85 165 139
scienceRed = [197, 22, 29]./255;

%% Calculate HC mean/std vs HD mean/std
% All
data = aggregate(aggregate.ROI=='leg_whole' & aggregate.Parameter=='relax2' & aggregate.AMPM=='AM',:); 
nRows = height(data); 
cdfAll = []; 
for i=1:nRows
    cdfAll = vertcat(cdfAll, data.cdf{i,1}); 
end

cdfAvg = mean(cdfAll,1); 
cdfStd = std(cdfAll,1); 
cdfBins = data.Bins{1,1}(2:end); 

%% Plot figure from paper (aggregate of all subjects and HD1, HD1b, HD2b overlaid)

% Prep data for plotting
x = [cdfBins, fliplr(cdfBins)]; 
y = [cdfAvg - cdfStd, fliplr(cdfAvg + cdfStd)]; 

% Plot
figure; hold on
fill(x, y, grey, 'FaceAlpha',0.4,'EdgeColor','none')
plot(cdfBins, cdfAvg, 'Color',grey, 'LineWidth',3)

plot(cdfBins, data(data.Subject=='HDMRI01',:).cdf{1,1},'LineWidth',3, 'Color',ishaRed)
plot(cdfBins, data(data.Subject=='HDMRI01b',:).cdf{1,1},'--','LineWidth',3, 'Color',ishaBlue)
plot(cdfBins, data(data.Subject=='HDMRI02b',:).cdf{1,1},'-.','LineWidth',3, 'Color',ishaGreen)

ylabel('Cumulative Probability')
xlabel('T_2 Relaxation Time (ms)')
title('T_2 Relaxation Time of the Whole Leg')
set(gcf, 'color', 'w'); set(gca, 'fontsize',21); box off
ylim([0 1])
legend('Standard Deviation', 'Mean', 'HD Subject 1','HD Subject 1b','HD Subject 2b', 'Location', 'Best')
legend boxoff
xlim([0 600])



%% ----- Plotting Separate HC vs HD as well as individual subjects ---------
aggregate2 = aggregate; 
% HC
data = aggregate2(aggregate2.ROI=='leg_whole' & aggregate2.Parameter=='relax2' & aggregate2.AMPM=='AM' & aggregate2.HD==0,:); 
nRows = height(data); 
cdfAll_HC = []; 
for i=1:nRows
    cdfAll_HC = vertcat(cdfAll_HC, data.cdf{i,1}); 
end

cdfAvg_HC = mean(cdfAll_HC,1); 
cdfStd_HC = std(cdfAll_HC,1); 
cdfBins = data.Bins{1,1}(2:end); 

% HD
data = aggregate2(aggregate2.ROI=='leg_whole' & aggregate2.Parameter=='relax2' & aggregate2.AMPM=='AM' & aggregate2.HD==1,:); 
nRows = height(data); 
cdfAll_HD = []; 
for i=1:nRows
    cdfAll_HD = vertcat(cdfAll_HD, data.cdf{i,1}); 
end

cdfAvg_HD = mean(cdfAll_HD,1); 
cdfStd_HD = std(cdfAll_HD,1); 

% prep
y_HC = [cdfAvg_HC - cdfStd_HC, fliplr(cdfAvg_HC + cdfStd_HC)]; 
y_HD = [cdfAvg_HD - cdfStd_HD, fliplr(cdfAvg_HD + cdfStd_HD)]; 

% plot
figure; hold on
fill(x, y_HC,grey,'FaceAlpha',0.4,'EdgeColor','none')
plot(cdfBins, cdfAvg_HC,'Color',grey,'LineWidth',2); hold on; 
fill(x, y_HD,keynoteBlue,'FaceAlpha',0.4,'EdgeColor','none')
plot(cdfBins, cdfAvg_HD,'Color',keynoteBlue,'LineWidth',2); hold on; 

set(gcf, 'color', 'w'); set(gca, 'fontsize',23); box off

% plot indiv
data = aggregate2(aggregate2.ROI=='leg_whole' & aggregate2.Parameter=='relax2' & aggregate2.AMPM=='AM',:); 
plot(cdfBins, data.cdf{8,1},'LineWidth',2)
plot(cdfBins, data.cdf{9,1},'LineWidth',2)
plot(cdfBins, data.cdf{10,1},'LineWidth',2)
plot(cdfBins, data.cdf{11,1},'LineWidth',2)
plot(cdfBins, data.cdf{12,1},'LineWidth',2)
plot(cdfBins, data.cdf{3,1},'LineWidth',2)
plot(cdfBins, data.cdf{7,1},'LineWidth',2)

legend('Avg HC', 'HC Std', 'Avg HD', 'HD Std','HD01','HD01b','HD02','HD02b','HD03','HC02')


