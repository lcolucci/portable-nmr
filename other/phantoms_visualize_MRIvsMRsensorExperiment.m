%% Phantom Result Visualization (figs. S8 and S9): MRI vs MR

%% Plot T2 Decays
%  Data Prep
clear; close all
doSave = 0; 
savePath = '/Users/linacolucci/Documents/Thesis/Thesis/Figures/MRIvsMR_Phantoms/T2decays'; 
load('/Volumes/CimaLabLargeFiles/Hydration/MRI_Study/ProcessedData/phantoms/phantoms_T2decays.mat')
data = master(master.Trial=='Average',:); % just average decays
data = data(data.nDummies=='0dummy' | data.nDummies=='shorterTE',:); % just 0dummy or shorterTE scans
%  Plot
nMeas = height(data); 
for iMeas = 1:nMeas
    figure; 
    plot(data.Time{iMeas,1}, data.T2Decay{iMeas,1},'LineWidth',2)
    set(gcf, 'color', 'w'); set(gca, 'fontsize',20,'linewidth',2); box off
    xlabel('Time (ms)')
    ylabel('Amplitude (a.u.)')
    title(sprintf('T2 Decay: %s %s %s %s', data.Sample(iMeas), data.Sensor(iMeas), data.nDummies(iMeas), data.Trial(iMeas))); 
    autoYlim = get(gca,'YLim'); 
    ylim([-0.5 autoYlim(2)])
    xlim([0 max(data.Time{iMeas,1})])
    
    %% Save
    if doSave == 1
        if ~(7==exist(savePath)) % exist = 7 when it's a folder
            mkdir(savePath)
        end

        saveas(gcf, fullfile(savePath, sprintf('T2decay_%s_%s.eps',data.Sample(iMeas), data.Sensor(iMeas))), 'epsc')
    end
    
end
clear data
    

%% Plot multi-exp results (MR Sensor and MRI overlaid) 

%% --- 1exp ---
clear; close all
doSave = 1; 
savePath = '/Users/linacolucci/Documents/Thesis/Thesis/Figures/MRIvsMR_Phantoms/MultiExp/Monoexponential'; 
load('/Volumes/CimaLabLargeFiles/Hydration/MRI_Study/ProcessedData/phantoms/phantoms_mri_pixelbypixel.mat')
load('/Volumes/CimaLabLargeFiles/Hydration/MRI_Study/ProcessedData/phantoms/phantoms_multiExp_1exp.mat')
data = result_1exp; 
pixeldata = pixelbypixel;
pixeldata.Sample(strcmp(pixeldata.Sample,'oil')) = {'oilPhantom'}; 
resultsTable = table(); 

% Plot
keynoteOrange = [255 147 0]./255;
keynoteBlue = [0 162 255]./255; 
grey = [150 150 150]./255; 
lineColors = {'k',keynoteOrange}; 
lineStyles = {'-','--'}; 
fnSamples = unique(pixeldata.Sample); 
nSamples = length(fnSamples); 
for iSample = 1:nSamples
    figure; 
    
    % plot pixel-by-pixel
    tempPixelData = pixeldata(strcmp(pixeldata.Sample,fnSamples{iSample}),:); 
    histogram(tempPixelData.Data{1,1}.relax1,'Normalization','probability','EdgeColor','none','FaceColor',grey); hold on
    getYLim = ylim; 
    
    tempResults = table(); 
    tempResults.Sample = fnSamples(iSample); 
    tempResults.MonoExp = 1; 
    data1 = tempPixelData.Data{1,1}.relax1; 
    data1(data1==999)=[]; 
    tempResults.MRI_MeanA = mean(data1); 
    tempResults.MRI_StdevA = std(data1); 
    tempResults.MRI_MeanB = nan; 
    tempResults.MRI_StdevB = nan;
    
    resultsTable = vertcat(resultsTable, tempResults);
    clear tempResults data1 data2
    
    % plot mri ROI 1-exp
    tempData = data(data.Sample== fnSamples{iSample},:); 
    tempData = sortrows(tempData,'Sensor','descend');
    nRows = height(tempData); 
    for j=1:nRows
        plot([tempData.relax1(j) tempData.relax1(j)], [0 getYLim(2)],'LineWidth',5,'Color',lineColors{j},'LineStyle',lineStyles{j})
    end
    set(gcf, 'color', 'w'); set(gca, 'fontsize',23,'linewidth',2); box off
    xlabel('T2 Relaxation Time (ms)')
    %legend(['MRI pixel-by-pixel';cellstr(tempData.Sensor)]); legend boxoff
    title(sprintf('%s 1-exp Results',fnSamples{iSample}))
    
    %% Save
    if doSave == 1
        if ~(7==exist(savePath)) % exist = 7 when it's a folder
            mkdir(savePath)
        end

        saveas(gcf, fullfile(savePath, sprintf('Monoexponential_%s.eps',fnSamples{iSample})), 'epsc')
    end
    
    % clear temp variables
    clear tempData tempPixelData getYLim nRows 
end

        
    

%% --- 2exp ---
clearvars -except resultsTable; close all
doSave = 1; 
savePath = '/Users/linacolucci/Documents/Thesis/Thesis/Figures/MRIvsMR_Phantoms/MultiExp/Biexponential'; 
load('/Volumes/CimaLabLargeFiles/Hydration/MRI_Study/ProcessedData/phantoms/phantoms_mri_pixelbypixel.mat')
load('/Volumes/CimaLabLargeFiles/Hydration/MRI_Study/ProcessedData/phantoms/phantoms_multiExp_2exp.mat')
data = result_2exp; 
pixeldata = pixelbypixel;
pixeldata.Sample(strcmp(pixeldata.Sample,'oil')) = {'oilPhantom'}; 
pixeldata = pixeldata(strcmp(pixeldata.Fit,'2exp'),:); 

% Plot
keynoteOrange = [255 147 0]./255;
keynoteBlue = [0 162 255]./255; 
grey = [150 150 150]./255; 
lineColors = {'k',keynoteOrange}; 
lineStyles = {'-','--'}; 
fnSamples = unique(pixeldata.Sample); 
nSamples = length(fnSamples); 
for iSample = 1:nSamples
    figure; 
    
    % plot pixel-by-pixel
    tempPixelData = pixeldata(strcmp(pixeldata.Sample,fnSamples{iSample}),:); 
    tempResults = vertcat(tempPixelData.Data{1,1}.relax1,tempPixelData.Data{1,1}.relax2); 
    tempResults(tempResults==999) = []; 
    histogram(tempResults,100,'Normalization','probability','EdgeColor','none','FaceColor',grey); hold on
    getYLim = ylim; 
    
    tempResults = table(); 
    tempResults.Sample = fnSamples(iSample); 
    tempResults.MonoExp = 0; 
    data1 = tempPixelData.Data{1,1}.relax1; 
    data1(data1==999)=[]; 
    tempResults.MRI_MeanA = mean(data1); 
    tempResults.MRI_StdevA = std(data1); 
    data2 = tempPixelData.Data{1,1}.relax2; 
    data2(data2==999)=[]; 
    tempResults.MRI_MeanB = mean(data2); 
    tempResults.MRI_StdevB = std(data2);
    
    resultsTable = vertcat(resultsTable, tempResults);
    clear tempResults data1 data2
    
    % plot mri ROI 2-exp
    tempData = data(data.Sample== fnSamples{iSample},:); 
    tempData = sortrows(tempData,'Sensor','descend');
    nRows = height(tempData); 
    for j=1:nRows
        plot([tempData.relax1(j) tempData.relax1(j)], [0 getYLim(2)],'LineWidth',5,'Color',lineColors{j},'LineStyle',lineStyles{j})
        plot([tempData.relax2(j) tempData.relax2(j)], [0 getYLim(2)],'LineWidth',5,'Color',lineColors{j},'LineStyle',lineStyles{j})
    end
    set(gcf, 'color', 'w'); set(gca, 'fontsize',23,'linewidth',2); box off
    xlabel('T2 Relaxation Time (ms)')
    %legend(['MRI pixel-by-pixel';cellstr(tempData.Sensor(1));cellstr(tempData.Sensor(1));cellstr(tempData.Sensor(2));cellstr(tempData.Sensor(2))]); legend boxoff
    title(sprintf('%s 2-exp Results',fnSamples{iSample}))
    
    %% Save
    if doSave == 1
        if ~(7==exist(savePath)) % exist = 7 when it's a folder
            mkdir(savePath)
        end

        saveas(gcf, fullfile(savePath, sprintf('Biexponential_%s.eps',fnSamples{iSample})), 'epsc')
    end
    
    % clear temp variables
    clear tempData tempPixelData getYLim nRows tempResults
end


%% Phantom Stability 
%% All of the PM measurements just have 999 in their fittings. 
clear; close all
doSave = 1; 
savePath = '/Users/linacolucci/Documents/Thesis/Thesis/Figures/MRIvsMR_Phantoms/Stability'; 
load('/Volumes/CimaLabLargeFiles/Hydration/MRI_Study/ProcessedData/phantoms/phantoms_mri_pixelbypixel.mat')
data =pixelbypixel;

whichPhantoms = {'agarShort','agarLong','legSensorPhantom','oil'}; 
nPhantoms = 1:length(whichPhantoms); 
for iPhantom = 1:3 %nPhantoms
    tempData = data(strcmp(data.Sample,whichPhantoms{iPhantom}),:); 
    fnTimePts = unique(tempData.TimePt); 
    nTimePts = length(fnTimePts); 
    timePtResults = table(); 
    for jTimePt = 1:nTimePts
        tempData2 = tempData(strcmp(tempData.TimePt, fnTimePts{jTimePt}),:); 
        nRows = height(tempData2); 
        tempResults = []; 
        for k = 1:nRows
            tempResults = vertcat(tempResults, tempData2.Data{k,1}.relax1); 
        end
        tempResults(tempResults==999)=[]; 
        timePtResults.timePt(jTimePt) = fnTimePts(jTimePt); 
        timePtResults.data(jTimePt) = {tempResults};
        timePtResults.group(jTimePt) = {jTimePt.*ones(size(tempResults))}; 
    end
    
    % Plot
    figure
    xData = []; 
    yData = []; 
    for kTimePt = 1:nTimePts
        xData = [xData; timePtResults.data{kTimePt}];
        yData = [yData; timePtResults.group{kTimePt}];
    end
    boxplot(xData, yData)
    %boxplot([timePtResults.data{1}; timePtResults.data{2}; timePtResults.data{3};timePtResults.data{4}], [timePtResults.group{1}; timePtResults.group{2}; timePtResults.group{3};timePtResults.group{4}])
    set(gca,'XTickLabel',timePtResults.timePt)
    set(gcf, 'color', 'w'); set(gca, 'fontsize',23,'linewidth',2); box off
    ylabel('T2 Relaxation Time (ms)')
    title(sprintf('MRI pixel-by-pixel: %s',whichPhantoms{iPhantom}))
    
    %% Save
    if doSave == 1
        if ~(7==exist(savePath)) % exist = 7 when it's a folder
            mkdir(savePath)
        end

        saveas(gcf, fullfile(savePath, sprintf('Stability_%s.eps',whichPhantoms{iPhantom})), 'epsc')
    end
    
    
end


