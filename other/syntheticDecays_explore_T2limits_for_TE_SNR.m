%% Exploration: For a given TE/SNR, what's the smallest/largest T2 you can measure? 
%  This script is aimed at answering the question: What's the lowest T2
%  that can be fit for a given TE/SNR? 
%  This will inform which relaxation times we throw away in the pixel-by-pixel 
%  fitting because they are too low. 
%       
%% -------------------------- NOTES ------------------------------------
% noise = (stdev of background)/sqrt(2-pi/2) = (stdev of background)/0.655
% SNR = maxSignal * 0.655/(stdev of background)

clear; close all

%% ----- USER INPUTS -----
% -- Decay Curve Definition 
TE = 8;  % what TE (echo time) spacing to use? (ms)
nPts = 32;  % How many points in T2 decay curve? 
time = TE:TE:TE*nPts;  % What time array to use? TE:TE:TE*nPts 
amp = 3000;  % What amplitude to use for the T2 arrays? (a.u.) Not sure that this matters to the results. 
% -- Skip first point?
skip1stPoint = 0;  %0 or 1. 0 = TE:TE:TE*nPts, 1 = 2*TE:TE:TE*nPts
additionalTitleText = ''; % Can be empty '' or '(skip 1st point)'
% -- T2 and SNR 
T2TeRatios = [.25 0.5 0.75 1 1.25 1.5 1.75 2 3 4 5 10 20 30 40 50 60 70 80 90 100 150 200]; % what T2/TE ratios to analyze? T2TeRatios * TE = T2 vector
targetSNR = [20, 35, 50, 100, 175, 250]; % Which SNRs do you want to analyze? 
% -- Save 
doSave = 0; % 1 to save
savePath = '/Users/linacolucci/Documents/Thesis/Figures/MRI/Pixel_by_Pixel/Criteria/Synthetic_Curves/Maximum_T2'; 
%% ---- END USER INPUTS ---
    
%% T2 Fits
lineStyles = repmat({'-','--',':','-.'},1,10); 
if skip1stPoint == 1
    time(1) = []; % ignore 1st point
end
results = table(); rowID = 0; 
nSNR = length(targetSNR); 
for iSNR=1:nSNR
    T2 = TE.*T2TeRatios; % generate T2 array 
    rowID = rowID+1; % which row of results table to populate? 
    fitresults = table(); figure
    
    for jT2=1:length(T2)
        
        % Generate raw T2 Decay 
        yRaw = amp.*exp(-time./T2(jT2)); 
        
        % Add wgn to the T2 decay to achieve target SNR
        [y, targetStd(jT2), actualStd(jT2), error_targetStd_actualStd(jT2), actualWgnPower(jT2)] = awgn_snr(yRaw, targetSNR(iSNR)); 
        snr2(jT2) = max(y) ./ actualStd(jT2); %calculate snr
        
        % Fit the noisy T2 decay
        fitresults(jT2,:) = t2_fitting_1exp(time, y);

        % Plot the noisy T2 decay 
        plot(time, y,'LineWidth',2,'LineStyle',lineStyles{jT2}); hold on
        set(gcf, 'color', 'w','Position',[64 49 1124 906]); set(gca, 'fontsize',16)
        xlabel('Time (ms)')
        ylabel('Amplitude (a.u.)')
        title(sprintf('Synthetically Generated T2 Decays with SNR ~%d %s', targetSNR(iSNR), additionalTitleText))

        % Clear temp variables
        clear y yRaw 
    end
    % Final plot formatting
    legend(cellstr(num2str(T2', 'T2=%.1fms'))); legend boxoff
    ylim([-200 amp+200]); xlim([0 max(time)+10])
    
    % Add additional paramters to the fitresults table
    fitresults.targetStd(:) = targetStd; 
    fitresults.actualStd(:) = actualStd; 
    fitresults.error_targetStd_actualStd(:) = error_targetStd_actualStd; 
    fitresults.actualWgnPower(:) = actualWgnPower; 

    fitresults.actualT2(:) = T2';
    fitresults.TE(:) = TE; 
    fitresults.t2_div_te = fitresults.actualT2 ./ fitresults.TE; 
    fitresults.SNR2(:) = snr2'; 
    fitresults.errorT2 = (fitresults.actualT2 - fitresults.relax1); 
    fitresults.errorT2_percent = fitresults.errorT2 ./ fitresults.actualT2 .*100; 

    %% Aggregate into massive data summary
    results.TE(rowID) = TE; 
    results.TargetSNR(rowID) = targetSNR(iSNR); 
    results.Skip1stPoint(rowID) = skip1stPoint; 
    results.Fits(rowID) = {fitresults}; 
    
    %% Save
    if doSave == 1
        if ~(7==exist(savePath)) % exist = 7 when it's a folder
                mkdir(savePath)
        end 
        saveas(gcf, fullfile(savePath, sprintf('syntheticDecays_TE%.1fms_SNR%.0f.png',TE, targetSNR(iSNR))) )
        saveas(gcf, fullfile(savePath, sprintf('syntheticDecays_TE%.1fms_SNR%.0f.eps',TE, targetSNR(iSNR))),'epsc' )
    end
end

%% Visualize Results: Both Percent Error (%) and Absolute Error (ms)    
% Loop through both error types
errorPlotYDataLoop = {'errorT2','errorT2_percent'}; % either 'errorT2' or 'errorT2_percent'
errorPlotYNameLoop = {'T_2 measurement error (ms)', 'T_2 measurement error (%)'}; % either 'T_2 measurement error (ms)' or 'T_2 measurement error (%)'

for iErrorType = 1:length(errorPlotYDataLoop)
    errorPlotYData = errorPlotYDataLoop{iErrorType}; 
    errorPlotYName = errorPlotYNameLoop{iErrorType}; 

    %% Visualize Results                                        
    dataStructColumn = 'Fits'; 
    xData = ''; % use default t2_div_te
    [figErrors, ax, textLabels] = plot_T2errors_vs_T2values_multipleSNRs(results, dataStructColumn, xData, errorPlotYData, errorPlotYName); 
    title(sprintf('Accuracy of T2 Measurements for TE=%.1fms based on SNR %s', results.TE(1), additionalTitleText))                                                    
    % Save
    if doSave == 1
        if ~(7==exist(savePath)) % exist = 7 when it's a folder
                mkdir(savePath)
        end 
        if strcmp(errorPlotYData,'errorT2_percent')
            saveas(gcf, fullfile(savePath, sprintf('syntheticDecays_ErrorBasedOnSNR_TE%.1fms_percent.png',TE)) )
            saveas(gcf, fullfile(savePath, sprintf('syntheticDecays_ErrorBasedOnSNR_TE%.1fms_percent.eps',TE)),'epsc' )
        end
        if strcmp(errorPlotYData,'errorT2')
            saveas(gcf, fullfile(savePath, sprintf('syntheticDecays_ErrorBasedOnSNR_TE%.1fms_absolute.png',TE)) )
            saveas(gcf, fullfile(savePath, sprintf('syntheticDecays_ErrorBasedOnSNR_TE%.1fms_absolute.eps',TE)),'epsc' )
        end
    end   

    %% Zoomed in Plot 
    xlim([0 10])
    delete(textLabels); % delete previous text labels
    xTicks = get(gca, 'xtick');
    yTicks = get(gca,'ytick');
    minY = min(yTicks);
    if strcmp(errorPlotYData,'errorT2_percent')
        HorizontalOffset = .05;
        VerticalOffset = 1;
    end
    if strcmp(errorPlotYData,'errorT2')
        HorizontalOffset = .07;
        VerticalOffset = 0.18;
    end
    for xx = 1:length(xTicks)
        h(xx) = text(xTicks(xx)-HorizontalOffset, minY - VerticalOffset, sprintf('%.0fms',TE.*xTicks(xx)), 'Fontsize',16); % Create a text box at every Tick label position
    end
    legend('Location','Best')
    % Save
    if doSave == 1
         if ~(7==exist(savePath)) % exist = 7 when it's a folder
                mkdir(savePath)
        end 
        if strcmp(errorPlotYData,'errorT2_percent')
            saveas(gcf, fullfile(savePath, sprintf('syntheticDecays_ErrorBasedOnSNR_TE%.1fms_percent_zoom.png',TE)) )
            saveas(gcf, fullfile(savePath, sprintf('syntheticDecays_ErrorBasedOnSNR_TE%.1fms_percent_zoom.eps',TE)),'epsc' )
        end
        if strcmp(errorPlotYData,'errorT2')
            saveas(gcf, fullfile(savePath, sprintf('syntheticDecays_ErrorBasedOnSNR_TE%.1fms_absolute_zoom.png',TE)) )
            saveas(gcf, fullfile(savePath, sprintf('syntheticDecays_ErrorBasedOnSNR_TE%.1fms_absolute_zoom.eps',TE)),'epsc' )
        end
    end  
end

%% Additional Code for Alternate Graphs
% line([0 10], [5 5],'Color','k','LineStyle','--','LineWidth',2)
line([0 200], [-10 -10],'Color','k','LineStyle','--','LineWidth',2)

%% Max. Measurable T2 vs SNR
% x = [ 20 35 50 100]; y = [ 51*8, 63*8, 73*8, 800]; <-- manually determined from plot
% cftool, best fit line = 4.789*x + 328.6, r2=0.9924, adjr2=0.9885
figure; scatter(x, y, 100, 'filled','MarkerFaceColor','k')
xline = [0:10:110]; yline = 4.789*xline + 328.6; 
hold on; plot(xline,yline,'k--') 
set(gcf, 'color', 'w'); set(gca, 'fontsize',16)
xlabel('SNR')
ylabel('Maximum Measurable T2 (ms)')
title('Maximum Measurable T2 vs SNR')
annotation('textbox', 'String','Max. Measurable T2 = 4.789*SNR + 328.6', 'FitBoxToText','on','LineStyle','none','FontSize',16)%,'Units','normalized','Position',[0.5 0.5 0.5 0.5])
    
    