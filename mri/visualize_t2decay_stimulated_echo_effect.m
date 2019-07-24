%% Thesis Figures


%% Show T2 Decay Curve w/ Stimulated Echo (Lower 1st point)
load('/Volumes/CimaLabLargeFiles/Hydration/MRI_Study/ProcessedData/mri/ROIs/ROImeans_table.mat')
TE=8; 
x = TE:TE:TE*32; y = ROImeansTable.T2Decay{120,1}  ; %ROImeansTable.T2Decay{10,1}; 
figure; plot(x, y,'ko-', 'MarkerFaceColor','k','MarkerSize',11)
set(gcf, 'color', 'w'); set(gca, 'fontsize',18)
xlabel('Time (ms)')
ylabel('Amplitude (a.u.)')
xlim([0 max(x)+10]); ylim([min(y)-100 max(y)+100])
savePath = '/Thesis/Figures/MRI/Pixel_by_Pixel/Criteria/Synthetic_Curves/png'; 
saveas(gcf, fullfile(savePath, 'decay_stimulatedEchoEffect.eps'),'epsc')
saveas(gcf, fullfile(savePath, 'decay_stimulatedEchoEffect.png'))