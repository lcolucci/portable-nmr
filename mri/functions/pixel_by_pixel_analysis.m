%% pixel_by_pixel_analysis
%  + This script does pixel-by-pixel fittings (e.g. 1-exp, 2-exp, ILT fits, etc.)
%  + You can use this script on a single patient, or put it inside a loop and use it to analyze many patients at once (this is how 'process_mri.m' uses this function)
%  + This function calls on fitting subfunctions like t2_fitting_2exp_table, t2_fitting_1exp_table, etc. You can add as many as you want
%  + This function gets rid of saturated pixel values before performing T2 fittings


function [pixelT2Decays, saturatedPixels, result_1exp, result_2exp, result_1expOffset, result_2expOffset,result_ILT] = pixel_by_pixel_analysis(image, mask, TE, ignore, doPixelT2Decays, doSaturatedPixels, do1ExpFit, do2ExpFit,do1ExpFitOffset,do2ExpFitOffset,doILT, startVals1Exp,startVals2Exp,startVals1ExpOffset,startVals2ExpOffset)

%% Do Checks and Define Constants
if length(size(image))==4
    nSlices = size(image, 3); 
    nTimePts = size(image,4); 
end
if length(size(image))==3 %nifti conversion not in freesurfer
    nSlices = 1; 
    nTimePts = size(image,3); 
    mask = mask(:,:,1); 
end
if (length(TE) ~= nTimePts)
    error('The time array (TE) has a length of %d. It should have a length of %d, the same length as the # of points on MRI data. Fix this in the USER INPUTS section (''timeLoop'').',length(TE), nTimePts)
end
saturatedIndex = 0; 

%% Loop through slices
for iSlice = 1:nSlices
    fprintf('Processing Slice %d of %d (%s) \n',iSlice,nSlices,datestr(now))
    
    %% Identify row and col IDs for mask
    [rowIdsInMask, colIdsInMask] = find(mask(:,:,iSlice)); % row and col IDs in the mask (i.e. where the 1's are located)
    if isempty(rowIdsInMask) %If there's no mask on this slice, go on to next slice
        continue
    end

    %% Pre-allocate Tables
    sliceName = strcat(sprintf('slice%d', iSlice));
    nMaskPixels = sum(sum(mask(:,:,iSlice))); 
    % pixelT2Decays
    if (doPixelT2Decays ==1)
        pixelT2Decays.(sliceName) = table(nan(nMaskPixels,1),nan(nMaskPixels,1),nan(nMaskPixels,1),nan(nMaskPixels,1),cell(nMaskPixels,1));
        pixelT2Decays.(sliceName).Properties.VariableNames = {'row','col','slice','pixelNum','decay'};
    end
    % saturatedPixels
    if (doSaturatedPixels==1)
        saturatedPixels.(sliceName) = table(nan(nMaskPixels,1),nan(nMaskPixels,1),nan(nMaskPixels,1),nan(nMaskPixels,1),nan(nMaskPixels,1),cell(nMaskPixels,1),cell(nMaskPixels,1),cell(nMaskPixels,1),...
                                      'VariableNames',{'row','col','slice','pixelNum','nSaturatedPts','idSaturatedPts','correspondingTE','decay'});
    end
    % Fittings: pre-allocate tables
    nCOLUMNS_EXTRA = 4; %slice, pixelNum, row, col
    if (do1ExpFit == 1)
        nCOLUMNS_1EXPFIT = 11; 
        result_1exp.(sliceName) = array2table(nan(nMaskPixels, nCOLUMNS_1EXPFIT+nCOLUMNS_EXTRA), 'VariableNames',{'relax1','relax1Std','amp1','amp1Std','r2','sse','rmse','maxVal','exitflag','startValAmp1','startValRelax1','slice','pixelNum','row','col'});
    end
    if (do2ExpFit == 1)
        nCOLUMNS_2EXPFIT = 17; 
        result_2exp.(sliceName) = array2table(nan(nMaskPixels, nCOLUMNS_2EXPFIT+nCOLUMNS_EXTRA), 'VariableNames',{'relax1','relax1Std','relax2','relax2Std','amp1','amp1Std','amp2','amp2Std','r2','sse','rmse','maxVal','exitflag','startValAmp1','startValRelax1','startValAmp2','startValRelax2','slice','pixelNum','row','col'});
    end
    if (do1ExpFitOffset == 1)
        nCOLUMNS_1EXPFITOFFSET = 14; 
        result_1expOffset.(sliceName) = array2table(nan(nMaskPixels, nCOLUMNS_1EXPFITOFFSET+nCOLUMNS_EXTRA), 'VariableNames',{'relax1','relax1Std','amp1','amp1Std','offset','offsetStd','r2','sse','rmse','maxVal','exitflag','startValAmp1','startValRelax1','startValOffset','slice','pixelNum','row','col'});
    end
    if (do2ExpFitOffset == 1)
        nCOLUMNS_2EXPFITOFFSET = 20; 
        result_2expOffset.(sliceName) = array2table(nan(nMaskPixels, nCOLUMNS_2EXPFITOFFSET+nCOLUMNS_EXTRA), 'VariableNames',{'relax1','relax1Std','relax2','relax2Std','amp1','amp1Std','amp2','amp2Std','offset','offsetStd','r2','sse','rmse','maxVal','exitflag','startValAmp1','startValRelax1','startValAmp2','startValRelax2','startValOffset','slice','pixelNum','row','col'});
    end
    if (doILT ==1)
        t2Min=0.01; t2Max=1000; num_guesses=100; t2Vector =logspace(log10(t2Min),log10(t2Max),num_guesses);
        lambdaValues = [0.001, 0.01, 0.1,0.11,0.12, 0.13];
        nLAMBDA_VALUES = length(lambdaValues); 
%         iltEmptyStruct = struct('lambda',nan,'xSpec',nan,'ySpec',nan,'ySpec_norm',nan,'ySpec_norm_integral',nan,'x_Lcurve',nan,'y_Lcurve',nan,'resnorm',nan,'residual',nan,'exitflag',nan,'output',struct(),'lagrangeMultipliers',nan,'date','');
        iltEmptyTable = table(nan(nLAMBDA_VALUES,1), repmat(struct(),nLAMBDA_VALUES,1),'VariableNames',{'lambda','ilt'}); % I don't think I need to pre-allocate the struct b/c that gets done in the ilt_multiplelambdas function. But maybe I'm wrong. 
        result_ILT.(sliceName) = table(nan(nMaskPixels,1),nan(nMaskPixels,1),nan(nMaskPixels,1),nan(nMaskPixels,1),repmat({iltEmptyTable},nMaskPixels,1),'VariableNames',{'slice','pixelNum','row','col','ilt'})
    end
        
        
    %% Loop through each pixel
    for jPixels = 1:nMaskPixels
        
        %% Calculate T2 Decay for this pixel
        if length(size(image))==3 %nifti conversion not in freesurfer
            tempPixelDecay(1:32) = image(rowIdsInMask(jPixels), colIdsInMask(jPixels), :); 
        else
            tempPixelDecay(1:32) = image(rowIdsInMask(jPixels), colIdsInMask(jPixels), iSlice,:); % 32 intensity values per pixel
        end
        %% Prepare x (time) and y (amplitude) data
        x_millisec = TE(:); % better name than TE?
        y = tempPixelDecay(:); 
        % initial points to ignore
        x_millisec = x_millisec(1+ignore:end); 
        y = y(1+ignore:end); 
        % check for any saturated values and remove them 
        [x_millisec, y] = remove_saturatedpixels(x_millisec, y); 
        
        %% pixelT2Decays: Generate table of T2 decays for each pixel 
        if (doPixelT2Decays ==1)
            pixelT2Decays.(sliceName).row(jPixels,1) = rowIdsInMask(jPixels); 
            pixelT2Decays.(sliceName).col(jPixels,1) = colIdsInMask(jPixels);
            pixelT2Decays.(sliceName).slice(jPixels,1) = iSlice;
            pixelT2Decays.(sliceName).pixelNum(jPixels,1) = jPixels;
            pixelT2Decays.(sliceName).decay(jPixels,1) = {tempPixelDecay}; 
        else
            pixelT2Decays = NaN;  % need to give it a value b/c I haven't designed this fxn with ability to do a variable # of outputs
        end
        
        %% saturated pixels: Check if there are any 4095 (maxed out pixel intensities) and put them in a structure so I know where they are
        if (doSaturatedPixels==1)
            SATURATED_PIXELVAL = 4095; %Scanner has 12-bit encoding, 2^12 = 4096-1 = 4095
            if sum(tempPixelDecay==SATURATED_PIXELVAL) > 0 
                saturatedIndex = saturatedIndex +1;
                saturatedPixels.(sliceName).row(saturatedIndex,1) = rowIdsInMask(jPixels); 
                saturatedPixels.(sliceName).col(saturatedIndex,1) = colIdsInMask(jPixels);
                saturatedPixels.(sliceName).slice(saturatedIndex,1) = iSlice;
                saturatedPixels.(sliceName).pixelNum(saturatedIndex,1) = jPixels;
                saturatedPixels.(sliceName).nSaturatedPts(saturatedIndex,1) = sum(tempPixelDecay==SATURATED_PIXELVAL); 
                saturatedPixels.(sliceName).idSaturatedPts(saturatedIndex,1) = {find(tempPixelDecay==SATURATED_PIXELVAL)}; 
                saturatedPixels.(sliceName).correspondingTE(saturatedIndex,1) = {TE(tempPixelDecay==SATURATED_PIXELVAL)}; 
                saturatedPixels.(sliceName).decay(saturatedIndex,1) = {tempPixelDecay};
            end
        else
            saturatedPixels=NaN; 
        end       
        
        %% Perform Fittings
        % 1-exp
        if (do1ExpFit ==1)
            try
                result_1exp.(sliceName)(jPixels, 1:nCOLUMNS_1EXPFIT) = t2_fitting_1exp(x_millisec, y,startVals1Exp); 
                result_1exp.(sliceName).slice(jPixels) = iSlice; 
                result_1exp.(sliceName).pixelNum(jPixels) = jPixels;
                result_1exp.(sliceName).row(jPixels) = rowIdsInMask(jPixels);
                result_1exp.(sliceName).col(jPixels) = colIdsInMask(jPixels);
            catch % if fit can't converge and function throws error, that pixel's fit columns will just stay NaN
                result_1exp.(sliceName).slice(jPixels) = iSlice; 
                result_1exp.(sliceName).pixelNum(jPixels) = jPixels;
                result_1exp.(sliceName).row(jPixels) = rowIdsInMask(jPixels);
                result_1exp.(sliceName).col(jPixels) = colIdsInMask(jPixels);
            end
        else
            result_1exp = NaN; %clear result_1exp %delete the empty initialized table
        end
        
        % 2-exp
        if (do2ExpFit == 1)
            try 
                result_2exp.(sliceName)(jPixels, 1:nCOLUMNS_2EXPFIT) = t2_fitting_2exp(x_millisec, y, startVals2Exp); 
                result_2exp.(sliceName).slice(jPixels) = iSlice; 
                result_2exp.(sliceName).pixelNum(jPixels) = jPixels;
                result_2exp.(sliceName).row(jPixels) = rowIdsInMask(jPixels);
                result_2exp.(sliceName).col(jPixels) = colIdsInMask(jPixels);
            catch % if fit can't converge and function throws error, that pixel's fit columns will just stay NaN
                result_2exp.(sliceName).slice(jPixels) = iSlice; 
                result_2exp.(sliceName).pixelNum(jPixels) = jPixels;
                result_2exp.(sliceName).row(jPixels) = rowIdsInMask(jPixels);
                result_2exp.(sliceName).col(jPixels) = colIdsInMask(jPixels);
            end
        else
            result_2exp=NaN; %delete the empty initialized table
        end
        
        % 1-exp with offset
        if (do1ExpFitOffset ==1)
            try
                result_1expOffset.(sliceName)(jPixels, 1:nCOLUMNS_1EXPFITOFFSET) = t2_fitting_1exp_offset(x_millisec, y,startVals1ExpOffset); 
                result_1expOffset.(sliceName).slice(jPixels) = iSlice; 
                result_1expOffset.(sliceName).pixelNum(jPixels) = jPixels;
                result_1expOffset.(sliceName).row(jPixels) = rowIdsInMask(jPixels);
                result_1expOffset.(sliceName).col(jPixels) = colIdsInMask(jPixels);
            catch % if fit can't converge and function throws error, that pixel's fit columns will just stay NaN
                result_1expOffset.(sliceName).slice(jPixels) = iSlice; 
                result_1expOffset.(sliceName).pixelNum(jPixels) = jPixels;
                result_1expOffset.(sliceName).row(jPixels) = rowIdsInMask(jPixels);
                result_1expOffset.(sliceName).col(jPixels) = colIdsInMask(jPixels);
            end
        else
            result_1expOffset = NaN; %clear result_1exp %delete the empty initialized table
        end
        
        % 2-exp with offset
         if (do2ExpFitOffset == 1)
            try 
                result_2expOffset.(sliceName)(jPixels, 1:nCOLUMNS_2EXPFIT) = t2_fitting_2exp_offset(x_millisec, y, startVals2ExpOffset); 
                result_2expOffset.(sliceName).slice(jPixels) = iSlice; 
                result_2expOffset.(sliceName).pixelNum(jPixels) = jPixels;
                result_2expOffset.(sliceName).row(jPixels) = rowIdsInMask(jPixels);
                result_2expOffset.(sliceName).col(jPixels) = colIdsInMask(jPixels);
            catch % if fit can't converge and function throws error, that pixel's fit columns will just stay NaN
                result_2expOffset.(sliceName).slice(jPixels) = iSlice; 
                result_2expOffset.(sliceName).pixelNum(jPixels) = jPixels;
                result_2expOffset.(sliceName).row(jPixels) = rowIdsInMask(jPixels);
                result_2expOffset.(sliceName).col(jPixels) = colIdsInMask(jPixels);
            end
        else
            result_2expOffset=NaN; %delete the empty initialized table
        end
        
        % 1-exp minus rician
        
        % 2-exp minus rician
        
        % ILT
        if (doILT ==1)
            try
                result_ILT.(sliceName).slice(jPixels) = iSlice; 
                result_ILT.(sliceName).pixelNum(jPixels) = jPixels;
                result_ILT.(sliceName).row(jPixels) = rowIdsInMask(jPixels);
                result_ILT.(sliceName).col(jPixels) = colIdsInMask(jPixels);
                result_ILT.(sliceName).ilt(jPixels,1) = {ilt_multiplelambdas(x_millisec, y, t2Vector, lambdaValues)}; %when to create t2Vector?
            catch
            end
        else
            result_ILT = NaN; 
        end
        %% <<<<<<<<<<<<<<<<<<<<<<<<<< Add other fit types later

        clear tempPixelDecay
    end
    
    % Delete any purely NaN rows of saturatedPixels.(sliceName)
    if (doSaturatedPixels==1)
        try
            saturatedPixels.(sliceName)= rmmissing(saturatedPixels.(sliceName)); %works for R2017a
        catch
            saturatedPixels.(sliceName)(isnan(saturatedPixels.slice1.row),:) = []; %works for R2016a
        end
    end
end


