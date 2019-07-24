%% Check fits
% 
%   Order of applying the criteria is important. RMSE cutoff must go first,
%   otherwise you'll have to delte

function result_checked = checkfits(result, TE)

%% Checks
if ~exist('TE')
    TE = 8; % Assume TE=8ms unless otherwise specified. 
            % If it's actually 25.5ms, then criteria will be even more restrictive, which is fine
end
if any(contains(result.Properties.VariableNames,'relax2'))
    is2exp=1; 
else
    is2exp=0; 
end

result_checked = result; 

%% RMSE cutoff
% rmseCutoff = 34.8760; %90th percentile of all RMSE values of 2-exp fits for all patients. prctile(allPixels,[80, 85, 90, 95, 97.5, 99]) = [31.7651, 33.0269, 34.8760, 38.5703, 43.8867, 54.3718]  
rmseCutoff = prctile(result.rmse, 99); % 99% RMSE value of this dataset                   
if any(result_checked.rmse > rmseCutoff)
    ids = result_checked.rmse > rmseCutoff;    
    result_checked.relax1(ids) = NaN; 
    result_checked.relax1Std(ids) = NaN; 
    result_checked.amp1(ids) = NaN; 
    result_checked.amp1Std(ids) = NaN; 
    try 
        result_checked.relax2(ids) = NaN; 
        result_checked.relax2Std(ids) = NaN; 
        result_checked.amp2(ids) = NaN; 
        result_checked.amp2Std(ids) = NaN; 
    catch 
        % if 1-exp fit and component 2 doesn't exist
    end
    result_checked.r2(ids) = 999; 
    result_checked.sse(ids) = 999; 
    result_checked.rmse(ids) = 999; 
    result_checked.criteria_rmseCutoff(ids) = 1; 
    clear ids
end 

%% Stdev = NaN 
if is2exp==1
    if any( isnan(result_checked.relax1Std) | isnan(result_checked.relax2Std) | isnan(result_checked.amp1Std) | isnan(result_checked.amp2Std) )
        ids = isnan(result_checked.relax1Std) | isnan(result_checked.relax2Std) | isnan(result_checked.amp1Std) | isnan(result_checked.amp2Std); 
        result_checked.relax1(ids) = NaN; 
        result_checked.relax1Std(ids) = NaN; 
        result_checked.amp1(ids) = NaN; 
        result_checked.amp1Std(ids) = NaN; 
        result_checked.relax2(ids) = NaN; 
        result_checked.relax2Std(ids) = NaN; 
        result_checked.amp2(ids) = NaN; 
        result_checked.amp2Std(ids) = NaN; 
        result_checked.r2(ids) = 999; 
        result_checked.sse(ids) = 999; 
        result_checked.rmse(ids) = 999; 
        result_checked.criteria_stdNaN(ids) = 1; 
        clear ids
    end
else
    if    any(isnan(result_checked.relax1Std) | isnan(result_checked.amp1Std) )
        ids = isnan(result_checked.relax1Std) | isnan(result_checked.amp1Std); 
        result_checked.relax1(ids) = NaN; 
        result_checked.relax1Std(ids) = NaN; 
        result_checked.amp1(ids) = NaN; 
        result_checked.amp1Std(ids) = NaN; 
        result_checked.r2(ids) = 999; 
        result_checked.sse(ids) = 999; 
        result_checked.rmse(ids) = 999; 
        result_checked.criteria_stdNaN(ids) = 1; 
        clear ids
    end
end
 
%% Minimum Criteria
minimumT2 = TE/2;
if is2exp==1
    if any( result_checked.relax1 < minimumT2 | result_checked.relax2 < minimumT2)
        ids = result_checked.relax1 < minimumT2 | result_checked.relax2 < minimumT2; 
        result_checked.relax1(ids) = NaN; 
        result_checked.relax1Std(ids) = NaN; 
        result_checked.amp1(ids) = NaN; 
        result_checked.amp1Std(ids) = NaN; 
        result_checked.relax2(ids) = NaN; 
        result_checked.relax2Std(ids) = NaN; 
        result_checked.amp2(ids) = NaN; 
        result_checked.amp2Std(ids) = NaN; 
        result_checked.r2(ids) = 999; 
        result_checked.sse(ids) = 999; 
        result_checked.rmse(ids) = 999; 
        result_checked.criteria_min(ids) = 1; 
        clear ids
    end
else
    if   any( result_checked.relax1 < minimumT2 )
        ids = result_checked.relax1 < minimumT2; 
        result_checked.relax1(ids) = NaN; 
        result_checked.relax1Std(ids) = NaN; 
        result_checked.amp1(ids) = NaN; 
        result_checked.amp1Std(ids) = NaN; 
        result_checked.r2(ids) = 999; 
        result_checked.sse(ids) = 999; 
        result_checked.rmse(ids) = 999; 
        result_checked.criteria_min(ids) = 1; 
        clear ids
    end
end

%% Maximum Criteria
% % 10ms error criteria
%     % result_checked.maxT2 = result_checked.SNR.*4.789 + 328.6; % relationship btw MaxT2 and SNR determined at very end of script syntheticDecays_explore_T2limits_for_TE_SNR.m
% 5% error criteria
    result_checked.maxT2 = result_checked.SNR .* 25.63 + 197.6; 
if is2exp==1
    if any(result_checked.relax1 > result_checked.maxT2 | result_checked.relax2 > result_checked.maxT2)
        ids = result_checked.relax1 > result_checked.maxT2 | result_checked.relax2 > result_checked.maxT2; 
        result_checked.relax1(ids) = NaN; 
        result_checked.relax1Std(ids) = NaN; 
        result_checked.amp1(ids) = NaN; 
        result_checked.amp1Std(ids) = NaN; 
        result_checked.relax2(ids) = NaN; 
        result_checked.relax2Std(ids) = NaN; 
        result_checked.amp2(ids) = NaN; 
        result_checked.amp2Std(ids) = NaN; 
        result_checked.r2(ids) = 999; 
        result_checked.sse(ids) = 999; 
        result_checked.rmse(ids) = 999; 
        result_checked.criteria_max(ids) = 1; 
        clear ids
    end
else
    if any(result_checked.relax1 > result_checked.maxT2 )
        ids = result_checked.relax1 > result_checked.maxT2 ; 
        result_checked.relax1(ids) = NaN; 
        result_checked.relax1Std(ids) = NaN; 
        result_checked.amp1(ids) = NaN; 
        result_checked.amp1Std(ids) = NaN; 
        result_checked.r2(ids) = 999; 
        result_checked.sse(ids) = 999; 
        result_checked.rmse(ids) = 999; 
        result_checked.criteria_max(ids) = 1; 
        clear ids
        end
end
   

%% Relaxation times are less than 10ms apart, i.e. abs(Relax1 - Relax 2) < 10ms 
closenessCutoff = 10; 
if is2exp==1
    if any( abs(result_checked.relax1 - result_checked.relax2) < closenessCutoff)
        ids = abs(result_checked.relax1 - result_checked.relax2) < closenessCutoff; 
        result_checked.relax1(ids) = NaN; 
        result_checked.relax1Std(ids) = NaN; 
        result_checked.amp1(ids) = NaN; 
        result_checked.amp1Std(ids) = NaN; 
        result_checked.relax2(ids) = NaN; 
        result_checked.relax2Std(ids) = NaN; 
        result_checked.amp2(ids) = NaN; 
        result_checked.amp2Std(ids) = NaN; 
        result_checked.r2(ids) = 999; 
        result_checked.sse(ids) = 999; 
        result_checked.rmse(ids) = 999; 
        result_checked.criteria_closeT2s(ids) = 1; 
        clear ids
    end
end


%% --------------- SCRAPS -------
%% SSE cutoff - same thing as RMSE cutoff
% sseCutoff = 32840; %90th percentile of all SSE values of 2-exp fits for all patients 
% if any(result_checked.sse > sseCutoff)
%     ids = result_checked.sse > sseCutoff;    
%     result_checked.relax1(ids) = NaN; 
%     result_checked.relax1Std(ids) = NaN; 
%     result_checked.amp1(ids) = NaN; 
%     result_checked.amp1Std(ids) = NaN; 
%     result_checked.relax2(ids) = NaN; 
%     result_checked.relax2Std(ids) = NaN; 
%     result_checked.amp2(ids) = NaN; 
%     result_checked.amp2Std(ids) = NaN; 
%     result_checked.r2(ids) = 999; 
%     result_checked.sse(ids) = 999; 
%     result_checked.rmse(ids) = 999; 
%     result_checked.criteria_sseCutoff(ids) = 1; 
%     clear ids
% end 

% %% Relax Stdev > Cutoff
% stdevCutoff = 1300; %99th percentile for all standard deviations for relaxation times across all 2-exponential fits
% if any( result_checked.relax1Std > stdevCutoff | result_checked.relax2Std > stdevCutoff )
%     ids = result_checked.relax1Std > stdevCutoff | result_checked.relax2Std > stdevCutoff; 
%     result_checked.relax1(ids) = NaN; 
%     result_checked.relax1Std(ids) = NaN; 
%     result_checked.amp1(ids) = NaN; 
%     result_checked.amp1Std(ids) = NaN; 
%     result_checked.relax2(ids) = NaN; 
%     result_checked.relax2Std(ids) = NaN; 
%     result_checked.amp2(ids) = NaN; 
%     result_checked.amp2Std(ids) = NaN; 
%     result_checked.r2(ids) = 999; 
%     result_checked.sse(ids) = 999; 
%     result_checked.rmse(ids) = 999; 
%     result_checked.criteria_stdCutoff(ids) = 1; 
%     clear ids
% end    
%     
% try 
% if any( result_checked.relax2Std > stdevCutoff )
%     ids = result_checked.relax2Std > stdevCutoff; 
%     result_checked.relax2(ids) = NaN; 
%     result_checked.relax2Std(ids) = NaN; 
%     result_checked.amp2(ids) = NaN; 
%     result_checked.amp2Std(ids) = NaN; 
%     result_checked.r2(ids) = 999; 
%     result_checked.sse(ids) = 999; 
%     result_checked.rmse(ids) = 999;
%     result_checked.criteria_stdCutoff(ids) = 1; 
%     clear ids
% end
% catch % if relax2 doesn't exist
% end    

%% AmpRatios > Cutoff
% try 
% result_checked.ampRatio = result_checked.amp1 ./ result_checked.amp2;         
% ampRatioCutoff = 100; 
% 
% % if AmpRatio > 100 (remove component 2, which has the smaller amplitude)
% if any( result_checked.ampRatio > ampRatioCutoff )
%     ids = result_checked.ampRatio > ampRatioCutoff; 
%     result_checked.relax2(ids) = NaN; 
%     result_checked.relax2Std(ids) = NaN; 
%     result_checked.amp2(ids) = NaN; 
%     result_checked.amp2Std(ids) = NaN; 
%     result_checked.r2(ids) = 999; 
%     result_checked.sse(ids) = 999; 
%     result_checked.rmse(ids) = 999; 
%     result_checked.criteria_ampRatio(ids) = 1; 
%     clear ids
% end
% 
%  % --- AmpRatio < 0.01 (remove component 1, which has the smaller amplitude)
% if any(result_checked.ampRatio < 1/ampRatioCutoff)
%     ids = result_checked.ampRatio < 1/ampRatioCutoff;    
%     result_checked.relax1(ids) = NaN; 
%     result_checked.relax1Std(ids) = NaN; 
%     result_checked.amp1(ids) = NaN; 
%     result_checked.amp1Std(ids) = NaN; 
%     result_checked.r2(ids) = 999; 
%     result_checked.sse(ids) = 999; 
%     result_checked.rmse(ids) = 999; 
%     result_checked.criteria_ampRatio(ids) = 1; 
%     clear ids
% end 
% 
% % --- If 'criteria_ampRatio' column was never created b/c neither 'if' statementwas triggered, populate it with 0's
% isTableCol = @(t, thisCol) ismember(thisCol, t.Properties.VariableNames);
% if ~isTableCol(result_checked, 'criteria_ampRatio')
%     result_checked.criteria_ampRatio(:) = 0; 
% end
%     
% catch
% end
    

