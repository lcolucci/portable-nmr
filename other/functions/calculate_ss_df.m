function [ss, df,result] = calculate_ss_df(time, data, whichModels, startVals)

% Checks
% lenght of time = data
% how many models? cell array 


% Go through each model type
if whichModels.do1exp == 1
    try
        [result.exp1,model] = t2_fitting_1exp(time, data, startVals.exp1); 
        ss.exp1 = result.exp1.sse; 
        nParams = length(coeffnames(model)); 
        df.exp1 = length(data) - nParams; 
        clear model nParams
    catch
        ss.exp1 = NaN; 
        df.exp1 = NaN; 
        disp('Fit: 1 exponential could not converge properly.')
    end
end

if whichModels.do1expoffset == 1
    try 
        [result.exp1offset,model] = t2_fitting_1exp_offset(time, data, startVals.exp1offset); 
        % ALTERNATIVE [~,model1o] = t2_fitting_1exp_offset(time, data, startVals.exp1offset); fit1o = fresult(time); ss1o = sum((data(:)-fit1o(:)).^2); % same thing as the sse of the model. So don't need these lines
        ss.exp1offset = result.exp1offset.sse; 
        nParams = length(coeffnames(model));  
        df.exp1offset = length(data) - nParams; 
        clear  model nParams
    catch
        ss.exp1offset = NaN; 
        df.exp1offset = NaN; 
        disp('Fit: 1 exponential w/ offset could not converge properly.')
    end
end
   
if whichModels.do2exp == 1
    try 
        [result.exp2,model] = t2_fitting_2exp(time, data, startVals.exp2); 
        ss.exp2 = result.exp2.sse; 
        nParams = length(coeffnames(model));  
        df.exp2 = length(data) - nParams; 
        clear  model nParams
    catch 
        ss.exp2 = NaN; 
        df.exp2 = NaN; 
        disp('Fit: 2 exponential could not converge properly.')
    end
end

if whichModels.do2expoffset == 1
    try 
        [result.exp2offset,model] = t2_fitting_2exp_offset(time, data, startVals.exp2offset); 
        ss.exp2offset = result.exp2offset.sse; 
        nParams = length(coeffnames(model)); 
        df.exp2offset = length(data) - nParams; 
        clear  model nParams
    catch 
        ss.exp2offset = NaN; 
        df.exp2offset = NaN; 
        disp('Fit: 2 exponential w/ offset could not converge properly.')
    end
end

if whichModels.do3exp == 1
    try 
        [result.exp3,model] = t2_fitting_3exp(time, data, startVals.exp3); 
        ss.exp3 = result.exp3.sse; 
        nParams = length(coeffnames(model)); 
        df.exp3 = length(data) - nParams; 
        clear  model nParams
    catch 
        ss.exp3 = NaN; 
        df.exp3 = NaN; 
        disp('Fit: 3 exponential could not converge properly.')
    end
end

if whichModels.do3expoffset == 1
    try 
        [result.exp3offset,model] = t2_fitting_3exp_offset(time, data, startVals.exp3offset); 
        ss.exp3offset = result.exp3offset.sse; 
        nParams = length(coeffnames(model)); 
        df.exp3offset = length(data) - nParams; 
        clear  model nParams
    catch 
        ss.exp3offset = NaN; 
        df.exp3offset = NaN;
        disp('Fit: 3 exponential w/ offset could not converge properly.')
    end
end

if whichModels.do4exp == 1
    try 
        [result.exp4,model] = t2_fitting_4exp(time, data, startVals.exp4); 
        ss.exp4 = result.exp4.sse; 
        nParams = length(coeffnames(model)); 
        df.exp4 = length(data) - nParams; 
        clear  model nParams
    catch 
        ss.exp4 = NaN; 
        df.exp4 = NaN;
        disp('Fit: 4 exponential could not converge properly.')
    end
end

if whichModels.do5exp == 1
    try 
        [result.exp5,model] = t2_fitting_5exp(time, data, startVals.exp5); 
        ss.exp5 = result.exp5.sse; 
        nParams = length(coeffnames(model)); 
        df.exp5 = length(data) - nParams; 
        clear  model nParams
    catch 
        ss.exp5 = NaN; 
        df.exp5 = NaN;
        disp('Fit: 5 exponential could not converge properly.')
    end
end

