%% Check if relax 1 > relax 2. Switch the entries (and their corresponding amplitudes) if so. 
%  Input
%       result - A table from 2exp fit (i.e. t2_fitting_2exp_table.m) - must have relax1,relax2,amp1,amp2 and their corresponding stdevs
%                Table can have 1 or many rows in it. 
%  Output
%       resultChecked - the same table with component 1 and component 2 in the right order (in case it wasn't already)

function resultChecked = checkorder_2expfit(result)
    resultChecked = result; 
    
    if sum(resultChecked.relax1 > resultChecked.relax2)>0
        ids = find(resultChecked.relax1 > resultChecked.relax2); 
        
        % for = 1:length(ids) <-- I think this is how I would add multi-row capability to this code
        
        % Temporary variables to store results of 1st component
        tempRelax1 = resultChecked.relax1(ids); 
        tempRelax1Std = resultChecked.relax1Std(ids); 
        tempAmp1 = resultChecked.amp1(ids); 
        tempAmp1Std  = resultChecked.amp1Std(ids); 

        % Replace component 1 with component 2 values
        resultChecked.relax1(ids) = resultChecked.relax2(ids); 
        resultChecked.relax1Std(ids) = resultChecked.relax2Std(ids); 
        resultChecked.amp1(ids) = resultChecked.amp2(ids); 
        resultChecked.amp1Std(ids) = resultChecked.amp2Std(ids); 

        % replace component 2 values with component 1 values
        resultChecked.relax2(ids) = tempRelax1; 
        resultChecked.relax2Std(ids) = tempRelax1Std; 
        resultChecked.amp2(ids) = tempAmp1; 
        resultChecked.amp2Std(ids) = tempAmp1Std; 
    end
end