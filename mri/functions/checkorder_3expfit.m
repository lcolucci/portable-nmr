%% Check if relax 1 > relax 2 > relax3. Switch the entries (and their corresponding amplitudes) if so. 
%  Input
%       result - A table from 2exp fit (i.e. t2_fitting_2exp_table.m) - must have relax1,relax2,amp1,amp2 and their corresponding stdevs
%                Table can have 1 or many rows in it. 
%  Output
%       resultChecked - the same table with component 1 and component 2 in the right order (in case it wasn't already)

% CAN ONLY HANDLE ONE ROW AT A TIME BECAUSE ORDER MIGHT BE DIFFERENT FOR
% EACH ONE. IF I WANT TO ADD MULTI-ROW CAPABILITY IN FUTURE, WILL PROB
% INVOLVE SETTING UP 6 DIFFERENT CASES FOR THE 6 POTENTIAL ORDERS OF THE 3
% COMPONENTS 

%% !!!!!!!! THIS FUNCTION DOESN'T WORK YET !!!!!!!!!!!!

function resultChecked = checkorder_3expfit(result)
    resultChecked = result; 
    
    % Determine sort order
    relaxArray = [resultChecked.relax1 resultChecked.relax2 resultChecked.relax3]; 
    [sortedRelax, sortedIDs] = sort(relaxArray);
    % if result is not 1,2,3 then sort them 
    
%     if sum(resultChecked.relax1 > resultChecked.relax2)>0 
%         ids = find(resultChecked.relax1 > resultChecked.relax2); 
        
        % Temporary variables to store results of 1st component
        tempRelax1 = resultChecked.relax1(ids); 
        tempRelax1Std = resultChecked.relax1Std(ids); 
        tempAmp1 = resultChecked.amp1(ids); 
        tempAmp1Std  = resultChecked.amp1Std(ids); 
        
        % Temp vars - 2nd component
        tempRelax2 = resultChecked.relax2(ids); 
        tempRelax2Std = resultChecked.relax2Std(ids); 
        tempAmp2 = resultChecked.amp2(ids); 
        tempAmp2Std  = resultChecked.amp2Std(ids); 
        
        % Temp vars - 3rd component
        tempRelax3 = resultChecked.relax3(ids); 
        tempRelax3Std = resultChecked.relax3Std(ids); 
        tempAmp3 = resultChecked.amp3(ids); 
        tempAmp3Std  = resultChecked.amp3Std(ids); 
        
        

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


