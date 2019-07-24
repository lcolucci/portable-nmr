%% check_cellarray_length
%   This function determines whether an input has length 1 (if so, we use
%   the same input across all patients) or a length > 1 (in which case, we
%   use a different input for each patient)
%
%  Inputs
%       input = a cell or cell array 
%       iLoopingVar = the iterating variable in the 'process_mri.m' script that loops through the patient ids (i.e. 'jPatient')
%  Output
%       checkedInput = either the input (if its length was 1) or the iLoopingVar-th element of the input (if its length was > 1)
%

function checkedInput = check_cellarray_length(iLoopingVar, input)

if ~iscell(input)
    error('The input is not a cell or cell array')
end

if length(input) == 1
    checkedInput = input{1};
elseif length(input) > 1
    checkedInput = input{iLoopingVar}; 
else
    error('Input does not have a length.')
end


%     
%         
%         if (isnumeric(timeLoop) && length(timeLoop)==32)
%             TE = timeLoop; 
%         elseif (iscell(timeLoop) && length(timeLoop)==1)
%             TE = timeLoop{1}; 
%         else
%             TE = timeLoop{jPatient}; %if timeLoop length > 1, it must be a cell array
%         end 
%         
%         if (ischar(whichImageToAnalyzeLoop))
%             whichImageToAnalyze = whichImageToAnalyzeLoop; 
%         elseif (iscell(whichImageToAnalyzeLoop) && length(whichImageToAnalyzeLoop)==1) 
%             whichImageToAnalyze = whichImageToAnalyzeLoop{1}; 
%         else
%             whichImageToAnalyze = whichImageToAnalyzeLoop{jPatient}; %if whichImageToAnalyzeLoop length > 1, it must be a cell array
%         end
%         
%         % <<<<<<<<<<<<<<<<<<<<<< add same thing for startVals 
%         if (isnumeric(startVals1ExpLoop))
%             startVals1Exp = startVals1ExpLoop; 
%         elseif (iscell(startVals1ExpLoop) && length(startVals1ExpLoop)==1)
%             startVals1Exp = startVals1ExpLoop{1}; 
%         else
%             startVals1Exp = startVals1ExpLoop{jPatient};
%         end
%         
%         if (isnumeric(startVals2ExpLoop))
%             startVals2Exp = startVals2ExpLoop; 
%         elseif (iscell(startVals2ExpLoop) && length(startVals2ExpLoop)==1)
%             startVals2Exp = startVals2ExpLoop{1}; 
%         else
%             startVals2Exp = startVals2ExpLoop{jPatient};
%         end