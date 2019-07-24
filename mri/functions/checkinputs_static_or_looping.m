%% checkinputs_static_or_looping
%  + This is a wrapper function around the base function 'check_cellarray_length.m', which determines 
%    whether a variable (an input to 'process_mri.m') should be used across multiple subjects (if length=1) or a different input for each subject (if length=# subject IDs). 
%  + This function just applies that base function to a variable number of inputs. 
%
%  Inputs
%       iLoopingVar = the iterating variable in the 'process_mri.m' script that loops through the patient ids (i.e. 'jPatient')
%       varargin = all the variable names that you want to go through this check (they must all be cells/cell arrays)
%
%  Outputs
%       varargout = all the variable names that were given as inputs. 
%                   Either the input itself (if its length was 1) or the iLoopingVar-th element of the input (if its length was > 1)

function varargout = checkinputs_static_or_looping(iLoopingVar, varargin)

nFIXED_INPUTS = 1; % change this # if more variables are added in front of 'varargin'

for iInputs = nFIXED_INPUTS+1:nargin %loop through varargin only
    
    iInputsVarargin = iInputs-nFIXED_INPUTS; % Use this looping variable when indexing the variable 'varargin' or 'varargout'
                                             % Use 'iInputs' for indexing everything else
                                             
    varargout{iInputsVarargin} = check_cellarray_length(iLoopingVar, varargin{iInputsVarargin}); 

end