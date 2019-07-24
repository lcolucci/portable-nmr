%% checkinput_cells
%  This function takes a variable # of inputs and checks if (1) they are cells, 
%  and (2) if the cell arrays are either length 1 or length equal to nIDsInLoop
%  
%  Inputs: 
%       nIDsInLoop: # of IDs that process_mri.m is analyzing (i.e. length of 'whichIdToAnalyzeLoop')
%       varargin: all the variable names that you want to go through this check
%
%       nFIXED_INPUTS: if you add any additional fixed inputs to the function, increase this number accordingly
%
%  Outputs:
%       Error messages if any of the 'varargin' do not make the checks

function checkinputs_cells(nIDsInLoop,varargin)

nFIXED_INPUTS = 1; % # of inputs before the 'varargin'

for iInputs = nFIXED_INPUTS+1:nargin %loop through varargin only
    
    iInputsVarargin = iInputs-nFIXED_INPUTS; % Use this looping variable when indexing the variable 'varargin' 
                                             % Use 'iInputs' for indexing everything else
    
    % 1. Check that input is a cell
    if ~iscell(varargin{iInputsVarargin})
        error('The input ''%s'' is a %s not a cell. Inputs must be single cells or cell arrays.', inputname(iInputs), class(varargin{iInputsVarargin}))
    end
    
    % 2. Check that the input cell array has a length of 1 or nIDsInLoop
    if ~(length(varargin{iInputsVarargin})==1 || length(varargin{iInputsVarargin})==nIDsInLoop)
        error('The input ''%s'' has length %d. The cell array must have a length of 1 or length equal to # of IDs (nIDs=%d).', inputname(iInputs), length(varargin{iInputsVarargin}), nIDsInLoop)
    end
    
    % <<<<<<<<<<<<<< ACTUALLY I DON'T THINK THIS SHOULD BE HERE. IT SHOULD
    % BE A CHECK DONE DURING THE FITTINGS, THAT X AND Y HAVE SAME LENGTH.
    % BECAUSE OTHERWISE I CAN'T CHECK STARTVALs
%     % 3. If cell array of is composed of numerics, check that the vectors have the proper length (i.e. REFERENCE_LENGTH)
%     if (isnumeric(varargin{iInputsVarargin}{1}))
%         nCellElements = length(varargin{iInputsVarargin}); 
%         for iCellElements = 1:nCellElements
%             if (length(varargin{iInputsVarargin}{iCellElements})~=REFERENCE_LENGTH)
%                 error('Element #%d of the input cell array ''%s'' has length %d. This vector must have a length equal to the REFERENCE_LENGTH = %d.', iCellElements, inputname(iInputs), length(varargin{iInputsVarargin}{iCellElements}), REFERENCE_LENGTH)
%             end
%         end
%     end
    
end
