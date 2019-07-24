%% Calculate ecdf (empirical cumulative distribution function) at specified values
%  This built-in Matlab function 'ecdf' cannot be calculated at specific
%  intervals (edges). This function does that. 
%
%   INPUTS
%       y - vector of values that you want to calculate CDF for
%       edges - vector of values that you want to evaluate the CDF at
%       startID - optional, index of where to start integration. default = 1.  
%       endID - optional, index of where to end integration. default = to the end of the vector.
%
%   OUTPUTS
%       cumDistVector - a vector with the CDF of 'y' evaluated at 'edges'. length = length(edges)-1

function cumDistVector = ecdf2(y, edges, startID, endID)

%% Checks
% if nargin 
if ~exist('startID') 
    startID = 1; 
end
if ~exist('endID')
    endID = length(edges); 
end
y = y(:); 
edges = edges(:); 

%% Calculate CDF
% Calculate bincounts, equivalent to bar heights of  "histogram(y,edges)" 
[bincounts, ~] = histcounts(y, edges);  % length(bincounts) = length(edges)-1
% Integrate across the bounds
cumDistVector = cumtrapz( edges(startID+1:endID), bincounts(startID:endID-1) ); 
% Normalize by max value so that max value of CDF is 1
cumDistVector = cumDistVector ./ max(cumDistVector); 
