%% rsquared_groups
%  This function loops through each group of a dataset (like gscatter does, for example) 
%  and finds the linear r2 value for each group 
%
%  see more properties under Matlab's 'LinearModel class' documentation page
%
%  Lina A. Colucci


function rSquared = rsquared_groups(x, y, groups)

%% Checks
if size(x)~=size(y)
    error('x and y must be the same size')
end

% If no 'groups' are entered as a function input, the script will just calculate 
% the r2 value for  the full x,y arrays
if ~exist('groups')
    groups = ones(size(x)); 
end

% If group is a cell array of grouping variables (i.e.{grouping1,grouping2}
% cell arrays must have same # of rows <-- write check for this later
if exist('groups') && max(size(groups))>1
    [groups,groupCombos,groupCombosMatrix,ignore2,maxgrp] = statslib.internal.mgrp2idx(groups,size(x,1),','); 
end

%% Identify Groups
listGroups = unique(groups); 
nGroups = size(listGroups,1); 

for iGroup = 1:nGroups
    whichRows = find(groups==listGroups(iGroup)); 
    tempX = x(whichRows);
    tempY = y(whichRows);
    
    rSquared{iGroup} = rsquared(tempX, tempY); 
end



