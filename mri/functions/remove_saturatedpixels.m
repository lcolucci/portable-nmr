%% remove_saturatedpixels
%  Input: 
%       x, y - vectors of the same length. x = time. y = T2 decay. 
%  Output: 
%       xChecked, yChecked - vectors equal or shorter in length than x,y
%
%  This code finds all saturated pixels present in the T2 decay (y), it
%  then deletes every point up until the last saturated value. It deletes
%  the coresponding pixels on the time (x) array so that these checked vectors can
%  be properly fed into any of the fitting functions. 

function [xChecked, yChecked] = remove_saturatedpixels(x, y)

%% Prep and Check Data
xChecked = x(:); yChecked = y(:); 
saturatedPixelVal = 4095; %Scanner has 12-bit encoding, 2^12 = 4096-1 = 4095. Can't trust pixel intensities at 4095.  

if length(x) ~= length(y)
    error('x and y vectors are not the same length')
end

%% Do the check 
if sum(yChecked==4095) > 0 
    %% This is the analysis to delete just the saturated values, but I think deleting all values up until the last saturated value is the correct thing to do. 
%     saturatedIDs = find(yChecked==4095);
%     xChecked(saturatedIDs) = []; 
%     yChecked(saturatedIDs) = []; 
    
    %% Delete all values up until the last saturated pixel
    maxSaturatedID = max(find(yChecked==4095));
    xChecked(1:maxSaturatedID) = []; 
    yChecked(1:maxSaturatedID) = []; 
end
    



