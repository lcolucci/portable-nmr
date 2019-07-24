%% awgn_snr
%  This function adds white gaussian noise to a curve to achieve a desired
%  SNR level, as defined by SNR = max(y)/stdev(noise floor of y)
%
%   The function 'awgn' exists in Matlab to add wgn to a curve and achieve a certain SNR level, 
%   but it does not calculate SNR the same way we do in the MR world. 
%   I needed a way to reliably add noise to synthetic decays to achieve a desired SNR level, so I wrote this function. 
%
%   INPUTS
%       y = the curve you want to add noise to (vector)
%       targetSNR = the desired SNR you want ySNR to have (numeric) 
%
%   OUTPUTS
%       ySNR = the original curve with added white gaussian noise that makes its SNR = targetSNR
%       --- optional ---
%       targetStd = the optimal stdev (noise level) based on max(y)/targetSNR
%       actualStd = the stdev (noise level) actually used in the calculation. Taken from noise lookup table. 
%       error_targetStd_actualStd = difference btw targetStd and actualStd
%       actualWgnPower = the power level that was used in the y = wgn(m,n,power) function to generate this SNR level
%
% 3/25/2018 LAC


function [ySNR, targetStd, actualStd, error_targetStd_actualStd, actualWgnPower] = awgn_snr(y, targetSNR) 

%% Load noise lookup table (the table should be located in same folder as this function)
filepath = mfilename('fullpath'); %filepath 
baseFolder = strrep(filepath, '/awgn_snr',''); 
load(fullfile(baseFolder, 'noiseLookupTable.mat')); %<-- if the error_targetStd_actualStd is too large, generate a new noiseLookupTable with better bounds (generate_noiselookuptable.m)

%% Generate Noisy Data
% Calculate Target Noise Level
yMax = max(y); 
targetStd = yMax / targetSNR; % calculate the target stdev (noise level) to achieve the desired SNR (based on formula SNR = maxVal/stdev_NoiseFloor)

% Find the Target Noise Level in the Noise Lookup Table
[error_targetStd_actualStd index] = min(abs(targetStd - noiseLookupTable.noiseStd)); % find the closest possible stdev match we have in the lookup table
actualStd = noiseLookupTable.noiseStd(index); % this is the closest noise level we could find in the lookup table
actualWgnPower = noiseLookupTable.wgnPower(index); % the 'power' level that corresponds to the 'actualStd' from the lookup table. 
                                               % This is the power level that will be entered into the 'wgn' function to obtain a decay with the desired SNR
if error_targetStd_actualStd > 1 % if stdError > 1, then noise lookup table needs to be re-generated with better noisePower values
    fprintf('\n*** NOTE *** \nThe error btw the target stdev and actual stdev is large. \nGenerate a new noiseLookupTable and \nchange bounds of the noisePower array to fix this.\n************\n')
end

% Generate Noise Vector
noiseVector = wgn(size(y,1),size(y,2), actualWgnPower,1,5); % 1 = default imp value. 5 = seed (can be anything but must be same # as used above in noise lookup table)

% Generate Noisy Data by adding data and noise vector
ySNR = y + noiseVector; 

%% SCRAPS
% awgn(yRaw, 22,'measured','linear'); %<-- This is the built-in Matlab
% function for doing this. supposed to lead to an SNR=22 but the way SNR is
% measured is different from how the MR world does it