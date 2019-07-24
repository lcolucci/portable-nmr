%% Generate Noise Lookup Table
%  This function generates a lookup table that allows you to go from a curve's noise floor stdev
%  (noiseStd), which is how the MR world calculates SNR, to a power level (noisePower),
%  which is how the function 'wgn' works. y = wgn(m,n,power)
%
%  INPUTS
%       minPower - minimum power level (i.e. from wgn(m,n,power))
%       maxPower - max power level 
%       nPts - # of points btw minPower and maxPower (linearly spaced). 
%              this corresponds to # of rows of lookup table
%
%  OUTPUTS
%       noiseLookupTable = table with 2 columns and # of rows = nPts 
%                           + wgnPower (the power level to use in wgn function)
%                           + noiseStd (denominator of SNR formula SNR = maxVal/noiseStd)
%                          

function noiseLookupTable = generate_noiselookuptable(minPower, maxPower, nPts)

noisePower = linspace(minPower, maxPower, nPts); % <--- If stdError is large, change the bounds of this array to reduce it. 
                                                 % amp=100 or less, linspace(-100, 10, 1000) 
                                                 % amp=3000 or so, linspace(-10, 100, 1000)
noiseLookupTable = table(nan(length(noisePower),1),nan(length(noisePower),1),'VariableNames',{'wgnPower','noiseStd'} ); %initialize empty table

for jj=1:length(noisePower)
    noiseVector = wgn(1000,1, noisePower(jj), 1, 5); % 1 = default imp value. 5 = seed (can be anything, jut make sure same # is used below)
    noiseStd = std(noiseVector); 
    
    noiseLookupTable.wgnPower(jj) = noisePower(jj); 
    noiseLookupTable.noiseStd(jj) = noiseStd; 
end    