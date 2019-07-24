%% Calculate proper leg segmental BI volumes (TBW, ECF, ICF)
%  This script loads in raw leg BI values and applies the proper segmental
%  equations to it. 
%
%  Originally written as script 'calculate_BI_segmental_leg.m' and later
%  turned into this function that's called by script 'process_bioimpedance.m'
%
%  Equations are taken from HYDRA ECF/ICF (Model 4200) OPERATING MANUAL
%  REVISION 1.03. Xitron Technologies, Inc. San Diego, CA, USA.
%
%  INPUTS
%       data - must include the following columns: 
%               wholeBody = double checks that all rows are 0
%               RHOem = resistivity of extracellular fluid (Ohm*cm).
%                       RHOem=males. RHOef=females. (all our subjects are male)
%               CalfElectrodeSpacing - (cm)
%               CalfLength - (cm)
%               C1 = leg circumference 1 (major)
%               C2 = leg circumference 2 (minor... but doesn't actually matter which is which)
%               Re = value from model fitting (Ohm)
%               Ri = value from model fitting (Ohm)
%               columns TBW-FFM1 = algorithm doesn't absolutely need these but code is written 
%                                  to clear out data from those columns
%
%  OUTPUTS
%       data - same output table with proper leg volumes calculated
%               additional columns are TBW,ECF,etc. calculated with *_calfLength (rather than electrode spacing) 
%                           since I'm not sure what the proper calf length to use is
%
%  Lina A. Colucci 2018

function data = calculateLegVolumesBI(data)

%% Calculate Bioimpedance: Segment Volume
% Check that no whole body data is present
data = data(data.WholeBody==0,:); 

% Clear out the TBW/ECF/etc. values for leg b/c automatically-calculated values are wrong
colsTBW = find(strcmpi(data.Properties.VariableNames, 'TBW'));
colsFFM1 = find(strcmpi(data.Properties.VariableNames,'FFM1')); 
data{:,colsTBW:colsFFM1} = NaN; % 

% Set Constants
rho_ecw = data.RHOem; %resistivity of extracellular fluid (Ohm*cm). RHOem=males. RHOef=females.
rho_icw = data.RHOim; % resistivity of intracellular fluid (Ohm*cm). RHOim=males. RHOif=females. 
L = data.CalfElectrodeSpacing; %segment length (cm)
L2 = data.CalfLength;
C1 = data.C1; %segment circumference (cm)
C2 = data.C2; %segment circumference (cm)
R_E = data.Re; % value from model fitting (Ohm)
R_I = data.Ri; % value from model fitting (Ohm)
k_p = rho_icw ./ rho_ecw; % equation B11, Xitron Manual

% Calculate leg ECW
data.ECF = (rho_ecw.^(2/3))./(3.*(4.*pi).^(1/3).*1000) .* ...
                L.*(C1.^2 + C2.^2 + C1.*C2) .* ...
                (L./(C1.*C2.*R_E)).^(2/3); % equation B9, Xitron Manual
data.ECF_calfLength = (rho_ecw.^(2/3))./(3.*(4.*pi).^(1/3).*1000) .* ...
                L2.*(C1.^2 + C2.^2 + C1.*C2) .* ...
                (L2./(C1.*C2.*R_E)).^(2/3); % equation B9, Xitron Manual

% Calculate Leg ICW. Can't isolate ICF from equation, so just solve graphically by finding intersection of the two sides.
for iRow = 1:height(data)
    plotData = data(iRow,:); 
    x = linspace(0,20,1000); % ICF must be a positive number
    y1 = ((plotData.Re + plotData.Ri)./plotData.Ri).*(1+ k_p(1).*x./plotData.ECF); % equation B10, Xitron Manual
    y1_calfLength = ((plotData.Re + plotData.Ri)./plotData.Ri).*(1+ k_p(1).*x./plotData.ECF_calfLength); % equation B10, Xitron Manual
    y2 = (1+x./plotData.ECF).^(5/2); % equation B10, Xitron Manual
    y2_calfLength = (1+x./plotData.ECF_calfLength).^(5/2); % equation B10, Xitron Manual
    try 
        [data.ICF(iRow),yout] = intersections(x,y1,x,y2,1); % equation B10, Xitron Manual
    catch 
        data.ICF(iRow)= NaN; 
    end
    try 
        [data.ICF_calfLength(iRow),yout] = intersections(x,y1_calfLength,x,y2_calfLength,1); % equation B10, Xitron Manual
    catch 
        data.ICF_calfLength(iRow) = NaN; 
    end
    
    clear x y1 y1_calfLength y2 y2_calfLength plotData
end

% Calculate Leg TBW
data.TBW = data.ECF + data.ICF; % equation B13, Xitron Manual
data.TBW_calfLength = data.ECF_calfLength + data.ICF_calfLength; % equation B13, Xitron Manual

% Calculate percentages 
data.TBW1 = data.TBW ./ data.Weight .*100; % TBW1 = TBW/Weight * 100%
data.ECF1 = data.ECF ./ data.TBW .*100; % ECF1 = ECF/TBW * 100%
data.ICF1 = data.ICF ./ data.TBW .*100; % ICF1 = ICF/TBW * 100%
data.ECF1_calfLength = data.ECF_calfLength ./ data.TBW_calfLength .*100; 
data.ICF1_calfLength = data.ICF_calfLength ./ data.TBW_calfLength .*100; 
data.TBW1_calfLength = data.TBW_calfLength ./ data.Weight .*100; 



