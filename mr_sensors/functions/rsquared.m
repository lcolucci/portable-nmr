%% rsquared
%  This function calculates the ordinary r2 value of a linear fit for an x and y vector
%
%  See more properties for which you could create functions like this under
%  Matlab's 'LinearModel class' documentation page (i.e. RMSE, residuals, etc.)
%  
%  Lina A. Colucci

function r2 = rsquared(x,y)

if size(x)~=size(y)
    error('x and y must be the same size')
end

mdl = fitlm(x, y); 
r2 = mdl.Rsquared.Ordinary; 
