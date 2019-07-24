function [result,fresult] = t2_fitting_1exp_offset(x,y,StartingPoints)
  % First order exponential fitting 
  % 
  % USAGE:
  %   [resultMatrix, gof] = exponential_fit_offset(data, StartPoint)
  %   
  %   resultMatrix = matrix with columns [Amp1, 95%CI_Amp1, Tau1, 95%CI_Tau1]
  %   gof = goodness of fit 
  %   data = [time in 1st column, amplitudes in 2nd column]
  %
  %   C. Vassiliou December 2006.
  %   L. Colucci September 2015
  %   2015-09-30 Edited to make it useful for MRI analysis
  
  % Load data
  x = x(:); y=y(:); 

  % Define the model and options
  model = fittype('a*exp(-x/b)+c', 'independent','x','dependent','y'); % This is the model 
  options = fitoptions('Method','NonLinearLeastSquares');
  Confidence = 0.95;
  if (~exist('StartingPoints','var') || isempty(StartingPoints))
      StartingPoints=[3,80,10]; % If no starting values entered, use default input
  end 
  set(options,'StartPoint',StartingPoints,'Lower',[0,0,0], 'Upper', [Inf,10000,10000]);

 % Perform the fit
  [fresult,gof,out] = fit(x,y,model,options);

  % Calculate error in the parameters given a confidence interval
  ci = confint(fresult,Confidence);

  % Put results in a table
  result = table(); 
  result.relax1 =  fresult.b;
  result.relax1Std = ci(2,2)'-fresult.b; 
  result.amp1 = fresult.a; 
  result.amp1Std = ci(2,1)'-fresult.a;
  result.offset = fresult.c; 
  result.offsetStd = ci(2,3)'-fresult.c; 
  result.r2 = gof.rsquare; 
  result.sse = gof.sse; 
  result.rmse = gof.rmse; 
  result.maxVal = max(y); 
  result.exitflag = out.exitflag; 
  result.startValAmp1 = StartingPoints(1); 
  result.startValRelax1 = StartingPoints(2); 
  result.startValOffset = StartingPoints(3); 
  
  