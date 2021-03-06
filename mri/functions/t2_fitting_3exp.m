function [result, fresult] = t2_fitting_3exp(x,y, StartingPoints)
  % L. Colucci September 2015

  % Load data
  x = x(:); y=y(:); 
  
  % Define the model and options
  model = fittype('a*exp(-x/b)+c*exp(-x/d)+e*exp(-x/f)');
  options = fitoptions(model);
  Confidence = 0.95;
  if (~exist('StartPoint','var') || isempty(StartingPoints))
      StartingPoints=[1,50,3,15,1,200]; % If no starting values entered, use this default 
  end 
  set(options,'StartPoint',StartingPoints,'Lower',[0,0,0,0,0,0], 'Upper', [10000,10000,10000,10000,10000,10000]);

  % Perform the fit
  [fresult,gof,out] = fit(x,y,model,options);

  % Calculate error in the parameters given a confidence interval
  ci = confint(fresult,Confidence);
  
  % Put results in a table
  result = table(); 
  result.amp1 = fresult.a;
  result.amp1Std = ci(2,1)'-fresult.a; 
  result.relax1 = fresult.b; 
  result.relax1Std = ci(2,2)'-fresult.b;
  result.amp2 = fresult.c; 
  result.amp2Std = ci(2,3)'-fresult.c;
  result.relax2 = fresult.d ;
  result.relax2Std = ci(2,4)'-fresult.d;
  result.amp3 = fresult.e ;
  result.amp3Std = ci(2,5)'-fresult.e;
  result.relax3 = fresult.f ;
  result.relax3Std = ci(2,6)'-fresult.f;
  result.r2 = gof.rsquare; 
  result.sse = gof.sse; 
  result.rmse = gof.rmse; 
  result.maxVal = max(y); 
  result.exitflag = out.exitflag; 
  result.startValAmp1 = StartingPoints(1); 
  result.startValRelax1 = StartingPoints(2); 
  result.startValAmp2 = StartingPoints(3); 
  result.startValRelax2 = StartingPoints(4); 
  result.startValAmp3 = StartingPoints(5); 
  result.startValRelax3 = StartingPoints(6); 
  
  % Make sure component 1,2 and 3 in the right order (comp1 must correspond to the smallest relaxation time. and comp3 to the largest relax time)
%   result = checkorder_3expfit(result); %<---- SCRIPT NOT READY YET
