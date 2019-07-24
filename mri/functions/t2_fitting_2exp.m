function [result, fresult] = t2_fitting_2exp(x,y, StartingPoints)
  % 2-exponential fit 
  % C. Vassiliou December 2006.
  % L. Colucci September 2015
  % 2017-11-19 LAC  Made output into a table format
 

  % Load data
  x = x(:); y=y(:); 

  % Define the model and options
  model = fittype('a*exp(-x/b)+c*exp(-x/d)');
  options = fitoptions(model);
  Confidence = 0.95;
  if (~exist('StartingPoints','var') || isempty(StartingPoints)) 
      StartingPoints=[1,50,3,100]; % If no starting values entered, use default input
  end 
  set(options,'StartPoint',StartingPoints,'Lower',[0,0,0,0], 'Upper', [10000,10000,10000,10000]);

  % Perform the fit
  [fresult,gof,out] = fit(x,y,model,options);

  % Calculate error in the parameters given a confidence interval
  ci = confint(fresult,Confidence);

  % Put results in a table
  result = table(); 
  result.relax1 =  fresult.b;
  result.relax1Std = ci(2,2)'-fresult.b; 
  result.relax2 =  fresult.d;
  result.relax2Std = ci(2,4)'-fresult.d; 
  result.amp1 = fresult.a; 
  result.amp1Std = ci(2,1)'-fresult.a; 
  result.amp2 = fresult.c; 
  result.amp2Std = ci(2,3)'-fresult.c; 
  result.r2 = gof.rsquare; 
  result.sse = gof.sse; 
  result.rmse = gof.rmse; 
  result.maxVal = max(y); 
  result.exitflag = out.exitflag; 
  result.startValAmp1 = StartingPoints(1); 
  result.startValRelax1 = StartingPoints(2); 
  result.startValAmp2 = StartingPoints(3); 
  result.startValRelax2 = StartingPoints(4); 

  % Make sure component 1 and 2 in the right order (comp1 must correspond to the smaller relaxation time)
  result = checkorder_2expfit(result); 
  
  
  
