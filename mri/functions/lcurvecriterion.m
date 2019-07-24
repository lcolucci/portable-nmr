function [lambdas, eta, rho, kappa, chosenLambda, chosenRho, chosenEta, s, betaCut, xi] = lcurvecriterion(time, data, t2Vector, plotLcurve,plotCurvature)
% find lambda
% plot L-curve and kappa curve (and where the ideal lamda value is located)



%% Double check that time and data are columns or rows (whatever it needs to be)
%% double check that time is in ms (or whatever it needs to be)
% maybe I just pass it the xSpec array already generated. NEed to figure
% that out. 


% Make time and data a column vector if not already
time = time(:); 
data = data(:); 

%% Construct an array of logarithmically spaced possible T2 values
% xSpec = logspace(log10(t2Min),log10(t2Max),nPts);

%% Generate A matrix (contains set of all possible T2 decays based on the min and max T2 defined as an input to the equation)
A = exp(-time*(1./t2Vector));

%% Perform singular value decomposition (SVD) on the A matrix
%   csvd is a function in the regtools (Hansen's matlab toolbox). It is an svd function (which is built in to Matlab) but with some checks and changes depending on the size of the matrix.
[U,s,V] = csvd(A,'full');

%% Intermediate Calculations
% Initialize variables and calculate lengths
b = data; % Let's call data 'b' to be consistent with Hansen's nomenclature
[m,n] = size(U); 
[p, ps] = size(s); 
% npoints = length(s);
npoints = 500; % # of points on L-curve. Can edit this # if you want more or less points on the curve.  
eta = zeros(npoints,1); rho=eta; lambdas = eta; % intialize vectors with 0's

% Calculate intermediate variables
beta = U'*b; % equation 15 in Hansen paper
betaCut = beta(1:p); % Trim beta vector to be the same size as 's' vector
xi = betaCut ./s; xi( isinf(xi) ) = 0; 
s2 = s.^2; 

%% Create vector of regularization parameters (lambda's)
smin_ratio = 16*eps;  % Smallest regularization parameter.
lambdas(npoints) = max([s(p), s(1)*smin_ratio]); % fill the last value in the lambda vector with the smallest possible lambda
ratio = (s(1)/lambdas(npoints))^(1/(npoints-1)); 
for i=npoints-1:-1:1
    lambdas(i) = ratio*lambdas(i+1); % fill in each successive value in the lambda vector by multiplying the smallest possible lambda by the ratio
end

%% Calculate eta, rho, and phi (i.e. eta prime) vectors
for i = 1:npoints
    f = s2 ./(s2 +lambdas(i)^2); % Equation 6 in Hansen's paper
    eta(i) = norm(f.*xi); % Equation 12 (left) in Hansen's paper
    rho(i) = norm((1-f).*betaCut); % Equation 12 (right) in Hansen's paper
    
    phi(i) = sum(-4.*f.*(1-f).*f.*xi.^2/lambdas(i));  % Equation 15 (left) in Hansen's paper. Phi (in code) = eta prime (in paper)
                     % Off by factor of 2 from lcfun's phi, but this expression matches eq. 15 in the paper. I chose to match what the paper says.   
                     % I DO NOT understand when to use (i) and when not to??                                    
end

%% Calculate Kappa (equation 18 in Hansen's paper)
for i = 1:npoints
    % Note: I broke up the expression for kappa into 3 parts (a,b,c), otherwise the line of code would have been too long and hard to understand
    a = eta.*rho./phi'; % left-most fraction, eq.18
    b = lambdas.^2.*rho.*phi' + 2.*lambdas(i).*rho.*eta + lambdas.^4.*eta.*phi'; % top of right-most fraction, eq.18
    c = (lambdas.^4.*eta.^2 + rho.^2).^1.5; % bottom of right-most fraction, eq.18
    %c2 = (reg_param.^2.*eta.^2 + rho.^2).^1.5; % This is the lambda^2 that Hansen's paper has for eq. 18 but I showed this is wrong by
                                                % calculating kappa with this c2 value. It doesn't even have the right shape (just looks like an exponential decay)

    kappa = 2.*a.*b./c; %  put parts of eq.18 together
end

%% Calculate the proper lambda value (corner of the L-curve)
[kappaMax, kappaMaxID] = max(kappa);
chosenLambda = lambdas(kappaMaxID);
chosenRho = rho(kappaMaxID); 
chosenEta = eta(kappaMaxID); 

%% Plot L-Curve
if plotLcurve==1
    figure; %plot(log10(rho), log10(eta), 'o-'); 
    loglog(rho, eta,'o-')
    xlabel('Residual Norm ||Ax-b|| (rho)'); ylabel('Solution Norm ||x|| (eta)')
    title('L-curve')
    hold on;
    loglog(chosenRho, chosenEta, '*','MarkerSize',20)
    %plot(log10(chosenRho), log10(chosenEta),'*','MarkerSize',20)
    set(gcf,'color','w'); set(gca,'fontsize',16)
    annotation('textbox','String',sprintf('Corner of L-curve: Lambda = %.3f',chosenLambda),'LineStyle','none','FontSize',16)
end

%% Plot Kappa (Curvature)
if plotCurvature==1
    figure; plot(kappa, 'LineWidth',3); 
    set(gcf,'color','w'); set(gca,'fontsize',16)
    hold on;
    plot(kappaMaxID, kappaMax,'*','MarkerSize',20)
    title('Curvature Plot of L-curve (kappa)')
    annotation('textbox','String',sprintf('Max kappa is %.2f at ID # %.0f, \nwhich corresponds to lambda = %.3f',kappaMax, kappaMaxID, chosenLambda),'LineStyle','none','FontSize',16)
end

%figure; plot(eta, rho,'o-')

% save rho, eta, kappa, lambda. L-curve plot with chosen point (both jpeg and fig). Kappa plot (both jpeg and fig). 
% create prettyfigure function that white background, big font, no line around legend 
% create save figure function that saves things as both png and fig ?

