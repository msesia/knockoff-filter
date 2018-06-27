% This demo illustrates the basic usage of the Knockoffs package 
% on a synthetic data set.

%% Synthetic problem parameters

n = 1000;         % Number of data points
p = 1000;         % Number of variables
k = 60;           % Number of variables with nonzero coefficients
amplitude = 4.5;  % Magnitude of nonzero coefficients
sigma = 1;        % Noise level
q = 0.10;         % Target false discovery rate (FDR)

rng(123);         % Random seed

%% Synthetic problem construction
% We generate the data by sampling the rows of X from a multivariate normal
% distribution with mean zero and identity covariance matrix.
% Conditional on X, the response y is drawn from a linear regression model
% with k non-zero coefficients

Sigma = eye(p);
mu = zeros(1,p);
X = mvnrnd(mu, Sigma, n);

S0 = randsample(p,k);
beta = zeros(p,1);
beta(S0) = amplitude/sqrt(n);
sampleY = @(X) X*beta + sigma .* randn(n,1);

trueDiscoveries = @(S) sum(beta(S) > 0);
power = @(S) trueDiscoveries(S)*100/k;
FDP = @(S) sum(beta(S) == 0) / max(1, length(S));
printSummary = @(S) fprintf(...
    ['%d true discoveries (Power = %2.2f%%)\n' ...
     'FDP = %2.2f%% (target FDR = %2.f%%)\n'], ...
    trueDiscoveries(S), power(S), 100*FDP(S), 100*q);

%% Running the knockoff filter
% Here we call the knockoff filter with all the default settings. We will
% explore some variations below.

y = sampleY(X);
Xmodel = {'gaussian', mu, Sigma};
S = knockoffs.filter(X, y, q, Xmodel);
printSummary(S);

%% Using a different method for creating knockoff variables
% By default, knockoff variables are created by solving a semi-definite 
% program (SDP). It is also possible to create equi-correlated knockoff
% variables at a lower computational cost.

S = knockoffs.filter(X, y, q, Xmodel, 'Method', 'equi');
printSummary(S);

%% Using yet a different method for creating knockoff variables
% SDP knockoffs may be too expensive to compute when the dimension p is large.
% It is also possible to create optimized knockoff variables that are the 
% solution to a simpler semi-definite programming (SDP) problem obtained by 
% approximating the covariance matrix as a block diagonal.

S = knockoffs.filter(X, y, q, Xmodel, 'Method', 'ASDP');
printSummary(S);


%% Using a different test statistic
% By default, a test statistic based on the lasso is used. Here we use
% a different statistic based on Lasso with cross-validation.

S = knockoffs.filter(X, y, q, Xmodel, 'Statistics', @knockoffs.stats.lassoLambdaSignedMax);
printSummary(S);

%% Using a custom test statistic
% It is also possible to define your own test statistic. To illustrate
% this, we implement a very simple statistic from the knockoff paper.

myKnockoffStatistic = @(X, X_ko, y) ...
    abs(X' * y) - abs(X_ko' * y);

S = knockoffs.filter(X, y, q, Xmodel, 'Statistics', myKnockoffStatistic);
printSummary(S);

%% Using another custom test statistic
% As another example, we show how to change the number of lambda values
% used to approximate the lasso path in the default test statistic
% (cf. the documentation for knockoffs.stats.lassoSignedMax).

myLassoStatistc = @(X, X_ko, y) ...
    knockoffs.stats.lassoLambdaSignedMax(X, X_ko, y, 10*p);

S = knockoffs.filter(X, y, q, Xmodel, 'Statistics', myLassoStatistc);
printSummary(S);

%% Using knockoffs for correlated Gaussian variables
% We generate the data by sampling the rows of X from a AR(1) process with 
% non-zero correlation coefficient.
% Conditional on X, the response y is drawn from a linear regression model
% with k non-zero coefficients.

rho = 0.1;
Sigma = toeplitz(rho.^(0:(p-1)));
Xmodel = {'gaussian', mu, Sigma};

beta(S0) = amplitude;
sampleY = @(X) X*beta + sigma .* randn(n,1);
X = mvnrnd(mu, Sigma, n);
y = sampleY(X);

S = knockoffs.filter(X, y, q, Xmodel);
printSummary(S);

%% Using second-order approximate knockoffs
% We generate the data by sampling the rows of X from non-Gaussian
% distribution.
% Conditional on X, the response y is drawn from a linear regression model
% with k non-zero coefficients.
% Knockoffs are created according to a multivariate normal distribution
% fitted to the observed covariate matrix X.

X = floor(mvnrnd(mu, Sigma, n)*2)/2;
y = sampleY(X);

muHat = mean(X);
SigmaHat = 0.75*eye(p) + 0.25*cov(X);
Xmodel = {'gaussian', muHat, SigmaHat};
X_k = knockoffs.create(X, Xmodel{:});

S = knockoffs.filter(X, y, q, Xmodel);
printSummary(S);
