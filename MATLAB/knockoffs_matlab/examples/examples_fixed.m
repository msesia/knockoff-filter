% This demo illustrates how to use the Knockoffs package 
% with Fixed-X knockoffs on a synthetic data set (i.e. therefore 
% assuming a linear regression model for the response).

%% Synthetic problem parameters

n = 2000;         % Number of data points
p = 700;          % Number of variables
k = 40;           % Number of variables with nonzero coefficients
amplitude = 4.5;  % Magnitude of nonzero coefficients
sigma = 1;        % Noise level
q = 0.10;         % Target false discovery rate (FDR)

rng(123);         % Random seed

%% Synthetic problem construction

X = randn(n,p) / sqrt(n);
S0 = randsample(p,k);
beta = zeros(p,1);
beta(S0) = amplitude;
sampleY = @() X*beta + sigma .* randn(n,1);

trueDiscoveries = @(S) sum(beta(S) > 0);
FDP = @(S) sum(beta(S) == 0) / max(1, length(S));
printSummary = @(S) fprintf(...
    ['%d true discoveries\n' ...
     'FDP = %2.2f%% (target FDR = %2.f%%)\n'], ...
    trueDiscoveries(S), 100*FDP(S), 100*q);

%% Running the knockoff filter

% Here we call the knockoff filter with the default arguments.
% By default, the statistics knockoffs.stats.lassoLambdaSignedMax
% are used with Fixed-X knockoffs.

y = sampleY();
S = knockoffs.filter(X, y, q, {'fixed'});
printSummary(S);

%% Using a different method for creating knockoff variables

% By default, knockoff variables are created by solving a semi-definite 
% programming (SDP) problem. If p is large, this can become computationally
% unfeasible. Therefore, one can alternatively create equi-correlated
% knockoffs without the need to solve an SDP.

S = knockoffs.filter(X, y, q, {'fixed'}, 'Method', 'equi');
printSummary(S);

%% Using yet a different method for creating knockoff variables
% SDP knockoffs may be too expensive to compute when the dimension p is
% large, but equicorrelated knockoffs may sometimes yield very low power.
% It is also possible to create optimized knockoff variables that are the 
% solution to a simpler semi-definite programming (SDP) problem obtained by 
% approximating the Gramm matrix as a block diagonal.

S = knockoffs.filter(X, y, q, {'fixed'}, 'Method', 'ASDP');
printSummary(S);

%% Using a different test statistic
% By default, a test statistic based on the lasso is used. Here we use
% a different statistic based on forward selection.

S = knockoffs.filter(X, y, q, {'fixed'}, 'Statistics', @knockoffs.stats.forwardSelection);
printSummary(S);

%% Using a custom test statistic

% It is also possible to define your own test statistic. To illustrate
% this, we implement a very simple statistic from the knockoff paper.

myKnockoffStatistic = @(X, X_ko, y) ...
    abs(X' * y) - abs(X_ko' * y);

S = knockoffs.filter(X, y, q, {'fixed'}, 'Statistics', myKnockoffStatistic);
printSummary(S);

% As another example, we show how to change the number of lambda values
% used to approximate the lasso path in the default test statistic
% (cf. the documentation for knockoffs.stats.lassoLambdaSignedMax).

myLassoStatistc = @(X, X_ko, y) ...
    knockoffs.stats.lassoLambdaSignedMax(X, X_ko, y, p);

S = knockoffs.filter(X, y, q, {'fixed'}, 'Statistics', myLassoStatistc);
printSummary(S);

%% Running the knockoff filter steps manually

% The main function 'knockoffs.filter' is a wrapper around simpler functions
% that create knockoffs, compute test statistics, and perform variable
% selection. When more control is necessary, these functions may be
% called directly. We demonstrate this below in reproducing the plot of
% Figure 1.

X_ko = knockoffs.create.fixed(X);
[W,Z] = knockoffs.stats.lassoLambdaSignedMax(X, X_ko, y);
t = knockoffs.threshold(W, q);

fig = figure();
hold on
set(fig, 'DefaultTextInterpreter', 'latex');
gscatter(Z(1:p), Z(p+1:2*p), ismember(1:p, S0), 'kr');
plot([t t 0], [0 t t], 'k');
hold off

xlabel('Value of $\lambda$ when $X_j$ enters model');
ylabel('Value of $\lambda$ when $\tilde X_j$ enters model');
limits = [0 ceil(max(Z))];
xlim(limits); ylim(limits);
title('Knockoff Filter with Lasso Statistic');
legend('Null feature', 'Non-null feature');
line = refline(1,0);
set(line, 'LineStyle', ':', 'Color', 'black');
