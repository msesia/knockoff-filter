% This demo illustrates some advanced usage of the Knockoffs package 
% on a synthetic data set.

%% Synthetic problem parameters

n = 1000;         % Number of data points
p = 1000;         % Number of variables
k = 60;           % Number of variables with nonzero coefficients
amplitude = 15;   % Magnitude of nonzero coefficients
sigma = 1;        % Noise level
q = 0.10;         % Target false discovery rate (FDR)

rng(123);         % Random seed

%% Synthetic problem construction
% We generate the data by sampling the rows of X from a multivariate normal
% distribution with mean zero and identity covariance matrix.
% Conditional on X, the response y is drawn from a logistic regression model
% with k non-zero coefficients

Sigma = eye(p);
mu = zeros(1,p);

S0 = randsample(p,k);
beta = zeros(p,1);
beta(S0) = amplitude/sqrt(n);
invlogit = @(x) exp(x)./(1+exp(x));
sampleY = @(X) binornd(1,invlogit(X*beta));

trueDiscoveries = @(S) sum(beta(S) > 0);
power = @(S) trueDiscoveries(S)*100/k;
FDP = @(S) sum(beta(S) == 0) / max(1, length(S));
FDPp = @(S) sum(beta(S) == 0) / (1/q + length(S));
printSummary = @(S) fprintf(...
    ['%d true discoveries (Power = %2.2f%%)\n' ...
     'FDP = %2.2f%% (target FDR = %2.f%%)\n'], ...
    trueDiscoveries(S), power(S), 100*FDP(S), 100*q);

%% Running the knockoff filter steps manually
% The main function 'knockoffs.filter' is a wrapper around simpler functions
% that create knockoffs, compute test statistics, and perform variable
% selection. When more control is necessary, these functions may be
% called directly. We demonstrate this below in reproducing the plot of
% Figure 1.

X = mvnrnd(mu, Sigma, n);
y = sampleY(X);
Xmodel = {'gaussian', mu, Sigma};
X_k = knockoffs.create(X, Xmodel{:});
[W,Z] = knockoffs.stats.lassoLambdaSignedMax_bin(X, X_k, y);
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
title('Knockoff Filter with Lasso Statistics');
legend('Null feature', 'Non-null feature');
line = refline(1,0);
set(line, 'LineStyle', ':', 'Color', 'black');


%% Running the knockoff filter manually multiple times
% The function 'knockoffs' is also a  wrapper around simpler 
% functions that create knockoffs for a specific model.
% When knockoffs for multiple identically distributed realizations of the 
% data matrix X are needed, one can save computation time by precomputing 
% the diagonal entries of 'diag_s' (e.g. by solving the SDP, for SDP
% knockoffs). We demonstrate this below.

diag_s = sparse(diag(knockoffs.create.solveSDP(Sigma)));
m=10;
[fdp, fdpp, pwr] = deal(nan(m,1));

for i = 1:m
    X = mvnrnd(mu, Sigma, n);
    y = sampleY(X);
    X_k = knockoffs.create.gaussian_sample(X, mu, Sigma, diag_s);
    W = knockoffs.stats.lassoCoefDiff_bin(X, X_k, y);
    S = knockoffs.select(W, q, 'knockoff');
    fdp(i) = FDP(S); fdpp(i) = FDPp(S); pwr(i) = trueDiscoveries(S)/k;
end
fprintf('Mean FDP: %2.2f, Mean FDPp: %2.2f, Mean Power: %2.2f%%\n', mean([fdp fdpp pwr*100]))
boxplot([fdp, fdpp, pwr], 'Labels',{'Fdp', 'Fdpp', 'Power'})
