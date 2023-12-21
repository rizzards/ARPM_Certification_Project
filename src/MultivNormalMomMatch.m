function [X_, p] = MultivNormalMomMatch(mu_, sigma2_, j_, method, d)
% This function generates antithetic normal simulations whose 
% moments match the theoretical moments
%  INPUTS
%   mu_     : [vector] (n_ x 1) vector of means
%   sigma2_ : [matrix] (n_ x n_) dispersion matrix
%   j_      : [scalar] (even) number of simulations
%   method  : [string] Riccati (default), CPCA, PCA, Cholesky-LDL, Gram-Schmidt
%   d       : [matrix] (k_ x n_) full rank constraints matrix for CPCA
%  OUTPUTS
%   X_      : [matrix] (n_ x j_) normal matrix of simulations
%   p       : [vector] (1 x j_) vector of Flexible Probabilities
% NOTE: Use always a large number of simulations j_ >> n_ to ensure that 
% MultivNormalMomMatch works properly

%% Code

if nargin<5 || isempty(d); d = []; end

if nargin<4 || isempty(method); method = 'Riccati'; end

n_ = length(mu_); % number of variables

% Step 1. t MC scenarios 
X_tilde = mvnrnd(zeros(1,n_), eye(n_), j_/2);
X_tilde = X_tilde';

% Step 2. Anthitetic (mean = 0)
X = [X_tilde, -X_tilde];
p = ones(1,j_)/j_; % flat probabilities

% Step 3. Twisted t MC scenarios
X_ = TwistScenMomMatch(X, p, mu_, sigma2_, method, d);