function [X_, p] = MultivtMomMatch(nu, mu_, sigma2_, j_, method, d)
% This function generates antithetic student t simulations whose
% moments match the theoretical moments mu_ nu/(nu-2)*sigma2_
%  INPUTS
%   nu      : [scalar] degrees of freedom
%   mu_     : [vector] (n_ x 1) vector of means
%   sigma2_ : [matrix] (n_ x n_) dispersion matrix
%   j_      : [scalar] (even) number of simulations
%   method  : [string] Riccati (default), CPCA, PCA, LDL-Cholesky, Gram-Schmidt
%   d       : [matrix] (k_ x n_) full rank constraints matrix for CPCA
%  OUTPUTS
%   X_      : [matrix] (n_ x j_) student t matrix of simulations
%   p       : [vector] (1 x j_) vector of Flexible Probabilities
% NOTE: Use always a large number of simulations j_ >> n_ to ensure that 
% MultivtMomMatch works properly

%% Code

if nargin<6 || isempty(d); d = []; end

if nargin<5 || isempty(method); method = 'Riccati'; end

% check definitiness covariance matrix
if nu < 2
    error('the covariance matrix is not defined')
end

n_ = length(mu_); % number of variables

% Step 1. t MC scenarios 
X_tilde = mvtrnd(eye(n_), nu, j_/2);
X_tilde = X_tilde';

% Step 2. Anthitetic (mean = 0)
X = [X_tilde, -X_tilde];
p = ones(1,j_)/j_; % flat probabilities

% Step 3. Twisted t MC scenarios
X_ = TwistScenMomMatch(X, p, mu_, nu/(nu-2)*sigma2_, method, d);
