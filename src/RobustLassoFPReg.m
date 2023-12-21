function [alpha_RMLFP,beta_RMLFP,sig2_RMLFP] = RobustLassoFPReg(X,Z,p,nu,tol,lambda_beta,lambda_phi,flag_rescale)
% Robust Regression - Max-Likelihood with Flexible Probabilites & Shrinkage
% (multivariate Student t distribution with given degrees of freedom = nu)
%  INPUTS
%   X             : [matrix] (n_ x t_ ) historical series of dependent variables
%   Z             : [matrix] (k_ x t_) historical series of independent variables
%   p             : [vector] flexible probabilities
%   nu            : [scalar] multivariate Student's t degrees of freedom
%   tol           : [scalar] or [vector] (3 x 1) tolerance, needed to check convergence
%   lambda_beta  : [scalar] lasso regression parameter
%   lambda_phi    : [scalar] graphical lasso parameter
%   flag_rescale  : [boolean flag] if 0 (default), the series is not rescaled
%
%  OUTPUTS
%   alpha_RMLFP   : [vector] (n_ x 1) shifting term
%   beta_RMLFP    : [matrix] (n_ x k_) optimal loadings
%   sig2_RMLFP    : [matrix] (n_ x n_) matrix of residuals' covariances 

%% Code
[n_, t_] = size(X);
k_ = size(Z,1);

% if flag_rescale is not provided it is set to 0
if nargin < 8 || isempty(flag_rescale)
    flag_rescale = 0;
end
% if lambda_lasso or lambda_phi are not provided they are set to 0
if nargin < 7 || isempty(lambda_phi)
    lambda_phi=0;
end
if nargin < 6 || isempty(lambda_beta)
    lambda_beta = 0;
end
% if FP are not provided, observations are equally weighted
if nargin<3 || isempty(p)
    p = ones(1,t_)/t_;
end
% adjust tolerance input
if numel(tol)==1 || numel(tol)==2
    tol = [tol(1) tol(1) tol(1)];
end

% rescale variables
if flag_rescale == 1
    [~,cov_Z]=FPmeancov(Z, p);
    sig_Z = sqrt(diag(cov_Z));
    [~,cov_X]=FPmeancov(X,p);
    sig_X = sqrt(diag(cov_X));
    Z = diag(1./sig_Z) * Z;
    X = diag(1./sig_X) * X;
end

% initialize variables
alpha = nan(n_,1);
beta = nan(n_,k_,1);
sig2 = nan(n_,n_,1);


% 0. Initialize
[alpha(:,1), beta(:,:,1), sig2(:,:,1), U] = LassoRegFP(X, Z, p, 0);

error = [10^6 10^6 10^6];
maxIter = 500;
i = 1;
while (sum(error>tol) >=1 && (i < maxIter))
    i = i+1;
    
    % 1. Update weights
    z2 = U'*(sig2(:,:,i-1)\U);
    w = (nu+n_)./(nu+diag(z2)');
    
    % 2. Update FP
    p_tilde = (p.*w) / sum(p.*w);
    
    % 3. Update output
    % Lasso regression
    [alpha(:,i),beta(:,:,i),sig2(:,:,i),U] = LassoRegFP(X,Z,p_tilde,lambda_beta);
    sig2(:,:,i) = sum(p.*w) * sig2(:,:,i);
    % Graphical lasso 
    if lambda_phi ~= 0
        [sig2(:,:,i),~,~,~,~]=GraphLasso(sig2(:,:,i),lambda_phi);
    end
    
    % 3. Check convergence
    error(1) = norm(alpha(:,i)-alpha(:,i-1))/norm(alpha(:,i-1));
    error(2) = norm(beta(:,:,i)-beta(:,:,i-1),'fro')/norm(beta(:,:,i-1),'fro');
    error(3) = norm(sig2(:,:,i)-sig2(:,:,i-1),'fro')/norm(sig2(:,:,i-1),'fro');
end

% Output
alpha_RMLFP = alpha(:,end);
beta_RMLFP = beta(:,:,end);
sig2_RMLFP = sig2(:,:,end);


% From rescaled variables to non-rescaled variables
if flag_rescale == 1
    alpha_RMLFP = diag(sig_X) * alpha_RMLFP;
    beta_RMLFP = diag(sig_X) * beta_RMLFP *diag(1./sig_Z);
    sig2_RMLFP = diag(sig_X) * sig2_RMLFP * diag(sig_X)';
end

