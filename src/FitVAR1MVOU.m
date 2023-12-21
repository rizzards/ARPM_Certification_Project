function [output1, output2, output3, output4] = FitVAR1MVOU(dX, X, tau, p, nu, lambda_beta, lambda_phi, flag_rescale, flag_filter, method)
% This function estimates the 1-step parameters of the MVOU process (levels) via lasso regression (on first differences)
%  INPUTS
%   dX               : [matrix] (n_ x t_ ) historical series of dependent variables
%   X                : [matrix] (n_ x t_) historical series of independent variables
%   tau              : [scalar] time between subsequent observations in the time series X, dX
%   p                : [vector] (1 x t_) flexible probabilities
%   nu               : [scalar] degrees of freedom of multivariate Student t
%   lambda_beta      : [scalar] lasso regression parameter for loadings
%   lambda_phi       : [scalar] lasso regression parameter for covariance matrix
%   flag_rescale     : [boolean flag] if 0 (default), the series is not rescaled before estimation
%   flag_filter      : [boolean flag] if 0 (default), no filter is imposed in
%                                    the theta-eigenvalues. If 1, the eigenvalues with negative real part
%                                    are set to 0
%   method           : [string] if method == MVOU (default), the function maps VAR1 parameters to MVOU parameters
%  OUTPUTS
%   output1          : if method == MVOU, [vector](n_ x 1) output1 = mu;  else if method == VAR1, [matrix](n_ x n_) output1 = beta
%   output2          : if method == MVOU, [matrix](n_ x n_) output2 = theta;  else if method == VAR1, [matrix](n_ x 1) output2 = alpha 
%   output3          : if method == MVOU, [matrix](n_ x n_) output3 = sig2;  else if method == VAR1, [matrix](n_ x n_) output3 = sigma2
%   output4          : if method == MVOU, [vector](1 x 2) relative error for computation of theta and sigma2; else if method == VAR1, output4 = []

%% Code

if nargin < 10 || isempty(method)
    method = 'MVOU';
end

if nargin < 9 || isempty(flag_filter)
    flag_filter = 0; % no filter on theta-eigenvalues
end

if nargin < 8 || isempty(flag_rescale)
    flag_rescale = 0; % no rescaling of the series
end

% input for lasso regression
if nargin < 7 || isempty(lambda_phi)
    lambda_phi = 0;
end

if nargin < 6 || isempty(lambda_beta)
    lambda_beta = 0;
end

if nargin < 5 || isempty(nu)
    nu = 10^9; % normal invariants
end

[n_, t_] = size(dX);
if nargin < 4 || isempty(p)
    p = ones(1,t_)/t_;
end

% robust lasso + glasso regression
[alpha, beta, sig2_U] = RobustLassoFPReg(dX, X, p, nu, 10^-6, lambda_beta, lambda_phi, flag_rescale);

switch method
    
    case 'VAR1'
        
        output1 = (eye(n_)+beta);
        output2 = alpha;
        output3 = sig2_U;
        output4 = [];
        
    case 'MVOU'
        
        % computation of theta
        Tol_eig = 10^-10;
        [e,Diag_lambda2] = eig(beta);
        beta_mod = beta;
        while ne(rank(beta_mod),rank(Diag_lambda2))
            beta_mod = beta + randn(n_)*Tol_eig;
            [e,Diag_lambda2] = eig(beta_mod);
        end
        lambda2_beta = diag(Diag_lambda2);% eigenvalues of beta
        % clean the small eigenvalues
        lambda2_beta(abs(lambda2_beta)<Tol_eig) = 0;
        Diag_lambda2 = diag(lambda2_beta);
        if norm(beta,'fro') == 0
            r(1) = 0;
        else
            r(1) = norm(e*Diag_lambda2/e-beta,'fro');
        end        
        theta_diag = -log(lambda2_beta + 1)/tau;
        % filter the negative eigenvalues of theta
        if (flag_filter == 1)&&(~isempty(find(real(theta_diag) < 0,1)))
            theta_diag(real(theta_diag) < 0) = 0;
        end
        theta = e*diag(theta_diag)/e;
        theta = real(theta);
        
        % computation of mu
        f = zeros(n_,1);
        f(theta_diag<=Tol_eig) = 1/tau;
        f(theta_diag>Tol_eig) = theta_diag(theta_diag>Tol_eig)./(1-exp(-theta_diag(theta_diag>Tol_eig)*tau));
        mu = e*diag(f)/e*alpha;
        mu = real(mu);
        
        % computation of sig2
        kronsum = kron(theta,eye(n_)) + kron(eye(n_),theta);
        [e_kron, Diag_lambda2_kron] = eig(kronsum);
        lambda2_kron = diag(Diag_lambda2_kron);
        Diag_lambda2_kron = diag(lambda2_kron);
        if norm(kronsum,'fro') == 0
            r(2) = 0;
        else
            r(2) = norm(e_kron*Diag_lambda2_kron/e_kron-kronsum,'fro');
        end        
        lambda_a = NaN(length(Diag_lambda2_kron),1);
        vec_sig2_U = reshape(sig2_U, n_^2,1);
        lambda_a((abs(lambda2_kron) <= Tol_eig)) = 1/tau;
        index = abs(lambda2_kron) > Tol_eig;
        lambda_a(index) = lambda2_kron(index)./(1-exp(-lambda2_kron(index)*tau));
        a = e_kron*diag(lambda_a)/e_kron;
        vec_sigma = a*vec_sig2_U;
        sig2 = reshape(vec_sigma,n_,n_);
        sig2 = real(sig2);
        
        % outputs
        output1 = mu;
        output2 = theta;
        output3 = sig2;
        output4 = r;
        
end

