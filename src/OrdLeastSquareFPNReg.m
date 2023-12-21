function [alpha_OLSFP, beta_OLSFP, s2_OLSFP, U] = OrdLeastSquareFPNReg(X, Z, p)
% This function computes the Ordinary Least Square with Flexible 
% Probabilities (OLSFP) estimator of loadings and dispersion of a 
% regression LFM.
%  INPUTS
%   X            :[matrix] (n_ x t_) time-series of target variables
%   Z            :[matrix] (k_ x t_) time-series of factors
%   p            :[vector] (1 x t_) flexible probabilities
%  OUTPUTS
%   beta_OLSFP   :[matrix] (n_ x k_) OLSFP estimator of loadings
%   s2_OLSFP     :[matrix] (n_ x n_) OLSFP estimator of dispersion of residuals
%   alpha_OLSFP  :[vector] (n_ x 1) OLSFP estimator of the shifting term
%   U            :[matrix] (n_ x t_) time-series of fitted residuals

%% code
[n_,t_] = size(X);
k_ = size(Z, 1);
[m_XZ, s2_XZ] = FPmeancov([X; Z], p);
s_XZ = s2_XZ(1 : n_, n_+1 : n_+k_);
s_ZZ = s2_XZ(n_+1 : n_+k_, n_+1 : n_+k_);

beta_OLSFP = s_XZ / s_ZZ;

alpha_OLSFP = m_XZ(1:n_) - beta_OLSFP*m_XZ(n_+1:n_+k_);

U = X - repmat(alpha_OLSFP,1,t_) - beta_OLSFP * Z;
[~, s2_OLSFP] = FPmeancov(U, p);