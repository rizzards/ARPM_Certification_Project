function [mu_MLFP,sigma2_MLFP,error] = MaxLikelihoodFPLocDispT(epsi,nu,p,threshold,last)
% This function estimates the Maximum Likelihood with Flexible Probabilities location and dispersion parameters of the invariants under the Student t distribution assumption by
% means of an iterative algorithm (MaxLikelihoodFPLocDispT routine)
%  INPUT
%  epsi            :[matrix](i_ x t_) timeseries of invariants
%  nu              :[scalar] degrees of freedom of the t-Student distribution
%  p               :[vector](1 x t_) flexible probabilities associated to the invariants
%  threshold       :[scalar] or [vector](1 x 2) convergence threshold
%  last  (optional):[scalar] if last~=0 only the last computed mean and covariance are returned
%  OUTPUT
%  mu_MLFP      :[matrix](i_ x k_) array containing the mean vectors computed at each iteration
%  sigma2_MLFP  :[array](i_ x i_ x k_) array containing the covariance matrices computed at each iteration
%  error       :[vector](2 x 1) vector containing the relative errors at the last iteration

%% Code

if numel(threshold) == 1
    threshold = [threshold threshold];
end

% if last == 0, all the intermediate steps estimates are returned (DEFAULT); if last~=0 only the final estimates are returned.
if nargin < 5 || isempty(last); last = 0; end

% initialize
[i_,t_] = size(epsi);
mu_MLFP(:,1) = epsi*p';
epsi_c = epsi-repmat(mu_MLFP(:,1),1,t_);
sigma2_MLFP(:,:,1) = epsi_c*diag(p)*epsi_c';

error = 10^6;
k = 1;
while sum(error>threshold) >= 1
    k = k+1;
    % update weigths
    epsi_c = epsi-repmat(mu_MLFP(:,k-1),1,t_);
    w_den = nu+sum(epsi_c.*(sigma2_MLFP(:,:,k-1)\epsi_c),1);
    w = (nu+i_)./w_den;
    % update output
    mu_MLFP(:,k) = sum(repmat(p.*w,i_,1).*epsi,2)/(p*w');
    epsi_c = epsi-repmat(mu_MLFP(:,k),1,t_);
    sigma2_MLFP(:,:,k) = epsi_c*diag(p.*w)*epsi_c';
    % convergence
    error(1) = norm(mu_MLFP(:,k)-mu_MLFP(:,k-1))/norm(mu_MLFP(:,k-1));
    error(2) = norm(sigma2_MLFP(:,:,k)-sigma2_MLFP(:,:,k-1),'fro')/norm(sigma2_MLFP(:,:,k-1),'fro');
end

if last ~= 0
    mu_MLFP = mu_MLFP(:,end);
    sigma2_MLFP = sigma2_MLFP(:,:,end);
end
