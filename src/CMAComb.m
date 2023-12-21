function X = CMAComb(x, u, U)
% This function computes the Combination step of Copula-Marginal Algorithm  
%  INPUTS
% x : [matrix] (n_ x k_) significant nodes: x(n, k) <= x(n, k+1)
% u : [matrix] (n_ x k_) cdf grid: u(n, k) = Fn(x(n, k))
% U : [matrix] (n_ x j_) FP-copula scenarios 
%  OUTPUTS
% X : [matrix] (n_ x j_) FP-joint scenarios

%% Code
[n_, j_] = size(U);

% joint scenarios by (linear) inter-/extra-polation of the grid (u, x)
X = NaN(n_, j_);
for n = 1 : n_
    X(n, :) = interp1(u(n, :), x(n, :), U(n, :), 'linear', 'extrap');
end
