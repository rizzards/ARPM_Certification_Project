function b2 = Riccati(phi2, sigma2)
% This function solves the algebraic Riccati equation s2 = b2 * phi2 * b2'
% by means of the Schur decomposition
%  INPUTS
%   sigma2 : [matrix] (n_ x n_) symmetric and positive (semi)definite matrix
%   phi2   : [matrix] (n_ x n_) symmetric and positive (semi)definite matrix
%  OUTPUTS
%   b2     : [matrix] (n_ x n_) symmetric matrix positive (semi)definite matrix

%% Code

n_ = length(sigma2);

% 1. Block matrix
h = [zeros(n_,n_), -phi2; -sigma2, zeros(n_,n_)];

% 2. Schur decomposition
[u, t] = schur(h, 'real');
u = ordschur(u, t, 'lhp');

% 3. Four n_ x n_ partitions
u_oneone = u(1 : n_, 1 : n_);
u_twoone = u(n_+1 : end, 1 : n_);

% Output
b2 = u_twoone / u_oneone;

