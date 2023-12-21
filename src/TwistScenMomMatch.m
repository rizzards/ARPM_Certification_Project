function x_ = TwistScenMomMatch(x, p, mu_, s2_, method, d)
% This function twists scenarios x to match arbitrary moments mu_ sigma2_
%  INPUTS
%   x      : [matrix] (n_ x j_) scenarios
%   p      : [vector] (1 x j_) flexible probabilities
%   mu_    : [vector] (n_ x 1) target means
%   s2_    : [matrix] (n_ x n_) target covariances
%   method : [string] Riccati (default), CPCA, PCA, LDL-Cholesky, Gram-Schmidt
%   d      : [matrix] (k_ x n_) full rank constraints matrix for CPCA
%  OUTPUTS
%   x_     : [matrix] (n_ x j_) twisted scenarios

%% Code

if nargin<6 || isempty(d); d = []; end

% default method: Riccati
if nargin<5 || isempty(method); method = 'Riccati'; end

% Step 1. Original moments
[mu_x, s2_x] = FPmeancov(x, p);

% Step 2. Transpose-square-root of s2_x
r_x = TransposeSquareRoot(s2_x, method, d);

% Step 3. Transpose-square-root of s2_
r_ = TransposeSquareRoot(s2_, method, d);

% Step 4. Twist factors
b = r_/r_x;

% Step 5. Shift factors
a = mu_ - b*mu_x;

% Step 6. Twisted scenarios
x_ = repmat(a, 1, size(x, 2)) + b*x;

