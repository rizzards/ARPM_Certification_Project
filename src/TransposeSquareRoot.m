function s = TransposeSquareRoot(sigma2, method, d)
% This function computes the transpose-square-root matrix s of a symmetric
% and positive (semi)definite matrix sigma2 such that sigma2 = s*s'
%  INPUTS
%   sigma2 : [matrix] (n_ x n_) positive definite matrix
%   method : [string] Riccati (default), CPCA, PCA, LDL-Cholesky, Gram-Schmidt
%   d      : [matrix] (k_ x n_) full rank constraints matrix for CPCA
%  OUTPUTS
%   s      : [matrix] (n_ x n_) transpose-square-root of sigma2

%% Code

n_ = length(sigma2);

if nargin<3 || isempty(d); d = []; end

% default method: Riccati
if nargin<2 || isempty(method); method = 'Riccati'; end 

switch method   
    case 'Riccati'
        s = Riccati(eye(n_), sigma2);
    case 'CPCA'
        [lambda2_d, e_d] = ConditionalPC(sigma2, d);
        s = (e_d')\diag(sqrt(lambda2_d));
    case 'PCA'
        [e, lambda2] = pcacov(sigma2);
        s = e*diag(sqrt(lambda2))*e';    
    case 'LDL-Cholesky'
        [l,delta_2] = ldl(sigma2);
        s = l*sqrt(delta_2); 
     case 'Gram-Schmidt'    
        g = GramSchmidt(sigma2);
        s = inv(g)';      
end  