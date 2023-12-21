function [x, u, U] = CMASep(X, p)
% This function computes the Separation step of Copula-Marginal Algorithm 
%  INPUTS
% X : [matrix] (n_ x j_) FP-joint scenarios
% p : [vector] (1 x j_) Flexible Probabilities 
%  OUTPUTS
% x : [matrix]  (n_ x j_) sorted scenarios: x(n, j) <= x(n, j+1)
% u : [matrix]  (n_ x j_) sorted cumulative probabilities: u(n, j) <= u(n, j+1)
% U : [matrix]  (n_ x j_) FP-copula scenarios 

% Preprocess variables
[n_, j_] = size(X);

 if size(p, 1) == n_
     for n = 1 : n_
         p(n, :) = max(p(n, :), 1/j_ * 10^(-8));
         p(n, :) = p(n, :)/sum(p(n, :));
     end
 else
    p = max(p, 1/j_ * 10^(-8)); 
    p = p/sum(p);
 end
 %% Code
 if j_<=10000
     l = j_/(j_ + 0.001); % to be fixed
 else
     l = j_/(j_ + 1);
 end

[x, Indx] = sort(X, 2); % sort the rows of X

u = NaN(n_, j_);
U = NaN(n_, j_);

for n = 1 : n_ 
    I = Indx(n, :); 
    cum_p = cumsum(p(I)); % cumulative probabilities
    u(n, :) = cum_p * l; % rescale to be <1 at the far right
    [~, Rnk] = sort(I); % compute ranking
    U(n, :) = cum_p(Rnk) * l; % grade scenarios
end