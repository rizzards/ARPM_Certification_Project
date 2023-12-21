function dep = SWDepMeasure(X, probs)
% This function estimates the Schweizer and Wolff measure of dependence
% between two random variables by means of Monte Carlo simulations
%  INPUTS
%   X     : [matrix]  (2 x j_) joint scenarios
%   probs : [vector]  (1 x j_) vector of Flexible probabilities
%  OUTPUTS
%   dep   : [scalar]  Schweizer-Wolff measure estimate

%% Code
[~, ~, U] = CMASep(X, probs); % grades scenarios

j_ = size(X, 2); % number of scenarios

g = NaN(j_, j_);
for i = 1 : j_
    for k = 1 : j_
        g(i, k) = abs(sum(probs .* (U(1, :) <= i/j_) .* (U(2, :) <= k/j_)) - (i * k)/j_^2);
    end
end

dep = 12 * mean(g(:));