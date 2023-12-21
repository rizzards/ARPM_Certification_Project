function q_HFP = HFPquantile(x,conf,p)
% This function computes the quantile of a Flexible probabilities
% distribution
%  INPUTS
%  x            :[vector](1 x t_) scenarios
%  conf         :[vector](1 x n_ql) confidence levels
%  p  (optional):[vector](1 x t_) Flexible Probabilities
%  OUTPUTS
%  q_HFP        :[vector](1 x n_ql) quantiles

%% Code

n_ql = length(conf);
% if the third argument is missing, the Flexible Probabilities are set to be uniform
if nargin < 3 || isempty(p); p = (1/length(x))*ones(1,length(x)); end

[x_sort,y] = sort(x);

p_sort = p(y);

q_HFP = zeros(1,n_ql);
cum = 0;
j = 1;
t = 1;
while j <= n_ql   
    while cum < conf(j) && t <= length(x)
        cum = cum + p_sort(t);
        t = t+1;
    end
    if t == 1
        q_HFP(j) = x_sort(1);
    else
        q_HFP(j) = x_sort(t-1);
    end
    j = j+1;
end
