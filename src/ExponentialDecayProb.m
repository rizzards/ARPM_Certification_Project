function p = ExponentialDecayProb(t_,tau_hl)
%Flexible probabilities computed via exponential decay with half-life
%tau_hl
%INPUTS
% t_ [scalar]: length of the desired vector of probabilities
% tau_hl [scalar]: half-life of the exponential decay
%OUTPUT
%p [vector]: (1 x t_) vector of flexible probabilities obtained via exponential decay
p=exp(-(log(2)/tau_hl)*(t_:-1:1)); 
p=p/sum(p);
end

