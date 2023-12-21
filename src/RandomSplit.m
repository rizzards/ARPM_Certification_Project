function [y1, y2] = RandomSplit(y)

% Splits vector y into two random mutually exclusive partitions of approximately same length.

t_ =length(y);
half_t_ = round(t_/2.0);

% making a random permutation of the given vector y

index = randperm(t_);
y_perm = y(index);

y1 = y_perm(1:half_t_);
y2 = y_perm(half_t_+1: end);