function pass_fin = KolmSTest(X,j,alpha)
% The code computes the percentage of the average over j simulations that 
% 2 random distribution are not equally distributed on alpha level of confidence. 
% kstest2 : the null hypothesis that the data in vectors x1 and x2 are from 
% the same continuous distribution, using the two-sample Kolmogorov-Smirnov test. 
% The alternative hypothesis is that x1 and x2 are from different continuous distributions. 
% The result h is 1 if the test rejects the null hypothesis at the alpha significance level,
% and 0 otherwise.


pass = zeros(1,j);


for op = 1:j
    [ks1_1, ks2_1] = RandomSplit(X);
    [pass(1,op),~] = kstest2(ks1_1, ks2_1,'Alpha',alpha);
end
pass_fin = sum(pass(1,:))/j;

