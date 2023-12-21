function [dep_Epsi,U] = CopulaInvariant(X,lag_)
%  INPUTS
%   X     : [matrix]  (2 x j_) joint scenarios

%% Code
% Estimate the SW measures of dependence over lag_
dep_Epsi = NaN(lag_, 1);
for l = 1 : lag_
    
    
    [~, j_] = size(X); % number of elements, number of scenarios
    probs = ones(1, j_ - l) / (j_ - l); % flat prob
    X = [X(1,1 + l : end); X(1,1 : end - l)];   
            
    % SW measures
    [~, ~,U] = CMASep(X, probs); % grades scenarios

    
    g = NaN(j_, j_);
    for i = 1 : j_
        for k = 1 : j_
            g(i, k) = abs(sum(probs .* (U(1, :) <= i/j_) .* (U(2, :) <= k/j_)) - (i * k)/j_^2);
        end
    end

    dep_Epsi(l) = 12 * mean(g(:));
end

U = U';


%% print 

j_ = size(U,1);

NumBins3d=round(sqrt(j_)/4);


% 'Position',[0.13 0.11 0.354945799457995 0.815]        
subplot(1,20,1:11)
e = hist3(U(:,[1,2]),[NumBins3d NumBins3d]);
hist3(U(:,[1,2]),[NumBins3d NumBins3d],'FaceColor',[0.87058824300766 0.921568632125854 0.980392158031464]);
Hsurface = get(gca,'children');
term = Hsurface.ZData/ mean(mean(e));
ydata1 = linspace(0,1,30);
xdata1 = linspace(0,1,30);
set(Hsurface,'ZData',term);
zticks([0 0.25 0.5 0.75 1]);
hold on
surf(xdata1,ydata1,ones(30,30),'FaceColor',[0.63 0.078 0.18])
xlabel('Grade Obs.'); ylabel('Grade Obs. with lag');
title('Copula Invariance Test')

% 'Position',[0.655067750677507 0.11 0.249932249322493 0.815]
subplot(1,20,15:19)
bar(dep_Epsi,0.45)
ylim([0 1]); title('SW measure of dependence');
xlabel('Lag'); ylabel('Dependence')


% ecd1 = ecdf(X(1,:));    inv1 = ecd1(1:end-1,1);
% ecd2 = ecdf(X(2,:));    inv2 = ecd2(1:end-1,1);
% inv = [inv1, inv2];
% 
% [NumObs,K]=size(inv);
% Copula = zeros(NumObs,K);
%         [~,C]=sort(inv);
%         for k=1:K
%             x=C(:,k);
%             y=[1:NumObs];
%             xi=[1:NumObs];
%             yi = interp1(x,y,xi);
%             Copula(:,k)=yi/(NumObs+1);
%         end 
% joint scenarios of a desired FP-distribution with given FP-copula 
% and FP-marginals  implementing the CMA

% Compute the grid of significant evaluation nodes and cdf grid
% eta = 0.06;
% k_ = j_-lag_;

% y = NaN(n_, k_);
% u = NaN(n_, k_);
% for n = 1 : n_
%     a = interp1(u_(n, :), y_(n, :), eta, 'linear', 'extrap'); % lower quantile
%     b = interp1(u_(n, :), y_(n, :), 1 - eta, 'linear', 'extrap'); % upper quantile
%     y(n, :) = linspace(a, b, k_); 
%     u(n, :) = interp1(y_(n, :), u_(n, :), y(n, :), 'linear', 'extrap');
% end

% Compute the joint scenarios through the CMA (combination) routine
% Copula = CMAComb(y_, u_, Cop');
