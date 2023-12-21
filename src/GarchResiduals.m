function [epsi, par] = GarchResiduals(x, t_garch, p_garch, g )
%  Computes the residuals of a Garch(1,1) fit on x.
%  
%     Args:
%     
%         x           (n_ x t_) or (t_ x 1) : time-series of observations.
%         t_garch                           : number of observations processed at every iteration.
%         p_garch     (1 x t_) or (t_ x 1)  : (optional) flexible probabilities.
%                                           : exponential decay with half-life of 6 months.
%         g                                 : stationarity constraint: a + b < g < 1 (default: g = 0.95).        
%         
%     Returns:
%     
%         epsi        (n_ x t_)  or (t_ x 1): residuals.
%         par          par[1]: c, par[2]: a, par[2]: b 

%% Code

% Initializations
% makes x a (n_ x t_) array also if n_ == 1

if size(x,2) == 1
    x = x';
end

% if FP are not provided, observations are exponentially weighted
if isempty(p_garch) 
    Lambda  = log(2.0)/180.0;
    p_garch = exp(-Lambda*(t_garch-1:-1:0));
    p_garch = p_garch/sum(p_garch);
    if size(p_garch,2) == 1
        p_garch = p_garch';
    end
end

% if g not provided
if isempty(g) 
    g = 0.95;
end


[n_, t_obs] = size(x);
p0 = [1.0, 0.45*g, 0.45*g];
lam = 0.7;

if t_garch == t_obs
    
    epsi = zeros(n_, t_obs);
    
    for n = 1: n_
        
        s2_0 = lam.*var(x(n,:)) + (1-lam).*sum(lam.^(0:1:t_obs-1).*(x(n,:).^2));
        [par, ~, epsi(n,:), ~] = FitGARCHFP(x(n,:), s2_0, p0, g, p_garch);
    end
    
else  % use rolling window
    
    t_   = round(t_obs - t_garch);        
    epsi = zeros(n_, t_);
    
    for t = 1:t_
        
        for n = 1:n_
            
            x_t  = x(n, t : t + t_garch - 1);
            s2_0 = lam.*var(x_t) + (1-lam).*sum(lam.^(0:1:t_garch-1).*(x_t.^2));
            
            [par, ~, e, ~] = FitGARCHFP(x_t, s2_0, p0, g, p_garch);
            
            epsi(n,t) = e(1,end-1);
        end
        
    end
end

if n_ == 1
    epsi = epsi';
end
            
    
