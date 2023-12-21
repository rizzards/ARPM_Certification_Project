function [par, sig2, epsi, lik] = FitGARCHFP(dx, s2_0, p0, g , FP)

%     Fit GARCH(1,1) model (WARNING: doesn't check for stationarity)
%  
%     Args:
%     
%         dx          (1 x t_) or (t_ x 1) : realized increments of the GARCH process x.
%         s2_0                 numpy.int64: initial value for sigma**2.
%         par0        dict: initial guess of parameters:
%                              par0[0]: c,
%                              par0[1]: a,
%                              par0[1]: b.
%         g           (len(g) x 1) : 
%                              stationarity constraint: a + b < g (default: g = 0.95).        
%         FP          (1 x t_) or (t_ x 1) : (optional) flexible probabilities.
%                              Default: eqully weighted observations.
%         
%     Returns:
%     
%         par         (3 x len(g))  : for each constraint g[i], par[:,g[i]] gives the estimates of the three parameters c,a,b.                             
%         sig2        (len(g) x t_) : estimated path of the squared scatter.
%         epsi        (len(g) x t_) : residuals.
%         lik         (len(g) x t_) : maximum likelihood achieved
%                              in correspondence of each value of g.

% Initializations
    
% makes dx a (1 x t_) array
if size(dx,2) == 1
    dx = dx';
end

t_ = size(dx,2);

% if FP are not provided, observations are equally weighted
if isempty(FP) 
    FP = ones(1,t_)*1.0/t_;
end

% makes FP a (1 x t_) array
if size(FP,2) == 1
    FP = FP';
end

% if g not provided
if isempty(g) 
    g = 0.95;
end

%reshape g (round to integer)
g = round(g);

par = zeros(3,g);

sig2 = zeros(g, t_);       
epsi = zeros(g, t_);       
lik  = zeros(1,g);

% no cycle for i = 1: g
[p, fval, ~, ~, Sig2, Epsi] = MLGarch(dx, s2_0, p0, g, FP);
par(:,1)  = p;
sig2(1,:) = Sig2;
epsi(1,:) = Epsi;
lik(1)    = -fval;





