function [likelihood, xp] = FilterStochasticVolatility(y, phi0, phi1, sQ, alpha, sR0, mu1, sR1)



% Filters stochastic volatility.
%  
%     Args:
%     
%         y                (t_ x 1) : time series of log(returns^2).
%         phi0,phi1,sQ,alpha,sR0,mu1,sR1      parameters of the stochastic volatility model.
%         
%         
%     Returns:
%     
%         likelihood       : -log(likelihood).
%         xp                (t_ x 1) : log of the squared-hidden volatility.


%% Code
% Initializations
t_ = length(y);
Q = sQ^2;

R0 = sR0^2;
R1 = sR1^2;

xf = 0.0;                              % <-- h_0**0   

Pf = max((sQ^2)/(1.0 - phi1), 0.0);    % <-- P0**0 

xp = zeros(1,t_);
Pp = zeros(1,t_);

pi1 = 0.5;
pi0 = 0.5;

fpi1 = 0.5;
fpi0 = 0.5;

likelihood = 0.0;           % -log(likelihood)

% Filtering

for i = 1 : t_
       
    xp(1,i) = phi1*xf + phi0;
    Pp(1,i) = phi1*Pf*phi1 + Q;

    sig1  = Pp(1,i) + R1;
    sig0  = Pp(1,i) + R0;

    k1 = Pp(1,i)/sig1;
    k0 = Pp(1,i)/sig0;

    e1 = y(1,i) - xp(1,i) - mu1 - alpha;
    e0 = y(1,i) - xp(1,i) - alpha;

    den1  = (1.0/sqrt(sig1))*exp(-0.5*(e1^2)/sig1);
    den0  = (1.0/sqrt(sig0))*exp(-0.5*(e0^2)/sig0);
    denom = pi1*den1 + pi0*den0;

    fpi1 = pi1*den1/denom;
    fpi0 = pi0*den0/denom;

    xf = xp(1,i) + fpi1*k1*e1 + fpi0*k0*e0;
    Pf = fpi1*(1.0 - k1)*Pp(1,i) + fpi0*(1.0 - k0)*Pp(1,i);

    likelihood = likelihood - 0.5*log(denom); 
    
        
end