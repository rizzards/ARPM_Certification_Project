function [par, fval, exitflag, output, sig2, eps] = MLGarch(Dx, S2_0, P0, G, fp)
% parameter bounds: c >= 0, 0 <= a <= 1, 0 <= b <= 1 

% stationarity constraint: a + b <= g < 1

T_ = length(Dx);

f = @(P0)LL(P0, T_, S2_0, Dx, fp);
options = optimoptions(@fmincon,'StepTolerance', 1e-9,'MaxIterations',2000);
[par,fval,exitflag,output] = fmincon(f, P0, [0,1,1], G, [], [], [0,0,0], [Inf,1,1],[], options);

sig2 = zeros(1,T_);
eps  = zeros(1,T_);
    
sig2(1,1) = S2_0;
eps(1,1)  = Dx(1,1)/sqrt(sig2(1,1));

for t = 2: T_
    sig2(1,t) = par(1) + par(3).*sig2(1,t-1) + par(2).*(Dx(1,t-1)^2);
    eps(1,t)  = Dx(1,t)/sqrt(sig2(t));
end
    
   
    