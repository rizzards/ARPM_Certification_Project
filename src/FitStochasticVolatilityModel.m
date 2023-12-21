function [P, fval, exitflag, output] = FitStochasticVolatilityModel(y, initpar)


%     Maximum likelihood estimation of the SV model parameters
%  
%     Args:
%     
%         y          (t_ x 1)   : time series.
%         initpar    (7 x 1)    : initial guess for parameters.
%         
%     Returns:
%     
%         P          (7 x 1)    : optimal parameters.
%         fval                  : maximum likelihood achieved.
%         exitflag   str        : exit status.
%         output     obj        : optimization exit object. 

%% Code
f = @(p)Linn(p,y);
options = optimoptions(@fminunc,'StepTolerance', 1e-8, 'Algorithm', 'quasi-newton', 'MaxIterations',3000);
[P, fval, exitflag, output] = fminunc(f, initpar, options);


