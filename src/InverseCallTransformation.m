function x= InverseCallTransformation(y, varargin)
% This function performes the Inverse Call Transformation, i.e. the inverse
% of the price profile of a Perpetual American Call option with Bachelier
% underlying (see A. Meucci, A. Loregian - "Neither "Normal" not "Lognormal": Modeling
% Interest Rates Across all Regimes" to appear (2013))
%
%INPUTS
% y [matrix]: (n_ x t_) time series to be transformed (e.g. yields)
% varargin are either
%   mu [scalar < 0]: drift parameter for the perpetual American call Bachelier underlying
%   sigma [scalar]: volatility parameter for the perpetual American call Bachelier underlying
%   k [scalar]: strike of the perpetual American call (default k=0) 
% or
%   eta [scalar]: inverse call transformation parameter (=-sigma^2/(2mu)). Then k=0 by default.
%OUTPUTS
% x [matrix]: (n_ x t_) output of the inverse-call transformation (e.g. shadow rates)

if length(varargin)==3;
    mu=varargin{1};
    sigma=varargin{2};
    k=varargin{3};
    
    gamma=(-2*mu)/(sigma^2);
    eta=1/gamma; %inverse call parameter
    
elseif length(varargin)==2;
    mu=varargin{1};
    sigma=varargin{2};
    k=0;
    
    gamma=(-2*mu)/(sigma^2);
    eta=1/gamma; %inverse call parameter

elseif length(varargin)==1;
    k=0;
    eta=varargin{1};
    gamma=1/eta;
    
end
       

x=zeros(size(y)); %storage

 boundary=k+eta; 
 const=eta*exp(-gamma*boundary); 

 x(y<eta)=log(y(y<eta)/const)/gamma;
 x(y>=eta)=y(y>=eta)+k; 
end

