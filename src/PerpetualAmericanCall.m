function [CallPrice]=PerpetualAmericanCall(x,varargin)
% Price of a perpetual American call option with Bachelier underlying
% INPUTS
% x: underlying [vector]
% varargin are either
%   mu [scalar < 0]: drift parameter for the underlying 
%   sigma [scalar]: volatility parameter for the underlying
%   k [scalar]: strike of the option (default k=0) 
% or
%   eta [scalar]: inverse-call transformation parameter (=-sigma^2/(2mu)). Then k=0 by default.
% OUTPUT
% CallPrice [vector]

if length(varargin)==3;
    mu=varargin{1};
    sigma=varargin{2};
    k=varargin{3};
    
    gamma=(-2*mu)/(sigma^2);
    eta=1/gamma; 
    
elseif length(varargin)==2;
    mu=varargin{1};
    sigma=varargin{2};
    k=0;
    
    gamma=(-2*mu)/(sigma^2);
    eta=1/gamma; 

elseif length(varargin)==1;
    k=0;
    eta=varargin{1};
    gamma=1/eta;
end


CallPrice=zeros(size(x));

boundary=k+eta; 

CallPrice(x>=boundary)=x(x>=boundary)-k;
CallPrice(x<boundary)=eta*exp(-gamma*boundary)*exp(gamma*x(x<boundary));

end