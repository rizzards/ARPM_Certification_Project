function [alpha, beta, s2_U, U] = LassoRegFP(X, Z, p, lambda)
% Weighted lasso regression function 
% Note: LassoRegFP includes the function solveLasso created by GAUTAM V. PENDSE, http://www.gautampendse.com)
% We made some changes in PENDSE's function in order to adapt it with SYMMYS's notation
% The changes are made in conformity with the Creative Commons Attribution 3.0 Unported License
%  INPUTS
% x       :[matrix](n_ x t_) time series of market observations;
% z       :[matrix](k_ x t_) time series of factors
% p       :[vector](t_ x 1) flexible probabilities
% lambda  :[vector](l_ x 1) vectot of penalties
%  OUTPUT
% alpha   :[matrix](n_ x l_) shifting term
% beta    :[array](n_ x k_ x l_) loadings
% s2_U    :[array](n_ x n_ x l_) covariance of residuals
% U       :[array](n_ x t_ x l_) residuals

%% Code

[n_,t_] = size(X);
k_ = size(Z,1);
l_ = length(lambda);

% if p are not provided, observations are equally weighted
if nargin<4 || isempty(p); p = (1/t_)*ones(1,t_); end


% solve optimization
if l_==1 && lambda==0
    [alpha, beta, s2_U, U] = OrdLeastSquareFPNReg(X,Z,p);
else
    % preliminary de-meaning of x and z
    [m_X,~] = FPmeancov(X,p);
    [m_Z,~] = FPmeancov(Z,p);
    X_c = X - repmat(m_X,1,t_);
    Z_c = Z - repmat(m_Z,1,t_);
    % trick to adapt function solveLasso to the FP framework
    X_p = X_c*sqrt(diag(p));
    Z_p = Z_c*sqrt(diag(p));
    % initialize variables
        beta = nan(n_,k_,l_);
        alpha = nan(n_,l_);
        s2_U = nan(n_,n_,l_);
        U = nan(n_,t_,l_);
    % solve lasso
    for l=1:l_
        for n=1:n_
            output = solveLasso(X_p(n,:),Z_p,lambda(l));
            beta(n,:,l) = output.beta;
        end
        alpha(:,l) = m_X - beta(:,:,l)*m_Z;
        U(:,:,l) = X - repmat(alpha(:,l),1,t_) - beta(:,:,l)*Z;
        [~,s2_U(:,:,l_)] = FPmeancov(U(:,:,l), p);
    end
end

end

function output = solveLasso(x, z, lambda)
%==========================================================================
%               AUTHOR: GAUTAM V. PENDSE                                  
%               DATE: 11 March 2011                                   
%==========================================================================
%               
%               PURPOSE:
%
%   Algorithm for solving the Lasso problem:
%
%           0.5 * (x - beta*z)*(x - beta*z)' + lambda * ||beta||_1
%                                               
%   where ||beta||_1 is the L_1 norm i.e., ||beta||_1 = sum(abs(beta))
%
%   We use the method proposed by Fu et. al based on single co-ordinate
%   descent. For more details see GP's notes or the following paper:
%
%   Penalized Regressions: The Bridge Versus the Lasso
%   Wenjiang J. FU, Journal of Computational and Graphical Statistics, 
%   Volume 7, Number 3, Pages 397?416, 1998
%   
%==========================================================================
%               
%               INPUTS:
%
%       =>      x = 1 by t_ response vector
%
%       =>      z = k_ by t_ design matrix
%
%       => lambda = regularization parameter for L1 penalty
%
%==========================================================================
%               
%               OUTPUTS:
%
%       => output.z = supplied design matrix
%
%       => output.x = supplied response vector
%
%       => output.lambda = supplied regularization parameter for L1 penalty
%
%       => output.beta = computed L1 regularized solution
%
%==========================================================================
%   
%       Copyright 2011 : Gautam V. Pendse
%
%               E-mail : gautam.pendse@gmail.com
%
%                  URL : http://www.gautampendse.com
%
%==========================================================================

%==========================================================================
%               check input args

    if (nargin ~= 3)
        disp('Usage: output = solveLasso(x, z, lambda)');
        output = [];
        return;
    end
    
    % check size of x
    [n_, t_] = size(x);
    
    % is x a row vector?
    if (n_ ~= 1)
        disp('x must be a 1 by t_ vector!!');
        output = [];
        return;
    end
    
    % check size of z
    [k_,t_1] = size(z);
    
    % does z have the same number of columns as x?
    if (t_1 ~= t_)
        disp('z must have the same number of rows as x!!!');
        output = [];
        return;
    end
    
    % make sure lambda > 0
    if (lambda < 0)
        disp('lambda must be >= 0!');
        output = [];
        return;
    end
    
%==========================================================================
%               initialize the Lasso solution

   % This assumes that the penalty is lambda * beta'*beta instead of lambda * ||beta||_1
   beta = (x*z') / (z*z' + 2*lambda);

%==========================================================================
%               start while loop

    % convergence flag
    found = 0;
    
    % convergence tolerance
    tol = 1e-6;
    
    while( found == 0 )
    
        % save current beta
        beta_old = beta;
        
        % optimize elements of beta one by one
        for k = 1:k_
            
            % optimize element i of beta
            
            % get ith col of z
            z_k = z(k,:);
            
            % get residual excluding ith col
            x_k = (x - beta*z) + beta(k)*z_k;           

            % calulate zi*xi' and see where it falls
            deltai = z_k*x_k'; % 1 by 1 scalar   
            if ( deltai < -lambda )        
                beta(k) = ( deltai + lambda )/(z_k*z_k');       
            elseif ( deltai > lambda )          
                beta(k) = ( deltai - lambda )/(z_k*z_k');
            else
                beta(k) = 0;
            end
        end
        
        % check difference between beta and beta_old
        if ( max(abs(beta - beta_old)) <= tol )
            found = 1;
        end
                      
    end
    
%==========================================================================
%   save outputs
    
    output.z = z;
    output.x = x;
    output.lambda = lambda;  
    output.beta = beta;
end
