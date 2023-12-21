clear; clc;

invariance_tests_e_k=0;  % if = 1 the script generates analyses performing tests for ellipsoid invariance test
                         % and numerical Kolmogorov-Smirnov test over j_ simulations
invariance_tests_c=1;    % if = 1 the script generates analyses performing tests for copula invariance test
invariance_graphs=0;     % if = 1 the script shows figures of the tests for invariance
hiddenVol = 0;           % if = 1 the script fit paramters to find hidden volatility for a stochastic volatility process
%% QUEST FOR INVARIANCE

load db_Stocks 



Stocks.n_= (5:13); % we consider the following stocks in the dataset

% Risk drivers: log- dividend adjusted prices

if ndims(Stocks.n_) == 1
    Stocks.x = log(AdjustedPrices.prices(1:Stocks.n_,:));
else
    Stocks.x = log(AdjustedPrices.prices(Stocks.n_,:));
end
    
Stocks.x_tnow = Stocks.x(:,end);

% Invariants: compounded returns (risk drivers' increments)
Stocks.epsi=diff(Stocks.x,1,2);

if ndims(Stocks.n_) == 1
    Stocks.i_ = Stocks.n_; %number of selected stocks
else
    Stocks.i_ = length(Stocks.n_);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%invariants with fit VAR1
% dt = 1;
% t_= size(Stocks.epsi,2); %length of the time series of the invariants
% p=ones(1,t_)/t_;  %flat flexible probabilities
% [b, a] = FitVAR1MVOU(Stocks.epsi, Stocks.x(:,1:end-1), dt, p, 100, 0, 0, 0, 0, 'VAR1');
% Stocks.epsi2=Stocks.x(:,2:end)-repmat(a,1,t_)-b*Stocks.x(:,1:end-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Invariance tests
if invariance_tests_e_k==1
    
    % input Kolmogorov-Smirnov test
    simu = 100; % number of simulations
    alpha = 0.01;
    pass_fin = zeros(1,Stocks.i_);
end   


if invariance_tests_c==1
    % input Copula based tests for invariance
    abs_epsi = abs(Stocks.epsi);
    lag_= 6;
    SW_stocks = zeros(lag_,Stocks.i_);
end
 
    
for i=1:Stocks.i_
    
    %  Ellipsoid Invariance Test
    if invariance_tests_e_k==1
        if invariance_graphs == 1
            figure()
            IIDAnalysis(Stocks.epsi(i,:))
        end

     %  Kolmogorov-Smirnov test
        pass = KolmSTest(Stocks.epsi(i,:),simu,alpha);
        pass_fin(1,i) = pass;
        
    end

    %  Copula based tests for invariance (The strongest)
    if invariance_tests_c==1
        if invariance_graphs == 1
            figure()
        end
        [dep_Epsi,~] = CopulaInvariant(abs_epsi(i,:),lag_);
        close
        SW_stocks(:,i) = dep_Epsi;               
    end
end
    
inv_.test = sum(SW_stocks > 0.1);
inv_.id = inv_.test > (lag_*0.4);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Invariance step 2

% Garch(1,1) residuals

[Stocks.epsiGarch, par] = GarchResiduals(Stocks.x(inv_.id,:), length(Stocks.x), [], 0.95 );

id = Stocks.n_(inv_.id)-Stocks.n_(1)+1;
[n,T]= size(Stocks.epsiGarch);

% exponential smoothing prob
lmd=0.0166;
p=exp(-lmd*(T-(1:T)'));
p=p/sum(p);

[u, sig] = FPmeancov(Stocks.x(id,:),p);

if invariance_tests_c==1
    % input Copula based tests for invariance
    abs_epsi = abs(Stocks.epsiGarch);
    abs_epsi = abs_epsi(:,2:end);
end

%  Copula based tests for invariance (The strongest)

if invariance_tests_c==1
    for z=1:sum(inv_.id)
    
        if invariance_graphs == 1
            figure()
        end
        [dep_Epsi,~] = CopulaInvariant(abs_epsi(z,:),lag_);
        close
        SW_stocks(:,id(z)) = dep_Epsi;               
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hidden Vol- Stoch Vol Computation

if hiddenVol == 1
    
    % Fit the stochastic volatility model

    % log(squared returns)
    if ndims(Stocks.n_) == 1
        Stocks.y = log(AdjustedPrices.prices(1:Stocks.n_,:).^2);
    else
        Stocks.y = log(AdjustedPrices.prices(Stocks.n_,:).^2);
    end

    %initial parameters
    phi0=0; phi1=.99; sQ=0.14;  sR0=0.9; mu1=-2; sR1=2;

    hidden_vol = zeros(Stocks.i_,length(Stocks.y(1,:)));

    for i=1:Stocks.i_

        % further initial parameters
        alpha=mean(Stocks.y(i,:));
        initpar=[phi0,phi1,sQ,alpha,sR0,mu1,sR1]; %group together 

        [param,~,~,~]=FitStochasticVolatilityModel(Stocks.y(i,:),initpar); 
        phi=param(1);
        phi1=param(2);   
        sQ=param(3);
        alpha=param(4);
        sR0=param(5);
        mu1=param(6);   
        sR1=param(7);
        [~,log_hiddenvol2] = FilterStochasticVolatility(Stocks.y(i,:),phi0,phi1,sQ,alpha,sR0,mu1,sR1);

        %hidden volatility
        hidden_vol(i,:) = sqrt(exp(log_hiddenvol2));

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Stocks.x_tnow1 = Stocks.x((inv_.test < (lag_*0.4)),end);
Stocks.x_tnow2 = Stocks.x((inv_.test > (lag_*0.4)),end);
Stocks.x_tnowSort = [Stocks.x_tnow1; Stocks.x_tnow2] ;
Stocks.epsiSort = Stocks.epsi((inv_.test < (lag_*0.4)),:);
Stocks.epsiSort = [Stocks.epsiSort ; Stocks.epsiGarch(:,2:end)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Projection

hist_approach = 0;       % if = 1 the historical approach (sampled sequences)
analytic_approach = 1;   % if = 1 analytical approach, normality of risk drivers)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Invariants (all) and Exponential decay probabilities
epsi=[Stocks.epsi];
[i_,t_]=size(epsi);
p=ExponentialDecayProb(t_,250);

% Current value of (all) the risk drivers
x_tnow=[Stocks.x_tnowSort];
d_=size(x_tnow,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% This script projects to an horizon t_hor = t_now + 18 days the risk drivers
tau_proj=18 ; %t_hor = tnow + 18 days 

% Path of the invariants: sampled sequences (bootstrap) approach
if hist_approach == 1
    
    j_=5000; %number of scenarios
    Epsi_path = zeros(Stocks.i_,tau_proj,j_);
    for tau=1:tau_proj
        Epsi_path(:,tau,:)=Bootstrap(epsi,p,j_);
    end
end

% Path of the invariants: sampled sequences (bootstrap) approach
% Estimation: normality assumption
if analytic_approach == 1
    [mu, s2]=MaxLikelihoodFPLocDispT(epsi,100,p,10^-6,1);
    j_=5000; %number of scenarios
    Epsi_path = zeros(Stocks.i_,tau_proj,j_);
    for tau=1:tau_proj
        Epsi_path(:,tau,:)=MultivNormalMomMatch(mu,s2,j_);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Path of the risk drivers
X_path=zeros(d_,tau_proj+1,j_); %initialization
X_path(:,1,:)=repmat(x_tnow,1,1,j_); %first node of the path: current value of the risk drivers

% Project stocks risk drivers according to a multivariate random walk
RandomWalk_idx=[1:Stocks.i_-sum(inv_.id)]; %position of the random walk entries in the risk drivers and invariants panels
for j=1:j_
X_path(RandomWalk_idx,2:end,j)= repmat(X_path(RandomWalk_idx,1,j),1,tau_proj)+cumsum(Epsi_path(RandomWalk_idx,:,j));
end

% Project stocks risk drivers according to a GARCH(1,1) model
Garch_idx=[Stocks.i_-sum(inv_.id)+1, Stocks.i_];
for z = 1 : Garch_idx(2)-Garch_idx(1)+1
    sig2 = sqrt(sig(z,z));
    for q = 3:tau_proj
        for j=1:j_
        
            X_path(z+Garch_idx(1)-1,q,j)= X_path(z+Garch_idx(1)-1,q-1,j)+ u(z) + sig2.*Epsi_path(z+Garch_idx(1)-1,q-1,j); 
            sig2 = par(1) + par(3).* sig(z,z) + par(2).* (X_path(z+Garch_idx(1)-1,q-1,j)-X_path(z+Garch_idx(1)-1,q-2,j)-u(z)).^2 ;
        end    
    end
end


% Probabilities associated to the projected paths
p=ones(1,j_)/j_;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Pricing 

% Stocks. Compute the scenarios of the ex-ante P&L of the stocks via exact pricing starting from the scenarios of the log-values at the horizon
% current values
    Stocks.v_tnow=exp(Stocks.x_tnowSort);   
% values at the horizon
    Stocks.V_thor=exp(squeeze(X_path(1:Stocks.i_,end,:)));
% P&L's
    Stocks.Pi=Stocks.V_thor-repmat(Stocks.v_tnow,1,j_); 