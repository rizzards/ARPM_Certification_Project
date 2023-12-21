function [ExpectedValue,Std_Deviation, StdHoldings] = EfficientFrontierQPReturn(NumPortf, Covariance, ExpectedValues,Current_Prices,Budget)
% this function computes the mean-variance efficient frontier for performance = return by quadratic programming
% (long only constraints)

warning off;
options=optimset('display','off');

NumAssets=size(Covariance,2);

% determine exp value of minimum-variance portfolio
FirstDegree=zeros(NumAssets,1);
SecondDegree=Covariance;
Aeq=[];
beq=[];
A=[-eye(NumAssets);Current_Prices'];
b=[zeros(NumAssets,1);1];
htilde_0=(Budget/NumAssets)./(Current_Prices.^2).*ones(NumAssets,1); %equally weighted portfolio
MinVol_htilde = quadprog(SecondDegree,FirstDegree,A,b,Aeq,beq,[],[],htilde_0,options);
MinVol_ExpVal=MinVol_htilde'*ExpectedValues;

% determine exp value of maximum-expected value portfolio

Ret=ExpectedValues./Current_Prices;
Max_ExpVal=ExpectedValues(Ret==max(Ret))/Current_Prices(Ret==max(Ret)); 

% slice efficient frontier in NumPortf equally thick horizontal sectors in the upper branch only
Target_ExpectedValues=linspace(MinVol_ExpVal, Max_ExpVal, NumPortf);

% compute the NumPortf compositions and risk-return coordinates
StdHoldings=[];
Std_Deviation=[];
ExpectedValue=[];

for i=1:NumPortf
    
    % determine least risky portfolio for given expected return
    AEq=[Aeq
        ExpectedValues'];
    bEq=[beq
        Target_ExpectedValues(i)];
    htilde = quadprog(SecondDegree,FirstDegree,A,b,AEq,bEq,[],[],htilde_0,options)';
    StdHoldings=[StdHoldings 
                    htilde];
    Std_Deviation=[Std_Deviation
                sqrt(htilde*Covariance*htilde')];
    ExpectedValue=[ExpectedValue
                    Target_ExpectedValues(i)];
end
