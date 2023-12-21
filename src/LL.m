function L = LL(p, T_, S2_0, Dx, fp)
    
t_   = T_;
s2_0 = S2_0;
dx   = Dx;
FP   = fp;

c = p(1);
a = p(2);
b = p(3);

s2 = zeros(1,t_);
B  = zeros(1,t_);

s2(1,1) = max(c + (a+b).*s2_0, 5e-6);
B(1,1)  = (dx(1,1)^2)/s2(1,1);

for t = 2: t_

    s2(1,t) = max(c + b.*s2(1,t-1) + a.*(dx(1,t-1)^2), 5e-6);
    B(1,t)  = (dx(1,t)^2)/s2(1,t);
    
end

% this is -L        
L = sum(FP.*log(s2)) + sum(FP.*B);