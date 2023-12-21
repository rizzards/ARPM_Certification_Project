function like = Linn(p,y)


phi0  = p(1);
phi1  = p(2);
sQ    = p(3);
alpha = p(4);
sR0   = p(5);
mu1   = p(6);
sR1   = p(7);

Y = y;

like = FilterStochasticVolatility(Y, phi0, phi1, sQ, alpha, sR0, mu1, sR1);