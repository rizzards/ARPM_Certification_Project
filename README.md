# ARPM Code Description

The main code "S_Main", which runs the entire project is divided in three main parts : quest for invariance, projection, pricing.

## 1) Quest for invariance

The first step in the quest for invariance is to retrieve the invariants from a group of stocks on our
choice taken from an available database 'db_Stocks' . The retrieved risk drivers are the log-dividend
ajusted prices. Then, under the assumption that the risk drivers follow a random walk I retrieve the
invariants. My quest for invariance is organized upon three test: the ellipsoid invariance test, a
numerical value based on a certain number of simulation of the Kolmogorov Smirnov test, and the
copula invariance test. The most strong one is the copula test.
By the means of the latter I organized a second step quest for invariance based on the value
obtained by the Schweizer and Wolff measure of dependence. I used this measure to filter the
acceptable invariants. All the risk drivers that had more than 0.1 of the SW measure of dependence
in more than the half of the lagged space value were checked again through a new quest.
The new assumption is that the risk drivers of these stocks follow a GARCH(1,1) process. The new
quest over the residuals of a GARCH process reveals better values for the SW measure.
The code offers, if requested, to fit the risk drivers of the second quest to a Stochastic Volatility
proccess and provides the values of the hidden volatility.
By setting some variables to 1 at the beginning of the code you have the possibility to choose if
skipping the quest for invariance or just to not plot the result graphs.

## 2) Value projection
The second step is about the projection of the invariants to the investment horizon, which is equal in this case to 18 business days. The code gives the opportunity to choose between the historical approach and the analytical approach in the projection of the invariants. On the basis of their initial assumption, consequently the invariants are projected with their respective stochastic process (random walk or GARCH(1,1)). 

## 3) Pricing
The last step is about the pricing. In this part are reported the scenarios of the ex-ante P&L of the stocks via exact pricing after having obtained the values of the risk drivers at the specified horizon.