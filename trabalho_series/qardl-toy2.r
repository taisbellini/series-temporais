# This code simulates a sample path of size T from a simple QARDL(1,1) model:
# Q(𝜏| ℱ[t-1]) = α₀(𝜏) + α₁(𝜏)*Y[t-1] + θ₁(𝜏)*Z[t-1]

# set.seed(1)
library(quantreg)
library(conquer)
library(forecast)
library(qrcmNP)
source('utils.R')

tau.grid = seq(from=.01,to=.99, by=.02)

# Value of the conditional quantile function at three
# vertices of the unit square (v10, v01, v00 and v11)
# Since v11 = v10 + v01 - v00 we can choose, say, v10 and v01
# freely and then choose v00 restricted with the restriction that
# v10 + v01 - v00 is non-decreasing and has its range contained
# in the unit interval [0,1]. If the functions are chosen 
# to be differentiable quantile functions, it is sufficient 
# then to choose v00 so that
# v10' + v01' \ge v00'
# or, equivalently, log(v10' + v01')-log(v00')\ge0
# The function qdbplot below allows us to investigate
# if this inequality holds when v00, v10 and v01 
# are chosen to be Beta quantile functions: 
qdbplot = function(a10,b10,a01,b01,a00,b00){
tau = seq(from=0.00001,to=.99999,by=.00001)
v00prime = 1/dbeta(qbeta(tau,a00,b00),a00,b00)
v10prime = 1/dbeta(qbeta(tau,a10,b10),a10,b10)
v01prime = 1/dbeta(qbeta(tau,a01,b01),a01,b01)
plot(tau, log(v10prime+v01prime)-log(v00prime), type = 'l'); abline(h=0, col='red')
}


# A function for visualizing the Beta densities:
dbplot = function(a,b){
 tau = seq(from=0,to=1,by=.01)
 plot(tau, dbeta(tau,a,b))
}

# Parameters for the Beta quantile functions v10, v01 and v00
{a10 = 1; b10 = 3; a01 = 5; b01 = 1; a00 = 3; b00 = 1}
qdbplot(a10,b10, a01,b01, a00,b00)

# Conditional quantile function of Y[t] given Y[t-1] = 1 and X[t-1]=0
v10 = function (tau) qbeta(tau, a10,b10)
# Conditional density of Y[t] given Y[t-1] = 1 and X[t-1]=0
dbplot(a10,b10)

# Conditional quantile function of Y[t] given Y[t-1] = 0 and X[t-1]=1
v01 = function(tau) qbeta(tau,a01,b01)
# Conditional density of Y[t] given Y[t-1] = 0 and X[t-1]=1
dbplot(a01,b01)

# Conditional quantile function of Y[t] given Y[t-1] = 0 and X[t-1]=0
v00 = function(tau) qbeta(tau,a00,b00)
# Conditional density of Y[t] given Y[t-1] = 0 and X[t-1]=0
dbplot(a00,b00)

# Conditional quantile function of Y[t] given Y[t-1] = 1 and X[t-1]=1
v11 = function(tau) v10(tau) + v01(tau) - v00(tau)
# Plot of v11 (must be nondecreasing and with range contained in [0,1])
plot(tau.grid, v11(tau.grid))
# Conditional density of Y[t] given Y[t-1] = 1 and X[t-1]=1
q11 = function(tau) {
1/dbeta(v10(tau),a10,b10) + 1/dbeta(v01(tau),a01,b01) - 1/dbeta(v00(tau),a00,b00)
}
plot(v11(tau.grid), 1/q11(tau.grid), xlim=c(0,1))

# Functional parameters for the quantile regression equation
alpha0 = v00
alpha1 = function(tau) v10(tau) - v00(tau)
theta1 = function(tau) v01(tau) - v00(tau)

# Quantile function of Y given ℱ[t-1]
Q = function(tau,Y.current,Z.current){
 alpha0(tau) + alpha1(tau)*Y.current + theta1(tau)*Z.current
}

# Simulating the sample paths:

# Arbitrary starting point (Y0,Z0)
Y0 = runif(1)
Z0 = runif(1)

Y = Z = numeric()
Z.current = Z0
Y.current = Y0
T = 10001
for (t in 1:T){
# Simulates Y[t] given ℱ[t-1] using the Fundamental Theorem of Simulation
 Y[t] = Q(runif(1), Y.current, Z.current)
 
 Z[t] = Q(runif(1), Z.current, Y.current)
 # Q(runif(1), Z.current, Y.current)
 # Z[t] = runif(1,max(Y[t]-.01,0),min(Y[t]+.01,1))
 
 Z.current = Z[t]
 Y.current = Y[t]
}

Y = Y[-1]
Z = Z[-1]
T = length(Y)
acf(Y, lwd=16, lend=3, col='gray')
pacf(Y, lwd=16, lend=3, col='gray')

plot(Y[2:T]~Y[1:(T-1)], pch=16, col=rgb(0,0,0,.4))
plot(Y[2:T]~Z[1:(T-1)], pch=16, col=rgb(0,0,0,.4))
hist(Y, border=NA, breaks="FD")
ts.plot(Y)

# True and estimated functional parameters

# quantile regression
Yvec = Y[2:T]
Xmat = cbind(Y[1:(T-1)], Z[1:(T-1)])
qrfit = rq(Yvec~Xmat, tau = tau.grid)
plot(tau.grid, alpha0(tau.grid), type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, coef(qrfit)[1,], lwd=2, col = rgb(0,0,0,.7))

plot(tau.grid, alpha1(tau.grid), type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, coef(qrfit)[2,], lwd=2, col = rgb(0,0,0,.7))

plot(tau.grid, theta1(tau.grid), type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, coef(qrfit)[3,], lwd=2, col = rgb(0,0,0,.7))

# kernel estimation
conquerfit = sapply(tau.grid, function(tau) conquer(Xmat,Yvec,tau=tau)$coeff)
plot(tau.grid, alpha0(tau.grid), type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, conquerfit[1,], lwd=2, col = rgb(0,0,0,.7))
# 
plot(tau.grid, alpha1(tau.grid), type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, conquerfit[2,], lwd=2, col = rgb(0,0,0,.7))
# 
plot(tau.grid, theta1(tau.grid), type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, conquerfit[3,], lwd=2, col = rgb(0,0,0,.7))


# global coefficients proposed
phi = phi_generator(3, tau.grid)
Xmat1 = cbind(rep(1,nrow(Xmat)), Xmat)
globalfit = global_qr(taus = tau.grid, phi = phi, X = Xmat1, y = Yvec, lambda=10)

plot(tau.grid, alpha0(tau.grid), type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, globalfit$bhat[1,], lwd=2, col = rgb(0,0,0,.7))

plot(tau.grid, alpha1(tau.grid), type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, globalfit$bhat[2,], lwd=2, col = rgb(0,0,0,.7))

plot(tau.grid, theta1(tau.grid), type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, globalfit$bhat[3,], lwd=2, col = rgb(0,0,0,.7))

# global coefficients Sottile
# Estimate qadl time series coefficients using qrcm
k.user = 3
p = length(tau.grid)
fo3o = piqr(Yvec~Xmat, formula.p = ~slp(p,k=k.user), lambda = 10)
fo3o$coefficients$lambda1
fo4o=slp(tau.grid,k=k.user)
PHI = cbind(1,fo4o)
BETA = fo3o$coef$lambda1%*%t(PHI)
Ypred = cbind(1,Xmat)%*%BETA

# True and estimated coefficients
plot(tau.grid, alpha0(tau.grid), type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, BETA[1,], lwd=2, col = rgb(0,0,0,.7))

plot(tau.grid, alpha1(tau.grid), type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, BETA[2,], lwd=2, col = rgb(0,0,0,.7))

plot(tau.grid, theta1(tau.grid), type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, BETA[3,], lwd=2, col = rgb(0,0,0,.7))


# ssr comparison
ssrqr = cbind((alpha0(tau.grid) - coef(qrfit)[1,])^2, (alpha1(tau.grid) - coef(qrfit)[2,])^2, (theta1(tau.grid) - coef(qrfit)[3,])^2)
ssrconquer = cbind((alpha0(tau.grid) - conquerfit[1,])^2, (alpha1(tau.grid) - conquerfit[2,])^2, (theta1(tau.grid) - conquerfit[3,])^2)
ssrglobal = cbind((alpha0(tau.grid) - globalfit$bhat[1,])^2, (alpha1(tau.grid) - globalfit$bhat[2,])^2, (theta1(tau.grid) - globalfit$bhat[3,])^2)
ssrpiqr = cbind((alpha0(tau.grid) - BETA[1,])^2, (alpha1(tau.grid) - BETA[2,])^2, (theta1(tau.grid) - BETA[3,])^2)

plot(tau.grid, ssrqr[,1], type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, ssrconquer[,1], lwd=2, col = rgb(0,1,0,.7))
lines(tau.grid, ssrglobal[,1], lwd=2, col = rgb(0,0,0,.7))

plot(tau.grid, ssrconquer[,2], type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, ssrglobal[,2], lwd=2, col = rgb(0,0,0,.7))

plot(tau.grid, ssrconquer[,3], type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, ssrglobal[,3], lwd=2, col = rgb(0,0,0,.7))

ssrqr_mean = apply(ssrqr, 1, mean)
ssrconquer_mean = apply(ssrconquer, 1, mean)
ssrglobal_mean = apply(ssrglobal, 1, mean)
plot(tau.grid, ssrconquer_mean, type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, ssrglobal_mean, lwd=2, col = rgb(0,0,0,.7))
lines(tau.grid, ssrqr_mean, lwd=2, col = rgb(0,1,0,.7))

cbind(mean(ssrqr),mean(ssrconquer), mean(ssrglobal), mean(ssrpiqr))
