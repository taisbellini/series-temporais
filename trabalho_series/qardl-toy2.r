# This code simulates a sample path of size T from a simple QARDL(1,1) model:
# Q(𝜏| ℱ[t-1]) = α₀(𝜏) + α₁(𝜏)*Y[t-1] + θ₁(𝜏)*Z[t-1]

# set.seed(1)
library(quantreg)
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

adf.test(Y)

# True and estimated functional parameters

qardl_sim = simulate_qardl(Q, Q)
Y = qardl_sim$Y
Z = qardl_sim$Z

T = 10001

# quantile regression
Yvec = Y[2:T]
Xmat = cbind(Y[1:(T-1)], Z[1:(T-1)])

Yvec_lags = Y[3:T]
Xmat_lags = cbind(Y[1:(T-1)], Z[1:(T-1)], c(NA, Y[1:(T-2)]), c(NA, Z[1:(T-2)]))[-1,]
colnames(Xmat_lags) = c("Yt-1", "Xt-1", "Yt-2", "Xt-2")

M = length(tau.grid)

real_coefs = cbind(alpha0(tau.grid), alpha1(tau.grid), theta1(tau.grid), rep(0,M), rep(0,M), rep(1, M))
colnames(real_coefs) = c("Intercept", "Yt-1", "Xt-1", "Yt-2", "Xt-2", "Method")

qrfit = rq(Yvec_lags~Xmat_lags, tau = tau.grid)
qrfit_coefs = cbind(t(coef(qrfit)), rep(2,M))
colnames(qrfit_coefs) = colnames(real_coefs)

plot(tau.grid, real_coefs[,1], type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, coef(qrfit)[1,], lwd=2, col = rgb(0,0,0,.7))

plot(tau.grid, real_coefs[,2], type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, coef(qrfit)[2,], lwd=2, col = rgb(0,0,0,.7))

plot(tau.grid, real_coefs[,3], type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, coef(qrfit)[3,], lwd=2, col = rgb(0,0,0,.7))

plot(tau.grid, real_coefs[,4], type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, coef(qrfit)[4,], lwd=2, col = rgb(0,0,0,.7))

plot(tau.grid, real_coefs[,5], type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, coef(qrfit)[5,], lwd=2, col = rgb(0,0,0,.7))

# kernel estimation
conquerfit = sapply(tau.grid, function(tau) conquer(Xmat_lags,Yvec_lags,tau=tau)$coeff)
plot(tau.grid, alpha0(tau.grid), type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, conquerfit[1,], lwd=2, col = rgb(0,0,0,.7))
# 
plot(tau.grid, alpha1(tau.grid), type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, conquerfit[2,], lwd=2, col = rgb(0,0,0,.7))
# 
plot(tau.grid, theta1(tau.grid), type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, conquerfit[3,], lwd=2, col = rgb(0,0,0,.7))

mean(conquerfit[4,])
mean(conquerfit[5,])

# global coefficients Sottile
# Estimate qadl time series coefficients using qrcm
k.user = 7
p = length(tau.grid)
fo3o = piqr(Yvec_lags~Xmat_lags, formula.p = ~slp(p,k=k.user), lambda = 10)
fo4o=slp(tau.grid,k=k.user)
PHI = cbind(1,fo4o)
BETA = fo3o$coef$lambda1%*%t(PHI)
piqr_coefs = cbind(t(BETA), rep(3,M))
colnames(piqr_coefs) = colnames(real_coefs)
Ypred = cbind(1,Xmat)%*%BETA

plot(tau.grid, alpha0(tau.grid), type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, BETA[1,], lwd=2, col = rgb(0,0,0,.7))

plot(tau.grid, alpha1(tau.grid), type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, BETA[2,], lwd=2, col = rgb(0,0,0,.7))

plot(tau.grid, theta1(tau.grid), type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, BETA[3,], lwd=2, col = rgb(0,0,0,.7))

mean(BETA[4,])
mean(BETA[5,])

# global coefficients proposed
phi = phi_generator(3, tau.grid)
Xmat_lags1 = cbind(rep(1,nrow(Xmat_lags)), Xmat_lags)
globalfit = global_qr(taus = tau.grid, phi = phi, X = Xmat_lags1, y = Yvec_lags, lambda=10, lags = 2)
global_coefs = cbind(t(globalfit$bhat), rep(4,M))
colnames(global_coefs) = colnames(real_coefs)


plot(tau.grid, alpha0(tau.grid), type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, globalfit$bhat[1,], lwd=2, col = rgb(0,0,0,.7))

plot(tau.grid, alpha1(tau.grid), type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, globalfit$bhat[2,], lwd=2, col = rgb(0,0,0,.7))

plot(tau.grid, theta1(tau.grid), type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, globalfit$bhat[3,], lwd=2, col = rgb(0,0,0,.7))

mean(globalfit$bhat[4,])
mean(globalfit$bhat[5,])

## Plot results 
coefs_results = data.frame(rbind(real_coefs, qrfit_coefs, piqr_coefs, global_coefs))
coefs_results$Method = as.factor(coefs_results$Method)
levels(coefs_results$Method) = c("Real", "QR", "piqr", "Global")


require(reshape2)
require(ggplot2)
coefs_results.m <- melt(coefs_results, id.var = "Method")
ggplot(data = coefs_results.m, aes(x=variable, y=value, fill=Method)) + geom_boxplot()


# Squared error 

sr_qr = (real_coefs - qrfit_coefs)^2
sr_qr[,6] = rep(2,M)
sr_piqr = (real_coefs - piqr_coefs)^2
sr_piqr[,6] = rep(3,M)
sr_global = (real_coefs - global_coefs)^2
sr_global[,6] = rep(4,M)

plot(tau.grid, sr_qr[,1], type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, sr_piqr[,1], lwd=2, col = rgb(0,1,0,.7))
lines(tau.grid, sr_global[,1], lwd=2, col = rgb(0,0,0,.7))

plot(tau.grid, sr_qr[,2], type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, sr_piqr[,2], lwd=2, col = rgb(0,1,0,.7))
lines(tau.grid, sr_global[,2], lwd=2, col = rgb(0,0,0,.7))

plot(tau.grid, sr_qr[,3], type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, sr_piqr[,3], lwd=2, col = rgb(0,1,0,.7))
lines(tau.grid, sr_global[,3], lwd=2, col = rgb(0,0,0,.7))

plot(tau.grid, sr_qr[,4], type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, sr_piqr[,4], lwd=2, col = rgb(0,1,0,.7))
lines(tau.grid, sr_global[,4], lwd=2, col = rgb(0,0,0,.7))

plot(tau.grid, sr_qr[,5], type = 'l', col='blue', lwd=2, lty='dotted')
lines(tau.grid, sr_piqr[,5], lwd=2, col = rgb(0,1,0,.7))
lines(tau.grid, sr_global[,5], lwd=2, col = rgb(0,0,0,.7))

sr_coefs = data.frame(rbind(sr_qr, sr_piqr, sr_global))
sr_coefs$Method = as.factor(sr_coefs$Method)

sr.m <- melt(sr_coefs, id.var = "Method")
ggplot(data = sr.m, aes(x=variable, y=value, fill=Method)) + geom_boxplot()

# Sum of squared errors
ssr_qr = apply(sr_qr, 2, sum)
ssr_qr[6] = 2
ssr_piqr = apply(sr_piqr, 2, sum)
ssr_piqr[6] = 3
ssr_global = apply(sr_global, 2, sum)
ssr_global[6] = 4

ssr_coefs = data.frame(rbind(ssr_qr, ssr_piqr, ssr_global))
ssr_coefs$Method = as.factor(ssr_coefs$Method)
levels(ssr_coefs$Method) = c("QR", "piqr", "Global")

ssr.m <- melt(ssr_coefs, id.var = "Method")
ggplot(data = ssr.m, aes(x=variable, y=value)) + geom_point(aes(colour=Method))

#means(TODO)

ssrqr_mean = apply(ssrqr, 2, mean)
ssrglobal_mean = apply(ssrglobal, 2, mean)
ssrpiqr_mean = apply(ssrpiqr, 2, mean)

plot(1:5, ssrglobal_mean, type = 'l', col='blue', lwd=2, lty='dotted')
lines(1:5, ssrqr_mean, lwd=2, col = rgb(0,1,0,.7))
lines(1:5, ssrconquer_mean, lwd=2, col = rgb(1,0,0,.7))
lines(1:5, ssrpiqr_mean, lwd=2, col = rgb(0,0,1,.7))


ssr_means = cbind(mean(ssrqr), mean(ssrconquer_mean), mean(ssrglobal), mean(ssrpiqr))
colnames(ssr_means) = c("QR", "Conquer", "Proposed", "PIQR")
ssr_means
         