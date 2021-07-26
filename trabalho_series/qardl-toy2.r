
# set.seed(1)
library(quantreg)
library(forecast)
library(qrcmNP)
source('utils.R')
source('qardl-gen.r')

tau.grid = seq(from=.01,to=.99, by=.02)

#### Simulating the sample paths Y-1 X-1 ####

Y1_X1 = simulate_qardl(N=500)

Y1_X1.Y = Y1_X1$Y
Y1_X1.Z = Y1_X1$Z

N = length(Y1_X1.Y)

acf(Y1_X1.Y, lwd=16, lend=3, col='gray')
pacf(Y1_X1.Y, lwd=16, lend=3, col='gray')

plot(Y1_X1.Y[2:N]~Y1_X1.Y[1:(N-1)], pch=16, col=rgb(0,0,0,.4))
plot(Y1_X1.Y[2:N]~Y1_X1.Z[1:(N-1)], pch=16, col=rgb(0,0,0,.4))

hist(Y1_X1.Y, border=NA, breaks="FD")
ts.plot(Y1_X1.Y)

#### Simulating the sample paths Y-2 X-2 ####

Y2_X2 = simulate_qardl(N=500, Ylag = 2, Xlag = 2)

Y2_X2.Y = Y2_X2$Y
Y2_X2.Z = Y2_X2$Z

N = length(Y2_X2.Y)

acf(Y2_X2.Y, lwd=16, lend=3, col='gray')
pacf(Y2_X2.Y, lwd=16, lend=3, col='gray')

# No corr
plot(Y2_X2.Y[2:N]~Y2_X2.Y[1:(N-1)], pch=16, col=rgb(0,0,0,.4))
# Corr
plot(Y2_X2.Y[3:N]~Y2_X2.Y[1:(N-2)], pch=16, col=rgb(0,0,0,.4))
#No corr
plot(Y2_X2.Y[2:N]~Y2_X2.Z[1:(N-1)], pch=16, col=rgb(0,0,0,.4))
plot(Y2_X2.Y[3:N]~Y2_X2.Z[1:(N-2)], pch=16, col=rgb(0,0,0,.4))

hist(Y2_X2.Y, border=NA, breaks="FD")
ts.plot(Y2_X2.Y)


#### Simulating the sample paths Y-5 X-1 ####

Y5_X1 = simulate_qardl(N=500, Ylag = 5, Xlag = 1)
Y5_X1.Y = Y5_X1$Y
Y5_X1.Z = Y5_X1$Z

N = length(Y5_X1.Y)

acf(Y5_X1.Y, lwd=16, lend=3, col='gray')
pacf(Y5_X1.Y, lwd=16, lend=3, col='gray')

# No corr
plot(Y5_X1.Y[2:N]~Y5_X1.Y[1:(N-1)], pch=16, col=rgb(0,0,0,.4))
# Corr
plot(Y5_X1.Y[6:N]~Y5_X1.Y[1:(N-5)], pch=16, col=rgb(0,0,0,.4))
#No corr
plot(Y5_X1.Y[2:N]~Y5_X1.Z[1:(N-1)], pch=16, col=rgb(0,0,0,.4))

hist(Y5_X1.Y, border=NA, breaks="FD")
ts.plot(Y5_X1.Y)

# quantile regression
#Yvec = Y
#Xmat = cbind(Y[1:(T-1)], Z[1:(T-1)])


Yvec_lags_10 = Y[11:N]
Xmat_lags_10 = cbind(Y[1:(N - 1)], Z[1:(N - 1)])
for (i in 2:10){
        Xmat_lags_10 = cbind(Xmat_lags_10, c(NA, Y[1:(N - i)]), c(NA, Z[1:(N -i)]))[-1, ] 
}
colnames(Xmat_lags_10) = c("Yt-1", "Xt-1", "Yt-2", "Xt-2", "Yt-3", "Xt-3", 
                           "Yt-4", "Xt-4", "Yt-5", "Xt-5", "Yt-6", "Xt-6", 
                           "Yt-7", "Xt-7", "Yt-8", "Xt-8", "Yt-9", "Xt-9", "Yt-10", "Xt-10")

M = length(tau.grid)

real_coefs = cbind(alpha0(tau.grid), rep(0,M), theta1(tau.grid), rep(0,M), rep(0,M), rep(0,M), rep(0,M), rep(0,M), rep(0,M), alpha1(tau.grid),rep(0,M))
for (i in 6:10){
        real_coefs=cbind(real_coefs, rep(0, M), rep(0, M))
}
colnames(real_coefs) = c("Intercept", "Yt-1", "Xt-1", "Yt-2", "Xt-2", "Yt-3", "Xt-3", 
                         "Yt-4", "Xt-4", "Yt-5", "Xt-5", "Yt-6", "Xt-6", 
                         "Yt-7", "Xt-7", "Yt-8", "Xt-8", "Yt-9", "Xt-9", "Yt-10", "Xt-10")

qrfit = rq(Yvec_lags_10~Xmat_lags_10, tau = tau.grid)
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
Xmat_lags1 = cbind(rep(1,nrow(Xmat_lags_10)), Xmat_lags_10)
globalfit = global_qr(taus = tau.grid, phi = phi, X = Xmat_lags1, y = Yvec_lags_10, lambda=10, lags = 10)
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
         