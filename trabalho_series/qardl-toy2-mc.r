
#### Function to plot generated data ####

# The function qdbplot below allows us to investigate
# if this inequality holds when v00, v10 and v01
# are chosen to be Beta quantile functions:
qdbplot = function(a10, b10, a01, b01, a00, b00) {
        tau = seq(from = 0.00001,
                  to = .99999,
                  by = .00001)
        v00prime = 1 / dbeta(qbeta(tau, a00, b00), a00, b00)
        v10prime = 1 / dbeta(qbeta(tau, a10, b10), a10, b10)
        v01prime = 1 / dbeta(qbeta(tau, a01, b01), a01, b01)
        plot(tau, log(v10prime + v01prime) - log(v00prime), type = 'l')
        abline(h = 0, col = 'red')
}

# A function for visualizing the Beta densities:
dbplot = function(a, b) {
        tau = seq(from = 0,
                  to = 1,
                  by = .01)
        plot(tau, dbeta(tau, a, b))
}

# Conditional density of Y[t] given Y[t-1] = 1 and X[t-1]=1
q11 = function(tau) {
        1 / dbeta(v10(tau), a10, b10) + 1 / dbeta(v01(tau), a01, b01) - 1 / dbeta(v00(tau), a00, b00)
}

#### Functions to generate QADL data ####

# This code simulates a sample path of size T from a simple QARDL(1,1) model:
# Q(𝜏| ℱ[t-1]) = α₀(𝜏) + α₁(𝜏)*Y[t-1] + θ₁(𝜏)*Z[t-1]

tau.grid = seq(from = .01, to = .99, by = .02)
#tau.grid = seq(from=.1, to=.9, by=.1)
M = length(tau.grid)

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


# Parameters for the Beta quantile functions v10, v01 and v00
{
        a10 = 1
        b10 = 3
        a01 = 5
        b01 = 1
        a00 = 3
        b00 = 1
}
qdbplot(a10, b10, a01, b01, a00, b00)

# Conditional quantile function of Y[t] given Y[t-1] = 1 and X[t-1]=0
v10 = function (tau)
        qbeta(tau, a10, b10)
# Conditional density of Y[t] given Y[t-1] = 1 and X[t-1]=0
dbplot(a10, b10)

# Conditional quantile function of Y[t] given Y[t-1] = 0 and X[t-1]=1
v01 = function(tau)
        qbeta(tau, a01, b01)
# Conditional density of Y[t] given Y[t-1] = 0 and X[t-1]=1
dbplot(a01, b01)

# Conditional quantile function of Y[t] given Y[t-1] = 0 and X[t-1]=0
v00 = function(tau)
        qbeta(tau, a00, b00)
# Conditional density of Y[t] given Y[t-1] = 0 and X[t-1]=0
dbplot(a00, b00)

# Conditional quantile function of Y[t] given Y[t-1] = 1 and X[t-1]=1
v11 = function(tau)
        v10(tau) + v01(tau) - v00(tau)
# Plot of v11 (must be nondecreasing and with range contained in [0,1])
plot(tau.grid, v11(tau.grid))

plot(v11(tau.grid), 1 / q11(tau.grid), xlim = c(0, 1))

# Functional parameters for the quantile regression equation
alpha0 = v00
alpha1 = function(tau)
        v10(tau) - v00(tau)
theta1 = function(tau)
        v01(tau) - v00(tau)

# Quantile function of Y given ℱ[t-1]
Q = function(tau, Y.current, Z.current) {
        alpha0(tau) + alpha1(tau) * Y.current + theta1(tau) * Z.current
}

#### Testing global QR equivalence with PIQR ####

library(qrcmNP)
source('utils.R')

### Simulate Data ###
# Arbitrary starting point (Y0,Z0)
Y0 = runif(1)
Z0 = runif(1)

Y = Z = numeric()
Z.current = Z0
Y.current = Y0
N = 1001
for (t in 1:N) {
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
N = length(Y)

Yvec = Y[2:N]
Xmat = cbind(Y[1:(N - 1)], Z[1:(N - 1)])
colnames(Xmat) = c("Yt-1", "Xt-1")

tau.grid = seq(from = .01, to = .99, by = .02)
tau.grid2 = seq(from = .05, to=.99, by=.05)
tau.grid3 = seq(from = .1, to=.9, by=.1)

# True parameters
# quantile regression
real_coefs = cbind(alpha0(tau.grid), alpha1(tau.grid), theta1(tau.grid))
colnames(real_coefs) = c("Intercept", "Yt-1", "Xt-1")

L = 5
M = length(tau.grid)
M2 = length(tau.grid2)
M3 = length(tau.grid3)

# global coefficients Sottile
# Estimate qadl time series coefficients using qrcm
k.user = L-1
p = M

# thinner grid
fo3o = piqr(Yvec ~ Xmat, formula.p = ~ slp(p, k = k.user), lambda = 1)
fo4o = slp(tau.grid, k = k.user)
PHI = cbind(1, fo4o)
BETA = fo3o$coef$lambda1 %*% t(PHI)
piqr_coefs = BETA

p = M2
fo3o = piqr(Yvec ~ Xmat, formula.p = ~ slp(p, k = k.user), lambda = 1)
fo4o = slp(tau.grid2, k = k.user)
PHI = cbind(1, fo4o)
BETA = fo3o$coef$lambda1 %*% t(PHI)
piqr_coefs2 = BETA

# thicker grid
p = M3
fo3o = piqr(Yvec ~ Xmat, formula.p = ~ slp(p, k = k.user), lambda = 1)
fo4o = slp(tau.grid3, k = k.user)
PHI = cbind(1, fo4o)
BETA = fo3o$coef$lambda1 %*% t(PHI)
piqr_coefs3 = BETA

# global coefficients equivalent to Sottile
phi = phi_generator(L, tau.grid)
phi2 = phi_generator(L, tau.grid2)
phi3 = phi_generator(L, tau.grid3)

Xmat_1 = cbind(rep(1, nrow(Xmat)), Xmat)
globalfit = global_qr(taus = tau.grid, phi = phi, X = Xmat_1, y = Yvec, lambda = 1, lags = 0, f = "piqr", w = F)
global_coefs = globalfit$bhat

globalfit = global_qr(taus = tau.grid2, phi = phi2, X = Xmat_1, y = Yvec, lambda = 1, lags = 0, f = "piqr", w = F)
global_coefs2 = globalfit$bhat

globalfit = global_qr(taus = tau.grid3, phi = phi3, X = Xmat_1, y = Yvec, lambda = 1, lags = 0, f = "piqr", w = F)
global_coefs3 = globalfit$bhat


plot(tau.grid, real_coefs[,1], type = 'l', col='black', lwd=2, lty='dotted', ylab = "Estimated alpha1")
lines(tau.grid, piqr_coefs[1,], lwd=2, col = rgb(0,0,1,.7))
lines(tau.grid, global_coefs[1,], lwd=2, col = rgb(1,0,0,.7))
legend(x = "bottomright",
       legend = c("Real", "Sottile Impl.", "Author Impl."),
       lty = c(2, 1, 1),
       col = c("black", rgb(0,0,1,.7), rgb(1,0,0,.7)),
       lwd = 2)

plot(tau.grid, real_coefs[,2], type = 'l', col='black', lwd=2, lty='dotted', ylab = "Estimated alpha1", ylim = c(-0.7, -0.1))
lines(tau.grid, piqr_coefs[2,], lwd=2, col = rgb(0,0,1,.7))
lines(tau.grid, global_coefs[2,], lwd=2, col = rgb(1,0,0,.7))
legend(x = "top",
       legend = c("Real", "Sottile Impl.", "Author Impl."),
       lty = c(2, 1, 1),
       col = c("black", rgb(0,0,1,.7), rgb(1,0,0,.7)),
       lwd = 2)


plot(tau.grid, real_coefs[,3], type = 'l', col='black', lwd=2, lty='dotted', ylab = "Estimated theta1")
lines(tau.grid, piqr_coefs[3,], lwd=2, col = rgb(0,0,1,.7))
lines(tau.grid, global_coefs[3,], lwd=2, col = rgb(1,0,0,.7))
legend(x = "topright",
       legend = c("Real", "Sottile Impl.", "Author Impl."),
       lty = c(2, 1, 1),
       col = c("black", rgb(0,0,1,.7), rgb(1,0,0,.7)),
       lwd = 2)


### Plots without real ###
plot(tau.grid, piqr_coefs[1,], type = 'l', col='black', lwd=2, lty='dotted', ylab = "Estimated alpha1")
lines(tau.grid, global_coefs[1,], lwd=2, col = rgb(1,0,0,.7))
legend(x = "bottomright",
       legend = c("Sottile Impl.", "Author Impl."),
       lty = c(2, 1),
       col = c("black", rgb(1,0,0,.7)),
       lwd = 2)

plot(tau.grid, piqr_coefs[2,], type = 'l', col='black', lwd=2, lty='dotted', ylab = "Estimated alpha1", ylim = c(-0.7, -0.1))
lines(tau.grid, global_coefs[2,], lwd=2, col = rgb(1,0,0,.7))
legend(x = "top",
       legend = c( "Sottile Impl.", "Author Impl."),
       lty = c(2, 1),
       col = c("black", rgb(1,0,0,.7)),
       lwd = 2)


plot(tau.grid, piqr_coefs[3,], type = 'l', col='black', lwd=2, lty='dotted', ylab = "Estimated theta1")
lines(tau.grid, global_coefs[3,], lwd=2, col = rgb(1,0,0,.7))
legend(x = "topright",
       legend = c("Sottile Impl.", "Author Impl."),
       lty = c(2, 1),
       col = c("black", rgb(1,0,0,.7)),
       lwd = 2)


### Lags = 1 ####
nrep = 200

tau.grid = seq(from = .1, to=.9, by=.1)
M = length(tau.grid)
L = 5

betas_qr_1 = betas_piqr_1 = betas_piqrW_1 = betas_piqrWL_1 = list()
betas_gLasso_1 = betas_gLassoW_1 = betas_gLassoWL_1 = list()
betas_qr_1[[1]] = betas_piqr_1[[1]] = betas_piqrW_1[[1]] = betas_piqrWL_1[[1]] = alpha0(tau.grid)
betas_gLasso_1[[1]] = betas_glassoW_1[[1]] = betas_glassoWL_1[[1]] = alpha0(tau.grid)
betas_qr_1[[2]] = betas_piqr_1[[2]] = betas_piqrW_1[[2]] = betas_piqrWL_1[[2]] = alpha1(tau.grid)
betas_gLasso_1[[2]] = betas_gLassoW_1[[2]] = betas_gLassoWL_1[[2]] = alpha1(tau.grid)
betas_qr_1[[3]] = betas_piqr_1[[3]] = betas_piqrW_1[[3]] = betas_piqrWL_1[[3]] = theta1(tau.grid)
betas_gLasso_1[[3]] = betas_gLassoW_1[[3]] = betas_gLassoWL_1[[3]] = theta1(tau.grid)

set.seed(205650)
for (i in 1:nrep) {
        # Arbitrary starting point (Y0,Z0)
        Y0 = runif(1)
        Z0 = runif(1)
        
        Y = Z = numeric()
        Z.current = Z0
        Y.current = Y0
        N = 501
        for (t in 1:N) {
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
        N = length(Y)
        
        Yvec = Y[2:N]
        Xmat = cbind(Y[1:(N - 1)], Z[1:(N - 1)])
        colnames(Xmat) = c("Yt-1", "Xt-1")
        
        
        qrfit = rq(Yvec ~ Xmat, tau = tau.grid)
        qrfit_coefs = coef(qrfit)
        
        # global coefficients proposed
        phi = phi_generator(L, tau.grid)
        Xmat_1 = cbind(rep(1, nrow(Xmat)), Xmat)
        
        piqr = global_qr(taus = tau.grid, phi = phi, X = Xmat_1, y = Yvec, lambda = 1, lags = 0, f = "piqr", w = F)
        piqr_coefs = piqr$bhat
        
        piqrW = global_qr(taus = tau.grid, phi = phi, X = Xmat_1, y = Yvec, lambda = 1, lags = 0, f = "piqr", w = T)
        piqrW_coefs = piqrW$bhat
        
        piqrWL = global_qr(taus = tau.grid, phi = phi, X = Xmat_1, y = Yvec, lambda = 1, lags = 1, f = "piqr", w = T)
        piqrWL_coefs = piqrWL$bhat
        
        gLasso = global_qr(taus = tau.grid, phi = phi, X = Xmat_1, y = Yvec, lambda = 1, lags = 0, f = "gLasso", w = F)
        gLasso_coefs = gLasso$bhat
        
        gLassoW = global_qr(taus = tau.grid, phi = phi, X = Xmat_1, y = Yvec, lambda = 1, lags = 0, f = "gLasso", w = T)
        gLassoW_coefs = gLassoW$bhat
        
        gLassoWL = global_qr(taus = tau.grid, phi = phi, X = Xmat_1, y = Yvec, lambda = 1, lags = 1, f = "gLasso", w = T)
        gLassoWL_coefs = gLassoWL$bhat
        
        if(piqr_coefs == "error" | piqrW_coefs == "error" | piqrWL_coefs == "error" | 
           gLasso_coefs == "error" | gLassoW_coefs == "error" | gLassoWL_coefs == "error") next
        
        
        for(i in 1:nrow(gLasso_coefs)){
                
                betas_qr_1[[i]] = rbind(betas_qr_1[[i]], qrfit_coefs[i,])
                betas_piqr_1[[i]] = rbind(betas_piqr_1[[i]], piqr_coefs[i,])
                betas_piqrW_1[[i]] = rbind(betas_piqrW_1[[i]], piqrW_coefs[i,])
                betas_piqrWL_1[[i]] = rbind(betas_piqrWL_1[[i]], piqrWL_coefs[i,])
                betas_gLasso_1[[i]] = rbind(betas_gLasso_1[[i]], gLasso_coefs[i,])
                betas_gLassoW_1[[i]] = rbind(betas_gLassoW_1[[i]], gLassoW_coefs[i,])
                betas_gLassoWL_1[[i]] = rbind(betas_gLassoWL_1[[i]], gLassoWL_coefs[i,])
        }
}



se_qr_1 = list()
se_piqr_1 = list()
se_piqrW_1 = list()
se_piqrWL_1 = list()
se_gLasso_1 = list()
se_gLassoW_1 = list()
se_gLassoWL_1 = list()


for (i in 1:length(betas_qr_1)){
        se_qr_1[[i]] = t(apply(betas_qr_1[[i]][-1,], 1, function(row) {(row - betas_qr_1[[i]][1,])^2}))
        se_piqr_1[[i]] = t(apply(betas_piqr_1[[i]][-1,], 1, function(row) {(row - betas_piqr_1[[i]][1,])^2}))
        se_piqrW_1[[i]] = t(apply(betas_piqrW_1[[i]][-1,], 1, function(row) {(row - betas_piqrW_1[[i]][1,])^2}))
        se_piqrWL_1[[i]] = t(apply(betas_piqrWL_1[[i]][-1,], 1, function(row) {(row - betas_piqrWL_1[[i]][1,])^2}))
        se_gLasso_1[[i]] = t(apply(betas_gLasso_1[[i]][-1,], 1, function(row) {(row - betas_gLasso_1[[i]][1,])^2}))
        se_gLassoW_1[[i]] = t(apply(betas_gLassoW_1[[i]][-1,], 1, function(row) {(row - betas_gLassoW_1[[i]][1,])^2}))
        se_gLassoWL_1[[i]] = t(apply(betas_gLassoWL_1[[i]][-1,], 1, function(row) {(row - betas_gLassoWL_1[[i]][1,])^2}))
}

mse_1 = list()
for (i in 1:length(betas_qr_1)){
        mse_1[[i]] = apply(se_qr_1[[i]], 2, mean)
        mse_1[[i]] = rbind(mse_1[[i]], apply(se_piqr_1[[i]], 2, mean))
        mse_1[[i]] = rbind(mse_1[[i]], apply(se_piqrW_1[[i]], 2, mean))
        mse_1[[i]] = rbind(mse_1[[i]], apply(se_piqrWL_1[[i]], 2, mean))
        mse_1[[i]] = rbind(mse_1[[i]], apply(se_gLasso_1[[i]], 2, mean))
        mse_1[[i]] = rbind(mse_1[[i]], apply(se_gLassoW_1[[i]], 2, mean))
        mse_1[[i]] = rbind(mse_1[[i]], apply(se_gLassoWL_1[[i]], 2, mean))
}

mse_taus_df_1 <- lapply(seq_along(mse_1), function(i) {
        df = data.frame(mse_1[[i]])
        df$var = rep(i, nrow(mse_1[[i]]))
        df$method = c(1,2,3,4,5,6,7)
        return(df)
})

mse_taus_df_1$method = as.factor(mse_taus_df_1$method)
mse_taus_df_1$var = as.factor(mse_taus_df_1$var)

mse_df_1 = data.frame("mse" = apply(mse_1[[1]], 1, mean), "var" = rep(1, 7), "method" = c(1,2,3,4,5,6,7))
mse_df_1 = rbind(mse_df_1, data.frame("mse" = apply(mse_1[[2]], 1, mean), "var" = rep(2, 7), "method"= c(1,2,3,4,5,6,7)))
mse_df_1 = rbind(mse_df_1, data.frame("mse" = apply(mse_1[[3]], 1, mean), "var" = rep(3, 7), "method"= c(1,2,3,4,5,6,7)))

mse_df_1$method = as.factor(mse_df_1$method)
mse_df_1$var = as.factor(mse_df_1$var)
levels(mse_df_1$method) = c("QR", "piqr", "piqrW", "piqrWL", "gLasso", "gLassoW", "gLassoWL")
levels(mse_df_1$var) = c("c", "Yt-1", "Xt-1")

library(reshape2)
mse.m_1 <- melt(mse_df_1, id.vars  = c("method", "var"))
ggplot(data = mse.m_1, aes(x=var, y=value)) + geom_point(aes(colour=method)) + xlab("Variable") + ylab("MSE")


#### Lags = 2 ####
### QR vs. PIQR vs. gLasso ###
nrep = 200
tau.grid = seq(from=.1, to=.9, by=.1)
M = length(tau.grid)
L = 5

betas_qr = betas_piqr = betas_piqrW = betas_piqrWL = list()
betas_gLasso = betas_gLassoW = betas_gLassoWL = list()
betas_qr[[1]] = betas_piqr[[1]] = betas_piqrW[[1]] = betas_piqrWL[[1]] = alpha0(tau.grid)
betas_gLasso[[1]] = betas_gLassoW[[1]] = betas_gLassoWL[[1]] = alpha0(tau.grid)
betas_qr[[2]] = betas_piqr[[2]] = betas_piqrW[[2]] = betas_piqrWL[[2]] = alpha1(tau.grid)
betas_gLasso[[2]] = betas_gLassoW[[2]] = betas_gLassoWL[[2]] = alpha1(tau.grid)
betas_qr[[3]] = betas_piqr[[3]] = betas_piqrW[[3]] = betas_piqrWL[[3]] = theta1(tau.grid)
betas_gLasso[[3]] = betas_gLassoW[[3]] = betas_gLassoWL[[3]] = theta1(tau.grid)
betas_qr[[4]] = betas_piqr[[4]] = betas_piqrW[[4]] = betas_piqrWL[[4]] = rep(0,M)
betas_gLasso[[4]] = betas_gLassoW[[4]] = betas_gLassoWL[[4]] = rep(0,M)
betas_qr[[5]] = betas_piqr[[5]] = betas_piqrW[[5]] = betas_piqrWL[[5]] = rep(0,M)
betas_gLasso[[5]] = betas_gLassoW[[5]] = betas_gLassoWL[[5]] = rep(0,M)


set.seed(205652)
for (i in 1:nrep) {
        # Arbitrary starting point (Y0,Z0)
        Y0 = runif(1)
        Z0 = runif(1)
        
        Y = Z = numeric()
        Z.current = Z0
        Y.current = Y0
        N = 501
        for (t in 1:N) {
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
        N = length(Y)
        
        Yvec_lags = Y[3:N]
        Xmat_lags = cbind(Y[1:(N - 1)], Z[1:(N - 1)], c(NA, Y[1:(N - 2)]), c(NA, Z[1:(N -2)]))[-1, ]
        colnames(Xmat_lags) = c("Yt-1", "Xt-1", "Yt-2", "Xt-2")
        
        # True and estimated functional parameters
        # quantile regression
        real_coefs = cbind(alpha0(tau.grid), alpha1(tau.grid), theta1(tau.grid), rep(0, M), rep(0, M))
        colnames(real_coefs) = c("Intercept", "Yt-1", "Xt-1", "Yt-2", "Xt-2")
        
        qrfit = rq(Yvec_lags ~ Xmat_lags, tau = tau.grid)
        qrfit_coefs = coef(qrfit)
        
        phi = phi_generator(L, tau.grid)
        Xmat_lags1 = cbind(rep(1, nrow(Xmat_lags)), Xmat_lags)
        
        piqr = global_qr(taus = tau.grid, phi = phi, X = Xmat_lags1, y = Yvec_lags, lambda = 1, lags = 0, f = "piqr", w = F)
        piqr_coefs = piqr$bhat
        
        piqrW = global_qr(taus = tau.grid, phi = phi, X = Xmat_lags1, y = Yvec_lags, lambda = 1, lags = 0, f = "piqr", w = T)
        piqrW_coefs = piqrW$bhat
        
        piqrWL = global_qr(taus = tau.grid, phi = phi, X = Xmat_lags1, y = Yvec_lags, lambda = 1, lags = 2, f = "piqr", w = T)
        piqrWL_coefs = piqrWL$bhat
        
        gLasso = global_qr(taus = tau.grid, phi = phi, X = Xmat_lags1, y = Yvec_lags, lambda = 1, lags = 0, f = "gLasso", w = F)
        gLasso_coefs = gLasso$bhat
        
        gLassoW = global_qr(taus = tau.grid, phi = phi, X = Xmat_lags1, y = Yvec_lags, lambda = 1, lags = 0, f = "gLasso", w = T)
        gLassoW_coefs = gLassoW$bhat
        
        gLassoWL = global_qr(taus = tau.grid, phi = phi, X = Xmat_lags1, y = Yvec_lags, lambda = 1, lags = 2, f = "gLasso", w = T)
        gLassoWL_coefs = gLassoWL$bhat
        
        if(piqr_coefs == "error" | piqrW_coefs == "error" | piqrWL_coefs == "error" | 
           gLasso_coefs == "error" | gLassoW_coefs == "error" | gLassoWL_coefs == "error") next
        
        for(i in 1:nrow(gLasso_coefs)){
                
                betas_qr[[i]] = rbind(betas_qr[[i]], qrfit_coefs[i,])
                betas_piqr[[i]] = rbind(betas_piqr[[i]], piqr_coefs[i,])
                betas_piqrW[[i]] = rbind(betas_piqrW[[i]], piqrW_coefs[i,])
                betas_piqrWL[[i]] = rbind(betas_piqrWL[[i]], piqrWL_coefs[i,])
                betas_gLasso[[i]] = rbind(betas_gLasso[[i]], gLasso_coefs[i,])
                betas_gLassoW[[i]] = rbind(betas_gLassoW[[i]], gLassoW_coefs[i,])
                betas_gLassoWL[[i]] = rbind(betas_gLassoWL[[i]], gLassoWL_coefs[i,])
        }
}

# Comparing PIQR and gLasso with real coefs
titles = c("Intercept", "Yt-1", "Xt-1", "Yt-2", "Xt-2")
par(mfrow=c(2,2))
for(var in 1:3){
        plot(tau.grid, betas_piqr[[var]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = titles[var], ylim = c(-1,1))
        for(i in 3:(nrep)){
                lines(tau.grid, betas_piqr[[var]][i,], lwd=2, col = rgb(0,0,0,.1))
        }
        lines(tau.grid, betas_piqr[[var]][1,], lwd=2, col = rgb(1,0,0,.7))
}
par(mfrow=c(2,2))
for(var in 4:5){
        plot(tau.grid, betas_gLassoW[[var]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = titles[var], ylim=c(-0.5,0.5))
        for(i in 3:nrep){
                lines(tau.grid, betas_gLassoW[[var]][i,], lwd=2, col = rgb(0,0,0,.1))
        }
        lines(tau.grid, betas_gLassoW[[var]][1,], lwd=2, col = rgb(1,0,0,.7))
}

titles = c("Intercept", "Yt-1", "Xt-1", "Yt-2", "Xt-2")

# Comparing piqr and gLasso with real coefs
par(mfrow=c(3,2))
## Var 1
plot(tau.grid, betas_piqr[[1]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[1], " PIQR"), ylim = c(0.2,1.2))
for(i in 3:nrep){
        lines(tau.grid, betas_piqr[[1]][i,], lwd=2, col = rgb(0,0,0,.1))
}
lines(tau.grid, betas_piqr[[1]][1,], lwd=2, col = rgb(1,0,0,.7))

plot(tau.grid, betas_gLasso[[1]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[1], " gLASSO"), ylim = c(0.2,1.2))
for(i in 3:nrep){
        lines(tau.grid, betas_gLasso[[1]][i,], lwd=2, col = rgb(0,0,0,.1))
}
lines(tau.grid, betas_gLasso[[1]][1,], lwd=2, col = rgb(1,0,0,.7))

# Var 2
plot(tau.grid, betas_piqr[[2]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[2], " PIQR"), ylim = c(-0.7, -0.2))
for(i in 3:nrep){
        lines(tau.grid, betas_piqr[[2]][i,], lwd=2, col = rgb(0,0,0,.1))
}
lines(tau.grid, betas_piqr[[2]][1,], lwd=2, col = rgb(1,0,0,.7))

plot(tau.grid, betas_gLasso[[2]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[2], " gLASSO"), ylim = c(-0.7,0.2))
for(i in 3:nrep){
        lines(tau.grid, betas_gLasso[[2]][i,], lwd=2, col = rgb(0,0,0,.1))
}
lines(tau.grid, betas_gLasso[[2]][1,], lwd=2, col = rgb(1,0,0,.7))

# Var 3
plot(tau.grid, betas_piqr[[3]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[3], " PIQR"), ylim = c(-0.5, 0.5))
for(i in 3:nrep){
        lines(tau.grid, betas_piqr[[3]][i,], lwd=2, col = rgb(0,0,0,.1))
}
lines(tau.grid, betas_piqr[[3]][1,], lwd=2, col = rgb(1,0,0,.7))

plot(tau.grid, betas_gLasso[[3]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[3], " gLASSO"), ylim = c(-0.5,0.5))
for(i in 3:nrep){
        lines(tau.grid, betas_gLasso[[3]][i,], lwd=2, col = rgb(0,0,0,.1))
}
lines(tau.grid, betas_gLasso[[3]][1,], lwd=2, col = rgb(1,0,0,.7))

#Vars 4 e 5
par(mfrow=c(2,2))
for(var in 4:5){
        plot(tau.grid, betas_piqr[[var]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[var], " PIQR"), ylim=c(-0.5,0.5))
        for(i in 3:nrep){
                lines(tau.grid, betas_piqr[[var]][i,], lwd=2, col = rgb(0,0,0,.1))
        }
        lines(tau.grid, betas_piqr[[var]][1,], lwd=2, col = rgb(1,0,0,.7))
        
        plot(tau.grid, betas_gLasso[[var]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[var], " gLASSO"), ylim=c(-0.5,0.5))
        for(i in 3:nrep){
                lines(tau.grid, betas_gLasso[[var]][i,], lwd=2, col = rgb(0,0,0,.1))
        }
        lines(tau.grid, betas_gLasso[[var]][1,], lwd=2, col = rgb(1,0,0,.7))
}


# Comparing piqr with different weight approaches

titles = c("Intercept", "Yt-1", "Xt-1", "Yt-2", "Xt-2")

par(mfrow=c(3,3))
## Var 1
plot(tau.grid, betas_piqr[[1]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[1], " PIQR"), ylim = c(0.2,1.2))
for(i in 3:nrep){
        lines(tau.grid, betas_piqr[[1]][i,], lwd=2, col = rgb(0,0,0,.1))
}
lines(tau.grid, betas_piqr[[1]][1,], lwd=2, col = rgb(1,0,0,.7))

plot(tau.grid, betas_piqrW[[1]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[1], "PIQR Weighted"), ylim = c(0.2,1.2))
for(i in 3:nrep){
        lines(tau.grid, betas_piqrW[[1]][i,], lwd=2, col = rgb(0,0,0,.1))
}
lines(tau.grid, betas_piqrW[[1]][1,], lwd=2, col = rgb(1,0,0,.7))

plot(tau.grid, betas_piqrWL[[1]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[1], "PIQR Weighted Lags"), ylim = c(0.2,1.2))
for(i in 3:nrep){
        lines(tau.grid, betas_piqrWL[[1]][i,], lwd=2, col = rgb(0,0,0,.1))
}
lines(tau.grid, betas_piqrWL[[1]][1,], lwd=2, col = rgb(1,0,0,.7))


# Var 2
plot(tau.grid, betas_piqr[[2]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[2], " PIQR"), ylim = c(-0.7, -0.2))
for(i in 3:nrep){
        lines(tau.grid, betas_piqr[[2]][i,], lwd=2, col = rgb(0,0,0,.1))
}
lines(tau.grid, betas_piqr[[2]][1,], lwd=2, col = rgb(1,0,0,.7))

plot(tau.grid, betas_piqrW[[2]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[2], " PIQR Weighted"), ylim = c(-0.7,0.2))
for(i in 3:nrep){
        lines(tau.grid, betas_piqrW[[2]][i,], lwd=2, col = rgb(0,0,0,.1))
}
lines(tau.grid, betas_piqrW[[2]][1,], lwd=2, col = rgb(1,0,0,.7))

plot(tau.grid, betas_piqrWL[[2]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[2], "PIQR Weighted Lags"), ylim = c(-0.7,0.2))
for(i in 3:nrep){
        lines(tau.grid, betas_piqrWL[[2]][i,], lwd=2, col = rgb(0,0,0,.1))
}
lines(tau.grid, betas_piqrWL[[2]][1,], lwd=2, col = rgb(1,0,0,.7))


# Var 3
plot(tau.grid, betas_piqr[[3]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[3], " PIQR"), ylim = c(-0.5, 0.5))
for(i in 3:nrep){
        lines(tau.grid, betas_piqr[[3]][i,], lwd=2, col = rgb(0,0,0,.1))
}
lines(tau.grid, betas_piqr[[3]][1,], lwd=2, col = rgb(1,0,0,.7))

plot(tau.grid, betas_piqrW[[3]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[3], "PIQR Weighted"), ylim = c(-0.5,0.5))
for(i in 3:nrep){
        lines(tau.grid, betas_piqrW[[3]][i,], lwd=2, col = rgb(0,0,0,.1))
}
lines(tau.grid, betas_piqrW[[3]][1,], lwd=2, col = rgb(1,0,0,.7))

plot(tau.grid, betas_piqrWL[[3]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[3], "PIQR Weighted Lags"), ylim = c(-0.5,0.5))
for(i in 3:nrep){
        lines(tau.grid, betas_piqrWL[[3]][i,], lwd=2, col = rgb(0,0,0,.1))
}
lines(tau.grid, betas_piqrWL[[3]][1,], lwd=2, col = rgb(1,0,0,.7))


#Vars 4 e 5
par(mfrow=c(2,3))
for(var in 4:5){
        plot(tau.grid, betas_piqr[[var]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[var], " PIQR"), ylim=c(-0.5,0.5))
        for(i in 3:nrep){
                lines(tau.grid, betas_piqr[[var]][i,], lwd=2, col = rgb(0,0,0,.1))
        }
        lines(tau.grid, betas_piqr[[var]][1,], lwd=2, col = rgb(1,0,0,.2))
        
        plot(tau.grid, betas_piqrW[[var]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[var], "PIQR Weighted"), ylim=c(-0.5,0.5))
        for(i in 3:nrep){
                lines(tau.grid, betas_piqrW[[var]][i,], lwd=2, col = rgb(0,0,0,.1))
        }
        lines(tau.grid, betas_piqrW[[var]][1,], lwd=2, col = rgb(1,0,0,.2))
        
        plot(tau.grid, betas_piqrWL[[var]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[var], "PIQR Weighted Lags"), ylim=c(-0.5,0.5))
        for(i in 3:nrep){
                lines(tau.grid, betas_piqrWL[[var]][i,], lwd=2, col = rgb(0,0,0,.1))
        }
        lines(tau.grid, betas_piqrWL[[var]][1,], lwd=2, col = rgb(1,0,0,.2))
}

# Comparing gLASSO with different weight approaches
# need to finish

titles = c("Intercept", "Yt-1", "Xt-1", "Yt-2", "Xt-2")

par(mfrow=c(3,3))
## Var 1
plot(tau.grid, betas_gLasso[[1]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[1], "gLASSO"), ylim = c(0.2,1.2))
for(i in 3:nrep){
        lines(tau.grid, betas_gLasso[[1]][i,], lwd=2, col = rgb(0,0,0,.1))
}
lines(tau.grid, betas_gLasso[[1]][1,], lwd=2, col = rgb(1,0,0,.7))

plot(tau.grid, betas_piqrW[[1]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[1], "PIQR Weighted"), ylim = c(0.2,1.2))
for(i in 3:nrep){
        lines(tau.grid, betas_piqrW[[1]][i,], lwd=2, col = rgb(0,0,0,.1))
}
lines(tau.grid, betas_piqrW[[1]][1,], lwd=2, col = rgb(1,0,0,.7))

plot(tau.grid, betas_piqrWL[[1]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[1], "PIQR Weighted Lags"), ylim = c(0.2,1.2))
for(i in 3:nrep){
        lines(tau.grid, betas_piqrWL[[1]][i,], lwd=2, col = rgb(0,0,0,.1))
}
lines(tau.grid, betas_piqrWL[[1]][1,], lwd=2, col = rgb(1,0,0,.7))


# Var 2
plot(tau.grid, betas_piqr[[2]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[2], " PIQR"), ylim = c(-0.7, -0.2))
for(i in 3:nrep){
        lines(tau.grid, betas_piqr[[2]][i,], lwd=2, col = rgb(0,0,0,.1))
}
lines(tau.grid, betas_piqr[[2]][1,], lwd=2, col = rgb(1,0,0,.7))

plot(tau.grid, betas_piqrW[[2]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[2], " PIQR Weighted"), ylim = c(-0.7,0.2))
for(i in 3:nrep){
        lines(tau.grid, betas_piqrW[[2]][i,], lwd=2, col = rgb(0,0,0,.1))
}
lines(tau.grid, betas_piqrW[[2]][1,], lwd=2, col = rgb(1,0,0,.7))

plot(tau.grid, betas_piqrWL[[2]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[2], "PIQR Weighted Lags"), ylim = c(-0.7,0.2))
for(i in 3:nrep){
        lines(tau.grid, betas_piqrWL[[2]][i,], lwd=2, col = rgb(0,0,0,.1))
}
lines(tau.grid, betas_piqrWL[[2]][1,], lwd=2, col = rgb(1,0,0,.7))


# Var 3
plot(tau.grid, betas_piqr[[3]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[3], " PIQR"), ylim = c(-0.5, 0.5))
for(i in 3:nrep){
        lines(tau.grid, betas_piqr[[3]][i,], lwd=2, col = rgb(0,0,0,.1))
}
lines(tau.grid, betas_piqr[[3]][1,], lwd=2, col = rgb(1,0,0,.7))

plot(tau.grid, betas_piqrW[[3]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[3], "PIQR Weighted"), ylim = c(-0.5,0.5))
for(i in 3:nrep){
        lines(tau.grid, betas_piqrW[[3]][i,], lwd=2, col = rgb(0,0,0,.1))
}
lines(tau.grid, betas_piqrW[[3]][1,], lwd=2, col = rgb(1,0,0,.7))

plot(tau.grid, betas_piqrWL[[3]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[3], "PIQR Weighted Lags"), ylim = c(-0.5,0.5))
for(i in 3:nrep){
        lines(tau.grid, betas_piqrWL[[3]][i,], lwd=2, col = rgb(0,0,0,.1))
}
lines(tau.grid, betas_piqrWL[[3]][1,], lwd=2, col = rgb(1,0,0,.7))


#Vars 4 e 5
par(mfrow=c(2,3))
for(var in 4:5){
        plot(tau.grid, betas_piqr[[var]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[var], " PIQR"), ylim=c(-0.5,0.5))
        for(i in 3:nrep){
                lines(tau.grid, betas_piqr[[var]][i,], lwd=2, col = rgb(0,0,0,.1))
        }
        lines(tau.grid, betas_piqr[[var]][1,], lwd=2, col = rgb(1,0,0,.2))
        
        plot(tau.grid, betas_piqrW[[var]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[var], "PIQR Weighted"), ylim=c(-0.5,0.5))
        for(i in 3:nrep){
                lines(tau.grid, betas_piqrW[[var]][i,], lwd=2, col = rgb(0,0,0,.1))
        }
        lines(tau.grid, betas_piqrW[[var]][1,], lwd=2, col = rgb(1,0,0,.2))
        
        plot(tau.grid, betas_piqrWL[[var]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[var], "PIQR Weighted Lags"), ylim=c(-0.5,0.5))
        for(i in 3:nrep){
                lines(tau.grid, betas_piqrWL[[var]][i,], lwd=2, col = rgb(0,0,0,.1))
        }
        lines(tau.grid, betas_piqrWL[[var]][1,], lwd=2, col = rgb(1,0,0,.2))
}






# Comparing betas of all methods with real coefs
titles = c("Intercept", "Yt-1", "Xt-1", "Yt-2", "Xt-2")
par(mfrow=c(3,3))
for(var in 1:3){
        plot(tau.grid, betas_piqrW[[var]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[var], " PIQR"), ylim = c(-1,1))
        for(i in 3:nrep){
                lines(tau.grid, betas_piqrW[[var]][i,], lwd=2, col = rgb(0,0,0,.1))
        }
        lines(tau.grid, betas_piqrW[[var]][1,], lwd=2, col = rgb(1,0,0,.7))
        
        plot(tau.grid, betas_gLassoW[[var]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[var], " Global"), ylim = c(-1,1))
        for(i in 3:nrep){
                lines(tau.grid, betas_gLassoW[[var]][i,], lwd=2, col = rgb(0,0,0,.1))
        }
        lines(tau.grid, betas_gLassoW[[var]][1,], lwd=2, col = rgb(1,0,0,.7))
}
par(mfrow=c(2,3))
for(var in 4:5){
        plot(tau.grid, betas_qr[[var]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[var], " QR"), ylim=c(-0.1,0.1))
        for(i in 3:nrep){
                lines(tau.grid, betas_qr[[var]][i,], lwd=2, col = rgb(0,0,0,.1))
        }
        lines(tau.grid, betas_qr[[var]][1,], lwd=2, col = rgb(1,0,0,.7))
        
        plot(tau.grid, betas_piqr[[var]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[var], " PIQR"), ylim=c(-0.1,0.1))
        for(i in 3:nrep){
                lines(tau.grid, betas_piqr[[var]][i,], lwd=2, col = rgb(0,0,0,.1))
        }
        lines(tau.grid, betas_piqr[[var]][1,], lwd=2, col = rgb(1,0,0,.7))
        
        plot(tau.grid, betas_global[[var]][2,], type = 'l', col=rgb(0,0,0,.1), lwd=2, ylab = "Coef", main = paste(titles[var], " Global"), ylim=c(-0.1,0.1))
        for(i in 3:nrep){
                lines(tau.grid, betas_global[[var]][i,], lwd=2, col = rgb(0,0,0,.1))
        }
        lines(tau.grid, betas_global[[var]][1,], lwd=2, col = rgb(1,0,0,.7))
}

se_qr = list()
se_piqr = list()
se_piqrW = list()
se_piqrWL = list()
se_gLasso = list()
se_gLassoW = list()
se_gLassoWL = list()


for (i in 1:length(betas_qr)){
        se_qr[[i]] = t(apply(betas_qr[[i]][-1,], 1, function(row) {(row - betas_qr[[i]][1,])^2}))
        se_piqr[[i]] = t(apply(betas_piqr[[i]][-1,], 1, function(row) {(row - betas_piqr[[i]][1,])^2}))
        se_piqrW[[i]] = t(apply(betas_piqrW[[i]][-1,], 1, function(row) {(row - betas_piqrW[[i]][1,])^2}))
        se_piqrWL[[i]] = t(apply(betas_piqrWL[[i]][-1,], 1, function(row) {(row - betas_piqrWL[[i]][1,])^2}))
        se_gLasso[[i]] = t(apply(betas_gLasso[[i]][-1,], 1, function(row) {(row - betas_gLasso[[i]][1,])^2}))
        se_gLassoW[[i]] = t(apply(betas_gLassoW[[i]][-1,], 1, function(row) {(row - betas_gLassoW[[i]][1,])^2}))
        se_gLassoWL[[i]] = t(apply(betas_gLassoWL[[i]][-1,], 1, function(row) {(row - betas_gLassoWL[[i]][1,])^2}))
}

mse = list()
for (i in 1:length(betas_qr)){
        mse[[i]] = apply(se_qr[[i]], 2, mean)
        mse[[i]] = rbind(mse[[i]], apply(se_piqr[[i]], 2, mean))
        mse[[i]] = rbind(mse[[i]], apply(se_piqrW[[i]], 2, mean))
        mse[[i]] = rbind(mse[[i]], apply(se_piqrWL[[i]], 2, mean))
        mse[[i]] = rbind(mse[[i]], apply(se_gLasso[[i]], 2, mean))
        mse[[i]] = rbind(mse[[i]], apply(se_gLassoW[[i]], 2, mean))
        mse[[i]] = rbind(mse[[i]], apply(se_gLassoWL[[i]], 2, mean))
}

mse_taus_df <- lapply(seq_along(mse), function(i) {
        df = data.frame(mse[[i]])
        df$var = rep(i, nrow(mse[[i]]))
        df$method = c(1,2,3,4,5,6,7)
        return(df)
})

mse_taus_df$method = as.factor(mse_taus_df$method)
mse_taus_df$var = as.factor(mse_taus_df$var)

mse_df = data.frame("mse" = apply(mse[[1]], 1, mean), "var" = rep(1, 7), "method" = c(1,2,3,4,5,6,7))
mse_df = rbind(mse_df, data.frame("mse" = apply(mse[[2]], 1, mean), "var" = rep(2, 7), "method"= c(1,2,3,4,5,6,7)))
mse_df = rbind(mse_df, data.frame("mse" = apply(mse[[3]], 1, mean), "var" = rep(3, 7), "method"= c(1,2,3,4,5,6,7)))
mse_df = rbind(mse_df, data.frame("mse" = apply(mse[[4]], 1, mean), "var" = rep(4, 7), "method"= c(1,2,3,4,5,6,7)))
mse_df = rbind(mse_df, data.frame("mse" = apply(mse[[5]], 1, mean), "var" = rep(5, 7), "method"= c(1,2,3,4,5,6,7)))

mse_df$method = as.factor(mse_df$method)
mse_df$var = as.factor(mse_df$var)
levels(mse_df$method) = c("QR", "piqr", "piqrW", "piqrWL", "gLasso", "gLassoW", "gLassoWL")
levels(mse_df$var) = c("c", "Yt-1", "Xt-1", "Yt-2", "Xt-2")

library(reshape2)
mse.m <- melt(mse_df, id.vars  = c("method", "var"))
ggplot(data = mse.m, aes(x=var, y=value)) + geom_point(aes(colour=method)) + xlab("Variable") + ylab("MSE")


# Lags = 10

nrep = 1

tau.grid = seq(from = .01, to = .99, by = .02)
M = length(tau.grid)

betas_qr_10 = betas_piqr_10 = betas_global_10 = list()
betas_qr_10[[1]] = betas_piqr_10[[1]] = betas_global_10[[1]] = alpha0(tau.grid)
betas_qr_10[[2]] = betas_piqr_10[[2]] = betas_global_10[[2]] = alpha1(tau.grid)
betas_qr_10[[3]] = betas_piqr_10[[3]] = betas_global_10[[3]] = theta1(tau.grid)

for (i in 4:10){
        betas_qr_10[[i]] = betas_piqr_10[[i]] = betas_global_10[[i]] = rep(0,M)
}

set.seed(205651)
for (i in 1:nrep) {
        # Arbitrary starting point (Y0,Z0)
        Y0 = runif(1)
        Z0 = runif(1)
        
        Y = Z = numeric()
        Z.current = Z0
        Y.current = Y0
        N = 10001
        for (t in 1:N) {
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
        N = length(Y)
        
        Yvec_lags_10 = Y[11:N]
        Xmat_lags_10 = cbind(Y[1:(N - 1)], Z[1:(N - 1)])
        for (i in 2:10){
                Xmat_lags_10 = cbind(Xmat_lags_10, c(NA, Y[1:(N - i)]), c(NA, Z[1:(N -i)]))[-1, ] 
        }
        colnames(Xmat_lags_10) = c("Yt-1", "Xt-1", "Yt-2", "Xt-2", "Yt-3", "Xt-3", 
                                   "Yt-4", "Xt-4", "Yt-5", "Xt-5", "Yt-6", "Xt-6", 
                                   "Yt-7", "Xt-7", "Yt-8", "Xt-8", "Yt-9", "Xt-9", "Yt-10", "Xt-10")
        
        # True and estimated functional parameters
        # quantile regression
        real_coefs = cbind(alpha0(tau.grid), alpha1(tau.grid), theta1(tau.grid))
        for (i in 2:10){
                real_coefs=cbind(real_coefs, rep(0, M), rep(0, M))
        }
        colnames(real_coefs) = c("Intercept", "Yt-1", "Xt-1", "Yt-2", "Xt-2", "Yt-3", "Xt-3", 
                                 "Yt-4", "Xt-4", "Yt-5", "Xt-5", "Yt-6", "Xt-6", 
                                 "Yt-7", "Xt-7", "Yt-8", "Xt-8", "Yt-9", "Xt-9", "Yt-10", "Xt-10")
        
        qrfit = rq(Yvec_lags_10 ~ Xmat_lags_10, tau = tau.grid)
        qrfit_coefs = coef(qrfit)
        
        # global coefficients Sottile
        # Estimate qadl time series coefficients using qrcm
        k.user = 3
        p = M
        fo3o = piqr(Yvec_lags_10 ~ Xmat_lags_10, formula.p = ~ slp(p, k = k.user), lambda = 10)
        fo4o = slp(tau.grid, k = k.user)
        PHI = cbind(1, fo4o)
        BETA = fo3o$coef$lambda1 %*% t(PHI)
        piqr_coefs = BETA
        
        # global coefficients proposed
        phi = phi_generator(3, tau.grid)
        Xmat_lags1_10 = cbind(rep(1, nrow(Xmat_lags_10)), Xmat_lags_10)
        globalfit = global_qr(taus = tau.grid, phi = phi, X = Xmat_lags1_10, y = Yvec_lags_10, lambda = 10, lags = 10)
        global_coefs = globalfit$bhat
        
        for(i in 1:nrow(global_coefs)){
                
                betas_qr_10[[i]] = rbind(betas_qr_10[[i]], qrfit_coefs[i,])
                betas_piqr_10[[i]] = rbind(betas_piqr_10[[i]], piqr_coefs[i,])
                betas_global_10[[i]] = rbind(betas_global_10[[i]], global_coefs[i,])
        }
}


se_qr_10 = list()
se_piqr_10 = list()
se_global_10 = list()

for (i in 1:length(betas_qr_10)){
        se_qr_10[[i]] = t(apply(betas_qr_10[[i]][-1,], 1, function(row) {(row - betas_qr_10[[i]][1,])^2}))
        se_piqr_10[[i]] = t(apply(betas_piqr_10[[i]][-1,], 1, function(row) {(row - betas_piqr_10[[i]][1,])^2}))
        se_global_10[[i]] = t(apply(betas_global_10[[i]][-1,], 1, function(row) {(row - betas_global_10[[i]][1,])^2}))
}

mse_10 = list()
for (i in 1:length(betas_global_10)){
        mse_10[[i]] = apply(se_qr_10[[i]], 2, mean)
        mse_10[[i]] = rbind(mse_10[[i]], apply(se_piqr_10[[i]], 2, mean))
        mse_10[[i]] = rbind(mse_10[[i]], apply(se_global_10[[i]], 2, mean))
}

for (i in 1:length(mse_10)){
        write.csv2(se_qr_10, file = paste("se_qr_10_", i, ".csv"))
        write.csv2(se_piqr_10, file = paste("se_piqr_10_", i, ".csv"))
        write.csv2(se_global_10, file = paste("se_global_10_", i, ".csv"))
        write.csv2(mse_10, file = paste("mse_10_", i, ".csv"))
}

mse_taus_df_10 <- lapply(seq_along(mse_10), function(i) {
        df = data.frame(mse_10[[i]])
        df$var = rep(i, nrow(mse_10[[i]]))
        df$method = c(1,2,3)
        return(df)
})

mse_taus_df_10$method = as.factor(mse_taus_df_10$method)
mse_taus_df_10$var = as.factor(mse_taus_df_10$var)

mse_df_10 = data.frame("mse" = apply(mse_10[[1]], 1, mean), "var" = rep(1, 3), "method" = c(1,2,3))
mse_df_10 = rbind(mse_df_10, data.frame("mse" = apply(mse_10[[2]], 1, mean), "var" = rep(2, 3), "method"= c(1,2,3)))
mse_df_10 = rbind(mse_df_10, data.frame("mse" = apply(mse_10[[3]], 1, mean), "var" = rep(3, 3), "method"= c(1,2,3)))

mse_df_10$method = as.factor(mse_df_10$method)
mse_df_10$var = as.factor(mse_df_10$var)
levels(mse_df_10$method) = c("QR", "piqr", "Global")
levels(mse_df_10$var) = c("c", "Yt-1", "Xt-1", "Yt-2", "Xt-2", "Yt-3", "Xt-3", 
                          "Yt-4", "Xt-4", "Yt-5", "Xt-5", "Yt-6", "Xt-6", 
                          "Yt-7", "Xt-7", "Yt-8", "Xt-8", "Yt-9", "Xt-9", "Yt-10", "Xt-10")

library(reshape2)
mse.m_10 <- melt(mse_df_10, id.vars  = c("method", "var"))
ggplot(data = mse.m_10, aes(x=var, y=value)) + geom_point(aes(colour=method)) + xlab("Variable") + ylab("MSE")

