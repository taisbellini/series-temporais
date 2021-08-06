
# set.seed(1)
library(quantreg)
library(forecast)
library(qrcmNP)
source('utils.R')
source('qardl-gen.r')

#### Simulating the sample paths Y-1 X-1 ####

Y1_X1 = simulate_qardl(N=500)

Y1_X1.Y = Y1_X1$Y
Y1_X1.Z = Y1_X1$Z

Y1_X1.N = length(Y1_X1.Y)

acf(Y1_X1.Y, lwd=16, lend=3, col='gray')
pacf(Y1_X1.Y, lwd=16, lend=3, col='gray')

plot(Y1_X1.Y[2:Y1_X1.N]~Y1_X1.Y[1:(Y1_X1.N-1)], pch=16, col=rgb(0,0,0,.4))
plot(Y1_X1.Y[2:Y1_X1.N]~Y1_X1.Z[1:(Y1_X1.N-1)], pch=16, col=rgb(0,0,0,.4))

hist(Y1_X1.Y, border=NA, breaks="FD")
ts.plot(Y1_X1.Y)

#### Simulating the sample paths Y-2 X-2 ####

Y2_X2 = simulate_qardl(N=500, Ylag = 2, Xlag = 2)

Y2_X2.Y = Y2_X2$Y
Y2_X2.Z = Y2_X2$Z

Y2_X2.N = length(Y2_X2.Y)

acf(Y2_X2.Y, lwd=16, lend=3, col='gray')
pacf(Y2_X2.Y, lwd=16, lend=3, col='gray')

# No corr
plot(Y2_X2.Y[2:Y2_X2.N]~Y2_X2.Y[1:(Y2_X2.N-1)], pch=16, col=rgb(0,0,0,.4))
# Corr
plot(Y2_X2.Y[3:Y2_X2.N]~Y2_X2.Y[1:(Y2_X2.N-2)], pch=16, col=rgb(0,0,0,.4))
#No corr
plot(Y2_X2.Y[2:Y2_X2.N]~Y2_X2.Z[1:(Y2_X2.N-1)], pch=16, col=rgb(0,0,0,.4))
plot(Y2_X2.Y[3:Y2_X2.N]~Y2_X2.Z[1:(Y2_X2.N-2)], pch=16, col=rgb(0,0,0,.4))

hist(Y2_X2.Y, border=NA, breaks="FD")
ts.plot(Y2_X2.Y)


#### Simulating the sample paths Y-5 X-1 ####

Y5_X1 = simulate_qardl(N=500, Ylag = 5, Xlag = 1)
Y5_X1.Y = Y5_X1$Y
Y5_X1.Z = Y5_X1$Z

Y5_X1.N = length(Y5_X1.Y)

acf(Y5_X1.Y, lwd=16, lend=3, col='gray')
pacf(Y5_X1.Y, lwd=16, lend=3, col='gray')

# No corr
plot(Y5_X1.Y[2:Y5_X1.N]~Y5_X1.Y[1:(Y5_X1.N-1)], pch=16, col=rgb(0,0,0,.4))
# Corr
plot(Y5_X1.Y[6:Y5_X1.N]~Y5_X1.Y[1:(Y5_X1.N-5)], pch=16, col=rgb(0,0,0,.4))
#No corr
plot(Y5_X1.Y[2:Y5_X1.N]~Y5_X1.Z[1:(Y5_X1.N-1)], pch=16, col=rgb(0,0,0,.4))

hist(Y5_X1.Y, border=NA, breaks="FD")
ts.plot(Y5_X1.Y)


#### Test all scenarios with a 10 lagged DB ####

tau.grid = seq(from=.05,to=.95, by=.05)

#### Test Yt-1 Xt-1 ####
Y1_X1.Yvec = Y1_X1.Y[11:Y1_X1.N]
Y1_X1.Xmat = cbind(Y1_X1.Y[1:(Y1_X1.N - 1)], Y1_X1.Z[1:(Y1_X1.N - 1)])
for (i in 2:10){
        Y1_X1.Xmat = cbind(Y1_X1.Xmat, c(NA, Y1_X1.Y[1:(Y1_X1.N - i)]), c(NA, Y1_X1.Z[1:(Y1_X1.N -i)]))[-1, ] 
}
colnames(Y1_X1.Xmat) = c("Yt-1", "Xt-1", "Yt-2", "Xt-2", "Yt-3", "Xt-3", 
                           "Yt-4", "Xt-4", "Yt-5", "Xt-5", "Yt-6", "Xt-6", 
                           "Yt-7", "Xt-7", "Yt-8", "Xt-8", "Yt-9", "Xt-9", "Yt-10", "Xt-10")

M = length(tau.grid)

real_coefs = cbind(alpha0(tau.grid), alpha1(tau.grid), theta1(tau.grid), rep(0,M), rep(0,M), rep(0,M), rep(0,M), rep(0,M), rep(0,M), rep(0,M),rep(0,M))
for (i in 6:10){
        real_coefs=cbind(real_coefs, rep(0, M), rep(0, M))
}
real_coefs = cbind(real_coefs, rep(1,M))
colnames(real_coefs) = c("Intercept", "Yt-1", "Xt-1", "Yt-2", "Xt-2", "Yt-3", "Xt-3", 
                         "Yt-4", "Xt-4", "Yt-5", "Xt-5", "Yt-6", "Xt-6", 
                         "Yt-7", "Xt-7", "Yt-8", "Xt-8", "Yt-9", "Xt-9", "Yt-10", "Xt-10", "method")

qrfit = rq(Y1_X1.Yvec ~ Y1_X1.Xmat, tau = tau.grid)
qrfit_coefs = cbind(t(coef(qrfit)), rep(2,M))
colnames(qrfit_coefs) = c(colnames(real_coefs))

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


#### Monte Carlo sim Y1 X1####

nrep = 500

tau.grid = seq(from = .1, to=.9, by=.1)
M = length(tau.grid)
L = 5

Y1_X1.betas_qr = Y1_X1.betas_piqr = Y1_X1.betas_gLasso = list()
Y1_X1.betas_piqrW = Y1_X1.betas_gLassoW = list()
Y1_X1.betas_piqrWL = Y1_X1.betas_gLassoWL = list()

for (i in 1:21){
        Y1_X1.betas_qr[[i]] = Y1_X1.betas_piqr[[i]] = Y1_X1.betas_gLasso[[i]] = 
                Y1_X1.betas_piqrW[[i]] = Y1_X1.betas_gLassoW[[i]] = Y1_X1.betas_piqrWL[[i]] = Y1_X1.betas_gLassoWL[[i]] = rep(0,M)
}

Y1_X1.betas_qr[[1]] = Y1_X1.betas_piqr[[1]] = Y1_X1.betas_gLasso[[1]] = alpha0(tau.grid)
Y1_X1.betas_piqrW[[1]] = Y1_X1.betas_gLassoW[[1]] = alpha0(tau.grid)
Y1_X1.betas_piqrWL[[1]] = Y1_X1.betas_gLassoWL[[1]] = alpha0(tau.grid)
Y1_X1.betas_qr[[2]] = Y1_X1.betas_piqr[[2]] = Y1_X1.betas_gLasso[[2]] = alpha1(tau.grid)
Y1_X1.betas_piqrW[[2]] = Y1_X1.betas_gLassoW[[2]] = alpha1(tau.grid)
Y1_X1.betas_piqrWL[[2]] = Y1_X1.betas_gLassoWL[[2]] = alpha1(tau.grid)
Y1_X1.betas_qr[[3]] = Y1_X1.betas_piqr[[3]] = Y1_X1.betas_gLasso[[3]] = theta1(tau.grid)
Y1_X1.betas_piqrW[[3]] = Y1_X1.betas_gLassoW[[3]] = theta1(tau.grid)
Y1_X1.betas_piqrWL[[3]] = Y1_X1.betas_gLassoWL[[3]] = theta1(tau.grid)

set.seed(205650)
for (i in 1:nrep) {
        Y1_X1 = simulate_qardl(N=500)
        
        Y1_X1.Y = Y1_X1$Y
        Y1_X1.Z = Y1_X1$Z
        
        Y1_X1.N = length(Y1_X1.Y)
        
        Y1_X1.Yvec = Y1_X1.Y[11:Y1_X1.N]
        Y1_X1.Xmat = cbind(Y1_X1.Y[1:(Y1_X1.N - 1)], Y1_X1.Z[1:(Y1_X1.N - 1)])
        for (i in 2:10){
                Y1_X1.Xmat = cbind(Y1_X1.Xmat, c(NA, Y1_X1.Y[1:(Y1_X1.N - i)]), c(NA, Y1_X1.Z[1:(Y1_X1.N -i)]))[-1, ] 
        }
        colnames(Y1_X1.Xmat) = c("Yt-1", "Xt-1", "Yt-2", "Xt-2", "Yt-3", "Xt-3", 
                                 "Yt-4", "Xt-4", "Yt-5", "Xt-5", "Yt-6", "Xt-6", 
                                 "Yt-7", "Xt-7", "Yt-8", "Xt-8", "Yt-9", "Xt-9", "Yt-10", "Xt-10")
        
        M = length(tau.grid)
        
        qrfit = rq(Y1_X1.Yvec ~ Y1_X1.Xmat, tau = tau.grid)
        qrfit_coefs = coef(qrfit)
        
        # global coefficients proposed
        phi = phi_generator(L, tau.grid)
        Y1_X1.Xmat1 = cbind(rep(1, nrow(Y1_X1.Xmat)), Y1_X1.Xmat)
        
        piqr = global_qr(taus = tau.grid, phi = phi, X = Y1_X1.Xmat1, y = Y1_X1.Yvec, lambda = 1, lags = 0, f = "piqr", w = F)
        piqr_coefs = piqr$bhat
        
        piqrW = global_qr(taus = tau.grid, phi = phi, X = Y1_X1.Xmat1, y = Y1_X1.Yvec, lambda = 1, lags = 0, f = "piqr", w = T)
        piqrW_coefs = piqrW$bhat
        
        piqrWL = global_qr(taus = tau.grid, phi = phi, X = Y1_X1.Xmat1, y = Y1_X1.Yvec, lambda = 1, lags = 10, f = "piqr", w = T)
        piqrWL_coefs = piqrWL$bhat
        
        gLasso = global_qr(taus = tau.grid, phi = phi, X = Y1_X1.Xmat1, y = Y1_X1.Yvec, lambda = 1, lags = 0, f = "gLasso", w = F)
        gLasso_coefs = gLasso$bhat
        
        gLassoW = global_qr(taus = tau.grid, phi = phi, X = Y1_X1.Xmat1, y = Y1_X1.Yvec, lambda = 1, lags = 0, f = "gLasso", w = T)
        gLassoW_coefs = gLassoW$bhat
        
        gLassoWL = global_qr(taus = tau.grid, phi = phi, X = Y1_X1.Xmat1, y = Y1_X1.Yvec, lambda = 1, lags = 10, f = "gLasso", w = T)
        gLassoWL_coefs = gLassoWL$bhat
        
        if(typeof(piqr_coefs) == "character" | typeof(piqrW_coefs) == "character" | typeof(piqrWL_coefs) == "character" | 
           typeof(gLasso_coefs) == "character" | typeof(gLassoW_coefs) == "character" | typeof(gLassoWL_coefs) == "character") next
        
        
        for(i in 1:nrow(gLasso_coefs)){
                
                Y1_X1.betas_qr[[i]] = rbind(Y1_X1.betas_qr[[i]], qrfit_coefs[i,])
                Y1_X1.betas_piqr[[i]] = rbind(Y1_X1.betas_piqr[[i]], piqr_coefs[i,])
                Y1_X1.betas_piqrW[[i]] = rbind(Y1_X1.betas_piqrW[[i]], piqrW_coefs[i,])
                Y1_X1.betas_piqrWL[[i]] = rbind(Y1_X1.betas_piqrWL[[i]], piqrWL_coefs[i,])
                Y1_X1.betas_gLasso[[i]] = rbind(Y1_X1.betas_gLasso[[i]], gLasso_coefs[i,])
                Y1_X1.betas_gLassoW[[i]] = rbind(Y1_X1.betas_gLassoW[[i]], gLassoW_coefs[i,])
                Y1_X1.betas_gLassoWL[[i]] = rbind(Y1_X1.betas_gLassoWL[[i]], gLassoWL_coefs[i,])
        }
}

#### Y1X1 SE ####

Y1_X1.se_qr = list()
Y1_X1.se_piqr = list()
Y1_X1.se_piqrW = list()
Y1_X1.se_piqrWL = list()
Y1_X1.se_gLasso = list()
Y1_X1.se_gLassoW = list()
Y1_X1.se_gLassoWL = list()


for (i in 1:length(Y1_X1.betas_qr)){
        Y1_X1.se_qr[[i]] = t(apply(Y1_X1.betas_qr[[i]][-1,], 1, function(row) {(row - Y1_X1.betas_qr[[i]][1,])^2}))
        Y1_X1.se_piqr[[i]] = t(apply(Y1_X1.betas_piqr[[i]][-1,], 1, function(row) {(row - Y1_X1.betas_piqr[[i]][1,])^2}))
        Y1_X1.se_piqrW[[i]] = t(apply(Y1_X1.betas_piqrW[[i]][-1,], 1, function(row) {(row - Y1_X1.betas_piqrW[[i]][1,])^2}))
        Y1_X1.se_piqrWL[[i]] = t(apply(Y1_X1.betas_piqrWL[[i]][-1,], 1, function(row) {(row - Y1_X1.betas_piqrWL[[i]][1,])^2}))
        Y1_X1.se_gLasso[[i]] = t(apply(Y1_X1.betas_gLasso[[i]][-1,], 1, function(row) {(row - Y1_X1.betas_gLasso[[i]][1,])^2}))
        Y1_X1.se_gLassoW[[i]] = t(apply(Y1_X1.betas_gLassoW[[i]][-1,], 1, function(row) {(row - Y1_X1.betas_gLassoW[[i]][1,])^2}))
        Y1_X1.se_gLassoWL[[i]] = t(apply(Y1_X1.betas_gLassoWL[[i]][-1,], 1, function(row) {(row - Y1_X1.betas_gLassoWL[[i]][1,])^2}))
}

#### Y1X1 MSE ####
Y1_X1.mse = list()
for (i in 1:length(Y1_X1.betas_qr)){
        Y1_X1.mse[[i]] = mean(apply(Y1_X1.se_qr[[i]], 2, mean))
        Y1_X1.mse[[i]] = c(Y1_X1.mse[[i]], mean(apply(Y1_X1.se_piqr[[i]], 2, mean)))
        Y1_X1.mse[[i]] = c(Y1_X1.mse[[i]], mean(apply(Y1_X1.se_piqrW[[i]], 2, mean)))
        Y1_X1.mse[[i]] = c(Y1_X1.mse[[i]], mean(apply(Y1_X1.se_piqrWL[[i]], 2, mean)))
        Y1_X1.mse[[i]] = c(Y1_X1.mse[[i]], mean(apply(Y1_X1.se_gLasso[[i]], 2, mean)))
        Y1_X1.mse[[i]] = c(Y1_X1.mse[[i]], mean(apply(Y1_X1.se_gLassoW[[i]], 2, mean)))
        Y1_X1.mse[[i]] = c(Y1_X1.mse[[i]], mean(apply(Y1_X1.se_gLassoWL[[i]], 2, mean)))
        names(Y1_X1.mse[[i]]) = c("QR","piqr","piqrW","piqrWL","gLasso","gLassoW","gLassoWL")
}

Y1_X1.mse.table = as.data.frame(do.call(rbind, Y1_X1.mse))
rownames(Y1_X1.mse.table) = c("c", "Yt-1", "Xt-1", "Yt-2", "Xt-2", "Yt-3", "Xt-3", 
                              "Yt-4", "Xt-4", "Yt-5", "Xt-5", "Yt-6", "Xt-6", 
                              "Yt-7", "Xt-7", "Yt-8", "Xt-8", "Yt-9", "Xt-9", "Yt-10", "Xt-10")
Y1_X1.mse.table.rel = Y1_X1.mse.table[c(1,2,3),]

png("img/Y1X1_mse_table.png", width=350,height=200,bg = "white")
grid.table(format(t(Y1_X1.mse.table.rel), scientific = T))
dev.off()


Y1_X1.mse_taus_df <- lapply(seq_along(Y1_X1.mse), function(i) {
        df = data.frame(Y1_X1.mse[[i]])
        df$var = rep(i, nrow(Y1_X1.mse[[i]]))
        df$method = c(1,2,3,4,5,6,7)
        return(df)
})

Y1_X1.mse_taus_df$method = as.factor(Y1_X1.mse_taus_df$method)
Y1_X1.mse_taus_df$var = as.factor(Y1_X1.mse_taus_df$var)

Y1_X1.mse_df = data.frame("mse" = apply(Y1_X1.mse[[1]], 1, mean), "var" = rep(1, 7), "method" = c(1,2,3,4,5,6,7))
for (i in 2:length(Y1_X1.mse)){
        Y1_X1.mse_df = rbind(Y1_X1.mse_df, data.frame("mse" = apply(Y1_X1.mse[[i]], 1, mean), "var" = rep(i, 7), "method"= c(1,2,3,4,5,6,7)))
}

Y1_X1.mse_df$method = as.factor(Y1_X1.mse_df$method)
Y1_X1.mse_df$var = as.factor(Y1_X1.mse_df$var)
levels(Y1_X1.mse_df$method) = c("QR", "piqr", "piqrW", "piqrWL", "gLasso", "gLassoW", "gLassoWL")
levels(Y1_X1.mse_df$var) = c("c", "Yt-1", "Xt-1", "Yt-2", "Xt-2", "Yt-3", "Xt-3", 
                             "Yt-4", "Xt-4", "Yt-5", "Xt-5", "Yt-6", "Xt-6", 
                             "Yt-7", "Xt-7", "Yt-8", "Xt-8", "Yt-9", "Xt-9", "Yt-10", "Xt-10")

library(reshape2)
X1_Y1.mse.m_1 <- melt(Y1_X1.mse_df, id.vars  = c("method", "var"))
ggplot(data = X1_Y1.mse.m_1, aes(x=var, y=value)) + geom_point(aes(colour=method), position=position_jitter(w=0.02)) + xlab("Variable") + ylab("MSE")


#### Variable Selection Y1X1 ####

Y1_X1.vs = list(length(Y1_X1.betas_qr))

for (i in 1:length(Y1_X1.betas_qr)){
        Y1_X1.vs[[i]] = ifelse(is.null(nrow(Y1_X1.betas_qr[[i]][rowSums(Y1_X1.betas_qr[[i]]) == 0,])), 0, nrow(Y1_X1.betas_qr[[i]][rowSums(Y1_X1.betas_qr[[i]]) == 0,])) 
        Y1_X1.vs[[i]] = c(Y1_X1.vs[[i]], ifelse(is.null(nrow(Y1_X1.betas_piqr[[i]][rowSums(Y1_X1.betas_piqr[[i]]) == 0,])), 0, nrow(Y1_X1.betas_piqr[[i]][rowSums(Y1_X1.betas_piqr[[i]]) == 0,])))
        Y1_X1.vs[[i]] = c(Y1_X1.vs[[i]], ifelse(is.null(nrow(Y1_X1.betas_piqrW[[i]][rowSums(Y1_X1.betas_piqrW[[i]]) == 0,])), 0, nrow(Y1_X1.betas_piqrW[[i]][rowSums(Y1_X1.betas_piqrW[[i]]) == 0,])))
        Y1_X1.vs[[i]] = c(Y1_X1.vs[[i]], ifelse(is.null(nrow(Y1_X1.betas_piqrWL[[i]][rowSums(Y1_X1.betas_piqrWL[[i]]) == 0,])), 0, nrow(Y1_X1.betas_piqrWL[[i]][rowSums(Y1_X1.betas_piqrWL[[i]]) == 0,])))
        Y1_X1.vs[[i]] = c(Y1_X1.vs[[i]], ifelse(is.null(nrow(Y1_X1.betas_gLasso[[i]][rowSums(Y1_X1.betas_gLasso[[i]]) == 0,])), 0, nrow(Y1_X1.betas_gLasso[[i]][rowSums(Y1_X1.betas_gLasso[[i]]) == 0,])))
        Y1_X1.vs[[i]] = c(Y1_X1.vs[[i]], ifelse(is.null(nrow(Y1_X1.betas_gLassoW[[i]][rowSums(Y1_X1.betas_gLassoW[[i]]) == 0,])), 0, nrow(Y1_X1.betas_gLassoW[[i]][rowSums(Y1_X1.betas_gLassoW[[i]]) == 0,])))
        Y1_X1.vs[[i]] = c(Y1_X1.vs[[i]], ifelse(is.null(nrow(Y1_X1.betas_gLassoWL[[i]][rowSums(Y1_X1.betas_gLassoWL[[i]]) == 0,])), 0, nrow(Y1_X1.betas_gLassoWL[[i]][rowSums(Y1_X1.betas_gLassoWL[[i]]) == 0,])))
        names(Y1_X1.vs[[i]]) = c("QR","piqr","piqrW","piqrWL","gLasso","gLassoW","gLassoWL")
}

Y1_X1.vs.table = as.data.frame(do.call(rbind, Y1_X1.vs))
rownames(Y1_X1.vs.table) = c("c", "Yt-1", "Xt-1", "Yt-2", "Xt-2", "Yt-3", "Xt-3", 
                             "Yt-4", "Xt-4", "Yt-5", "Xt-5", "Yt-6", "Xt-6", 
                             "Yt-7", "Xt-7", "Yt-8", "Xt-8", "Yt-9", "Xt-9", "Yt-10", "Xt-10")
library(gridExtra)

png("Y1X1_vs_table.png", width=480,height=480,bg = "white")
grid.table(Y1_X1.vs.table)
dev.off()

Y1_X1.vs_taus_df <- lapply(seq_along(Y1_X1.vs), function(i) {
        df = data.frame(Y1_X1.vs[[i]])
        df$var = i
        df$method = c(1,2,3,4,5,6,7)
        return(df)
})

Y1_X1.vs_taus_df$method = as.factor(Y1_X1.vs_taus_df$method)
Y1_X1.vs_taus_df$var = as.factor(Y1_X1.vs_taus_df$var)

Y1_X1.vs_df = data.frame("vs" = Y1_X1.vs[[1]], "var" = rep(1, 7), "method" = c(1,2,3,4,5,6,7))
for (i in 2:length(Y1_X1.vs)){
        Y1_X1.vs_df = rbind(Y1_X1.vs_df, data.frame("vs" = Y1_X1.vs[[i]], "var" = rep(i, 7), "method"= c(1,2,3,4,5,6,7)))
}

Y1_X1.vs_df$method = as.factor(Y1_X1.vs_df$method)
Y1_X1.vs_df$var = as.factor(Y1_X1.vs_df$var)
levels(Y1_X1.vs_df$method) = c("QR", "piqr", "piqrW", "piqrWL", "gLasso", "gLassoW", "gLassoWL")
levels(Y1_X1.vs_df$var) = c("c", "Yt-1", "Xt-1", "Yt-2", "Xt-2", "Yt-3", "Xt-3", 
                            "Yt-4", "Xt-4", "Yt-5", "Xt-5", "Yt-6", "Xt-6", 
                            "Yt-7", "Xt-7", "Yt-8", "Xt-8", "Yt-9", "Xt-9", "Yt-10", "Xt-10")

library(reshape2)
Y1_X1.vs.m_1 <- melt(Y1_X1.vs_df, id.vars  = c("method", "var"))
ggplot(data = Y1_X1.vs.m_1, aes(x=var, y=value,)) + geom_point(aes(colour=method), position=position_jitter(w=0.02)) + 
        xlab("Variable") + ylab("# of times var was set to zero")


#### Monte Carlo sim Y5 X1####

nrep = 500

tau.grid = seq(from = .1, to=.9, by=.1)
M = length(tau.grid)
L = 5

Y5_X1.betas_qr = Y5_X1.betas_piqr = Y5_X1.betas_gLasso = list()
Y5_X1.betas_piqrW = Y5_X1.betas_gLassoW = list()
Y5_X1.betas_piqrWL = Y5_X1.betas_gLassoWL = list()

for (i in 1:21){
        Y5_X1.betas_qr[[i]] = Y5_X1.betas_piqr[[i]] = Y5_X1.betas_gLasso[[i]] = 
                Y5_X1.betas_piqrW[[i]] = Y5_X1.betas_gLassoW[[i]] = Y5_X1.betas_piqrWL[[i]] = Y5_X1.betas_gLassoWL[[i]] = rep(0,M)
}

Y5_X1.betas_qr[[1]] = Y5_X1.betas_piqr[[1]] = Y5_X1.betas_gLasso[[1]] = alpha0(tau.grid)
Y5_X1.betas_piqrW[[1]] = Y5_X1.betas_gLassoW[[1]] = alpha0(tau.grid)
Y5_X1.betas_piqrWL[[1]] = Y5_X1.betas_gLassoWL[[1]] = alpha0(tau.grid)
Y5_X1.betas_qr[[10]] = Y5_X1.betas_piqr[[10]] = Y5_X1.betas_gLasso[[10]] = alpha1(tau.grid)
Y5_X1.betas_piqrW[[10]] = Y5_X1.betas_gLassoW[[10]] = alpha1(tau.grid)
Y5_X1.betas_piqrWL[[10]] = Y5_X1.betas_gLassoWL[[10]] = alpha1(tau.grid)
Y5_X1.betas_qr[[3]] = Y5_X1.betas_piqr[[3]] = Y5_X1.betas_gLasso[[3]] = theta1(tau.grid)
Y5_X1.betas_piqrW[[3]] = Y5_X1.betas_gLassoW[[3]] = theta1(tau.grid)
Y5_X1.betas_piqrWL[[3]] = Y5_X1.betas_gLassoWL[[3]] = theta1(tau.grid)

set.seed(205650)
for (i in 1:nrep) {
        Y5_X1 = simulate_qardl(N=500, Ylag = 5, Xlag = 1)
        
        Y5_X1.Y = Y5_X1$Y
        Y5_X1.Z = Y5_X1$Z
        
        Y5_X1.N = length(Y5_X1.Y)
        
        Y5_X1.Yvec = Y5_X1.Y[11:Y5_X1.N]
        Y5_X1.Xmat = cbind(Y5_X1.Y[1:(Y5_X1.N - 1)], Y5_X1.Z[1:(Y5_X1.N - 1)])
        for (i in 2:10){
                Y5_X1.Xmat = cbind(Y5_X1.Xmat, c(NA, Y5_X1.Y[1:(Y5_X1.N - i)]), c(NA, Y5_X1.Z[1:(Y5_X1.N -i)]))[-1, ] 
        }
        colnames(Y5_X1.Xmat) = c("Yt-1", "Xt-1", "Yt-2", "Xt-2", "Yt-3", "Xt-3", 
                                 "Yt-4", "Xt-4", "Yt-5", "Xt-5", "Yt-6", "Xt-6", 
                                 "Yt-7", "Xt-7", "Yt-8", "Xt-8", "Yt-9", "Xt-9", "Yt-10", "Xt-10")
        
        M = length(tau.grid)
        
        qrfit = rq(Y5_X1.Yvec ~ Y5_X1.Xmat, tau = tau.grid)
        qrfit_coefs = coef(qrfit)
        
        # global coefficients proposed
        phi = phi_generator(L, tau.grid)
        Y5_X1.Xmat1 = cbind(rep(1, nrow(Y5_X1.Xmat)), Y5_X1.Xmat)
        
        piqr = global_qr(taus = tau.grid, phi = phi, X = Y5_X1.Xmat1, y = Y5_X1.Yvec, lambda = 1, lags = 0, f = "piqr", w = F)
        piqr_coefs = piqr$bhat
        
        piqrW = global_qr(taus = tau.grid, phi = phi, X = Y5_X1.Xmat1, y = Y5_X1.Yvec, lambda = 1, lags = 0, f = "piqr", w = T)
        piqrW_coefs = piqrW$bhat
        
        piqrWL = global_qr(taus = tau.grid, phi = phi, X = Y5_X1.Xmat1, y = Y5_X1.Yvec, lambda = 1, lags = 10, f = "piqr", w = T)
        piqrWL_coefs = piqrWL$bhat
        
        gLasso = global_qr(taus = tau.grid, phi = phi, X = Y5_X1.Xmat1, y = Y5_X1.Yvec, lambda = 1, lags = 0, f = "gLasso", w = F)
        gLasso_coefs = gLasso$bhat
        
        gLassoW = global_qr(taus = tau.grid, phi = phi, X = Y5_X1.Xmat1, y = Y5_X1.Yvec, lambda = 1, lags = 0, f = "gLasso", w = T)
        gLassoW_coefs = gLassoW$bhat
        
        gLassoWL = global_qr(taus = tau.grid, phi = phi, X = Y5_X1.Xmat1, y = Y5_X1.Yvec, lambda = 1, lags = 10, f = "gLasso", w = T)
        gLassoWL_coefs = gLassoWL$bhat
        
        if(typeof(piqr_coefs) == "character" | typeof(piqrW_coefs) == "character" | typeof(piqrWL_coefs) == "character" | 
           typeof(gLasso_coefs) == "character" | typeof(gLassoW_coefs) == "character" | typeof(gLassoWL_coefs) == "character") next
        
        
        for(i in 1:nrow(gLasso_coefs)){
                
                Y5_X1.betas_qr[[i]] = rbind(Y5_X1.betas_qr[[i]], qrfit_coefs[i,])
                Y5_X1.betas_piqr[[i]] = rbind(Y5_X1.betas_piqr[[i]], piqr_coefs[i,])
                Y5_X1.betas_piqrW[[i]] = rbind(Y5_X1.betas_piqrW[[i]], piqrW_coefs[i,])
                Y5_X1.betas_piqrWL[[i]] = rbind(Y5_X1.betas_piqrWL[[i]], piqrWL_coefs[i,])
                Y5_X1.betas_gLasso[[i]] = rbind(Y5_X1.betas_gLasso[[i]], gLasso_coefs[i,])
                Y5_X1.betas_gLassoW[[i]] = rbind(Y5_X1.betas_gLassoW[[i]], gLassoW_coefs[i,])
                Y5_X1.betas_gLassoWL[[i]] = rbind(Y5_X1.betas_gLassoWL[[i]], gLassoWL_coefs[i,])
        }
}

#### SE ####
## difference of estimated coefficients with real ones for each row (nrep) ## 
Y5_X1.se_qr = list()
Y5_X1.se_piqr = list()
Y5_X1.se_piqrW = list()
Y5_X1.se_piqrWL = list()
Y5_X1.se_gLasso = list()
Y5_X1.se_gLassoW = list()
Y5_X1.se_gLassoWL = list()


for (i in 1:length(Y5_X1.betas_qr)){
        Y5_X1.se_qr[[i]] = t(apply(Y5_X1.betas_qr[[i]][-1,], 1, function(row) {(row - Y5_X1.betas_qr[[i]][1,])^2}))
        Y5_X1.se_piqr[[i]] = t(apply(Y5_X1.betas_piqr[[i]][-1,], 1, function(row) {(row - Y5_X1.betas_piqr[[i]][1,])^2}))
        Y5_X1.se_piqrW[[i]] = t(apply(Y5_X1.betas_piqrW[[i]][-1,], 1, function(row) {(row - Y5_X1.betas_piqrW[[i]][1,])^2}))
        Y5_X1.se_piqrWL[[i]] = t(apply(Y5_X1.betas_piqrWL[[i]][-1,], 1, function(row) {(row - Y5_X1.betas_piqrWL[[i]][1,])^2}))
        Y5_X1.se_gLasso[[i]] = t(apply(Y5_X1.betas_gLasso[[i]][-1,], 1, function(row) {(row - Y5_X1.betas_gLasso[[i]][1,])^2}))
        Y5_X1.se_gLassoW[[i]] = t(apply(Y5_X1.betas_gLassoW[[i]][-1,], 1, function(row) {(row - Y5_X1.betas_gLassoW[[i]][1,])^2}))
        Y5_X1.se_gLassoWL[[i]] = t(apply(Y5_X1.betas_gLassoWL[[i]][-1,], 1, function(row) {(row - Y5_X1.betas_gLassoWL[[i]][1,])^2}))
}

#### Y5X1 MSE ####
Y5_X1.mse = list()
for (i in 1:length(Y5_X1.betas_qr)){
        Y5_X1.mse[[i]] = mean(apply(Y5_X1.se_qr[[i]], 2, mean))
        Y5_X1.mse[[i]] = c(Y5_X1.mse[[i]], mean(apply(Y5_X1.se_piqr[[i]], 2, mean)))
        Y5_X1.mse[[i]] = c(Y5_X1.mse[[i]], mean(apply(Y5_X1.se_piqrW[[i]], 2, mean)))
        Y5_X1.mse[[i]] = c(Y5_X1.mse[[i]], mean(apply(Y5_X1.se_piqrWL[[i]], 2, mean)))
        Y5_X1.mse[[i]] = c(Y5_X1.mse[[i]], mean(apply(Y5_X1.se_gLasso[[i]], 2, mean)))
        Y5_X1.mse[[i]] = c(Y5_X1.mse[[i]], mean(apply(Y5_X1.se_gLassoW[[i]], 2, mean)))
        Y5_X1.mse[[i]] = c(Y5_X1.mse[[i]], mean(apply(Y5_X1.se_gLassoWL[[i]], 2, mean)))
        names(Y5_X1.mse[[i]]) = c("QR","piqr","piqrW","piqrWL","gLasso","gLassoW","gLassoWL")
}

Y5_X1.mse.table = as.data.frame(do.call(rbind, Y5_X1.mse))
rownames(Y5_X1.mse.table) = c("c", "Yt-1", "Xt-1", "Yt-2", "Xt-2", "Yt-3", "Xt-3", 
                             "Yt-4", "Xt-4", "Yt-5", "Xt-5", "Yt-6", "Xt-6", 
                             "Yt-7", "Xt-7", "Yt-8", "Xt-8", "Yt-9", "Xt-9", "Yt-10", "Xt-10")
Y5_X1.mse.table.rel = Y5_X1.mse.table[c(1,3,10),]

png("img/Y5X1_mse_table.png", width=350,height=200,bg = "white")
grid.table(format(t(Y5_X1.mse.table.rel), scientific = T))
dev.off()


Y5_X1.mse_taus_df <- lapply(seq_along(Y5_X1.mse), function(i) {
        df = data.frame(Y5_X1.mse[[i]])
        df$var = rep(i, nrow(Y5_X1.mse[[i]]))
        df$method = c(1,2,3,4,5,6,7)
        return(df)
})

Y5_X1.mse_taus_df$method = as.factor(Y5_X1.mse_taus_df$method)
Y5_X1.mse_taus_df$var = as.factor(Y5_X1.mse_taus_df$var)

Y5_X1.mse_df = data.frame("mse" = Y5_X1.mse[[1]], "var" = rep(1, 7), "method" = c(1,2,3,4,5,6,7))
for (i in c(3,10)){
        Y5_X1.mse_df = rbind(Y5_X1.mse_df, data.frame("mse" = Y5_X1.mse[[i]], "var" = rep(i, 7), "method"= c(1,2,3,4,5,6,7)))
}

Y5_X1.mse_df$method = as.factor(Y5_X1.mse_df$method)
Y5_X1.mse_df$var = as.factor(Y5_X1.mse_df$var)
levels(Y5_X1.mse_df$method) = c("QR", "piqr", "piqrW", "piqrWL", "gLasso", "gLassoW", "gLassoWL")
levels(Y5_X1.mse_df$var) = c("c", "Xt-1", "Yt-5")

library(reshape2)
Y5_X1.mse.m_1 <- melt(Y5_X1.mse_df, id.vars  = c("method", "var"))
ggplot(data = Y5_X1.mse.m_1, aes(x=var, y=value)) + geom_point(aes(colour=method), position=position_jitter(w=0.02)) + xlab("Variable") + ylab("MSE")


#### SSE ####
Y5_X1.sse = list()
for (i in 1:length(Y5_X1.betas_qr)){
        Y5_X1.sse[[i]] = apply(Y5_X1.se_qr[[i]], 2, mean)
        Y5_X1.sse[[i]] = rbind(Y5_X1.sse[[i]], apply(Y5_X1.se_piqr[[i]], 2, sum))
        Y5_X1.sse[[i]] = rbind(Y5_X1.sse[[i]], apply(Y5_X1.se_piqrW[[i]], 2, sum))
        Y5_X1.sse[[i]] = rbind(Y5_X1.sse[[i]], apply(Y5_X1.se_piqrWL[[i]], 2, sum))
        Y5_X1.sse[[i]] = rbind(Y5_X1.sse[[i]], apply(Y5_X1.se_gLasso[[i]], 2, sum))
        Y5_X1.sse[[i]] = rbind(Y5_X1.sse[[i]], apply(Y5_X1.se_gLassoW[[i]], 2, sum))
        Y5_X1.sse[[i]] = rbind(Y5_X1.sse[[i]], apply(Y5_X1.se_gLassoW[[i]], 2, sum))
}

Y5_X1.sse_taus_df <- lapply(seq_along(Y5_X1.sse), function(i) {
        df = data.frame(Y5_X1.sse[[i]])
        df$var = rep(i, nrow(Y5_X1.sse[[i]]))
        df$method = c(1,2,3,4,5,6,7)
        return(df)
})

Y5_X1.sse_taus_df$method = as.factor(Y5_X1.sse_taus_df$method)
Y5_X1.sse_taus_df$var = as.factor(Y5_X1.sse_taus_df$var)

Y5_X1.sse_df = data.frame("mse" = apply(Y5_X1.sse[[1]], 1, sum), "var" = rep(1, 7), "method" = c(1,2,3,4,5,6,7))
for (i in 2:length(Y5_X1.sse)){
        Y5_X1.sse_df = rbind(Y5_X1.sse_df, data.frame("mse" = apply(Y5_X1.sse[[i]], 1, sum), "var" = rep(i, 7), "method"= c(1,2,3,4,5,6,7)))
}

Y5_X1.sse_df$method = as.factor(Y5_X1.sse_df$method)
Y5_X1.sse_df$var = as.factor(Y5_X1.sse_df$var)
levels(Y5_X1.sse_df$method) = c("QR", "piqr", "piqrW", "piqrWL", "gLasso", "gLassoW", "gLassoWL")
levels(Y5_X1.sse_df$var) = c("c", "Yt-1", "Xt-1", "Yt-2", "Xt-2", "Yt-3", "Xt-3", 
                             "Yt-4", "Xt-4", "Yt-5", "Xt-5", "Yt-6", "Xt-6", 
                             "Yt-7", "Xt-7", "Yt-8", "Xt-8", "Yt-9", "Xt-9", "Yt-10", "Xt-10")

library(reshape2)
Y5_X1.sse.m_1 <- melt(Y5_X1.sse_df, id.vars  = c("method", "var"))
ggplot(data = Y5_X1.sse.m_1, aes(x=var, y=value)) + geom_point(aes(colour=method), position=position_jitter(w=0.02)) + xlab("Variable") + ylab("SSE")

#### Variable Selection ####

Y5_X1.vs = list(length(Y5_X1.betas_qr))

for (i in 1:length(Y5_X1.betas_qr)){
        Y5_X1.vs[[i]] = ifelse(is.null(nrow(Y5_X1.betas_qr[[i]][rowSums(Y5_X1.betas_qr[[i]]) == 0,])), 0, nrow(Y5_X1.betas_qr[[i]][rowSums(Y5_X1.betas_qr[[i]]) == 0,])) 
        Y5_X1.vs[[i]] = c(Y5_X1.vs[[i]], ifelse(is.null(nrow(Y5_X1.betas_piqr[[i]][rowSums(Y5_X1.betas_piqr[[i]]) == 0,])), 0, nrow(Y5_X1.betas_piqr[[i]][rowSums(Y5_X1.betas_piqr[[i]]) == 0,])))
        Y5_X1.vs[[i]] = c(Y5_X1.vs[[i]], ifelse(is.null(nrow(Y5_X1.betas_piqrW[[i]][rowSums(Y5_X1.betas_piqrW[[i]]) == 0,])), 0, nrow(Y5_X1.betas_piqrW[[i]][rowSums(Y5_X1.betas_piqrW[[i]]) == 0,])))
        Y5_X1.vs[[i]] = c(Y5_X1.vs[[i]], ifelse(is.null(nrow(Y5_X1.betas_piqrWL[[i]][rowSums(Y5_X1.betas_piqrWL[[i]]) == 0,])), 0, nrow(Y5_X1.betas_piqrWL[[i]][rowSums(Y5_X1.betas_piqrWL[[i]]) == 0,])))
        Y5_X1.vs[[i]] = c(Y5_X1.vs[[i]], ifelse(is.null(nrow(Y5_X1.betas_gLasso[[i]][rowSums(Y5_X1.betas_gLasso[[i]]) == 0,])), 0, nrow(Y5_X1.betas_gLasso[[i]][rowSums(Y5_X1.betas_gLasso[[i]]) == 0,])))
        Y5_X1.vs[[i]] = c(Y5_X1.vs[[i]], ifelse(is.null(nrow(Y5_X1.betas_gLassoW[[i]][rowSums(Y5_X1.betas_gLassoW[[i]]) == 0,])), 0, nrow(Y5_X1.betas_gLassoW[[i]][rowSums(Y5_X1.betas_gLassoW[[i]]) == 0,])))
        Y5_X1.vs[[i]] = c(Y5_X1.vs[[i]], ifelse(is.null(nrow(Y5_X1.betas_gLassoWL[[i]][rowSums(Y5_X1.betas_gLassoWL[[i]]) == 0,])), 0, nrow(Y5_X1.betas_gLassoWL[[i]][rowSums(Y5_X1.betas_gLassoWL[[i]]) == 0,])))
        names(Y5_X1.vs[[i]]) = c("QR","piqr","piqrW","piqrWL","gLasso","gLassoW","gLassoWL")
}

Y5_X1.vs.table = as.data.frame(do.call(rbind, Y5_X1.vs))
rownames(Y5_X1.vs.table) = c("c", "Yt-1", "Xt-1", "Yt-2", "Xt-2", "Yt-3", "Xt-3", 
                             "Yt-4", "Xt-4", "Yt-5", "Xt-5", "Yt-6", "Xt-6", 
                             "Yt-7", "Xt-7", "Yt-8", "Xt-8", "Yt-9", "Xt-9", "Yt-10", "Xt-10")
library(gridExtra)

png("Y5X1_vs_table.png", width=480,height=480,bg = "white")
grid.table(Y5_X1.vs.table)
dev.off()

Y5_X1.vs_taus_df <- lapply(seq_along(Y5_X1.vs), function(i) {
        df = data.frame(Y5_X1.vs[[i]])
        df$var = i
        df$method = c(1,2,3,4,5,6,7)
        return(df)
})

Y5_X1.vs_taus_df$method = as.factor(Y5_X1.vs_taus_df$method)
Y5_X1.vs_taus_df$var = as.factor(Y5_X1.vs_taus_df$var)

Y5_X1.vs_df = data.frame("vs" = Y5_X1.vs[[1]], "var" = rep(1, 7), "method" = c(1,2,3,4,5,6,7))
for (i in 2:length(Y5_X1.vs)){
        Y5_X1.vs_df = rbind(Y5_X1.vs_df, data.frame("vs" = Y5_X1.vs[[i]], "var" = rep(i, 7), "method"= c(1,2,3,4,5,6,7)))
}

Y5_X1.vs_df$method = as.factor(Y5_X1.vs_df$method)
Y5_X1.vs_df$var = as.factor(Y5_X1.vs_df$var)
levels(Y5_X1.vs_df$method) = c("QR", "piqr", "piqrW", "piqrWL", "gLasso", "gLassoW", "gLassoWL")
levels(Y5_X1.vs_df$var) = c("c", "Yt-1", "Xt-1", "Yt-2", "Xt-2", "Yt-3", "Xt-3", 
                             "Yt-4", "Xt-4", "Yt-5", "Xt-5", "Yt-6", "Xt-6", 
                             "Yt-7", "Xt-7", "Yt-8", "Xt-8", "Yt-9", "Xt-9", "Yt-10", "Xt-10")

library(reshape2)
Y5_X1.vs.m_1 <- melt(Y5_X1.vs_df, id.vars  = c("method", "var"))
ggplot(data = Y5_X1.vs.m_1, aes(x=var, y=value)) + geom_point(aes(colour=method), position=position_jitter(w=0.02)) + xlab("Variable") + ylab("# of times var was set to zero")















#### Y10X1 ####
nrep = 500

tau.grid = seq(from = .1, to=.9, by=.1)
M = length(tau.grid)
L = 5

Y10X1.betas_qr = Y10X1.betas_piqr = Y10X1.betas_gLasso = list()
Y10X1.betas_piqrW = Y10X1.betas_gLassoW = list()
Y10X1.betas_piqrWL = Y10X1.betas_gLassoWL = list()

for (i in 1:21){
        Y10X1.betas_qr[[i]] = Y10X1.betas_piqr[[i]] = Y10X1.betas_gLasso[[i]] = 
                Y10X1.betas_piqrW[[i]] = Y10X1.betas_gLassoW[[i]] = Y10X1.betas_piqrWL[[i]] = Y10X1.betas_gLassoWL[[i]] = rep(0,M)
}

Y10X1.betas_qr[[1]] = Y10X1.betas_piqr[[1]] = Y10X1.betas_gLasso[[1]] = alpha0(tau.grid)
Y10X1.betas_piqrW[[1]] = Y10X1.betas_gLassoW[[1]] = alpha0(tau.grid)
Y10X1.betas_piqrWL[[1]] = Y10X1.betas_gLassoWL[[1]] = alpha0(tau.grid)
Y10X1.betas_qr[[20]] = Y10X1.betas_piqr[[20]] = Y10X1.betas_gLasso[[20]] = alpha1(tau.grid)
Y10X1.betas_piqrW[[20]] = Y10X1.betas_gLassoW[[20]] = alpha1(tau.grid)
Y10X1.betas_piqrWL[[20]] = Y10X1.betas_gLassoWL[[20]] = alpha1(tau.grid)
Y10X1.betas_qr[[3]] = Y10X1.betas_piqr[[3]] = Y10X1.betas_gLasso[[3]] = theta1(tau.grid)
Y10X1.betas_piqrW[[3]] = Y10X1.betas_gLassoW[[3]] = theta1(tau.grid)
Y10X1.betas_piqrWL[[3]] = Y10X1.betas_gLassoWL[[3]] = theta1(tau.grid)

set.seed(205650)
for (i in 1:nrep) {
        Y10X1 = simulate_qardl(N=500, Ylag = 10, Xlag = 1)
        
        Y10X1.Y = Y10X1$Y
        Y10X1.Z = Y10X1$Z
        
        Y10X1.N = length(Y10X1.Y)
        
        Y10X1.Yvec = Y10X1.Y[11:Y10X1.N]
        Y10X1.Xmat = cbind(Y10X1.Y[1:(Y10X1.N - 1)], Y10X1.Z[1:(Y10X1.N - 1)])
        for (i in 2:10){
                Y10X1.Xmat = cbind(Y10X1.Xmat, c(NA, Y10X1.Y[1:(Y10X1.N - i)]), c(NA, Y10X1.Z[1:(Y10X1.N -i)]))[-1, ] 
        }
        colnames(Y10X1.Xmat) = c("Yt-1", "Xt-1", "Yt-2", "Xt-2", "Yt-3", "Xt-3", 
                                 "Yt-4", "Xt-4", "Yt-5", "Xt-5", "Yt-6", "Xt-6", 
                                 "Yt-7", "Xt-7", "Yt-8", "Xt-8", "Yt-9", "Xt-9", "Yt-10", "Xt-10")
        
        M = length(tau.grid)
        
        qrfit = rq(Y10X1.Yvec ~ Y10X1.Xmat, tau = tau.grid)
        qrfit_coefs = coef(qrfit)
        
        # global coefficients proposed
        phi = phi_generator(L, tau.grid)
        Y10X1.Xmat1 = cbind(rep(1, nrow(Y10X1.Xmat)), Y10X1.Xmat)
        
        piqr = global_qr(taus = tau.grid, phi = phi, X = Y10X1.Xmat1, y = Y10X1.Yvec, lambda = 1, lags = 0, f = "piqr", w = F)
        piqr_coefs = piqr$bhat
        
        piqrW = global_qr(taus = tau.grid, phi = phi, X = Y10X1.Xmat1, y = Y10X1.Yvec, lambda = 1, lags = 0, f = "piqr", w = T)
        piqrW_coefs = piqrW$bhat
        
        piqrWL = global_qr(taus = tau.grid, phi = phi, X = Y10X1.Xmat1, y = Y10X1.Yvec, lambda = 1, lags = 10, f = "piqr", w = T)
        piqrWL_coefs = piqrWL$bhat
        
        gLasso = global_qr(taus = tau.grid, phi = phi, X = Y10X1.Xmat1, y = Y10X1.Yvec, lambda = 1, lags = 0, f = "gLasso", w = F)
        gLasso_coefs = gLasso$bhat
        
        gLassoW = global_qr(taus = tau.grid, phi = phi, X = Y10X1.Xmat1, y = Y10X1.Yvec, lambda = 1, lags = 0, f = "gLasso", w = T)
        gLassoW_coefs = gLassoW$bhat
        
        gLassoWL = global_qr(taus = tau.grid, phi = phi, X = Y10X1.Xmat1, y = Y10X1.Yvec, lambda = 1, lags = 10, f = "gLasso", w = T)
        gLassoWL_coefs = gLassoWL$bhat
        
        if(typeof(piqr_coefs) == "character" | typeof(piqrW_coefs) == "character" | typeof(piqrWL_coefs) == "character" | 
           typeof(gLasso_coefs) == "character" | typeof(gLassoW_coefs) == "character" | typeof(gLassoWL_coefs) == "character") next
        
        
        for(i in 1:nrow(gLasso_coefs)){
                
                Y10X1.betas_qr[[i]] = rbind(Y10X1.betas_qr[[i]], qrfit_coefs[i,])
                Y10X1.betas_piqr[[i]] = rbind(Y10X1.betas_piqr[[i]], piqr_coefs[i,])
                Y10X1.betas_piqrW[[i]] = rbind(Y10X1.betas_piqrW[[i]], piqrW_coefs[i,])
                Y10X1.betas_piqrWL[[i]] = rbind(Y10X1.betas_piqrWL[[i]], piqrWL_coefs[i,])
                Y10X1.betas_gLasso[[i]] = rbind(Y10X1.betas_gLasso[[i]], gLasso_coefs[i,])
                Y10X1.betas_gLassoW[[i]] = rbind(Y10X1.betas_gLassoW[[i]], gLassoW_coefs[i,])
                Y10X1.betas_gLassoWL[[i]] = rbind(Y10X1.betas_gLassoWL[[i]], gLassoWL_coefs[i,])
        }
}

#### SE ####
## difference of estimated coefficients with real ones for each row (nrep) ## 
Y10X1.se_qr = list()
Y10X1.se_piqr = list()
Y10X1.se_piqrW = list()
Y10X1.se_piqrWL = list()
Y10X1.se_gLasso = list()
Y10X1.se_gLassoW = list()
Y10X1.se_gLassoWL = list()


for (i in 1:length(Y10X1.betas_qr)){
        Y10X1.se_qr[[i]] = t(apply(Y10X1.betas_qr[[i]][-1,], 1, function(row) {(row - Y10X1.betas_qr[[i]][1,])^2}))
        Y10X1.se_piqr[[i]] = t(apply(Y10X1.betas_piqr[[i]][-1,], 1, function(row) {(row - Y10X1.betas_piqr[[i]][1,])^2}))
        Y10X1.se_piqrW[[i]] = t(apply(Y10X1.betas_piqrW[[i]][-1,], 1, function(row) {(row - Y10X1.betas_piqrW[[i]][1,])^2}))
        Y10X1.se_piqrWL[[i]] = t(apply(Y10X1.betas_piqrWL[[i]][-1,], 1, function(row) {(row - Y10X1.betas_piqrWL[[i]][1,])^2}))
        Y10X1.se_gLasso[[i]] = t(apply(Y10X1.betas_gLasso[[i]][-1,], 1, function(row) {(row - Y10X1.betas_gLasso[[i]][1,])^2}))
        Y10X1.se_gLassoW[[i]] = t(apply(Y10X1.betas_gLassoW[[i]][-1,], 1, function(row) {(row - Y10X1.betas_gLassoW[[i]][1,])^2}))
        Y10X1.se_gLassoWL[[i]] = t(apply(Y10X1.betas_gLassoWL[[i]][-1,], 1, function(row) {(row - Y10X1.betas_gLassoWL[[i]][1,])^2}))
}

#### Y10X1 MSE ####
Y10X1.mse = list()
for (i in 1:length(Y10X1.betas_qr)){
        Y10X1.mse[[i]] = mean(apply(Y10X1.se_qr[[i]], 2, mean))
        Y10X1.mse[[i]] = c(Y10X1.mse[[i]], mean(apply(Y10X1.se_piqr[[i]], 2, mean)))
        Y10X1.mse[[i]] = c(Y10X1.mse[[i]], mean(apply(Y10X1.se_piqrW[[i]], 2, mean)))
        Y10X1.mse[[i]] = c(Y10X1.mse[[i]], mean(apply(Y10X1.se_piqrWL[[i]], 2, mean)))
        Y10X1.mse[[i]] = c(Y10X1.mse[[i]], mean(apply(Y10X1.se_gLasso[[i]], 2, mean)))
        Y10X1.mse[[i]] = c(Y10X1.mse[[i]], mean(apply(Y10X1.se_gLassoW[[i]], 2, mean)))
        Y10X1.mse[[i]] = c(Y10X1.mse[[i]], mean(apply(Y10X1.se_gLassoWL[[i]], 2, mean)))
        names(Y10X1.mse[[i]]) = c("QR","piqr","piqrW","piqrWL","gLasso","gLassoW","gLassoWL")
}

Y10X1.mse.table = as.data.frame(do.call(rbind, Y10X1.mse))
rownames(Y10X1.mse.table) = c("c", "Yt-1", "Xt-1", "Yt-2", "Xt-2", "Yt-3", "Xt-3", 
                              "Yt-4", "Xt-4", "Yt-5", "Xt-5", "Yt-6", "Xt-6", 
                              "Yt-7", "Xt-7", "Yt-8", "Xt-8", "Yt-9", "Xt-9", "Yt-10", "Xt-10")
Y10X1.mse.table.rel = Y10X1.mse.table[c(1,3,20),]

png("img/Y10X1_mse_table.png", width=350,height=200,bg = "white")
grid.table(format(t(Y10X1.mse.table.rel), scientific = T))
dev.off()

Y10X1.mse_df = data.frame("mse" = Y10X1.mse[[1]], "var" = rep(1, 7), "method" = c(1,2,3,4,5,6,7))
for (i in c(3,20)){
        Y10X1.mse_df = rbind(Y10X1.mse_df, data.frame("mse" = Y10X1.mse[[i]], "var" = rep(i, 7), "method"= c(1,2,3,4,5,6,7)))
}

Y10X1.mse_df$method = as.factor(Y10X1.mse_df$method)
Y10X1.mse_df$var = as.factor(Y10X1.mse_df$var)
levels(Y10X1.mse_df$method) = c("QR", "piqr", "piqrW", "piqrWL", "gLasso", "gLassoW", "gLassoWL")
levels(Y10X1.mse_df$var) = c("c", "Xt-1", "Yt-10")

library(reshape2)
Y10X1.mse.m_1 <- melt(Y10X1.mse_df, id.vars  = c("method", "var"))
ggplot(data = Y10X1.mse.m_1, aes(x=var, y=value)) + geom_point(aes(colour=method), position=position_jitter(w=0.02)) + xlab("Variable") + ylab("MSE")


#### SSE ####
Y10X1.sse = list()
for (i in 1:length(Y10X1.betas_qr)){
        Y10X1.sse[[i]] = apply(Y10X1.se_qr[[i]], 2, mean)
        Y10X1.sse[[i]] = rbind(Y10X1.sse[[i]], apply(Y10X1.se_piqr[[i]], 2, sum))
        Y10X1.sse[[i]] = rbind(Y10X1.sse[[i]], apply(Y10X1.se_piqrW[[i]], 2, sum))
        Y10X1.sse[[i]] = rbind(Y10X1.sse[[i]], apply(Y10X1.se_piqrWL[[i]], 2, sum))
        Y10X1.sse[[i]] = rbind(Y10X1.sse[[i]], apply(Y10X1.se_gLasso[[i]], 2, sum))
        Y10X1.sse[[i]] = rbind(Y10X1.sse[[i]], apply(Y10X1.se_gLassoW[[i]], 2, sum))
        Y10X1.sse[[i]] = rbind(Y10X1.sse[[i]], apply(Y10X1.se_gLassoW[[i]], 2, sum))
}

Y10X1.sse_taus_df <- lapply(seq_along(Y10X1.sse), function(i) {
        df = data.frame(Y10X1.sse[[i]])
        df$var = rep(i, nrow(Y10X1.sse[[i]]))
        df$method = c(1,2,3,4,5,6,7)
        return(df)
})

Y10X1.sse_taus_df$method = as.factor(Y10X1.sse_taus_df$method)
Y10X1.sse_taus_df$var = as.factor(Y10X1.sse_taus_df$var)

Y10X1.sse_df = data.frame("mse" = apply(Y10X1.sse[[1]], 1, sum), "var" = rep(1, 7), "method" = c(1,2,3,4,5,6,7))
for (i in 2:length(Y10X1.sse)){
        Y10X1.sse_df = rbind(Y10X1.sse_df, data.frame("mse" = apply(Y10X1.sse[[i]], 1, sum), "var" = rep(i, 7), "method"= c(1,2,3,4,5,6,7)))
}

Y10X1.sse_df$method = as.factor(Y10X1.sse_df$method)
Y10X1.sse_df$var = as.factor(Y10X1.sse_df$var)
levels(Y10X1.sse_df$method) = c("QR", "piqr", "piqrW", "piqrWL", "gLasso", "gLassoW", "gLassoWL")
levels(Y10X1.sse_df$var) = c("c", "Yt-1", "Xt-1", "Yt-2", "Xt-2", "Yt-3", "Xt-3", 
                             "Yt-4", "Xt-4", "Yt-5", "Xt-5", "Yt-6", "Xt-6", 
                             "Yt-7", "Xt-7", "Yt-8", "Xt-8", "Yt-9", "Xt-9", "Yt-10", "Xt-10")

library(reshape2)
Y10X1.sse.m_1 <- melt(Y10X1.sse_df, id.vars  = c("method", "var"))
ggplot(data = Y10X1.sse.m_1, aes(x=var, y=value)) + geom_point(aes(colour=method), position=position_jitter(w=0.02)) + xlab("Variable") + ylab("SSE")

#### Variable Selection ####

Y10X1.vs = list(length(Y10X1.betas_qr))

for (i in 1:length(Y10X1.betas_qr)){
        Y10X1.vs[[i]] = ifelse(is.null(nrow(Y10X1.betas_qr[[i]][rowSums(Y10X1.betas_qr[[i]]) == 0,])), 0, nrow(Y10X1.betas_qr[[i]][rowSums(Y10X1.betas_qr[[i]]) == 0,])) 
        Y10X1.vs[[i]] = c(Y10X1.vs[[i]], ifelse(is.null(nrow(Y10X1.betas_piqr[[i]][rowSums(Y10X1.betas_piqr[[i]]) == 0,])), 0, nrow(Y10X1.betas_piqr[[i]][rowSums(Y10X1.betas_piqr[[i]]) == 0,])))
        Y10X1.vs[[i]] = c(Y10X1.vs[[i]], ifelse(is.null(nrow(Y10X1.betas_piqrW[[i]][rowSums(Y10X1.betas_piqrW[[i]]) == 0,])), 0, nrow(Y10X1.betas_piqrW[[i]][rowSums(Y10X1.betas_piqrW[[i]]) == 0,])))
        Y10X1.vs[[i]] = c(Y10X1.vs[[i]], ifelse(is.null(nrow(Y10X1.betas_piqrWL[[i]][rowSums(Y10X1.betas_piqrWL[[i]]) == 0,])), 0, nrow(Y10X1.betas_piqrWL[[i]][rowSums(Y10X1.betas_piqrWL[[i]]) == 0,])))
        Y10X1.vs[[i]] = c(Y10X1.vs[[i]], ifelse(is.null(nrow(Y10X1.betas_gLasso[[i]][rowSums(Y10X1.betas_gLasso[[i]]) == 0,])), 0, nrow(Y10X1.betas_gLasso[[i]][rowSums(Y10X1.betas_gLasso[[i]]) == 0,])))
        Y10X1.vs[[i]] = c(Y10X1.vs[[i]], ifelse(is.null(nrow(Y10X1.betas_gLassoW[[i]][rowSums(Y10X1.betas_gLassoW[[i]]) == 0,])), 0, nrow(Y10X1.betas_gLassoW[[i]][rowSums(Y10X1.betas_gLassoW[[i]]) == 0,])))
        Y10X1.vs[[i]] = c(Y10X1.vs[[i]], ifelse(is.null(nrow(Y10X1.betas_gLassoWL[[i]][rowSums(Y10X1.betas_gLassoWL[[i]]) == 0,])), 0, nrow(Y10X1.betas_gLassoWL[[i]][rowSums(Y10X1.betas_gLassoWL[[i]]) == 0,])))
        names(Y10X1.vs[[i]]) = c("QR","piqr","piqrW","piqrWL","gLasso","gLassoW","gLassoWL")
}

Y10X1.vs.table = as.data.frame(do.call(rbind, Y10X1.vs))
rownames(Y10X1.vs.table) = c("c", "Yt-1", "Xt-1", "Yt-2", "Xt-2", "Yt-3", "Xt-3", 
                             "Yt-4", "Xt-4", "Yt-5", "Xt-5", "Yt-6", "Xt-6", 
                             "Yt-7", "Xt-7", "Yt-8", "Xt-8", "Yt-9", "Xt-9", "Yt-10", "Xt-10")
library(gridExtra)

png("Y10X1_vs_table.png", width=480,height=480,bg = "white")
grid.table(Y10X1.vs.table)
dev.off()

Y10X1.vs_df$method = as.factor(Y10X1.vs_df$method)
Y10X1.vs_df$var = as.factor(Y10X1.vs_df$var)
levels(Y10X1.vs_df$method) = c("QR", "piqr", "piqrW", "piqrWL", "gLasso", "gLassoW", "gLassoWL")
levels(Y10X1.vs_df$var) = c("c", "Yt-1", "Xt-1", "Yt-2", "Xt-2", "Yt-3", "Xt-3", 
                            "Yt-4", "Xt-4", "Yt-5", "Xt-5", "Yt-6", "Xt-6", 
                            "Yt-7", "Xt-7", "Yt-8", "Xt-8", "Yt-9", "Xt-9", "Yt-10", "Xt-10")

library(reshape2)
Y10X1.vs.m_1 <- melt(Y10X1.vs_df, id.vars  = c("method", "var"))
ggplot(data = Y10X1.vs.m_1, aes(x=var, y=value)) + geom_point(aes(colour=method), position=position_jitter(w=0.02)) + xlab("Variable") + ylab("# of times var was set to zero")













