#### Utils trabalho series ####
#### Loading required packages ####

if(!require(quantreg)) install.packages('quantreg')
if(!require(CVXR)) install.packages('CVXR')

library(quantreg)
library(CVXR)


#### Objective function ####

R = function(a, Y, X, TAUS, phi){
  B = a%*%phi # candidate optimizer
  R.eval = (2*TAUS-1)*(Y-X%*%B) + abs(Y-X%*%B)
  return(sum(R.eval)) # mean of the preceding objective functions
}


#### Simulation function #### 
### Run one experiment estimating using the traditional qr estimator and the proposed one
### Save ahat and bhat found

penalty = function(a, lambda=0){
  w = weights(a)
  a_norm = numeric()
  for (d in 1:nrow(a)){
    a_norm[d] = value(cvxr_norm(a[d,]))
  }
  pen = lambda*sum(w*a_norm)
  return(pen)
}

#TODO test the norm for aj here as well
weights = function(a){
  w = numeric()
  for (j in 1:nrow(a)){
    aj = value(cvxr_norm(a[j,]))
    w[j] = 1/aj*(exp(-0.5))
  }
  return(w)
}

#todo - optimize via optim
global_qr = function(taus = c(0.5), phi = matrix(), X = matrix(), y = c(), lambda = 0){
  # Get parameters
  M = length(taus)
  L = nrow(phi)
  D = ncol(X)
  N = nrow(X)
  
  Y = matrix(rep(y,M),N,M)
  TAUS = matrix(rep(taus,N),N,M, byrow = TRUE)
  
  a = Variable(D,L)
  objective = R(a, Y, X, TAUS, phi) + penalty(a, lambda)
  problem = Problem(Minimize(objective))
  result = solve(problem)
  ahat = result$getValue(a)
  bhat = ahat%*%phi
  
  return(list(
    "ahat" = ahat,
    "bhat" = bhat
  ))
}

bhatQR = function(taus = c(0.5), X = matrix(), y = c()){
  bhat = rq(y~X[,-1],taus)$coef
  return(list("bhat" = bhat))
} 

R_univ = function(a = a, taus = taus, y = y){ # objective function
  N = length(y)
  M = length(taus)
  
  b = phi%*%a # candidate optimizer
  foo = sapply(1:M, function(m){
    argeval = y-b[m]
    # objective function corresponding to quantile level = taus[m]:
    sum((2*taus[m]-1)*(argeval) + abs(argeval))/N
    # sum((taus[m]-1)*(argeval) + argeval*pnorm(argeval,sd=h) + h^2*dnorm(argeval,sd=h))/N
  }
  )
  return(sum(foo)/M) # mean of the preceding objective functions
}

global_qr_uni = function(phi = phi){
  L = nrow(phi)
  ahat = optim(rep(1,L),R)$par
  bhat = phi%*%ahat
  return (list(
    "ahat" = ahat,
    "bhat" = bhat
  ))
}





#### MC function ####
### Run a MC experiment estimating using the traditional qr estimator and the proposed one
### Save MSE for each iteration and evaluate results
### Save Bhat found for each

# nrepl: number of monte carlo replications
# N: sample size
# D: covariables 
# M: tau grid size

run_mc = function(nrepl = 200, N = 100, taus = c(0.5), phi = matrix(), beta = c(1), lambda = 0, sd_error = 0.1){
  
  # Get parameters
  M = length(taus)
  D = length(beta)
  L = nrow(phi)
  
  # Initialize response variables
  Bhat_list = list()
  BhatQR_list = list()
  mse_Bhat_list = list()
  mse_BhatQR_list = list()
  
  for (i in 1:D){
    Bhat_list[[i]] = matrix(0,nrepl,M)
    BhatQR_list[[i]] = matrix(0,nrepl,M)
    mse_Bhat_list[[i]] = matrix(0,nrepl,M)
    mse_BhatQR_list[[i]] = matrix(0,nrepl,M)
  }
  
  # Run simulation 
  for (i in 1:nrepl){
    
    X = matrix(runif(N*D),N,D)
    X[,1] = 1
    Xbar = colMeans(X)
    
    y = X%*%beta + rnorm(N,sd=sd_error)
    
    Y = matrix(rep(y,M),N,M)
    TAUS = matrix(rep(taus,N),N,M, byrow = TRUE)
    
    bhatQR = rq(y~X[,-1],taus)$coef
    
    a = Variable(D,L)
    penalty = function(a) sum(abs(a))
    objective = R(a, Y, X, TAUS) + lambda*penalty(a)
    problem = Problem(Minimize(objective))
    result = solve(problem)
    if(result$status == 'solver_error') next
    ahat = result$getValue(a)
    bhat = ahat%*%phi
    
    for(k in 1:D){
      Bhat_list[[k]][i,] = bhat[k,]
      BhatQR_list[[k]][i,] = bhatQR[k,]
      mse_Bhat_list[[k]][i,] = (beta[k] - bhat[k,])^2
      mse_BhatQR_list[[k]][i,] = (beta[k] - bhatQR[k,])^2
    }
  }
  
  return(list(
    "bhat" = Bhat_list,
    "bhatQR" = BhatQR_list,
    "mse_bhat" = mse_Bhat_list,
    "mse_bhatQR" = mse_BhatQR_list
  ))
}


#### Multiplot ####
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#### Helper functions ####

# Generate the phi matrix based on the taus grid and number of basis vectors
phi_generator = function(L, taus){
  M = length(taus)
  phi = matrix(0,L,M) # matrix to store basis vectors
  colnames(phi) = as.character(taus)
  
  # finite dimensional analogue of Legendre polynomials
  phi[1,] = rep(1,M)
  phi[1,] = phi[1,]/sqrt(sum(phi[1,]^2))
  phi[2,] = lm(taus~0+phi[1,])$res
  phi[2,] = phi[2,]/sqrt(sum(phi[2,]^2))
  if (L>=3) {
    for (ell in 3:L){
      phi[ell,] = lm(taus^(ell-1)~t(phi[2:ell,]))$res
      phi[ell,] = phi[ell,]/sqrt(sum(phi[ell,]^2))
    }
  }
  return(phi)
}

