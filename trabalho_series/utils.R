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
### Save ahat and bhat found

penalty = function(a, f, w){
  if(f == "piqr"){
    penalty = w[1]*sum(abs(a[1,]))
    for(d in 2:nrow(a)){
      penalty = penalty + w[d]*sum(abs(a[d,]))
    }
  }else if (f=="gLasso"){
    penalty = w[1]*cvxr_norm(a[1,])
    for (d in 2:nrow(a)){
      penalty  = penalty + w[d]*cvxr_norm(a[d,])
    }
  }
  return(penalty)
}

# Function to calculate weights based on ahat with no penalty
weights = function(a, lags, Y, X, TAUS, phi){
  w = numeric()
  if (lags != 0){
    lag_penalty = numeric()
    lag_penalty[1] = 1
    for(l in 1:lags){
      lag_penalty = c(lag_penalty, l, l)
    }
  }
  else {
    lag_penalty = rep(1,nrow(a))
  }
  
  objective = R(a, Y, X, TAUS, phi)
  problem = Problem(Minimize(objective))
  result = solve(problem)
  if(result$status == 'solver_error'){
    return("error")   
  }
  ahat = result$getValue(a)
  
  for (j in 1:nrow(ahat)){
    aj = sum(abs(ahat[j,]))
    w[j] = 1/(aj*(exp(-0.5*lag_penalty[j])))
  }
  return(w)
}

global_qr = function(taus = c(0.5), phi = matrix(), X = matrix(), y = c(), lambda = 0, lags=0, f = "piqr", w = F){
  # Get parameters
  M = length(taus)
  L = nrow(phi)
  D = ncol(X)
  N = nrow(X)
  Y = matrix(rep(y,M),N,M)
  TAUS = matrix(rep(taus,N),N,M, byrow = TRUE)
  
  tol = 1e-6
  
  a = Variable(D,L)
  
  if(w){
    we = weights(a, lags, Y, X, TAUS, phi)
  } else {
    we = rep(1, nrow(a))
  }
  
  if(we == "error") return (list("bhat" = "error"))
  
  objective = R(a, Y, X, TAUS, phi) + lambda*penalty(a, f, we)
  problem = Problem(Minimize(objective))
  result = solve(problem)
  if(result$status == 'solver_error') return (list("bhat" = "error"))
  ahat = result$getValue(a)
  ahat_tol = apply(ahat, 1, function(row) {
    if (sum(abs(row))< tol){
      print(sum(abs(row)))
      return(rep(0, length(row)))
    }
    return(row)
  })
  bhat = t(ahat_tol)%*%phi
  
  return(list(
    "ahat" = ahat,
    "bhat" = bhat
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

