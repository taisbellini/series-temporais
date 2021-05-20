# This code simulates a sample path of size T from a simple QARDL(1,1) model:
# Q(𝜏| ℱ[t-1]) = α₀(𝜏) + α₁(𝜏)*Y[t-1] + θ₁(𝜏)*Z[t-1]

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
dbplot = function(a,b, title="", x_lab="tau.grid", y_lab=""){
 tau = seq(from=0,to=1,by=.01)
 plot(tau, dbeta(tau,a,b), type = "l", main = title, xlab=x_lab, ylab=y_lab)
}

# Parameters for the Beta quantile functions v10, v01 and v00
{a10 = 1; b10 = 3; a01 = 5; b01 = 1; a00 = 3; b00 = 1}

# Conditional quantile function of Y[t] given Y[t-1] = 1 and X[t-1]=0
v10 = function (tau) qbeta(tau, a10,b10)

# Conditional quantile function of Y[t] given Y[t-1] = 0 and X[t-1]=1
v01 = function(tau) qbeta(tau,a01,b01)

# Conditional quantile function of Y[t] given Y[t-1] = 0 and X[t-1]=0
v00 = function(tau) qbeta(tau,a00,b00)

# Conditional quantile function of Y[t] given Y[t-1] = 1 and X[t-1]=1
v11 = function(tau) v10(tau) + v01(tau) - v00(tau)

# Conditional density of Y[t] given Y[t-1] = 1 and X[t-1]=1
q11 = function(tau) {
1/dbeta(v10(tau),a10,b10) + 1/dbeta(v01(tau),a01,b01) - 1/dbeta(v00(tau),a00,b00)
}

# Functional parameters for the quantile regression equation
alpha0 = v00
alpha1 = function(tau) v10(tau) - v00(tau)
theta1 = function(tau) v01(tau) - v00(tau)

# Quantile function of Y given ℱ[t-1]
Q = function(tau,Y.current,Z.current){
 alpha0(tau) + alpha1(tau)*Y.current + theta1(tau)*Z.current
}

# Simulating the sample paths:

simulate_qardl = function(N = 10001){
        
        # Arbitrary starting point (Y0,Z0)
        Y0 = runif(1)
        Z0 = runif(1)
        
        Y = Z = numeric()
        Z.current = Z0
        Y.current = Y0

        for (t in 1:N){
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
        return(list(
                "Y" = Y,
                "Z" = Z
        ))
}
