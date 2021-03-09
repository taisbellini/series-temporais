# Motivação
n = 1000
b0 = 20
b1 = 1
x = rep(0, n)
y = z = x
dx = 1
dy = 1
dz = 1

set.seed(864)
ex= rnorm(n, 0, dx)
ey = rnorm(n, 0, dy)
ez = rnorm(n, 0, dz)

# 2 random walks
for (i in 2:n){
  x[i] = x[i-1] + ex[i]
  y[i] = y[i-1] + ey[i]
  z[i] = b0 + b1*x[i] + ez[i]
}

x = ts(x)
y = ts(y)
z = ts(z)

ts.plot(x,y, z, col = c("red", "blue", "green"))
# x e z são correlacionadas ao longo do tempo -> raiz unitaria comum
# (z foi gerada a partir de x) 
# cointegração: existe uma combinação linear entre as séries que é estacionária

plot(x,y)
m1 = lm(y~x)
summary(m1)
abline(m1, col = "blue")
ts.plot(residuals(m1))
acf(residuals(m1))
# regressão espúria


#No caso das séries que são correlacionadas - relação de cointegração
# Aqui, a regressão é "verdadeira" (sabemos pois geramos os dados por simulação)
plot(x,z)
m2 = lm(z~x)
summary(m2)
abline(m2, col = "blue")
ts.plot(residuals(m2)) # estacionario e zero -> cointegram
acf(residuals(m2))

# Metodologia para identificar cointegração:
# regredir uma série na outra, avaliar os resíduos
# se zero, cointegram, há uma relação entre elas
# significa que há uma raiz unitária única - caso de x e z
# raiz unitária de z é consequência da raiz unitária de x
# no caso de x e y há duas raizes unitárias independentes

