y = AirPassengers
ts.plot(y)
# observa-se pelo grafico que a variância aumenta ao longo do tempo, sendo não-estacionario


ly = log(y)
ts.plot(ly)

acf(ly, lag=48)

# tomamos a diferença para que fique estacionario
dly = diff(ly)
acf(dly, lag=48)
pacf(dly, lag=48)

# uma opção de estimação: modelo ARIMA
m1 = arima(ly, order = c(0,1,0), seasonal=list(order = c(0,1,0), frequency=12))
r1 = residuals(m1)
pacf(r1)

# pelos resultados, identifico que possivelmente AR(1) no componente sazonal se ajusta melhor
m2 = arima(ly, order = c(0,1,0), seasonal=list(order = c(1,1,0), frequency=12))
r2 = residuals(m2)
acf(r2, lag = 48)

#MA(2)
m3 = arima(ly, order = c(0,1,0), seasonal=list(order = c(1,1,3), frequency=12))
r3 = residuals(m3)
acf(r3, lag = 48)
m3


m4 = arima(ly, order = c(0,1,0), seasonal=list(order = c(1,1,3), frequency=12), fixed = c(NA,NA,0,NA))
r4 = residuals(m4)
pacf(r4, lag = 48)
m4
