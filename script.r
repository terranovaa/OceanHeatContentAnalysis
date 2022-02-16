# Il dataset in formato JSON è stato opportunamente convertito in CSV e il formato della data è stato modificato
table = read.csv("tabella.csv")
table <- data.frame(table[,-1], row.names = table[,1])
# Modifico il valore nullo iniziale, settandolo pari a quello del mese successivo
# per evitare problemi nel modello moltiplicativo di HoltWinters che richiede valori non nulli
table[1,1] = -0.4499674
# Modifico il formato della data

# Esplorazione
T = ts(table)
start(T)
end(T)
plot(T, ylab="Zettajoule", main="Andamento serie")
plot(diff(T), ylab="Serie al netto del trend")
acf(T,50)
acf(diff(T),50)

# Sovrappongo diversi periodi
par(bg = "black")
m_table = matrix(table[179:286,], 12, 9)
ts.plot(m_table, col = heat.colors(9), main="Andamento annuale")
ts.plot(scale(m_table, scale = F), col = heat.colors(9), main="Andamento annuale", xlab="Mese")
lines(rowMeans(scale(m_table, scale = F)), lwd = 3, col = "white")
par(bg = "white")

# Bande di confidenza empiriche
nt.sd = vector("numeric", 12)
for (i in 1:12) {
  nt.sd[i] = sd(m_table[i, ])
}
m_table = scale(m_table, scale=F) 
nt.m = rowMeans(m_table)
plot(nt.m, main="Bande di confidenza empiriche", xlab="Mese", pch = 20, type = "b", ylim = range(c(nt.m - 3 * nt.sd, nt.m + 3 * nt.sd)))
arrows(1:12, nt.m - nt.sd, 1:12, nt.m + nt.sd, length = 0.02, angle = 90, 
       code = 3, col = "green3")
points(nt.m + nt.sd, type = "b", pch = 20, col = "gray")
points(nt.m - nt.sd, type = "b", pch = 20, col = "gray")
boxplot(t(m_table), pch = "*")

# Decomposizione
T = ts(table, start=c(1992,1), end=c(2017,12), frequency = 12)
# Additiva
T.d = decompose(T)
plot(T.d)
# Moltiplicativa
T.dm = decompose(T, "multiplicative")
plot(T.dm)
# Confronto ACF
acf(na.omit(T.d$random),100, main="ACF Decomposizione Additiva")
acf(na.omit(T.dm$random),100, main="ACF Decomposizione Moltiplicativa")

# Confronto residui
T.dr = as.vector(window(T.d$random, c(1992,7), c(2017,6)))
plot(T.dr, pch = 20, main="Residui decomposizione additiva")
T.dmr = as.vector(window(T.dm$random, c(1992,7), c(2017,6)))
T.dmrl = log(T.dmr)
plot(T.dmrl, pch=20, main="Residui decomposizione moltiplicativa")
# Zoom nella zona centrale
plot(T.dmrl[-which(T.dmr > 1 | T.dmr < 0.7)],pch=20) 

var(T.dr)/var(window(T, c(1992, 7), c(2017, 6)))
var(na.omit(T.dmrl))/var(na.omit(window(log(T), c(1992, 7), c(2017, 6))))

hist(T.dr, 20,freq=F, main="Istogramma Decomposizione Additiva")
lines(density(T.dr),col="blue")
lines(sort(T.dr), dnorm(sort(T.dr),mean(T.dr),sd(T.dr)),col="red")
hist(T.dmrl, 20,freq=F, main="Istogramma Decomposizione Moltiplicativa")
lines(density(na.omit(T.dmrl)),col="blue")
lines(sort(na.omit(T.dmrl)), dnorm(sort(na.omit(T.dmrl)),mean(na.omit(T.dmrl)),sd(na.omit(T.dmrl))),col="red")

qqnorm(T.dr, main="Q-Q Plot Decomposizione Additiva")
qqline(T.dr)
qqnorm(T.dmrl, main="Q-Q Plot Decomposizione Moltiplicativa")
qqline(T.dmrl)

shapiro.test(T.dr)
shapiro.test(T.dmrl)

sd(acf(T.dr, plot = F)$acf)
sd(acf(T.dmr, plot = F)$acf)

# stl
plot(stl(T[,1],7))
plot(stl(T[,1],9))

# Analisi residui stl (window = 7)
T.stl = stl(T[,1],7)
T.stlr = T.stl$time.series[,3]
acf(na.omit(T.stlr),100)
plot(T.stlr, pch = 20)
var(T.stlr)/var(window(T, c(1992, 7), c(2017, 6)))
hist(T.stlr, 20,freq=F)
lines(density(T.stlr),col="blue")
lines(sort(T.stlr), dnorm(sort(T.stlr),mean(T.stlr),sd(T.stlr)),col="red")
qqnorm(T.stlr)
qqline(T.stlr)
shapiro.test(T.stlr)

# Analisi residui stl (window = 9)
T.stl = stl(T[,1],9)
T.stlr = T.stl$time.series[,3]
acf(na.omit(T.stlr),100)
plot(T.stlr, pch = 20)
var(T.stlr)/var(window(T, c(1992, 7), c(2017, 6)))
hist(T.stlr, 20,freq=F)
lines(density(T.stlr),col="blue")
lines(sort(T.stlr), dnorm(sort(T.stlr),mean(T.stlr),sd(T.stlr)),col="red")
qqnorm(T.stlr)
qqline(T.stlr)
shapiro.test(T.stlr)

# Sovrappongo trend e stagionalità con quelli di decompose
plot(decompose(T)$trend, col="blue")
lines(stl(T[,1],7)$time.series[,2],col="red")
plot(decompose(T)$seasonal, col="blue")
lines(stl(T[,1],7)$time.series[,1],col="red")

# HoltWinters
T.hw = HoltWinters(T, seasonal="additive")
plot(T.hw)
T.hw$alpha
T.hw$beta
T.hw$gamma

# Determino le condizioni iniziali
x = 1:24
coefficients(lm(T[x]~x))
plot(HoltWinters(T, l.start = -1.06, b.start = -0.04))

# Esploriamo una zona di parametri nell'intorno di quelli forniti dal software
for(alpha in 7:9){
  for(beta in 0:2){
    for(gamma in 8:10){
      plot(HoltWinters(T,alpha=alpha/10,beta=beta/10,gamma=gamma/10,l.start = -1.06, b.start = -0.04),xlab=paste("alpha:",alpha/10," beta:",beta/10," gamma:",gamma/10))
    }
  }
}

# Sovrappongo i risultati a quelli ottenuti dalla decomposizione
ts.plot(T.d$trend, T.hw$fitted[,2],col=c("black","red"), main="Sovrapposizione trend")
legend("bottomright",legend = c("Holt-Winters", "Decomposizione additiva"),col = c("red","black"),lwd=2)
ts.plot(T.d$seasonal, T.hw$fitted[,4],col=c("black","red"), main="Sovrapposizione stagionalità")
legend("bottomright",legend = c("Holt-Winters", "Decomposizione additiva"),col = c("red","black"),lwd=2)

# Previsione
plot(T.hw, predict(T.hw,12), main="Previsione a 12 mesi")
legend("bottomright",legend = c("Holt-Winters", "Serie Originale"),col = c("red","black"),lwd=2)

# Incertezza per via non parametrica
T.hw.r = residuals(T.hw)
plot(T.hw, predict(T.hw, 12))
lines(predict(T.hw, 12) + quantile(T.hw.r, 0.05), col = "green3")
lines(predict(T.hw, 12) + quantile(T.hw.r, 0.95), col = "green3")

# Confronto con modello moltiplicativo 
T.hwm = HoltWinters(T, seasonal="multiplicative")
ts.plot(T, T.hw$fitted[,1], T.hwm$fitted[,1],col=c("black","blue","red"))

# Confronto residui
T.hw.r = residuals(T.hw)
T.hwm.r = residuals(T.hwm)

plot(T.hw.r, type = "p", pch = 20)
plot(T.hwm.r, type = "p", pch = 20)
plot(T.hw$fitted[, 1], T.hw.r, pch = 20)
plot(T.hwm$fitted[, 1], T.hwm.r, pch = 20)

var(T.hw.r)/var(window(T, c(1992, 7), c(2017, 6)))
var(T.hwm.r)/var(window(T, c(1992, 7), c(2017, 6)))

acf(T.hw.r,100)
acf(T.hwm.r,100)

hist(T.hw.r, 20, freq = F)
lines(density(T.hw.r))
lines(sort(T.hw.r), dnorm(sort(T.hw.r), mean(T.hw.r), sd(T.hw.r)), col = "red")
hist(T.hwm.r, 20, freq = F)
lines(density(T.hwm.r))
lines(sort(T.hwm.r), dnorm(sort(T.hwm.r), mean(T.hwm.r), sd(T.hwm.r)), col = "red")

qqnorm(T.hw.r, pch = 20)
qqline(T.hw.r)
qqnorm(T.hwm.r, pch = 20)
qqline(T.hwm.r)

shapiro.test(T.hw.r)
shapiro.test(T.hwm.r)

# Autovalidazione
l = length(T)
res.hw = rep(0, 24)
res.hwm = rep(0, 24)
j = 1
for (i in (l - 24):(l - 1)) {
  T_cv = ts(T[1:i], frequency = 12, start = c(1, 1))
  T.hw = HoltWinters(T_cv, seasonal = "additive")
  T.hwm = HoltWinters(T_cv, seasonal = "multiplicative")
  T.hw.p = predict(T.hw, 1)
  T.hwm.p = predict(T.hwm, 1)
  res.hw[j] = T.hw.p - T[i + 1]
  res.hwm[j] = T.hwm.p - T[i + 1]
  j = j + 1
}
sqrt(mean(res.hw^2))
sqrt(mean(res.hwm^2))
plot(res.hw, type = "b", pch = 20, col = "blue", main="Confronto Capacità di Predizione", ylim=c(-0.6,0.75))
lines(res.hwm, type = "b", pch = 20, col = "green3")
legend("bottomleft",legend = c("Holt-Winters Additivo", "Holt-Winters moltiplicativo"),col = c("blue","green3"),lwd=2)

# Incertezza per via parametrica modello additivo
plot(predict(T.hw, 12), ylim=c(18.2,19.2), col="red", main="Incertezza parametrica previsione")
lines(predict(T.hw, 12) + qnorm(0.05, mean(T.hw.r), sd(T.hw.r)), col = "green")
lines(predict(T.hw, 12) + qnorm(0.95, mean(T.hw.r), sd(T.hw.r)), col = "green")

# Confronto modello con parametri scelti manualmente
T.hw2 = HoltWinters(T, l.start = -1.06, b.start = -0.04, alpha = 0.8, beta = 0.2, gamma = 0.8)
plot(T.hw2, predict(T.hw2,12), main="Previsione a 12 mesi")
legend("bottomright",legend = c("Holt-Winters", "Serie Originale"),col = c("red","black"),lwd=2)

# Confronto residui
T.hw.r = residuals(T.hw)
T.hw2.r = residuals(T.hw2)

plot(T.hw.r, type = "p", pch = 20)
plot(T.hw2.r, type = "p", pch = 20)
plot(T.hw$fitted[, 1], T.hw.r, pch = 20)
plot(T.hw2$fitted[, 1], T.hw2.r, pch = 20)

var(T.hw.r)/var(window(T, c(1992, 7), c(2017, 6)))
var(T.hw2.r)/var(window(T, c(1992, 7), c(2017, 6)))

acf(T.hw.r,100)
acf(T.hw2.r,100)

hist(T.hw.r, 20, freq = F)
lines(density(T.hw.r))
lines(sort(T.hw.r), dnorm(sort(T.hw.r), mean(T.hw.r), sd(T.hw.r)), col = "red")
hist(T.hw2.r, 20, freq = F)
lines(density(T.hw2.r))
lines(sort(T.hw2.r), dnorm(sort(T.hw2.r), mean(T.hw2.r), sd(T.hw2.r)), col = "red")

qqnorm(T.hw.r, pch = 20)
qqline(T.hw.r)
qqnorm(T.hw2.r, pch = 20)
qqline(T.hw2.r)

shapiro.test(T.hw.r)
shapiro.test(T.hw2.r)

# Autovalidazione
l = length(T)
res.hw = rep(0, 24)
res.hw2 = rep(0, 24)
j = 1
for (i in (l - 24):(l - 1)) {
  T_cv = ts(T[1:i], frequency = 12, start = c(1, 1))
  T.hw = HoltWinters(T_cv, l.start = -1.06, b.start = -0.04)
  T.hw2 = HoltWinters(T_cv,l.start = -1.06, b.start = -0.04, alpha = 0.8, beta = 0.2, gamma = 0.8)
  T.hw.p = predict(T.hw, 1)
  T.hw2.p = predict(T.hw2, 1)
  res.hw[j] = T.hw.p - T[i + 1]
  res.hw2[j] = T.hw2.p - T[i + 1]
  j = j + 1
}
sqrt(mean(res.hw^2))
sqrt(mean(res.hw2^2))
plot(res.hw, type = "b", pch = 20, col = "blue", main="Confronto Capacità di Predizione", ylim=c(-0.6,0.75))
lines(res.hw2, type = "b", pch = 20, col = "green3")
legend("bottomleft",legend = c("Holt-Winters Additivo", "Holt-Winters moltiplicativo"),col = c("blue","green3"),lwd=2)

# Metodi regressivi
pacf(T,60, main="Funzione di autocorrelazione parziale")
pacf(diff(T),60, main="Funzione di autocorrelazione parziale al netto del trend")

# Metodo Yule-Walker
T.ar = ar(T)
T.ar
ts.plot(T, T - T.ar$resid, col = c("black", "red"))

# Predizione
T.ar.pt = predict(T.ar, n.ahead = 12, se.fit = FALSE)
plot(T.ar.pt)
ts.plot(T,  T - T.ar$resid, T.ar.pt,col=c("black","red","red"), main="Previsione Metodo Yule-Walker")
legend("bottomright",legend = c("Yule-Walker", "Serie Originale"),col = c("red","black"),lwd=2)

# Incertezza non parametrica previsione
T.ar.r = T.ar$resid[2:312]
y.max = max(T.ar.pt + quantile(T.ar.r, 0.975))
y.min = min(T.ar.pt + quantile(T.ar.r, 0.025))
ts.plot(T.ar.pt, ylim = c(y.min, y.max))
lines(T.ar.pt + quantile(T.ar.r, 0.975), col = "red")
lines(T.ar.pt + quantile(T.ar.r, 0.025), col = "red")

# Analisi residui
plot(T.ar.r, type = "p", pch = 20, main="Residui Yule-Walker")
T.ar.fitted = as.double(na.omit(T - T.ar$resid))
plot(T.ar.fitted, T.ar.r, pch = 20)
l = length(T)
v = var(T[2:l])
var(na.omit(T.ar$resid))/v
acf(T.ar.r,100, main="Autocorrelazione residui Yule-Walker")
pacf(T.ar.r)
hist(T.ar.r, 20, freq = F)
lines(density(T.ar.r))
lines(sort(T.ar.r), dnorm(sort(T.ar.r), mean(T.ar.r), sd(T.ar.r)), col = "red")
qqnorm(T.ar.r, pch = 20)
qqline(T.ar.r)
shapiro.test(T.ar.r)

# Metodo dei minimi quadrati
T.ls = ar(T, method="ols")
T.ls
ts.plot(T, T - T.ls$resid, col = c("black", "red"))

# Predizione
T.ls.pt = predict(T.ls, n.ahead = 12, se.fit = FALSE)
plot(T.ls.pt)
ts.plot(T, T - T.ls$resid, T.ls.pt,col=c("black","red","red"), main="Predizione metodo dei minimi quadrati")

# Incertezza non parametrica previsione
T.ls.r = as.double(na.omit(T.ls$resid))
y.max = max(T.ls.pt + quantile(T.ls.r, 0.975))
y.min = min(T.ls.pt + quantile(T.ls.r, 0.025))
ts.plot(T.ls.pt, ylim = c(y.min, y.max))
lines(T.ls.pt + quantile(T.ls.r, 0.975), col = "green3")
lines(T.ls.pt + quantile(T.ls.r, 0.025), col = "green3")

# Analisi residui
plot(T.ls.r, type = "p", pch = 20)
T.ls.fitted = as.double(na.omit(T - T.ls$resid))
plot(T.ls.fitted, T.ls.r, pch = 20)
var(T.ls.r)/var(T.ls.r + T.ls.fitted)
acf(T.ls.r,100)
pacf(T.ls.r)
hist(T.ls.r, 20, freq = F)
lines(density(T.ls.r))
lines(sort(T.ls.r), dnorm(sort(T.ls.r), mean(T.ls.r), sd(T.ls.r)), col = "red")
qqnorm(T.ls.r, pch = 20)
qqline(T.ls.r)
shapiro.test(T.ls.r)

# Incertezza parametrica previsione
plot(T.ls.pt, ylim=c(18.4,19.2), col="red",  main="Incertezza parametrica previsione")
lines(T.ls.pt + qnorm(0.05, mean(T.ls.r), sd(T.ls.r)), col = "green")
lines(T.ls.pt + qnorm(0.95, mean(T.ls.r), sd(T.ls.r)), col = "green")

# Autovalidazione
res.hw = rep(0, 11)
res.ls = rep(0, 11)
for (i in 1:11) {
  train = window(T, end = c(2017, i))
  test = window(T, start = c(2017, i + 1))
  res.hw[i] = predict(HoltWinters(train), 1)
  res.ls[i] = predict(ar(train, method="ols"), n.ahead = 1, se.fit = F)
}
test = window(T, start = c(2017, 2))
sqrt(mean((test - res.hw)^2))
sqrt(mean((test - res.ls)^2))
plot(res.hw, type = "b", pch = 20, col = "blue", main="Confronto Capacità di Predizione", ylim=c(18,19), xlab="Mese")
lines(res.ls, type = "b", pch = 20, col = "green3")
legend("bottomleft",legend = c("Holt-Winters Additivo", "Metodo dei minimi quadrati"),col = c("blue","green"),lwd=2)

