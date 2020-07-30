John <- read_excel("St.John MonthlySnow.xlsx")
StJohnData <- John[,4]
JohnSnow <- ts(StJohnData, start=1942,frequency=3)
win.graph(width=10, height=10,pointsize=12)
plot.ts(JohnSnow, type="o", pch = 8, cex = 0.75, main="Snow fall(cm) for the St.John winter Quarter Starting in 1942")

win.graph(width=10, height=10,pointsize=12)
par(mfrow=c(2,1))        # set up the graphics 
acf(JohnSnow)
pacf(JohnSnow)
eacf(JohnSnow)

#-----------------------

fit1=arima(JohnSnow,order=c(0,0,0),seasonal=list(order=c(0,1,1), period=3))
fit1$coef #(parameter estimates)
sqrt(diag(fit1$var.coef))# standard errors
fit1$sigma2 # noise variance

res1=residuals(fit1)
win.graph(width=9.5, height=7,pointsize=12)
par(mfrow=c(2,1))        # set up the graphics 
hist(res1, prob=TRUE, 12)   # histogram    
lines(density(res1))     # smooth it - ?density for details 
qqnorm(res1)             # normal Q-Q plot  
qqline(res1)             # add a line 
win.graph(width=9.5, height=12,pointsize=12)
par(mfrow=c(2,1))
acf(res1, lag=40)
pacf(res1, lag=40)

win.graph(width=9.5, height=7,pointsize=12)
tsdiag(fit1,gof=20,omit.initial=T)
detectAO(fit1) # detect outliers 
detectAO(fit1, robust=F)
detectIO(fit1)

diag1=function(xdata, fitit)
{
  n = length(xdata)
  k = length(fitit$coef)
  BIC = -2*fitit$loglik+k*log(n)
  
  AIC = -2*fitit$loglik+2*k
  AICc = AIC+2*(1+k)*(k+2)/(n-k-2)
  list(AIC=AIC, AICc=AICc, BIC=BIC)
}

diag1(JohnSnow,fit1)

#---------------------------
fit2=arima(JohnSnow,order=c(0,0,1),seasonal=list(order=c(0,1,1), period=3))
fit2$coef #(parameter estimates)
sqrt(diag(fit2$var.coef))# standard errors
fit2$sigma2 # noise variance

res2=residuals(fit2)
win.graph(width=9.5, height=7,pointsize=12)
par(mfrow=c(2,1))        # set up the graphics 
hist(res2, prob=TRUE, 12)   # histogram    
lines(density(res2))     # smooth it - ?density for details 
qqnorm(res2)             # normal Q-Q plot  
qqline(res2)             # add a line 
win.graph(width=9.5, height=12,pointsize=12)
par(mfrow=c(2,1))
acf(res2, lag=40)
pacf(res2, lag=40)

win.graph(width=9.5, height=7,pointsize=12)
tsdiag(fit2,gof=20,omit.initial=T)
detectAO(fit2) # detect outliers 
detectAO(fit2, robust=F)
detectIO(fit2)

diag1=function(xdata, fitit)
{
  n = length(xdata)
  k = length(fitit$coef)
  BIC = -2*fitit$loglik+k*log(n)
  
  AIC = -2*fitit$loglik+2*k
  AICc = AIC+2*(1+k)*(k+2)/(n-k-2)
  list(AIC=AIC, AICc=AICc, BIC=BIC)
}

diag1(JohnSnow,fit2)
