#------------------------------------
# Stat 370
# Name: Daneille Gagne and David Cao
# Instructor: Cristina Anton 
# Project
#------------------------------------
library(readxl)
require(TSA)

lag.plot1=function(data1,max.lag=1,corr=TRUE,smooth=FALSE){ 
  name1=paste(deparse(substitute(data1)),"(t-",sep="")
  name2=paste(deparse(substitute(data1)),"(t)",sep="")
  data1=as.ts(data1)
  max.lag=as.integer(max.lag)
  prow=ceiling(sqrt(max.lag))
  pcol=ceiling(max.lag/prow)
  a=acf(data1,max.lag,plot=FALSE)$acf[-1]
  par(mfrow=c(prow,pcol), mar=c(2.5, 4, 2.5, 1), cex.main=1.1, font.main=1)
  for(h in 1:max.lag){                       
    plot(lag(data1,-h), data1, xy.labels=FALSE, main=paste(name1,h,")",sep=""), ylab=name2, xlab="") 
    if (smooth==TRUE) 
      lines(lowess(ts.intersect(lag(data1,-h),data1)[,1],
                   ts.intersect(lag(data1,-h),data1)[,2]), col="red")
    if (corr==TRUE)
      legend("topright", legend=round(a[h], digits=2), text.col ="blue", bg="white", x.intersp=0)
  }
}

diag1=function(xdata, fitit)
{
  n = length(xdata)
  k = length(fitit$coef)
  BIC = -2*fitit$loglik+k*log(n)
  
  AIC = -2*fitit$loglik+2*k
  AICc = AIC+2*(1+k)*(k+2)/(n-k-2)
  list(AIC=AIC, AICc=AICc, BIC=BIC)
}

erro=function(xdata, preddata)
{
  #xdata are the real values
  # preddata are the predicted values
  # xdata and preddata should have the same dimension
  n = length(preddata)
  m=length(xdata)
  e=xdata-preddata# the error
  MSE=mean(e*e)
  MPE=mean(e/xdata)
  MAE=mean(abs(e))
  MAPE=mean(abs(e/xdata))   
  list(MPE=MPE, MSE=MSE, MAE=MAE, MAPE=MAPE)
}

#------------------------------------

Edmonton <- read_excel("Edmonton.xlsx")
EdmontonData <- Edmonton[1:195,4] # test data set
EdmontonFull <- Edmonton[,4] # full data set
ESnowFull <- ts(EdmontonFull, start = 1937, frequency = 3)

ESnow <- ts(EdmontonData, start=1937,frequency=3)
win.graph(width=10, height=10,pointsize=12)
plot.ts(ESnow, type="o", pch = 8, cex = 0.75, main="Snow fall(cm) for the winter Quarter Starting in 1937")

summary(ESnow)
win.graph(width=10, height=10,pointsize=12)
par(mfrow=c(2,1))        # set up the graphics 
hist(ESnow, prob=TRUE, 12)   # histogram    
lines(density(ESnow))     # smooth it - ?density for details 
qqnorm(ESnow)             # normal Q-Q plot  
qqline(ESnow)             # add a line    

win.graph(width=10, height=10,pointsize=12)
par(mfrow=c(2,1))        # set up the graphics 
acf(ESnow)
pacf(ESnow)
eacf(ESnow)

win.graph(width=15, height=15,pointsize=12)
lag.plot1(ESnow,12,corr=TRUE,smooth=TRUE)

win.graph(width=15, height=15,pointsize=12)
bct<-BoxCox.ar(ESnow + 0.01,method='yw')
bct$ci
bct$mle

DESnow <- diff(ESnow,3)

win.graph(width=10, height=10,pointsize=12)
par(mfrow=c(2,1))        # set up the graphics 
acf(DESnow)
pacf(DESnow)

fit1=arima(ESnow,order=c(0,0,0),seasonal=list(order=c(0,1,1), period=3))
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


#Model with purely seasonal MA(1) and outlier
fit2=arima(ESnow,order=c(0,0,0),seasonal=list(order=c(0,1,1), period=3), xreg=as.numeric (seq (ESnow) == 136))
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

FullEdmontonSnow = diff(ESnowFull,3)


#Model with (0,0,1)x(0,1,1)3
fit3=arima(ESnow, order=c(0,0,1), seasonal = list(order=c(0,1,1), period=3), xreg=as.numeric(seq (ESnow)== 136))
fit3$coef #(parameter estimates)
sqrt(diag(fit3$var.coef))# standard errors
fit3$sigma2 # noise variance

res3=residuals(fit3)
win.graph(width=9.5, height=7,pointsize=12)
par(mfrow=c(2,1))        # set up the graphics 
hist(res3, prob=TRUE, 12)   # histogram    
lines(density(res3))     # smooth it - ?density for details 
qqnorm(res3)             # normal Q-Q plot  
qqline(res3)             # add a line 
win.graph(width=9.5, height=12,pointsize=12)
par(mfrow=c(2,1))
acf(res3, lag=40)
pacf(res3, lag=40)

win.graph(width=9.5, height=7,pointsize=12)
par(mfrow=c(2,1)) 
plot(fit3,n.ahead=12,type='b', pch=16,newxreg=0)
plot(fit3,n.ahead=12,type='b', pch=16, xlim=c(2000,2005), newxreg=0)
points(time(ESnowFull), ESnowFull, col = 2, pch = 16, type = 'b')

#AR coefficient is not significant, and predictions do not improve in terms of including all true values

###################  Predictions
win.graph(width=9.5, height=7,pointsize=12)
par(mfrow=c(2,1)) 
plot(fit2,n.ahead=12,type='b', pch=16,newxreg=0)
plot(fit2,n.ahead=12,type='b', pch=16, xlim=c(2000,2005), newxreg=0)
points(time(ESnowFull), ESnowFull, col = 2, pch = 16, type = 'b')

#Prediction errors 

xESnow=ESnowFull[196:205]
prem1=predict(fit2, n.ahead = 10, newxreg = 0)
preddata1=(prem1$pred)
preddata1
erro(xESnow, preddata1)


diag1(ESnowFull,fit2)


#------------------------------------------------
#------------------------------------------------
John <- read_excel("St.John MonthlySnow.xlsx")
StJohnData <- John[,4]
JohnSnow <- ts(StJohnData, start=1942,frequency=3)
win.graph(width=10, height=10,pointsize=12)
plot.ts(JohnSnow, type="o", pch = 8, cex = 0.75, main="Snow fall(cm) for the St.John winter Quarter Starting in 1942")

JohnTest=John[1:201, 4]
JohnSnowtest=ts(JohnTest, start=1942, frequency=3)

summary(JohnSnowtest)
win.graph(width=10, height=10,pointsize=12)
par(mfrow=c(2,1))        # set up the graphics 
hist(JohnSnowtest, prob=TRUE, 12)   # histogram    
lines(density(JohnSnowtest))     # smooth it - ?density for details 
qqnorm(JohnSnowtest)             # normal Q-Q plot  
qqline(JohnSnowtest)             # add a line  

win.graph(width=15, height=15,pointsize=12)
bct<-BoxCox.ar(JohnSnowtest+0.01,method='yw')
bct$ci
bct$mle

#squareroot transformation applied

sqrtJohnSnowtest <- (JohnSnowtest)^(0.5)

win.graph(width=14, height=10,pointsize=12)
plot.ts(sqrtJohnSnowtest, type="o", pch = 8, cex = 0.75, main="Snow fall(cm) for the St.John winter Quarter Starting in 1942, sqrt transformed")

win.graph(width=14, height=10,pointsize=12)
par(mfrow=c(2,1))        # set up the graphics 
hist(sqrtJohnSnowtest, prob=TRUE, 12)   # histogram    
lines(density(sqrtJohnSnowtest))     # smooth it - ?density for details 
qqnorm(sqrtJohnSnowtest)             # normal Q-Q plot  
qqline(sqrtJohnSnowtest)             # add a line 

win.graph(width=10, height=10,pointsize=12)
par(mfrow=c(2,1))        # set up the graphics 
acf(sqrtJohnSnowtest)
pacf(sqrtJohnSnowtest)
eacf(sqrtJohnSnowtest)

win.graph(width=15, height=15,pointsize=12)
lag.plot1(sqrtJohnSnowtest,12,corr=TRUE,smooth=TRUE)

#-----------------------

diff.sqrt.stjt <- diff(sqrtJohnSnowtest, 3)
win.graph(width=10, height=10,pointsize=12)
par(mfrow=c(2,1))        # set up the graphics 
acf(diff.sqrt.stjt)
pacf(diff.sqrt.stjt)

fit1.1=arima(sqrtJohnSnowtest,order=c(0,0,0),seasonal=list(order=c(0,1,1), period=3))
fit1.1$coef #(parameter estimates)
sqrt(diag(fit1.1$var.coef))# standard errors
fit1.1$sigma2 # noise variance

res1.1=residuals(fit1.1)
win.graph(width=9.5, height=7,pointsize=12)
par(mfrow=c(2,1))        # set up the graphics 
hist(res1.1, prob=TRUE, 12)   # histogram    
lines(density(res1.1))     # smooth it - ?density for details 
qqnorm(res1.1)             # normal Q-Q plot  
qqline(res1.1)             # add a line 
win.graph(width=9.5, height=12,pointsize=12)
par(mfrow=c(2,1))
acf(res1.1, lag=40)
pacf(res1.1, lag=40)

win.graph(width=9.5, height=7,pointsize=12)
tsdiag(fit1.1,gof=20,omit.initial=T)
detectAO(fit1.1) # detect outliers 
detectAO(fit1.1, robust=F)
detectIO(fit1.1)

diag1(sqrtJohnSnowtest,fit1.1)

#------------------------
#Predicting on JohnSnowtest data using the first fitted model
#-----------------------
fit1.1=arima(sqrtJohnSnowtest,order=c(0,0,0),seasonal=list(order=c(0,1,1), period=3))

sqrtJohnSnow=(JohnSnow^(0.5))

win.graph(width=9.5, height=7, pointsize=12)
par(mfrow=c(2,1))

plot(fit1.1, n.ahead=10, type='b', pch=16, main="Sqrt Transformed Series")
plot(fit1.1, n.ahead=10, type='b', pch=16, xlim=c(2008,2012))
grid(NULL, NA, lwd=2)
points(time(sqrtJohnSnow), sqrtJohnSnow, col="red", pch=16)


win.graph(width=9.5, height=7, pointsize = 12)
par(mfrow=c(2,1))
plot.ts(JohnSnow, col=1, type='o', xlim=c(1942, 2012), pch=16, main="Initial Series")
lines((prem3$pred)^2, type="p", col=2, pch=16)
plot.ts(JohnSnow, col=1, type="o", xlim=c(2005, 2012), pch=16)
lines((prem3$pred)^2, type="p", col=2, pch=16)
grid(NULL, NA, lwd=2)


#---- prediction errors-----
xsqrtJohnSnow=sqrtJohnSnow[202:211]
xdata1=JohnSnow[202:211]

prem2=predict(fit1.1, n.ahead = 10)
preddata2=(prem2$pred)
preddata2


prem3=predict(fit1.1, n.ahead=10)
preddata3=(prem3$pred)^2
preddata3

erro(xsqrtJohnSnow, preddata2)
erro(xdata1, preddata3)


prem2=predict(fit2, n.ahead=10)

preddata2=(prem2$pred)^2
erro(xdata1, preddata2)

