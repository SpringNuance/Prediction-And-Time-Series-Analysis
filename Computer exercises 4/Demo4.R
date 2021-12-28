
# 4.1

# install.packages("forecast") # Run only once
library(forecast) # Run every time you want to use functions of the forecast package

INTEL <- read.table("INTEL.txt",header=T)
SUNSPOT <- read.table("SUNSPOT.txt",header=T,row.names=1)
SALES <- read.table("SALES.txt",header=T)

Intel_Close <- ts(INTEL$Intel_Close)
Spots <- ts(SUNSPOT,start=1749)
Sales <- ts(SALES$Sales,frequency=12,start=1970)

par(mfrow=c(1,2),mar=c(2.5,2.5,3.5,1.5))
acf(Intel_Close,main="ACF") #decays exponentially?
pacf(Intel_Close,main="PACF") #cuts off after lag 2?

model=Arima(Intel_Close, order=c(2,0,0))
par(mfrow=c(1,2),mar=c(2.5,2.5,3.5,1.5))
acf(model$res,main="ACF")
pacf(model$res,main="PACF")


#Perform Ljung-Box test with lags that are greater than the number of estimated parameters
ljung_box <- c(rep(NA,17))
k <- 2
for(i in 1:17){
  ljung_box[i] <- Box.test(model$res,lag=(i+k),fitdf=k,
                           type="Ljung-Box")$p.value
}

ljung_box

# reset the canvas
dev.off()

# lags for which the null hypothesis is accepted
which(ljung_box > 0.05) + k

fit <- fitted(model)

plot(fit,type="b",col="blue",ylim=c(60,68), ylab="Price",xlab="Time")

lines(Intel_Close,col="red",type="b")

legend(16,68, legend=c("time series", "fit"), col=c("red","blue"),lty=c(2,2),cex=0.8)

model_ver <- Arima(Intel_Close[1:16],order=c(2,0,0))
prediction <- forecast(model_ver,h=4,level=FALSE)$mean
#level=FALSE, omits confidence intervals

plot(Intel_Close,col="red",type="b",ylim=c(60,68),ylab="Price",xlab="Time")
lines(prediction,col="blue",type="b")
legend(16,68, legend=c("time series", "prediction"),col=c("red", "blue"), lty=c(1,1), cex=0.8)


#SPOTS

par(mfrow=c(1,2),mar=c(2.5,2.5,3.5,1.5))
acf(Spots,main="ACF", lag.max=50)
pacf(Spots,main="PACF", lag.max=50)

model2 <- Arima(Spots,order=c(2,0,0))

acf(model2$res,main="AR(2)-model - ACF of residuals", lag.max=50)
pacf(model2$res, main="AR(2)-model - PACF of residuals",lag.max=50)

k <- 2
spots_lb <- rep(NA,47)
for (i in 1:47)
{
  spots_lb[i]=Box.test(model2$res,lag=(i+k),fitdf=k,
                       type="Ljung-Box")$p.value
}
round(spots_lb,3)
which(spots_lb > 0.05)+k

model.auto <- auto.arima(Spots)
model.auto

acf(model.auto$res,lag.max=50,main="ARMA(3,1)-model - ACF of residuals")
pacf(model.auto$res,lag.max=50,main="ARMA(3,1)-model - PACF of residuals") 


k <- 4
spots_lb <- rep(NA,47)
for (i in 1:47)
{
  spots_lb[i]=Box.test(model.auto$res,lag=(i+k),fitdf=k,
                       type="Ljung-Box")$p.value
}
round(spots_lb,3)

which(spots_lb > 0.05)+k

fit.sun <- fitted(model2)
fit.auto <- fitted(model.auto)

dev.off()

plot(fit.sun,type="b",col="blue",ylab="Spots",xlab="Time")

lines(Spots,col="red",type="b")
lines(fit.auto,col="green",type="b")

legend("topleft", legend=c("Spots time series", "AR(2)-fit",
  "ARMA(3,1)-fit"),col=c("red","blue","green"),lty=c(1,1),cex=0.8)

#predict the end of the time series by estimating a model from the beginning of the time series 
model_ver <- Arima(Spots[1:172],order=c(2,0,0))
model_ver1 <- Arima(Spots[1:210],order=c(2,0,0))
# Set the starting year correctly
prediction <- ts(forecast(model_ver,h=43,level=FALSE)$mean,start=1921)
prediction1 <- ts(forecast(model_ver1,h=5,level=FALSE)$mean,start=1959)

par(mfrow=c(1,2),mar=c(2.5,2.5,3.5,1.5))
plot(Spots,col="red",type="l",ylim=c(0,200),
     ylab="Lkm.",xlab="Time",main="43 step prediction")
lines(prediction,col="blue",type="l")

plot(Spots,col="red",type="l",ylim=c(0,200),
     ylab="Lkm.",xlab="Time",main="5 step prediction")
lines(prediction1,col="blue",type="l")

dev.off()

#Below the same is done for the model given by auto.arima

model_vera <- Arima(Spots[1:172],order=c(3,0,1))
model_vera1 <- Arima(Spots[1:210],order=c(3,0,1))
# Set the starting year correctly
predictiona <- ts(forecast(model_vera,h=43,level=FALSE)$mean,start=1921)
predictiona1 <- ts(forecast(model_vera1,h=5,level=FALSE)$mean,start=1959)

par(mfrow=c(1,2),mar=c(2.5,2.5,3.5,1.5))
plot(Spots,col="red",type="l",ylim=c(0,200),
     ylab="Lkm.",xlab="Time",main="43 step prediction")
lines(predictiona,col="blue",type="l")

plot(Spots,col="red",type="l",ylim=c(0,200),
     ylab="Lkm.",xlab="Time",main="5 step prediction")
lines(predictiona1,col="blue",type="l")


#SALES

model_sales <- auto.arima(Sales)
acf(model_sales$res,main="ACF of residuals",lag.max=50)
pacf(model_sales$res,main="PACF of residuals",lag.max=50)

# fitted parameters 2+2+1+2 = 7
k <- 7
sales_lb <- rep(NA,47)
for (i in 1:47)
{
  sales_lb[i]=Box.test(model_sales$res,lag=(i+k),fitdf=k,
                       type="Ljung-Box")$p.value
}
round(sales_lb,3)
which(sales_lb > 0.05)+k # only accepted for lag 9


model_sales2 <- Arima(Sales,order=c(2,1,0),season=c(1,1,0))
acf(model_sales2$res,main="ACF of residuals",lag.max=50)
pacf(model_sales2$res, main="PACF of residuals", lag.max=50)


# fitted parameters 2+1=3
k <- 3
sales_lb2 <- rep(NA,47)
for (i in 1:47)
{
  sales_lb2[i]=Box.test(model_sales2$res,lag=(i+k),fitdf=k,
                       type="Ljung-Box")$p.value
}
round(sales_lb2,3)
which(sales_lb2 > 0.05)+k

dev.off()

fit.sales <- fitted(model_sales2)

plot(fit.sales,type="b",col="blue", ylab="Sales",xlab="Time")

lines(Sales,col="red",type="b")

prediction_sales <- forecast(model_sales2,h=48,level=FALSE)$mean

plot(Sales,col="red",type="b",ylim=c(100,340),
     xlim=c(1970,1987),ylab="Sales",xlab="Time")
lines(prediction_sales,col="blue",type="b")

