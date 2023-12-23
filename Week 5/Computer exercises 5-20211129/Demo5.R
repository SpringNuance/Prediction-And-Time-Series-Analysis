
#install.packages("lmtest")
#install.packages("car")
#install.packages("forecast")

library(car) # vif() is contained in the package car
library(forecast)
library(lmtest)
alko<-read.table("alcohol.txt",header=T,sep="\t")
LR1C<-ts(alko$LR1C,start=1950)
LQ1CPC<-ts(alko$LQ1CPC,start=1950)
LQTOTALPC<-ts(alko$LQTOTALPC,start=1950)

model1<-lm(LQ1CPC~LR1C+LQTOTALPC)
summary(model1)


#pdf("kx1.pdf",width=8,height=8)
hist(model1$residuals,xlab="Residuals",ylab="Frequency",main=" ")
#dev.off()

#pdf("kx2.pdf",width=8,height=8)
acf(model1$residuals,main="")
#dev.off()

#pdf("k1.pdf",width=8,height=8)
qqnorm(model1$residuals,pch=16)
qqline(model1$residuals,col="red",lwd=2)
#dev.off()

#pdf("k2.pdf",width=8,height=8)
plot(model1$residuals,type="p",ylab="Residuals",xlab="Year",pch=16,xaxt="n")
axis(1,at=seq(from=1,to=32,by=3),labels=seq(from=1950,to=1981,by=3))
abline(0,0)
#dev.off()

fit <- ts(predict(model1),start=1950)

plot(model1$fitted.values)
#does the same as
points(predict(model1))


# Depending on your resolution, the legend in the plots
# might not look nice
# If you uncomment the pdf()-function and the
# dev.off()-function, the corresponding pdf will be generated
# into your current working directory

# pdf("k3.pdf",width=8,height=8)
plot(LQ1CPC,col="red",xlab="Time",ylab="")
lines(fit,col="blue")
legend("topleft", legend=c("LQ1CPC", "fit"),
       col=c("red","blue"),lty=c(1,1),cex=1.8)
# dev.off()

#pdf("k4.pdf",width=8,height=8)
plot(cooks.distance(model1),ylab="Cook's distances",xlab="Year",pch=16,xaxt="n")
axis(1,at=seq(from=1,to=32,by=3),labels=seq(from=1950,to=1981,by=3))
#dev.off()

vif(model1)

zeros <- rep(0,19)
ones <- rep(1,13)
LAW=ts(c(zeros,ones),start=1950)
model2=lm(LQ1CPC~LR1C+LQTOTALPC+LAW)

summary(model2)

#pdf("kx3.pdf",width=8,height=8)
hist(model2$residuals,breaks=seq(from=-0.2,to=0.1,by=0.02),xlab="Residuals",ylab="Frequency",main=" ")
#dev.off()

#pdf("kx4.pdf",width=8,height=8)
acf(model2$residuals,main="")
#dev.off()


#pdf("k5.pdf",width=8,height=8)
qqnorm(model2$residuals,pch=16)
qqline(model2$residuals)
#dev.off()

#pdf("k6.pdf",width=8,height=8)
plot(model2$residuals,type="p",ylab="Residuals",xlab="Year",pch=16,xaxt="n")
axis(1,at=seq(from=1,to=32,by=3),labels=seq(from=1950,to=1981,by=3))
abline(0,0)
#dev.off()

fit2 <- ts(predict(model2),start=1950)
#pdf("k7.pdf",width=8,height=8)
plot(LQ1CPC,col="red",xlab="Time",ylab="")
lines(fit2,col="blue")
legend("topleft", legend=c("LQ1CPC", "fit"),
       col=c("red","blue"),lty=c(1,1),cex=1.8)
#dev.off()

#pdf("k8.pdf",width=8,height=8)
plot(cooks.distance(model2),ylab="Cook's distances",xlab="Index",pch=16,xaxt="n")
axis(1,at=seq(from=1,to=32,by=3),labels=seq(from=1950,to=1981,by=3))
#dev.off()

vif(model2)


DLQ1CPC <- diff(LQ1CPC)
DLR1C <- diff(LR1C)
DLQTOTALPC <- diff(LQTOTALPC)
DLAW <- diff(LAW)


model3<-lm(DLQ1CPC~DLR1C+DLQTOTALPC+DLAW)

summary(model3)


#pdf("kx5.pdf",width=8,height=8)
hist(model3$residuals,breaks=seq(from=-0.1,to=0.1,by=0.02),xlab="Residuals",ylab="Frequency",main=" ")
#dev.off()

#pdf("kx6.pdf",width=8,height=8)
acf(model3$residuals,main="")
#dev.off()


#pdf("k9.pdf",width=8,height=8)
qqnorm(model3$residuals,pch=16)
qqline(model3$residuals)
#dev.off()

#pdf("k10.pdf",width=8,height=8)
plot(model3$residuals,type="p",ylab="Residuals",xlab="Year",pch=16,xaxt="n")
axis(1,at=seq(from=1,to=32,by=3),labels=seq(from=1950,to=1981,by=3))
abline(0,0)
#dev.off()

plot(model3$fitted.values, model3$residuals,type="p",ylab="Residuals", xlab = "Fitted values", pch=16)
abline(0,0)

fit3 <- ts(predict(model3),start=1951)
#pdf("k11.pdf",width=8,height=8)
plot(diff(LQ1CPC),col="red",xlab="Time",ylab="")
lines(fit3,col="blue")
legend("topleft", legend=c("DLQ1CPC", "fit"),
       col=c("red","blue"),lty=c(1,1),cex=1.8)
#dev.off()

#pdf("k12.pdf",width=8,height=8)
plot(cooks.distance(model3),ylab="Cook's distances",xlab="Index",pch=16,xaxt="n")
axis(1,at=seq(from=1,to=32,by=3),labels=seq(from=1950,to=1981,by=3))
#dev.off()

vif(model3)

model3_bg <- rep(NA,27)

# Breusch-Godfrey can be performed up to order: 
# (sample size) - (number of estimated parameters) = 31-4 = 27

for (i in 1:27)
{
  model3_bg[i]= bgtest(model3, order=i)$p.value
}

which(model3_bg > 0.05) 
# Null hypothesis of no autocorrelation accepted with all lags


# Model 4
n = nrow(alko)
model4 <- lm(LQ1CPC[-1]~ LQ1CPC[-n] + LR1C[-1] +LR1C[-n]+
               LQTOTALPC[-1] + LQTOTALPC[-n]+ LAW[-1] +LAW[-n])


summary(model4)

#pdf("kx7.pdf",width=8,height=8)
hist(model4$residuals,breaks=seq(from=-0.1,to=0.1,by=0.02),xlab="Residuals",ylab="Frequency",main=" ")
#dev.off()

#pdf("kx8.pdf",width=8,height=8)
acf(model4$residuals,main="")
#dev.off()

#pdf("k13.pdf",width=8,height=8)
qqnorm(model4$residuals,pch=16)
qqline(model4$residuals)
#dev.off()

#pdf("k15.pdf",width=8,height=8)
plot(model4$residuals,type="p",ylab="Residuals",xlab="Year",pch=16,xaxt="n")
axis(1,at=seq(from=1,to=32,by=3),labels=seq(from=1950,to=1981,by=3))
abline(0,0)
#dev.off()

plot(model3$fitted.values, model3$residuals,type="p",ylab="Residuals", xlab = "Fitted values", pch=16)
abline(0,0)

fit4 <- ts(predict(model4),start=1951)
#pdf("k14.pdf",width=8,height=8)
plot(ts(LQ1CPC,start=1950),col="red",xlab="Time",ylab="",type="l")
lines(fit4,col="blue")
legend("topleft", legend=c("LQ1CPC", "fit"),
       col=c("red","blue"),lty=c(1,1),cex=1.8)
#dev.off()

#pdf("k16.pdf",width=8,height=8)
plot(cooks.distance(model4),ylab="Cook's distances",xlab="Index",pch=16,xaxt="n")
axis(1,at=seq(from=1,to=32,by=3),labels=seq(from=1950,to=1981,by=3))
#dev.off()

model4_bg <- rep(NA,23)

# Breusch-Godfrey can be performed up to order: 
# (sample size) - (number of estimated parameters) = 31 - 8 = 23

for (i in 1:23)
{
  model4_bg[i]= bgtest(model4, order=i)$p.value
}

which(model4_bg > 0.05)
# Null hypothesis of no autocorrelation accepted with all lags

vif(model4)

names(model4)
coef <- model4$coefficients

elasticity1 <- (coef[3]+coef[4])/(1-coef[2])
elasticity2 <- (coef[5]+coef[6])/(1-coef[2])
elasticity3 <- (coef[7]+coef[8])/(1-coef[2])


