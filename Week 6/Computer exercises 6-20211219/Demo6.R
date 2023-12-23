######
# 6.1
######

CONST<- read.table("Const.txt",header=T,sep=",",row.names=1)

const <- ts(CONST,start=1966, frequency=12)
const.stl <- stl(const[,1], s.window="periodic")

# with s.window-parameter it is possible to set the method 
# for estimating the seasonal component

plot(const.stl)

const.filt <- filter(const, c(1,rep(2,11),1)/24 )
trend <- const.stl$time.series[,2]
plot(const, lty=3)
lines(trend, col="blue")
lines(const.filt, lty=2, col="red")
legend("topleft", legend=c("Time series","Filter","STL"),col=c(1,"red","blue"),
       lty=c(3,2,1))


const.diff <- diff(diff(const, lag=12))  

const.diff.stl <- stl(const.diff[,1], s.window="periodic")

plot(const.diff.stl)

plot(const.diff,lty=3,ylab=" ")
trend <- const.diff.stl$time.series[,2]
lines(trend, col="blue")
legend("topleft", legend=c("Time series","STL trend"),col=c(1,"blue"), lty=c(3,1))

######
# 6.2
######

#install.packages("forecast")
library(forecast)

simu <- read.table("arsimulation.txt",header=TRUE)
X <- simu[,1]
Y <- simu[,2]
Z <- simu[,3]

ts.plot(X, main="Time series (1)",ylab="x")
ts.plot(Y, main="Time series (2)",ylab="y")
ts.plot(Z, main="Time series (3)",ylab="z")

X[43] # The peak visible in the first figure

fitX <- Arima(X,order=c(1,0,0), include.mean = FALSE, method="ML")
fitY <- Arima(Y,order=c(1,0,0), include.mean = FALSE, method="ML")
fitZ <- Arima(Z,order=c(1,0,0), include.mean = FALSE, method="ML")

fitX$coef # -0.499
fitY$coef # -0.494
fitZ$coef # -0.492

# We bootstrap the confidence intervals according to steps 1-5,
# given in the lecture slides

set.seed(3141)

# Choose number of repetitions
m <- 5000 # if slow, make this smaller
w <- 10 # Theoretical minimum would be w =2
# However, estimating the AR(1)-parameter from 1 observation
# makes the ML-estimation procedure unstable
# Thus, we set w = 10 

n <- length(X) # The length of the time series is 100
resX <- rep(NA,m) # Initialize an empty vector for the results

for(i in 1:m){
  # Step 1: Select two time points s and u, 0 < s < u <= n, u-s >= w 
  #         uniformly
  
  # Keep choosing the time points, until u-s >= w is satisfied
  u <- 0 # initialize u and s
  s <- 0
  while(u-s < w){
    # The same point cannot be chosen twice
    su <- sample(1:n,2,replace=FALSE) 
    s <- min(su)
    u <- max(su)
  }# Keep repeating until u-s >= w = 10 
  
  # Step 2: Calculate a new parameter vector esimate from the series 
  #         x_s, x_(s+1),..., x_u
  
  # print(i) #uncomment this, if you want to track progress
  
  resX[i] <- Arima(X[s:u],order=c(1,0,0),include.mean = FALSE,
                   method="ML")$coef[1]
} # Step 3: Repeat m-1 times 

# Step 4: order the obtained m estimates from the smallest to the 
#         largest

resXsort <- sort(resX)

# Step 5: Set lower end of the boostrap confidence interval to be 
#         smaller than or equal to the 125th ordered estimate and
#         set the upper end of the bootstrap confidence interval 
#         to be larger than or equal to the 4875th ordered estimate.

confintX <- c(resXsort[125],resXsort[4875])

# Repeat the same for Y and Z

resY <- rep(NA,m) # Initialize an empty vector for the results

for(i in 1:m){
  # Step 1: Select two time points s and u, 0 < s < u <= n, u-s >= w 
  #         uniformly
  
  # Keep choosing the time points, until u-s >= w is satisfied
  u <- 0 # initialize u and s
  s <- 0
  while(u-s < w){
    # The same point cannot be chosen twice
    su <- sample(1:n,2,replace=FALSE) 
    s <- min(su)
    u <- max(su)
  }# Keep repeating until u-s >= 2 
  
  # Step 2: Calculate a new parameter vector esimate from the series 
  #         x_s, x_(s+1),..., x_u
  
  # print(i) #uncomment this, if you want to track progress
  
  resY[i] <- Arima(Y[s:u],order=c(1,0,0),include.mean = FALSE,
                   method="ML")$coef[1]
} # Step 3: Repeat m-1 times 

# Step 4: Order the obtained m estimates from the smallest to the largest

resYsort <- sort(resY)

# Step 5: Set lower end of the boostrap confidence interval to be smaller than
#         or equal to the 125th ordered estimate and set the upper end of the 
#         bootstrap confidence interval to be larger than or equal to the
#         4875th ordered estimate.

confintY <- c(resYsort[125],resYsort[4875])

resZ <- rep(NA,m) # Initialize an empty vector for the results

for(i in 1:m){
  # Step 1: Select two time points s and u, 0 < s < u <= n, u-s >= w 
  #         uniformly
  
  # Keep choosing the time points, until u-s >= w is satisfied
  u <- 0 # initialize u and s
  s <- 0
  while(u-s < w){
    # The same point cannot be chosen twice
    su <- sample(1:n,2,replace=FALSE) 
    s <- min(su)
    u <- max(su)
  }# Keep repeating until u-s >= 2 
  
  # Step 2: Calculate a new parameter vector esimate from the series 
  #         x_s, x_(s+1),..., x_u
  
  # print(i) #uncomment this, if you want to track progress
  
  resZ[i] <- Arima(Z[s:u],order=c(1,0,0),include.mean = FALSE,method="ML")$coef[1]
} # Step 3: Repeat m-1 times 

# Step 4: order the obtained m estimates from the smallest to the largest
# Note that Step 4 is not compulsory in R 

resZsort <- sort(resZ)

# Step 5: Set lower end of the boostrap confidence interval to be smaller than
#         or equal to the 125th ordered estimate and set the upper end of the 
#         bootstrap confidence interval to be larger than or equal to the
#         4875th ordered estimate.

confintZ <- c(resZsort[125],resZsort[4875])

round(confintX,2) # -0.94 and -0.17
round(confintY,2) # -0.75 and -0.26
round(confintZ,2) # -0.75 and -0.23

######
# 6.3
######

#install.packages("KFAS")
library(KFAS)
alko <-ts(read.table("alcoholdeaths.txt"),start=1969)

a1 <- c(0,0) # Initial guess for mu and nu
Zt <- matrix(c(1, 0), 1, 2) 
Ht <- matrix(NA) 
Tt <- matrix(c(1, 0, 1, 1), 2, 2) 
Rt <- matrix(c(1, 0), 2, 1) 
Qt <- matrix(NA) 
P1 <- matrix(0, 2, 2)
P1inf <- diag(2)

# -1 sets that no constant is estimated in the model
model_gaussian <- SSModel(alko~ -1 + SSMcustom(a1=a1,Z = Zt, T = Tt, 
                                               R = Rt, Q = Qt, 
                                               P1 = P1,P1inf = P1inf),H = Ht)

fit_gaussian <- fitSSM(model_gaussian, inits = c(0, 0))
# above inits-parameter is related to the estimation procedure of the
# unknown variances

fit_gaussian$model$Q
fit_gaussian$model$H
  
out_gaussian <- KFS(fit_gaussian$model)

plot(alko)
lines(out_gaussian$a[,1],col="red")

out_gaussian$a[,1]
out_gaussian$a[,2]