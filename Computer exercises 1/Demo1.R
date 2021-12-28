

# 1.1
# (a)
# setwd("path")

emis <- read.table("emissions.txt",header=T,sep="\t")
# Note that the imported data is of the data type called "data.frame"
# data.frame behaves usually as matrix, but some matrix operations are
# not available for data frames, e.g. matrix multiplication

emis.matrix <- as.matrix(emis)
# More about different data types on Introduction to R-programming course

set.seed(123)
# set.seed fixes the random number generator in a certain "position"

colnames(emis) #returns the names of the variables
emis$NOx
emis[,"NOx"]
emis[,2] 
# All 3 commands above give the observations of NOx

# With View() it is possible to examine the imported data
View(emis)

# summary() can be used to extract useful information from different objects in R
summary(emis)

# Note that histogram is not the same thing as barplot
hist(emis[,"NOx"],main="NOx", xlab="kg/cm^2")
hist(emis[,"Humidity"],main="Humidity",xlab="unit")
hist(emis[,"Pressure"],main="Pressure",xlab="unit")
hist(emis$Temp,main="Temperature",xlab="unit")

# The number of bins in a histogram can be changed. For example
b1 <- c(seq(min(emis$NOx),max(emis$NOx),length.out=10))
# The command above divides the interval between
# the smallest and largest value of NOx with 10 points.

hist(emis[,"NOx"],breaks=b1)

# Try also
hist(emis[,"NOx"],main="NOx", xlab=expression(kg/cm^2))

# By the histograms, the variables do not seem to be normally distributed

# (b), (c)

fit1 <- lm(NOx~Humidity+Temp+Pressure, data=emis) #linear model = lm

summary(fit1) 

summary(fit1)$r.squared # coefficient of determination

# (d)
# Remember to run set.seed() every time before running the permutation-loop to
# obtain identical results

# The results may vary depending on the operating 
# system and the version of R
set.seed(123)

k <- 2000 # Number of repetitions
y.mean <- mean(emis$NOx) 

# Mathematical formulations of the formulas below can be found from the lecture slides
SST <- sum((emis$NOx-y.mean)^2)
SSE <- sum((fit1$res)^2)
Rsquare1 <- 1 - SSE/SST

# Create an empty matrix, NA = not available.
# The number of rows is k  (number of permutations)
# The number of columns is 3 (number of explanatory variables)
perm<- matrix(NA,nrow=k,ncol=3)

for(i in 1:k){
  # Create 3 temporary variables to permute one explanatory variable at a time
  # Could be written more generally for m explanatory variables
  tmp1 <- emis
  tmp2 <- emis
  tmp3 <- emis
  
  # By default sample() draws a sample of the size of the input without replacement
  # try help(sample) 
  tmp1$Pressure <- sample(tmp1$Pressure)
  tmp2$Humidity <- sample(tmp2$Humidity)
  tmp3$Temp <- sample(tmp3$Temp)
  
  # LS-estimates
  tmpfit1 <- lm(NOx ~ Humidity+Temp+Pressure, data=tmp1)
  tmpfit2 <- lm(NOx ~ Humidity+Temp+Pressure, data=tmp2)
  tmpfit3 <- lm(NOx ~ Humidity+Temp+Pressure, data=tmp3)

  perm[i,1] <- summary(tmpfit1)$r.squared
  perm[i,2] <- summary(tmpfit2)$r.squared
  perm[i,3] <- summary(tmpfit3)$r.squared
}
# Note that it is not feasible to calculate all possible permutations,
# since 20! is of the magnitude 10^18:
factorial(20)

# The perm-matrix now contains three columns with k observations, such that the elements of the
# perm-matrix are values for the coefficient of determination (R-squared values).
# The first column is formed such that the variables Humidity and Temp are exactly as in the original
# data set in every k iterations. The variable Pressure is given a random permutation in every iteration
# and the R-squared value is calculated again in every iteration.

# The second column is formed such that the variables Pressure and Temp are left untouched and
# the variable Humidity is permutated in every iteration.

# Similarly, the third column is formed such that the variables Pressure and Humidity are left untouched
# and the variable Temp is permutated in every iteration.

# Choose the significance level alpha = 5% and order the observations:
preord <- sort(perm[,1])
humord <- sort(perm[,2])
tempord <- sort(perm[,3])

# (note that in R, the sorting is not necessary for the next step)

# The null hypothesis (H_0) of beta_j = 0 is rejected if the 
# original R-squared is larger than the calculated 95th percentile
Rsquare1 > quantile(preord,0.95)  # H_0 Accepted (at a significance level of 5%)
Rsquare1 > quantile(humord,0.95)  # H_0 Rejected
Rsquare1 > quantile(tempord,0.95) # H_0 Rejected

# Alternatively, you get the same results from:
pre <- 1-sum(perm[,1] < Rsquare1)/k #0.4635, H_0 Accepted
hum <- 1-sum(perm[,2] < Rsquare1)/k #0, H_0 Rejected
temp <- 1-sum(perm[,3] < Rsquare1)/k #0.0115, H_0 Rejected


# (e)
fit2 <- lm(NOx~Humidity+Temp, data=emis)
summary(fit2)

# (f)
colnames(emis) #The variables Humidity and Temp are in the third and fourth columns
tmp <-as.matrix(emis[,c(3,4)])

# dim, nrow and ncol commands can be used to find out the dimensions of data
n <- nrow(emis)

Intercept <- rep(1,n) #vector of ones

emis2 <- cbind(Intercept,tmp)
p <- ncol(emis2)

res <- fit2$residuals

s2 <- sum(res^2)/(n-p)

vari <- s2*diag( solve( t(emis2)%*%emis2) )
stdev <- sqrt(vari)

# (g)

confint(fit2,level=0.95)

# (h)
coef <- fit2$coef
up <- coef + stdev * qt(0.975,n-p)
down <- coef - stdev * qt(0.975,n-p)

# (i)
# The logic in bootstrapping is mainly the same as in permutation test

# Below are steps 1-5 from the lecture slides

k <- 2000
bootmat <- matrix(NA,nrow=k,ncol=3) # Empty matrix for the results

y <- emis$NOx
X <- emis2

set.seed(123)

for(i in 1:(k-1)){
  # 1. Select n data points randomly with replacement from the 
  # original observations
  # (Note that in the permutation test the data points are selected
  # without replacement)
  ind <- sample(1:n,replace = TRUE) # Indices for the random rows
  
  Xtmp <- X[ind,]
  ytmp <- y[ind]
  
  # 2. Calculate a new parameter vector estimate from the new sample
  btmp <- solve(t(Xtmp)%*%Xtmp)%*%t(Xtmp)%*%ytmp
  
  # 3. Save the results to bootmat and repeat k-1 times
  bootmat[i,] <- t(btmp)
}

# 4. Include the original estimate and order the obtained estimates 
# from the smallest to the largest
boriginal <- solve(t(X)%*%X)%*%t(X)%*%y

bootmat[k,] <- t(boriginal)

boot1 <- sort(bootmat[,1])
boot2 <- sort(bootmat[,2])
boot3 <- sort(bootmat[,3])

# 5. Set alpha=5% and set the lower end of the confidence interval
# to be smaller or equal to the [0.025*2000] ordered estimate.
# Set the upper end to be larger or equal to the [0.975*2000] ordered estimate.

qconst <- quantile(boot1, probs = c(0.025,0.975))
qhum <- quantile(boot2, probs = c(0.025,0.975))
qtemp <- quantile(boot3, probs = c(0.025,0.975))

# Note that step 4 is not necessary in R, 
# you would get the same results by inputting the
# columns of bootmat into the quantile-function. 

# Here, we have wanted to emphasize the steps of the lecture slides and thus 
# have included step 4 separately
  
# Plot the histograms of the bootstrap estimates
hist(bootmat[,1])
hist(bootmat[,2])
hist(bootmat[,3])

