# install.packages("ggplot2") in console to install library


library(ggplot2)
y <- 1
y = 1

# Prefer = to <- when using functions due to scoping rules
?mean
mean(x=2)
mean(x <- 2)

# vectors and matrices
?c # help(c) to see its documentations
z <- c(1,2,3,4)
x <- c(5,6,7,8)
?cbind
m <- cbind(z,x)

#indexing (Starts at 1)
?'['
z[1]
m[1,1]
m[1,]
m[,'x']
z[z>2]

?read.table
emission_data_og <- read.table('emissions.txt', header=T, sep = '\t')
emission_data_og$Obsno

# Let's drop the first column
emission_data <- emission_data_og[,-1]

# Set seed to get deterministic results
set.seed(123)

# a)
?hist
?par

par(mfrow = c(2,2))

hist(emission_data$NOx)
hist(emission_data$Humidity)
hist(emission_data$Temp)
hist(emission_data$Pressure)

?apply #very useful
apply(emission_data, 2, hist) # not that pretty

?lapply # same but with lists
var_names <- colnames(emission_data)
var_names
?colnames
sapply(var_names, function(x) hist(emission_data[,x], main = x, xlab = 'Unit'))

# all vars seem to be non-normal (Skewness, general asymmetry)

# b)
# y = w^Tx
?lm #fitting linear models
# ~ symbol means the prev is response variable and subsequent vars are explanatory vars
model_b <- lm(NOx ~ Humidity + Temp + Pressure, data = emission_data)
lm(NOx ~ .,data = emission_data) #shorthand for above

# Useful function all around
?summary
summary(model_b)

# Elements of model summary
# Residuals: certain stats of estimated residuals

# Coefficients:
# estimate: the least squares (LS) estimates for explanatory vars and the intercept
# standard error: Estimate for the standard deviations of the LS estimator
# t value: Parametric t-test for the hypothesis: H_0: beta_i = 0 H_1: beta_i != 0
# Is reasonable if the residuals are normally distributed
# Pr(>|t|): Associated two-tailed p-value

# Residual standard error: Estimate for the residuals standard deviation
# Multiple R-Squared: Coefficient of Determination
# Adjusted R-squared: COeff.of Deter. adjusted for model degrees of freedom (num of explanatory vars)
# F-stats: A parametric test for the entire model against a base model that only has the intercept
# Test is: H_(): beta_i = 0 for all i (excluding the intercept) H_1: beta_i != 0 for some i
# is reasonable if residuals are normally distributed

# c)
summary(model_b)$r.squared

# Calculate R^2 by hand 
# total sum of squares
SST <- sum((emission_data$NOx - mean(emission_data$NOx))^2)
# error sum of squares 
SSE <- sum(resid(model_b)^2)
# R^2 (R^2 is coefficient of determination) = 1 - SSE/SST 
1 - SSE/SST # The proportion of variation that can be explained by the model
# model sum of squares
SSM <- sum((fitted(model_b) - mean(emission_data$NOx))^2)
SSE + SSM # SST = SSE + SSM
SST

# alternatively
r_squared_og <- cor(fitted(model_b), emission_data$NOx)^2
r_squared_og

# d) 
# Permutation test: H_0: beta_i = 0, H_1: beta_i != 0
# High level idea: Permute the values of a given explanatory variable. This breaks the dependency structure
# Using many different permutations, we can estimate the distribution of R^2 under the null hypothesis
# Then the p-value is the prob. of observing at least as good of an R^2 as that computed from the original data
# under the null hypothesis

# X: Data Matrix, y: Response
fit_helper_d <- function(X,y,perm_var) {
  X[,perm_var] <- sample(X[,perm_var])
  # LS estimate
  beta <- solve(t(X) %*% X) %*% t(X) %*% y
  # fitted values
  y_hat <- X %*% beta
  # R^2
  cor(y_hat,y)^2
}

perm_replicator_d <- function(n_perm, X, y, var_name) {
  replicate(n_perm, fit_helper_d(X,y,var_name))
}

n <- nrow(emission_data)
Intercept <- rep(1,n)
X <- cbind(as.matrix(emission_data[,-1]), Intercept)
y <- emission_data$NOx

n_perm <- 2000
alpha <- 0.05

expl_var <- var_names[-1]
r_squares <- sapply(expl_var, function(name) perm_replicator_d(n_perm, X, y, name))

p_values_perm <- apply(r_squares, 2, function(x) sum(x > r_squared_og)/length(x))
p_values_perm
p_values_perm < alpha # The null is rejected for humidity and Temperature

# e)
model_e <- lm(NOx ~ . - Pressure, data = emission_data)
summary(model_e)
summary(model_b)

# f) 
X_f <- X[, -3]

# Residual variance
res_e <- resid(model_e)
res_var <- 1/(n - ncol(X_f)) * sum(res_e^2)
sqrt(res_var)

# Standard deviation of LS estimators
Ls_std <- sqrt(diag(res_var * solve(t(X_f) %*% X_f)))
Ls_std

Ls_std <- Ls_std[c(3,1,2)]

# g) 

confint(model_e, level = 1 - alpha)

#) h)

beta <- coef(model_e)
?qt
beta - qt(1-alpha/2, n - ncol(X_f)) * Ls_std
beta + qt(1-alpha/2, n - ncol(X_f)) * Ls_std

# i) 
# Bootstrapping. The high level idea: Sample with replacement from the original data set and calculate bootstrap estimate
# of the statistic of interest (in the case LS estimate). Using many bootstrap samples, we can estimate the true sampling
# distribution of the statistic

fit_helper_i <- function(X,y) {
  # Sample with replacement
  inds <- sample(1:nrow(X), replace = T)
  X <- X[inds,]
  y <- y[inds]
  # Bootstrap LS estimate
  solve((t(X) %*% X)) %*% t(X) %*% y
}

n_boot <- 2000

boot_samples <- replicate(n_boot, fit_helper_i(X_f,y), simplify = 'matrix')
boot_samples <- cbind(boot_samples, beta[c(2,3,1)])
apply(boot_samples, 1, function(x) quantile(x, probs = c(alpha/2, 1- alpha/2, 1-alpha/2)))

par(mfrow = c(2,2))
var_names <- rownames(boot_samples)
sapply(var_names, function(x) hist(boot_samples[x,], main =x, xlab = 'Unit'))

# For HW
?plot
?abline

