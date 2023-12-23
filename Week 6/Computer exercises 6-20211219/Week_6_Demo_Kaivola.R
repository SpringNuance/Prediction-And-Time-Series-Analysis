setwd('~/Prediction/Week_6')

# 6.1

# a)

CONST<- read.table("Const.txt",header=T,sep=",",row.names=1)

const <- ts(CONST$V1,start=1966, frequency=12)

ts_plotter <- function(ts, lag_max) {
  par(mfrow = c(2,2))
  plot(ts, main = 'Run Plot')
  acf(ts, main = 'Autocorrelation', lag.max = lag_max)
  pacf(ts, main = 'Partial Autocorrelation', lag.max = lag_max)
  spec <- spectrum(ts, plot = F)
  plot(spec$freq, spec$spec, main = 'Spectrum' ,
       xlab = 'Frequency', ylab = 'Energy density')
}

ts_plotter(const,50)

# Some sort of seasonal component modulates the shape of the time series
# The expected value also seems to evolve as a function of time, though
# the trend is not monotonic (i.e. increasing or decreasing over the range)
# Perhaps some kind of polynomial trend would suffice for us.
# In any case, the time series does not seem to be stationary.

# We can think of the time series as being composed of three components
# x_t = m_t + s_t + e_t
# where m_t is the deterministic trend, s_t the deterministic seasonal and e_t the stochastic
# component.

# b)

?stl
# The basic structure of this decomposition method:
# (1) Extract the seasonal component by taking the mean of all seasonal sub-series
# (when s.window = periodic).
# (2) Remove the extracted seasonal component from the series and use local polynomial regression
# i.e. fit low-degree polynomials using a small window of the data to extract the trend.
# (LOESS = "Locally Estimated Scatter-Plot Smoothing")
# (3) Remove the extracted trend from the seasonal component and add to the trend component
# (4) Repeat the above for a couple of iterations. What remains is the stochastic component.

stl_decomp <- stl(const, s.window = "per")
plot(stl_decomp)
# The bar on the right-hand side provides a reference range to evaluate the relative contribution
# of each component. The smaller the component (relative to the remainder), the higher its contribution is
# relative to noise. So in this case the seasonal component is not as substantial in explaining the 
# the variation of the time series as the trend component.

# c) The linear filter given in the exercise is the Savitzky-Golay filter
# https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter
# using a window of 13 observations and some low-degree polynomial.
# The coefficients are obtained by minimizing a least squares loss function.
?filter
const_filt <- filter(const, c(1,rep(2,11),1) * 1/24)
stl_trend <- stl_decomp$time.series[,'trend']

plot(const,lty=3)
lines(stl_trend, col="blue")
lines(const_filt, lty=2, col="red")
legend("topleft", legend=c("Time series","STL","Filter"),
       col=c(1,"blue",'red'), lty=c(3,1,2))

# The overall behavior of the trend estimates are similar but the filter trend is somewhat rougher
# than the STL trend. Perhaps caused by the usage of higher order polynomials whereas STL uses
# first order polynomials to perform the smoothing.

# d)

const_diff <- diff(diff(const, lag = 12))

const_diff_stl <- stl(const_diff, s.window = 'per')
plot(const_diff_stl)

# After differencing, the contribution of the seasonal and trend components have become minor
# relative to the stochastic component. Such behavior indicates that the differenced time
# series could be stationary. This kind of descriptive analysis is probably the most reasonable
# way to use these model-less smoothing methods since it is very difficult to infer anything
# about the precision or significance of the components.

# 6.2

# a)

library(forecast)

simu <- read.table("arsimulation.txt",header=TRUE)

data_6_2 <- apply(simu,2,ts, simplify = F)

lapply(data_6_2, function(x) ts_plotter(x,50))

# The kurtosis of the distribution decreases when degrees of freedom increases
# Thus large deviations are substantially more often observed when dof = 1.

# b)

models <- lapply(data_6_2, function(x) arima(x, order = c(1,0,0), include.mean = F, method = 'ML'))
sapply(models, coef) # The estimated AR(1) parameters are quite close to the true parameter value
# in all cases

# c)

set.seed(3141)

interval_boot_helper <- function(data, min_len) {
  start = sample(1:(length(data) - min_len + 1), size = 1)
  end = sample((start + min_len - 1):length(data), size = 1)
  
  coef(arima(data[start:end], order = c(1,0,0), include.mean = F, method = 'ML'))
  
}

boot_replicator <- function(n_boot ,data, min_len) {
  replicate(n_boot,interval_boot_helper(data, min_len))
}

n_boot <- 5000
alpha <- 0.05

X_boot <- boot_replicator(n_boot, data_6_2$X, 10)
Y_boot <- boot_replicator(n_boot, data_6_2$Y, 10)
Z_boot <- boot_replicator(n_boot, data_6_2$Z, 10)
boot_samples <- cbind(X_boot, Y_boot, Z_boot)

apply(boot_samples, 2, function(x) quantile(x, probs = c(alpha/2, 1 - alpha/2) ))

#boot_samples <- sapply(data_6_2, function(x) boot_replicator(n_boot, x, 10))

# d)
# An important observation with regards to the residual processes is that the moments of the t-distribution
# with n degree of freedom are defined up to order n-1. Thus process (1) has residuals with an undefined
# expected value and variance. Since the variance has to be well-defined (i.e. finite), the process (1)
# is not weakly stationary.

# The remaining processes have finite variance. Since the AR(1) coefficient is |.| < 1 and the residuals
# are mutually independent w.r.t. themselves and past elements of the time series, they are weakly stationary.

# 6.3
library(KFAS)

# Linear State Space Model (LSSM)
# A flexible generalization of the SARIMA models (i.e. all SARIMA models can be represented as
# linear state-space models). LSSM processes do not need to be stationarity.
# Basic definition:
# Hidden state x_t, observation y_t
# x_{t+1} = Tx_t + Rn_t, n_t ~ N(0,Q) # State Process
# y_t = Zx_t + e_t, e_t ~ N(0,H) # Observation Process
# x_1 ~ N(x_1, P_1)
# Noise terms and the initial state distribution are mutually independent.

# In the lecture slides, we are concerned with inferring properties of the hidden
# states x_t given the model parameters, mainly its expected value and variance. 

# There are various different ways of incorporating the observations to the inference.
# The two main ways are: Filtering Distribution P(x_t | y_t, y_{t-1}, ..., y_1)
# Smoothing Distribution P(x_t | y_T, y_{T-1}, ..., y_1) where T is the final observed time point.
# We also consider the Predictive Distribution P(x_t | y_{t-1}, ..., y_1)

# Which one of them is useful to you depends on what the use case is. Smoothing is probably useful
# if you are performing retrospective analysis i.e. all data has already been collected.
# Filtering is more useful in online cases where observations are received sequentially, especially
# if there are control inputs to the system.

# Kalman Filter is an algorithm which can be used to calculate the distributions of the hidden state
# given above. The recursive update equations show that we only need to keep track of the current
# state prediction and its variance which is again useful for online applications.
# Check the Wikipedia article, you can find derivations and other interesting material about the
# update equations. 
# Roughly speaking, the Filter incorporates the current observation to the
# running state estimate using the optimal Kalman Gain which minimizes the mean square error
# related to the state estimate.

# The estimation of the model parameters is not really discussed in this course.
# The approach is maximum likelihood but it is complicated by the fact that
# the model involves hidden states. The traditional way of dealing with this
# is the Expectation-Maximization (EM) algorithm which is discussed at least in
# the course Advanced Probabilistic Methods.

alko <-ts(read.table("alcoholdeaths.txt"),start=1969)

# Our model:
# mu_{t+1} = mu_t + v + n_t, n_t ~ N(0,Q)
# y_t = mu_t + e_t, e_t ~ N(0,H)

# Construct the hidden state from mu_t and v.

a1 <- c(0,0) # Mean vector for the initial state
Zt <- matrix(c(1, 0), 1, 2) # Emission matrix
Ht <- matrix(NA) # Observation noise variance
Tt <- matrix(c(1, 0, 1, 1), 2, 2) # Transition matrix
Rt <- matrix(c(1, 0), 2, 1) # State noise matrix
Qt <- matrix(NA) # State noise variance
P1 <- matrix(0, 2, 2) # Non-diffuse state covariance matrix
P1inf <- diag(2) # Diffuse state indicator matrix

?SSModel # Specify the model

model_gaussian <- SSModel(alko~-1+SSMcustom(a1=a1,Z=Zt,T=Tt,R = Rt,Q=Qt,
                                            P1=P1,P1inf=P1inf),H=Ht)

?fitSSM # ML estimation of the variance parameters

fit_gaussian <- fitSSM(model_gaussian, inits = c(0, 0))

# ML estimates
fit_gaussian$model$Q # State noise variance
fit_gaussian$model$H # Observation noise variance

?KFS # Filtering, smoothing, whatever distribution you want to know.

hidden_state_infer <- KFS(fit_gaussian$model)
plot(alko, ylab = 'Deaths per 100 000')
lines(hidden_state_infer$a[,1],col="red") # One-step-ahead predictions of state
lines(hidden_state_infer$alphahat[,1],col='blue') # Smoothed states
legend('topleft', legend = c('Original', 'Prediction', 'Smoothed'), 
       col = c('black', 'red','blue'), lty = 1)
pred_states <- hidden_state_infer$a[,1]
std_pred_states <- sqrt(hidden_state_infer$P[1,1,])

arrows(x0=time(alko), y0=pred_states - std_pred_states,
       x1=time(alko), y1=pred_states + std_pred_states, code=3, angle=90, length=0.1)

# For the drift, we can check how the predictions evolve
hidden_state_infer$a[,2] # The last value is 0.84
hidden_state_infer$alphahat[,2] # Equal everywhere as expected
