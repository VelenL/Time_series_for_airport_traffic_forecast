###############################
# Analysis of Marseille airport
###############################
library(readxl)
Mar <- read_excel("marseille-1982-2016.xls")

traf <- ts(Mar[,3], start=c(1982,1), frequency=12)

##########################
# Analysis of stationarity
##########################
plot.ts(traf)
# non stationarity, increasing variance, non constant trend, potential seasonal effects
acf(ts(traf, frequency=1))
# non stationarity, persistence of significant coefficients, no fast decay to 0
# may have seasonal effects
pacf(ts(traf, frequency=1))
# The main series is non stationary

##################################################################
# 1st transformation: log to remove for increasing variance effect
##################################################################

ltraf <- log(traf)
plot.ts(ltraf)
# reduce increasing variance effect
# still: non constant trend, potential seasonal effect

acf(ts(ltraf, frequency=1), main="Autocorrelogram  log(series)")
# same patterns of non stationarity:
# persistence of significant coefficients, no fast decay to 0

pacf(ts(ltraf, frequency=1), main="Partial autocorrelogram log(series)")

# The log(series) is still non stationary

##############################################################
# 2nd transformation: 1st order difference to remove the trend
##############################################################

dltraf <- diff(ltraf, 1)

plot.ts(dltraf)
# no more increasing variance and no more trend
# still potential seasonal effects

acf(ts(dltraf, frequency=1), main="Autocorrelogram  1st difference of log(series)")
# persistence of significant coefficients with seasonality s=12

pacf(ts(dltraf, frequency=1), main="Partial autocorrelogram 1st difference of log(series)")
# strong significant coefficients around lag 7 and 12

# The series is still non stationary


##################################################################
# 3rd transformation: difference of order 12 to remove seasonality
##################################################################

dltraf_12 <- diff(dltraf, 12)

plot.ts(dltraf_12)
# no more seasonal effects
# We can identify some perturbations in the trafic before 1990

acf(ts(dltraf_12, frequency=1), main="Autocorrelogram  1st difference of log(series) wo seasonality")
# significant coefficients at lags 1, 12

pacf(ts(dltraf_12, frequency=1), main="Partial autocorrelogram 1st difference of log(series) wo seasonality")
# significant coefficients at lags 1, 2, 3, 11, 12, 13, 23, 24


# We can try and fit a multiplicative SARIMA on log(traffic)

#---------------------------------
# Identification of orders p and q
#---------------------------------
# d=D=1, s=12
# p=1, q=1
# P=Q=1 as a starting point
# Maybe then try P=2 to identify effects at lag 24

#------------------------------------------------------------------------
# Estimation of a multiplicative SARIMA(1,1,1)(1,1,1) with seasonality 12 
#------------------------------------------------------------------------

mod1 <- arima(ltraf, c(1,1,1), seasonal=list(order=c(1,1,1), period=12), method='ML')
mod1
#aic = -1056.55

#plot of the fitted value
library("TSA")
fit1 <- fitted(mod1)
fit1
plot.ts(cbind(ltraf,fit1), plot.type='single', col=c('black','red'))

#--------------------
# Validation of mod1
#--------------------

#-------------------------------------------------------------------
# Significance of coefficients: pvalue of Student test on each coeff
#-------------------------------------------------------------------

#-----------------------------------------------------
# Student test to check for the significance of coeffs
#-----------------------------------------------------
# Ho: coeff=0 against H1: coeff#0
# pvalue< 5%, we accept H1

mod1$coef  #value of coefficients
mod1$var.coef #variance of coefficients
tstat <- mod1$coef/sqrt(diag(mod1$var.coef)) # test statistic of the Student test
pvalue <- 2*(1-pnorm(abs(tstat))) #pvalue of the test 
tstat
pvalue

# the pvalues of AR1 and SAR1 are larger than 5%
# So we could remove these coefficients from the model

#-------------------------------------------------------
# Residuals analysis: assumption of Gaussian White Noise
#-------------------------------------------------------
res1 <- mod1$residuals
plot(res1)

# We can identify two extreme value: the residuals are less than -0.2

# General overview of residuals
library(forecast)
checkresiduals(mod1)

#---------------------------------
# Autocorrelation of the residuals
#---------------------------------

acf(ts(res1,frequency=1))
# three significant coeff at lag 6, 11, 21

# Box-Pierce test ou Ljung-Box test
# Ho: no autocorrelation against H1: autocorrelation

Box.test(res1, lag=20, type="Ljung-Box")
# Ho: all correlations up to lag 20 are =0 
# H1: at least one is different from 0
# pvalue= 7% > 5%, we accept HO: no autocorrelation

#-------------------------------
# Normal Distribution assumption
##------------------------------

res_norm <- res1/sqrt(mod1$sigma2) # normalized residuals
summary(res_norm)

# If the residuals follow a Normal distribution, the values of res_norm
# should lie in between -2 and 2, with 95% of chance

plot(res_norm)
abline(a=2,b=0,col="red")
abline(a=-2, b=0, col="red")

out1 <- which(res_norm < -4) # identification number of the outlier
out1 # the outlier corresponds ot observation n°17 and n°50

library("zoo")
index(res_norm)[out1] #date of the outlier
# the outlier occurs at date 1983.333 1986.083 (1 month= 1/12)
traf[out1] # value of the outlier

#QQplot 
# check whether the points are close to the line
qqnorm(res1)
qqline(res1)

#Shapiro test
shapiro.test(res1)
# pvalue < 5% so we reject the normality assumption

# Jarque-Bera test
install.packages('tseries')
library(tseries)
jarque.bera.test(res1)

# constant variance hypothesis for the residuals
sq.res <- mod1$residuals^2
acf(ts(sq.res, frequency = 1), main=" ")
# 3 highly significant coeffs: lag 1, lag 12 and lag 21
vartest <- McLeod.Li.test(mod1, plot=FALSE)
vartest
# The analysis of the 1st pvalue is enough
# pvalue << 5% so we reject constant variance assumption


#-------------------------------------------------------
# Check of the quality of the fit wrt confidence bounds
#-------------------------------------------------------
cb80 <- mod1$sigma2^.5*qnorm(0.9) # confidence bound of security 80%
plot(cbind(ltraf, fit1-cb80, fit1+cb80), plot.type='single', 
     lty=c(1,2,2), xlim=c(2000,2016))

# Proportion of points in the confidence interval
indi <- (ltraf-(fit1-cb80))>0&(fit1+cb80-ltraf)>0
prop <- 100*sum(indi)/length(indi)
prop
# prop = 89%

# prop > 80%, then the fit is considered good

#------------------------
# Validation set approach
#------------------------

data.train <- window(ltraf, start=c(2000,1), end=c(2015,12))
data.test <- window(ltraf, start=c(2016,1), end=c(2016,12))

mod1.train <- arima(data.train, c(0,1,1), seasonal=list(order=c(0,1,1), period=12), method='ML')
pred1.test <- predict(mod1.train, n.ahead=36)

#install.packages("forecast")
library("forecast")
accuracy(pred1.test$pred,data.test)

#                  ME       RMSE        MAE       MPE      MAPE       ACF1 Theil's U
# Test set 0.04964222 0.06032309 0.04964222 0.7649353 0.7649353 -0.2478596 0.5468511

# plot comparing observed values and prediction of the traffic
ts.plot(traf, xlim=c(2014,2019))
lines(2.718^(pred1.test$pred), col="red")
lines(2.718^(pred1.test$pred-1.96*pred1.test$se), col=4, lty=2)
lines(2.718^(pred1.test$pred+1.96*pred1.test$se), col=4, lty=2)

#--------------------------
# Estimation of a 2nd model
#--------------------------

# 1st idea: remove the non significant coeff AR1 and SAR1

# Without AR1 and SAR1
mod2 <- arima(ltraf, c(0,1,1), seasonal=list(order=c(0,1,1), period=12), method='ML')
mod2

#AIC_mod1 = -1056.55
#AIC_mod2 = -1058.51, mod2 is a little bit better sice AIC_mod2 < AIC_mod1

mod2$coef  #value of coefficients
mod2$var.coef #variance of coefficients
tstat <- mod2$coef/sqrt(diag(mod2$var.coef)) # test statistic of the Student test
pvalue <- 2*(1-pnorm(abs(tstat))) #pvalue of the test 
tstat
pvalue
# both MA1 and SMA1 are highly significant!

# Check for residuals autocorrelation for model2
res2 <- mod2$residuals
plot.ts(res2)
acf(ts(res2,frequency=1)) 

# still no significant coefficients before lag 6
Box.test(res2, lag=20, type="Ljung-Box")

# Check for residuals normality

res_norm2 <- res2/sqrt(mod2$sigma2) # normalized residuals
summary(res_norm2)

# If the residuals follow a Normal distribution, the values of res_norm
# should lie in between -2 and 2, with 95% of chance

plot(res_norm2)
abline(a=2,b=0,col="red")
abline(a=-2, b=0, col="red")

#Shapiro test
shapiro.test(res2)

# Jarque-Bera test
jarque.bera.test(res2)

# Check for residuals constant variance
sq.res2 <- (res2)^2
acf(ts(sq.res2, frequency=1))
# 2 highly significant coeffs: lag 12 and lag 21
# The issue of non constant variance remains
vartest2 <- McLeod.Li.test(mod2, plot=FALSE)
vartest2

# Check of the quality of the fit wrt confidence bounds
library("TSA")
fit2 <- fitted(mod2)
fit2
# if (!is.ts(fit2)) fit2 <- ts(fit2)

cb80 <- mod2$sigma2^.5*qnorm(0.9) # confidence bound of security 80%
plot(cbind(ltraf, fit2-cb80, fit2+cb80), plot.type='single', lty=c(1,2,2), xlim=c(2000,2016))

# Proportion of points in the confidence interval
indi <- (ltraf-(fit2-cb80))>0&(fit2+cb80-ltraf)>0
prop <- 100*sum(indi)/length(indi)
prop
#prop_mod1 = 89%
#prop_mod2 = 89%
# no improvement

#------------------------------------------------------
# Let's add dummy variable to correct for the bad fit
#------------------------------------------------------
# Adding an external variable X = fitting a SARIMAX model

out2 <- which(res_norm2 < -4) # identification number of the outlier
out2 # the outlier corresponds ot observation n°50
library("zoo")
index(res_norm2)[out2] #date of the outlier
# the outlier occurs at date 1983.333 and 1986.083 (1 month= 1/12)

# Create Dummy variable at date May 1983 and February 1986
Mar$dum1 <- 0
Mar$dum1[out2] <- 1

mod3 <- arima(ltraf, c(0,1,1), seasonal=list(order=c(0,1,1), period=12), method='ML', xreg=Mar$dum1)
mod3

# AIC_mod3 = -1273.14
# Big improvement wrt the previous models

res3 <- mod3$residuals
par(mfrow=c(2,1))
plot(res2)
plot(res3)
par(mfrow=c(1,1))

# We can identify all value: the residuals are between -0.2 and 0.2

#Shapiro test
shapiro.test(res3)

# Jarque-Bera test
jarque.bera.test(res3)

# Check for residuals constant variance
sq.res3 <- (res3)^2
acf(ts(sq.res3, frequency=1))
vartest3 <- McLeod.Li.test(mod3, plot=FALSE)
vartest3

#------------------------------------------------
# Remark in order to do predictions using model 3
#------------------------------------------------

pred3 <- predict(mod3, n.ahead=36, newxreg=0)
ts.plot(traf, 2.718^pred3$pred,log="y", lty=c(1,3))


#########################
# Change the time period#
#########################
# Because the Pvalues of shapiro test of previous models is less than 5%, 
# which means that models do not have the normality. 
# So we select data from 2001 to 2016 to create another model

sub.traf <- window(traf, start = c(2001, 1), end = c(2016, 12))
plot(sub.traf)
lstraf <- log(sub.traf)
dlstraf <- diff(lstraf, 1)
dlstraf_12 <- diff(dlstraf, 12)
plot.ts(dlstraf_12)
par(mfrow = c(2,1))
acf(ts(dlstraf_12, frequency = 1), main ="sub-period")
# significant coefficients at lags 1，12, 22
pacf(ts(dlstraf_12, frequency = 1), main ="sub-period")
# significant coefficients at lags 1, 12，13
par(mfrow = c(1,1))

#---------------------------------
# Identification of orders p and q
#---------------------------------
# d=D=1, s=12
# p=q=1
# P=Q=1 as a starting point


mod4 <- arima(lstraf, c(1,1,1), 
              seasonal = list(order = c(1,1,1), period = 12),
              method='ML')
mod4
#AIC = -579.28

#plot of the fitted value
library(TSA)
fit4 <- fitted(mod4)
fit4
plot.ts(cbind(lstraf, fit4), 
        plot.type = 'single', 
        col = c('black', 'red'))

#-------------------------------------
# check for the significance of coeffs
#-------------------------------------
mod4$coef  #value of coefficients
mod4$var.coef #variance of coefficients
tstat <- mod4$coef/sqrt(diag(mod4$var.coef)) # test statistic of the Student test
pvalue <- 2*(1-pnorm(abs(tstat))) #pvalue of the test 
tstat
pvalue

# the pvalues of MA1 and SAR1 are much larger than 5%
# So we could remove these coefficients from the model and create a new model

#---------------------------------
# Identification of orders p and q
#---------------------------------
# d=D=1, s=12
# p=1, q=0
# P=0, Q=1 


mod5 <- arima(lstraf, c(1,1,0), 
              seasonal = list(order = c(0,1,1), period = 12),
              method='ML')
mod5
#AIC = -585.22

#plot of the fitted value
library(TSA)
fit5 <- fitted(mod5)
fit5
plot.ts(cbind(lstraf, fit5), 
        plot.type = 'single', 
        col = c('black', 'red'))

#-------------------------------------
# check for the significance of coeffs
#-------------------------------------
mod5$coef  #value of coefficients
mod5$var.coef #variance of coefficients
tstat <- mod5$coef/sqrt(diag(mod5$var.coef)) # test statistic of the Student test
pvalue <- 2*(1-pnorm(abs(tstat))) #pvalue of the test 
tstat
pvalue

# the pvalues of AR1 and SAR1 is less than 5%
# so those coefficients are significant

#--------------------------------------------------------
# Residuals analysis: assumption of Gaussian White Noise
#--------------------------------------------------------

res5 <- mod5$residuals
plot(res5)
# We can identify all residuals are between -0.2 and 0.2

# General overview of residuals
library("forecast")
checkresiduals(mod5)

#---------------------------------
# Autocorrelation of the residuals
#---------------------------------

acf(ts(res5,frequency=1))

# Box-Pierce test ou Ljung-Box test
# Ho: no autocorrelation against H1: autocorrelation

Box.test(res5, lag=20, type="Ljung-Box")
# Ho: all correlations up to lag 20 are =0 
# H1: at least one is different from 0
# pvalue= 87,57% > 5%, we accept HO: no autocorrelation

#-------------------------------
# Normal Distribution assumption
#-------------------------------

res_norm5 <- res5/sqrt(mod5$sigma2) # normalized residuals
summary(res_norm5)

# If the residuals follow a Normal distribution, the values of res_norm
# should lie in between -2 and 2, with 95% of chance

plot(res_norm5)
abline(a=2,b=0,col="red")
abline(a=-2, b=0, col="red")
# QQplot 
# check whether the points are close to the line
qqnorm(res5)
qqline(res5)

# Shapiro test
shapiro.test(res5)
# pvalue=5,74% > 5% so we accept the normality assumption

# Jarque-Bera test
#install.packages('tseries')
library(tseries)
jarque.bera.test(res5)
# pvalue=13,46% > 5% so we accept the normality assumption

# constant variance hypothesis for the residuals
sq.res5 <- mod5$residuals^2
acf(ts(sq.res5, frequency = 1), main=" ")
# 1 highly significant coefficient: lag 1

vartest5 <- McLeod.Li.test(mod5, plot=FALSE)
vartest5
# The analysis of the 1st pvalue is enough
# pvalue=10% > 5% so we accept constant variance assumption


#-------------------------------------------------------
# Check of the quality of the fit wrt confidence bounds
#-------------------------------------------------------
cb80 <- mod5$sigma2^.5*qnorm(0.9) # confidence bound of security 80%
plot(cbind(lstraf, fit5-cb80, fit5+cb80), plot.type='single', 
     lty=c(1,2,2), xlim=c(2000,2016))

# Proportion of points in the confidence interval
indi <- (lstraf-(fit5-cb80))>0&(fit5+cb80-lstraf)>0
prop <- 100*sum(indi)/length(indi)
prop

#Pop 79.17% ~ 80%,this means the fit is considered good



###################
#### Validation ###
###################

data.train5 <- window(lstraf, start=c(2001,1), end=c(2015,12))
data.test5 <- window(lstraf, start=c(2016,1), end=c(2016,12))

mod5.train <- arima(data.train5, c(1,1,0), seasonal=list(order=c(0,1,1), period=12), method='ML')
pred5.test <- predict(mod5.train, n.ahead=36)

accuracy(pred5.test$pred,data.test5)

#                ME       RMSE        MAE       MPE      MAPE       ACF1 Theil's U
#Test set 0.04847248 0.05931827 0.04847248 0.7468494 0.7468494 -0.2703915 0.5398999

# plot comparing observed values and prediction of the traffic
ts.plot(traf, xlim=c(2014,2019))
lines(2.718^(pred5.test$pred), col="red")
lines(2.718^(pred5.test$pred-1.96*pred5.test$se), col=4, lty=2)
lines(2.718^(pred5.test$pred+1.96*pred5.test$se), col=4, lty=2)

# In the end, we got 2 good models : mod3 and mod5
# mod3 has better quality of the fit white noise confidence bounds than mod5 : 
# prop_mod3=89% > prop_mod5 = 79%
# and better AIC_mod3 = -1273.14 <  AIC_mod5= -585.22
# but mod3 doesn't have normality and constant variance.
# However, mod5 has normality and constant variance.