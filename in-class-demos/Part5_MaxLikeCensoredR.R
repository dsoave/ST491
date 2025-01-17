rm(list=ls()); options(scipen=999)
wdata = read.table("https://raw.githubusercontent.com/dsoave/Datasets/refs/heads/main/Weibull.data2.txt")
dim(wdata)
head(wdata)
summary(wdata)
Time = wdata$Time; Uncensored = wdata$Uncensored # Avoiding the attach() function

# Find MLE numerically

mloglike = function(theta,t,delta)
     { # Minus log likelihood function
     alpha = theta[1]; lambda = theta[2]
     # logf and logS will be of length n
     logf = log(alpha)+log(lambda)+(alpha-1)*log(lambda*t) + -(lambda*t)^alpha
     logS = -(lambda*t)^alpha
     value = -sum(logf*delta) - sum(logS*(1-delta))
     return(value)
     } # End of function mloglike

# Testing
mloglike(c(3,0.2),t=Time,delta=Uncensored)

yes = Time[Uncensored==1]; no = Time[Uncensored==0]
-sum(dweibull(yes,shape=3, scale=5,log=TRUE)) - sum(pweibull(no, shape=3, scale=5, lower.tail = FALSE, log.p = TRUE))



################################################################################
# Find MLE
################################################################################

startvals = c(1,1/2) # I tried a few values

search1 = optim(par=startvals, fn=mloglike, t=Time,delta=Uncensored,
                     hessian=TRUE, lower=c(0,0), method='L-BFGS-B')


search1


# If both eigenvalues of the Hessian are positive, minus LL is concave up.
H = search1$hessian
eigen(H)$values

alphahat = search1$par[1]; lambdahat = search1$par[2]

# These data were simulated, so I know the true parameter values.
truealpha = 3; truelambda = 1/5
# Compare:
c(alphahat,truealpha)
c(lambdahat,truelambda)


################################################################################
# Calculate the estimated asymptotic covariance matrix of the MLEs.
################################################################################

Vhat = solve(H); Vhat # Solve returns the inverse.

################################################################################
# Point estimate and confidence interval for the median
# Median = log(2)^(1/alpha) / lambda
################################################################################

# Point estimate of median
medhat = 1/lambdahat * log(2)^(1/alphahat); medhat
# Compare the truth
truemedian = log(2)^(1/truealpha) / truelambda; truemedian

median(Time) # Sample median is way off, because it ignores censoring



################################################################################
# Plot hazard function h(t) = alpha*lambda^alpha * t^(alpha-1)
################################################################################

x = seq(from=0,to=5,length=101)
esthazard = alphahat*lambdahat^alphahat * x^(alphahat-1)
truehazard = truealpha*truelambda^truealpha * x^(truealpha-1)

plot(x,truehazard,type='l',xlab='Time',ylab='Hazard',
 main='Hazard Function for the Weibull Data')
lines(x,esthazard,lty=2)
# Annotate the plot (Make the legend)
x1 = c(0.5,1.5); y1 = c(0.5,0.5)
lines(x1,y1,lty=1)
text(2,0.5,'True h(t)')
x2 = c(0.5,1.5); y2 = c(0.45,0.45)
lines(x2,y2,lty=2)
text(2.2,0.45,'Estimated h(t)')


################################################################################
# What if the (Weibull) model is wrong? Try bowl shaped hazard. See A#3 Q5
################################################################################
t = seq(from=0,to=10,length=101)
Density = exp(-8/3)* (t-2)^2 * exp(-(t-2)^3/3)
plot(t,Density,type='l',main='Strange Density with h(t) = (t-2)^2')
bowldat = read.table("https://raw.githubusercontent.com/dsoave/Datasets/refs/heads/main/bowlhaz.data.txt")
head(bowldat); summary(bowldat)

Time = bowldat$Time; Uncensored = bowldat$Uncensored # Writing over earlier vars
# Weibull minus log likelihood again
mloglike = function(theta,t,delta)
     { # Minus log likelihood function for Weibull
     alpha = theta[1]; lambda = theta[2]
     # logf and logS will be of length n
     logf = log(alpha)+log(lambda)+(alpha-1)*log(lambda*t) + -(lambda*t)^alpha
     logS = -(lambda*t)^alpha
     value = -sum(logf*delta) - sum(logS*(1-delta))
     return(value)
     } # End of function mloglike

################################################################################
# Find MLE
################################################################################

startvals = c(1,1/2)
search = optim(par=startvals, fn=mloglike, t=Time,delta=Uncensored,
                     hessian=TRUE, lower=c(0,0), method='L-BFGS-B')
search

alphahat = search$par[1]; lambdahat = search$par[2]

################################################################################
# Plot hazard function h(t) = alpha*lambda^alpha * t^(alpha-1)
################################################################################

x = seq(from=0,to=5,length=101)
esthazard = alphahat*lambdahat^alphahat * x^(alphahat-1)
truehazard = (x-2)^2

plot(x,truehazard,type='l',xlab='Time',ylab='Hazard', main='Hazard Function for the Bowl Data')
lines(x,esthazard,lty=2)
# Annotate the plot (Make the legend)
x1 = c(0.5,1.5); y1 = c(7,7)
lines(x1,y1,lty=1)
text(2,7,'True h(t)')
x2 = c(0.5,1.5); y2 = c(6,6)
lines(x2,y2,lty=2)
text(2.2,6,'Estimated h(t)')






# This is bad, but the worst part is outside the range of the data: Max = 3.72
# 75th percentile is 0.38
# From HW4, S(t) = exp(-1/3 ((t-2)^3 + 8)
# At t=2, where hazard starts increasing,
exp(-8/3)
# So for 93% of the distribution, the hazard function is decreasing.

# Try estimating the survival function
x = seq(from=0,to=5,length=101)
Shat = exp(-(lambdahat*x)^alphahat)
trueS = exp( -1/3*((x-2)^3 + 8) )
tstring = 'Survival Function for the Bowl Data'
plot(x,trueS,type='l',xlab='Time',ylab='Survival',ylim=c(0,1), main=tstring)
lines(x,Shat,lty=2)
# Annotate the plot (Make the legend)
x1 = c(2.5,3.5); y1 = c(0.8,0.8)
lines(x1,y1,lty=1)
text(4,0.8,'True S(t)')
x2 = x1; y2 = c(0.7,0.7)
lines(x2,y2,lty=2)
text(4.2,0.7,'Estimated S(t)')
