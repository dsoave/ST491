rm(list=ls()); options(scipen=999)
wdata = read.table("https://raw.githubusercontent.com/dsoave/Datasets/refs/heads/main/Weibull.data2.txt")
# head(wdata)
Time = wdata$Time; Uncensored = wdata$Uncensored # Avoiding the attach() function

# install.packages("survival",dependencies=TRUE) # Only need to do this once
library(survival) # Do this every time


################################################################################
# (a)
################################################################################

y = Surv(Time,Uncensored); y[1:20] # A pre-processing step
km1 = survfit(y ~ 1) # Like a regression model with just an intercept: No x values
km1
# For comparison, MLE of the median was 4.466, 95% CI = (4.199 4.733)

summary(km1) # Returns a matrix, first column t_j etc.

phat1=261/262; phat2=251/252; phat3=239/240
Shat = phat1*phat2*phat3; Shat # Compare S-hat(1.07) = 0.9881

# Get SE of S-hat(1.07):
estvarlog = 1/(262*261) + 1/(252*251) + 1/(239*240)
estvarShat = estvarlog * Shat^2 # One-var delta method
seShat = sqrt(estvarShat); seShat # Compare 0.00684
Shat + 1.96*seShat
# Upper confidence limit was truncated to one.

################################################################################
# (b)
################################################################################

plot(km1)
title('Kaplan-Meier Estimate for the Weibull Data')


# These data were simulated, so I know the true parameter values.
# Add true S(t) to the plot
truealpha = 3; truelambda = 1/5
x = seq(from=0,to=8,length=101)
trueS = exp(-(truelambda*x)^truealpha)
lines(x,trueS,lty=1, col = "red1")

# Parametric (Weibull) estimate of S(t) and median, for comparison.

mloglike = function(theta,t,delta)
     { # Minus log likelihood function
     alpha = theta[1]; lambda = theta[2]
     # logf and logS will be of length n
     logf = log(alpha)+log(lambda)+(alpha-1)*log(lambda*t) + -(lambda*t)^alpha
     logS = -(lambda*t)^alpha
     value = -sum(logf*delta) - sum(logS*(1-delta))
     return(value)
     } # End of function mloglike

############
# (c) Find MLE #
############

startvals = c(1,1/2) # I tried a few values
search1 = optim(par=startvals, fn=mloglike, t=Time,delta=Uncensored,
                     hessian=TRUE, lower=c(0,0), method='L-BFGS-B')
# search1
alphahat = search1$par[1]; lambdahat = search1$par[2]

# Compare true and estimated median
truealpha = 3; truelambda = 1/5
H = search1$hessian
Vhat = solve(H) # Solve returns the inverse.

################################################################################
# (d) Point estimate for the median
# Median = log(2)^(1/alpha) / lambda
################################################################################

# Point estimate of median
medhat = 1/lambdahat * log(2)^(1/alphahat); medhat
# Compare the truth
truemedian = log(2)^(1/truealpha) / truelambda; truemedian

################################################################################
# (e)
################################################################################

# Estimate the survival function
x = seq(from=0,to=10,length=101)
Shat = exp(-(lambdahat*x)^alphahat)
trueS = exp(-(truelambda*x)^truealpha)
tstring = 'Survival Function for the Weibull Data'
plot(x,trueS,type='l',xlab='Time',ylab='Survival',ylim=c(0,1), main=tstring, col = "red1")
lines(x,Shat,lty=2, col = "blue1")
# Annotate the plot (Make the legend)
x1 = c(6,7.5); y1 = c(0.8,0.8)
lines(x1,y1,lty=1)
text(8.5,0.8,'True S(t)', col = "red1")
x2 = x1; y2 = c(0.7,0.7)
lines(x2,y2,lty=2)
text(8.9,0.7,'Estimated S(t)', col = "blue1")

################################################################################
# (f) Add MLE to Kaplan-Meier plot
################################################################################

plot(km1)
title("Kaplan-Meier Estimate for the Weibull Data")
lines(x,trueS,lty=1, col = "red1")
lines(x,Shat,lty=2, col = "blue1")

################################################################################
# (g) What if the (Weibull) model is wrong? Bowl (shaped hazard) data
################################################################################

rm(list=ls()); options(scipen=999)
bowldat = read.table("http://www.utstat.toronto.edu/brunner/data/legal/bowlhaz.data.txt")
head(bowldat); summary(bowldat); attach(bowldat)


# install.packages("survival",dependencies=TRUE) # Only need to do this once
library(survival) # Do this every time

################################################################################
# (g) (i)
################################################################################

y = Surv(Time,Uncensored) # A pre-processing step
km2 = survfit(y ~ 1) # A regression model with just an intercept: No x values
km2

# What is the true median?
# From HW3, S(t) = exp(-1/3 ((t-2)^3 + 8)
# Set exp(-1/3 ((t-2)^3 + 8) = 1/2, get t = 0.190935

# Plot true S(t) and MLE for comparison to K-M


################################################################################
# (g) (ii)Find MLE
################################################################################
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

startvals = c(1,1/2)
search = optim(par=startvals, fn=mloglike, t=Time,delta=Uncensored,
                     hessian=TRUE, lower=c(0,0), method='L-BFGS-B')
alphahat = search$par[1]; lambdahat = search$par[2]

################################################################################
# (g) (iii)
################################################################################
# MLE of median
medhat = 1/lambdahat * log(2)^(1/alphahat); medhat
# Compared to truth of 0.191 and K-M estimate of 0.178

0.203-0.191 # Error of MLE

0.178-0.191 # Error of Kaplan-Meier

################################################################################
# (g) (iv)
################################################################################

#True S(t)
x = seq(from=0,to=8,length=101)
trueS = exp(-1/3*((x-2)^3 + 8))

plot(km2)
# Add title and other estimates to Kaplan-Meier plot
title('Kaplan-Meier Estimate for the Bowl Data')
lines(x,trueS,lty=1, col = "red1")

################################################################################
# (g) (v)
################################################################################
Shat = exp(-(lambdahat*x)^alphahat)
lines(x,Shat,lty=2, col = "blue1")

# Annotate the plot (Make the legend)
x1 = c(2,2.75); y1 = c(0.8,0.8)
lines(x1,y1,lty=1, col = "red1")
text(3,0.8,'True S(t)', col = "red1")
x2 = x1; y2 = c(0.7,0.7)
lines(x2,y2,lty=1)
text(3.15,0.7,'Kaplan-Meier')
x3 = x1; y3 = c(0.6,0.6)
lines(x3,y3,lty=2,  col = "blue1")
text(3.15,0.6,'Weibull MLE', col = "blue1")
