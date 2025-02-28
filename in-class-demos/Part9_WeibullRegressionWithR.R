# Part9_WeibullRegressionWithR.R

rm(list=ls()); # options(scipen=999)
# install.packages("survival",dependencies=TRUE) # Only need to do this once
library(survival) # Do this every time
# install.packages("asaur",dependencies=TRUE) # Only need to do this once
library(asaur)
# help(pharmacoSmoking)

head(pharmacoSmoking)
summary(pharmacoSmoking)

quit = within(pharmacoSmoking, {TimeToRelapse = Surv(ttr,relapse)} )
sort(quit$TimeToRelapse)

# Fit a Weibull model with just treatment group, which was randomly assigned.
Model0 = survreg(TimeToRelapse~grp,dist='weibull', data=quit)
#Error in survreg(TimeToRelapse ~ grp, dist = "weibull", data = quit) :
#  Invalid survival times for this distribution
# Likely choking on the zeros.
# Day of relapse starts with one.
quit = within(quit, {DayOfRelapse = Surv(ttr+1,relapse)} )
Model0 = survreg(DayOfRelapse~grp,dist='weibull', data=quit)
summary(Model0)

#Conclusion is that combination therapy is more effective. 
# But the alphabetical order of treatments makes combination the reference category, and this is clumsy. 
# Make patch-only the reference category and re-run. Write down the table of population means for the two groups. Now write the one we want.

quit = within(quit,{
 contrasts(grp) = contr.treatment(2,base=2)
 colnames(contrasts(grp)) = c('Combo') # Names of dummy vars -- just one
 })
Model1 = survreg(DayOfRelapse~grp,dist='weibull', data=quit)
summary(Model1)

betahat = Model1$coefficients; betahat
betahat0 = betahat[1]; betahat1 = betahat[2]
sigmahat = Model1$scale; sigmahat
Vhat = vcov(Model1); Vhat

# Asymptotic covariance matrix comes out in terms of Log(scale), which is
# a minor pain.

# 1) When patients receive the combination drug therapy rather than nicotine patch only, expected relapse time is multiplied by _______ .
    # a) Give an estimate
    # b) Modify the CI for beta1 to get a 95% confidence interval (don't use the delta method).

# a) Give an estimate
exp(betahat1)
grpCombo
3.196004

# b) Modify the CI for beta1 to get a 95% confidence interval (don't use the delta method).
L = 1.162 - 1.96*0.3999; U = 1.162 + 1.96*0.3999
c(exp(L),exp(U))

summary(Model1) # Repeating

# 2) Estimate and plot the density of relapse time for the two experimental conditions.

# Okay, lambda = exp(-mu), alpha = 1/sigma
alpha = 1/sigmahat
lambda0 = exp(-betahat0); lambda1 = exp(-betahat0-betahat1)
t = seq(from=0,to=300,length=101)
f0 = alpha*lambda0^alpha * t^(alpha-1) * exp(-(lambda0*t)^alpha)
f1 = alpha*lambda1^alpha * t^(alpha-1) * exp(-(lambda1*t)^alpha)
plot(t,f0,pch=' ',xlab='Day of Relapse',ylab='Density') # Empty plot
title('Estimated Density of Relapse Time')
lines(t,f0,lty=2); lines(t,f1,lty=1)
# Annotate the plot
x0=c(200,250); y0 = c(0.015,0.015); lines(x0,y0,lty=2)
text(300,0.015,'Patch Only')
x1=c(200,250); y1 = c(0.017,0.0171); lines(x1,y1,lty=1)
text(300,0.017,'Combination')

max(quit$ttr) # Maximum time value (censored)

# 3) Estimate median time to relapse for the 2 groups, with Cis
