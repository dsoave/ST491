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

# Asymptotic covariance matrix comes out in terms of Log(scale), which is
# unfortunate.
# Denote log(sigma) by s, and re-write g(theta) = exp(beta0) * log(2)^sigma
# as g(theta) = exp(beta0) * log(2)^exp(s)
shat = log(sigmahat)

# Patch Only
medianhat0 = exp(betahat0)*log(2)^sigmahat
# gdot will be a 1 x 3 matrix.
gdot0 = cbind( exp(betahat0)*log(2)^exp(shat), 0, exp(betahat0) * log(2)^exp(shat) * log(log(2)) * exp(shat) )
se0 = sqrt( as.numeric(gdot0 %*% Vhat %*% t(gdot0)) ); se0
lower0 = medianhat0 - 1.96*se0; upper0 = medianhat0 + 1.96*se0
patchonly = c(medianhat0,lower0,upper0)
names(patchonly) = c('Median','Lower95','Upper95')
patchonly

# Combination drug treatment
medianhat1 = exp(betahat0+betahat1)*log(2)^sigmahat
gdot1 = cbind( exp(betahat0+betahat1)*log(2)^exp(shat), exp(betahat0+betahat1)*log(2)^exp(shat),
                exp(betahat0+betahat1) * log(2)^exp(shat) * log(log(2)) * exp(shat) )
se1 = sqrt( as.numeric(gdot1 %*% Vhat %*% t(gdot1)) ); se1
lower1 = medianhat1 - 1.96*se1; upper1 = medianhat1 + 1.96*se1
combination = c(medianhat1,lower1,upper1)
names(combination) = c('Median','Lower95','Upper95')
combination

# There is an easier way to get these numbers
# help(predict.survreg)
Justpatch = data.frame(grp='patchOnly')
# A data frame with just one case and one variable.
Combination = data.frame(grp='combination')
treatments = rbind(Justpatch,Combination); treatments

# The 0.5 quantile is the median
medians = predict(Model1,newdata=treatments,type='quantile',p=0.5,se=TRUE)
medians

cbind(medians$fit,medians$se)
rbind(c(medianhat0,se0),
       c(medianhat1,se1) )

# 4) Plot the Kaplan-Meier estimates and MLEs of S(t)

KM = survfit(DayOfRelapse~grp, type="kaplan-meier", data=quit)  
# Kaplan-Meier is the default estimation method anyway.
summary(KM)

# Look at K-M estimates of medians
KM[1] # Combination
KM[2] # Patch Only

# Repeat MLEs for comparison
combination
patchonly

plot(KM,  xlab='t', ylab='Survival Probability', lwd=2, col=1:2)
# 1 is black,  2 is red
legend(x=125,y=1.0, col=1:2, lwd=2, legend=c('Combination','Patch Only'))
title(expression(paste(hat(S)(t),': Kaplan-Meier and Maximum Likelihood Estimates')))
# MLEs

x = 1:185
lambda0 = exp(-betahat0); lambda1 = exp(-betahat0-betahat1); alpha=1/sigmahat
Shat0 = exp(-(lambda0*x)^alpha); Shat1 = exp(-(lambda1*x)^alpha)
lines(x,Shat0,lty=2,col=2) # Patch only is red
lines(x,Shat1,lty=2)       # Combination is black (default)

# Non-parametric rank test of equal survival functions
# See http://dwoll.de/rexrepos/posts/survivalKM.html.
survdiff(DayOfRelapse~grpgrp, data=quit)

# Compare p = 0.00367 from Z-test of H0: beta1=0
twoby2
prop.table(twoby2,1) # Proportions of row totals
chisq.test(twoby2)

fisher.test(twoby2) # p = 0.01713






