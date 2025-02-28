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

# Conclusion is that combination therapy is more effective. 
# But the alphabetical order of treatments makes combination 
# the reference category, and this is clumsy. 
# Make patch-only the reference category and re-run. 
# Write down the table of population means for the two groups. 
# Now write the one we want.

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

# 1) When patients receive the combination drug therapy rather than 
#    nicotine patch only, expected relapse time is multiplied by _______ .
    # a) Give an estimate
    # b) Modify the CI for beta1 to get a 95% confidence interval 
    #    (don't use the delta method).

# a) Give an estimate
exp(betahat1)

# b) Modify the CI for beta1 to get a 95% confidence interval 
#    (don't use the delta method).
L = 1.162 - 1.96*0.3999; U = 1.162 + 1.96*0.3999
c(exp(L),exp(U))

summary(Model1) # Repeating

# 2) Estimate and plot the density of relapse time for the two 
#    experimental conditions.

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
medianhat0 = exp(betahat0)*log(2)^sigmahat; medianhat0
# CI requires SE estimate using multivariate delta method (not this course!)

# Combination drug treatment
medianhat1 = exp(betahat0+betahat1)*log(2)^sigmahat;medianhat1

# There is an easier way to get these numbers using R (no need for the delta method by hand)
# help(predict.survreg)
Justpatch = data.frame(grp='patchOnly')
# A data frame with just one case and one variable.
Combination = data.frame(grp='combination')
treatments = rbind(Justpatch,Combination); treatments

# The 0.5 quantile is the median
medians = predict(Model1,newdata=treatments,type='quantile',p=0.5,se=TRUE)
medians

cbind(medians$fit,medians$se,medians$fit-1.96*medians$se,medians$fit+1.96*medians$se)

# 4) Plot the Kaplan-Meier estimates and MLEs of S(t)

KM = survfit(DayOfRelapse~grp, type="kaplan-meier", data=quit)  
# Kaplan-Meier is the default estimation method anyway.
summary(KM)

# Look at K-M estimates of medians
km_med1<-summary(KM[1])$table["median"];km_med1 # Combination
km_med2<-summary(KM[2])$table["median"];km_med2 # Patch Only

# Repeat MLEs for comparison
cbind(medians$fit,medians$se,medians$fit-1.96*medians$se,medians$fit+1.96*medians$se)


plot(KM,  xlab='t', ylab='Survival Probability', lwd=2, col=1:2)
# 1 is black,  2 is red
legend(x=125,y=1.0, col=1:2, lwd=2, legend=c('Combination','Patch Only'))
title(expression(paste(hat(S)(t),': Kaplan-Meier and Maximum Likelihood Estimates')))
#medians
abline(h=0.5)
points(c(km_med1,km_med2),c(0.5,0.5),pch=16) #medians


# MLEs

x = 1:185
lambda0 = exp(-betahat0); lambda1 = exp(-betahat0-betahat1); alpha=1/sigmahat
Shat0 = exp(-(lambda0*x)^alpha); Shat1 = exp(-(lambda1*x)^alpha)
lines(x,Shat0,lty=2,col=2) # Patch only is red
lines(x,Shat1,lty=2)       # Combination is black (default)
points(c(medians$fit),c(0.5,0.5),pch=16) #medians

# Non-parametric rank test of equal survival functions
# See http://dwoll.de/rexrepos/posts/survivalKM.html.
survdiff(DayOfRelapse~grp, data=quit)

# Compare p = 0.00367 from Z-test of H0: beta1=0

twoby2<-table(quit$grp,quit$relapse);twoby2
prop.table(twoby2,1) # Proportions of row totals
chisq.test(twoby2)

fisher.test(twoby2) # p = 0.01713


#####################
# Next we are going to look at some other covariates in larger
# multiple variabel model

# Make fixed-up data frame called quit
quit = within(pharmacoSmoking,{
 DayOfRelapse = Surv(ttr+1,relapse)
 contrasts(grp) = contr.treatment(2,base=2) # Patch only is reference category
 colnames(contrasts(grp)) = c('Combo') # Names of dummy vars -- just one
 # Collapse race categories
 Race = as.character(race) # Small r race is a factor. This is easier to modify.
 Race[Race!='white'] = 'blackOther'; Race=factor(Race)
 }) # Finished making data frame quit
 with(quit, table(race,Race) )

full = survreg(DayOfRelapse ~ grp + age + gender + Race + employment
        + yearsSmoking + levelSmoking + priorAttempts, dist='weibull', data=quit)
summary(full)
 
# I am thinking about dropping Race, yearsSmoking, levelSmoking and priorAttempts. 
# The last 3 variables all represent smoking history and could be correlated 
# highly enough to wash out each other's effects. Test them simultaneously.

# Fit the restricted model: Restricted by H0
rest1 =  survreg(DayOfRelapse ~ grp + age + gender + Race + employment ,
                  dist='weibull', data=quit)
anova(rest1,full) # LR test


# Is Race significant with those variables dropped?
summary(rest1)

# Decision: Drop race and gender.
full2 =  survreg(DayOfRelapse ~ grp + age + employment , dist='weibull', data=quit)
summary(full2)

# Test employment status controlling for age and experimental treatment.
rest2 = survreg(DayOfRelapse ~ grp + age , dist='weibull', data=quit)
anova(rest2,full2) # LR test

# Predict the day of relapse for a 50 year old patient who is employed 
# full time and gets the patch-only treatment.
thetahat<-full2$coefficients;thetahat
sigmahat<-full2$scale
  
x = c(1,0,50,0,0)
xb = sum(x*thetahat)


# a) The estimated mean
exp(xb) * gamma(sigmahat+1)

# b) The estimated mean
exp(xb) * log(2)^sigmahat

# I think the median is preferable to mean because the Weibull distribution
# is skewed. Also, the predict function for Weibull regression works as expected
# for medians (but not means).

oldguy = data.frame(grp='patchOnly',age=50,employment='ft')
predict(full2,newdata=oldguy,type='quantile',p=0.5,se=TRUE)

# The 0.5 quantile is the median. se is from the delta method.

# Estimate and plot S(t) for the old guy.

# The se of S-hat(t) is straightforward in theory, but messy in practice. 
# e.g. requires multivariate delta method

t = 1:365
Shat = exp( -(exp(-xb/sigmahat)*t^(1/sigmahat)) )

plot(t,Shat,type='l',ylim=c(0,1),xlab='t',ylab='Probability that Relapse is After Day t')
tstring = expression(paste(hat(S)(t), " = Probability Relapsing After Day t"))
title(tstring)

# Plot estimated hazard function for that 50 year old patient who is employed 
# full time and gets the patch-only treatment.

h = 1/sigmahat * exp(-xb/sigmahat) * t^(1/sigmahat - 1)
plot(t,h,type='l',xlab='Day',ylab='Risk',main='Estimated Risk of Relapse')

