############
# Sample Questions for Longnormal Regression
# #15 Log-normal Regression with R*


rm(list=ls());  options(scipen=999)
# install.packages("survival",dependencies=TRUE) # Only need to do this once
library(survival) # Do this every time
# install.packages("asaur",dependencies=TRUE) # Only need to do this once
library(asaur)
# summary(pharmacoSmoking)
# Make fixed-up data frame called quit
quit = within(pharmacoSmoking,{
 DayOfRelapse = Surv(ttr+1,relapse)
 contrasts(grp) = contr.treatment(2,base=2) # Patch only is reference category
 colnames(contrasts(grp)) = c('Combo') # Names of dummy vars -- just one
 # Collapse race categories
 Race = as.character(race) # Small r race is a factor. This is easier to modify.
 Race[Race!='white'] = 'blackOther'; Race=factor(Race)
 }) # Finished making data frame quit
# with(quit, table(race,Race) )

#################################
# (a)
wmod =  survreg(DayOfRelapse ~ grp + age + employment , dist='weibull', data=quit)
summary(wmod) # This was model full2 in an earlier analysis

#################################
# (b)
lognorm = survreg(DayOfRelapse ~ grp + age + employment,dist='lognormal', data=quit)
summary(lognorm)

#################################
# (c)
# Now predict the day of relapse for a 50-year-old in the patch-only condition
# who is working part-time. Use the estimated median exp(muhat) as a prediction.

###############################
# First do it the hard way

betahat = lognorm$coefficients; betahat

xnplus1 = cbind(1,0,50,0,1); xnplus1
sigmahat = lognorm$scale
Vn = vcov(lognorm); round(Vn,4)

Cn = Vn[(1:5),(1:5)]; round(Cn,4)

yhat = sum(xnplus1*betahat); yhat # This is estimated x'beta
exp(yhat) # This is predicted number of days
se = sqrt( sigmahat^2 + as.numeric(xnplus1 %*% Cn %*% t(xnplus1)) ); se
A = yhat-1.96*se; B = yhat+1.96*se; c(A,B)
c( exp(A), exp(B)) # Prediction interval
exp(A)*24  # Lower limit in hours
exp(B)/365 # Upper limit in years

# So with 95% confidence, the old guy will be able to hold out between 
# 5.2 hours and 3.9 years. Thank you very much.

# Now the easy way
oldguy = data.frame(grp='patchOnly',age=50,employment='pt')
pred1 = predict(lognorm,newdata=oldguy,type='linear',se=TRUE) ; pred1

# Construct prediction interval
se = sqrt(sigmahat^2+pred1$se^2)
L = yhat - 1.96*se; U = yhat + 1.96*se
t_hat= exp(yhat)
lower95 = exp(L); upper95 = exp(U)
pi = c(t_hat,lower95,upper95)
names(pi) = c('t-hat','lower95','upper95')

pi



#Comparing Survival Curves (Weibull VS Lognormal)


t = 1:365

xnplus1 = cbind(1,0,50,0,0); xnplus1
sigmahat = lognorm$scale
sigmahatW = wmod$scale
betahat = lognorm$coefficients; betahat
betahatW = wmod$coefficients; betahatW

yhat = sum(xnplus1*betahatW); yhat # This is estimated x'beta
exp(yhat)

Shat = exp( -(exp(-yhat/sigmahatW)*t^(1/sigmahatW)) )
plot(t,Shat,type='l',ylim=c(0,1),xlab='t',ylab='Probability that Relapse is After Day t')
tstring = expression(paste(hat(S)(t), " = Probability Relapsing After Day t"))
title(tstring)
lines(t,Shat,col=1)
points(exp(yhat) * log(2)^sigmahatW,0.5)

yhat = sum(xnplus1*betahat); yhat # This is estimated x'beta
ShatLN = 1-pnorm((log(t)-yhat)/sigmahat)
lines(t,ShatLN,lty=2)
points(exp(yhat),0.5)

xnplus1 = cbind(1,1,50,0,0); xnplus1
yhat = sum(xnplus1*betahatW); yhat # This is estimated x'beta
exp(yhat)

Shat = exp( -(exp(-yhat/sigmahatW)*t^(1/sigmahatW)) )
title(tstring)
lines(t,Shat,col=2)
points(exp(yhat) * log(2)^sigmahatW,0.5)

yhat = sum(xnplus1*betahat); yhat # This is estimated x'beta
ShatLN = 1-pnorm((log(t)-yhat)/sigmahat)
lines(t,ShatLN,lty=2,col=2)
points(exp(yhat),0.5)
