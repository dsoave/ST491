rm(list=ls())
# install.packages("survival",dependencies=TRUE) # Only need to do this once
library(survival) # Do this every time
mux = 50; sdx = 12; n = 500
beta0 = 2.3; beta1 = 0.05

# Simulate one data set
set.seed(9999); options(scipen=999)
x = rnorm(n,mean=mux,sd=sdx)
hist(x)
epsilon = rexp(n)
hist(epsilon)
tstar = exp(beta0 + beta1*x) * epsilon
hist(tstar)
U = abs(rnorm(n, mean=0, sd=100)) # Censoring times
hist(U)
t = pmin(tstar,U) # Vector of minima, n values

delta = as.numeric(tstar<U); table(delta)

# Try exponential regression and proportional hazards regression
stime = Surv(t,delta)
plot(survfit(stime~1))
plot(survfit(stime~(x>50)))

expo = survreg(stime ~ x , dist='exponential')
summary(expo)
beta0;beta1


phaz = coxph(stime ~ x)
summary(phaz)

expo = survreg(stime ~ x , dist='weibull')
summary(expo)




# Pharmaco-smoking
rm(list=ls());  options(scipen=999)
# install.packages("survival",dependencies=TRUE) # Only need to do this once
library(survival) # Do this every time
# install.packages("asaur",dependencies=TRUE) # Only need to do this once
library(asaur)
# Make fixed-up data frame called quit
quit = within(pharmacoSmoking,{
 DayOfRelapse = Surv(ttr+1,relapse)
 contrasts(grp) = contr.treatment(2,base=2) # Patch only is reference category
 colnames(contrasts(grp)) = c('Combo') # Names of dummy vars -- just one
 # Collapse race categories
 Race = as.character(race) # Small r race is a factor. This is easier to modify.
 Race[Race!='white'] = 'blackOther'; Race=factor(Race)
 }) # Finished making data frame quit

head(quit)
w_All = survreg(DayOfRelapse ~ grp + age +  employment + yearsSmoking + levelSmoking + priorAttempts, dist='weibull', data=quit); summary(w_All)

ph_All = coxph(DayOfRelapse ~ grp + age +  employment + yearsSmoking + levelSmoking + priorAttempts, data=quit); summary(ph_All)

# Checking the re-parameterization, a Weibull -beta_j/sigma should be close to a proportional hazards beta_j. Look at Racewhite.

-1.092580/1.72  # Compare  -0.1394286 from PH
# Not bad. This is n=125

wfull =  survreg(DayOfRelapse ~ grp + age + employment , dist='weibull', data=quit)
summary(wfull)

phfull =  coxph(DayOfRelapse ~ grp + age + employment, data=quit)
summary(phfull)

# How are they getting the confidence intervals for those hazard ratios?


# Try Partial Likelihood and Wald tests for employment, controlling for age and
# experimental treatment.

# Partial Likelihood Ratio test
nojob = coxph(DayOfRelapse ~ grp + age, data=quit)
anova(nojob,phfull) # LR test
