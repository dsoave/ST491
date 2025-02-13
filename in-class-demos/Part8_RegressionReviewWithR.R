kars = read.table("https://raw.githubusercontent.com/dsoave/Datasets/refs/heads/main/mcars4.data")
dim(kars)
head(kars)
summary(kars)

n = dim(kars)[1]; n
# Modifying the data frame, not the original data file
kars = within(kars,{
+ # Make dummy variables
+ c1 = numeric(n); c1[Cntry=='Europ'] = 1
+ c2 = numeric(n); c2[Cntry=='Japan'] = 1
+ c3 = numeric(n); c3[Cntry=='US'] = 1
+ # Make Cntry a factor
+ Cntry = factor(Cntry)
+ }) # End of within kars

head(kars)
summary(kars)

# Checking dummy variables
with(kars,table(c1,Cntry))
with(kars,table(c2,Cntry))
with(kars,table(c3,Cntry))

# Take a look at mean fuel consumption for each country
with(kars, aggregate(lper100k,by=list(Cntry),FUN=mean) )

# Must specify a LIST of grouping factors


# H0: mu1=mu2=mu3
justcountry = lm(lper100k ~ c1+c2, data=kars)
summary(justcountry)

# Which means are different?
# Have t-tests. What about Europe vs. Japan?
# test H0: beta1 = beta2
# A cheap way is to use a different reference category.
# R can make the dummy variables for you
is.factor(kars$Cntry)
# The factor Cntry has dummy vars built in. What are they?
contrasts(kars$Cntry) # Note alphabetical order

jc2 = lm(lper100k~Cntry, data=kars); summary(jc2)


# You can select the dummy variable coding scheme and the reference category.
contr.treatment(3,base=2) # Category 2 is the reference category


# U.S. as reference category again
kars$Country = kars$Cntry
contrasts(kars$Country) = contr.treatment(3,base=3)
summary( lm(lper100k~Country, data=kars) )


# Names of dummy variables 1=Europe, 2=Japan could be nicer
colnames(contrasts(kars$Country)) = c("Europe","Japan")
contrasts(kars$Country)

summary( lm(lper100k~Country, data=kars) )
# Include covariates
fullmodel = lm(lper100k ~ weight+length+Country, data=kars)
summary(fullmodel) # Look carefully at the signs!


# Test country controlling for size, using full versus reduced approach.

justsize = lm(lper100k ~ weight+length, data=kars); summary(justsize)

with(kars, cor(weight,length)  )


# I advise using anova ONLY to compare full and reduced models
anova(justsize,fullmodel) # Full vs reduced

# Test car size controlling for country too -- why not?
anova(justcountry,fullmodel)

###### Predictions, confidence intervals and prediction intervals ######

# Predict litres per 100 km for a Japanese car weighing
# 1295kg, 4.52m long (1990 Toyota Camry)

b = fullmodel$coefficients; b

ell = c(1,1295,4.52,1,0)
yhat = sum(ell*b); # ell-prime b
yhat

# Confidence interval for E(ell-prime beta)
# First the hard way

tcrit = qt(0.975,df=fullmodel$df.residual) # t_alpha/2
MSE.XpXinv = vcov(fullmodel)
ell = as.matrix(ell) # Now it's a column vector
me95 = tcrit * sqrt( as.numeric(t(ell) %*% MSE.XpXinv %*% ell) )
lower95 = yhat - me95; upper95 = yhat + me95
c(lower95, upper95) # 95% Confidence interval for ell-prime beta

# Use the predict function
# help(predict.lm)

camry1990 = data.frame(weight=1295,length=4.52,Cntry='Japan')
camry1990
predict(fullmodel,newdata=camry1990) # Compare yhat = 12.38739
predict(fullmodel,newdata=camry1990, interval='confidence')

# With 95 percent prediction interval (95 is default)

predict(fullmodel,newdata=camry1990, interval='prediction')

# Multiple predictions
cadillac1990 = data.frame(weight=1800,length=5.22,Cntry='US')
volvo1990 = data.frame(weight=1371,length=4.823,Cntry='Europ')
newcars = rbind(camry1990,cadillac1990,volvo1990); newcars

is.data.frame(newcars)
predict(fullmodel,newdata=newcars, interval='prediction')







