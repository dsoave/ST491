# Title: Review: Normal Regression with R*

# Step 1: Load the Data
# Load the dataset directly from GitHub
kars <- read.table("https://raw.githubusercontent.com/dsoave/Datasets/refs/heads/main/mcars4.data")
dim(kars)  # Check dimensions
head(kars) # View first few rows
summary(kars) # Summary statistics

# Step 2: Create Dummy Variables
# Create dummy variables for the 'Cntry' column
n <- dim(kars)[1]
kars <- within(kars, {
  c1 <- numeric(n); c1[Cntry == 'Europ'] <- 1
  c2 <- numeric(n); c2[Cntry == 'Japan'] <- 1
  c3 <- numeric(n); c3[Cntry == 'US'] <- 1
  Cntry <- factor(Cntry)
})

# Validate the results
head(kars)
summary(kars)

# Step 3: Validate Dummy Variables
# Check the dummy variables
with(kars, table(c1, Cntry))
with(kars, table(c2, Cntry))
with(kars, table(c3, Cntry))

# Step 4: Analyze Mean Fuel Consumption
# Mean fuel consumption by country
with(kars, aggregate(lper100k, by = list(Cntry), FUN = mean))

# Step 5: Fit a Linear Regression Model
# Test if fuel consumption differs by country
justcountry <- lm(lper100k ~ c1 + c2, data = kars)
summary(justcountry)

# Step 6: Change Reference Categories
# Change reference category to "US"
kars$Country <- kars$Cntry
contrasts(kars$Country) <- contr.treatment(3, base = 3)
summary(lm(lper100k ~ Country, data = kars))

# Step 7: Add Covariates to the Model
# Full model with covariates
fullmodel <- lm(lper100k ~ weight + length + Country, data = kars)
summary(fullmodel)

# Step 8: Compare Models
# Reduced model for comparison
justsize <- lm(lper100k ~ weight + length, data = kars)
anova(justsize, fullmodel)  # Full vs reduced

# Compare country effects controlling for size
anova(justcountry, fullmodel)

……………. Remove!!!!!
# General linear test approach. H0: L beta = h

source("http://www.utstat.utoronto.ca/brunner/Rfunctions/ftest.txt")

# Test country controlling for size: Compare F = 6.8999
# Full model again for comparison
summary(fullmodel)

# Now the F-test of couuntry controlling for size
L0 = rbind(c(0,0,0,1,0),
            c(0,0,0,0,1))
L0

ftest(fullmodel,L0)

# As before, t-tests give comparison of U.S with Europe and Japan.
# Test Europe vs. Japan controlling for size.
L1 = cbind(0,0,0,1,-1) # One row, 5 columns
ftest(fullmodel,L1)

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




