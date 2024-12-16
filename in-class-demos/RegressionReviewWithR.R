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
