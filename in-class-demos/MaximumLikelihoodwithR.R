################################################################################
# Question 2: Normal random sample
################################################################################

rm(list=ls()); options(scipen=999) # Option to suppress scientific notation
# Read normal data
normaldata = scan("http://www.utstat.toronto.edu/brunner/data/legal/normal.data.txt")

################################################################################
# (a) Find the maximum likelihood estimates of mu and sigma-squared numerically.
################################################################################
mloglike = function(theta,x) # theta is (mu,sigmasq)
     {
     mu = theta[1]; sigmasq = theta[2]
     n = length(x)
     value = n/2*log(sigmasq) + n/2*log(2*pi) + sum((x-mu)^2)/(2*sigmasq)
     return(value)
     } # End function mloglike

# Test the function -- easy in this case
mloglike(c(90,220),normaldata)

-sum(dnorm(normaldata,90,sqrt(220),log=TRUE)) # Get log of density as an option

# help(optim) # Minimize a function
startvals = c(0,1) # Of course the sample mean and variance would be better.

normalsearch = optim(par=startvals, fn=mloglike, x=normaldata,
                      hessian=TRUE, lower=c(-Inf,0), method='L-BFGS-B')


normalsearch


# If you really have no idea what the parameter values are, it's wise to try
# several sets of starting values. Here's a really lucky guess.
optim(par=c(mean(normaldata),var(normaldata)), fn=mloglike, x=normaldata,
             hessian=TRUE, lower=c(-Inf,0), method='L-BFGS-B')


# We will stay with the first solution.

################################################################################
# Compare the answer to your closed-form solution.
################################################################################

muhat = mean(normaldata)
n = length(normaldata); sigmasqhat = (n-1)*var(normaldata)/n
c(muhat,sigmasqhat) # The exact MLE

################################################################################
# (b) Show that the minus log likelihood is indeed minimized at the MLE for
#     this data set.
################################################################################

# If both eigenvalues of the Hessian are positive, it's concave up.
H = normalsearch$hessian
eigen(H)$values

################################################################################
# (c) Calculate the estimated asymptotic covariance matrix of the MLEs.
################################################################################

Vhat = solve(H); Vhat # Solve returns the inverse.

################################################################################
# (d) Give a ``better" estimated asymptotic covariance matrix based on your
#     closed-form solution.
################################################################################

betterVhat = rbind(c(sigmasqhat/n, 0               ),
                    c(0           , 2*sigmasqhat^2/n) )
betterVhat     

################################################################################
# (e) Calculate a large-sample 95\% confidence interval for sigma-squared.
################################################################################

se = sqrt(Vhat[2,2])
sigmasqhat = normalsearch$par[2] # Replacing the exact answer. It's almost
                                   # the same.
lower95 = sigmasqhat - 1.96*se; upper95 = sigmasqhat + 1.96*se
c(lower95,upper95)



################################################################################
# (f) Test H0: mu = 103 with
    # (i) Z test
    # (ii) Likelihood ratio chi-squared test. Compare the closed-form version.
    # (iii) Wald chi-squared test.
# Give the test statistic and the p-value for each test.
################################################################################

################################################################################
# (i) Z test
################################################################################
se2 = sqrt(Vhat[1,1]); Z = (muhat-103)/se2
pval = 2 * (1-pnorm(abs(Z)))
c(Z,pval)

################################################################################
# (ii) Likelihood ratio chi-squared test
################################################################################

# LR test could be done in closed form, but do it numerically
# Function for constrained minimization
mloglike0 = function(theta,x) # theta is just sigmasq
     {
     mu = 103; sigmasq = theta
     n = length(x)
     value = n/2*log(sigmasq) + n/2*log(2*pi) + sum((x-mu)^2)/(2*sigmasq)
     return(value)
     } # End function mloglike0

# Estimate. Don't ask for the Hessian.
restricted1 = optim(par=sigmasqhat, fn=mloglike0, x=normaldata,
               lower=0, method='L-BFGS-B')
restricted1


Gsq = 2 * (restricted1$value - normalsearch$value)
pval = 1-pchisq(Gsq,1)
c(Gsq,pval)
Gsq = 2 * (restricted1$value - normalsearch$value)
pval = 1-pchisq(Gsq,1)
c(Gsq,pval)


################################################################################
# Compare the closed-form version.
################################################################################

sigmasqhat0 = sum( (normaldata-103)^2 )/n
n * ( log(sigmasqhat0) - log(sigmasqhat) ) # G-squared in closed form

