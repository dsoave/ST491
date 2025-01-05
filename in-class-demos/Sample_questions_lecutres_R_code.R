##### Part 1 - Sample Questions: Maximum Likelihood Part One
rm(list=ls())
x = scan("https://raw.githubusercontent.com/dsoave/Datasets/refs/heads/main/pareto.data.txt")
x

thetahat = 1/mean(log(x)); thetahat
n = length(x); vhat = thetahat^2/n; se = sqrt(vhat); se
low95 = thetahat - 1.96*se; up95 = thetahat + 1.96*se
c(low95,up95) # 95% CI

# Critical value(s). Just say plus and minus 1.96, or ...
c(qnorm(0.025),qnorm(0.975))
Z = (thetahat-1.16)/se; Z
pvalue = 2 * (1-pnorm(abs(Z))); pvalue

