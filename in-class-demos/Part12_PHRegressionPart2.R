rm(list=ls());  # options(scipen=999)
# install.packages("survival",dependencies=TRUE) # Only need to do this once
library(survival) # Do this every time

help(kidney) ### A description of the data set

# help(kidney)
head(kidney)
summary(kidney)

dim(kidney)

table(kidney$disease)

contrasts(kidney$disease)


# Make a new data frame with
     # 1=F, 0=M
     # age and frailty centered

Kidney = within(kidney,{
 sex = sex-1 # Indicator for female
 # Centering age and frailty
 age = age-mean(age)
 frail = frail-mean(frail)
 })
with(Kidney,cor(age,frail))

kmod1 = coxph( Surv(time,status) ~ age + sex + disease, data=Kidney)
summary(kmod1)
# Are se(coef) labelled correctly?
se = sqrt(diag(vcov(kmod1))); se # Yes

# CI for the hazard ratio exp(beta1)
betahat = coef(kmod1); betahat

CIbeta1 = c(betahat[1]-1.96*se[1], betahat[1]+1.96*se[1]); CIbeta1

exp(CIbeta1)
# So summary is giving us confidence intervals for the hazard ratios,
# not the coefficients.


summary(kmod1)

# (1) Estimated hazard of infection is ____ times as great for female as male
# (2) Estimated hazard of infection is ____ times as great for disease type AN as it is for Other.
# (3) Estimated hazard of infection is ____ times as great for disease type AN as it is for disease type PKD.
betahat = coef(kmod1); betahat

exp(betahat[4]-betahat[5]) # Hazard ratio of AN/PKD


# Test disease type with  a partial likeihood ratio test
k2 = coxph( Surv(time,status) ~ age + sex, data=Kidney)
anova(k2,kmod1)


# Comparing survival functions for males and females
male = data.frame(age=0, sex=0, disease="Other") # An average guy
female = data.frame(age=0, sex=1, disease="Other") # An average gal
sexcomp = rbind(male,female); sexcomp

rownames(sexcomp) = c("M","F"); sexcomp

s1 = survfit(kmod1,newdata=sexcomp)
s1


ls(s1)

summary(s1)

head(s1$cumhaz)


S = s1$surv[1:10,1] # Col 1 is males
H = s1$cumhaz[1:10,1]
Q = exp(-H) # Question: Is this the survival function?
cbind(H,S,Q)


#Check that Shat(t) for females is S0hat(t)^exp(beta1hat)
cbind(s1$surv[1:10,1],s1$surv[1:10,2],s1$surv[1:10,1]^exp(-1.483137252 ))


plot(s1)


# Try to locate the medians
xx = c(0,400); yy = c(.5,0.5)
lines(xx,yy,col="red")
# Median for M = 26, F = 141 ?

xm = c(26,26); ym = c(0,1); lines(xm,ym,col="red")
xf = c(141,141); yf = c(0,1); lines(xf,yf,col="blue")

# How about a nicer plot?
plot(s1,lty = c(1,2),xlab="Days", ylab="Probability")
title('Estimated "Survival" Probabilities for the Catheter')

xm = c(350,450); ym = c(0.9,0.9); lines(xm,ym,lty=1)
text(500,0.9,"Males  ")

xf = c(350,450); yf = c(0.8,0.8); lines(xf,yf,lty=2)
text(500,0.8,"Females")



# Compare disease types, just for women
table(Kidney$disease)

Other = data.frame(age=0, sex=1, disease="Other", frail=0)
GN = data.frame(age=0, sex=1, disease="GN", frail=0)
AN = data.frame(age=0, sex=1, disease="AN", frail=0)
PKD = data.frame(age=0, sex=1, disease="PKD", frail=0)
discomp = rbind(Other, GN, AN, PKD)
rownames(discomp) = c("Other", "GN", "AN", "PKD")
s2 = survfit(kmod1,newdata=discomp); s2

summary(kmod1)

# Catheters for patients with PKD stay in longest.
# How about some pairwise comparisons?

s2

plot(s2,lty = 1:4,xlab="Days", ylab="Probability")
title('Estimated "Survival" Probabilities by Disease Type')
### Which disease type has the obvious favourable prognosis?


### Compare the estimated survival curves using coxph and 
### the Breslow estimate of H0(t)
### VS the survival curves estimate from the KM approach
### Both approaches stratify by sex only
kmod1 = coxph( Surv(time,status) ~ sex, data=Kidney)
plot(survfit(kmod1,newdata=sexcomp))
# KM Curves for the two sex categories.
lines(survfit(Surv(time,status) ~ sex ,data=Kidney),lty=2)


kmod1 = coxph( Surv(time,status) ~ 1, data=Kidney)
plot(survfit(kmod1),conf.int = FALSE)
# KM Curves for the two sex categories.
lines(survfit(Surv(time,status) ~ 1 ,data=Kidney),lty=2,conf.int = F)
