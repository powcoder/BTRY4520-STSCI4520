https://powcoder.com
代写代考加微信 powcoder
Assignment Project Exam Help
Add WeChat powcoder
https://powcoder.com
代写代考加微信 powcoder
Assignment Project Exam Help
Add WeChat powcoder
#########################################################
#############   BTRY/STSCI 4520  ########################
#############      Homework 4    ########################
############# Due: April 13, 2018 ########################
#########################################################


# Instructions: save this file in the format <NetID>_HW4.R. 
# Complete each question using code below the question number.
# You need only upload this file to CMS. 

# Note, we assume your working directory contains any files
# that accompany this one. 

# Further note: 10% will be deducted if your file produces
# an error when run. If your code produces an error and you 
# cannot find it, comment out that portion of the code and we
# will give partial credit for it. 

# Do not use either the function set.seed() or rm(list=ls())
# in your code. 

# The file 1DOptimizers.r contains the functions 
# GoldenSection and NewtonRaphson from Lecture 11. 
# If you wish to use them, you can load these up with

source('1DOptimizers.R')

# if you wish to use them.  DO NOT modify these -- when
# we look at your code, we will have this file in our
# working directory as it currently is. We will not 
# replace it with a modified version. 

# A note on these, in many cases you will want to optimize
# functions that have more than one argument, but neither
# NewtonRaphson nor GoldenSection have the option to add
# more inputs.  

# You can either choose to re-write these functions, or 
# an alternative is to define a new function that inherits 
# the other inputs as constants. 

# Eg in 1a below, within Mixed.GS, I can define a function
#
# newfn = function(x){ fn(x,mean1,mean2,sd1,sd2) }
#
# to use within GoldenSection. Here newfn will not see 
# mean1, mean2, sd1, sd2, but will look in the environment
# in Mixed.GS to find these quantities, rather than in the
# global environment. 

# These sort of issues come up all the time when using
# code that others have written. Not allowing additional inputs
# is poor form (sorry!), but there are many other times
# when you have to find ways around in-built (often unintentional)
# restrictions. 


######################################
# Question 1: Domains of Convergence #
######################################


# In class we saw that even finding the mode (the 
# most probable value) of a mixture distribution can
# be difficult. In particular, we looked at 
#
# f(x) = 0.5*dnorm(x, mean = 0, sd = 1) + 0.5*dnorm(x,mean=2,sd=2)
#
# Here we will observe that transforming f can have
# an important effect. 

# In particular, we will examine what happens when we
# try to optimize f(x), versus log( f(x) ) when we 
# examine different starting points

# a) Write a function to find the maximum of f(x) 
# using Golden Section search starting from an interval 
# of [-m,m]. Your function should have inputs that are 
#
#  - mean1, mean2: the means of the two components
#  - sd1, sd2: the standard deviations of the two components
#  - m - to define the search interval [-m, m]
#  - tol - the tolerance criteria
#  - maxit -- maximum number of iterations
#
# It should return the optimizing value of x in argmax,
# the value of f in value, and the number of iterations.


Mixed.GS = function(mean1,mean2,sd1,sd2,m,tol=1e-6,maxit=1000)
{

 return(list( argmax = , value = , niter = )
}

# Using the values

m = 2:10

# along with settings

mean1 = 0
mean2 = 2
sd1 = 1
sd2 = 2

# record the argmax values found and the number of iterations taken

opt.x.gs = 
niter.gs = 


# b) We will compare this to using Newton-Raphson. Write a function
# to obtain the NewtonRaphson estimate starting from m with 
# the same inputs and outputs as part a:

Mixed.NR = function(mean1,mean2,sd1,sd2,m,tol=1e-6,maxit=1000)
{

 return(list( argmax = , value = , niter = )
}

# Again, store the values that you found but this time for

m2 = -10:10

# along with the means and variances above. 

opt.x.nr = 
niter.nr = 

# c) Repeat a and b using log(f) instead of f

Mixed.GS.log = function(mean1,mean2,sd1,sd2,m,tol=1e-6,maxit=1000)
{

 return(list( argmax = , value = , niter = )
}


Mixed.NR.log = function(mean1,mean2,sd1,sd2,m,tol=1e-6,maxit=1000)
{

 return(list( argmax = , value = , niter = )
}


# How do the convergence properties change? Why?

opt.x.gs.log = 
niter.gs.log = 

opt.x.nr.log = 
niter.nr.log = 



# d) Repeat the experiments above but with 
# params = (mean1,mean2,sd1,sd2) given by

params2 = c(-2,2,2,1.5)

# to find

opt.x.gs2 = 
niter.gs2 = 

opt.x.nr2 = 
niter.nr2 = 

opt.x.gs.log2 = 
niter.gs.log2 = 

opt.x.nr.log2 = 
niter.nr.log2 = 

# and

params3 = c(-2,2,1,0.8)

# to find

opt.x.gs3 = 
niter.gs3 = 

opt.x.nr3 = 
niter.nr3 = 

opt.x.gs.log3 = 
niter.gs.log3 = 

opt.x.nr.log3 = 
niter.nr.log3 = 


# Describe how your results change in this case. Why?





############################
# Question 2: Optimization #
############################

# In particle physics, when muons decay, they emit electrons
# in a direction with distribution relative to their orientation. 
# This distribution has been calculated to have density
#
# f(x;alpha) = (1+alpha x)/2
#
# where x is cos(theta) and theta is the angle of the direction, so
# x lies in [-1, 1]. Here alpha also lies in [-1,1].
#
# A sample of such x's is in 

muon = read.table('muon.txt')$V1


# A standard way to estimate alpha is to maximize the log 
# likelihood. That is, we seek to maximize
#
# l(alpha) = \sum_i log[ (1 + alpha x_i)/2 ]
#
# a) Write a function to apply a Golden section search 
# (you may copy from notes) to find the maximum likelihood estimate
# for alpha, up to some tolerance.  
#
# It should return
#
# - opt.alpha: the value that maximizes the log likelihood
# - value: the value of the log likelihood
# - niter: the number of iterations taken


muon.golden.section = function(data, tol=1e-6,maxit=1000)
{



  return( list( opt.alpha = , value = , niter = ) )
}



# b) Write a function to provide a bias-corrected normal-theory
# bootstrap confidence interval for your estimate. Test this using 
# nboot = 200 bootstrap samples. You should return the bias-corrected
# estimate and the confidence interval

boot.muon.ci = function(data,tol=1e-6,nboot=200)
{


 return( list(corrected.estimate = , confint = ) )
}


# c) If amax is the value of alpha that maximizes l(alpha), an 
# alternative confidence interval is to find values alpha such
# that 
#
# l(amax) - l(alpha) = 3.84
#
# there should be one point  smaller than amax and one larger than it. 
# To obtain these we will try a secant method. Since amax should
# be at the top of the likelihood function, we should be able to 
# head left or right to start heading downhill. 

# Your function should implement the secant method (you may code
# this up as a separate function) and then start from initial 
# points amax and -1 to head left, and amax and 1
# to head right. 

# When you carry out the iteration, you should make sure that
# you never have an estimate outside the interval [-1,1], by
# taking any such value back to the nearest endpoint. 

# Return the confidence interval.

likelihood.ratio.ci = function(data,amax,tol=1e-6)
{


  return(confint)
}




#########################################################
# Question 3: Rejection Sampling for a Truncated Normal #
#########################################################

# In this question we will look at generating data at one
# end of a normal distribution. Specifically, we consider
# trying to simulate X as N(0,1), conditional on X > a
# where a is positive. Formally, we can write the density in 
# R as
#
# f(x) =  dnorm(x)*(x>a)/(1-pnorm(a))
#
# This can be quite challenging when a is very large as 
# the approximations in pnorm and in qnorm become unstable. 


# a) The first option is to simply generate normal data and
# throw away the data for which x < a.  Write a function
# to generate n samples this way. Return the samples X and the 
# number of random numbers N needed to generate them.

rtruncnorm1 = function(n,a){

  return(list( X = , N = ) )
}

# What percentage of samples do  you keep? (Bonus 1 point
# for being as efficient as possible. 




# b) An alternative is to perform rejection sampling. To
# do this, we will use a shifted exponential distribution. 
# That is, we will look at Z=Y+a where Y ~ Exp(r). Since
# the Y's are all positive, Z is necessarily greater 
# than a.

# The density of Z can be written down as 
#
#  g(z) = r*exp(-r*(z-a))
# 
# and we can obtain a Z from a uniform U by
#
# Z = a - log(1-U)/r

# An important note about rejection sampling; we only need 
# to know the form of a density, not the constants. That is
# if we can generate
#
# X ~ g(x)
# 
# and want Y ~ C*f(x) where C is a number we can't calculate
# easily, we can still perform rejection sampling to generate
# y by finding k such that  k*g(x) > f(x) for all x. Then 
# the recipe is the same:
#
#  1. Simulate X ~ g(x)
#  2. Simulate U ~ [0, k*g(x)]
#  3. Keep X if  U < f(x)

# i) For a fixed a and r, find the smallest value of k such that
# 
#  k * dexp(x-a,r) > dnorm(x) 
#
# for all x > a (the x>a restriction may be important).    You may 
# submit scanned hand-calculations if you wish, or give the math 
# in comments below. 



# ii) For a fixed value of a, what is the optimal exponential rate 
# r to choose? What value of k does this yield?


# c) Write a function to derive truncated normal random
# variables using this rejection sampler. As in part a, return
# n samples as well as a count of the number of random numbers
# used.  

rtruncnorm2 = function(n,a){

  return(list( X = , N = ) )
}



# d) Plot the effort to generate 1000 samples for a = 1,2,3,4  
# for both methods above along with the expected effort for 
# rtruncnorm1. 
#
# Store the effort in a 4-by-3 table

effort = 



# BONUS: Find a way to generate truncated normals with no 
# waste, assuming you can use qnorm and pnorm. How does this
# perform at a = 2? Compare all three methods at a=10. 




##################################
# Question 4: Poisson Regression #
##################################

# The Poisson distribution is often used for modeling count
# data (eg, the number of trucks going past Giles' house between
# 5 and 6 in the morning). It has density parameterized by the 
# rate r with formula:
#
# P(Y = k; r) = r^k * exp(-r)/factorial(k)
#
# Here we will examine relating a count to a covariate (day of 
# the week in the case of Giles' ability to sleep in, but see
# examples below).  
#
# The idea is to model the Poisson rate as changing with X and
# in particular
#
# log(r) = b0 + b1 X
#
# So that each Yi has its own rate Ri that is determined 
# by Xi.   This means that the probability that we see 
# the data Yi at covariate Xi is
#
# P(Y = Yi| Xi) = exp(Yi*(b0+b1*Xi)) * exp(-exp(b0+b1*Xi))/factorial(Yi)
#
# so that the likelihood is the product of these over i. 
#
# Here, we will maximize the log likelihood (ie, the sum of the
# log probabilities)
#
# l(b0,b1|Y,X) = sum log( P(Yi | Xi; b0,b1) )
#         = sum  Y_i*(b0 + b1*Xi) - exp(b0 + b1*Xi) - log factorial(Yi)
#
# Since the last term does not change with b0 and b1, we usually drop it
# and just try to maximize
#
# l(b0,b1|Y,X) = sum  Y_i*(b0 + b1*Xi) - exp(b0 + b1*Xi)

# As an example, the data in 

Flu = read.csv('Influenza.txt')

# gives the number of Influenza cases reported for high school
# students along with the number of days since the outbreak. 
#
# We will model Students as Y and Days as X above.  

# You may want to write functions evaluating the first derivative 
# vector and matrix of second derviatives to use in parts 
# a and b below. 

# a) Write a co-ordinate ascent algorithm to maximize l(b0,b1|Y,X)
# using a one-dimensional Newton-Raphson method by first maximizing
# for b0, then b1, then b0 and so on until convergence. You may copy 
# the Newton-Raphson algorithm from notes if you wish. 

# Your function should return the optimum values of b0 and b1
# along with the value of the likelihood at the optimum and the 
# number of iterations taken to achieve it. How many iterations
# did you need to fit the Influenza data?

# Start from b0=0, b1=0


PoissonReg1 = function(b0,b1,X,Y,tol=1e-6,maxit=1000)
{

	return( list( b0 = , b1 = , likelihood = , num.iter = ) )
}


# b) Compare this to applying a two-dimensional Newton-Raphson 
# method to update b0 and b1 simultaneously.  Write a function to 
# perform this method and report the same quantities along with the 
# iteration history. 


PoissonReg2 = function(beta,X,Y,tol=1e-6,maxit=1000)
{

	return( list( beta = , likelihood = , num.iter = , iterhist = ) )
}

# Produce a contour plot of the likelihood over the range of the iteration 
# history and add the path taken by the Newton steps. To view this 
# reasonably, you should set any likelihood value less than -800 to be 800. 


# C) Carry out a parametric bootstrap based on your fitted parameters
# (you can generate Poisson random variables with the rpois() function).



# Provide confidence intervals for b0, b1 and b0+b1*5 -- the last of these
# is the log of the expected count at 5 days.  Examine the bootstrap
# distribution for exp(b0+b1*5) -- would a symmetric confidence interval
# be appropriate here?

Flu.confint.b0 = 
Flu.confint.b1 = 
Flu.confint.5 = 




# BONUS: In the linear regression function lm, you can
# specify a vector of weights as an input. If you do,
# the coefficients are calculated to minimize
#
#  sum  Wi (Yi - b0 - b1 Xi)^2
#
# this minimum can be calculated in matrix terms to 
# satisfy (in R code):
#
# b = solve(t(X) %*% diag(W) %*% X,  t(X) %*% diag(W) %*% Y)
#
# Re-write your function in part b to repeatedly use the
# lm function and weights (which change each iteration). 

# This is one example of a general procedure called
# iteratively reweighted least squares (IRWLS) that is used
# throughout Generalized Linear Models. 