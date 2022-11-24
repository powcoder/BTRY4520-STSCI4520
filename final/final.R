#########################################################
#############   BTRY/STSCI 4520  ########################
#############      Final Exam    ########################
############# Due: May 20, 2018  ########################
#########################################################


# Instructions: save this file in the format <NetID>_Final.R. 
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


#### IMPORTANT INSTRUCTIONS FOR THE FINAL

## The final is to be completed as though it were in-class. 
## That means 
##
## 1. You must complete this work independently, without 
## collaboration or external assistance. Violations will
## be treated under the academic code. 
##
## 2. We will not provide office hours. We will monitor 
## Piazza and provide clarifications where questions are
## unclear. 
##
## 3. While we will not, in general, debug your work for you,
## we are happy to try to explain error messages if you can
## isolate the line of code that produces them and, preferably,
## reproduce them with a few-line script.
##
## [Example: R tries to interpret dimensions, but doesn't 
## always get it right.  So the code
##  
##   b = matrix(2,1,1)
##   t(1:3)/b
##
##  returns an error, but (1:3)/b does not.  We think it reasonable
##  to explain something like this to you if you get a message
##
##   Error in t(1:3)/b : non-conformable arrays
##
##  along with other instances of "I know where the error is, I just
##  don't know why R has a problem or how to fix it"]


################################################
# Question 1: Control Variates and Functionals #
################################################

# This question covers some recent ideas in Monte Carlo
# integration first published just two years ago in 
#
# https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/rssb.12185

# Recall that if we are interested in estimating
#
#  mu = E g(X)   where  X ~ f(x) 
#
# a Monte Carlo estimate is given by 
#
#  gbar =  (1/n) sum  g(Xi)   where  Xi ~ f(x)
#
# Control Variates use a function h(x) where we
# know that 
#
#  E h(X) = 0  if  X ~ f(x)
#
# (If we know Eh(X) = m, then use h*(x) = h(x)-m in place of h(x).)
#
# Then if h(x) looks like g(x) we can expect that  
#
#  hbar =  1/n sum h(Xi)    
#
# (which "should" be zero) might tell us something about
# how far  gbar is from mu.

# That is, if hbar is too high, it is likely that gbar is also
# too high.  In particular, we can improve gbar by
#
#  gbar2 =  gbar - alpha hbar
#
# Note that since hbar has mean zero, the expectation of gbar2 is
# still mu.  We found in class that the optimal alpha was 
#
#  alpha = cov(g(X),h(X))/var(h(X))
#
# which we can calculate from the g(Xi) and h(Xi) values. 

# a) Write a function to calculate Vanilla and Control Variate
# estimates of gbar from a sample X. It should take arguments
# X (a sample of data)  and functions g and h, and should
# return a 2-vector giving the Vanilla and Control Variate estimates.


MCfunc1 = function(X,g,h){
 
}

# b) i) Use this to write a function to conduct a simulation
# study in which you set g to be  sin(x+0.1) and h to be 
# -exp(x) + exp(1/2)  and generate X from a standard normal
# (the expectation of exp(X) is exp(1/2 when X is N(0,1)). 
# You should have length(X) be N = 100 and conduct R = 50 
# replications and return the percentage variance reduction 
# due to control variates

MCsim1 = function(N,R){
	

}


# ii) What variance reduction do you achieve?



# c) We can combine both control variates and antithetic
# sampling by looking at 
#
#  (1/n) sum (g(Xi) + g(-Xi))/2 + alpha (h(Xi) + h(-Xi))/2
#
# (note that if f(x) is symmetric about zero, X and -X
# have the same distribution). 

# i) Write a function to conduct a simulation as above 
# but with antithetic sampling, too.  You should return
# the Vanilla Monte Carlo variance and the variance
# reduction due to each of antithetic sampling, control
# variates and both. 
#
# It makes a big difference that you do the antithetic
# average before the control variate calculation.
#
# In this case use the same number N of random samples
# for antithetic sampling rather than fixing the number
# of function evaluations. 


MCsim2 = function(N,R){

	
  return( list( Vanilla.Var = , 
                Anti.Reduction = ,
                Control.Reduction = ,
                Both.Reduction =  ) )
}


# ii) Which adds most reduction? Does the joint reduction
# exceed what you would expect after observing 
# the reduction for each individual strategy?





# Bonus:  The original question had you set g(x) = sin(x)
# and h(x) = x  (we'll use this with Control Functionals below). 
# Why is this an uninformative choice of functions for the
# simulation you just conducted?




# d) The better the agreement between h and g, the more
# effective the control variates strategy will be.  Eg, 
# compare the h we used in parts b and c with using
# h(x) = x. 

# Oates et. al. 2016 propose to approximate g non-parametrically.
#
# Specifically, we'll let
#
#  h(x) = sum_j d_j psi_j(x)
#
# To do this, we need to know the integral of h(x). To accomplish
# this, the authors propose using a modified form of a kernel. 
# Specifically, they use the first X_1,...,X_(N/2) samples and 
# define  
#
#   psi_j(x) =  phi'(x-X_j;s) + phi(x-X_j;s) f'(x)/f(x)
#
# where phi(x;s) is a normal density with standard deviation s and ' 
# represents a derivative with respect to x. 
#
# The reason for this choice is that 
#
#  int psi_j(x)f(x)dx = int [phi'(x-X_j;s) f(x) + phi(x-X_j;s) f'(x)]dx
#                     = int [phi(x-X_j;s) f(x)]' dx 
#                     = phi(Inf,s)f(Inf) - phi(-Inf,s)f(-Inf)
#                     = 0
#
# With this choice of psi, we get the d_j by simple least-squares regression
# of  g(Xi) on  psi_j(Xi) for  i = 1,...,N/2. 
#
# To carry this out you should 
#    - evaluate psi_j(X_i) for for the X_i in your data
#    - estimate d_j by using the matrix of phi_j(X_i) to predict g(X_i)
#      using lm.  Some coefficients will be returned as NA -- set them to 0. 
#    - then use these d_j to give you sum_j d_j phi_j(x) for any x. 

# i) Write a function to evaluate this h(x) for a vector of values x. It 
# should have inputs
#
# x -- the vector of values at which to evaluate h(x) 
# X -- a vector of draws from f(x)
# s -- the bandwidth
# g -- the function we want to approximate
#
# In this case you can assume that f(x) is a standard normal. It should
# return a vector of values of length(x).  You may use 

hfunc = function(x,X,s,g){

  
}

# ii) Using g(x) = sin(x), and s = 0.2, set 

xvec = seq(-2,2,by = 0.1)

# iii) plot h(x) based on 100 standard normals. 




# e) i) Write a function to implement
#
#  1.  Vanilla Monte Carlo
#  2.  Control Variates, using h(x) = x
#  3.  Control Functionals
#        - Use the first half of the values X_1,...,X_N
#          to produce hfunc as in 1d
#        - Use the second half to obtain the control variate
#          estimate with this hfunc
# 
# These should be based on X being N(0,1)
#
# Your function should take arguments N, g and s and return
# the three estimates above as a 3-vector

ControlFunctionals = function(N,g,s){
  
}

# ii) When g(x) = sin(x) we know that E g(X) = 0. Based on 
# 50 replicates at N = 100, using s = 0.2, what is the 
# expected squared error of each of your estimates?





# f) One of the claims of the paper is that we can 
# achieve a faster convergence rate than (1/sqrt(N))
# because as we get more data, h(x) gets closer to g(x). 
# 
# i) Write a function to conduct a simulation study so that
# for each value of 

Ns = c(100,200,500,1000)

# you estimate the root-average-squared-error (RMSE) of
# each of the three estimates in 1e based on R = 20 
# simulations using values of s in 

svec = 2/sqrt(Ns)

# that correspond to each value in Ns.  (Here we have
# smaller s so that we can approximate g better when
# we have more data points). 

# You should return the matrix of mean squared errors.

# CAUTION: this may take a few minutes to run. 

ControlFuncSim = function(Ns,R,s){

}

# ii) Plot log(RMSE) versus log(N) for the three estimates.
# This should give you approximately a straight line with
# the order of convergence given by the slope. 





# iii) Confirm that vanilla and control variates approaches have
# MSE = O(1/Ns); what is the order for Control Functionals?








####################################################
# Question 2: Kernel Density Estimation and        # 
#             Minimum Hellinger Distance estimates #
####################################################

# In this question we will use Kernel Density Estimation
# as an intermediate step in statistical procedures. 
#
# But we will first need to deal with KDE's. Recall that
# we can estimate a nonparametric density by 
#
#  fhat(x,s) =  (1/n) sum phi(x-Xi;s)
#
# (using the same notation as in Question 1). 

# a) i) Write a function to evaluate fhat for a vector
# of evaluation points x, using data X with bandwidth s. 
# You should not need to use either for loops or apply 
# statements. 

kde = function(x,X,s){
   

}

# The data we will use for this question comes from a study of
# the effectiveness of a drug pyrantel used to treat parasites
# in horses.  These data record the logit of the ratio of 
# the number of parasite eggs found in horse feces before
# and after treatment.

eggrate = c(-1.1653, -0.7538, -1.3218, -2.3394, -1.9766, -1.8718, -1.5041)

# ii) Using s = 0.2 produce a plot of the kernel density estimate
# on values of x from -2.5 to -0.5.





# b) In order to apply this well, we need to choose a 
# bandwidth. To do this we'll use cross-validation. 
# Specifically we'll define a score for each h as
#
# CV(s) =  sum log  fhat^(-i)(Xi;s)
#
# that is. For each Xi, we leave Xi out of the data
# and estimate fhat^(-i) without Xi, then we see how high
# fhat^(-i) thinks the density is at Xi.  The higher
# the density, the better we are at predicting where future
# data might be.  Hnce the optimal s is the one that maximizes
# CV(s).
#
# i) Write a function to evaluate CV(s) for each s. Bonus 2 points
# for avoiding for loops and apply statements. 

kde.cv = function(X,s){

}


# ii) Use this function to decide on the optimal s for the eggrate
# data above with possible bandwidths given in 

ss = seq(0.1,1,by=0.1)

# to give the value sbest that maximizes the cross-validated
# likelihood. 

sbest = 




# c) We will also need to be able to simulate data from fhat. 
# To do this, we observe that we can express fhat as being exactly
# the density that corresponds to the following
#
#   1. Choose one of the Xi at random
#   2. Simulate a new z from a N(Xi,s) distribution. 
#
# i) Write a function to simulate N observations from fhat given 
# X and s. 2 bonus points for avoiding for loops 

kde.sim = function(N,X,s){

}

# ii) Use this to simulate 1000 data points from your estimate
# above (with optimal bandwidth). Draw a histogram of these points
# and add the original data as vertical red lines. 








# d) One of the advantages of having a density estimate is
# that there are more ways that you can compare your data to 
# a parametric family of models. 
#
# In particular, one measure of how different two distributions
# are is Hellinger distance
#
#    HD(f,g) = int ( sqrt(f(x)) - sqrt(g(x)) )^2 dx
#            = 2 - 2 int sqrt(f(x)*g(x)) dx
#
# For a family of densities f(x,theta) and an estimate 
# fhat(x), the minimum Hellinger Distance estimate (MHDE)
# is obtained by maximizing
#
#    A(theta) =  \int sqrt( f(x,theta)*fhat(x) ) dx
#
# In our case, we will use a normal family for the mean
#
#  f(x,theta) = dnorm(x,mean=theta)
#
# and our KDE. 
#
# We still need to approximate the integral. Here we will 
# sample from fhat and use the approximation
#
#  A(theta,fhat) = int sqrt( f(x,theta)/fhat(x) ) fhat(x) dx
#           ~=  (1/N) sum  sqrt( f(Xi,theta)/fhat(Xi) )
#
# where the Xi are sampled from fhat as in part c. 

# i) Write a function to calculate A(theta), for f(x,theta) given
# by dnorm(x,mean=theta) using the original data X and a sample
# Xsamp generated from fhat with bandwidth s. 



Afn = function(theta,X,Xsamp,s){

}

# ii) And a further function to find the optimum value based on 
# Data X, bandwidth s and sample size N. 
#
# You may use any optimization function you like (Giles
# used 'optimize') and can assume the maximum lies in 
# (-10,10). 


HD.opt = function(X,s,N){

}


# iii) Use this to obtain the value that minimizes Hellinger
# distance for the data using the optimal bandwidth
# above based on 1000 Monte Carlo samples. 


hellinger.opt = 




# e) i) One reason to look for the MHDE is robustness
# to outliers. To see this add an additional data point
# to your data set with values in 

O = seq(-1,10,by=0.2)

# plot the value of the MHDE and the average 
# value in your data versus the values in O above. Keep
# the optimal bandwidth you calculated in part c. 



# ii) Which of these (average or MHDE) do you think is 
# more reasonable?







# f) One of the reasons that Hellinger distance has
# received some attention is that it is supposed to be
# as precise as the mean, at least asymptotically. 
#
# However, for a given sample, it can be difficult to
# work out how precise it is.  One suggestion was
# to use it within a Bayesian analysis. See
#
#    https://link.springer.com/article/10.1007%2Fs11749-014-0360-z
#
# The idea is to replace the posterior with
#
#    exp( n*A(theta) ) pi(theta)
#
# of pi(theta) is the prior.  
#
# i) Using a N(0,10) prior, write a function to conduct a random 
# walk MCMC algorithm. The function should take in
#
#    X -- your data
#    s -- bandwidth for the KDE
#    N -- number of Monte Carlo samples to evaluate A
#    nmc -- length of the MCMC chain
#    sig -- random walk variance
#
# and return a set of posterior samples and the acceptance
# rate. You should keep the same Monte Carlo samples for
# each evaluation. 


HD.MCMC = function(X,s,N,nmc,sig){
 
 
 
 return( list( samples = ,  acc =  ) )
}

# ii) Experimentally choose sig to give you between 30% and 40% 
# acceptance rate and produce a histogram of the posterior 
# based on 1000 samples after 1000 samples burn-in and thinning
# to every 5th sample.  




# Bonus: in fact, MCMC can be used with a stochastic likelihood.
# That is, it still works if you use new Monte Carlo samples
# each time.  Apply this using N = 100, but draw new samples
# for each step in the MCMC. How high can you get your 
# acceptance rate?



HD.MCMC2 = function(X,s,nsamples,nmc,sig){
 

 return( list( samples = ,  acc =  ) )
}





######################################
# Question 3: Random Effects Methods #
######################################

# In the additional exercises after Lecture 18, we
# saw how Poisson random effects models can be fit
# with MCMC. Here we'll look at extensions of
# this and Frequentist alternatives. 

# In this question we will consider a logistic model
# with random effects. That is, we observe outcomes that
# are either 0 or 1 and write out
#
#  P(Y=1|X,Z) = plogis( beta0 + beta1*x + z)
#
# where x is a covariate and z is a subject-specific
# random effect. 
#
# We can turn this into an over-all probability by
#
# P(Y|X,z) = P(Y=1|X,Z)^Y*(1-P(Y=1|X,Z))^(1-Y)
#
# which evaluates to P(Y=1|X,Z) when Y = 1 and 
# (1-P(Y=1|X,Z)) when Y = 0. 

# Here we regard Z as a latent normal random variable that
# applies to several data points. By way of example, the
# data in 

toenail = read.table('toenail.txt',head=TRUE)
toenail$Subject = as.factor(toenail$Subject)

# gives an indicator of nail health (good or bad) for 
# treatment for foot fungus. This is recorded for 7 visits each
# for 12 Subjects at different times given by number of Months. 
#
# Here we will write Y_ij for the i-th measurement from the j-th 
# patient and say
#
#  P(Y_ij=1|X_ij,Z_j,beta0,beta1) = plogis( beta0 + beta1*X_ij + Z_j)
#
# where $X_ij$ is the time of the i-th visit of subject j. Z_j
# is a different value for each subject so that some subjects 
# have consistently worse feet than others depending on the size 
# of Z_j. 
#
# We don't get to see the Z_j, so we will specify that 
#
#   Z_j  ~ N(0, \sigma^2) 
#
# The likelihood of the observed data in group j is
# then 
#
# P(Y_1j,..,Y_7j|X_j) = \int prod( P(Y_ij|X_ij,Z_j,beta0,beta) ) phi(Z_j;0,sigma)dZ_j
#
# and the likelihood for the parameters beta0, beta1, sigma is
#
# l(beta0,beta1,sigma) = sum_j log P(Y_1j,...,Y_7j|X_j)


# a) When only one random effect applies to each notation we can
# evaluate the integral above numerically.  To do this, we will
# use a Gauss-Hermite approximation. This can be obtained through 
# the ecoreg package  (you will need to install it). 

library(ecoreg)

# The function gauss.hermite produce Gauss-Hermite points

gh = gauss.hermite(21)

# Where we can approximate  \int f(x) phi(x) dx by 
#
#   sum(gh[,'Weights']*f(gh[,'Points'])
#
# Note that this is for a standard N(0,1), to make this N(0,s^2)
# you need to multiply gh[,'Points] by s. 

# i) Write a function to evaluate the negative of the log likelihood 
# of the toenail data for any beta0, beta1 and sigma given 
# the vector theta = c(beta0,beta1,sigma)

logistic.nll = function(theta,data){

	
}	
	


# ii) Hence find the values of beta0, beta1 and sigma that
# maximize the log likelihood. You may use optim or any other
# optimizer that you wish. 

# Good starting values are 

theta = c(-0.2,-0.2,1.2)   







# b) As an alternative, we can use Monte Carlo integration for 
# each random effect. To do this, we'll replace the Gauss-Hermite
# points and weights with a sample of 1000 N(0,1) random values
# and weights 1/1000.  To change this to integrating with 
# respect to a N(0,s^2) you can make the same changes as for 
# Gauss-Hermite quadrature.  Use the same 1000 points for each
# group. 
#
# i) Write a function to calculate the negative log likelihood based 
# on this approximation. You may hard-code the use of N = 1000 
# Monte Carlo points. 

logistic.nll.mc = function(theta,data){
  	
}

# ii) How would you go about optimizing this likelihood? Do so. 






# c) Given a likelihood, we can always carrying out an MCMC
# analysis to find a posterior. In this case we will
# fix sigma to be 1.2 to make the computation simpler. 
# Setting your prior on beta0 and beta1 to be 
#
# P(beta0,beta1) = N(0,10)*N(0,10)
#
#
# i) Carry out a random walk MCMC for 5000 steps using the 
# likelihood you found in part a using a standard deviation
# of

sds = c(1,0.2) 

# for the steps in each dimension of your random walk. Use the same
# initial value for beta0 and beta1 as in part a. 




# ii) How far is your expected a-posteriori estimator from the maximum 
# likelihood estimate in part a?




# iii) Produce 95% credible intervals for beta0 and beta1












# d) Rather than doing the integration numerically, MCMC also
# allows you to simply include the random effects as additional 
# parameters. That is, we can write the joint distribution of
# all of the data, the Z's and and beta0 and beta1 (still keeping
# sigma fixed at 1.2) as 
#
#  P(Y_ij|X,Z_j,beta0,beta1)P(Z_j)P(beta0)P(beta1)
#
# And conducting a random walk on beta0,beta1,Z_1,..,Z_K

# i) Write a function to calculate this log posterior using the same
# prior as as in part c. 

joint.posterior = function(theta,Z,data){

}

# and carry out a random walk MCMC using step sizes with variance

sds2 = sds/4

# for theta and steps of standard deviation 0.5 for Z.
# Initialize Z to be all zero. 





# ii) Have your expected a posteriori estimates for bet0
# and beta1 moved outside the credible interval found 
# in part d?



# iii) Plot a histogram of the posterior distribution of Z_1
# based on every 5th draw from the chain after dropping the first 
# 1000 draws






# e) i) Can we calculate the conditional distribution of Z_1
# as a frequentist? Here is one last way of generating random
# variables that takes inspiration from importance sampling
# in Monte Carlo integration. 
#
# Specifically, suppose we want to sample from g(x), but can
# readily draw from f(x). The idea is to draw X_1,...,X_N but 
# to then re-sample these X's with weights 
#
#  W_i = g(X_i)/f(X_i)
#
# (the sample() function allows you to specify a vector 'prob'
# of weights). This means that the density of the resulting
# X's is 
#
#   P(X selected)*P(X in dx) = W f(x) = g(x)/f(x) * f(x) = g(x)
#
# To get a posterior distribution of Z1, first, use the beta0,
# beta1 and sigma from part a.  Then 
#
#   1. Generate 1000 samples of Z from N(0,sigma^2) (you may use rnorm)
#   2. Resample the Z with replacement, but with weights given by 
#        prod(  P(Y_i1 | X_i1, Z,beta0,beta1) )
#

# Write a function to obtain this samples using beta0, beta1
# and sigma calculated in part a but using the data in 


Subj1 = toenail[1:7,]

ranef.post = function(N,theta,Subj1){

}



# ii) Produce a histogram based on your 1000 resulting samples. How does
# this compare to the posterior distribution in part b?





# iii) This approach to generating random numbers is known as sequential
# Monte Carlo (SMC) and is particularly useful when you want to keep
# updating distributions -- you can take a sample from f_1(x) to f_2(x) 
# and then update again to f_3(x).  This arises, for example, if you
# are tracking a noisily observed process. 
#
# Here, however, we are only generating from one distribution. Is it a 
# good or a bad idea? Why?





# f)  Describe how you would go about constructing a bootstrap distribution
# for this model, accounting for the fact that you might have different
# subjects next time.  You do not need to carry this out. 





#######################
# BONUS: Optimization #
#######################

# This repeats the final exercise in Lab 7 (which the labs
# did not get to).  Additional 10 marks for an answer within
# 10% of the true values


# The data in FhN.csv come from an example that Giles was 
# working on for his most recent book. They arise from
# a model of the way neurons transmit signals and are given
# by an ordinary differential equation. 
#
# This is a simplification of a more complex model than won 
# Alan Hodgkin and Andrew Huxley a Nobel prize in 1963 (so
# it's fitting that Giles is writing this lab while in 
# Stockholm). 
#
# The data is in 

FHN = read.csv('fhn.csv',head=FALSE)

# It contains a set of times

t = FHN[,1]

# and two-columns of values

V = FHN[,2:3]

# so we can look at

matplot(t,V,type='l')

# These values were produced using the following function, 
# which relies on the deSolve package to solve ODE equations

library(deSolve)


FHN.fn = function(t,p){

	x0 = p[4:5]
	
	fhnfunode = function(t,x,p)
	{
	r = x;
	r[1] = p[3]*(x[1] - x[1]^3/3 + x[2])
	r[2] = -(x[1] -p[1] + p[2]*x[2])/p[3]
	return(list(dx=r))
	}
	
	res = lsoda(x0,t,fhnfunode,p)

	return( res[,2:3] )
}

# This depends on a vector p of 5 parameters, 
# for instance, if you set

p = c(1,1,3,0.5,-1.2)

# you get 

matplot(t, FHN.fn(t,p),type='l')

# The game here is to find the values of p that generated
# the data. There is no noise in the data, so the right
# p should produce V exactly. 

# To do this, use the optim function in R, and experiment
# with different optimization solvers and/or starting points. 
# You will need to define an objective function -- I recommend
# squared error and trying to find where it is zero. 

# Two things might help:
#
#  1. I will tell you the true parameters are in the ranges
#      [-1,1],  [-1,1],   [0.5,5],  [-2,2],  [-2, 2]
#
#  2. You will find that at bad values of p, FHN.fn will
#  produce an error. This will cause your optimization to
#  to stop. However you can get around this with the
#  try-catch construction. 
#
#   Here's an example

for(i in 1:3){
  try({
      m = b/8 
	  })
  print(i)
}

# Although there is an error (we have never defined b) the for 
# loop does not terminate.  You can set a number of lines of 
# code inside the {} in try.   
#
# You might use this to try evaluating at some p, for example,
# and return Inf, or a really big number if that p produces
# an error. 
#
# You may search for good starting values, but describe how you
# decided to start from there, and include code to get from them
# to your final estimate. 



