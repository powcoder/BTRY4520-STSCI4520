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
#############      Homework 5    ########################
############# Due: May 4, 2018 ########################
#########################################################


# Instructions: save this file in the format <NetID>_HW5.R. 
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
# In this homework we have provided 

source('HW5Optimizers.R')

# implementing 1 dimensional Newton Raphson and a multi-dimensional
# Levenberg-Marquardt scheme if you wish to use them. 


#############################
# Question 1: Poisson LASSO #
#############################

# We saw previously that the LASSO for least-squares regression
# can be formulated as two steps, each applied to every 
# covariate until convergence:
#
# 1. Test is the minimum at zero (is the derivative of squared
#    error smaller than the derivative of absolute value at zero?)
# 2. If not, find the minimum of the penalized criterion. 
# 
# Here we will do the same thing but for Poisson regression. 
#
# For Poisson regression, we have a set of covariates X_i1,...,X_ip 
# and a response Y that takes integer values. If we assume that 
#
#     Y_i ~  Poisson( beta0 + beta1*X_i1 + ... + betap*X_ip)
# 
#  or in vector form  
#
#     Y_i ~ Poisson( beta0 + X_i beta )
#
# then we can choose beta0, beta to minimize the negative log likelihood
#
# l(b|Y,X) = sum  - log( P(Yi | Xi; beta0,beta) )
#         = sum  [- Y_i*(beta0 + Xi beta) + exp(beta0 + Xi beta) + log factorial(Yi)]
#
# where we drop the final factorial(Yi) because it doesn't change
# with beta. 
#
# Here we will consider trying to set some of the bj exactly 
# to zero by adding on a LASSO penalty and minimizing
#
#    pen.l(beta0,beta,lambda) = l(beta0,beta|X,Y) + lambda * sum_j |beta_j|
#
# However, we will not add a penalty for beta0. 



# For this problem, we will use data on school absenteeism 
# from Walgett, NSW -- a town with a large Aboriginal population
# that exhibits the social disfunctions that are common to 
# dispossessed indigenous groups around the world.  (The actual
# data is from 1978, but there are few statistical indicators
# that suggest much has changed in the past 40 years). 

Walgett = read.table('Walgett.csv',sep=',',head=TRUE)

# The data cover 148 children with columns
#
# Days -- number of days missed (our Poisson-distributed response)
# EthN -- Aboriginal or Not (0 = Aboriginal)
# SexM -- 0 = Female
# AgeF1 -- First grade indicator
# AgeF2 -- Second grade indicator
# AgeF3 -- Third grade indicator (note there is an F0 that is treated
#          as the default) 
# LrnSL -- Slow learner indicator
#
# Here we will try and work out which of these is predictive of
# the number of days missed. 
#
# It will be useful to take about 14 observation (10% of the data)
# at random. Here we will do this as 

test.ind = c(135,124,114,30,88,110,131,113,141,63,126,94,82,68)

# and we will create a set

test.Walgett = Walgett[test.ind,]

# and remove these from the original data

Walgett = Walgett[-test.ind,]

# a) Write out the derivative of l(b|X,Y) with respect to an 
# individual b_j. Hence determine a rule to decide whether the 
# optimum b_j (keeping all the others fixed) is at zero. 




# Hence write a function to provide a function to minimize 
# pen.l(b,lambda)  for a fixed b_j. You may use the NewtonRaphson
# code provided in HW5.Optimizers.R if you wish. 

pois.update = function(j,beta0,beta,lambda,X,Y){



}


# A useful starting point for beta0 is log(mean(Y)) with the
# remaining beta_j set to zero. Using the Walgett data set above
# examine which beta_j would move off zero with lambda = 100 (don't
# do this in sequence over j for this part of the question; for each
# beta_j leave all the others at zero). 
#
# Report the indices 

which.nonzero = 


# b) Just as in least squares LASSO, we can now perform
# coordinate descent by repeatedly optimizing each coordinte
# until convergence. 
#
# Write a function to carry this out for Poisson regression. 
#
# Remember that we don't penalize b0, but we do update it along
# with the others. A good starting point is that b0 = log(mean(Y)) 
# with all the other b_i being zero. 
#
# Your function should return an estimated beta vector (including
# the values that are zero) along with the number of optimization
# rounds, niter. 
#
# Remember to optimize beta0 in your update (but don't penalize it!)
# on each round through the co-ordinate descent algorithm: the initial
# guess at beta0 = log(mean(Y)) was just that -- a guess. You should
# also include beta0 in your convergence criterion. 
#
# Your function should allow a vector of lambda values, in which
# case you should return a vector of beta0's (one for each lambda)
# and a matrix of betas, where the rows of the matrix correspond
# to values of lambda. 


Pois.LASSO = function(X,Y,lambda,tol=1e-6,maxit=1000){



 return( list(beta0 = , beta = , niter = ) )

}

# Use this to obtain estimated beta for the sequence 

lambdas = seq(0,300,by=10)

# Produce a plot of the values of beta (without b0) as lambda
# changes. Put the value of each beta on the y axis
# with the x-axis given by  sum_(j >= 1) |beta_j|. You can plot
# all the betas on one plot. 



# d) We still need to choose lambda. One way to do that is to 
# use the test-set we left out in Walgett.test. Using your
# output from part c, calculate the likelihood for the predictions
# for the test set. That is obtain
#
#   sum_i  log P( Y_i| X_i beta )
#
# Write a function that takes in our output from c, and the test
# set and computes the log likelihood for each row of the 
# sequence of betas. Which coefficients are used in your final model?


test.likelihood = function(LASSO.result, testY, testX){

  

  return( loglik = )
} 





# BONUS: Rather than relying on one test set, or on leaving-one
# observation out at a time, 10-fold cross validation divides the 
# data into 10 approximately equal pieces. It then leaves one
# piece out in turn: estimates a sequence of betas using the remaining
# 9 pieces and obtains the log likelihood on the left-out set. 
# Then you average over each of the 10 left cross-validated likelihoods. 
#
# To do this, you should split the data at random. 
#
# Carry this out: does it change the value of lambda that you 
# select? 






######################################################################
# Question 2: Partially Linear Models + Uniform Confidence Intervals #
######################################################################

# This function will consider a variation on local linear regression.
# Here we will instead look at what is described as a varying
# coefficient model.

# To describe this model this, the data in 

ethanol = read.csv('ethanol.csv')

# which contains readings of Nitrous Oxide (NOx) produced by a 
# gasoline engine as well as the compression of the gas (C) in
# the engine and E -- a measure of the air/ethanol mix. 

# For these data, the relationship between NOx and C changes
# depending on E. If we examine

plot(ethanol$C[ethanol$E<0.8],ethanol$NOx[ethanol$E<0.8] )
plot(ethanol$C[ethanol$E>1],ethanol$NOx[ethanol$E>1] )

# We see somewhat different relationships. For that reason
# we consider a model of the form
#
#  NOx = beta0(E) + beta1(E) C
#
# That is, NOx has a linear relationship with C, but that 
# relationship depends in a non-parametric way on E. 

# a) We will be interested in being able to find beta0(e) and
# beta1(e) for any value e in the range of E. To do this, we 
# can perform a linear regression for  Y~C  (using the lm()
# function) with weights given by dnorm(e-E,sd=h) for bandwidth h. 
#
# Write a function to compute a varying coefficient estimate 
# of beta0(e) and beta1(e) at a given e and return the two-vector
# of coefficients. 

VarCoef = function(e,data,h)
{




}


# Use this to plot the estimated varying coefficient 
# estimates for the ethanol data over the range of 
# E, using h = 0.1.



# b) To choose a bandwidth, we can again apply cross
# validation. It's possible to represent our predicted
# values as a linear smoother  yhat = S(h)Y, but here
# you can just do it manually with a for loop. 
#
# Write a function to calculate the cross-validated
# squared error for predicting Y from these data. 
#
# In this case, unlike in local linear regression,
# we haven't centered by C and do need to use the whole
# model to predict each Y, rather than just an intercept. 


VarCoefCV = function(h,data){


}

# use this to choose an optimal bandwidth from among those in

hs = seq(0.01,0.2,by = 0.01)





# c) Fixing the bandwidth at the optimum you found in b,
# perform  a residual bootstrap for the values of
# beta0(e) and beta1(e) for each of the 61 values of e
# in

e = seq(0.6,1.2,by=0.01)

# (you should calculate all 61 values for each bootstrap
# re-sample). 
#
# Write a function to use this sequence e, the data 
# the bandwidth h, and produce a length(e)-by-3 
# matrices called conf.band.b0 and conf.band.b1 with 
# columns giving
#
# estimate - 1.96 boostrap.se,   estimate,  estimate+1.96 bootstrap se
#
# for each value in e where bootstrap.se is the standard
# deviation over nboot = 200 bootstrap samples. 
#
# You should also calculate the percentage of bootstrap 
# samples for which the bootstrap beta0(e) is inside the 
# confidence interval for ALL of the values of e (this is 
# like a multiple testing problem). 
 
VarCoefBoot = function(e,data,h,nboot)
{


  return( conf.band.b0 =  , conf.band.b1 = ,  
          numinside.b0 = ,  numinside.b1 = )

}

# Use the result to plot our estimated beta0(e) and beta1(e) with
# confidence bands around them.  Are these a good representation
# of the variability of your estimates?







# d) The intervals in Part c are described as being POINTWISE
# confidence intervals because they cover the value of beta0(e)
# (say) at each e. But it's possible that every individual 
# bootstrap curve goes outside the bands at some value of e, even
# if at any given e, 95% of the bootstrap curves are contained
# in the interval.  (eg, bootstrap 1 is outside at e[45], bootstrap 2
# at e[22] etc).

# Instead, UNIFORM confidence bands are intended to ensure that 
# 95% of bootstrap curves are completely contained within the band. 
# To produce these, we take an idea from Lab 5 and do the following. 
#
# 1.  From your bootstrap, calculate the t-statistics 
#
#  beta0.t(e) =  (beta0.boot(e) - beta0.est(e))/sd(beta0.boot(e))
#
# where beta0.est is the original estimate, and beta0.boot is 
# the estimate for each bootstrap sample. That is you should have
# length(e)-by-nboot t statistics. 
# 
# 2. Find the maximum value 
#
#    max.beta0.t = max |beta0.t(e)| 
#
# over e for each bootstrap to give  nboot values of  max.beta0.t.
#
# 3. Obtain Q as the 95th quantile of max.beta0.t -- that is 
# at most 5% of the curves are ever more that Q*sd(beta0.boot(e)) 
# away from beta0.est(e). 
#
# 4. Return the bands  
#
#    [beta0.est(e) - Q*sd(beta0.boot(e)),  beta0.est(e) + beta0.boot(e)]
#
#
# Write a function to repeat your calcuations in part c, but
# replace 1.96 with Q calculated as above.  


VarCoefBoot2 = function(e,data,h,nboot)
{


  return( conf.band.b0 =  , conf.band.b1 = ,  
          numinside.b0 = ,  numinside.b1 = )

}

# Do 95% of your bootstrap curves fall inside these bands?
# Produce a plot of your estimate and uniform confidence bands. 



# Small bonus: Why do we bootstrap t-statistics instead of just looking
# at a band of fixed width C (ie   [beta0(e) - C,  beta0(e) + C]
# with C chosen so that 95% of curves lie in this band?



# Big Bonus: Describe how to estimate this model with splines
# instead of kernels? What would you do to estimate the 
# partially linear model
#
# Y = beta0(E) + beta1*C?



#############################################
# Question 3: Log-Spline Density Estimation #
#############################################

# This question will examine maximum likelihood in a flexible 
# class of densities. To do this, we'll use a data set giving 
# the area burned in 270 forest fires. Because there is a lot
# of skew in this distribution, we'll look at the log of the 
# area. 

# This is given in the following data set; we've included a 
# histogram so you can get a sense of the shape of the distribution. 

fire = as.numeric(read.table('ForestFires.csv',head=FALSE)[[1]])

hist(fire,20,prob=TRUE)

# This looks normal-ish, but might be a bit skewed, and maybe 
# a bit light-tailed. We'll try to find a more flexible 
# set of densities to use to model this. 

# The idea is that the log of a normal density is quadratic
#
# log( dnorm(x) ) = [-mu^2/(2*sigma^2) - 1/2*log(2*pi*sigma)] 
#                         + (mu/sigma^2)*x - x^2/(2*sigma^2)   
#                 = a + b*x + c*x^2
# 
# Maybe we could add some additional terms. Say
#
# log( density(x)) = c0 + c1*x + c2*x^2 + c3*x^3 + c4*x^4

# Doing this runs us into a few difficulties: we know that 
# densities must be positive and integrate to one. To see how 
# this affects what we do, let's look at
# 
# density(x) = exp(c1*x + c2*x^2 + c3*x^3 + c4*x^4)/C
#
# Here it is natural to think of C = exp(-c0) as the number
# that makes the density integrate to 1. 

# Of course, we have to work out what C is. In this case, we
# will assume that the density is only defined on the interval
# [-3, 7].  Thus for any [c1,c2,c3,c4] we can obtain C by using 
# Simpson's rule to approximate
#
# int_(-3)^7 exp(c1*x + c2*x^2 + c3*x^3 + c4*x^4)

# a) Write a function to approximate C via Simpson's rule using
# a spacing of h = 0.25 between quadrature points. It should 
# taken in a vector c and return the number C

NormConst = function(c){
  
}

# Try this out on 

coefs = c(0.5,-0.125,0,0)

# to get 

Const1


# Use this to plot your (appropriately normalized) density overlayed on the
# histogram of log fire areas.  You should plot this at each of the 
# evaluation points used to evaluate the constant. 






# b) For real-world data, we need to maximize the log likelihood.
# We will use a nonlinear optimization algorithm to do this. 

# Recalling that the log likelihood is given as 
#
# sum_i log(density(x_i; c))
#
# write (separate) functions to evaluate the log likelihood along
# with its gradient with respect to c and its Hessian. Remember that
# your normalizing constant will change with C as well. 
#
# You may find the diag function useful, but note that if you wish to
# create a diagonal matrix with vector v on the diagonal, diag(v) only
# works if v is of mode vector. 

loglik = function(coefs,X){
   
  
}


ll.gradient = function(coefs,X){
   
 
}

ll.Hess = function(coefs,X){

}



# c) The file LevenbergMarquardt.R in HW5Optimizers.R implements 
# the Levenberg-Marquardt optimizer from Lecture 10. It assumes 
# inputs of the form above. 

# You may choose to write your own. In which case it should be pasted in here. 

# Use this to provide maximum likelihood estimates of your coefficients


opt.result = 


## PARTIAL CREDIT if you use the built-in function 'optim' to 
## optimize the log likelihood. (Only applies if your opt.result
## gives the wrong answer). You may complete the rest of the 
## questions with the coefs that you get this way. 
##
## Remember that 'optim' tries to minimize rather than maximize. 


# What is the new normalizing constant for these coefficients?

newconst = 


# Produce a plot of the density overlayed on the histogram of fire values
# using the same points as in Part a).  



# d)  Another direction the researchers consider is to have a more flexible
# model. As a way to manage this, we'll replace the 
#
# log(density(x)) = c_0 + c_1*x  + c2*x^2 + c3*x^3 + c4*x^4
#
# with a B-spline basis.  The idea is that with a set of functions 
# phi1(t),...,phiK(t) we can approximate a function
#
# f(t) = c1*phi1(t)+c2*phi2(t)+...+cK*phiK(t). Then we will write
#
#
# density(x) =  exp( f(x) )/ int_(-3)^7 exp(f(x))

# For this question we'll use some functions already coded
# up in the 'splines' library. You may need to install this library.  

library('splines')

# For this problem we will use 
#
#   - order 6 b-splines
#   - knots at the integers from -3 to 7
#
# Re-write your loglikelihood function, gradient and Hessian
# from part B to use the B-spline basis rather than the four
# terms x, x^2, x^3, x^4. 

# When you do this, you will need to drop the first B-spline.
# This is analogous to the way we removed the term c0 and 
# incorporated it into the constant.

# You should also write these so that they accept an input
# X. Because you may need to provide more input than just the data
# this can be a list of objects that you can then access
# inside your functions. Include whatever you think to be appropriate
# but aim to do as little extra computational work as possible.   

# We will assume whatever you are putting in X is now in an 
# object called

Xlist = 
	
# And produce the functions to evaluate the log likelihood
# gradient and hessian:
	
spline.loglik = function(coefs,X){
   
   
}


spline.ll.gr = function(coefs,X){
  
}


spline.ll.hess = function(coefs,X){
   
}

# We will check the evaluation of these using a simple set of coefficients

spline.coefs = rep(0,14)
spline.coefs[6:7] = 1


spline.val = spline.loglik(spline.coefs,Xlist)

spline.gr = spline.ll.gr(spline.coefs,Xlist)

spline.Hess = spline.ll.hess(spline.coefs,Xlist)


# e) Using the LevenbergMarquardt function provided, find the optimal values for these
# coefficients and plot the resulting density estimate along with the histogram (don't 
# forget to normalize).

# Note that you may find that LevenbergMarquardt does not converge in the default number
# of iterations. That is fine -- it gets us close enough. 

result2.opt = 





# f) We can penalize this estimate for smoothness, too. If we write
#
# density(x) =  exp(f(x))/int( f(x) ) then we can add a penalty
# to the log likelihood
#
# log likelihood =  sum_i log(density(x_i);c) - lambda*int (f'''(x))^2 dx
#
# here we penalize the third derivative. The reason for that is that 
# the log of a Gaussian is quadratic and its third derivative is 0. So large
# lambda will make the density look Gaussian. 
#
# Importantly, note that we SUBTRACT lambda*int (f'''(x))^2 dx; this is
# because we are maximizing the log likelihood rather than minimizing 
# a criterion. 

# Modify your functions above to include a this penalty.  

penspline.loglik = function(coefs,X){
   
   
}


penspline.ll.gr = function(coefs,X){
  
}


penspline.ll.hess = function(coefs,X){
   
}


# Verify that at lambda = 0, you get the same result as you did for the 
# un-penalized estimate.  Plot your result for lambda = 0.05 (no relationship
# to standard alpha levels). 




# BONUS: in part f, how would you go about choosing lambda?
# You do not need to implement this. 









############################################
# Question 4: Some Monte Carlo Integration #
############################################

# Based on JMR 18.6.18. You may not use the rcauchy, pcauchy or qcauchy functions 
# in this question. 

# The Cauchy distribution with parameter alpha has density
#
# f(x; alpha) = alpha/[pi*(alpha^2 + x^2)]
#
# with quantile (inverse CDF) function
#
# q(u;alpha) = alpha*tan(pi*(u - 0.5))

# a) Write a function to evaluate the inverse CDF for the Cauchy distribution. 

iCauchyCDF = function(u,alpha){


}


# Use this to sample 1,000 Cauchy-distributed samples with 
# alpha = 0.01 and produce a histogram of these with 100 bins. 
# Limit your x-axis to the range [-2, 2].
# 
# Overlay a plot of the Cauchy density to verify that you are 
# getting sensible answers.



# b) i) Write a function to estimate the expectation of
#
# g(x) = exp[(X-5)/50]/(1 + exp[(X-5)/50])
#
# when X is a Cauchy(2) random variable. g(x) can be calculated 
# with the plogis function. Use 1,000 Monte Carlo sample points.  
#
# Note that you may need to re-write g(x) to avoid numerical 
# overflow -- try not to exponentiate the positive value of a Cauchy 
# random variable. 

Exp.g.vanilla = function(N){


  return( estimate )
}




# ii) Repeat the evaluation using antithetic sampling. Here
# generate only N/2 random number so you evaluate the function
# a total of N times.

Exp.g.antithetic = function(N){


  return( estimate )
}

# iii) Estimate the variance of each of these estimates based 
# on R = 100 replicates and return the variance reduction (as a 
# percentage of the vanilla Monte Carlo estimate) due to antithetic
# sampling?

var.reduction = function(R,N){


   return( reduction )
}

# c) i) It is well-known that the Cauchy distribution does not 
# have a well-defined mean. We can see this in simulation by 
# estimating EX when X is Cauchy(2). Using your inverse sampler, 
# estimate the variance a vanilla Monte Carlo evaluation of X based 
# on sample size in a vector N. In this case we will use

Ns = c(10, 100, 10000)

# samples. Use 50 replicated Monte Carlo expectations 
# for each sample size. Store your variances in the 3-vector cvars


cauchy.mean.var = function(R,N){


   return( cvars =  )
}




# ii) Does the variance of your estimate decrease at rate 1/sqrt(n)? 
# Plot the relationship between your variance and n. 






